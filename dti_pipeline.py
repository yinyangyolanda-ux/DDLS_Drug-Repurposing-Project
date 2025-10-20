import os
import random
import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from chembl_webresource_client.new_client import new_client
from DeepPurpose import utils, DTI
from oddt.metrics import enrichment_factor

# -----------------------------------------------------------
# SECTION 0: Enforce Determinism and Setup
# -----------------------------------------------------------
SEED = 42

def set_seed(seed):
    """Set seeds for reproducibility across all libraries."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

set_seed(SEED)
print(f"Random seed set to {SEED} for full pipeline determinism.") [2]

# Placeholder for complex GDSC feature extraction/distillation
def load_and_distill_gdsc_features():
    """Simulate loading distilled GDSC genomic features (VAE/Autoencoder step)."""
    # This step represents reducing 57k+ omics features to a dense, stable vector [3]
    print("Simulating distillation of GDSC multi-omics features...") 
    return "GDSC_Embedding_Model.h5"

GDSC_CONTEXT_MODEL = load_and_distill_gdsc_features()

# -----------------------------------------------------------
# SECTION I: Data Engine Establishment and Curation Protocol
# -----------------------------------------------------------

# I.A. Data Acquisition and Curation (ChEMBL)
target_api = new_client.target
activity_api = new_client.activity

# Example Target: Human Cyclin-dependent kinase 2 (CDK2) 
CDK2_ID = 'P24941' 
targets = target_api.get(target_components__accession=CDK2_ID).only(['target_chembl_id', 'organism', 'pref_name', 'target_components'])

if targets:
    target_chembl_id = targets['target_chembl_id']
    target_name = targets['pref_name']
    print(f"Found target: {target_name} ({target_chembl_id})")
else:
    raise ValueError(f"Target with UniProt ID {CDK2_ID} not found.")

# Filter for high-quality, quantitative dose-response data [1, 4]
activities = activity_api.filter(target_chembl_id=target_chembl_id).filter(
    pchembl_value__isnull=False # Must have standardized pChEMBL value
).filter(
    assay_type__in= # Binding or Functional assays
).filter(
    standard_type__in=['IC50', 'EC50', 'Ki', 'Kd'] # Dose-response measurements
).only(['molecule_chembl_id', 'canonical_smiles', 'target_chembl_id', 'pchembl_value'])

df_raw = pd.DataFrame(list(activities))
print(f"--- Raw Filtered ChEMBL Data (N={len(df_raw)}) ---")

# --- Curation and Binary Label Generation ---
ACTIVITY_THRESHOLD = 6.0 # pChEMBL value >= 6.0 (1 uM affinity)
df_curated = df_raw[['canonical_smiles', 'target_chembl_id', 'pchembl_value']].copy()

# Deduplication and standardization (using median pChEMBL per SMILES) [4]
df_curated = df_curated.groupby('canonical_smiles').median().reset_index()

# Generate binary label (Interactor/Non-Interactor)
df_curated['label'] = (df_curated['pchembl_value'] >= ACTIVITY_THRESHOLD).astype(int)
df_curated['target_sequence'] = target_name 

drugs = df_curated['canonical_smiles'].values
targets = df_curated['target_sequence'].values
labels = df_curated['label'].values

print(f"\n--- Final Curated Data (N={len(df_curated)}) ---")
print(f"Positive Interactions (Active, label=1): {df_curated['label'].sum()}")


# -----------------------------------------------------------
# SECTION II: Enforcing Scaffold-Aware Generalization
# -----------------------------------------------------------

def generate_scaffold(smiles):
    """Generates the Bemis-Murcko Scaffold SMILES for a given molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 'Invalid_SMILES'
    # Uses the RDKit implementation for Murcko Scaffold [5, 6]
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) 

def scaffold_split_data(X_drugs, Y_labels, frac_test=0.2):
    """Splits data based on Bemis-Murcko scaffolds (out-of-distribution split)."""
    
    df = pd.DataFrame({'smiles': X_drugs, 'label': Y_labels})
    df['scaffold'] = df['smiles'].apply(generate_scaffold)
    
    # 1. Group data by unique scaffold
    scaffold_groups = df.groupby('scaffold')['smiles'].apply(list).to_dict()
    
    # 2. Sort groups by size (largest first) [7]
    sorted_scaffolds = sorted(scaffold_groups.items(), key=lambda x: len(x[1]), reverse=True) 
    
    test_size_needed = int(len(df) * frac_test)
    test_smiles =
    
    # 3. Assign full scaffolds to the test set until the required size is met [5, 8, 9]
    for _, smiles_list in sorted_scaffolds:
        if len(test_smiles) < test_size_needed:
            test_smiles.extend(smiles_list)
        else:
            break
            
    test_set_smiles = set(test_smiles)
    df_test = df[df['smiles'].isin(test_set_smiles)]
    df_train = df[~df['smiles'].isin(test_set_smiles)]
    
    # Validation check: essential for a rigorous OOD test [5]
    shared_scaffolds = set(df_train['scaffold'].unique()) & set(df_test['scaffold'].unique())
    print(f"\nScaffold Split Results:")
    print(f"Test Set Size: {len(df_test)} ({len(df_test)/len(df)*100:.2f}%)")
    print(f"Shared Scaffolds (must be 0): {len(shared_scaffolds)}")
    assert len(shared_scaffolds) == 0, "Scaffold splitting failed: Scaffolds are shared."
    
    return df_train['smiles'].values, df_train['target_sequence'].values, df_train['label'].values, \
           df_test['smiles'].values, df_test['target_sequence'].values, df_test['label'].values

# Execute the scaffold split
X_train_drug, X_train_target, Y_train, \
X_test_drug, X_test_target, Y_test = scaffold_split_data(drugs, targets, frac_test=0.2)


# -----------------------------------------------------------
# SECTION III: Model Construction and Training
# -----------------------------------------------------------

# Define the Hybrid Encoders [10, 11]
drug_encoding = 'AttentiveFP' 
target_encoding = 'Transformer' # DeepPurpose uses Transformer for contextual sequence embeddings

# 1. Prepare Data for DeepPurpose
print(f"\nPreparing data for Drug Encoder: {drug_encoding}, Target Encoder: {target_encoding}")

# Note: DeepPurpose will handle featurization (GNN for drug, Transformer for target)
train, val, _ = utils.data_process(X_drug=X_train_drug, X_target=X_train_target, Y=Y_train, 
                                   drug_encoding=drug_encoding, target_encoding=target_encoding, 
                                   split_method='random', frac=[0.8, 0.2, 0.0], random_seed=SEED)

test_set = utils.data_process(X_drug=X_test_drug, X_target=X_test_target, Y=Y_test, 
                              drug_encoding=drug_encoding, target_encoding=target_encoding, 
                              split_method='random', frac=[1.0, 0.0, 0.0], random_seed=SEED)
test_data = test_set


# 2. Model Initialization and Training (Tier II Hybrid)
# The multi-modal fusion layer (MLP with Attention in the proposal) is implicitly handled by 
# DeepPurpose's multi-modal classifier when both encoders are specified. [10, 12]
config = utils.generate_config(drug_encoding=drug_encoding, 
                               target_encoding=target_encoding, 
                               cls_hidden_dims = , 
                               train_epoch=10, 
                               LR=0.001, 
                               batch_size=128,
                               binary_thre=ACTIVITY_THRESHOLD) 

# Check if GPU is available and set device
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

model = DTI.model_initialize(**config).to(device)

print("\n--- Starting Training on Scaffold Split Data ---")
# Training monitors PR-AUC, preferred for imbalanced DTI data [13]
model.train(train, val, test_data, 
            training_loss = 'BCEWithLogitLoss', 
            model_path = 'multi_modal_dti_model.pth')


# -----------------------------------------------------------
# SECTION IV: Evaluation and Translational Metrics
# -----------------------------------------------------------

print("\n--- Evaluating on Rigorous Scaffold-Split Test Set ---")

result, _, pred_probas = model.test(test_data)

roc_auc = result['roc_auc']
pr_auc = result['precision_recall_auc']

print(f"\n*** Generalization Results (Scaffold-Split) ***")
print(f"ROC-AUC (Target >= 0.75): {roc_auc:.4f}")
print(f"PR-AUC (Preferred Metric for Imbalance): {pr_auc:.4f}")

# 2. Translational Metric: Enrichment Factor (EF@1%)
y_true = Y_test 
y_score = pred_probas

# Calculate EF at 1% cutoff using ODDT [14, 15]
EF_1_percent = enrichment_factor(y_true=y_true, 
                                y_score=y_score, 
                                percentage=1.0, 
                                pos_label=1, 
                                kind='fold')

print(f"\n*** Virtual Screening Performance ***")
print(f"Enrichment Factor @ 1% (Target >= 5): {EF_1_percent:.2f} Fold")

if roc_auc >= 0.75 and EF_1_percent >= 5.0:
    print("STATUS: PROJECT GOALS MET - Superior generalization and strong enrichment demonstrated.")
else:
    print("STATUS: OPTIMIZATION REQUIRED - Results did not meet all target thresholds.")

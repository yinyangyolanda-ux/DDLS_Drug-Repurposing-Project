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
print(f"Random seed set to {SEED} for full pipeline determinism.") [1]

# Placeholder for complex GDSC feature extraction/distillation
def load_and_distill_gdsc_features():
    """Simulate loading distilled GDSC genomic features (VAE/Autoencoder step)."""
    # This represents reducing >57k multi-omics features (GDSC/CCLE) to a dense vector
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

# Filter for high-quality, quantitative dose-response data
activities = activity_api.filter(target_chembl_id=target_chembl_id).filter(
    pchembl_value__isnull=False # Must have standardized pChEMBL value
).filter(
    # CORRECTION: Specify acceptable assay types (B: Binding, F: Functional)
    assay_type__in=
).filter(
    standard_type__in=['IC50', 'EC50', 'Ki', 'Kd'] # Only dose-response measurements 
).only(['molecule_chembl_id', 'canonical_smiles', 'target_chembl_id', 'pchembl_value'])

df_raw = pd.DataFrame(list(activities))
print(f"--- Raw Filtered ChEMBL Data (N={len(df_raw)}) ---")

# --- Curation and Binary Label Generation ---
ACTIVITY_THRESHOLD = 6.0 # pChEMBL value >= 6.0 (1 uM affinity)
df_curated = df_raw[['canonical_smiles', 'target_chembl_id', 'pchembl_value']].copy()

# Deduplication and standardization (using median pChEMBL per SMILES)
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
    # Uses RDKit's MurckoScaffoldSmiles function
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) 

def scaffold_split_data(X_drugs, Y_labels, frac_test=0.2):
    """Splits data based on Bemis-Murcko scaffolds (out-of-distribution split)."""
    
    df = pd.DataFrame({'smiles': X_drugs, 'label': Y_labels})
    df['scaffold'] = df['smiles'].apply(generate_scaffold)
    
    # 1. Group data by unique scaffold
    scaffold_groups = df.groupby('scaffold')['smiles'].apply(list).to_dict()
    
    # 2. Sort groups by size (largest first for deterministic, efficient splitting) [S23, S27

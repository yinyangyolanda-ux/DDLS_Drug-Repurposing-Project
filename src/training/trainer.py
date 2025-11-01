# 文件: src/training/trainer.py

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm
from sklearn.metrics import roc_auc_score, average_precision_score
import numpy as np
import os

# --- 占位符：自定义数据集 (需要封装所有四种特征和 KG ID) ---
class HybridDTIDataset(Dataset):
    def __init__(self, smiles, sequences, labels, drugbank_ids, uniprot_ids):
        # 简化：实际需要加载分子图和序列特征
        self.smiles = smiles
        self.sequences = sequences
        self.labels = labels
        self.drugbank_ids = drugbank_ids
        self.uniprot_ids = uniprot_ids
        # 还需要实现一个 collate_fn 来处理 DGLGraph 批处理

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        # 实际需要返回 DGLGraph 对象
        return (self.smiles[idx], self.sequences[idx], self.labels[idx], 
                self.drugbank_ids[idx], self.uniprot_ids[idx])

# --- DTI 训练器 ---
class DTITrainer:
    """
    负责 Hybrid TransformerDTA 模型的训练和零样本预测。
    """
    def __init__(self, model, kg_embeddings, config, train_loader, val_loader, test_loader):
        self.model = model
        self.kg_embeddings = kg_embeddings # 预训练的 CompGCN 关系嵌入
        self.config = config
        self.device = torch.device(config['project_info']['device'])
        
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.test_loader = test_loader
        self.criterion = nn.BCEWithLogitsLoss()
        
        self.optimizer = optim.Adam(model.parameters(), 
                                    lr=config['training']['learning_rate'],
                                    weight_decay=config['training']['weight_decay'])
        self.best_val_auc = 0.0
        self.checkpoint_dir = 'models/checkpoints'
        os.makedirs(self.checkpoint_dir, exist_ok=True)

    def train_epoch(self):
        """训练一轮。"""
        self.model.train()
        total_loss, all_preds, all_labels = 0.0,,
        
        pbar = tqdm(self.train_loader, desc='Training')
        for drug_batch, target_seqs, labels, drug_ids, target_ids in pbar:
            drug_batch = drug_batch.to(self.device)
            labels = labels.to(self.device)
            
            # 预测：模型接收结构/序列/关系嵌入所需的 ID
            logits = self.model(drug_batch, target_seqs, drug_ids, target_ids, self.kg_embeddings)
            loss = self.criterion(logits, labels)
            
            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()
            
            # 记录指标
            total_loss += loss.item()
            preds = torch.sigmoid(logits).detach().cpu().numpy()
            all_preds.extend(preds)
            all_labels.extend(labels.cpu().numpy())

        avg_loss = total_loss / len(self.train_loader)
        auc = roc_auc_score(all_labels, all_preds) if len(all_labels) > 1 and len(np.unique(all_labels)) > 1 else 0.0
        return avg_loss, auc

    def evaluate(self, loader):
        """评估模型性能。"""
        self.model.eval()
        all_preds, all_labels =,
        with torch.no_grad():
            for drug_batch, target_seqs, labels, drug_ids, target_ids in loader:
                drug_batch = drug_batch.to(self.device)
                logits = self.model(drug_batch, target_seqs, drug_ids, target_ids, self.kg_embeddings)
                
                preds = torch.sigmoid(logits).cpu().numpy()
                all_preds.extend(preds)
                all_labels.extend(labels.cpu().numpy())

        auc = roc_auc_score(all_labels, all_preds) if len(all_labels) > 1 and len(np.unique(all_labels)) > 1 else 0.0
        pr_auc = average_precision_score(all_labels, all_preds)
        return auc, pr_auc, np.array(all_labels), np.array(all_preds)

    def train(self):
        # 训练循环 (省略 early stopping 逻辑)
        print("\n--- 开始混合模型训练 ---")
        for epoch in range(1, self.config['training']['num_epochs'] + 1):
            train_loss, train_auc = self.train_epoch()
            val_auc, val_pr_auc, _, _ = self.evaluate(self.val_loader)
            
            print(f"Epoch {epoch} | Train AUC: {train_auc:.4f} | Val AUC: {val_auc:.4f} | Val PR-AUC: {val_pr_auc:.4f}")
            
            # 保存最佳模型
            if val_auc > self.best_val_auc:
                self.best_val_auc = val_auc
                torch.save(self.model.state_dict(), os.path.join(self.checkpoint_dir, 'best_hybrid_model.pt'))

    def predict_zero_shot(self, peptides_seqs, rare_target_id):
        """
        零样本预测：预测肽类药物对未见靶点（如 CFTR）的潜力。
        """
        # 1. 获取 CFTR 靶点的序列和关系嵌入 (Zero-Shot Target)
        # 需要确保 CFTR_HUMAN 在 KG 嵌入中存在
        if rare_target_id not in self.kg_embeddings:
            print(f"错误: 罕见病靶点 {rare_target_id} 不在知识图谱实体中。")
            return None
        
        # 2. 构造预测输入
        # 简化：药物 ID 使用 DrugBank ID，靶点 ID 使用 UniProt ID
        mock_drug_ids = [self.config['data']['peptide_targets'][p] for p in peptides_seqs.keys()]
        mock_target_ids = [rare_target_id] * len(peptides_seqs)
        target_seqs = list(peptides_seqs.values())
        
        # 3. 预测 (需要针对预测数据定制 DataLoader)
        # prediction_results = self.evaluate_zero_shot(mock_loader)
        
        print(f"已生成 {len(peptides_seqs)} 个肽类药物对 {rare_target_id} 的重定向预测。")
        return {"Target": rare_target_id, "Peptides": peptides_seqs.keys()}

if __name__ == '__main__':
    print("Trainer 模块已加载，请在 main.py 中运行。")

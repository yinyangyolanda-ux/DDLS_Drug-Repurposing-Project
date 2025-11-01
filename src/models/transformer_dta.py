# 文件: src/models/transformer_dta.py

import torch
import torch.nn as nn
import torch.nn.functional as F
from dgl.nn.pytorch import AttentiveFPGNN, AttentiveFPReadout
from dgllife.utils import CanonicalAtomFeaturizer, CanonicalBondFeaturizer, mol_to_bigraph
from transformers import AutoTokenizer, AutoModel
from rdkit import Chem
import dgl
import numpy as np

# --- 辅助函数：DGL 图批处理 ---
def collate_graphs(samples):
    graphs, target_sequences, labels, drug_ids = map(list, zip(*samples))
    batched_graph = dgl.batch(graphs)
    labels = torch.tensor(labels, dtype=torch.float32).unsqueeze(1)
    return batched_graph, target_sequences, labels, drug_ids

# --- 1. 结构编码器：AttentiveFP GNN [2, 3] ---
class AttentiveFPDrugEncoder(nn.Module):
    """拓扑感知药物编码器，输出固定维度嵌入。"""
    def __init__(self, node_feat_size=39, edge_feat_size=10, num_layers=2, 
                 num_timesteps=2, graph_feat_size=200, dropout=0.2):
        super().__init__()
        self.graph_feat_size = graph_feat_size
        self.gnn = AttentiveFPGNN(node_feat_size=node_feat_size, edge_feat_size=edge_feat_size,
                                  num_layers=num_layers, graph_feat_size=graph_feat_size, dropout=dropout)
        self.readout = AttentiveFPReadout(feat_size=graph_feat_size, num_timesteps=num_timesteps, dropout=dropout)

    def forward(self, g):
        node_feats = g.ndata['h'].float()
        edge_feats = g.edata['e'].float()
        
        atom_feats = self.gnn(g, node_feats, edge_feats)
        graph_emb = self.readout(g, atom_feats, edge_feats)
        return graph_emb

# --- 2. 序列编码器：ProtBERT/ESM Transformer [4, 3] ---
class ProteinTransformerEncoder(nn.Module):
    """使用 ProtBERT/ESM 预训练模型的蛋白质序列编码器，骨干层冻结。"""
    def __init__(self, model_name='Rostlab/prot_bert_bfd', embedding_dim=1024, freeze_backbone=True):
        super().__init__()
        self.model_name = model_name
        self.embedding_dim = embedding_dim
        
        # 加载分词器和模型 (ProtBERT 适用于序列编码)
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModel.from_pretrained(model_name)
        
        # 冻结骨干层 (符合项目要求) 
        if freeze_backbone:
            for param in self.model.parameters():
                param.requires_grad = False

        # 池化层，将序列输出降维到所需的嵌入维度
        self.pooling = nn.Linear(self.model.config.hidden_size, embedding_dim)
    
    def forward(self, sequences):
        # 准备序列 (ProtBERT 需要在氨基酸之间加入空格)
        sequences = [' '.join(list(seq)) for seq in sequences]
        
        inputs = self.tokenizer(sequences, return_tensors='pt', padding=True, truncation=True, max_length=1024)
        device = next(self.model.parameters()).device
        inputs = {k: v.to(device) for k, v in inputs.items()}
        
        with torch.no_grad(): # 在冻结层上禁用梯度
            outputs = self.model(**inputs)
            
        # 提取 CLS token (表示整个序列)
        cls_output = outputs.last_hidden_state[:, 0, :]
        
        # 降维到目标嵌入维度
        return self.pooling(cls_output)


# --- 3. 混合融合模块：Cross-Attention Hybrid Fusion [5] ---
class CrossAttentionFusion(nn.Module):
    """
    将四种嵌入（药物结构、药物关系、靶点序列、靶点关系）通过 Cross-Attention 融合。
    """
    def __init__(self, structural_dim, relational_dim, fusion_hidden_dim, num_heads):
        super().__init__()
        self.fusion_hidden_dim = fusion_hidden_dim
        
        # 线性投影，将所有输入维度统一到 fusion_hidden_dim
        self.proj_d_struct = nn.Linear(structural_dim, fusion_hidden_dim)
        self.proj_d_rel = nn.Linear(relational_dim, fusion_hidden_dim)
        self.proj_t_seq = nn.Linear(structural_dim, fusion_hidden_dim)
        self.proj_t_rel = nn.Linear(relational_dim, fusion_hidden_dim)

        # 多头交叉注意力 (用于融合结构和关系)
        self.cross_attention = nn.MultiheadAttention(fusion_hidden_dim, num_heads, batch_first=True)
        
        # 融合后的全连接层 (FFN)
        self.ffn = nn.Sequential(
            nn.Linear(fusion_hidden_dim, fusion_hidden_dim),
            nn.ReLU(),
            nn.Linear(fusion_hidden_dim, fusion_hidden_dim)
        )
        self.norm = nn.LayerNorm(fusion_hidden_dim)

    def forward(self, d_struct, d_rel, t_seq, t_rel):
        # 1. 投影到公共空间 (hidden_dim)
        d_h = self.proj_d_struct(d_struct) + self.proj_d_rel(d_rel)
        t_h = self.proj_t_seq(t_seq) + self.proj_t_rel(t_rel)

        d_h = d_h.unsqueeze(1) # (B, 1, H)
        t_h = t_h.unsqueeze(1) # (B, 1, H)

        # 2. Cross-Attention: Drug attends to Target (Q=D, K=V=T)
        d_attended, _ = self.cross_attention(d_h, t_h, t_h)
        
        # 3. Self-Attention: Target attends to Drug (Q=T, K=V=D)
        t_attended, _ = self.cross_attention(t_h, d_h, d_h)
        
        # 4. 残差连接与 FFN
        d_out = self.norm(d_h + d_attended)
        t_out = self.norm(t_h + t_attended)
        
        # 5. 拼接输出 (B, 2*H)
        fused = torch.cat([d_out.squeeze(1), t_out.squeeze(1)], dim=-1)
        return fused

# --- 4. 完整的 Hybrid DTA 模型 (Prediction Head) ---
class HybridTransformerDTA(nn.Module):
    """
    完整的混合 DTA 预测模型，整合 GNN, Transformer, CompGCN 嵌入。
    """
    def __init__(self, drug_encoder, target_encoder, kg_encoder, config):
        super().__init__()
        self.drug_encoder = drug_encoder
        self.target_encoder = target_encoder
        self.kg_encoder = kg_encoder # CompGCN 关系编码器
        
        # 确保结构和关系嵌入维度匹配
        structural_dim = config['model']['drug_encoder']['graph_feat_size']
        relational_dim = config['model']['kg_encoder']['embedding_dim']
        fusion_h_dim = config['model']['fusion']['hidden_dim']
        
        # 核心融合模块
        self.fusion = CrossAttentionFusion(
            structural_dim=structural_dim,
            relational_dim=relational_dim,
            fusion_hidden_dim=fusion_h_dim,
            num_heads=config['model']['fusion']['num_heads']
        )
        
        # 预测 MLP
        fusion_output_dim = fusion_h_dim * 2
        mlp_layers =
        in_dim = fusion_output_dim
        
        for hidden_dim in config['model']['decoder']['mlp_hidden_dims']:
            mlp_layers.extend(['decoder']['dropout'])
            ])
            in_dim = hidden_dim
            
        mlp_layers.append(nn.Linear(in_dim, 1)) # 输出为单个 DTI 分数
        self.mlp = nn.Sequential(*mlp_layers)
        
    def forward(self, drug_graphs, target_sequences, drug_ids, target_ids, kg_embeddings):
        # 1. 结构和序列编码
        d_struct = self.drug_encoder(drug_graphs)
        t_seq = self.target_encoder(target_sequences)
        
        # 2. 关系嵌入 (从预训练的 KG 嵌入中查询)
        # 假设 kg_embeddings 是一个包含所有实体嵌入的字典/映射
        d_rel = torch.tensor([kg_embeddings[d_id] for d_id in drug_ids], dtype=torch.float32).to(d_struct.device)
        t_rel = torch.tensor([kg_embeddings[t_id] for t_id in target_ids], dtype=torch.float32).to(d_struct.device)

        # 3. 混合融合
        fused = self.fusion(d_struct, d_rel, t_seq, t_rel)
        
        # 4. 最终预测
        return self.mlp(fused)

if __name__ == '__main__':
    # 简单的模型初始化演示
    class MockConfig:
        def __init__(self):
            self.model = {
                'drug_encoder': {'graph_feat_size': 200},
                'target_encoder': {'embedding_dim': 256},
                'kg_encoder': {'embedding_dim': 128},
                'fusion': {'hidden_dim': 256, 'num_heads': 4},
                'decoder': {'mlp_hidden_dims': , 'dropout': 0.2}
            }
    
    mock_config = MockConfig()
    
    # 初始化子模块
    drug_enc = AttentiveFPDrugEncoder(graph_feat_size=200)
    target_enc = ProteinTransformerEncoder(embedding_dim=256)
    
    # Mock KG 编码器 (用于初始化 Hybrid DTA 模型，实际嵌入来自预训练的 KG)
    class MockKGEncoder(nn.Module):
        def __init__(self, embed_dim):
            super().__init__()
            self.embedding_dim = embed_dim
            
    kg_enc = MockKGEncoder(mock_config.model['kg_encoder']['embedding_dim'])
            
    # 初始化 Hybrid DTA 模型
    hybrid_model = HybridTransformerDTA(drug_enc, target_enc, kg_enc, mock_config)
    print("Hybrid TransformerDTA Model Successfully Initialized.")
    print(f"总参数量: {sum(p.numel() for p in hybrid_model.parameters()):,}")

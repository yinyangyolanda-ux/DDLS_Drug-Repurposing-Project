DDLS_Drug-Repurposing-Project: 基于知识图谱和混合深度学习的零样本肽类药物重定向
项目概述
本项目旨在构建一个Hybrid Zero-Shot Repurposing Framework，通过整合多模态生物医学大数据（DrugBank XML, BindingDB, GDSC MoBEM），预测肽类药物（如 Tirzepatide, Semaglutide）对**罕见病靶点（如 CFTR）**的重定向潜力。

核心方法论是融合结构特征和关系特征，实现对训练集中未见的新靶点（Zero-Shot Targets）的精确泛化。

核心技术路线
多模态编码： 采用 AttentiveFP GNN 编码小分子拓扑结构；采用 ProtBERT/ESM Transformer 编码肽链和蛋白质序列的上下文信息。

知识图谱（DRKG）构建： 将药物、靶点、基因、疾病作为实体，从 DrugBank、BindingDB 和 DisGeNET 提取关系，构建异构图。

零样本学习： 使用 R-GCN/CompGCN 模型训练 DRKG，获取实体关系嵌入（Relational Embeddings）。

混合融合（Hybrid Fusion）： 在 Cross-Attention 机制中融合 结构/序列嵌入 和 关系嵌入，作为最终 DTI 预测模型的输入。

肽类药物重定向用例
本项目将使用以下 GLP-1 类似物作为案例，预测它们在胰岛功能障碍和炎症通路中的潜在重定向靶点，例如囊性纤维化相关的靶点：

Tirzepatide (DB15171): GIP/GLP-1 双重激动剂

Semaglutide (DB13928): GLP-1 激动剂

Pegloxenatide (DB19230)

2. 项目目录结构
DDLS_Drug-Repurposing-Project/ ├── src/ │ ├── data/ │ │ ├── drugbank_loader.py # 【实现】解析 DrugBank XML，提取序列、SMILES、DGI │ │ ├── bindingdb_loader.py # 【实现】解析 BindingDB TSV/Dump，提取定量亲和力数据 │ │ ├── gdsc_loader.py # 【实现】GDSC MoBEM 数据加载及基因组特征降维（VAE/Autoencoder） │ │ ├── uniprot_loader.py # 【实现】通过 UniProt API 批量获取蛋白质序列 │ │ ├── disgenet_loader.py # 【新增】获取基因-疾病关联数据 (G-D)，强化零样本路径 │ │ ├── kg_builder.py # 【新增】协调数据源，构建 DGL/PyG 异构知识图谱 (DRKG) │ │ └── featurize.py # 【新增】分子图/肽链序列的特征化工具 │ ├── models/ │ │ ├── gnn_models.py # 【实现】AttentiveFP GNN 编码器，R-GCN/CompGCN 关系编码器 │ │ ├── sequence_encoders.py # 【实现】ProtBERT/ESM Transformer 序列编码器 │ │ ├── fusion_models.py # 【实现】Cross-Attention 混合融合模块 │ │ └── transformer_dta.py # 【实现】总混合模型架构 (DTA Prediction Head) │ ├── training/ │ │ ├── data_splits.py # 【实现】支架感知分割 (Scaffold Split) │ │ ├── trainer.py # 【实现】训练循环及零样本预测逻辑 │ │ └── hyperparameter_tuning.py │ └── evaluation/ │ ├── metrics.py # 【实现】ROC-AUC, EF@1%, PR-AUC 计算 │ ├── visualization.py # 【实现】分子结构 Grad-CAM 热力图 │ └── repurposing_predictor.py # 【新增】零样本预测与结果解释 (SHAP) ├── config/ │ ├── base_config.yaml # 【实现】项目配置，数据路径，超参数 │ └── model_configs.yaml ├── data/ │ ├── raw/ │ ├── processed/ │ └── knowledge_graph/ # 存放 DRKG 文件和实体嵌入 ├── notebooks/ └── main.py # 【实现】主执行脚本

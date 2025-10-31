# 文件: src/data/bindingdb_loader.py

import pandas as pd
import os
import zipfile
from io import StringIO

BINDINGDB_PATH = '/content/gdrive/MyDrive/Colab_Projects/BDB-mySQL_All_202511_dmp.zip'
TSV_FILE_NAME = 'BindingDB_All.tsv' # 假设核心数据文件名为TSV

class BindingDBLoader:
    """
    用于加载 BindingDB TSV 数据，专注于定量亲和力数据。
    """
    def __init__(self, zip_path, cache_dir='data/processed'):
        self.zip_path = zip_path
        self.cache_dir = cache_dir

    def load_and_preprocess(self):
        """从TSV加载并筛选定量结合数据（Ki, IC50, Kd, EC50）。"""
        print(f"正在读取 BindingDB 文件: {self.zip_path}")
        
        # 假设 TSV 在 ZIP 文件内部
        with zipfile.ZipFile(self.zip_path, 'r') as zf:
            # 找到包含 BindingDB 数据的 TSV 文件
            tsv_data_name = next((name for name in zf.namelist() if name.endswith('.tsv') or name.endswith('.txt')), TSV_FILE_NAME)
            
            with zf.open(tsv_data_name) as f:
                # 使用 StringIO 处理读取的数据，并指定分隔符为 Tab 
                # 注意：需要跳过可能的注释行（# 或其他元数据行）
                data = f.read().decode('utf-8')
                df = pd.read_csv(StringIO(data), sep='\t', skipinitialspace=True, low_memory=False)

        print(f"原始 BindingDB 记录数: {len(df)}")

        # 筛选关键列 
        key_cols =
        
        df_clean = df[key_cols].copy()
        
        # 筛选至少具有一种定量亲和力测量的记录
        affinity_cols = ['Ki (nM)', 'IC50 (nM)', 'Kd (nM)', 'EC50 (nM)']
        df_clean = df_clean.dropna(subset=affinity_cols, how='all')
        
        # 清理数据：移除非数值符号（如 '>'）并转换为 pAffinity
        for col in affinity_cols:
            df_clean[col] = df_clean[col].astype(str).str.replace('[^0-9\.]', '', regex=True)
            df_clean[col] = pd.to_numeric(df_clean[col], errors='coerce')
        
        # 进一步的 pAffinity 标准化和合并将是下一步骤

        df_clean.rename(columns={'Ligand SMILES': 'smiles', 'Target Name': 'target_name'}, inplace=True)
        df_clean.to_csv(os.path.join(self.cache_dir, 'bindingdb_processed_dti.csv'), index=False)
        print(f"BindingDB 数据处理完成，有效DTI记录数: {len(df_clean)}")
        return df_clean

if __name__ == '__main__':
    # 示例运行
    loader = BindingDBLoader(zip_path=BINDINGDB_PATH)
    loader.load_and_preprocess()

# 文件: src/data/bindingdb_loader.py

import pandas as pd
import os
import zipfile
from io import StringIO
from tqdm import tqdm

BINDINGDB_FILE_NAME = 'BindingDB_All.tsv' # 假设核心数据文件名为 TSV 或类似格式

class BindingDBLoader:
    """
    加载和预处理 BindingDB 数据，提取定量 DTI 记录。
    数据源: BDB-mySQL_All_202511_dmp.zip
    """
    def __init__(self, zip_path, cache_dir='data/processed'):
        self.zip_path = zip_path
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)
        # 确定需要提取的定量亲和力列 
        self.affinity_cols = ['Ki (nM)', 'IC50 (nM)', 'Kd (nM)', 'EC50 (nM)']
        self.key_cols = + self.affinity_cols

    def load_and_preprocess(self):
        """从 ZIP 文件中加载 TSV 数据，并进行初步清洗和筛选。"""
        print(f"正在读取 BindingDB 文件: {self.zip_path}")
        
        with zipfile.ZipFile(self.zip_path, 'r') as zf:
            # 尝试找到核心 TSV 文件
            tsv_data_name = next((name for name in zf.namelist() if name.endswith('.tsv') or name.endswith('.txt')), None)
            
            if tsv_data_name is None:
                raise FileNotFoundError(f"在 {self.zip_path} 中未找到 BindingDB 核心 TSV/TXT 数据文件。")
            
            with zf.open(tsv_data_name) as f:
                # 使用 StringIO 处理读取的数据，分隔符为 Tab 
                data = f.read().decode('utf-8', errors='ignore')
                df = pd.read_csv(StringIO(data), sep='\t', skipinitialspace=True, low_memory=False)

        print(f"原始 BindingDB 记录数: {len(df)}")
        
        # 1. 筛选关键列 (确保列名存在，BindingDB 版本可能略有不同)
        df_clean = df.rename(columns={'Ligand SMILES': 'smiles', 
                                      'UniProt (SwissProt) Primary ID of Target Chain': 'uniprot_id'}).copy()
        
        # 2. 筛选至少具有一种定量亲和力测量的记录 
        df_clean = df_clean.dropna(subset=self.affinity_cols, how='all')
        
        # 3. 清理亲和力值：移除非数值符号 (如 ">", "<") 并转换为数值
        for col in tqdm(self.affinity_cols, desc="Cleaning Affinity Values"):
            if col in df_clean.columns:
                # 移除所有非数字和小数点的字符
                df_clean[col] = df_clean[col].astype(str).str.replace('[^0-9\.]', '', regex=True)
                df_clean[col] = pd.to_numeric(df_clean[col], errors='coerce')
        
        # 4. 筛选人类靶点，并确保具有 SMILES 和 UniProt ID
        df_clean = df_clean.str.contains('Human', na=False)]
        df_clean = df_clean.dropna(subset=['smiles', 'uniprot_id'])

        df_clean.to_csv(os.path.join(self.cache_dir, 'bindingdb_processed_dti.csv'), index=False)
        print(f"BindingDB 数据处理完成，有效DTI记录数: {len(df_clean)}")
        return df_clean

if __name__ == '__main__':
    # 示例路径，在 Colab 中执行时需确认 Drive 挂载
    TEST_PATH = '../../Colab_Projects/BDB-mySQL_All_202511_dmp.zip'
    if os.path.exists(TEST_PATH):
        loader = BindingDBLoader(zip_path=TEST_PATH)
        loader.load_and_preprocess()
    else:
        print(f"警告: 绑定数据库测试文件 {TEST_PATH} 不存在。")

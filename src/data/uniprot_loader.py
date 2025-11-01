# 文件: src/data/uniprot_loader.py

import requests
from Bio import SeqIO
from io import StringIO
import os
import time
from tqdm import tqdm

class UniProtFetcher:
    """
    通过 UniProt REST API 批量获取蛋白质序列，并支持缓存。
    """
    def __init__(self, cache_dir='data/raw/uniprot'):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.base_url = 'https://rest.uniprot.org/uniprotkb/'

    def fetch_sequence(self, uniprot_id, retry=3):
        """获取单个 UniProt ID 的序列。"""
        cache_file = os.path.join(self.cache_dir, f'{uniprot_id}.fasta')
        
        # 1. 检查缓存
        if os.path.exists(cache_file):
            with open(cache_file, 'r') as f:
                try:
                    record = next(SeqIO.parse(f, 'fasta'))
                    return str(record.seq)
                except StopIteration:
                    return None
        
        # 2. 从 UniProt API 获取
        url = f"{self.base_url}{uniprot_id}.fasta"
        for attempt in range(retry):
            try:
                response = requests.get(url, timeout=10)
                response.raise_for_status() # 检查 HTTP 错误
                
                # 解析 FASTA 响应
                record = next(SeqIO.parse(StringIO(response.text), 'fasta'))
                sequence = str(record.seq)
                
                # 保存到缓存
                with open(cache_file, 'w') as f:
                    f.write(response.text)
                return sequence

            except requests.exceptions.RequestException as e:
                print(f"[{uniprot_id}] 第 {attempt+1} 次请求失败: {e}")
                if attempt < retry - 1:
                    time.sleep(2 ** attempt) # 指数退避
                else:
                    return None
        return None

    def fetch_batch(self, uniprot_ids):
        """批量获取蛋白质序列。"""
        sequences = {}
        print(f"开始批量获取 {len(uniprot_ids)} 个序列...")
        for uid in tqdm(uniprot_ids, desc="Fetching UniProt Sequences"):
            seq = self.fetch_sequence(uid)
            if seq:
                sequences[uid] = seq
        print(f"成功获取 {len(sequences)} 个序列。")
        return sequences

if __name__ == '__main__':
    fetcher = UniProtFetcher()
    # 案例靶点：CFTR, GLP1R (Semaglutide/Tirzepatide靶点), GIPR (Tirzepatide靶点)
    test_ids = ['P13569', 'P43220', 'P48039'] 
    
    sequences = fetcher.fetch_batch(test_ids)
    
    # 验证 CFTR 序列
    if 'P13569' in sequences:
        print(f"CFTR (P13569) 序列长度: {len(sequences['P13569'])}")

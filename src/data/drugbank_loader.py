# 文件: src/data/drugbank_loader.py

import zipfile
import xml.etree.ElementTree as ET
import pandas as pd
import os
from tqdm import tqdm

# 注意：需要将 DrugBank XML 文件上传到您的 Google Drive
DB_ZIP_PATH = '/content/gdrive/MyDrive/Colab_Projects/drugbank_all_full_database.xml.zip'
XML_FILE_NAME = 'full_database.xml' # 假设XML文件名

class DrugBankXMLLoader:
    """
    用于解析 DrugBank 完整 XML 数据库的类，提取药物、靶点、肽链序列和DGI数据。
    """
    def __init__(self, zip_path, cache_dir='data/processed'):
        self.zip_path = zip_path
        self.cache_dir = cache_dir
        self.namespace = '{http://www.drugbank.ca}' # DrugBank XML 命名空间

    def extract_xml(self):
        """解压 XML 文件并返回根元素。"""
        print(f"正在读取 DrugBank 文件: {self.zip_path}")
        with zipfile.ZipFile(self.zip_path, 'r') as zf:
            with zf.open(XML_FILE_NAME) as f:
                tree = ET.parse(f)
                return tree.getroot()

    def parse_drugs(self):
        """解析 XML 并提取核心药物数据，包括肽链序列和靶点信息。"""
        root = self.extract_xml()
        data =

        for drug in tqdm(root, desc="Parsing DrugBank XML"):
            db_id = drug.find(f'{self.namespace}drugbank-id').text
            name = drug.find(f'{self.namespace}name').text
            mol_type = drug.find(f'{self.namespace}kind').text # e.g., Small molecule drug, Polypeptide
            
            # 提取 SMILES
            smiles = ''
            calculated_props = drug.find(f'{self.namespace}calculated-properties')
            if calculated_props is not None:
                for prop in calculated_props:
                    if prop.find(f'{self.namespace}kind').text == 'SMILES':
                        smiles = prop.find(f'{self.namespace}value').text
                        break

            # 提取肽链/蛋白序列
            sequence = ''
            polypeptide = drug.find(f'{self.namespace}polypeptide')
            if polypeptide is not None:
                sequence = polypeptide.find(f'{self.namespace}sequence').text
                
            # 提取靶点（D-T / D-G 关系）
            targets =
            targets_node = drug.find(f'{self.namespace}targets')
            if targets_node is not None:
                for target in targets_node:
                    gene_name = target.find(f'{self.namespace}name').text
                    uniprot_id = target.find(f'{self.namespace}polypeptide/{self.namespace}external-identifiers')
                    uniprot_acc = ''
                    if uniprot_id is not None:
                        for ext_id in uniprot_id:
                            if ext_id.find(f'{self.namespace}resource').text == 'UniProtKB':
                                uniprot_acc = ext_id.find(f'{self.namespace}identifier').text
                                break
                    targets.append({'gene': gene_name, 'uniprot': uniprot_acc})

            data.append({
                'drugbank_id': db_id,
                'name': name,
                'kind': mol_type,
                'smiles': smiles,
                'sequence': sequence,
                'targets': targets # D-T/D-G关系，用于KG构建
            })

        df = pd.DataFrame(data)
        df.to_csv(os.path.join(self.cache_dir, 'drugbank_processed_base.csv'), index=False)
        print("\nDrugBank XML解析完成，基础数据已保存。")
        return df

if __name__ == '__main__':
    # 示例运行：在Colab中运行此脚本前，请确保挂载了Google Drive
    loader = DrugBankXMLLoader(zip_path=DB_ZIP_PATH)
    loader.parse_drugs()

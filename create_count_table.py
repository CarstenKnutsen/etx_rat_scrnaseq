import os
import scanpy as sc
import numpy as np
import pandas as pd
import string
import anndata
from gtfparse import read_gtf
from anndata import AnnData
gtf_fn = '/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/Rattus_norvegicus.mRatBN7.2.105.filtered.gtf'
data = '/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/soupx'
output = '/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/processed/single_cell_files'
os.makedirs(output, exist_ok=True)
def read_adata(folder):
    adata = sc.read_mtx(f'{folder}/matrix.mtx').T
    features = pd.read_csv(f'{folder}/genes.tsv',
                          sep = '\t',
                          header = None)
    bc = pd.read_csv(f'{folder}/barcodes.tsv',
                          sep = '\t',
                    header = None)
    features.rename(columns={0:'gene_id',
                            1: 'gene_symbol',
                            2: 'category'},
                   inplace = True)

    adata.var = features
    adata.obs_names = bc[0]
    adata.var_names = adata.var['gene_id'].values
    return adata
if __name__ == '__main__':
    runs = os.listdir(data)
    adatas = []
    gtf = read_gtf(gtf_fn)
    gtf['gene_id'] = gtf['gene_id'].str.split('.').str[0]
    gene_name_dict = pd.Series(gtf.gene_name.values, index=gtf.gene_id).to_dict()
    for x in gene_name_dict.keys():
        if gene_name_dict[x] == '':
            gene_name_dict[x] = x
    for run in runs:
        print(run)
        folder = f'{data}/{run}/outs/soupx_filt'
        adata = read_adata(folder)
        adata.var_names = [gene_name_dict[x] for x in adata.var_names]
        df = pd.DataFrame(adata.X.todense(),
                          columns=adata.var_names.tolist(),
                          index=adata.obs_names.tolist())
        df = df.sum(axis=1, level=0)
        adata2 = AnnData(df)
        adata2.obs = adata.obs
        adata = adata2.copy()
        adata.obs_names = run + '_' + adata.obs_names
        adata.obs['Treatment'] = run.split('_')[0]
        adata.obs['Timepoint'] = run.split('_')[1]
        adata.obs['Individual'] = '_'.join(run.split('_')[0:3])
        adata.obs['Run'] = run
        adatas.append(adata.copy())
    adata = anndata.concat(adatas)
    for column in ['gene_id', 'gene_biotype', 'seqname', 'transcript_name', 'protein_id']:
        temp_dict = pd.Series(gtf[column].values, index=gtf['gene_name']).to_dict()
        temp_dict.update(pd.Series(gtf[column].values, index=gtf['gene_id']).to_dict())
        adata.var[column] = [temp_dict[x] for x in adata.var.index]
    adata.var['mt'] = [True if x == 'MT' else False for x in adata.var['seqname']]
    adata.write(f'{output}/rat_etx_all_cells_raw.gz.h5ad', compression='gzip')
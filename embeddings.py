import pandas as pd
import os
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import leidenalg as la
import northstar
import itertools
import anndata

figures = "/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/figures/multiple_embedding"
data = "/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/processed/single_cell_files"
os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi=300, dpi_save=300, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    # read in processed data
    adata = sc.read(f"{data}/etx_rat_processed_cell_typed_raw.gz.h5ad")
    adata.obs['Timepoint'] = pd.Categorical(adata.obs['Timepoint'], categories= ['E20','D0','D14'], ordered=True)
    ct_order = []
    ct_number = 0
    ct_number_dict = {}
    for lin in adata.obs['lineage'].cat.categories:
        for ct in adata[adata.obs['lineage'] == lin].obs['Cell Subtype'].cat.categories:
            ct_order.append(ct)
            ct_number += 1
            ct_number_dict[ct] = str(ct_number)
    adata.obs['ct_number_name'] = [f'{ct_number_dict[x]}. {x} ' for x in adata.obs['Cell Subtype']]
    adata.obs['ct_number'] = [f'{ct_number_dict[x]}' for x in adata.obs['Cell Subtype']]
    adata.obs['Cell Subtype'] = pd.Categorical(adata.obs['Cell Subtype'], categories=ct_order)


    sc.pp.normalize_total(adata, target_sum = 1e6)
    sc.pp.log1p(adata, base =10)
    print(adata)
    # marker gene dicts for each lineage
    gene_dict = {
        "mesenchymal": [
            "Col3a1",
            "G0s2",
            "Limch1",
            "Col13a1",
            "Col14a1",
            "Serpinf1",
            "Pdgfra",
            "Scara5",
            "Acta2",
            "Hhip",
            "Fgf18",
            "Wif1",
            "Tgfbi",
            "Tagln",
            "Mustn1",
            "Aard",
            "Pdgfrb",
            "Cox4i2",
            # 'Higd1b',
            "Myh11",
            'Msln',
            "Wt1",
            "Lrrn4",
            "Upk3b",
            "Mki67",
            "Acta1",
            # 'Crh',
            'Lgr6',
            'Cd34',
            'Tnc',
            'Rbp1',
            'Macf1',
            'Aspn',
            'Cthrc1',
            "Epcam",
            "Ptprc",
            "Pecam1",
        ],
        "endothelial": [
            "Gja5",
            "Bmx",
            "Fn1",
            "Ctsh",
            # 'Kcne3',
            "Cdh13",
            # 'Car8',
            "Mmp16",
            "Slc6a2",
            "Thy1",
            "Mmrn1",
            # 'Ccl21a',
            "Reln",
            "Neil3",
            "Mki67",
            "Aurkb",
            "Depp1",
            "Ankrd37",
            "Peg3",
            "Mest",
            "Hpgd",
            # 'Cd36',
            'Vwf',
            'Ca4',
            "Sirpa",
            "Fibin",
            "Col1a1",
            "Epcam",
            "Ptprc",
            "Pecam1",
        ],
        "immune": [
            "Cd68",
            "Gal",
            "Itgax",
            'Ca4',
            "C1qa",
            'Ccl3',
            'Nkg7',
            "Plac8",
            "Batf3",
            'Rora',
            "Itgae",
            # 'Cd209a',
            "Mreg",
            # 'Mcpt8',
            "Retnlg",
            "Ms4a1",
            "Gzma",
            "Cd3e",
            "Areg",
            "Mki67",
            "Col1a1",
            "Epcam",
            "Ptprc",
            "Pecam1",
        ],
        "epithelial": [
            "Scg5",
            "Ascl1",
            # 'Lyz1',
            "Lyz2",
            "Sftpc",
            "Slc34a2",
            "S100g",
            "Sftpa1",
            "Akap5",
            "Hopx",
            # 'Col4a4',
            "Vegfa",
            "Tmem212",
            "Dynlrb2",
            "Cdkn1c",
            "Tppp3",
            'Krt17',
            'Sfn',
            'Krt5',
            # 'Scgb3a2',
            # 'Cyp2f2',
            "Scgb1a1",
            # 'Reg3g',
            "Scgb3a1",
            "Mki67",
            "Col1a1",
            "Epcam",
            "Ptprc",
            "Pecam1",
        ],
    }
    metadata = [
                 'leiden',
                 'lineage',
                 'n_genes_by_counts',
                 'Run',
                 'Treatment',
                 'log1p_total_counts',
                 'Timepoint',
        'Cell Subtype'
                 ]
    b_dict = {"mesenchymal": 1, "endothelial": 0.7, "epithelial": 1, "immune": 1.2}
    figures_all = f"{figures}/all"
    figures_all_genes = f"{figures_all}/all/genes"
    figures_all_meta = f"{figures_all}/all/meta"
    gene_ls = []
    for k in gene_dict.keys():
        gene_ls = gene_ls + gene_dict[k]
    sc.settings.figdir = figures_all_genes
    for gene in gene_ls:
        sc.pl.umap(adata, color = gene, alpha=0.5,show=False, save = f'_{gene}.png')
    sc.settings.figdir = figures_all_meta
    for meta in metadata:
        sc.pl.umap(adata, color = meta, alpha=0.5,show=False, save = f'_{meta}.png')
    for tp in adata.obs['Timepoint'].cat.categories:
        figures_tp_gene = f"{figures_all}/{tp}/genes"
        figures_tp_meta = f"{figures_all}/{tp}/meta"
        tp_adata = adata[adata.obs['Timepoint'] == tp]
        sc.pp.highly_variable_genes(tp_adata, n_top_genes=500,batch_key='Run')
        sc.pp.pca(tp_adata)
        sc.pp.neighbors(tp_adata)
        sc.tl.umap(tp_adata)
        sc.settings.figdir = figures_tp_gene
        for gene in gene_ls:
            sc.pl.umap(tp_adata, color=gene, alpha=0.5, show=False,save=f'_{tp}_{gene}.png')
        sc.settings.figdir = figures_tp_meta
        for meta in metadata:
            if meta == 'leiden':
                sc.pl.umap(tp_adata, color=meta, alpha=0.8,legend_loc='on data', show=False, save=f'_{tp}_{meta}.png')
            else:
                sc.pl.umap(tp_adata, color=meta, alpha=0.8, show=False,save=f'_{tp}_{meta}.png')
    for lineage in adata.obs['lineage'].cat.categories:
        metadata_lin = metadata +[f'leiden_{lineage}']
        figures_lineage_gene = f"{figures}/{lineage}/all/genes"
        figures_lineage_meta = f"{figures}/{lineage}/all/meta"
        lineage_adata = adata[adata.obs['lineage'] == lineage]
        sc.settings.figdir = figures_lineage_gene
        for gene in gene_dict[lineage]:
            sc.pl.embedding(lineage_adata,basis=f"X_umap_{lineage}", color=gene, alpha=0.5, show=False,save=f'_{gene}.png')
        sc.settings.figdir = figures_lineage_meta
        for meta in metadata_lin:
            if meta in ['leiden',f'leiden_{lineage}']:
                sc.pl.embedding(lineage_adata,basis= f"X_umap_{lineage}", color=meta, alpha=0.8,legend_loc='on data', show=False, save=f'_{meta}.png')
            else:
                sc.pl.embedding(lineage_adata,basis= f"X_umap_{lineage}", color=meta, alpha=0.8, show=False,save=f'_{meta}.png')
        for tp in adata.obs['Timepoint'].cat.categories:
            figures_tp_genes = f"{figures}/{lineage}/{tp}/genes"
            figures_tp_meta = f"{figures}/{lineage}/{tp}/meta"
            tp_adata = lineage_adata[lineage_adata.obs['Timepoint'] == tp]
            sc.pp.highly_variable_genes(tp_adata, n_top_genes=500, batch_key='Run')
            sc.pp.pca(tp_adata)
            sc.pp.neighbors(tp_adata)
            sc.tl.umap(tp_adata, a=2, b=b_dict[lineage])
            sc.settings.figdir = figures_tp_genes
            for gene in gene_dict[lineage]:
                sc.pl.umap(tp_adata, color=gene, alpha=0.5, show=False,save=f'_{lineage}_{tp}_{gene}.png')
            sc.settings.figdir = figures_tp_meta
            for meta in metadata_lin:
                if meta in ['leiden',f'leiden_{lineage}']:
                    sc.pl.umap(tp_adata, color=meta, alpha=0.8, legend_loc='on data', show=False, save=f'_{lineage}_{tp}_{meta}.png')
                else:
                    sc.pl.umap(tp_adata, color=meta, alpha=0.8, show=False, save=f'_{lineage}_{tp}_{meta}.png')
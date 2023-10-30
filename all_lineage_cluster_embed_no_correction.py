import pandas as pd
import os
import scanpy as sc
import leidenalg as la
import itertools
import scanpy.external as sce

figures = '/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/figures/tissue_embedding_base'
data = '/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/processed/single_cell_files'
os.makedirs(figures, exist_ok = True)
sc.set_figure_params(dpi = 300, dpi_save = 300, format = 'png')
sc.settings.figdir = figures
if __name__ == '__main__':
    adata = sc.read(f'{data}/rat_etx_qc_trimmed_raw.gz.h5ad')
    print(adata)
    adata.layers['raw'] = adata.X.copy()
    sc.pp.normalize_total(adata, key_added=None, target_sum=1e4)
    adata.layers['cp10k'] = adata.X.copy()
    sc.pp.log1p(adata)
    adata.layers['log10'] = adata.X.copy()
    sc.pp.highly_variable_genes(adata,
                                n_top_genes = 2000,
                                batch_key='Run'
                                )
    sc.pp.pca(adata, use_highly_variable=True)
    sc.pp.neighbors(adata)
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
    sc.tl.umap(adata, min_dist=0.5)
    sc.tl.tsne(adata)
    sc.tl.leiden(adata, resolution=0.0005, partition_type=la.CPMVertexPartition, key_added='leiden')
    print(len(adata.obs['leiden'].unique()))
    # sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    # sc.pl.rank_genes_groups_dotplot(adata, groupby='leiden',
    #                                 n_genes=int(150 / len(adata.obs['leiden'].unique())), show=False,
    #                                 save=f'leiden_markers.png')
    adata = adata[~adata.obs['leiden'].isin(['39','72','77']),:] # Cluster 40/43 look low quality
    # sc.pl.pca_overview(adata, color='leiden', show=False, save = True)
    # sc.pl.pca_variance_ratio(adata, show = False, save ='variance_ratio')
    # sc.pl.pca_loadings(adata, components = ','.join([str(x) for x in range(1,10)]), show = False, save = True)
    ## Assign lineage from gene expression
    genes = ['Col1a1', 'Cdh5', 'Ptprc', 'Epcam', 'leiden']
    genedf = sc.get.obs_df(adata, keys=genes)
    grouped = genedf.groupby("leiden")
    mean = grouped.mean()
    mean_t = mean.T
    mean_t.to_csv(f'{figures}/leiden_lineage_scores.csv')
    lineage_dict = {}
    for cluster in mean_t.columns:
        gene = mean_t[cluster].idxmax()
        if gene == 'Cdh5':
            lineage_dict[cluster] = 'endothelial'
        elif gene == 'Ptprc':
            lineage_dict[cluster] = 'immune'
        elif gene == 'Epcam':
            lineage_dict[cluster] = 'epithelial'
        elif gene == 'Col1a1':
            lineage_dict[cluster] = 'mesenchymal'
    adata.obs['lineage'] = [lineage_dict[x] for x in adata.obs['leiden']]
    adata.write(f'{data}/rat_etx_processed_log_no_batch.gz.h5ad', compression='gzip')
    print(adata)
    for gene in ['Col1a1',
                 'Epcam',
                 'Ptprc',
                 'Pecam1',
                 'Cdh5',
                 'Acta2',
                 'Sirpa',
         'Ca4',
        'Fibin',
        'Mmrn1',
                 # 'Higd1b',
        'Peg3',
                 'Sftpc',
                 'Ager',
                 'Acta1',
                 'Tubb3',
                 'Mki67',
                 'leiden',
                 'lineage',
                 'n_genes_by_counts',
                 'Run',
                 'Treatment',
                 'total_counts',
                 'Timepoint',
                 'doublet',
                 'doublet_score',
                 'pct_counts_mt'
                 ]:
        if gene == 'leiden':
            sc.pl.umap(adata, legend_loc='on data', color=gene, alpha=0.5, show = False,save=f'_{gene}')
        else:
            sc.pl.umap(adata, color=gene, alpha=0.5, show = False, save=f'_{gene}')
    pd.DataFrame(index=adata.obs_names, data = adata.obsm['X_umap']).to_csv(f'{figures}/allcells_umap.csv')
    with pd.ExcelWriter(f'{figures}/metadata_counts.xlsx', engine='xlsxwriter') as writer:
        obs_list = ['lineage', 'Treatment','Timepoint']
        num_obs = len(obs_list) + 1
        for ind in range(0, num_obs):
            for subset in itertools.combinations(obs_list, ind):
                if len(subset) != 0:
                    subset = list(subset)
                    if len(subset) == 1:
                        key = subset[0]
                        adata.obs[key].value_counts().to_excel(writer, sheet_name=key)
                    else:
                        key = "_".join(subset)
                        adata.obs.groupby(subset[:-1])[subset[-1]].value_counts().to_excel(writer, sheet_name=key[:31])



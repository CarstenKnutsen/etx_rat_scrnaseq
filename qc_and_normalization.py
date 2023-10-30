import pandas as pd
import os
import scanpy as sc
from anndata import AnnData
import seaborn as sns
import matplotlib.pylab as plt
import doubletdetection

figures = '/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/figures/qc'
data = '/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/processed/single_cell_files'
gtf_fn = '/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/Rattus_norvegicus.mRatBN7.2.105.filtered.gtf'
os.makedirs(figures, exist_ok = True)
os.makedirs(data, exist_ok = True)
sc.set_figure_params(dpi = 300, format = 'png')
sc.settings.figdir = figures

if __name__ == '__main__':
    ##Convert gene_id to gene_name and add metadata from gtf file
    adata = sc.read(f'{data}/rat_etx_all_cells_raw.gz.h5ad')
    print(adata)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
    sc.pl.scatter(adata, x='log1p_total_counts', y='n_genes_by_counts', show=False, save='genes_by_counts_pretrim_log')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=False, save='genes_by_counts_pretrim_log')
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, save=f'_pretrim')
    sc.pl.violin(adata, ['log1p_total_counts'], save = 'counts_pretrim')
    sc.pl.violin(adata, ['n_genes_by_counts'], save = 'genes_pretrim')
    sc.pp.filter_cells(adata, min_counts=500)
    sc.pp.filter_cells(adata, min_genes=200)
    print(adata)
    sc.pp.filter_genes(adata, min_cells=10)
    adata = adata[adata.obs.pct_counts_mt < 20, :]
    print(adata)
    sc.pl.scatter(adata, x='log1p_total_counts', y='n_genes_by_counts', show = False, save = 'genes_by_counts_posttrim_log')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show = False, save = 'genes_by_counts_posttrim')
    sc.pl.highest_expr_genes(adata, n_top=20, show = False, save = f'_posttrim')
    sc.pl.violin(adata, ['log1p_total_counts'], save='counts_posttrim')
    sc.pl.violin(adata, ['n_genes_by_counts'], save='genes_posttrim')
    # doubletdetection
    clf = doubletdetection.BoostClassifier(
        n_iters=20,
        clustering_algorithm="leiden",
        standard_scaling=True,
        pseudocount=0.1,
        n_jobs=-1,
    )
    doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
    doublet_score = clf.doublet_score()
    adata.obs["doublet"] = doublets
    adata.obs["doublet_score"] = doublet_score
    sc.pl.violin(adata, ['doublet'], save='doublet.png')
    sc.pl.violin(adata, ['doublet_score'], save='doublet_score.png')
    f = doubletdetection.plot.convergence(clf, save=f'{figures}/convergence_test.pdf', show=True, p_thresh=1e-16, voter_thresh=0.5)
    f3 = doubletdetection.plot.threshold(clf, save='threshold_test.pdf', show=True, p_step=6)
    adata = adata[~(adata.obs['doublet'] !=0)]
    print(adata)
    adata.write(f'{data}/rat_etx_qc_trimmed_raw.gz.h5ad', compression='gzip')
    sns.jointplot(
        data=adata.obs,
        x="log1p_total_counts",
        y="log1p_n_genes_by_counts",
        kind="hex",
    )
    plt.savefig(os.path.join(figures, 'genes_by_counts.png'))

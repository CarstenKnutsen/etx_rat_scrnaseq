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

figures = "/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/figures/tissue_embedding"
data = "/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/processed/single_cell_files"
os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi=300, dpi_save=300, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    # read in processed data
    adata = sc.read(f"{data}/rat_etx_processed_log.gz.h5ad")
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
            'Aspn',
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
    # create empty cell type columns for both northstar and cluster identities
    adata.obs["Cell Subtype"] = pd.Series(index=adata.obs.index, data=None, dtype="str")
    # do this for every lineage
    for lineage in adata.obs["lineage"].unique():
        print(lineage)
        genes = gene_dict[lineage]
        figures_lin = f"{figures}/{lineage}"
        os.makedirs(figures_lin, exist_ok=True)
        sc.settings.figdir = figures_lin
        lin_adata = anndata.AnnData(
            adata[adata.obs["lineage"] == lineage].X,
            obs=adata[adata.obs["lineage"] == lineage].obs,
            var=adata[adata.obs["lineage"] == lineage].var,
        )
        del lin_adata.var["highly_variable"]

        sc.pp.highly_variable_genes(lin_adata, n_top_genes=2000, batch_key="Run")
        sc.pl.highly_variable_genes(lin_adata, save="highly_variable", show=False)
        sc.pp.pca(lin_adata, use_highly_variable=True)
        sc.pp.neighbors(lin_adata)
        print(lineage)
        print(lin_adata.obs["Individual"].value_counts())
        sce.pp.bbknn(
            lin_adata,
            batch_key="Run",
            n_pcs=10,
            neighbors_within_batch=2,
            set_op_mix_ratio=0.5,
        )
        lin_adata.uns["neighbors"]["connectivities"] = lin_adata.obsp["connectivities"]
        leiden_res_dict = {
            "mesenchymal": 0.0002,
            "endothelial": 0.0004,
            "epithelial": 0.0004,
            "immune": 0.0005,
        }
        sc.tl.leiden(
            lin_adata,
            resolution=leiden_res_dict[lineage],
            partition_type=la.CPMVertexPartition,
            key_added=f"leiden_{lineage}",
        )
        # drop low quality clusters
        if lineage == "endothelial":
            lin_adata = lin_adata[
                lin_adata.obs[f"leiden_{lineage}"] != '8'
            ]  # low quality cells with no endo markers, lots of Rpl/Rps and some epi markers
        elif lineage == 'mesenchymal':
            lin_adata = lin_adata[
                lin_adata.obs[f"leiden"] != '52'
                ]  # Low quality E20 cells after iteration, fall in center of umap and have no markers outside of ribosomal genes
        # leiden cluster identity based off of cell type markers
        leiden_ct_dict = {
            "mesenchymal": {
                "0": "Myofibroblast",
                "1": "Alveolar fibroblast",
                "2": "Adventitial fibroblast",
                "3": "Alveolar fibroblast",
                "4": "ASM MyoF precursor",
                "5": "Alveolar fibroblast",
                "6": "Airway smooth muscle",
                "7": "Alveolar fibroblast",
                "8": "Airway smooth muscle",
                "9": "Alveolar fibroblast",
                "10": "Adventitial fibroblast",
                "11": "Alveolar fibroblast",
                "12": "Proliferating fibroblast",
                "13": "Alveolar fibroblast",
                "14": "Alveolar fibroblast",
                "15": "Vascular smooth muscle",
                "16": "Pericyte",
                "17": "Proliferating myofibroblast",
                "18": "Proliferating pericyte",
                "19": "Alveolar fibroblast",
            },
            "endothelial": {
                "0": "Capillary-g",
                "1": "Capillary-g",
                "2": "Capillary-g",
                "3": "Capillary-a",
                "4": "Proliferateive Cap",
                "5": "Arterial EC",
                "6": "Venous EC",
                "7": "Capillary-g",
                "8": "low-quality",
                "9": "Lymphatic EC",
                "10": "Arterial EC",
            },
            "immune": {
                "0": "Alveolar macrophage",
                "1": "Alveolar macrophage",
                "2": "Monocyte",
                "3": "Alveolar macrophage",
                "4": "Neutrophil",
                "5": "Alveolar macrophage",
                "6": "B cell",
                "7": "Intermediate monocyte",
                "8": "Alveolar macrophage",
                "9": "Interstitial macrophage",
                "10": "T cell",
                "11": "Alveolar macrophage",
                "12": "NK cell",
                "13": "Neutrophil",
                "14": "Proliferating macrophage",
                "15": "B cell",
                "16": "Intermediate monocyte",
            },
            "epithelial": {
                "0": "AT1",
                "1": "AT1_AT2",
                "2": "AT2",
                "3": "AT1",
                "4": "AT2",
                "5": "AT2",
                "6": "AT1",
                "7": "AT2",
                "8": "Club",
                "9": "AT1_AT2",
                "10": "AT1",
                "11": "AT1",
                "12": "Club",
                "13": "AT2",
                "14": "AT1",
                "15": "AT1",
                "16": "Ciliated",
                "17": "AT1",
                "18": "AT1",
                "19": "AT2",
                "20": "Club",
                "21": "Basal",
            },
        }
        lin_adata.obs["Cell Subtype"] = (
            lin_adata.obs[f"leiden_{lineage}"].map(leiden_ct_dict[lineage]).astype("str")
        )
        # Clean up some mesenchymal annotation
        if lineage =='mesenchymal':
            cts = []
            for x,y in zip(lin_adata.obs['leiden'], lin_adata.obs['Cell Subtype']):
                if x == '62':
                    cts.append('Mesothelial')
                elif x == '29' and y == 'Proliferating pericyte':
                    cts.append('Proliferating vascular smooth muscle')
                else:
                    cts.append(y)
            lin_adata.obs['Cell Subtype'] = cts
        sc.tl.dendrogram(lin_adata, groupby=f"leiden_{lineage}")
        sc.tl.rank_genes_groups(lin_adata, f"leiden_{lineage}", method="wilcoxon")
        print(lin_adata.obs[f"leiden_{lineage}"].cat.categories)
        sc.pl.rank_genes_groups_dotplot(
            lin_adata,
            groupby=f"leiden_{lineage}",
            n_genes=int(150 / len(lin_adata.obs[f"leiden_{lineage}"].unique())),
            show=False,
            save=f"{lineage}_leiden_markers.png",
        )
        sc.pl.DotPlot(lin_adata, genes, groupby=f"leiden_{lineage}",).style(
            cmap="Reds"
        ).add_totals().savefig(
            os.path.join(figures_lin, f"dotplot_{lineage}_leiden.png"), dpi=300
        )
        plt.clf()
        adata.obs["Cell Subtype"].loc[lin_adata.obs.index] = lin_adata.obs[
            "Cell Subtype"
        ]
        # Visualize cell type connections, generate umap for lineage, save lineage umap to main object
        sc.tl.paga(lin_adata, groups="Cell Subtype")
        sc.pl.paga(
            lin_adata,
            color="Cell Subtype",
            show=False,
            save=f"_{lineage}_Cell Subtype.png",
        )
        b_dict = {"mesenchymal": 1, "endothelial": 0.7, "epithelial": 1, "immune": 1.2}
        sc.tl.umap(lin_adata, a=2, b=b_dict[lineage])
        # Add lineage umaps and leiden clusters to top level
        adata.obs[f"umap_{lineage}_1"] = np.nan
        adata.obs[f"umap_{lineage}_2"] = np.nan
        lin_adata.obs[f"umap_{lineage}_1"] = [x[0] for x in lin_adata.obsm["X_umap"]]
        lin_adata.obs[f"umap_{lineage}_2"] = [x[1] for x in lin_adata.obsm["X_umap"]]
        adata.obs[f"umap_{lineage}_1"].loc[lin_adata.obs.index] = lin_adata.obs[
            f"umap_{lineage}_1"
        ]
        adata.obs[f"umap_{lineage}_2"].loc[lin_adata.obs.index] = lin_adata.obs[
            f"umap_{lineage}_2"
        ]
        adata.obs[f"leiden_{lineage}"] = np.nan
        adata.obs[f"leiden_{lineage}"].loc[lin_adata.obs.index] = lin_adata.obs[
            f"leiden_{lineage}"
        ]
        adata.obsm[f"X_umap_{lineage}"] = adata.obs[
            [f"umap_{lineage}_1", f"umap_{lineage}_2"]
        ].to_numpy()
        del adata.obs[f"umap_{lineage}_1"]
        del adata.obs[f"umap_{lineage}_2"]
        sc.pl.pca_overview(lin_adata, color=f"leiden_{lineage}", show=False, save=True)
        sc.pl.pca_variance_ratio(lin_adata, show=False, save="variance_ratio")
        sc.pl.pca_loadings(
            lin_adata, components="1,2,3,4,5", show=False, save="loadings"
        )
        plt.close()
        plt.clf()
        metadata = f"{figures_lin}/metadata"
        os.makedirs(metadata, exist_ok=True)
        sc.settings.figdir = metadata
        for color in [
            "leiden",
            f"leiden_{lineage}",
            "lineage",
            "Treatment",
            "Run",
            "Individual",
            "n_genes_by_counts",
            "Timepoint",
            "Cell Subtype",
        ]:
            if color in ["leiden",f"leiden_{lineage}"]:
                sc.pl.umap(
                    lin_adata,
                    color=color,
                    alpha=0.5,
                    legend_loc="on data",
                    show=False,
                    save=f"{color}",
                )
            else:
                sc.pl.umap(
                    lin_adata, color=color, alpha=0.5, show=False, save=f"{color}"
                )
        gene_fn = f"{figures_lin}/genes"
        os.makedirs(gene_fn, exist_ok=True)
        sc.settings.figdir = gene_fn
        for gene in genes:
            sc.pl.umap(lin_adata, alpha=0.5, color=gene, show=False, save=f"{gene}")
        sc.settings.figdir = figures_lin
        sc.pl.DotPlot(lin_adata, genes, groupby="Cell Subtype",).add_dendrogram(
            show=True
        ).style(cmap="Reds").add_totals().savefig(
            os.path.join(figures_lin, f"dotplot_{lineage}_Cell Subtype.png"), dpi=300
        )
        plt.clf()

        pd.DataFrame(index=lin_adata.obs_names, data=lin_adata.obsm["X_umap"]).to_csv(
            f"{figures_lin}/{lineage}_umap.csv"
        )
        lin_adata = lin_adata[~lin_adata.obs["Cell Subtype"].isna(), :]
        lin_adata.obs = lin_adata.obs.dropna(axis=1, how="all")
        lin_adata.write(f"{data}/{lineage}_processed_log.gz.h5ad", compression="gzip")
        with pd.ExcelWriter(
            f"{figures_lin}/{lineage}_metadata_counts.xlsx", engine="xlsxwriter"
        ) as writer:
            obs_list = ["Cell Subtype", "Treatment", "Timepoint"]
            num_obs = len(obs_list) + 1
            for ind in range(0, num_obs):
                for subset in itertools.combinations(obs_list, ind):
                    if len(subset) != 0:
                        subset = list(subset)
                        if len(subset) == 1:
                            key = subset[0]
                            lin_adata.obs[key].value_counts().to_excel(
                                writer, sheet_name=key
                            )
                        else:
                            key = "_".join(subset)
                            lin_adata.obs.groupby(subset[:-1])[
                                subset[-1]
                            ].value_counts().to_excel(writer, sheet_name=key[:31])
        sc.tl.rank_genes_groups(lin_adata, "Cell Subtype", method="wilcoxon")
        sc.pl.rank_genes_groups_dotplot(
            lin_adata,
            groupby="Cell Subtype",
            n_genes=int(100 / len(lin_adata.obs["leiden"].unique())),
            show=False,
            save=f"_{lineage}_Cell Subtype_markers.png",
        )
    print('Done witn lineages')
    adata.X = adata.layers["raw"].copy()
    del adata.layers
    adata = adata[~adata.obs['Cell Subtype'].isna()] # drop the low quality clusters from the lineage portion
    print('trimmed')
    print(adata)
    adata.write(f"{data}/etx_rat_processed_cell_typed_raw.gz.h5ad", compression="gzip")
    sc.pp.normalize_total(adata, target_sum = 1e4)
    sc.pp.log1p(adata)
    sc.settings.figdir = f"{figures}/tissue_embedding"


    gene_ls = []
    for k in gene_dict.keys():
        gene_ls = gene_ls + gene_dict[k]
    gene_ls = gene_ls + [
                 'leiden',
                 'lineage',
                 'n_genes_by_counts',
                 'Run',
                 'Treatment',
                 'total_counts',
                 'Timepoint',
                 'doublet',
                 'doublet_score',
                 'pct_counts_mt',
        'Cell Subtype'
                 ]
    for gene in gene_ls:
        if gene == 'leiden':
            sc.pl.umap(adata, legend_loc='on data', color=gene, alpha=0.5, show = False,save=f'_{gene}')
        else:
            sc.pl.umap(adata, color=gene, alpha=0.5, show = False, save=f'_{gene}')
    with pd.ExcelWriter(f'{figures}/tissue_embedding/metadata_counts.xlsx', engine='xlsxwriter') as writer:
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



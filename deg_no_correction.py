import pandas as pd
import os
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import numpy as np
import gseapy as gp

figures = "/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/figures/degs_base"
data = "/home/carsten/alvira_bioinformatics/etx_rat_scrnaseq/data/processed/single_cell_files"
os.makedirs(figures, exist_ok=True)
if __name__ == "__main__":
    # read in processed data
    adata = sc.read(f"{data}/etx_rat_processed_cell_typed_raw.gz.h5ad")
    print(adata)
    adata.var['exclude'] = adata.var_names.str.startswith(("Rps","Rpl",'Mt-'))
    adata =adata[:,adata.var['exclude']==False]
    print(adata)
    # normalize data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # do this for every lineage
    no_comparison = []
    for lineage in adata.obs["lineage"].unique():
        print(lineage)
        figures_lin = f"{figures}/{lineage}"
        os.makedirs(figures_lin, exist_ok=True)
        sc.settings.figdir = figures_lin
        lin_adata = adata[adata.obs["lineage"] == lineage]
        ## broad ct markers
        sc.tl.rank_genes_groups(
            lin_adata,
            "Cell Subtype",
            method="wilcoxon",
            pts=True,
            key_added="rank_genes_groups_Cell Subtype",
        )
        with pd.ExcelWriter(
            f"{figures_lin}/{lineage}_Cell Subtype_markers.xlsx", engine="xlsxwriter"
        ) as writer:
            for ct in lin_adata.obs["Cell Subtype"].cat.categories:
                sc.get.rank_genes_groups_df(
                    lin_adata, key="rank_genes_groups_Cell Subtype", group=ct
                ).to_excel(writer, sheet_name=f"{ct} v rest"[:31])
        ##  ct v every other ct
        figures_lin_comp = f"{figures_lin}/cell_type_comparisons"
        os.makedirs(figures_lin_comp, exist_ok=True)
        for ct in lin_adata.obs["Cell Subtype"].cat.categories:
            print(ct)
            with pd.ExcelWriter(
                f"{figures_lin_comp}/{ct}_deg.xlsx", engine="xlsxwriter"
            ) as writer:
                for ct2 in lin_adata.obs["Cell Subtype"].cat.categories:
                    cts_adata = lin_adata[lin_adata.obs["Cell Subtype"].isin([ct, ct2])]
                    sc.tl.rank_genes_groups(
                        cts_adata,
                        "Cell Subtype",
                        groups=[ct, ct2],
                        method="wilcoxon",
                        pts=True,
                        key_added="rank_genes_groups_Cell Subtype",
                    )
                    sc.get.rank_genes_groups_df(
                        cts_adata, key="rank_genes_groups_Cell Subtype", group=ct
                    ).to_excel(writer, sheet_name=f"{ct} v {ct2}"[:31])

        ##  etx v control
        figures_lin_etx = f"{figures_lin}/etx_v_ctl"
        os.makedirs(figures_lin_etx, exist_ok=True)
        for ct in lin_adata.obs["Cell Subtype"].cat.categories:
            with pd.ExcelWriter(
                f"{figures_lin_etx}/{ct}_deg.xlsx", engine="xlsxwriter"
            ) as writer:
                with pd.ExcelWriter(
                    f"{figures_lin_etx}/{ct}_pathway.xlsx", engine="xlsxwriter"
                ) as writer2:
                    for tp in ["D0", "D14"]:
                        ct_tp_adata = lin_adata[
                            (lin_adata.obs["Timepoint"] == tp)
                            & (lin_adata.obs["Cell Subtype"] == ct)
                        ]
                        try:

                            sc.tl.rank_genes_groups(
                                ct_tp_adata,
                                "Treatment",
                                method="wilcoxon",
                                pts=True,
                                key_added="rank_genes_groups_Treatment",
                            )
                            df = sc.get.rank_genes_groups_df(
                                ct_tp_adata,
                                key="rank_genes_groups_Treatment",
                                group="ETX",
                            )
                            df.to_excel(writer, sheet_name=f"{tp}")
                            path = sc.queries.enrich(
                                {
                                    "up": df["names"].values.tolist()[:100],
                                    "down": df["names"].values.tolist()[-100:],
                                },
                                org="rnorvegicus",
                                gprofiler_kwargs={"no_evidences": False},
                            )
                            path.to_excel(writer2, sheet_name=f"{tp}")
                        except:
                            no_comparison.append(f"{lineage}_{ct}_treatment_{tp}")
                            continue

        ##  timepoints
        figures_lin_etx = f"{figures_lin}/timepoints"
        os.makedirs(figures_lin_etx, exist_ok=True)
        for ct in lin_adata.obs["Cell Subtype"].cat.categories:
            ct_adata = lin_adata[lin_adata.obs["Cell Subtype"] == ct]
            with pd.ExcelWriter(
                f"{figures_lin_etx}/{ct}_deg.xlsx", engine="xlsxwriter"
            ) as writer:
                with pd.ExcelWriter(
                    f"{figures_lin_etx}/{ct}_pathway.xlsx", engine="xlsxwriter"
                ) as writer2:
                    for treat in ["CTL", "ETX"]:
                        ct_lin = ct_adata[ct_adata.obs["Treatment"] == treat]

                        e20_d0 = ct_lin[ct_lin.obs["Timepoint"].isin(["E20", "D0"])]
                        d0_d14 = ct_lin[ct_lin.obs["Timepoint"].isin(["D14", "D0"])]
                        if treat == 'CTL':
                            try:
                                sc.tl.rank_genes_groups(
                                    e20_d0,
                                    "Timepoint",
                                    method="wilcoxon",
                                    pts=True,
                                    key_added="rank_genes_groups_Timepoint",
                                )
                                df = sc.get.rank_genes_groups_df(
                                    e20_d0, key="rank_genes_groups_Timepoint", group="D0"
                                )
                                df.to_excel(writer, sheet_name=f"E20 v D0_{treat}"[:31])
                                path = sc.queries.enrich(
                                    {
                                        "up": df["names"].values.tolist()[:100],
                                        "down": df["names"].values.tolist()[-100:],
                                    },
                                    org="rnorvegicus",
                                    gprofiler_kwargs={"no_evidences": False},
                                )
                                path.to_excel(writer2, sheet_name=f"E20 v D0_{treat}"[:31])
                            except:
                                no_comparison.append(f"{lineage}_{ct}_e20vd0_{treat}")
                                continue
                        try:
                            sc.tl.rank_genes_groups(
                                d0_d14,
                                "Timepoint",
                                method="wilcoxon",
                                pts=True,
                                key_added="rank_genes_groups_Timepoint",
                            )
                            df = sc.get.rank_genes_groups_df(
                                d0_d14, key="rank_genes_groups_Timepoint", group="D14"
                            )
                            df.to_excel(writer, sheet_name=f"D0 v D14_{treat}"[:31])
                            path = sc.queries.enrich(
                                {
                                    "up": df["names"].values.tolist()[:100],
                                    "down": df["names"].values.tolist()[-100:],
                                },
                                org="rnorvegicus",
                                gprofiler_kwargs={"no_evidences": False},
                            )
                            path.to_excel(writer2, sheet_name=f"D0 v D14_{treat}"[:31])
                        except:
                            no_comparison.append(f"{lineage}_{ct}_d0vd14_{treat}")
                            continue
    pd.Series(no_comparison).to_csv(f"{figures}/no comparison.csv")

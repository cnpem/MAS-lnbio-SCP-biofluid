# ;==========================================
# ; Title:  Figure 3 - PG
# ; Author: Fábio Patroni
# ; Date:   03 Fev 2025
# ;==========================================

from itertools import product
from glob import glob
from os import makedirs
from os.path import join as pjoin
from re import sub

from matplotlib.patches import Ellipse
from sklearn.decomposition import PCA
from umap import UMAP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


def preprocess(file):
    data = pd.read_csv(file, sep="\t")
    new_columns = [sub(r"\[\d{1,2}\]\s{1}", "", x)
                   for x in data.columns[1:].values]
    new_columns = [sub(r"^X\.\d+\.{2}", "", x) for x in new_columns]
    new_columns = [sub(".htrms", "", x) for x in new_columns]
    data.rename(columns={old: new for old, new in zip(
        data.columns[1:], new_columns)}, inplace=True)

    gene_names = data["protein"].map(
        lambda x: x.lstrip(";") if x else pd.NA)
    gene_names = gene_names.apply(lambda x: x.split(";")[0])
    data["protein"] = gene_names

    return data


def export_reductions(data, name, col_map, diameter, pg=None,
                      pg_layer=False, base_saida="/tmp"):
    if pg_layer:
        camada_cor = pg_number
    else:
        camada_cor = data.loc[pg, :]

    sizes = np.array(diameter)*10
    min_val = np.min(sizes)
    max_val = np.max(sizes)
    range_size = (max_val - min_val) / 3.5
    boundary1 = round(min_val + range_size)
    boundary2 = round(boundary1 + range_size)
    boundary3 = round(boundary2 + range_size)
    sns.set_theme(style="white", font="Arial", font_scale=2)
    plt.figure(figsize=(12, 8))
    plt.rc('font', family='Arial')
    grafico = plt.scatter(umap_result[:, 0], umap_result[:, 1],
                          s=sizes.tolist(), c=camada_cor, alpha=0.7,
                          edgecolors='black', linewidth=0.5, cmap=col_map)
    cbar = plt.colorbar(grafico, label=pg, shrink=0.7, aspect=20,
                        fraction=0.046, pad=0.12)
    cbar.ax.set_position([0.85, 0.05, 0.05, 0.6])
    plt.title(f'UMAP Projection colored by {pg}')
    plt.xlabel('UMAP Component 1', labelpad=10)
    plt.ylabel('UMAP Component 2', labelpad=10)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    legend_elements = [
        plt.scatter([], [], s=boundary1, c='grey',
                    label=f'Small ({round(boundary1/10)} μm)'),
        plt.scatter([], [], s=boundary2, c='grey',
                    label=f'Medium ({round(boundary2/10)} μm)'),
        plt.scatter([], [], s=boundary3, c='grey',
                    label=f'Large ({round(boundary3/10)} μm)')
    ]
    plt.legend(handles=legend_elements, title='Cell diameter',
               labelspacing=1, title_fontsize=16, bbox_to_anchor=(1.01, 0.9),
               loc='center left', frameon=False, fontsize=14)
    plt.savefig(pjoin(base_saida, f"{name}_umap_norm_fillbefore_{pg}.svg"),
                bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return True


def draw_confidence_ellipse(ax, x, y, group, n_std=2.0, **kwargs):
    """Draws a confidence ellipse based on the PCA data."""
    if len(x) < 2:
        return
    mean_x, mean_y = np.mean(x), np.mean(y)
    lambda_, v = np.linalg.eig(np.cov(x, y))
    lambda_ = np.sqrt(lambda_)
    ellipse = Ellipse(xy=(mean_x, mean_y),
                      width=lambda_[0]*n_std*2, height=lambda_[1]*n_std*2,
                      angle=np.rad2deg(np.arccos(v[0, 0])), **kwargs)
    ax.add_patch(ellipse)


def check_batch(annotation, df, name, hue_order, pca_result=None):
    adata = sc.AnnData(X=df)
    in_common_cols = df.index.intersection(annotation["Run Label"])
    adata.obs['batch'] = [str(x)
                          for x in annotation.loc[in_common_cols, "batch"]]
    biofluid_labels = ["_".join(x.split("_")[0:2])
                       for x in df.index]
    biofluid_labels = [sub(r"^\[\d+\]\s", "", x) for x in biofluid_labels]
    adata.obs['cell_type'] = biofluid_labels

    adata.var = df.columns.to_frame()
    adata.var_names = df.columns

    # calculate lovain
    sc.pp.neighbors(adata, n_pcs=30, random_state=random_state)
    adata.obsm["X_umap"] = umap_result
    sc.tl.leiden(adata, key_added="leiden_res0_25",
                 resolution=0.25, random_state=random_state)
    sc.tl.leiden(adata, key_added="leiden_res0_5",
                 resolution=0.5, random_state=random_state)
    sc.tl.leiden(adata, key_added="leiden_res1",
                 resolution=1.0, random_state=random_state)

    if name == "batch_fluid":
        batch_palette = {keys: values for keys, values in zip(
            hue_order, ["#feedde", "#fdbe85", "#fd8d3c", "#d94701",
                        "#eff3ff", "#bdd7e7", "#6baed6", "#2171b5"])}
    else:
        batch_palette = {keys: values for keys, values in zip(
            adata.obs['batch'].unique(), ["#8fd9c8", "#e8adc4", "#d5ce9c",
                                          "#9cc1ec"])}
    explained_variance_ratio = pca.explained_variance_ratio_[:2] * 100
    pca_result = pca_result
    groups = adata.obs['batch']

    # NOTE Batch effect check
    sns.set_theme(style="white", font="Arial", font_scale=2)
    _, ax = plt.subplots(figsize=(12, 8))
    for group in hue_order:
        idx = groups == group
        x = pca_result[idx, 0]
        y = pca_result[idx, 1]
        ax.scatter(x, y, label=group, color=batch_palette[group],
                   alpha=0.6, s=50)
        draw_confidence_ellipse(ax, x, y, group=group,
                                edgecolor=batch_palette[group],
                                facecolor='none')
    ax.set_xlabel(f'PCA 1 ({explained_variance_ratio[0]:.2f}%)')
    ax.set_ylabel(f'PCA 2 ({explained_variance_ratio[1]:.2f}%)')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.legend(bbox_to_anchor=(1.01, 0.5),
               loc='center left', frameon=False, fontsize=12)
    plt.savefig(
        pjoin(base_output, f"batch_eval{name}.svg"), bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return adata


def plot_leiden(adata, name, col_map):
    sns.set_theme(style="white", font="Arial", font_scale=2)
    plt.figure(figsize=(12, 8))
    sns.scatterplot(x=umap_result[:, 0], y=umap_result[:, 1],
                    hue=adata.obs[name], palette=col_map)
    sns.despine(top=True, right=True)
    plt.title(name)
    plt.xlabel('UMAP Component 1', labelpad=10)
    plt.ylabel('UMAP Component 2', labelpad=10)
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)
    plt.savefig(
        pjoin(base_output, f"leiden{name}.svg"),
        bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()


def export_ranks(tmp_adata, col_map, leiden_name):
    sc.tl.rank_genes_groups(
        tmp_adata, groupby=leiden_name,
        method="wilcoxon", key_added="dea_leiden"
    )
    sc.pl.rank_genes_groups_dotplot(
        tmp_adata, groupby=leiden_name, standard_scale="var",
        n_genes=5, key="dea_leiden", show=False, dendrogram=True,
        save="top5_by_cluster.svg"
    )
    sc.get.rank_genes_groups_df(
        tmp_adata, None,
        key="dea_leiden",
    ).to_csv(pjoin(base_output, "PG_DE.tsv"), sep="\t", index=None)
    sc.pl.rank_genes_groups(tmp_adata, key='dea_leiden', show=False,
                            save="rankPG_by_cluster.svg")

    # filter the DE genes to select for more cluster-specific DE genes
    sc.tl.filter_rank_genes_groups(
        tmp_adata,
        min_in_group_fraction=0.2,
        max_out_group_fraction=0.2,
        key="dea_leiden",
        key_added="dea_leiden_filtered",
    )
    pg_filtered = sc.get.rank_genes_groups_df(
        tmp_adata, None,
        key="dea_leiden_filtered",
    )
    pg_filtered = pg_filtered.dropna().reset_index(drop=True)
    pg_filtered.sort_values(["group", "pvals_adj"], inplace=True)
    pg_filtered.to_csv(pjoin(base_output, "PG_DE_filtered.tsv"),
                       sep="\t", index=None)
    pg_filtered.query("pvals_adj <= 0.05").to_csv(
        pjoin(base_output, "PG_DE_filtered_fdr.tsv"),
        sep="\t", index=None)
    sc.pl.rank_genes_groups_dotplot(
        tmp_adata,
        groupby=leiden_name,
        standard_scale="var", n_genes=5, show=False, dendrogram=True,
        save="top5_by_cluster_filtered.svg", key="dea_leiden_filtered",
    )
    sc.pl.rank_genes_groups(tmp_adata, key='dea_leiden_filtered', show=False,
                            save="rankPG_by_cluster_filtered.svg")

    PGs = []
    for g in pg_filtered.group.unique():
        PGs.append(pg_filtered.query(
            f"group  == '{g}'").reset_index(drop=True).names[0])
    sc.pl.umap(
        tmp_adata, color=[*PGs, leiden_name],
        vmax="p99", legend_loc="on data", frameon=False, cmap=col_map,
        save="umap_topDE_by_cluster.svg", show=False, palette="Set2"
    )

    return PGs


base_path = "/home/fabio/Desktop/output/spectronaut/figura3"
random_state = 2025

base_output = "/tmp/figura3"
sc.settings.figdir = base_output
makedirs(base_output, exist_ok=True)

morfo_path = "/home/fabio/Desktop/output/morphometrics"
morfo = pd.read_csv(pjoin(
    morfo_path, "Morphometrics_tears_saliva.csv"), sep=";")
morfo = morfo.drop("Lagrima", axis=1).set_index(
    "Saliva")
morfo.index = [f"{x[:1]}{int(x[1:]):02d}" for x in morfo.index]

pgs2remove = pd.read_csv((pjoin(base_path, "exclude_pgs_116.csv")))[
    "remove"].tolist()
final_df = preprocess(pjoin(base_path, "saliva-tears_normalized.tsv"))
final_df.set_index("protein", drop=False, inplace=True)
final_df.drop([*pgs2remove], inplace=True)
final_df.index[final_df.index.duplicated()]
final_df = final_df.loc[~final_df.index.duplicated(keep='first')]

hue_order = ["Biofluid_Saliva", "Biofluid_Tears"]
custom_palette = {key: value for key, value in zip(
    hue_order, ["#D9762B", "#3466AC"])}

# NOTE boxplot/violin PG per cell/control
final_melted = final_df.melt(id_vars="protein",
                             var_name="sample", value_name="ms1")
labels_melted = ["_".join(x.split("_")[0:2]) for x in final_melted["sample"]]
labels_melted = [sub(r"^\[\d+\]\s", "", x) for x in labels_melted]
final_melted["groups"] = [sub("Biofluid_", "", x) for x in labels_melted]
final_melted.dropna(axis=0, inplace=True)
final_melted.reset_index(drop=True, inplace=True)
final_melted["well_line"] = [x.split("_")[2].split(".")[
    0] for x in final_melted["sample"]]
final_melted["well_line"] = [
    f"{x[:1]}{int(x[1:]):02d}" for x in final_melted["well_line"]]
final_melted["diameter"] = pd.NA

for i, (poco, tipo) in enumerate(zip(
        final_melted.well_line, final_melted.groups)):
    if tipo == "Saliva":
        col_touse = "Diameter"
    elif tipo == "Tears":
        col_touse = "Diameter.1"
    final_melted.loc[i, "diameter"] = morfo.loc[poco, col_touse]

final_melted["diameter"] = final_melted["diameter"].astype(float)
final_melted.to_csv(pjoin(base_output, "PG_melted.tsv"),
                    sep="\t", index=False)

# NOTE export PCA, UMAP, t-sne add PG number as coloring option
# fill than norm
pgs_de = pjoin(base_output, "pgs_de")
makedirs(pgs_de, exist_ok=True)

# NOTE Batch effect check
annot_glob = glob(
    "/home/fabio/Desktop/output/spectronaut/figura2/annotation_*")
annotation_saliva = pd.read_csv(annot_glob[0], sep="\t")
annotation_saliva["Run Label"].replace(".htrms$", "", regex=True, inplace=True)
annotation_saliva.index = annotation_saliva["Run Label"]
annotation_tears = pd.read_csv(annot_glob[1], sep="\t")
annotation_tears["batch"] = annotation_tears["batch"].replace(
    "^batch", "", regex=True).astype(int)
annotation_tears.index = annotation_tears["Run Label"]
annotation = pd.concat([annotation_saliva, annotation_tears])
annotation["batch"] = pd.Categorical(
    annotation["batch"], categories=range(1, 5, 1))
annotation["fluid"] = [x.split('_')[1] for x in annotation["Run Label"]]
annotation["well_line"] = [x.split("_")[2].split(".")[
    0] for x in annotation["Run Label"]]
annotation["well_line"] = [
    f"{x[:1]}{int(x[1:]):02d}" for x in annotation["well_line"]]
annotation["diameter"] = pd.NA

for sample, poco, tipo in zip(annotation["Run Label"],
                              annotation.well_line, annotation.fluid):
    if tipo == "Saliva":
        col_touse = "Diameter"
    elif tipo == "Tears":
        col_touse = "Diameter.1"
    annotation.loc[sample, "diameter"] = morfo.loc[poco, col_touse]
annotation

annotation2 = annotation.copy()
annotation2["batch"] = ["_".join([str(y), str(x)]) for x, y in zip(
    annotation["batch"], annotation["fluid"])]
tmp_batch_lvl = ["_".join([str(x), str(y)])
                 for x, y in product(["Saliva", 'Tears'], range(1, 5, 1))]
annotation2["batch"] = pd.Categorical(
    annotation2["batch"], categories=tmp_batch_lvl)

data = np.log1p(final_df.drop("protein", axis=1))

pca = PCA(n_components=2, random_state=random_state)
pca_result = pca.fit_transform(data.T)

umap = UMAP(random_state=random_state)
umap_result = umap.fit_transform(data.T)

adata = check_batch(annotation2, data.T, "batch_fluid",
                    tmp_batch_lvl, pca_result)
_ = check_batch(annotation, data.T, "batch", [
    str(x) for x in range(1, 5, 1)], pca_result)

# NOTE plot leiden
for col in ["leiden_res0_25", "leiden_res0_5", "leiden_res1"]:
    pg_number = adata.obs[col].astype(int)
    plot_leiden(adata, col, col_map="Set2")

final_palette = {keys: values for keys, values in zip(
    tmp_batch_lvl, ["#feedde", "#fdbe85", "#fd8d3c", "#d94701",
                    "#eff3ff", "#bdd7e7", "#6baed6", "#2171b5"])}
tmp_incommon = list(set(annotation["Run Label"]).intersection(data.columns))

diameter = annotation.loc[tmp_incommon, "diameter"]
pgs = export_ranks(adata, "Greens", "leiden_res1")
for pg in pgs:
    if len(data.index.intersection([pg])) == 0:
        continue
    export_reductions(data, "Green", "Greens", diameter, pg, base_saida=pgs_de)

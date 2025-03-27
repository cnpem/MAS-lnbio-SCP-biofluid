# ;==========================================
# ; Title:  Figure 2 (PG)
# ; Author: Fábio Patroni
# ; Date:   23 Jan 2025
# ;==========================================

from glob import glob
from os import makedirs
from os.path import basename
from os.path import join as pjoin
from re import sub

from matplotlib.patches import Ellipse
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.manifold import TSNE
from umap import UMAP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


def visualize_high_dim_clusters(df, labels, random_state):
    """
    Visualize clustering results for high-dimensional data using multiple

    Parameters:
    data: array-like of shape (n_samples, n_features)
    labels: array-like of shape (n_samples,)
    random_state: int, random seed for reproducibility
    """

    # color_labels = pd.Series(labels).map(labels_mapper)

    # data = RobustScaler().fit_transform(df)
    data = df

    # Create figure with subplots
    plt.figure(figsize=(30, 15))

    # 1. PCA
    pca = PCA(n_components=2, random_state=random_state)
    pca_result = pca.fit_transform(data)

    plt.subplot(231)
    plt.scatter(pca_result[:, 0], pca_result[:, 1], c="#7398C4")
    plt.title('PCA Projection')

    # plt.legend(handles=[mpatches.Patch(color=y, label=x)
    #            for x, y in labels_mapper.items()])
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')

    # 2. t-SNE
    tsne = TSNE(n_components=2, random_state=random_state)
    tsne_result = tsne.fit_transform(data)

    plt.subplot(232)
    plt.scatter(tsne_result[:, 0], tsne_result[:, 1], c="#7398C4")
    plt.title('t-SNE Projection')

    # 3. UMAP
    umap = UMAP(random_state=random_state)
    umap_result = umap.fit_transform(data)

    plt.subplot(233)
    plt.scatter(umap_result[:, 0], umap_result[:, 1], c="#7398C4")
    plt.title('UMAP Projection')

    # 4. Feature Importance Heatmap (using PCA components)
    plt.subplot(234)
    pca_full = PCA(random_state=random_state)
    pca_full.fit(data.T)

    loadings = pd.DataFrame(
        pca_full.components_[:5, :],  # First 5 components
        columns=range(len(labels))
        # columns=labels
    )
    sns.heatmap(loadings, cmap='RdBu_r', center=0, vmax=1, vmin=-1)
    plt.title('Top PCA Feature Loadings')
    plt.xlabel('Features')
    plt.ylabel('PCA Components')

    # 5. Cluster Size Distribution
    plt.subplot(235)
    unique_labels, counts = np.unique(labels, return_counts=True)
    plt.bar(unique_labels, counts, color="#7398C4")
    plt.title('Cluster Size Distribution')
    plt.xlabel('Cluster Label')
    plt.xticks(rotation=0, ha='center', rotation_mode='anchor')
    plt.ylabel('Number of Samples')

    # 6. Hierarchical Clustering Dendrogram
    plt.subplot(236)
    linkage_matrix = linkage(data, method='ward')
    dendrogram(linkage_matrix, labels=labels)
    plt.title('Hierarchical Clustering Dendrogram\n(sample of data)')
    plt.xlabel('Sample Index')
    plt.ylabel('Distance')
    return True


def preprocess(file):
    data = pd.read_csv(file, sep="\t")
    crap = np.array(data['PG.FastaFiles'].map(lambda x: 'crap' in x.lower()))

    crap_df = data[crap]
    crap_df.to_csv(
        pjoin(base_output,
              f"{crap_df.shape[0]}_crap_removed_{basename(file)}"),
        sep="\t", index=False)

    data = data[~crap]
    data.drop(['PG.FastaFiles'], axis=1, inplace=True)
    try:
        data.drop(['PG.Genes'], axis=1, inplace=True)
        data.drop(['PG.ProteinDescriptions'], axis=1, inplace=True)
    except KeyError:
        print("col not present")
    new_columns = [sub(r"\[\d{1,2}\]\s{1}", "", x)
                   for x in data.columns[1:].values]
    new_columns = [sub(r"^X\.\d+\.{2}", "", x) for x in new_columns]
    new_columns = [sub(".htrms", "", x) for x in new_columns]
    data.rename(columns={old: new for old, new in zip(
        data.columns[1:], new_columns)}, inplace=True)

    gene_names = data["PG.ProteinGroups"].map(
        lambda x: x.lstrip(";") if x else pd.NA)
    gene_names = gene_names.apply(lambda x: x.split(";")[0])
    data["PG.ProteinGroups"] = gene_names

    # remove library
    library_col = np.array(data.columns.map(
        lambda x: 'library' in x.lower()), dtype=bool)

    return data[data.columns[~library_col]]


def boxplot_outliers(df, name):
    sns.set_theme(style="white", font="Arial", font_scale=2)
    df.plot(kind="box", figsize=(18, 8), ylabel='MS1', logy=True,
            color=custom_palette[name])
    plt.xticks(rotation=60, ha='right', rotation_mode='anchor')
    sns.despine(top=True, right=True)
    # plt.legend(bbox_to_anchor=(1.05, 1),
    #            loc='upper left', title="concentration")
    # plt.savefig(pjoin(base_output, f"outliers_{name}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output, f"outliers_{name}.svg"),
                bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()
    return True


def pg_plots(df, x_axis, y_axis, cor_col, name, hue_order,
             base_output, y_label, fluid, caixa=True):
    sns.set_theme(style="white", font="Arial", font_scale=2)
    plt.figure(figsize=(12, 8))
    if caixa:
        sns.boxplot(x=x_axis, y=y_axis, data=df,
                    hue=cor_col, order=hue_order,
                    flierprops={
                        'marker': 'o',
                        'markerfacecolor': 'white',
                        'markeredgecolor': 'black'},
                    palette=custom_palette, dodge=False)
    else:
        sns.violinplot(data=df, x=x_axis, y=y_axis, hue=x_axis,
                       clip_on=False, palette=custom_palette, dodge=False,
                       fill=True, alpha=1, linewidth=1.5, legend="brief",
                       order=hue_order)
    sns.despine(top=True, right=True)
    plt.ylabel(y_label)
    plt.xlabel(fluid)
    plt.xticks(rotation=0, ha='center', rotation_mode='anchor')
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)
    # plt.savefig(pjoin(base_output, f"{name}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output, f"{name}.svg"), bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return True


def export_reductions(data, name, col_map, diameter, pg=None,
                      pg_layer=False, base_saida="/tmp"):
    if pg_layer == "#D9762B" or pg_layer == "#3466AC":
        camada_cor = pg_layer
    elif pg_layer:
        camada_cor = pg_number
    else:
        camada_cor = data.loc[pg, :]

    sns.set_theme(style="white", font="Arial", font_scale=2)
    plt.figure(figsize=(12, 8))
    plt.rc('font', family='Arial')
    grafico = plt.scatter(pca_result[:, 0], pca_result[:, 1],
                          c=camada_cor, alpha=0.7, edgecolors='black',
                          linewidth=0.5, cmap=col_map)
    if pg_layer == "#D9762B" or pg_layer == "#3466AC":
        pass
    else:
        plt.colorbar(grafico, label=pg)
    plt.title(f'PCA Projection colored by {pg}')
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    # plt.savefig(pjoin(base_output, f"{name}_pca_norm_fillbefore_{pg}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_saida, f"{name}_pca_norm_fillbefore_{pg}.svg"),
                bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    # 2D tsne
    sns.set_theme(style="white", font="Arial", font_scale=1.5)
    plt.figure(figsize=(12, 8))
    plt.rc('font', family='Arial')
    grafico = plt.scatter(tsne_result[:, 0], tsne_result[:, 1],
                          c=camada_cor, alpha=0.7, edgecolors='black',
                          linewidth=0.5, cmap=col_map)
    if pg_layer == "#D9762B" or pg_layer == "#3466AC":
        pass
    else:
        plt.colorbar(grafico, label=pg, shrink=0.7, aspect=20)
    plt.title(f't-SNE Projection colored by {pg}')
    plt.xlabel('t-SNE Component 1', labelpad=10)
    plt.ylabel('t-SNE Component 2', labelpad=10)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    # plt.savefig(pjoin(base_output, f"{name}_tsne_norm_fillbefore_{pg}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_saida, f"{name}_tsne_norm_fillbefore_{pg}.svg"),
                bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    # 3D tsne
    sns.set_theme(style="white", font="Arial", font_scale=1.5)
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    plt.rc('font', family='Arial')
    grafico = plt.scatter(tsne_result[:, 0], tsne_result[:, 1],
                          tsne_result[:, 2], c=camada_cor, alpha=0.7,
                          edgecolors='black', linewidth=0.5, cmap=col_map)
    if pg_layer == "#D9762B" or pg_layer == "#3466AC":
        pass
    else:
        plt.colorbar(grafico, label=pg, shrink=0.7, aspect=40,
                     orientation='horizontal', pad=0.1)
    plt.title(f't-SNE 3D Projection colored by {pg}')
    plt.xlabel('t-SNE Component 1', labelpad=10)
    plt.ylabel('t-SNE Component 2', labelpad=10)
    ax.set_zlabel('t-SNE Component 3', labelpad=10)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    ax.invert_yaxis()
    ax.view_init(elev=20, azim=-35)
    plt.savefig(pjoin(base_saida, f"{name}_tsne3D_norm_fillbefore_{pg}.svg"),
                bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    sizes = np.array(diameter)*10
    min_val = np.min(sizes)
    max_val = np.max(sizes)
    range_size = (max_val - min_val) / 3.5
    boundary1 = round(min_val + range_size)
    boundary2 = round(boundary1 + range_size)
    boundary3 = round(boundary2 + range_size)
    sns.set_theme(style="white", font="Arial", font_scale=1.5)
    plt.figure(figsize=(12, 8))
    plt.rc('font', family='Arial')
    grafico = plt.scatter(umap_result[:, 0], umap_result[:, 1],
                          s=sizes.tolist(), c=camada_cor, alpha=0.7,
                          edgecolors='black', linewidth=0.5, cmap=col_map)
    if pg_layer == "#D9762B" or pg_layer == "#3466AC":
        pass
    else:
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


def check_batch(annotation, df, fluid, pca_result=None):
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
    # sc.tl.umap(adata, random_state=random_state)
    adata.obsm["X_umap"] = umap_result
    sc.tl.leiden(adata, key_added="leiden_res0_25",
                 resolution=0.25, random_state=random_state)
    sc.tl.leiden(adata, key_added="leiden_res0_5",
                 resolution=0.5, random_state=random_state)
    sc.tl.leiden(adata, key_added="leiden_res1",
                 resolution=1.0, random_state=random_state)

    batch_palette = {keys: values for keys, values in zip(
        adata.obs['batch'].unique(), ["#8fd9c8", "#e8adc4", "#d5ce9c",
                                      "#9cc1ec"])}
    explained_variance_ratio = pca.explained_variance_ratio_[:2] * 100
    pca_result = pca_result
    name = ""

    groups = adata.obs['batch']

    # NOTE Batch effect check
    sns.set_theme(style="white", font="Arial", font_scale=2)
    _, ax = plt.subplots(figsize=(12, 8))
    for group in groups.unique():
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
        pjoin(base_output, f"batch_eval{name}_{fluid}.svg"),
        bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return adata


def plot_leiden(adata, name, fluid, col_map):
    sns.set_theme(style="white", font="Arial", font_scale=2)
    plt.figure(figsize=(12, 8))
    sns.scatterplot(x=umap_result[:, 0], y=umap_result[:, 1],
                    hue=adata.obs[name], palette=col_map)
    sns.despine(top=True, right=True)
    plt.title(f"{name} - {fluid}")
    plt.xlabel('UMAP Component 1', labelpad=10)
    plt.ylabel('UMAP Component 2', labelpad=10)
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)
    plt.savefig(
        pjoin(base_output, f"leiden{name}_{fluid}.svg"),
        bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()


def export_ranks(tmp_adata, key, col_map):
    sc.tl.rank_genes_groups(
        tmp_adata, groupby="leiden_res0_5",
        method="wilcoxon", key_added="dea_leiden"
    )
    sc.pl.rank_genes_groups_dotplot(
        tmp_adata, groupby="leiden_res0_5", standard_scale="var",
        n_genes=5, key="dea_leiden", show=False, dendrogram=True,
        save=f"top5_by_cluster_{key}.svg"
    )
    sc.get.rank_genes_groups_df(
        tmp_adata, None,
        key="dea_leiden",
    ).to_csv(pjoin(base_output, f"PG_DE_{key}.tsv"), sep="\t", index=None)
    sc.pl.rank_genes_groups(tmp_adata, key='dea_leiden', show=False,
                            save=f"rankPG_by_cluster_{key}.svg")

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
    pg_filtered.to_csv(pjoin(base_output, f"PG_DE_filtered_{key}..tsv"),
                       sep="\t", index=None)
    pg_filtered.query("pvals_adj <= 0.05").to_csv(
        pjoin(base_output, f"PG_DE_filtered_fdr_{key}..tsv"),
        sep="\t", index=None)
    sc.pl.rank_genes_groups_dotplot(
        tmp_adata,
        groupby="leiden_res0_5",
        standard_scale="var", n_genes=5, show=False, dendrogram=True,
        save=f"top5_by_cluster_{key}_filtered.svg", key="dea_leiden_filtered",
    )
    sc.pl.rank_genes_groups(tmp_adata, key='dea_leiden_filtered', show=False,
                            save=f"rankPG_by_cluster_{key}_filtered.svg")

    PGs = []
    for g in pg_filtered.group.unique():
        PGs.append(pg_filtered.query(
            f"group  == '{g}'").reset_index(drop=True).names[0])
    sc.pl.umap(
        tmp_adata, color=[*PGs, "leiden_res0_5"],
        vmax="p99", legend_loc="on data", frameon=False, cmap=col_map,
        save=f"umap_{key}_topDE_by_cluster.svg", show=False, palette="Set2"
    )

    return PGs


base_path = "/home/fabio/Desktop/output/spectronaut/figura2/pg"
# base_path = "/home/fabio.patroni/Downloads/spectronaut/figura_2/pg"
df_globs = glob(pjoin(base_path, "*.tsv"))
random_state = 2025
# np.random.RandomState(random_state)

base_output = "/tmp/figura2"
sc.settings.figdir = base_output
makedirs(base_output, exist_ok=True)

remove_pg_path = pjoin(
    "/home/fabio/Desktop/output/spectronaut", "figura3/exclude_pgs_116.csv")
pgs2remove = pd.read_csv(remove_pg_path)["remove"].tolist()

df_dict = {}
for file in df_globs:
    tmp_df = preprocess(file)
    tmp_df = tmp_df[~tmp_df["PG.ProteinGroups"].isin(pgs2remove)]
    tmp_df.columns = [tmp_df.columns[0]] + \
        [sub(".PG.MS1Quantity", "", x) for x in tmp_df.columns[1:]]
    df_dict[basename(file).replace(".tsv", "")] = tmp_df
df_dict.keys()

# NOTE PCA, UMAP, t-sne - Saliva + Tears + Blank

for i_dataset, (key, df) in enumerate(df_dict.items()):
    final_df = df.copy(deep=True)
    final_df.set_index("PG.ProteinGroups", drop=True, inplace=True)
    final_df.sort_index(inplace=True)
    labels = ["_".join(x.split("_")[0:2]) for x in final_df.columns]
    labels = [sub(r"^\[\d+\]\s", "", x) for x in labels]

    visualize_high_dim_clusters(final_df.fillna(
        0).T, labels, random_state=random_state)
    # plt.savefig(pjoin(base_output, f"dimension_reduction_{key}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output, f"dimension_reduction_{key}.svg"),
                bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

# custom_palette = {key: value for key, value in zip(
#     df_dict.keys(),
#     ["#ABCBFF", "#8EAFE4", "#7398C4", "#5E74A6", "#404D73"])}

# NOTE boxplot/violin PG per cell/control
for i_dataset, (key, df) in enumerate(df_dict.items()):
    if i_dataset == 0:
        final_melted = df.melt(id_vars="PG.ProteinGroups",
                               var_name="sample", value_name="ms1")
        final_melted["cond"] = key
    else:
        tmp_df = df.melt(id_vars="PG.ProteinGroups",
                         var_name="sample", value_name="ms1")
        tmp_df["cond"] = key
        final_melted = pd.concat([final_melted, tmp_df])
labels_melted = ["_".join(x.split("_")[0:2]) for x in final_melted["sample"]]
labels_melted = [sub(r"^\[\d+\]\s", "", x) for x in labels_melted]
final_melted["groups"] = [sub("Biofluid_", "", x) for x in labels_melted]
final_melted.dropna(axis=0, inplace=True)
final_melted.reset_index(drop=True, inplace=True)
final_melted["preparo"] = final_melted["cond"].map(lambda x: x.split("_")[1])
final_melted["preparo"] = final_melted["preparo"].replace(
    "ALL", "before filtering")
final_melted["preparo"] = final_melted["preparo"].replace(
    "filteredwarnings", "after filtering")
hue_order = list(reversed(final_melted["preparo"].unique()))

unicos = final_melted.groupby(["groups", "preparo", "sample"])[
    'PG.ProteinGroups'].nunique().reset_index(name='nunique')
unicos.to_csv(pjoin(base_output, "PG_unicos.tsv"), sep="\t", index=None)

custom_palette = {key: value for key, value in zip(
    hue_order, ["#96C5DE", "#3466AC"])}
pg_plots(unicos.query("groups == 'Tears'"), 'preparo', 'nunique',
         'preparo', 'tears_PG', hue_order,
         base_output, "Avg of Protein Groups", 'Tears')
pg_plots(unicos.query("groups == 'Tears'"), 'preparo', 'nunique',
         'preparo', 'tears_PG_violin', hue_order,
         base_output, "Avg of Protein Groups", 'Tears', caixa=False)

custom_palette = {key: value for key, value in zip(
    hue_order, ["#E1AF7A", "#D9762B"])}
pg_plots(unicos.query("groups == 'Saliva'"), 'preparo', 'nunique', 'preparo',
         'saliva_PG', hue_order, base_output, "Avg of Protein Groups",
         'Saliva')
pg_plots(unicos.query("groups == 'Saliva'"), 'preparo', 'nunique', 'preparo',
         'saliva_PG_violin', hue_order, base_output, "Avg of Protein Groups",
         'Saliva', caixa=False)

# NOTE Corr Diametro das células vs protein groups (lm)

# morfo_path="/home/fabio.patroni/Downloads/scp_jackinho_janeiro/morphometrics"
morfo_path = "/home/fabio/Desktop/output/morphometrics"
morfo = pd.read_csv(pjoin(
    morfo_path, "Morphometrics_tears_saliva.csv"), sep=";")
morfo = morfo.drop("Lagrima", axis=1).set_index(
    "Saliva")
morfo.index = [f"{x[:1]}{int(x[1:]):02d}" for x in morfo.index]

unicos["well_line"] = [x.split("_")[2].split(".")[0] for x in unicos["sample"]]
unicos["well_line"] = [f"{x[:1]}{int(x[1:]):02d}" for x in unicos["well_line"]]
unicos["diameter"] = pd.NA

for i, (poco, tipo) in enumerate(zip(
        unicos.well_line, unicos.groups)):
    if tipo == "Saliva":
        col_touse = "Diameter"
    elif tipo == "Tears":
        col_touse = "Diameter.1"
    unicos.loc[i, "diameter"] = morfo.loc[poco, col_touse]

unicos["diameter"] = unicos["diameter"].astype(float)

model = LinearRegression()
model.fit(np.array(unicos["nunique"]).reshape(-1, 1),
          unicos["diameter"])
y_pred = model.predict(np.array(unicos["nunique"]).reshape(-1, 1))
n = len(unicos["nunique"])
p = 1
residuals = unicos["diameter"] - y_pred
mse = np.sum(residuals**2) / (n - p - 1)
var_x = np.sum(
    (unicos["nunique"] - np.mean(unicos["nunique"]))**2)
se_slope = np.sqrt(mse / var_x)
t_stat = model.coef_[0] / se_slope
p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=n-p-1))
r_sqr = model.score(pd.DataFrame(
    unicos["nunique"]), pd.DataFrame(unicos["diameter"]))

unicos_saliva = unicos.query("groups == 'Saliva'")
unicos_tear = unicos.query("groups == 'Tears'")
sns.set_theme(style="white", font="Arial", font_scale=2)
plt.figure(figsize=(12, 8))
plt.scatter(unicos_saliva["nunique"],
            unicos_saliva["diameter"],
            color='#D9762B', s=100, label='Saliva')
plt.scatter(unicos_tear["nunique"], unicos_tear["diameter"],
            color='#3466AC', s=100, label='Tear')
plt.plot(unicos["nunique"], y_pred, color='black', label=f'R² = {r_sqr:.3f}')
plt.xlabel('Protein Groups')
plt.ylabel('Diameter (µm)')
sns.despine(top=True, right=True)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)
equation = f'y = {model.coef_[0]:.4f}x + {model.intercept_:.2f}'
plt.text(0.8, 0.05, equation, transform=plt.gca().transAxes, fontsize=13,
         bbox=dict(facecolor='white', alpha=0.8), fontname='Arial')
plt.title(
    f"Saliva and Tears Linear Regression (p-value={p_value: .4f})",
    fontname='Arial')
# plt.savefig(pjoin(base_output, "lm_diamter_PG.png"),
#             dpi=300, bbox_inches='tight')
plt.savefig(pjoin(base_output, "lm_diamter_PG.svg"), bbox_inches='tight')
plt.close("all")
sns.reset_defaults()

unicos.to_csv(pjoin(base_output, "LM_PG_unicos.tsv"), sep="\t", index=False)

sns.set_theme(style="white", font="Arial", font_scale=2)
plt.figure(figsize=(12, 8))
sns.heatmap(unicos[["nunique", "diameter"]].corr(), annot=True,
            vmax=1, vmin=-1, cmap="coolwarm", center=0)
plt.title("pairwise correlation method = 'pearson'")
# plt.savefig(pjoin(base_output, "lm_diamter_PG_corr.png"),
#             dpi=300, bbox_inches='tight')
plt.savefig(pjoin(base_output, "lm_diamter_PG_corr.svg"), bbox_inches='tight')
plt.close("all")
sns.reset_defaults()

# NOTE export PCA, UMAP, t-sne add PG number as coloring option
# fill than norm
biomarkers = ["P49748", "P02649", "P00966", "P06276", "P04626", "P11413",
              "P11414", "P11415", "P11416", "P11417", "P11418", "P11419",
              "P11420", "P11421", "P06280", "P04439", "P01911", "O75874",
              "P48735", "P02545", "P06748", "P63151", "Q9Y6L6", "P04637",
              "O14773", "P02766", "P10253"]
biomarkers_saida = pjoin(base_output, "biomarkers_dimention")
makedirs(biomarkers_saida, exist_ok=True)

pgs_de = pjoin(base_output, "pgs_de")
makedirs(pgs_de, exist_ok=True)

abs_saida = pjoin(base_output, "abs")
makedirs(abs_saida, exist_ok=True)

# NOTE Batch effect check

# annot_glob = glob("/home/fabio.patroni/Downloads/spectronaut/annotation_*")
annot_glob = glob(
    "/home/fabio/Desktop/output/spectronaut/figura2/annotation_*")
annotation_saliva = pd.read_csv(annot_glob[0], sep="\t")
annotation_saliva["Run Label"].replace(".htrms$", "", regex=True, inplace=True)
annotation_saliva.index = annotation_saliva["Run Label"]
annotation_tears = pd.read_csv(annot_glob[1], sep="\t")
annotation_tears["batch"] = annotation_tears["batch"].replace(
    "^batch", "", regex=True).astype(int)
annotation_tears.index = annotation_tears["Run Label"]

annotation_dict = {key: value for key, value in zip(
    ['tears_filteredwarnings', 'saliva_filteredwarnings'],
    [annotation_tears, annotation_saliva])}

adata_dict = {key: None for key in [
    'tears_filteredwarnings', 'saliva_filteredwarnings']}

for key in ['tears_filteredwarnings', 'saliva_filteredwarnings']:
    final_df = df_dict[key].copy(deep=True)
    final_df.set_index("PG.ProteinGroups", drop=True, inplace=True)
    final_df.sort_index(inplace=True)

    # data = pd.DataFrame(RobustScaler().fit_transform(
    #     final_df.fillna(0)), columns=final_df.columns, index=final_df.index)
    data = pd.DataFrame(final_df.fillna(
        0), columns=final_df.columns, index=final_df.index)
    data = np.log1p(data)

    pca = PCA(n_components=2, random_state=random_state)
    pca_result = pca.fit_transform(data.T)

    tsne = TSNE(n_components=3, random_state=random_state)
    tsne_result = tsne.fit_transform(data.T)

    umap = UMAP(random_state=random_state)
    umap_result = umap.fit_transform(data.T)

    adata_dict[key] = check_batch(
        annotation_dict[key], data.T, key, pca_result)

    abs_PG = ["P01730", "P01732"]

    if 'tears' in key:
        col_map = 'Blues'
        pg_number_tmp = unicos.query(
            "groups == 'Tears' and preparo == 'after filtering'")
        pg_number = pg_number_tmp.copy(deep=True)
        diameter = pg_number_tmp['diameter']
        pg_number_tmp.set_index(pg_number_tmp["sample"], inplace=True)
        export_reductions(data, key, pg="Tears", diameter=diameter,
                          col_map=None, pg_layer="#3466AC",
                          base_saida=base_output)
        tmp_abs_PG = ["P08575", "P07766", "P15391"]

    else:
        col_map = 'Oranges'
        pg_number_tmp = unicos.query(
            "groups == 'Saliva' and preparo == 'after filtering'")
        pg_number = pg_number_tmp.copy(deep=True)
        diameter = pg_number_tmp['diameter']
        pg_number_tmp.set_index(pg_number_tmp["sample"], inplace=True)
        export_reductions(data, key, pg="Saliva", diameter=diameter,
                          col_map=None, pg_layer="#D9762B",
                          base_saida=base_output)
        tmp_abs_PG = ["P07766", "P15391"]

    pg_number = pg_number_tmp['nunique']
    export_reductions(data, key, pg="Protein Groups", diameter=diameter,
                      col_map=col_map, pg_layer=True, base_saida=base_output)

    # NOTE plot leiden
    for col in ["leiden_res0_25", "leiden_res0_5", "leiden_res1"]:
        pg_number = adata_dict[key].obs[col].astype(int)
        plot_leiden(adata_dict[key], col, key, col_map="Set2")

    for pg in biomarkers:
        if len(data.index.intersection([pg])) == 0:
            print(f"{pg} not in data")
            continue
        export_reductions(data, key, pg=pg, col_map=col_map,
                          diameter=diameter, base_saida=biomarkers_saida)

    tmp_abs_PG = list(set(tmp_abs_PG + abs_PG))
    for pg in tmp_abs_PG:
        if len(data.index.intersection([pg])) == 0:
            print(f"{pg} not in data")
            continue
        export_reductions(data, key, pg=pg, col_map=col_map,
                          diameter=diameter, base_saida=abs_saida)

    pgs = export_ranks(adata_dict[key], key, col_map)
    for pg in pgs:
        if len(data.index.intersection([pg])) == 0:
            continue
        export_reductions(data, key, pg=pg, col_map=col_map,
                          diameter=diameter, base_saida=pgs_de)

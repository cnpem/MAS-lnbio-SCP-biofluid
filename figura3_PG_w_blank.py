# ;==========================================
# ; Title:  Figure 3 - PG w/ blank
# ; Author: Fábio Patroni
# ; Date:   23 Jan 2025
# ;==========================================

from glob import glob
from os import makedirs
from os.path import basename
from os.path import join as pjoin
from re import sub

from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.manifold import TSNE
from sklearn.preprocessing import RobustScaler
from umap import UMAP
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def visualize_high_dim_clusters(df, labels, random_state, labels_mapper):
    """
    Visualize clustering results for high-dimensional data using multiple

    Parameters:
    data: array-like of shape (n_samples, n_features)
    labels: array-like of shape (n_samples,)
    random_state: int, random seed for reproducibility
    """

    color_labels = pd.Series(labels).map(labels_mapper)

    data = RobustScaler().fit_transform(df)
    top_features = pd.Series(
        data.mean(axis=0)).sort_values(ascending=False)[:10]

    # Create figure with subplots
    plt.figure(figsize=(30, 15))

    # 1. PCA
    pca = PCA(n_components=2, random_state=random_state)
    pca_result = pca.fit_transform(data)

    plt.subplot(231)
    plt.scatter(pca_result[:, 0], pca_result[:, 1],
                c=color_labels, cmap='Set1')
    plt.title('PCA Projection')

    plt.legend(handles=[mpatches.Patch(color=y, label=x)
               for x, y in labels_mapper.items()])
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')

    # 2. t-SNE
    tsne = TSNE(n_components=2, random_state=random_state)
    tsne_result = tsne.fit_transform(data)

    plt.subplot(232)
    plt.scatter(tsne_result[:, 0], tsne_result[:, 1],
                c=color_labels, cmap='Set1')
    plt.title('t-SNE Projection')

    # 3. UMAP
    umap = UMAP(random_state=random_state)
    umap_result = umap.fit_transform(data)

    plt.subplot(233)
    plt.scatter(umap_result[:, 0], umap_result[:, 1],
                c=color_labels, cmap='Set1')
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
    sns.heatmap(loadings, cmap='RdBu', center=0)
    plt.title('Top PCA Feature Loadings')
    plt.xlabel('Features')
    plt.ylabel('PCA Components')

    # 5. Cluster Size Distribution
    plt.subplot(235)
    unique_labels, counts = np.unique(labels, return_counts=True)
    plt.bar(unique_labels, counts)
    plt.title('Cluster Size Distribution')
    plt.xlabel('Cluster Label')
    plt.xticks(rotation=30, ha='right', rotation_mode='anchor')
    plt.ylabel('Number of Samples')

    # 6. Hierarchical Clustering Dendrogram
    plt.subplot(236)
    linkage_matrix = linkage(data, method='ward')
    dendrogram(linkage_matrix, labels=labels)
    plt.title('Hierarchical Clustering Dendrogram\n(sample of data)')
    plt.xlabel('Sample Index')
    plt.ylabel('Distance')
    # return fig
    return top_features


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
        data.drop(['PG.Genes', 'PG.ProteinDescriptions'], axis=1, inplace=True)
    except KeyError:
        print("col not present")
    new_columns = [sub(r"\[\d{1,2}\]\s{1}", "", x)
                   for x in data.columns[1:].values]
    data.rename(columns={old: new for old, new in zip(
        data.columns[1:], new_columns)}, inplace=True)

    gene_names = data["PG.ProteinGroups"].map(
        lambda x: x.lstrip(";") if x else pd.NA)
    gene_names = gene_names.apply(lambda x: x.split(";")[0])
    data["PG.ProteinGroups"] = gene_names

    # remove library
    library_col = np.array(data.columns.map(
        lambda x: 'library' in x.lower()))

    return data[data.columns[~library_col]]


def boxplot_outliers(df, name):
    sns.set_theme(style="white", font="Arial", font_scale=2)
    df.plot(kind="box", figsize=(18, 8), ylabel='MS1', logy=True,
            color=custom_palette[name])
    plt.xticks(rotation=60, ha='right', rotation_mode='anchor')
    sns.despine(top=True, right=True)
    # plt.legend(bbox_to_anchor=(1.05, 1),
    #            loc='upper left', title="concentration")
    plt.savefig(pjoin(base_output, f"outliers_{name}.png"),
                dpi=300, bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()
    return True


def pg_plots(df, x_axis, y_axis, hue, name, c_order,
             base_output, y_label, caixa=True):
    sns.set_theme(style="white", font="Arial", font_scale=2)
    plt.figure(figsize=(12, 8))
    if caixa:
        sns.boxplot(x=x_axis, y=y_axis, data=df, order=c_order, legend='brief',
                    hue=hue, palette=custom_palette, dodge=False)
    else:
        sns.violinplot(data=df, x=x_axis, y=y_axis, hue=x_axis,
                       clip_on=False, palette=custom_palette, order=col_order,
                       fill=True, alpha=1, linewidth=1.5, legend="brief")
    sns.despine(top=True, right=True)
    plt.ylabel(y_label)
    plt.xlabel("")
    plt.xticks(rotation=30, ha='right', rotation_mode='anchor')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left', title=x_axis)
    plt.savefig(pjoin(base_output, f"{name}.png"),
                dpi=300, bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return True


base_path = "/home/fabio/Desktop/output/spectronaut/figura2/pg"
df_globs = glob(pjoin(base_path, "*_Report.tsv"))

random_state = 2025
seed = 30
rng = np.random.RandomState(seed)
plot_num = 1

base_output = "/tmp/figura2"
makedirs(base_output, exist_ok=True)

df_dict = {}
for file in df_globs:
    tmp_df = preprocess(file)
    df_dict[basename(file).replace("_Report.tsv", "")] = tmp_df
df_dict.keys()

# blank_1_tears
# blank_1_saliva
# blank do cellenone somente blank 2
# blank_10samples: blank 3

custom_palette = {key: value for key, value in zip(
    ['Biofluid_Saliva', 'blank_1_saliva', 'Biofluid_Tears',
     'blank_3', 'blank_1_tears', 'blank_2'],
    ['#D9762B', '#E1AF7A', '#3466AC',
     '#A6DBA0', '#96C5DE', '#C294CF'])}

# NOTE PCA, UMAP, t-sne - Saliva and Tears
for key, df in df_dict.items():
    # get outliers
    boxplot_outliers(df, key)
    if key == "Biofluid_Saliva" or key == "Biofluid_Tears":
        tmp_df = df.copy(deep=True)
        tmp_df.set_index(tmp_df["PG.ProteinGroups"], drop=True, inplace=True)
        tmp_df.drop("PG.ProteinGroups", axis=1, inplace=True)
        tmp_df.sort_index(inplace=True)
        tmp_labels = ["_".join(x.split("_")[0:2]) for x in tmp_df.columns]
        tmp_labels = [sub(r"^\[\d+\]\s", "", x) for x in tmp_labels]
        visualize_high_dim_clusters(tmp_df.fillna(
            0).T, tmp_labels, random_state=random_state,
            labels_mapper=custom_palette)
        plt.savefig(pjoin(base_output, f"dimension_reduction_{key}.png"),
                    dpi=300, bbox_inches='tight')
        plt.close("all")
        sns.reset_defaults()

# NOTE PCA, UMAP, t-sne - Saliva + Tears

biofluid_df = df_dict['Biofluid_Saliva'].merge(
    df_dict['Biofluid_Tears'], on="PG.ProteinGroups", how="outer")
biofluid_df.set_index(biofluid_df["PG.ProteinGroups"], drop=True, inplace=True)
biofluid_df.drop("PG.ProteinGroups", axis=1, inplace=True)
biofluid_df.sort_index(inplace=True)
biofluid_labels = ["_".join(x.split("_")[0:2]) for x in biofluid_df.columns]
biofluid_labels = [sub(r"^\[\d+\]\s", "", x) for x in biofluid_labels]
visualize_high_dim_clusters(biofluid_df.fillna(
    0).T, biofluid_labels, random_state=random_state)
plt.savefig(pjoin(base_output, "dimension_reduction_biofluid.png"),
            dpi=300, bbox_inches='tight')
plt.close("all")
sns.reset_defaults()

# NOTE PCA, UMAP, t-sne - Blank

tmp_n_file = 0
blank_df = pd.DataFrame()
for key, df in df_dict.items():
    if "lank" in key and tmp_n_file == 1:
        blank_df = blank_df.merge(df, on="PG.ProteinGroups", how="outer")
    if "lank" in key and tmp_n_file == 0:
        blank_df = df.copy(deep=True)
        tmp_n_file = 1
blank_df.set_index(blank_df["PG.ProteinGroups"], drop=True, inplace=True)
blank_df.drop("PG.ProteinGroups", axis=1, inplace=True)
blank_df.sort_index(inplace=True)
blank_labels = ["_".join(x.split("_")[0:2]) for x in blank_df.columns]
blank_labels = [sub(r"^\[\d+\]\s", "", x) for x in blank_labels]
visualize_high_dim_clusters(blank_df.fillna(
    0).T, blank_labels, random_state=random_state)
plt.savefig(pjoin(base_output, "dimension_reduction_blank.png"),
            dpi=300, bbox_inches='tight')
plt.close("all")
sns.reset_defaults()

for i_dataset, (key, df) in enumerate(df_dict.items()):
    if i_dataset == 0:
        final_df = df.copy(deep=True)
    else:
        final_df = final_df.merge(df, on="PG.ProteinGroups", how="outer")

# NOTE boxplot/violin PG per cell/control
final_melted = final_df.melt(
    id_vars="PG.ProteinGroups", var_name="sample", value_name="ms1")
labels_melted = ["_".join(x.split("_")[0:2]) for x in final_melted["sample"]]
labels_melted = [sub(r"^\[\d+\]\s", "", x) for x in labels_melted]
final_melted["groups"] = labels_melted
col_order = pd.Categorical(pd.Series(labels_melted).unique(),
                           categories=['Biofluid_Saliva', 'blank_1_saliva',
                                       'Biofluid_Tears', 'blank_1_tears',
                                       'blank_2', 'blank_3'], ordered=True)
final_melted.dropna(axis=0, inplace=True)
final_melted.reset_index(drop=True, inplace=True)

unicos = final_melted.groupby(["groups", "sample"])[
    'PG.ProteinGroups'].nunique().reset_index(name='nunique')
unicos.to_csv(
    pjoin(base_output, "PG_unicos.tsv"), sep="\t")

pg_plots(unicos, 'groups', 'nunique', 'groups',
         'unique_PG_caixa', col_order, base_output, "N# Unique PG")
pg_plots(unicos, 'groups', 'nunique', 'groups',
         'unique_PG_violin', col_order, base_output, "N# Unique PG",
         caixa=False)

# NOTE Corr Diametro das células vs protein groups (lm)

morfo_path = "/home/fabio.patroni/Downloads/scp_jackinho_janeiro/morphometrics"
morfo = pd.read_csv(pjoin(
    morfo_path, "Morphometrics_tears_saliva.csv"), sep=";")
morfo = morfo.drop("Lagrima", axis=1).set_index(
    "Saliva")
morfo.index = [f"{x[:1]}{int(x[1:]):02d}" for x in morfo.index]
morfo.head()

biofluid_df = df_dict['Biofluid_Saliva'].merge(
    df_dict['Biofluid_Tears'], on="PG.ProteinGroups", how="outer")
biofluid_melted = biofluid_df.melt(
    id_vars="PG.ProteinGroups", var_name="sample", value_name="ms1")
labels_biofluid = ["_".join(x.split("_")[0:2])
                   for x in biofluid_melted["sample"]]
labels_biofluid = [sub(r"^\[\d+\]\s", "", x) for x in labels_biofluid]
biofluid_melted["groups"] = labels_biofluid
biofluid_order = pd.Categorical(pd.Series(labels_biofluid).unique(),
                                categories=[
                                'Biofluid_Saliva', 'Biofluid_Tears'],
                                ordered=True)
biofluid_melted.dropna(axis=0, inplace=True)
biofluid_melted.reset_index(drop=True, inplace=True)

unicos_biofluid = biofluid_melted.groupby(["groups", "sample"])[
    'PG.ProteinGroups'].nunique().reset_index(name='nunique')
unicos_biofluid["well_line"] = [x.split("_")[2].split(
    ".")[0] for x in unicos_biofluid["sample"]]
unicos_biofluid["well_line"] = [
    f"{x[:1]}{int(x[1:]):02d}" for x in unicos_biofluid["well_line"]]
unicos_biofluid["diameter"] = pd.NA

for i, (poco, tipo) in enumerate(zip(
        unicos_biofluid.well_line, unicos_biofluid.groups)):
    if tipo == "Biofluid_Saliva":
        col_touse = "Diameter"
    elif tipo == "Biofluid_Tears":
        col_touse = "Diameter.1"
    unicos_biofluid.loc[i, "diameter"] = morfo.loc[poco, col_touse]
unicos_biofluid.head()

model = LinearRegression()
model.fit(np.array(unicos_biofluid["nunique"]).reshape(-1, 1),
          unicos_biofluid["diameter"])
y_line = model.predict(np.array(unicos_biofluid["nunique"]).reshape(-1, 1))

unicos_biofluid_saliva = unicos_biofluid.query("groups == 'Biofluid_Saliva'")
unicos_biofluid_tear = unicos_biofluid.query("groups == 'Biofluid_Tears'")
sns.set_theme(style="white", font="Arial", font_scale=2)
plt.figure(figsize=(12, 8))
plt.scatter(unicos_biofluid_saliva["nunique"],
            unicos_biofluid_saliva["diameter"],
            color='#b2182b', s=100, label='Saliva')
plt.scatter(unicos_biofluid_tear["nunique"], unicos_biofluid_tear["diameter"],
            color='#2166ac', s=100, label='Tear')
plt.plot(unicos_biofluid["nunique"], y_line, color='black',
         label=f'R² = {model.score(pd.DataFrame(unicos_biofluid["nunique"]), pd.DataFrame(unicos_biofluid["diameter"])):.3f}')
plt.xlabel('Unique PG')
plt.ylabel('Diameter')
sns.despine(top=True, right=True)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)
equation = f'y = {model.coef_[0]:.4f}x + {model.intercept_:.2f}'
plt.title(f'Linear Regression equation {equation}')
plt.savefig(pjoin(base_output, "lm_diamter_PG.png"),
            dpi=300, bbox_inches='tight')
plt.close("all")
sns.reset_defaults()

unicos_biofluid.to_csv(
    pjoin(base_output, "LM_PG_unicos_biofluid.tsv"), sep="\t", index=False)

sns.set_theme(style="white", font="Arial", font_scale=2)
plt.figure(figsize=(12, 8))
sns.heatmap(unicos_biofluid[["nunique", "diameter"]].corr(),
            vmax=1, vmin=-1, cmap="coolwarm", center=0)
plt.title("pairwise correlation method = 'pearson'")
plt.savefig(pjoin(base_output, "lm_diamter_PG_corr.png"),
            dpi=300, bbox_inches='tight')
plt.close("all")
sns.reset_defaults()


# NOTE PCA, UMAP, t-sne - Saliva + Tears + Blank

final_df.set_index(final_df["PG.ProteinGroups"], drop=True, inplace=True)
final_df.drop("PG.ProteinGroups", axis=1, inplace=True)
final_df.sort_index(inplace=True)
labels = ["_".join(x.split("_")[0:2]) for x in final_df.columns]
labels = [sub(r"^\[\d+\]\s", "", x) for x in labels]

visualize_high_dim_clusters(final_df.fillna(
    0).T, labels, random_state=random_state)
plt.savefig(pjoin(base_output, "dimension_reduction_finaldf.png"),
            dpi=300, bbox_inches='tight')
plt.close("all")
sns.reset_defaults()

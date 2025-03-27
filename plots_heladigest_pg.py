# ;==========================================
# ; Title:  Hela digest Analysis (PG)
# ; Author: Fábio Patroni
# ; Date:   08 Jan 2025
# ;==========================================

from os import makedirs
from os.path import basename
from os.path import join as pjoin

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np


def plot_barplot(df):

    ms_groupby = df.groupby(["concentration", "R.Replicate", "ms_type"],
                            observed=False)[
        'value'].agg(["mean", "min", "max", "std"])
    ms_groupby.to_csv(pjoin(base_output, "ms_by_rep.tsv"), sep="\t")

    sns.set_theme(style="white", font="Arial", font_scale=2)
    _, axs = plt.subplots(1, 2, sharex=True, figsize=(24, 8))
    for coluna, ax in zip(["MS1", "MS2"], (axs[0], axs[1])):
        tmp_df = df.query(f"ms_type == '{coluna}'")
        sns.barplot(x="concentration", y='value', data=tmp_df, errorbar='se',
                    hue="concentration", palette=custom_palette,
                    order=col_order, dodge=False, ax=ax, capsize=0.1)
        sns.despine(top=True, right=True, ax=ax)
        ax.set_title("errorbar='se'")
        ax.set_xlabel("")
        ax.set_ylabel(coluna)
    # plt.savefig(pjoin(base_output, "barplot.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output, "barplot.svg"), bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return True


def plot_barplot2(df, y_label):
    plt.figure(figsize=(12, 8))
    sns.set_theme(style="white", font="Arial", font_scale=2)
    sns.barplot(x="concentration", y='nunique', data=df, errorbar=None,
                hue="R.Replicate", palette=rep_palette,
                order=col_order, dodge=True)
    sns.despine(top=True, right=True)
    plt.ylabel(y_label)
    plt.xlabel("")
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left',
               title="Replicate", frameon=False)
    # plt.savefig(pjoin(base_output, "barplot2.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output, "barplot2.svg"), bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return True


def single_barplot(df, y_axis, hue, name, c_order, cmap, ebar=None):
    sns.set_theme(style="white", font="Arial", font_scale=2)
    plt.figure(figsize=(12, 8))
    sns.barplot(x="concentration", y=y_axis, data=df, order=c_order,
                hue=hue, palette=cmap, errorbar=ebar, dodge=False,
                legend=None, capsize=0.1,)
    sns.despine(top=True, right=True)
    plt.ylabel("Avg of Peptide")
    plt.xlabel("")
    # plt.savefig(pjoin(base_output, f"{name}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output, f"{name}.svg"), bbox_inches='tight')
    plt.close("all")
    return True


def caixaplot(df, y_axis, hue, name, c_order,
              base_output, paleta, y_label="CV(%)", lado=False, log_y=False):
    sns.set_theme(style="white", font="Arial", font_scale=2)
    plt.figure(figsize=(12, 8))
    g = sns.boxplot(x="concentration", y=y_axis, data=df, order=c_order,
                    hue=hue, palette=paleta, dodge=lado, legend="brief")
    sns.despine(top=True, right=True)
    plt.ylabel(y_label)
    plt.xlabel('')
    if log_y:
        g.set_yscale("log")
    plt.legend(bbox_to_anchor=(1.01, 1), frameon=False, loc='upper left')
    # plt.savefig(pjoin(base_output, f"{name}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output, f"{name}.svg"), bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return True


def export_std(df, cond_order, nivel_output, nivel="PG.Genes"):
    ms_groupby = df.groupby(["concentration", nivel, "ms_type"],
                            observed=False)[
        'value'].agg(["mean", "min", "max", "std"]).reset_index()
    ms_groupby["CVP"] = round((ms_groupby['std']/ms_groupby['mean'])*100, 3)
    ms_groupby.to_csv(
        pjoin(nivel_output, "ms_key_by_cond.tsv"), sep="\t")

    caixaplot(ms_groupby, 'CVP', 'ms_type', 'MS1_MS2', cond_order,
              base_output, paleta=ms_palette, lado=True, log_y=False)

    ms_groupby.dropna(inplace=True)

    return ms_groupby


def cv_below20(ms1_ms2_df, datatype, gkey='PG.Genes'):

    for level in ["MS1", "MS2"]:
        tmp_df = ms1_ms2_df.query(f"ms_type == '{level}'")
        total_pep = tmp_df[gkey].nunique()
        below20 = tmp_df.query("CVP <= 20")[gkey].nunique()
        type_list = tmp_df.query("CVP <= 20")[gkey].unique().tolist()

        with open(pjoin(
                base_output, f"{datatype}_<=20CV_{level}.txt"), 'w') as f:
            f.write(f"N# {datatype} below 20% CV = {below20}" + '\n')
            f.write(
                f"{round(below20/total_pep, 3)*100}% from total" + '\n')
            f.write("="*79 + '\n')
            f.write("="*39 + f"Full {datatype} list" + '\n')
            f.write("="*79 + '\n')
            f.writelines(str(item) + '\n' for item in type_list)

    return True


def lineplot_concentration(ms1_ms2_df, groupkey, top=True):

    for level in ["MS1", "MS2"]:
        tmp_df = ms1_ms2_df.query(f"ms_type == '{level}'")
        first_df = tmp_df.query("concentration == '50pg'")
        sobe = False if top else True
        first_df = first_df.sort_values("mean", ascending=sobe).iloc[:5, :]
        peps = first_df[groupkey].tolist()
        bool_string = f"' | `{groupkey}` == '".join(peps)
        downstream = tmp_df.query(f"`{groupkey}` == '{bool_string}'")[
            ["concentration", groupkey, "mean"]]
        downstream.reset_index(drop=True, inplace=True)

        file_name = "top5" if top else 'bottom5'

        plt.figure(figsize=(12, 8))
        sns.set_theme(style="white", font="Arial", font_scale=2)
        g = sns.lineplot(data=downstream, x='concentration', y='mean',
                         hue='PG.Genes', palette="Blues")
        g.set_yscale("log")
        sns.despine(top=True, right=True)
        plt.ylabel(f"{level} Mean")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
        # plt.savefig(pjoin(base_output, f"{file_name}_{level}.png"),
        #             dpi=300, bbox_inches='tight')
        plt.savefig(
            pjoin(base_output, f"{file_name}_{level}.svg"),
            bbox_inches='tight')
        plt.close("all")
        sns.reset_defaults()
    return True


def dynamic_range(df, ms_col, groupingkey, name, cor, titulo, together=False):

    if together:
        # mean by replicate by pep
        dinamic_df_tmp = df.groupby([groupingkey, "R.Replicate", "ms_type"])[
            "value"].mean().reset_index()

        # mean by pep
        dinamic_df_final = dinamic_df_tmp.groupby([groupingkey, "ms_type"])[
            "value"].agg(["mean", "std"]).reset_index()

        # MS1
        dinamic_df_final_ms1 = dinamic_df_final.query(
            f"ms_type == '{ms_col[0]}'")
        dinamic_df_final_ms1 = dinamic_df_final_ms1.sort_values(
            "mean", ascending=False).reset_index(drop=True)
        dinamic_df_final_ms1['rank'] = np.arange(
            1, dinamic_df_final_ms1.shape[0] + 1)

        # MS2
        dinamic_df_final_ms2 = dinamic_df_final.query(
            f"ms_type == '{ms_col[1]}'")
        dinamic_df_final_ms2 = dinamic_df_final_ms2.sort_values(
            "mean", ascending=False).reset_index(drop=True)
        dinamic_df_final_ms2['rank'] = np.arange(
            1, dinamic_df_final_ms2.shape[0] + 1)

        sns.set_theme(style="white", font="Arial", font_scale=2)
        plt.figure(figsize=(12, 8))
        dinamic_df_final_ms1 = dinamic_df_final_ms1.sort_values(
            'mean', ascending=False).reset_index(drop=True)
        dinamic_df_final_ms1['rank'] = np.arange(
            1, dinamic_df_final_ms1.shape[0] + 1)
        g = sns.scatterplot(x="rank", y='mean', data=dinamic_df_final_ms1,
                            label="MS1", edgecolor=None,
                            alpha=0.4, color=cor[0])
        g.set_yscale("log")
        dinamic_df_final_ms2 = dinamic_df_final_ms2.sort_values(
            'mean', ascending=False).reset_index(drop=True)
        dinamic_df_final_ms2['rank'] = np.arange(
            1, dinamic_df_final_ms2.shape[0] + 1)
        h = sns.scatterplot(x="rank", y='mean', data=dinamic_df_final_ms2,
                            label="MS2", edgecolor=None,
                            alpha=0.4, color=cor[1])
        h.set_yscale("log")
        sns.despine(top=True, right=True)
        plt.xlabel("Rank")
        plt.ylabel("Mean LFQ intensity")
        plt.title(titulo)
        plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)
        # plt.savefig(pjoin(base_output, "dynamic_range_concat.png"),
        #             dpi=300, bbox_inches='tight')
        plt.savefig(pjoin(base_output, "dynamic_range_concat.svg"),
                    bbox_inches='tight')
        plt.close("all")
        sns.reset_defaults()
    else:
        tmp_df = df.query(f"ms_type == '{ms_col}'")

        # mean by replicate by pep
        dinamic_df_tmp = tmp_df.groupby([groupingkey, "R.Replicate"])[
            "value"].mean().reset_index()

        # mean by pep
        dinamic_df_final = dinamic_df_tmp.groupby([groupingkey])[
            "value"].agg(["mean", "std"]).reset_index()

        dinamic_df_final = dinamic_df_final.sort_values(
            "mean", ascending=False).reset_index(drop=True)
        dinamic_df_final['rank'] = np.arange(1, dinamic_df_final.shape[0] + 1)

        sns.set_theme(style="white", font="Arial", font_scale=2)
        plt.figure(figsize=(12, 8))
        g = sns.scatterplot(x="rank", y='mean', data=dinamic_df_final,
                            edgecolor=None, alpha=0.4, color=cor)
        g.set_yscale("log")
        sns.despine(top=True, right=True)
        plt.xlabel("Rank")
        plt.ylabel(f"Mean LFQ intensity at {name.upper()} level")
        # plt.savefig(pjoin(base_output, f"dynamic_range_{name}.png"),
        #             dpi=300, bbox_inches='tight')
        plt.title(titulo)
        plt.savefig(
            pjoin(base_output, f"dynamic_range_{name}.svg"),
            bbox_inches='tight')
        plt.close("all")
        sns.reset_defaults()

        rank_mapper = {key: value for key, value in zip(
            dinamic_df_final[groupingkey],
            np.arange(1, dinamic_df_final.shape[0] + 1))}
        dinamic_df_tmp["rank"] = dinamic_df_tmp[groupingkey].map(rank_mapper)
        sns.set_theme(style="white", font="Arial", font_scale=2)
        plt.figure(figsize=(12, 8))
        g = sns.scatterplot(x="rank", y='value', data=dinamic_df_tmp,
                            edgecolor=None, alpha=0.4, color=cor)
        sns.despine(top=True, right=True)
        plt.xlabel("Rank")
        plt.ylabel(f"LFQ intensity at {name.upper()} level")
        g.set_yscale("log")
        # plt.savefig(pjoin(
        #     base_output, f"dynamic_range_withReplicates_{name}.png"),
        #     dpi=300, bbox_inches='tight')
        plt.title(titulo)
        plt.savefig(pjoin(
            base_output, f"dynamic_range_withReplicates_{name}.svg"),
            bbox_inches='tight')
        plt.close("all")
        sns.reset_defaults()

    return True


base_path = "/home/fabio/Desktop/output/Hela_digest/combined"
rep_palette = {key: value for key, value in zip(
    range(1, 4, 1), ["#ABCBFF", "#7398C4", "#404D73"])}
# range(1, 4, 1), ["#001e72", "#4ba52f", "#ff6656"])}

file = pjoin(base_path, "PG_combined.tsv")
base_output = pjoin("/tmp", basename(file).replace(".tsv", ""))

makedirs(base_output, exist_ok=True)

df = pd.read_csv(file, sep="\t")
df = df[df["PG.Genes"].notna()].reset_index(drop=True)
df.concentration.replace('1ng', '1000pg', inplace=True)
df.concentration.replace('500ng', '500pg', inplace=True)
df.concentration.replace('50ng', '50pg', inplace=True)
gene_names = df["PG.Genes"].map(lambda x: str(x).lstrip(";"))
df["PG.Genes"] = gene_names.apply(lambda x: x.split(";")[0])
df = df.melt(id_vars=('PG.Genes', 'concentration'), var_name="sample")
df["R.Replicate"] = [int(x.split(" ")[0].strip("[").strip("]"))
                     for x in df["sample"]]
df["ms_type"] = [x.split(".")[-1].replace("Quantity", "")
                 for x in df["sample"]]
df = df.dropna().reset_index(drop=True)
col_order = sorted(df["concentration"].unique(), key=lambda x: int(
    ''.join(filter(str.isdigit, x))))
custom_palette = {key: value for key, value in zip(
    col_order, ["#ABCBFF", "#8EAFE4", "#7398C4", "#5E74A6", "#404D73"])}
df["concentration"] = pd.Categorical(
    df["concentration"], categories=col_order, ordered=True)

ms_palette = {key: value for key, value in zip(
    ['MS1', 'MS2'], ["#ABCBFF", "#404D73"])}
# ['MS1', 'MS2'], ["#73E07E", "#E05C87"])}

# NOTE -Boxplot de cv ms1 e ms2


ms1_ms2_df = export_std(df, col_order, base_output, nivel='PG.Genes')
ms1_ms2_df.query("concentration == '50pg'").groupby(["ms_type"], observed=False)['CVP'].agg(
    ["median", "mean", "min", "max", "std"]).reset_index().to_csv(pjoin(base_output, "ms_metrics.tsv"), sep="\t", index=None)

# NOTE - barplot MS1 e MS2

plot_barplot(df)

# NOTE -Cv below 20 % (quantas proteínas foram quantificadas abaixo de 20%)

cv_below20(ms1_ms2_df, "PG")

# NOTE -Olhar os peptídeos tb nessa curva de concentração(intersecção entre
# as concentrações) - será que é linear? top5 up/down grafico de linha

lineplot_concentration(ms1_ms2_df, 'PG.Genes', top=True)
lineplot_concentration(ms1_ms2_df, 'PG.Genes', top=False)


# NOTE -Numero de PGs por concentração e peptídeos(methdo eval e todos
# juntos - duas abordagens)

unicos = df.groupby(["concentration", "R.Replicate"])[
    'PG.Genes'].nunique().reset_index(name='nunique')
unicos.to_csv(
    pjoin(base_output, "Pep_unicos.tsv"), sep="\t")
single_barplot(unicos, 'nunique', 'concentration',
               'PG_unicos', col_order, custom_palette, ebar='sd')

# NOTE -Stack plot por unique identification de replicata e concentração
plot_barplot2(unicos, "Avg of PG Genes")

custom_cmap = mcolors.LinearSegmentedColormap.from_list(
    'rep_palette', list(rep_palette.values()))
df_pivot = unicos.pivot(index='concentration',
                        columns='R.Replicate', values='nunique')
df_pivot.index = pd.Categorical(
    df_pivot.index, categories=col_order, ordered=True)
df_pivot.sort_index(inplace=True)
sns.set_theme(style="white", font="Arial", font_scale=1.5)
df_pivot.plot(kind='bar', stacked=True,
              ylabel="Cumulative N# Unique PG Genes", colormap=custom_cmap)
plt.xticks(rotation=0)
plt.xlabel("")
sns.despine(top=True, right=True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
           title="Replicate", frameon=False)
# plt.savefig(pjoin(base_output, "stack_count.png"),
#             dpi=300, bbox_inches='tight')
plt.savefig(pjoin(base_output, "stack_count.svg"), bbox_inches='tight')
plt.close("all")
sns.reset_defaults()

distribution = df_pivot.div(np.array(df_pivot.sum(axis=1)).reshape(-1, 1))
distribution.to_csv(pjoin(base_output, "distribution_stack.tsv"), sep="\t")
to_plot = distribution.iloc[:, ::-1].cumsum(axis=1).stack()
to_plot.index.set_names(["concentration", 'R.Replicate'], inplace=True)
sns.set_theme(style="white", font="Arial", font_scale=1.5)
sns.barplot(data=to_plot.reset_index(name='Dist'),
            x='concentration', y='Dist', hue='R.Replicate', dodge=False,
            order=col_order, hue_order=distribution.columns,
            palette=rep_palette)
plt.ylabel("Proportion Unique PG Genes")
plt.xlabel("")
sns.despine(top=True, right=True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
           title="Replicate", frameon=False)
# plt.savefig(pjoin(base_output, "stack_percent.png"),
#             dpi=300, bbox_inches='tight')
plt.savefig(pjoin(base_output, "stack_percent.svg"), bbox_inches='tight')
plt.close("all")
sns.reset_defaults()

# NOTE dinamic range

titulo = "50 pg HeLa digest (directDIA+, three replicates)"
dinamic_df = df.query("concentration == '50pg'").reset_index(drop=True)
dynamic_range(dinamic_df, "MS1", "PG.Genes", "MS1", "#ABCBFF", titulo)
dynamic_range(dinamic_df, "MS2", "PG.Genes", "MS2", "#404D73", titulo)

# dinamic_df = df.query("concentration == '50pg'").reset_index(drop=True)
dynamic_range(dinamic_df, ['MS1', 'MS2'], "PG.Genes",
              "combined", ("#ABCBFF", "#404D73"), titulo, together=True)

# ;==========================================
# ; Title:  Hela digest Analysis (Pep)
# ; Author: Fábio Patroni
# ; Date:   16 Dez 2024
# ;==========================================

from os import makedirs
from os.path import basename
from os.path import join as pjoin

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np


def plot_kde(df, value):
    sns.set_theme(style="white", font="Arial", rc={
                  "axes.facecolor": (0, 0, 0, 0)})
    g = sns.FacetGrid(df, row="concentration", col="R.Replicate",
                      hue="concentration", dropna=True, aspect=12, height=2,
                      palette=custom_palette, margin_titles=True,
                      row_order=col_order)
    #   row_order=row_order)
    g.map(sns.kdeplot, value, clip_on=False, fill=True, alpha=1, linewidth=1.5)
    g.map(sns.kdeplot, value, clip_on=False, color="w", lw=2)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color, fontsize=18,
                ha="left", va="center", transform=ax.transAxes)

    g.map(label, value)
    g.figure.subplots_adjust(hspace=-.6)
    g.set_titles("{col_name}")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    for ax in g.axes.flat:
        ax.axvline(x=df[value].mean(), color='r',
                   linestyle='--', linewidth=0.3)
    # plt.savefig(pjoin(base_output, f"{value}_kde.png"), dpi=300)
    plt.savefig(pjoin(base_output, f"{value}_kde.svg"), bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return True


def plot_violin(df, value):

    sns.set_theme(style="white", font="Arial", rc={
                  "axes.facecolor": (0, 0, 0, 0)}, font_scale=2)
    g = sns.FacetGrid(df, row="concentration", col="R.Replicate",
                      hue="concentration", dropna=True, aspect=8, height=4,
                      palette=custom_palette, margin_titles=True,
                      row_order=col_order)
    g.map(sns.violinplot, value, clip_on=False, order=col_order,
          fill=True, alpha=1, linewidth=1.5)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .5, label, fontweight="bold", color=color, fontsize=23,
                ha="left", va="center", transform=ax.transAxes)

    g.map(label, value)
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    # plt.savefig(pjoin(base_output, f"{value}_violin.png"), dpi=300)
    plt.savefig(pjoin(base_output, f"{value}_violin.svg"), bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return True


def plot_barplot(df, cols, ms1_col="PG.MS1Quantity", ms2_col="PG.MS2Quantity"):

    ms1_groupby = df.groupby(["concentration", "R.Replicate"])[
        ms1_col].agg(["mean", "min", "max", "std"])
    ms1_groupby.to_csv(pjoin(base_output, "ms1_by_rep.tsv"), sep="\t")
    ms2_groupby = df.groupby(["concentration", "R.Replicate"])[
        ms2_col].agg(["mean", "min", "max", "std"])
    ms2_groupby.to_csv(pjoin(base_output, "ms2_by_rep.tsv"), sep="\t")

    sns.set_theme(style="white", font="Arial", font_scale=1.5)
    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(14, 8))
    for coluna, ax in zip(cols, (axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1])):
        sns.barplot(x="concentration", y=coluna, data=df, errorbar='sd',
                    hue="concentration", palette=custom_palette,
                    order=col_order, dodge=False, ax=ax, capsize=0.1)
        sns.despine(top=True, right=True, ax=ax)
        ax.set_xlabel("")
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
              base_output, paleta, y_label="CV(%)", lado=False):
    sns.set_theme(style="white", font="Arial", font_scale=2)
    plt.figure(figsize=(12, 8))
    sns.boxplot(x="concentration", y=y_axis, data=df, order=c_order,
                hue=hue, palette=paleta, dodge=lado, legend="brief")
    sns.despine(top=True, right=True)
    plt.ylabel(y_label)
    plt.xlabel('')
    plt.legend(bbox_to_anchor=(1.05, 1), frameon=False, loc='upper left')
    # plt.savefig(pjoin(base_output, f"{name}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output, f"{name}.svg"), bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return True


def export_std(df, cond_order, nivel_output, nivel='PG.ProteinGroups',
               ms1_col="PG.MS1Quantity", ms2_col="PG.MS2Quantity"):
    ms1_groupby = df.groupby(["concentration", nivel])[ms1_col].agg(
        ["mean", "min", "max", "std"]).reset_index()
    ms1_groupby["CVP"] = round((ms1_groupby['std']/ms1_groupby['mean'])*100, 3)
    ms1_groupby.to_csv(
        pjoin(nivel_output, "ms1_key_by_cond.tsv"), sep="\t")
    ms1_groupby.groupby(["concentration"])["CVP"].describe().to_csv(
        pjoin(nivel_output, "ms1_cvp_cond.tsv"), sep="\t")
    caixaplot(ms1_groupby, 'CVP', 'concentration',
              'MS1', cond_order, nivel_output, paleta=custom_palette)

    ms2_groupby = df.groupby(["concentration", nivel])[ms2_col].agg(
        ["mean", "min", "max", "std"]).reset_index()
    ms2_groupby["CVP"] = round((ms2_groupby['std']/ms2_groupby['mean'])*100, 3)
    ms2_groupby.to_csv(
        pjoin(nivel_output, "ms2_key_by_cond.tsv"), sep="\t")
    ms2_groupby.groupby(["concentration"])["CVP"].describe().to_csv(
        pjoin(nivel_output, "ms2_cvp_cond.tsv"), sep="\t")
    caixaplot(ms2_groupby, 'CVP', 'concentration',
              'MS2', cond_order, nivel_output, paleta=custom_palette)

    ms1_groupby.dropna(inplace=True)
    ms2_groupby.dropna(inplace=True)

    ms1_groupby["cond"] = 'ms1'
    ms2_groupby["cond"] = 'ms2'

    final_df = pd.concat([ms1_groupby, ms2_groupby], axis=0)

    return final_df


def cv_below20(ms1_ms2_df, datatype, level="ms1", gkey='PEP.GroupingKey'):
    tmp_df = ms1_ms2_df.query(f"cond == '{level}'")
    total_pep = tmp_df[gkey].nunique()
    below20 = tmp_df.query("CVP <= 20")[gkey].nunique()
    type_list = tmp_df.query("CVP <= 20")[gkey].unique().tolist()

    with open(pjoin(base_output, f"{datatype}_<=20CV_{level}.txt"), 'w') as f:
        f.write(f"N# {datatype} below 20% CV = {below20}" + '\n')
        f.write(
            f"{round(below20/total_pep, 3)*100}% from total" + '\n')
        f.write("="*79 + '\n')
        f.write("="*39 + f"Full {datatype} list" + '\n')
        f.write("="*79 + '\n')
        f.writelines(str(item) + '\n' for item in type_list)

    return True


def lineplot_concentration(ms1_ms2_df, groupkey, level, top=True):
    tmp_df = ms1_ms2_df.query(f"cond == '{level}'")
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
                     hue='PEP.GroupingKey', palette="muted")
    g.set_yscale("log")
    sns.despine(top=True, right=True)
    plt.ylabel(f"log({level} Mean)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
    # plt.savefig(pjoin(base_output, f"{file_name}_{level}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(
        pjoin(base_output, f"{file_name}_{level}.svg"), bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()
    return True


def dynamic_range(df, ms_col, groupingkey, name, cor, titulo, together=False):

    # mean by replicate by pep
    dinamic_df = df.groupby([groupingkey, "R.Replicate"])[
        ms_col].mean().reset_index()

    # mean by pep
    dinamic_df_final = dinamic_df.groupby([groupingkey])[
        ms_col].agg(["mean", "std"]).reset_index()

    if together:
        sns.set_theme(style="white", font="Arial", font_scale=2)
        plt.figure(figsize=(12, 8))
        dinamic_df_final = dinamic_df_final.sort_values(
            (ms_col[0], 'mean'), ascending=False).reset_index(drop=True)
        dinamic_df_final['rank'] = np.arange(1, dinamic_df_final.shape[0] + 1)
        g = sns.scatterplot(x="rank", y=(
            ms_col[0], 'mean'), data=dinamic_df_final, label="MS1",
            edgecolor=None, alpha=0.4, color=cor[0])
        g.set_yscale("log")
        dinamic_df_final = dinamic_df_final.sort_values(
            (ms_col[1], 'mean'), ascending=False).reset_index(drop=True)
        dinamic_df_final['rank'] = np.arange(1, dinamic_df_final.shape[0] + 1)
        h = sns.scatterplot(x="rank", y=(
            ms_col[1], 'mean'), data=dinamic_df_final, label="MS2",
            edgecolor=None, alpha=0.4, color=cor[1])
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
        dinamic_df["rank"] = dinamic_df[groupingkey].map(rank_mapper)
        sns.set_theme(style="white", font="Arial", font_scale=2)
        plt.figure(figsize=(12, 8))
        g = sns.scatterplot(x="rank", y=ms_col, data=dinamic_df,
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

file = pjoin(base_path, "Pep_combined.tsv")
base_output = pjoin("/tmp", basename(file).replace(".tsv", ""))

makedirs(base_output, exist_ok=True)

df = pd.read_csv(file, sep="\t")
df.concentration.replace('1ng', '1000pg', inplace=True)
df.concentration.value_counts()

df.groupby(["concentration", "R.Replicate"])[
    'PEP.GroupingKey'].nunique().reset_index(name='nunique')

cols_plot = ["PEP.MS1Quantity", 'PEP.MS2Quantity',
             'EG.DatapointsPerPeak (MS1)', 'EG.DatapointsPerPeak']

# pep repetindo
tmp = df.groupby(["R.Replicate", "concentration", 'PEP.GroupingKey'])
tmp = tmp['PEP.GroupingKey'].agg(["count"]).sort_values(
    'count', ascending=False).reset_index()
tmp.sort_values(["R.Replicate", "count"], ascending=[True, False])

plt.figure(figsize=(12, 8))
sns.heatmap(df[["R.Replicate", "EG.DatapointsPerPeak",
            "EG.DatapointsPerPeak (MS1)", "EG.FWHM", "PEP.MS1Quantity",
                "PEP.MS2Quantity"]].corr(),
            vmax=1, vmin=-1, cmap="coolwarm", center=0)
plt.xticks(rotation=30, ha='right', rotation_mode='anchor')
# plt.savefig(pjoin(base_output, "cor_heatmap.png"),
#             dpi=300, bbox_inches='tight')
plt.savefig(pjoin(base_output, "cor_heatmap.svg"), bbox_inches='tight')
plt.close("all")

# NOTE -Boxplot de cv ms1 e ms2

# PEP
col_order = sorted(df["concentration"].unique(), key=lambda x: int(
    ''.join(filter(str.isdigit, x))))
custom_palette = {key: value for key, value in zip(
    col_order, ["#ABCBFF", "#8EAFE4", "#7398C4", "#5E74A6", "#404D73"])}

df["concentration"] = pd.Categorical(
    df["concentration"], categories=col_order, ordered=True)

ms1_ms2_df = export_std(df, col_order, base_output,
                        nivel='PEP.GroupingKey',
                        ms1_col="PEP.MS1Quantity",
                        ms2_col="PEP.MS2Quantity")

ms_palette = {key: value for key, value in zip(
    ['ms1', 'ms2'], ["#ABCBFF", "#404D73"])}
# ['ms1', 'ms2'], ["#73E07E", "#E05C87"])}
caixaplot(ms1_ms2_df, 'CVP', 'cond',
          'MS1_MS2', col_order, base_output, paleta=ms_palette, lado=True)

# NOTE -Data point per peak, MS1 e MS2
plot_kde(df, 'EG.DatapointsPerPeak (MS1)')
caixaplot(df, 'EG.DatapointsPerPeak (MS1)', 'concentration',
          'MS1_peak', col_order, base_output, y_label="PpP",
          paleta=custom_palette)
plot_violin(df, 'EG.DatapointsPerPeak (MS1)')

plot_kde(df, 'EG.DatapointsPerPeak')
caixaplot(df, 'EG.DatapointsPerPeak', 'concentration',
          'MS2_peak', col_order, base_output, y_label="PpP",
          paleta=custom_palette)
plot_violin(df, 'EG.DatapointsPerPeak')

plot_barplot(df, cols_plot, ms1_col="PEP.MS1Quantity",
             ms2_col="PEP.MS2Quantity")

# NOTE -Cv below 20 % (quantas proteínas foram quantificadas abaixo de 20%)

cv_below20(ms1_ms2_df, "pep", "ms1")
cv_below20(ms1_ms2_df, "pep", "ms2")

# NOTE -Olhar os peptídeos tb nessa curva de concentração(intersecção entre
# as concentrações) - será que é linear? top5 up/down grafico de linha

for level in ['ms1', 'ms2']:
    lineplot_concentration(
        ms1_ms2_df, 'PEP.GroupingKey', level=level, top=True)
    lineplot_concentration(
        ms1_ms2_df, 'PEP.GroupingKey', level=level, top=False)

# NOTE -Numero de PGs por concentração e peptídeos(methdo eval e todos
# juntos - duas abordagens)

unicos = df.groupby(["concentration", "R.Replicate"])[
    'PEP.GroupingKey'].nunique().reset_index(name='nunique')
unicos.to_csv(
    pjoin(base_output, "Pep_unicos.tsv"), sep="\t")
single_barplot(unicos, 'nunique', 'concentration',
               'Pep_unicos', col_order, custom_palette, ebar='sd')

# NOTE -Stack plot por unique identification de replicata e concentração
plot_barplot2(unicos, "Avg of Peptide")

distribution = pd.crosstab(
    unicos.concentration, unicos["R.Replicate"], normalize='index')
to_plot = distribution.iloc[:, ::-1].cumsum(axis=1).stack()
sns.set_theme(style="white", font="Arial", font_scale=1.5)
sns.barplot(data=to_plot.reset_index(name='Dist'),
            x='concentration', y='Dist', hue='R.Replicate', dodge=False,
            order=col_order, hue_order=distribution.columns,
            palette=rep_palette)
plt.ylabel("Proportion Unique Peptide")
plt.xlabel("")
sns.despine(top=True, right=True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
           title="Replicate", frameon=False)
# plt.savefig(pjoin(base_output, "stack_percent.png"),
#             dpi=300, bbox_inches='tight')
plt.savefig(pjoin(base_output, "stack_percent.svg"), bbox_inches='tight')
plt.close("all")
sns.reset_defaults()

custom_cmap = mcolors.LinearSegmentedColormap.from_list(
    'rep_palette', list(rep_palette.values()))
df_pivot = unicos.pivot(index='concentration',
                        columns='R.Replicate', values='nunique')
df_pivot.index = pd.Categorical(
    df_pivot.index, categories=col_order, ordered=True)
df_pivot.sort_index(inplace=True)
sns.set_theme(style="white", font="Arial", font_scale=1.5)
df_pivot.plot(kind='bar', stacked=True,
              ylabel="Cumulative N# Unique Peptide", colormap=custom_cmap)
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

# NOTE dinamic range

titulo = "50 pg HeLa digest (directDIA+, three replicates)"
dinamic_df = df.query("concentration == '50pg'").reset_index(drop=True)

dynamic_range(dinamic_df, "PEP.MS1Quantity",
              "PEP.GroupingKey", "ms1", "#ABCBFF", titulo)
dynamic_range(dinamic_df, "PEP.MS2Quantity",
              "PEP.GroupingKey", "ms2", "#404D73", titulo)

dynamic_range(dinamic_df, ["PEP.MS1Quantity", "PEP.MS2Quantity"],
              "PEP.GroupingKey", "combined", ("#ABCBFF", "#404D73"), titulo,
              together=True)

# dinamic_df.groupby(["PEP.GroupingKey"])[
# ["PEP.MS1Quantity", "PEP.MS2Quantity"]].agg(["min", "max", "std"]).describe()

# ;==========================================
# ; Title:  Hela Single Cell Analysis
# ; Author: Fábio Patroni
# ; Date:   15 Jan 2025
# ;==========================================

from glob import glob
from json import dump
from os import makedirs
from os.path import join as pjoin
from re import sub

from matplotlib_venn import venn2
from scipy import stats
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def caixa_id(df, y_axis, y_label, font_size, data_type="PG",
             swarm=True, nocrap="", pad=11):
    plt.figure(figsize=(12, 8))
    sns.set_theme(style="white", font="Arial", font_scale=2)
    if swarm:
        name = "swarm"
        sns.swarmplot(x="cond", y=y_axis, data=df, dodge=True,
                      edgecolor="grey", hue="RunNr", palette=custom_palette,
                      alpha=0.50, legend=False, linewidth=0.8)
    else:
        name = ""
    g = sns.barplot(x="cond", y=y_axis, data=df, errorbar="sd",
                    hue="RunNr", palette=custom_palette,
                    dodge=True, capsize=0.1)
    for i in g.containers:
        g.bar_label(i, fmt='%.0f', padding=pad, fontsize=font_size)
    sns.despine(top=True, right=True)
    plt.ylabel(y_label)
    plt.xlabel("")
    plt.legend(bbox_to_anchor=(1.01, 1),
               loc='upper left', frameon=False)
    # plt.savefig(pjoin(base_output, f"ID_{data_type}_{name}{nocrap}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output,
                      f"ID_{data_type}_{name}{nocrap}_{font_size}.svg"),
                bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()
    return True


def preprocess(file, name):
    data = pd.read_csv(file, sep="\t")
    crap = np.array(data['PG.FastaFiles'].map(lambda x: 'crap' in x.lower()))

    crap_df = data[crap]
    crap_df.to_csv(
        pjoin(base_output,
              f"{crap_df.shape[0]}_crap_removed_{name}.tsv"),
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


def plot_venn(set1, set2, title, name,
              labels_set=('Before', 'After'),
              color_list=["#af8dc3", "#7fbf7b"]):
    sns.set_theme(style="white", font="Arial", font_scale=2.5)
    plt.figure(figsize=(12, 8))
    venn2([set(set1), set(set2)], set_labels=labels_set,
          set_colors=color_list)
    plt.title(title)
    # plt.savefig(pjoin(base_output, f"venn_{name}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output, f"venn_{name}.svg"), bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    esquerda = list(set(set1).difference(set(set2)))
    meio = list(set(set1).intersection(set(set2)))
    direita = list(set(set2).difference(set(set1)))

    return {"esquerda": esquerda, "meio": meio, "direita": direita}


def plot_lm(df, title_name, otimizado='after'):
    # check if data after removing library for each PG
    assert (df.isna().sum(axis=1) == df.shape[0]).sum() == 0

    biofluid_melted = df.melt(
        id_vars="PG.ProteinGroups", var_name="sample", value_name="ms1")
    biofluid_melted.dropna(axis=0, inplace=True)
    biofluid_melted.reset_index(drop=True, inplace=True)

    unicos_biofluid = biofluid_melted.groupby("sample")[
        'PG.ProteinGroups'].nunique().reset_index(name='nunique')
    unicos_biofluid["well_line"] = [x.split("_")[2].split(
        ".")[0] for x in unicos_biofluid["sample"]]
    unicos_biofluid["well_line"] = [
        f"{x[:1]}{int(x[1:]):02d}" for x in unicos_biofluid["well_line"]]
    unicos_biofluid["diameter"] = pd.NA

    for i, poco in enumerate(unicos_biofluid.well_line):
        unicos_biofluid.loc[i, "diameter"] = morfo_df.loc[poco, 'Diameter']

    model = LinearRegression()
    model.fit(np.array(unicos_biofluid["nunique"]).reshape(-1, 1),
              unicos_biofluid["diameter"])
    y_pred = model.predict(np.array(unicos_biofluid["nunique"]).reshape(-1, 1))
    n = len(unicos_biofluid["nunique"])
    p = 1
    residuals = unicos_biofluid["diameter"] - y_pred
    mse = np.sum(residuals**2) / (n - p - 1)
    var_x = np.sum(
        (unicos_biofluid["nunique"] - np.mean(unicos_biofluid["nunique"]))**2)
    se_slope = np.sqrt(mse / var_x)
    t_stat = model.coef_[0] / se_slope
    p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=n-p-1))
    r_sqr = model.score(pd.DataFrame(
        unicos_biofluid["nunique"]), pd.DataFrame(unicos_biofluid["diameter"]))

    sns.set_theme(style="white", font="Arial", font_scale=2)
    plt.figure(figsize=(12, 8))
    plt.scatter(unicos_biofluid["nunique"],
                unicos_biofluid["diameter"], color='#404D73', s=100)
    plt.plot(unicos_biofluid["nunique"], y_pred,
             color='black', label=f'R² = {r_sqr:.3f}')
    plt.xlabel('Protein Groups')
    plt.ylabel(' Diameter (µm)')
    sns.despine(top=True, right=True)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)
    plt.title(f"{title_name} (p-value={p_value: .4f})", fontname='Arial')
    equation = f'y = {model.coef_[0]:.4f}x + {model.intercept_:.2f}'
    plt.text(0.65, 0.05, equation, transform=plt.gca().transAxes, fontsize=13,
             bbox=dict(facecolor='white', alpha=0.8), fontname='Arial')
    # plt.savefig(pjoin(base_output, f"lm_diamter_ME_PG_{otimizado}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output, f"lm_diamter_ME_PG_{otimizado}.svg"),
                bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    unicos_biofluid.to_csv(
        pjoin(base_output, f"LM_ME_PG_unicos_biofluid_{otimizado}.tsv"),
        sep="\t", index=False)

    sns.set_theme(style="white", font="Arial", font_scale=2)
    plt.figure(figsize=(12, 8))
    sns.heatmap(unicos_biofluid[["nunique", "diameter"]].corr(),
                vmax=1, vmin=-1, cmap="coolwarm", center=0, annot=True)
    plt.title("pairwise correlation method = 'pearson'", fontname='Arial')
    # plt.savefig(pjoin(base_output, f"lm_diamter_PG_corr_ME_{otimizado}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output, f"lm_diamter_PG_corr_ME_{otimizado}.svg"),
                bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return True


def caixaplot(df, y_axis, col_cor, name, c_order, h_order,
              base_output, paleta, y_label="CV(%)", lado=False, y_log=False):
    sns.set_theme(style="white", font="Arial", font_scale=2)
    plt.figure(figsize=(12, 8))
    g = sns.boxplot(x="cond", y=y_axis, data=df, order=c_order,
                    hue=col_cor, hue_order=h_order, palette=paleta,
                    dodge=lado, legend="brief")
    sns.despine(top=True, right=True)
    plt.ylabel(y_label)
    plt.xlabel('')
    if y_log:
        g.set_yscale("log")
    else:
        ax = plt.gca()
        ax.yaxis.set_major_locator(plt.MultipleLocator(50))
    plt.legend(bbox_to_anchor=(1.01, 1), frameon=False, loc='upper left')
    # plt.savefig(pjoin(base_output, f"{name}.png"),
    #             dpi=300, bbox_inches='tight')
    plt.savefig(pjoin(base_output, f"{name}.svg"), bbox_inches='tight')
    plt.close("all")
    sns.reset_defaults()

    return True


def remove_ms2(df):
    tmp_df = df.columns[1:].map(lambda x: "ms1quantity" in x.lower()).to_list()
    return df.loc[:, [True]+tmp_df]


def dynamic_range(df, ms_col, groupingkey, name, cor, titulo, together=False):

    if together:
        # mean by replicate by pep
        dinamic_df_tmp = df.groupby([groupingkey, "well_line", "ms_type"])[
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
        dinamic_df_tmp = tmp_df.groupby([groupingkey, "well_line"])[
            "value"].mean().reset_index()

        # mean by pep
        dinamic_df_final = dinamic_df_tmp.groupby([groupingkey])[
            "value"].agg(["mean"]).reset_index()

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
        plt.savefig(pjoin(base_output, f"dynamic_range_{name}.svg"),
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


base_path = "/home/fabio/Desktop/output/Hela_single-cell"
seed = 2025
appro = ["Blanks_after_opt", "Blanks_before_opt"]
after_glob = sorted(glob(pjoin(base_path, appro[0], "*.tsv")))
before_glob = sorted(glob(pjoin(base_path, appro[1], "*.tsv")))
base_output = "/tmp/cluster"
df_globs = before_glob + after_glob

ms_palette = {key: value for key, value in zip(
    ['MS1', 'MS2'], ["#ABCBFF", "#404D73"])}
# ['MS1', 'MS2'], ["#73E07E", "#E05C87"])}

makedirs(base_output, exist_ok=True)

df_dict = {}
for file in df_globs:
    key_name = "_".join(file.split("/")[-2:]).replace(".tsv", "")
    tmp_df = preprocess(file, key_name)
    df_dict[key_name] = tmp_df
df_dict.keys()

# NOTE -Abordagem de busca entre os 4 métodos possíveis no cenário branco sujo
# e branco limpo (comparar PGs e peptídeos)

# com crap
blank_before = pd.read_csv(
    pjoin(base_path, appro[1], "ID_compiladoCRAPinclude.csv"), sep="\t")
blank_before["cond"] = "before"
blank_after = pd.read_csv(
    pjoin(base_path, appro[0], "ID_compilado_CRAPinclude.csv"), sep="\t")
blank_after["cond"] = "after"

all_ids = pd.concat([blank_before, blank_after])
all_ids["cond"] = pd.Categorical(
    all_ids["cond"], categories=["before", "after"], ordered=True)
all_ids['RunNr'].replace("DirectDIA", "directDIA+", inplace=True)
all_ids['RunNr'].replace("library", "library search (10 cells)", inplace=True)
all_ids['RunNr'].replace(
    "matching enhancer", "matching enhancer (10 cells)", inplace=True)

col_order = ['blank', 'method evaluation', 'directDIA+',
             'library search (10 cells)', 'matching enhancer (10 cells)']
all_ids["RunNr"] = pd.Categorical(
    all_ids["RunNr"], categories=col_order, ordered=True)

custom_palette = {key: value for key, value in zip(
    col_order, ["#ABCBFF", "#8EAFE4", "#7398C4", "#5E74A6", "#404D73"])}
# col_order, ["#B0B0AC", "#A8C570", "#7EA6C4", "#DB993B", "#DAA8E3"])}

for i in range(16, 26, 1):
    caixa_id(all_ids, 'ProteinGroups', "Avg of Protein Groups", i, pad=70)
    caixa_id(all_ids, "Peptides", "Avg of Peptides",
             i, data_type="Pep", pad=70)

# sem crap

for i_dataset, (key, df) in enumerate(df_dict.items()):
    # add key to columns
    tmp_df = df.copy(deep=True)
    tmp_df.rename(columns=lambda x: x + "|" +
                  key if x != "PG.ProteinGroups" else x, inplace=True)
    if i_dataset == 0:
        final_df = tmp_df
    else:
        final_df = final_df.merge(
            tmp_df, on="PG.ProteinGroups", how="outer", suffixes=(None, key))
final_melted = final_df.melt(id_vars="PG.ProteinGroups", var_name="sample")
labels_melted = [x.split("|")[1] for x in final_melted["sample"]]
final_melted["RunNr"] = [x.split("_")[3] for x in labels_melted]
final_melted["RunNr"].replace(
    "ME", "matching enhancer (10 cells)", inplace=True)
final_melted["RunNr"].replace("Methodeval", "method evaluation", inplace=True)
final_melted['RunNr'].replace("directDIA", "directDIA+", inplace=True)
final_melted['RunNr'].replace("Blanks", "blank", inplace=True)
final_melted['RunNr'].replace(
    "Library", "library search (10 cells)", inplace=True)
no_crap_order = ['blank', 'method evaluation', 'directDIA+',
                 'library search (10 cells)', 'matching enhancer (10 cells)']
final_melted["RunNr"] = pd.Categorical(
    final_melted["RunNr"], categories=no_crap_order, ordered=True)
final_melted["cond"] = [x.split("_")[1] for x in labels_melted]
final_melted.dropna(axis=0, inplace=True)
final_melted.reset_index(drop=True, inplace=True)
final_melted["ms_type"] = [
    "MS1" if 'ms1quan' in x.lower() else "MS2" for x in final_melted["sample"]]

custom_palette = {key: value for key, value in zip(
    no_crap_order, ["#ABCBFF", "#8EAFE4", "#7398C4", "#5E74A6", "#404D73"])}
# no_crap_order, ["#B0B0AC", "#A8C570", "#7EA6C4", "#DB993B", "#DAA8E3"])}

unicos = final_melted.query("ms_type == 'MS1'")
unicos = unicos.groupby(["cond", "RunNr", "sample"], observed=False)[
    'PG.ProteinGroups'].count().reset_index(name='count')
unicos = unicos.query("count > 0")
unicos["cond"] = pd.Categorical(
    unicos["cond"], categories=["before", "after"], ordered=True)
unicos.to_csv(
    pjoin(base_output, "PG_unicos_noCrap.tsv"), sep="\t", index=None)

for i in range(16, 26, 1):
    caixa_id(unicos, 'count', "Avg of Protein Groups",
             i, nocrap="_noCrap", pad=70)

# NOTE -CV entre as células(MS1) branco sujo(before PO - HIGH) vs limpo(after
# PO - LOW) - Uma abordagem só(talvez a que identificou mais PGs) - ME

me_df = final_melted.query(
    "RunNr == 'matching enhancer (10 cells)'").reset_index(drop=True)
ms_groupby = me_df.groupby(["cond", 'PG.ProteinGroups', 'ms_type'])[
    'value'].agg(["mean", "min", "max", "std"]).reset_index()
ms_groupby["CVP"] = round((ms_groupby['std']/ms_groupby['mean'])*100, 3)
ms_groupby.to_csv(
    pjoin(base_output, "ms_key_by_cond.tsv"), sep="\t")

ms_groupby.groupby(["cond", "ms_type"], observed=False)['CVP'].agg(
    ["median", "mean", "min", "max", "std"]).reset_index().to_csv(pjoin(base_output, "ms_metrics.tsv"), sep="\t", index=False)

caixaplot(ms_groupby, 'CVP', 'ms_type', 'MS1_MS2', ['before', 'after'], [
          'MS1', 'MS2'], base_output, paleta=ms_palette, lado=True,
          y_log=False)

# NOTE -Diametro das células vs protein groups(tanto no branco sujo quanto
# branco limpo, lm) - tabela de diametro está em outro arquivo - indexar pela
# posição do poço

morfo_path = "/home/fabio/Desktop/output/morphometrics"
morfo_df = pd.DataFrame()
for i, file in enumerate(glob(pjoin(morfo_path, "*.tsv"))):
    tmp_df = pd.read_csv(file, sep="\t")
    if i == 1:
        morfo_df = pd.concat([morfo_df, tmp_df])
    if i == 0:
        morfo_df = tmp_df.copy(deep=True)
morfo_df = morfo_df.set_index("Position")
morfo_df.index = [f"{x[:1]}{int(x[1:]):02d}" for x in morfo_df.index]

# remove ms2 before
plot_lm(remove_ms2(df_dict["Blanks_after_opt_ME_PGs"]),
        'Blank after protocol optimization')
plot_lm(remove_ms2(df_dict["Blanks_before_opt_ME_PGs"]),
        'Blank before protocol optimization', "before")


# NOTE -Venn do branco  do que está identificado no branco e nas células,
# intersecção entre o branco sujo e o limpo intercção das PGs identificadas das
# células na condição branco sujo e branco limpo

# all vs all
before_glob_keys = [
    "_".join(x.split("/")[-2:]).replace(".tsv", "") for x in before_glob]
after_glob_keys = ["_".join(x.split("/")[-2:]).replace(".tsv", "")
                   for x in after_glob]
for antes, depois in zip(before_glob_keys, after_glob_keys):
    tmp_before = df_dict[antes]['PG.ProteinGroups']
    tmp_after = df_dict[depois]['PG.ProteinGroups']
    name = antes.split("_")[-2]
    intersection_dict = plot_venn(
        tmp_before, tmp_after, f"{name} vs {name}", f"{name}_{name}",
        color_list=["#ABCBFF", "#5E74A6"])
    with open(pjoin(base_output, f"{name}_{name}.json"), "w") as json_file:
        dump(intersection_dict, json_file, indent=2)
df_dict.keys()

# blank vs me
tmp_before = df_dict["Blanks_after_opt_Blanks_Report"]['PG.ProteinGroups']
tmp_after = df_dict["Blanks_after_opt_ME_PGs"]['PG.ProteinGroups']
blank_ME_dict = plot_venn(tmp_before, tmp_after, "Blank vs ME", "blank_me",
                          labels_set=("Blank", "ME"),
                          color_list=["#ABCBFF", "#5E74A6"])
# color_list = ["#B0B0AC", "#DAA8E3"])

with open(pjoin(base_output, "blank_ME.json"), "w") as json_file:
    dump(blank_ME_dict, json_file, indent=2)

# NOTE dinamic range

me_df = me_df.query("cond == 'after'").reset_index(drop=True)
me_df["well_line"] = [x.split("_")[2].split(".")[0] for x in me_df["sample"]]
me_df["well_line"] = [f"{x[0]}{int(x[1:]):02d}" for x in me_df["well_line"]]
me_df["diameter"] = me_df["well_line"].map(morfo_df.to_dict()['Diameter'])

titulo = "Single HeLa cells after protocol optimization"
dynamic_range(me_df, "MS1", "PG.ProteinGroups", "MS1", "#ABCBFF", titulo)
dynamic_range(me_df, "MS2", "PG.ProteinGroups", "MS2", "#404D73", titulo)

# dinamic_df = df.query("concentration == '50pg'").reset_index(drop=True)
dynamic_range(me_df, ['MS1', 'MS2'], "PG.ProteinGroups",
              "combined", ("#ABCBFF", "#404D73"), titulo, together=True)

tmp_df = me_df.query("ms_type == 'MS1'")
tmp_df = tmp_df.groupby(["sample", "diameter"])[
    'PG.ProteinGroups'].nunique().reset_index(name='nunique')
tmp_df = tmp_df.sort_values(
    "nunique", ascending=False).reset_index(drop=True)
tmp_df['rank'] = np.arange(1, tmp_df.shape[0] + 1).astype(int)
sns.set_theme(style="white", font="Arial", font_scale=2.5)
plt.figure(figsize=(12, 8))
g = sns.scatterplot(x="rank", y='nunique', data=tmp_df, sizes=(69, 1050),
                    edgecolor=None, alpha=0.4, size="diameter")
plt.xlabel("Cell Rank Order")
plt.ylabel("N# PG")
sns.despine(top=True, right=True)
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left',
           frameon=False, title="HeLa cell diameter")
# plt.savefig(pjoin(base_output, "diameter_PG_rank.png"),
#             dpi=300, bbox_inches='tight')
plt.savefig(pjoin(base_output, "diameter_PG_rank.svg"), bbox_inches='tight')
plt.close("all")
sns.reset_defaults()

# Ms1 total

tmp_df = me_df.query("ms_type == 'MS1'")
tmp_df = tmp_df.groupby(["sample", "diameter"])['value'].mean().reset_index()
tmp_df = tmp_df.sort_values("value", ascending=False).reset_index(drop=True)
tmp_df['rank'] = np.arange(1, tmp_df.shape[0] + 1).astype(int)
sns.set_theme(style="white", font="Arial", font_scale=2.5)
plt.figure(figsize=(12, 8))
g = sns.scatterplot(x="rank", y='value', data=tmp_df, sizes=(69, 1050),
                    edgecolor=None, alpha=0.4, size="diameter")
plt.xlabel("Cell Rank Order ")
plt.ylabel("Mean MS1 Quantity")
sns.despine(top=True, right=True)
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left',
           frameon=False, title="HeLa cell diameter")
# plt.savefig(pjoin(base_output, "diameter_MS1_rank.png"),
#             dpi=300, bbox_inches='tight')
plt.savefig(pjoin(base_output, "diameter_MS1_rank.svg"), bbox_inches='tight')
plt.close("all")
sns.reset_defaults()

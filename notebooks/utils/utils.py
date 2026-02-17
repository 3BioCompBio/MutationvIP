import math
from .constants import DF_ORDER, TYPE_ORDER
import polars as pl
import pandas as pd
from viper import Viper

from scipy import stats
import numpy as np
from tabulate import tabulate
import matplotlib.pyplot as plt
import seaborn as sns


VIPER = Viper()


def get_pyrimidine_strand(seq):
    mapper = {"A": "T", "C": "G", "G": "C", "T": "A"}
    if seq[math.floor(len(seq) / 2)] in ["T", "C"]:
        return seq
    new_seq = ""
    for nucl in seq[::-1]:
        new_seq += mapper.get(nucl)
    return new_seq


def get_table_of_corr(
    data,
    x="vip",
    y="F",
    with_len=False,
    categories=("db", DF_ORDER),
    sub_categories=("type", TYPE_ORDER),
    corr_method=stats.pearsonr,
    early_return=False,
):
    """
    categories : tuple
                      key on which split the df
                      sorted list of categories (col for output)
    sub_categories : tuple
                      key on which split the df
                      sorted list of categories (row for output)
    """

    data = data.drop_nulls()
    sort_key = {key: i for i, key in enumerate(sub_categories[1])}

    sort_cat_key = {cat: i for i, cat in enumerate(categories[1])}

    summary_data = {
        categories[0]: [],
        "correlation": [],
        "pvalue": [],
        sub_categories[0]: [],
        "count": [],
    }

    for cat in categories[1]:
        for to_check in sub_categories[1]:
            tmp_cat = data.filter(
                pl.col(categories[0]) == cat, pl.col(sub_categories[0]) == to_check
            )
            try:
                s = corr_method(tmp_cat[y], tmp_cat[x])
            except Exception as e:
                print(f"Error calculating correlation for {cat} and {to_check}: {e}")
                continue
            summary_data[categories[0]].append(cat)
            summary_data["correlation"].append(s.statistic)
            summary_data["pvalue"].append(s.pvalue)
            summary_data["count"].append(tmp_cat["count"].sum())
            summary_data[sub_categories[0]].append(to_check)

    if early_return:
        return summary_data
    unique_cat = [sub_categories[0]] + sorted(
        categories[1], key=lambda x: sort_cat_key.get(x)
    )
    df = (
        pl.DataFrame(summary_data)
        .drop_nulls()
        .pivot(values="correlation", index=[sub_categories[0]], on=categories[0])
        .select(unique_cat)
        .with_columns(
            pl.col(sub_categories[0]).replace_strict(sort_key).alias("sort_key")
        )
        .sort(by="sort_key")
        .drop("sort_key")
    )
    if with_len:
        df.drop("count")
    return df


def aggregate_data(df, kmer_context, on=["type", "context", "db"], over=["type", "db"]):
    """
    kmer_context : df whith at least 'context' and 'all_f' columns
    on = on chich columns we group the values
    over = which columns are use to compute the frequency of occurence of the motifs
    """
    return (
        df.with_columns(
            pl.col("small_context")
            .map_elements(get_pyrimidine_strand, return_dtype=pl.String)
            .alias("context"),
        )
        .group_by(on)
        .agg(pl.col("len").sum().alias("count"))
        .with_columns(
            pl.col("context")
            .map_elements(VIPER.compute_double_strand_score, return_dtype=pl.Float64)
            .alias("vip"),
            (pl.col("count") / pl.col("count").sum().over(over)).alias("freq"),
        )
        .join(kmer_context["context", "all_f"], on="context")
        .with_columns(np.log(pl.col("freq") / pl.col("all_f")).alias("F"))
    )


def print_markdown(df, latex=False, decimals=2):
    """
    Print a dataframe as a markdown or latex table.

    Parameters:
    -----------
    df : polars.DataFrame
        The dataframe to print
    latex : bool, default=False
        If True, output in latex format, otherwise markdown
    decimals : int, default=2
        Number of decimal places for floating point columns
    """
    # Convert to pandas and format numeric columns
    df_pandas = df.to_pandas()

    # Format all numeric columns to specified decimal places
    for col in df_pandas.columns:
        if df_pandas[col].dtype in ["float64", "float32", "float16"]:
            df_pandas[col] = df_pandas[col].apply(
                lambda x: f"{x:.{decimals}f}" if not pd.isna(x) else x
            )

    md_table = tabulate(
        df_pandas,
        headers="keys",
        tablefmt="latex" if latex else "github",
        showindex=False,
    )

    print(md_table)


def plot_regression_plot(data, x, y, **kwargs):
    hue = kwargs.get("hue", None)
    fig, ax = plt.subplots()
    palette = kwargs.get("palette", sns.color_palette())
    if hue is not None:
        unique_hues = data[hue].unique()
        colors = dict(zip(unique_hues, palette))
        for i, val in enumerate(unique_hues):
            subset = data[data[hue] == val]
            plt.scatter(subset[x], subset[y], label=str(val), color=colors[val])
            if len(subset[x]) > 1:
                slope, intercept = np.polyfit(subset[x], subset[y], 1)
                regression_line = slope * subset[x] + intercept
                plt.plot(subset[x], regression_line, color=colors[val])
                r, p = stats.pearsonr(subset[x], subset[y])
                plt.text(
                    0.95,
                    0.95 - i * 0.12,
                    f"{val}: ρ={r:.2f}, p={p:.1e}, b={slope:.2f}",
                    verticalalignment="top",
                    transform=plt.gca().transAxes,
                    bbox=dict(facecolor="white", alpha=0.5, edgecolor="black"),
                    color=colors[val],
                )
        plt.legend(title=hue)
    else:
        plt.scatter(data[x], data[y])
        r, p = stats.pearsonr(data[x], data[y])
        slope, intercept = np.polyfit(data[x], data[y], 1)
        regression_line = slope * data[x] + intercept
        plt.plot(data[x], regression_line, color=palette[1])
        plt.text(
            0.95,
            0.95,
            f"ρ = {r:.2f}\np = {p:.1e}\nb = {slope:.2f}",
            verticalalignment="top",
            transform=plt.gca().transAxes,
            bbox=dict(facecolor="white", alpha=0.5, edgecolor="black"),
        )
    return fig, ax
    # plt.scatter(data[x], data[y])
    # r, p = stats.pearsonr(data[x], data[y])

    # slope, intercept = np.polyfit(data[x], data[y], 1)
    # regression_line = slope * data[x] + intercept

    # # Plot the regression line
    # plt.plot(data[x], regression_line, color=sns.color_palette()[1])

    # plt.text(
    #     0.95,
    #     0.95,
    #     f"ρ = {r:.2f}\np = {p:.1e}\nb = {slope:.2f}",
    #     verticalalignment="top",
    #     transform=plt.gca().transAxes,
    #     bbox=dict(facecolor="white", alpha=0.5, edgecolor="black"),
    # )

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


def partial_correlation(df, x_col, y_col, control_col, stat_method=stats.pearsonr):
    """
    Partial correlation between x and y controlling for control.
    Formula: r(x,y|z) = (r(x,y) - r(x,z)*r(y,z)) / sqrt((1-r(x,z)²)(1-r(y,z)²))
    stat_method : callable with signature f(x, y) -> (statistic, pvalue), e.g. stats.pearsonr or stats.spearmanr
    """
    from scipy.stats import t as t_dist

    x = df[x_col].to_numpy()
    y = df[y_col].to_numpy()
    z = df[control_col].to_numpy()

    r_xy, _ = stat_method(x, y)
    r_xz, _ = stat_method(x, z)
    r_yz, _ = stat_method(y, z)

    numerator = r_xy - r_xz * r_yz
    denominator = np.sqrt((1 - r_xz**2) * (1 - r_yz**2))

    if denominator == 0:
        return np.nan, 1.0

    r_partial = numerator / denominator
    n = len(x)
    t_stat = r_partial * np.sqrt(n - 3) / np.sqrt(1 - r_partial**2)
    p_value = 2 * (1 - t_dist.cdf(abs(t_stat), n - 3))

    return r_partial, p_value


def aggregate_data(
    df,
    kmer_context,
    on=["type", "context", "db"],
    over=["type", "db"],
    context_col="context",
):
    """
    df : dataframe with at least a context column and a 'len' column
    kmer_context : dataframe with at least 'context' and 'frequencies' columns
    on : columns to group by when aggregating counts
    over : columns used to compute the frequency of occurrence of the motifs
    context_col : column in df used as the raw context (default: 'context')
    """
    return (
        df.with_columns(
            pl.col(context_col)
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
        .join(kmer_context["context", "frequencies"], on="context")
        .with_columns(np.log(pl.col("freq") / pl.col("frequencies")).alias("F"))
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

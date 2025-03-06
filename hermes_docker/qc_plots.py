import pandas as pd
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from qmplot import manhattanplot, qqplot


def manhattan_plot(datafile, outputdir, filename="manhattan.png"):
    columns_and_types = {
        'varId': 'str',
        'chromosome': 'str',
        'position': 'float64',
        'pValue': 'float64'
    }
    df = pd.read_table(datafile, sep="\t", usecols=columns_and_types.keys(), dtype=columns_and_types)
    df = df.dropna(how="any", axis=0)
    chromosome_order = {str(i): i for i in range(1, 23)}
    chromosome_order.update({
        'X': 23,
        'Y': 24,
        'XY': 25,
        'MT': 26
    })

    def chromosome_sort_key(chromosome):
        return chromosome_order.get(chromosome, float('inf'))

    df['chromosome_sort'] = df['chromosome'].apply(chromosome_sort_key)
    df_sorted = df.sort_values(by=['chromosome_sort', 'position'])

    f, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    manhattanplot(data=df_sorted, chrom="chromosome", pos="position", pv="pValue", snp="varId", ax=ax, is_annotate_topsnp=True)
    plt.savefig(f'{outputdir}/{filename}.png', dpi=300)
    plt.close()


def qq_plot(datafile, outputdir, filename="qqplot"):
    df = pd.read_table(datafile, sep="\t", usecols=['pValue'])
    df = df.dropna(how="any", axis=0)
    f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
    qqplot(df["pValue"], ax=ax)
    plt.savefig(f'{outputdir}/{filename}.png', dpi=300)
    plt.close()


def pz_plot(datafile, outputdir, filename='pz_plot'):
    df = pd.read_table(datafile, sep="\t", usecols=['beta', 'pValue', 'stdErr'])
    df = df.dropna(how="any", axis=0)
    df['observed'] = -np.log10(df['pValue'])
    expected_p = 2 * norm.sf(np.abs(df['beta'] / df['stdErr']))
    df["expected"] = -np.log10(expected_p)
    df = df.sort_values("expected")

    plt.figure(figsize=(6, 6))
    plt.scatter(df['expected'], df['observed'], c='darkblue', alpha=0.75, s=10)
    plt.plot([0, df["expected"].max()], [0, df["expected"].max()], linestyle="--", color="red", label="y=x")
    plt.xlabel("-log10(p.ztest)")
    plt.ylabel("-log10(p)")
    plt.savefig(f'{outputdir}/{filename}.png', dpi=300)
    plt.close()


def eaf_plot(datafile, outputdir, filename='eaf_plot'):
    df = pd.read_table(datafile, sep="\t", usecols=['eaf', 'g1000_eaf'])
    df = df.dropna(how="any", axis=0)

    plt.figure(figsize=(6, 6))
    plt.scatter(df['eaf'], df['g1000_eaf'], alpha=0.5, s=10, color='darkblue')
    plt.plot([0, 1], [0, 1], color='red', linestyle='--')
    plt.xlabel('Cohort Frequency')
    plt.ylabel('1000G Allele Frequency')
    plt.savefig(f'{outputdir}/{filename}.png', dpi=300)
    plt.close()

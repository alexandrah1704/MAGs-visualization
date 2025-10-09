import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
import numpy as np
import os
import pandas as pd
import seaborn as sns

prefix_map = {
    'domain': 'd',
    'phylum': 'p',
    'class': 'c',
    'order': 'o',
    'family': 'f',
    'genus': 'g',
    'species': 's',
}

def extract_rank(classification, rank):
        prefix = prefix_map.get(rank)
        try:
            return next(x for x in classification.split(';') if x.startswith(f'{prefix}__')).replace(f'{prefix}__', '').strip()
        except StopIteration:
            return None

def completeness_contamination_plot(checkm, output_path):
    df = checkm.loc[:,['Completeness', 'Contamination']]

    df.rename(columns={'Completeness': 'Completeness (%)', 'Contamination': 'Contamination (%)'}, inplace=True)

    conditions = [(df['Completeness (%)'] >= 90.0) & (df['Contamination (%)'] <= 5.0), (df['Completeness (%)'] >= 90.0) & (df['Contamination (%)'] > 5.0)]

    choices = ['red', 'blue']

    colors = np.select(conditions, choices, default='grey')

    plt.scatter(x=df['Completeness (%)'], y=df['Contamination (%)'], c=colors, edgecolor='black')

    plt.xlim(max(math.floor((df['Completeness (%)'].min() - 10) / 10) * 10 ,0), 100)
    plt.ylim(0, min(math.ceil(df['Contamination (%)'].max() + 1), 100))

    plt.xlabel('Completeness (%)')
    plt.ylabel('Contamination (%)')

    legend_patches = [
        mpatches.Patch(color='red', label='High Quality (≥90% comp, ≤5% cont)'),
        mpatches.Patch(color='blue', label='Contaminated (≥90% comp, >5% cont)'),
        mpatches.Patch(color='grey', label='Others (<90% comp)')
    ]
    plt.legend(handles=legend_patches, loc='best')

    plt.savefig(os.path.join(output_path,"comp_conta.png"))

def rank_completeness_contamination_plot(checkm, gtdb_bac, gtdb_ar, rank, output_path, n):

    checkm.index = checkm.index.str.replace('.', '_', regex=False)

    gtdb_df = pd.concat([gtdb_ar, gtdb_bac], ignore_index=False)
    
    gtdb_df[f'{rank.capitalize()}'] = gtdb_df['classification'].apply(lambda x: extract_rank(x, rank))

    # ---- Join GTDB phylum into CheckM based on index ----
    merged_df = checkm.join(gtdb_df[[f'{rank.capitalize()}']])

    # ---- Create a "Phylum (n)" column for labeling ----
    rank_counts = merged_df[f'{rank.capitalize()}'].value_counts()
    top_ranks = rank_counts.nlargest(n).index.tolist()

    other_count = rank_counts[~rank_counts.index.isin(top_ranks)].sum()
    
    def map_rank(r):
        if pd.isnull(r):
            return f"Unknown {rank.capitalize()}"
        elif r in top_ranks:
            return f"{r} ({rank_counts[r]})"
        else:
            return f"Other ({other_count})"

    merged_df[f'{rank.capitalize()} (n)'] = merged_df[f'{rank.capitalize()}'].map(map_rank)

    # ---- Plot: Completeness vs Contamination with counts in legend ----
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 6))
    sns.scatterplot(
        data=merged_df,
        y='Contamination',
        x='Completeness',
        hue=f'{rank.capitalize()} (n)',
        palette='tab20',
        s=80,
        alpha=0.8
    )

    plt.title(f'CheckM: Completeness vs Contamination (Colored by {rank.capitalize()})', fontsize=14, weight='bold')
    plt.ylabel('Contamination (%)')
    plt.xlabel('Completeness (%)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title=f'{rank.capitalize()} (n)')
    plt.tight_layout()

    plt.savefig(os.path.join(output_path, "comp_conta_by_rank.png"))
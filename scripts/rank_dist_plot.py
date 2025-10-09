import matplotlib.pyplot as plt
import os

prefix_map = {
    'domain': ('d', 0),
    'phylum': ('p', 1),
    'class': ('c', 2),
    'order': ('o', 3),
    'family': ('f', 4),
    'genus': ('g', 5),
    'species': ('s', 6),
}

def rank_distribution_pie(gtdb, output_path, rank, n):
    df = gtdb.reset_index().loc[:,['user_genome', 'classification']]

    df.rename(columns={'user_genome': 'Bin', 'classification': f'{rank.capitalize()}'}, inplace=True)

    taxa_level, pos = prefix_map.get(rank, ("taxa",-1))
    df[f'{rank.capitalize()}'] = (
        df[f'{rank.capitalize()}']
        .str.split(';')
        .str[pos]
        .str.replace(f'{taxa_level}__', '', regex=False)
    )

    df[f'{rank.capitalize()}'] = df[f'{rank.capitalize()}'].replace('', f'Unknow {rank.capitalize()}')

    rank_counts = df[f'{rank.capitalize()}'].value_counts()

    top_counts = rank_counts.head(n)

    plt.figure(figsize=(8,8))
    top_counts.plot(kind="pie", autopct='%1.1f%%')
    plt.ylabel("")
    plt.title(f"{rank.capitalize()}-level distribution of MAGs")

    plt.savefig(os.path.join(output_path, "rank_dist_pie.png"))
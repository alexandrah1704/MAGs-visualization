import os
import matplotlib.pyplot as plt

def create_n50_histogram(checkm2, output_path):
    df = checkm2.loc[:,['Contig_N50']]

    df.rename(columns={'Contig_N50': 'N50(kbp)'}, inplace=True)

    df = df / 1000

    plt.figure()
    plt.hist(df, bins='auto', color='skyblue', edgecolor='black')

    plt.xlim(0, df['N50(kbp)'].max())

    plt.xlabel('N50 (kbp)')
    plt.ylabel('Frequency')

    plt.xticks(plt.xticks()[0], [f'{int(x)} kb' for x in plt.xticks()[0]])

    plt.savefig(os.path.join(output_path, "n50_histogram.png"))

def number_of_contigs(checkm2, output_path):
    df = checkm2.loc[:,['Total_Contigs']]

    df.rename(columns={'Total_Contigs': '# Contig'}, inplace=True)
    
    plt.figure()
    plt.hist(df, bins='auto', color='skyblue', edgecolor='black')

    plt.xlim(0, df['# Contig'].max())

    plt.xlabel('Number of contigs per genome')
    plt.ylabel('Frequency')

    plt.savefig(os.path.join(output_path, "number_of_contig_his.png"))

def create_assambly_info_histo(checkm2, output_path):
    df = checkm2.reset_index().loc[:,['Name', 'Contig_N50', 'Genome_Size', 'Max_Contig_Length', 'Coding_Density']]

    df.rename(columns={'Name': 'Bin', 'Genome_Size': 'Total length Assembly', 'Contig_N50' : 'N50', 'Max_Contig_Length': 'Longest Contig', 'Coding_Density': 'Coding Density'}, inplace=True)

    metrics = ['N50','Total length Assembly','Longest Contig','Coding Density']

    fig, axes = plt.subplots(2, 2, figsize=(12,8))
    axes = axes.flatten()

    for i, metric in enumerate(metrics):
        axes[i].hist(df[metric], bins='auto', color='skyblue', edgecolor='black')
        axes[i].set_title(f'Distribution of {metric}')
        axes[i].set_xlabel(metric)
        axes[i].set_ylabel('Count')

    plt.tight_layout()

    plt.savefig(os.path.join(output_path, "assambly_info_histo.png"))
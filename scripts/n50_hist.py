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

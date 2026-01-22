import matplotlib.pyplot as plt
import seaborn as sns
import os

def mag_detection_heatmap(coverm, output_path):

    num_bins, num_samples = coverm.shape

    plt.figure(figsize=(max(10, num_samples*1.2), max(6, num_bins*0.3)))

    sns.heatmap(coverm, cmap= sns.diverging_palette(240, 10, as_cmap=True), annot=True, cbar_kws={'label': 'Abundance'}, linewidths=0.5,linecolor='gray')

    plt.xticks(rotation=45, ha="right", fontsize=10)
    plt.yticks(rotation=0, fontsize=8)
    plt.title("MAG Detection Across Samples")
    plt.xlabel("Samples")
    plt.ylabel("MAGs")
    plt.tight_layout()
    
    plt.savefig(os.path.join(output_path, "mag_detection_heatmap.png"))
import matplotlib.pyplot as plt
import os
import seaborn as sns

def binner_plot(amber, output_path):

    df = amber[amber['binning type'] == 'genome']

    sns.set_theme(style="whitegrid")

    # Create barplot
    plt.figure(figsize=(12, 6))
    ax = sns.barplot(
        data=df,
        x='Tool',
        y='f1_score_per_bp',
        hue='Tool',
        legend=False,
        palette='crest'
    )

    # Customize plot
    ax.set_title('F1 Score per Base Pair (Genome Binning)', fontsize=14, weight='bold')
    ax.set_xlabel('Binner', fontsize=12)
    ax.set_ylabel('F1 Score per bp', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    plt.savefig(os.path.join(output_path, "binner_compare.png"))
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
import numpy as np
import os

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
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

def rarefaction_curve(df, n_iter=100, step=1):
    genomes = df["Bin"].tolist()
    species_map = dict(zip(df["Bin"], df["Cluster"]))
    
    max_depth = len(genomes)
    depths = range(step, max_depth+1, step)
    results = []
    
    for depth in depths:
        species_counts = []
        for _ in range(n_iter):
            subsample = np.random.choice(genomes, size=depth, replace=False)
            clusters = set(species_map[g] for g in subsample)
            species_counts.append(len(clusters))
        results.append({
            "depth": depth,
            "mean_species": np.mean(species_counts),
            "std_species": np.std(species_counts)
        })
    
    return pd.DataFrame(results)

def create_n50_plot(drep, output_path):
    df = drep.reset_index().loc[:,['genome', 'secondary_cluster']]

    df.rename(columns={'genome': 'Bin', 'secondary_cluster': 'Cluster'}, inplace=True)

    results = rarefaction_curve(df, n_iter=200, step=1)


    plt.figure(figsize=(7,5))
    plt.plot(results["depth"], results["mean_species"], label="Mean species richness", color="black")
    plt.fill_between(results["depth"],
                    results["mean_species"]-results["std_species"],
                    results["mean_species"]+results["std_species"],
                    color="gray", alpha=0.3, label="±1 SD")
    plt.xlabel("Number of genomes sampled")
    plt.ylabel("Number of species clusters (>95% ANI)")
    plt.title("Species-level Rarefaction Curve")
    plt.legend()

    plt.savefig(os.path.join(output_path, "n50_plot.png"))
import pandas as pd
import plotly.graph_objects as go
import os
import re

prefix_map = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species"
}

def generate_taxa_sanky(gtdb, output_path):
    tax_split = gtdb["classification"].str.split(";", expand=True)
    tax_split.columns = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    tax_split = tax_split.replace({"": None, " ": None})

    links = []
    for i in range(len(tax_split.columns) - 1):
        pairs = (
            tax_split[[tax_split.columns[i], tax_split.columns[i+1]]]
            .dropna()
            .value_counts()
            .reset_index(name="count")
        )
        pairs.columns = ["source", "target", "count"]
        links.append(pairs)

    links_df = pd.concat(links, ignore_index=True)

    nodes = pd.Index(pd.concat([links_df["source"], links_df["target"]]).unique())
    node_map = {name: i for i, name in enumerate(nodes)}

    links_df["source_idx"] = links_df["source"].map(node_map)
    links_df["target_idx"] = links_df["target"].map(node_map)

    rank_colors = {
        "domain": "#1f77b4",   # blue
        "phylum": "#ff7f0e",   # orange
        "class": "#2ca02c",    # green
        "order": "#d62728",    # red
        "family": "#9467bd",   # purple
        "genus": "#8c564b",    # brown
        "species": "#e377c2"   # pink
    }

    node_colors = []
    for node in nodes:
        assigned = False
        for rank in tax_split.columns:
            if node in tax_split[rank].values:
                node_colors.append(rank_colors[rank])
                assigned = True
                break
        if not assigned:
            node_colors.append("lightgray")

    clean_labels = []
    for label in nodes:
        if not label:
            clean_labels.append("")
        else:
            match = re.match(r"^([a-zA-Z]+)_", label)
            prefix = match.group(1)
            taxa_level = prefix_map.get(prefix, "taxa")
            cleaned = re.sub(r"^[a-z]__*", "", label)
            clean_labels.append(cleaned if cleaned else f"Unknown {taxa_level.capitalize()}")


    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=clean_labels,
            color=node_colors
        ),
        link=dict(
            source=links_df["source_idx"],
            target=links_df["target_idx"],
            value=links_df["count"]
        )
    )])

    fig.update_layout(
        title_text="Taxonomic Classification Sankey",
        font=dict(size=10),
        width=1600,
        height=900
    )

    fig.write_html(os.path.join(output_path,"sankey_plot.html"))
    
    
    #fig.write_image(os.path.join(output_path,"sankey_plot.png")) --> Possible but there are a lot of libraries needed to make this work so if this is wanted i can add them all as requirements

def taxa_sanky_rank(gtdb, output_path, rank):

    tax_split = gtdb.reset_index()["classification"].str.split(";", expand=True)
    tax_split.columns = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    tax_split = tax_split.replace({"": None, " ": None})

    df_sankey = gtdb.reset_index()[["user_genome"]].copy()
    df_sankey[rank] = tax_split[rank]

    links_df = df_sankey.groupby(["user_genome", rank]).size().reset_index(name="count")
    links_df.columns = ["source", "target", "count"]

    nodes = pd.Index(pd.concat([links_df["source"], links_df["target"]]).unique())
    node_map = {name: i for i, name in enumerate(nodes)}
    links_df["source_idx"] = links_df["source"].map(node_map)
    links_df["target_idx"] = links_df["target"].map(node_map)

    rank_colors = {
        "genome": "#7f7f7f",  
        "domain": "#1f77b4",
        "phylum": "#ff7f0e",
        "class": "#2ca02c",
        "order": "#d62728",
        "family": "#9467bd",
        "genus": "#8c564b",
        "species": "#e377c2"
    }

    node_colors = []
    for node in nodes:
        if node in df_sankey["user_genome"].values:
            node_colors.append(rank_colors["genome"])
        else:
            node_colors.append(rank_colors.get(rank, "lightgray"))

    clean_labels = []
    for label in nodes:
        if label.split('_')[-1] == 'fasta':
            clean_labels.append(label)
        elif label.split('__')[-1] == '':
            clean_labels.append(f'Unknow {rank.capitalize()}')
        else:
            clean_labels.append(label.split('__')[-1])

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=15,
            line=dict(color="black", width=0.5),
            label=clean_labels,
            color=node_colors
        ),
        link=dict(
            source=links_df["source_idx"],
            target=links_df["target_idx"],
            value=links_df["count"]
        )
    )])

    fig.update_layout(
        title_text=f"Genome → {rank.capitalize()} Sankey",
        font=dict(size=10),
        width=1200,
        height=800
    )

    fig.write_html(os.path.join(output_path,"sankey_plot_rank_filtered.html"))
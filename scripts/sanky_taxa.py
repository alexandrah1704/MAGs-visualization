import pandas as pd
import plotly.graph_objects as go
import os

def generate_taxa_sanky(gtdb, output_path):
    tax_split = gtdb["classification"].str.split(";", expand=True)
    tax_split.columns = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    # Clean up missing/empty values
    tax_split = tax_split.replace({"": None, " ": None})

    # --- Build source-target pairs for Sankey ---
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

    # --- Build mapping of node names to indices ---
    nodes = pd.Index(pd.concat([links_df["source"], links_df["target"]]).unique())
    node_map = {name: i for i, name in enumerate(nodes)}

    # Map source/target to indices
    links_df["source_idx"] = links_df["source"].map(node_map)
    links_df["target_idx"] = links_df["target"].map(node_map)

    # --- Assign colors based on taxonomy rank ---
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
            node_colors.append("lightgray")  # fallback color

    # --- Make labels nicer (add line breaks to avoid overlap) ---
    clean_labels = [label.replace(";", ";\n") if label else "" for label in nodes]

    # --- Create Sankey diagram ---
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

    # --- Layout adjustments ---
    fig.update_layout(
        title_text="Taxonomic Classification Sankey",
        font=dict(size=10),
        width=1600,   # wider to reduce overlap
        height=900
    )

    fig.write_html(os.path.join(output_path,"sankey_plot.html"))
    
    
    #fig.write_image(os.path.join(output_path,"sankey_plot.png")) --> Possible but there are a lot of libraries needed to make this work so if this is wanted i can add them all as requirements
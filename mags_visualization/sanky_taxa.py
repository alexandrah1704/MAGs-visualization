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

def normalize_rank(rank: str) -> str:
    """Accepts short codes (d/p/c/o/f/g/s) or full names; returns full rank name."""
    if not isinstance(rank, str):
        raise ValueError("rank must be a string")
    r = rank.strip().lower()
    # map short code -> full
    short2full = {
        "d": "domain", "p": "phylum", "c": "class",
        "o": "order",  "f": "family", "g": "genus", "s": "species"
    }
    return short2full.get(r, r)  # if already 'phylum', returns 'phylum'


def clean_tax_label(label: str, rank_name: str) -> str:
    """
    Convert gtdb labels ("g_Gilliamella") into display labels.
    Missing values = unkown.
    """
    if pd.isna(label) or not str(label).strip():
        return f"Unknown {rank_name.capitalize()}"
    
    label = str(label).strip()
    cleaned = re.sub(r"^[a-z]__*", "", label).strip()
    
    if cleaned:
        return cleaned
    return f"Unknown {rank_name.capitalize()}"


def unique_tax_id(row, rank_idx, tax_columns):
    """"
    Create internal node ID.
    Known taxa grouped by rank + taxon value.
    Missing taxa grouped by their parent taxon.
    """
    parts = []
    for i in range(rank_idx + 1):
        rank = tax_columns[i]
        value = row[rank]

        if pd.isna(value) or not str(value).strip():
            value = f"unknown_{rank}"

        parts.append(f"{rank}={value}")

    return "||".join(parts)


def is_missing_taxon(value) -> bool:
    """
    Return True if a taxonomy entry is missing or unclassified
    """
    if pd.isna(value):
        return True

    value = str(value).strip()
    if not value:
        return True

    # Remove GTDB prefix like d__, p__, g__
    cleaned = re.sub(r"^[a-z]__*", "", value).strip()

    if cleaned == "":
        return True

    cleaned_lower = cleaned.lower()

    # Treat generic unknown/unclassified labels as missing
    if cleaned_lower.startswith("unclassified"):
        return True
    if cleaned_lower.startswith("unknown"):
        return True
    if cleaned_lower in {"uncultured", "unidentified", "na", "n/a", "none"}:
        return True

    return False

def generate_taxa_sanky(gtdb, output_path, rank):
    if rank is not None:
        rank = normalize_rank(rank)

    tax_columns = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    tax_split = gtdb.reset_index()["classification"].str.split(";", expand=True)
    tax_split.columns = tax_columns[:tax_split.shape[1]]
    tax_split = tax_split.reindex(columns=tax_columns)
    tax_split = tax_split.replace({"": None, " ": None})

    if rank is not None and rank not in tax_columns:
        raise ValueError(f"Requested rank '{rank}' is not in parsed columns {tax_columns}")

    links = []

    for _, row in tax_split.iterrows():
        valid_ranks = []

        # collect only real classified ranks
        for current_rank in tax_columns:
            current_value = row[current_rank]
            if not is_missing_taxon(current_value):
                valid_ranks.append((current_rank, str(current_value).strip()))

        # skip completely unclassified rows
        if not valid_ranks:
            continue

        # build links between known ranks
        previous_node = None
        for current_rank, current_value in valid_ranks:
            current_node = f"{current_rank}|{current_value}"

            if previous_node is not None:
                links.append((previous_node, current_node))

            previous_node = current_node

        # if species is missing, add one terminal Unknown Species node
        species_value = row["species"]
        if is_missing_taxon(species_value):
            last_rank, last_value = valid_ranks[-1]

            if last_rank == "genus":
                unknown_species_node = f"species|unknown_from|{last_rank}|{last_value}"
                links.append((previous_node, unknown_species_node))

    if not links:
        raise ValueError("No classified taxonomy paths were found for Sankey plotting.")

    links_df = (
        pd.DataFrame(links, columns=["source", "target"])
        .value_counts()
        .reset_index(name="count")
    )

    nodes = pd.Index(pd.concat([links_df["source"], links_df["target"]]).unique())
    node_map = {name: i for i, name in enumerate(nodes)}

    links_df["source_idx"] = links_df["source"].map(node_map)
    links_df["target_idx"] = links_df["target"].map(node_map)

    rank_colors = {
        "domain": "#1f77b4",
        "phylum": "#ff7f0e",
        "class": "#2ca02c",
        "order": "#d62728",
        "family": "#9467bd",
        "genus": "#8c564b",
        "species": "#e377c2"
    }

    clean_labels = []
    node_colors = []

    for node in nodes:
        parts = str(node).split("|")
        node_rank = parts[0]
        node_colors.append(rank_colors.get(node_rank, "lightgray"))

        # normal node: "rank|value"
        if len(parts) == 2:
            raw_value = parts[1]
            clean_labels.append(clean_tax_label(raw_value, node_rank))

        # unknown species: "species|unknown_from|genus|g__Gilliamella"
        elif len(parts) == 4 and parts[1] == "unknown_from":
            parent_rank = parts[2]
            parent_value = parts[3]
            parent_label = clean_tax_label(parent_value, parent_rank)
            clean_labels.append(f"Unknown Species")

        else:
            clean_labels.append(node)

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

    os.makedirs(output_path, exist_ok=True)
    fig.write_html(os.path.join(output_path, "sankey_plot.html"))
    
    
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
        label_str = str(label)

        if label_str.endswith(".fasta"):
            clean_labels.append(label_str)
        elif "__" in label_str:
            val = label_str.split("__")[-1]
            clean_labels.append(val if val else f"Unknown {rank.capitalize()}")
        else:
            clean_labels.append(label_str)

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
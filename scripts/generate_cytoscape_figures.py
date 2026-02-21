#!/usr/bin/env python3
"""Generate publication-ready network figures using Cytoscape via py4cytoscape.

Produces 3 figures from existing KG data in research-results/knowledge_graph/:
  1. Hub Gene Subnetwork — 20 hubs + top neighbors, colored by community
  2. Drug-Hub Target Map — hub genes + FDA-approved/clinical-stage drugs
  3. Community Overview — full 336-node network with community-aware layout

Prerequisites:
  - Cytoscape desktop running (CyREST on localhost:1234)
  - py4cytoscape installed in pathwayenv
  - All data files present in research-results/knowledge_graph/

Usage:
  source pathwayenv/bin/activate
  python scripts/generate_cytoscape_figures.py
"""

import json
import math
import sys
import time
from pathlib import Path

import networkx as nx
import pandas as pd
import py4cytoscape as p4c

# ---------------------------------------------------------------------------
# Phase A: Setup & Data Loading
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path("research-results/knowledge_graph")

# Color palettes
COMMUNITY_COLORS = {
    1: "#E74C3C",  # Red - Dopamine/Wnt
    2: "#3498DB",  # Blue - Synaptic/Glutamate/GABA
    3: "#2ECC71",  # Green - Chromatin/Wnt/mTOR
    4: "#F39C12",  # Orange - Calcium/Ion channels
    5: "#9B59B6",  # Purple - Circadian/GABA
    6: "#1ABC9C",  # Teal - Neurodevelopment/Axon
    7: "#E67E22",  # Dark orange - Immune/Oxidative
}
COMMUNITY_COLORS_DEFAULT = "#CCCCCC"  # For singleton communities (0, 8-13)

DISEASE_COLORS = {
    "ASD": "#4A90D9",
    "SCZ": "#D94A4A",
    "Both": "#9B59B6",
}
DISEASE_COLORS_DEFAULT = "#999999"


def verify_cytoscape():
    """Verify Cytoscape is running and accessible."""
    try:
        p4c.cytoscape_ping()
        info = p4c.cytoscape_version_info()
        print(f"Connected to Cytoscape {info['cytoscapeVersion']}")
        return True
    except Exception as e:
        print(f"ERROR: Cannot connect to Cytoscape: {e}")
        print("Make sure Cytoscape desktop is running.")
        sys.exit(1)


def load_data():
    """Load all input data files."""
    print("\n--- Loading data ---")
    data = {}

    # Full PPI network (GraphML)
    data["G_full"] = nx.read_graphml(OUTPUT_DIR / "cross_disease_network.graphml")
    print(
        f"Full network: {data['G_full'].number_of_nodes()} nodes, "
        f"{data['G_full'].number_of_edges()} edges"
    )

    # Hub genes
    data["hub_df"] = pd.read_csv(OUTPUT_DIR / "hub_genes_ranked.csv")
    data["hub_genes"] = set(data["hub_df"]["gene"].tolist())
    print(f"Hub genes: {len(data['hub_genes'])}")

    # Community analysis
    data["community_df"] = pd.read_csv(OUTPUT_DIR / "community_analysis.csv")

    # PPI edges with proper scores
    data["ppi_df"] = pd.read_csv(OUTPUT_DIR / "string_ppi_edges.csv")
    print(f"PPI edges: {len(data['ppi_df'])}")

    # Drug interactions
    data["drug_df"] = pd.read_csv(OUTPUT_DIR / "dgidb_interactions.csv")
    print(f"Drug interactions: {len(data['drug_df'])}")

    # Drug repurposing table
    data["drug_repurp_df"] = pd.read_csv(OUTPUT_DIR / "drug_repurposing_table.csv")

    return data


def clear_cytoscape_session():
    """Delete all existing networks in Cytoscape."""
    try:
        p4c.close_session(save_before_closing=False)
        time.sleep(1)
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Phase B: Figure 1 — Hub Gene Subnetwork
# ---------------------------------------------------------------------------


def build_hub_subnetwork(data):
    """Build a subnetwork of hub genes + their top PPI neighbors."""
    print("\n--- Figure 1: Building hub subnetwork ---")

    G_full = data["G_full"]
    hub_genes = data["hub_genes"]
    ppi_df = data["ppi_df"]

    # Build a PPI graph with proper scores from the CSV
    G_ppi = nx.Graph()
    for _, row in ppi_df.iterrows():
        g1, g2, score = row["gene1"], row["gene2"], row["combined_score"]
        if g1 in G_full and g2 in G_full:
            G_ppi.add_edge(g1, g2, combined_score=score)

    # Get hub gene neighbors with high-confidence edges (score > 0.6)
    neighbor_degrees = {}
    for hub in hub_genes:
        if hub not in G_ppi:
            continue
        for neighbor in G_ppi.neighbors(hub):
            if neighbor in hub_genes:
                continue
            score = G_ppi[hub][neighbor]["combined_score"]
            if score > 0.6:
                if neighbor not in neighbor_degrees:
                    neighbor_degrees[neighbor] = 0
                neighbor_degrees[neighbor] += 1

    # Take top 30 neighbors by number of hub connections
    top_neighbors = sorted(neighbor_degrees, key=neighbor_degrees.get, reverse=True)[
        :30
    ]
    selected_nodes = hub_genes | set(top_neighbors)

    # Build subgraph
    G_hub = nx.Graph()
    for node in selected_nodes:
        if node in G_full:
            attrs = dict(G_full.nodes[node])
            comm = int(attrs.get("community", 0))
            disease = attrs.get("disease", "Both")
            is_hub = node in hub_genes
            G_hub.add_node(
                node,
                community=str(comm),
                community_color=COMMUNITY_COLORS.get(comm, COMMUNITY_COLORS_DEFAULT),
                disease=disease,
                disease_color=DISEASE_COLORS.get(disease, DISEASE_COLORS_DEFAULT),
                is_hub=is_hub,
                degree=G_ppi.degree(node) if node in G_ppi else 0,
            )

    # Add edges between selected nodes using PPI scores
    for u, v, d in G_ppi.edges(data=True):
        if u in selected_nodes and v in selected_nodes:
            G_hub.add_edge(u, v, combined_score=d["combined_score"])

    print(f"Hub subnetwork: {G_hub.number_of_nodes()} nodes, {G_hub.number_of_edges()} edges")
    return G_hub


def render_hub_subnetwork(G_hub, data):
    """Load hub subnetwork into Cytoscape, style, and export."""
    print("Loading into Cytoscape...")

    hub_genes = data["hub_genes"]
    hub_df = data["hub_df"]
    style_name = "HubSubnetwork"

    # Create network in Cytoscape
    suid = p4c.create_network_from_networkx(G_hub, title="Hub Gene Subnetwork")
    time.sleep(2)

    # Create visual style
    defaults = {
        "NODE_SHAPE": "ELLIPSE",
        "NODE_SIZE": 25,
        "NODE_FILL_COLOR": "#CCCCCC",
        "NODE_BORDER_WIDTH": 1,
        "NODE_BORDER_PAINT": "#666666",
        "NODE_LABEL_FONT_SIZE": 10,
        "NODE_LABEL_COLOR": "#000000",
        "NODE_TRANSPARENCY": 230,
        "EDGE_WIDTH": 1.0,
        "EDGE_TRANSPARENCY": 60,
        "EDGE_STROKE_UNSELECTED_PAINT": "#AAAAAA",
        "NETWORK_BACKGROUND_PAINT": "#FFFFFF",
    }
    p4c.create_visual_style(style_name, defaults=defaults)
    p4c.set_visual_style(style_name)
    time.sleep(1)

    # Map node fill color by community
    comm_values = sorted(set(nx.get_node_attributes(G_hub, "community").values()))
    comm_colors = [COMMUNITY_COLORS.get(int(c), COMMUNITY_COLORS_DEFAULT) for c in comm_values]
    p4c.set_node_color_mapping(
        "community",
        table_column_values=comm_values,
        colors=comm_colors,
        mapping_type="d",
        style_name=style_name,
    )

    # Map node border color by disease
    disease_vals = ["ASD", "SCZ", "Both"]
    disease_cols = [DISEASE_COLORS[d] for d in disease_vals]
    p4c.set_node_border_color_mapping(
        "disease",
        table_column_values=disease_vals,
        colors=disease_cols,
        mapping_type="d",
        default_color="#666666",
        style_name=style_name,
    )

    # Map node border width: hub=4, non-hub=1
    p4c.set_node_border_width_mapping(
        "is_hub",
        table_column_values=["True", "False"],
        widths=[4, 1],
        mapping_type="d",
        style_name=style_name,
    )

    # Map node size by degree (continuous mapping)
    degrees = nx.get_node_attributes(G_hub, "degree")
    min_deg = min(degrees.values()) if degrees else 1
    max_deg = max(degrees.values()) if degrees else 10
    p4c.set_node_size_mapping(
        "degree",
        table_column_values=[min_deg, max_deg],
        sizes=[20, 80],
        mapping_type="c",
        style_name=style_name,
    )

    # Set labels only for hub genes using bypass
    all_nodes = list(G_hub.nodes())
    hub_node_list = [n for n in all_nodes if n in hub_genes]
    non_hub_list = [n for n in all_nodes if n not in hub_genes]

    # Set all labels first via mapping, then clear non-hub labels
    p4c.set_node_label_mapping("name", style_name=style_name)
    time.sleep(0.5)

    # Hide non-hub labels by setting them to empty string
    if non_hub_list:
        p4c.set_node_label_bypass(non_hub_list, [""] * len(non_hub_list))

    # Make hub gene labels bold and larger
    if hub_node_list:
        p4c.set_node_font_size_bypass(hub_node_list, [14] * len(hub_node_list))

    # Map edge width by combined_score
    p4c.set_edge_line_width_mapping(
        "combined_score",
        table_column_values=[0.4, 1.0],
        widths=[0.5, 3.0],
        mapping_type="c",
        style_name=style_name,
    )

    # Layout
    print("Applying force-directed layout...")
    p4c.layout_network("force-directed")
    time.sleep(3)

    # Export (Cytoscape 3.10+ API: use zoom instead of deprecated resolution/units/height/width)
    p4c.fit_content()
    time.sleep(1)
    out_path = str(OUTPUT_DIR / "hub_subnetwork_figure.png")
    print(f"Exporting to {out_path}")
    p4c.export_image(
        out_path,
        type="PNG",
        zoom=400,
        all_graphics_details=True,
        overwrite_file=True,
    )
    time.sleep(1)
    print(f"Figure 1 exported: {out_path}")
    return suid


# ---------------------------------------------------------------------------
# Phase C: Figure 2 — Drug-Hub Target Map
# ---------------------------------------------------------------------------


def build_drug_hub_network(data):
    """Build bipartite network of hub genes and their drug targets."""
    print("\n--- Figure 2: Building drug-hub target map ---")

    hub_genes = data["hub_genes"]
    hub_df = data["hub_df"]
    drug_df = data["drug_df"]
    G_full = data["G_full"]

    # Filter drugs for hub genes
    hub_drugs = drug_df[drug_df["gene"].isin(hub_genes)].copy()

    # Curated key drugs (case-insensitive match)
    curated_drugs = {
        "EVEROLIMUS", "LITHIUM", "MEMANTINE HYDROCHLORIDE",
        "IPATASERTIB", "CAPIVASERTIB", "GANAXOLONE", "INFLIXIMAB",
        "TIDEGLUSIB", "SIROLIMUS", "TEMSIROLIMUS",
    }

    # Select drugs: approved + high-evidence + curated
    selected_drugs = set()

    # All approved drugs for hub genes
    approved = hub_drugs[hub_drugs["approved"] == True]  # noqa: E712
    for _, row in approved.iterrows():
        selected_drugs.add(row["drug"])

    # High-evidence drugs (n_sources >= 3)
    high_evidence = hub_drugs[hub_drugs["n_sources"] >= 3]
    for _, row in high_evidence.iterrows():
        selected_drugs.add(row["drug"])

    # Add curated drugs
    for _, row in hub_drugs.iterrows():
        if row["drug"].upper() in curated_drugs or row["drug"] in curated_drugs:
            selected_drugs.add(row["drug"])

    # Prioritize by n_sources and cap at 25
    drug_priority = (
        hub_drugs[hub_drugs["drug"].isin(selected_drugs)]
        .groupby("drug")
        .agg({"n_sources": "max", "approved": "first", "interaction_type": "first"})
        .sort_values("n_sources", ascending=False)
    )
    top_drugs = set(drug_priority.head(25).index)

    # Build bipartite graph
    G_drug = nx.Graph()

    # Add hub gene nodes
    for gene in hub_genes:
        attrs = dict(G_full.nodes[gene]) if gene in G_full else {}
        comm = int(attrs.get("community", 0))
        disease = attrs.get("disease", "Both")
        ppi_degree = G_full.degree(gene) if gene in G_full else 0
        G_drug.add_node(
            gene,
            node_type="gene",
            community=str(comm),
            community_color=COMMUNITY_COLORS.get(comm, COMMUNITY_COLORS_DEFAULT),
            disease=disease,
            disease_color=DISEASE_COLORS.get(disease, DISEASE_COLORS_DEFAULT),
            ppi_degree=ppi_degree,
        )

    # Add drug nodes and edges
    drug_interactions = hub_drugs[hub_drugs["drug"].isin(top_drugs)]
    for _, row in drug_interactions.iterrows():
        drug_name = row["drug"]
        gene = row["gene"]
        approved_status = row["approved"]
        interaction = row["interaction_type"] if pd.notna(row["interaction_type"]) else "unknown"

        if drug_name not in G_drug:
            G_drug.add_node(
                drug_name,
                node_type="drug",
                approved=str(approved_status),
                community="drug",
                disease="drug",
            )

        # Categorize interaction type for coloring
        if interaction in ("inhibitor", "blocker", "negative modulator"):
            edge_class = "inhibitor"
        elif interaction in ("agonist", "activator", "positive modulator"):
            edge_class = "agonist"
        elif interaction in ("modulator", "antibody", "binder"):
            edge_class = "modulator"
        else:
            edge_class = "unknown"

        G_drug.add_edge(
            drug_name,
            gene,
            interaction_type=interaction,
            edge_class=edge_class,
        )

    n_genes = sum(1 for _, d in G_drug.nodes(data=True) if d.get("node_type") == "gene")
    n_drugs = sum(1 for _, d in G_drug.nodes(data=True) if d.get("node_type") == "drug")
    print(f"Drug-hub network: {n_genes} genes + {n_drugs} drugs = {G_drug.number_of_nodes()} nodes, {G_drug.number_of_edges()} edges")
    return G_drug


def render_drug_hub_network(G_drug, data):
    """Load drug-hub network into Cytoscape, style, and export."""
    print("Loading into Cytoscape...")

    hub_genes = data["hub_genes"]
    style_name = "DrugHubTarget"

    suid = p4c.create_network_from_networkx(G_drug, title="Drug-Hub Target Map")
    time.sleep(2)

    # Create visual style
    defaults = {
        "NODE_SHAPE": "ELLIPSE",
        "NODE_SIZE": 35,
        "NODE_FILL_COLOR": "#CCCCCC",
        "NODE_BORDER_WIDTH": 2,
        "NODE_BORDER_PAINT": "#333333",
        "NODE_LABEL_FONT_SIZE": 11,
        "NODE_LABEL_COLOR": "#000000",
        "NODE_TRANSPARENCY": 240,
        "EDGE_WIDTH": 2.0,
        "EDGE_TRANSPARENCY": 150,
        "EDGE_STROKE_UNSELECTED_PAINT": "#AAAAAA",
        "NETWORK_BACKGROUND_PAINT": "#FFFFFF",
    }
    p4c.create_visual_style(style_name, defaults=defaults)
    p4c.set_visual_style(style_name)
    time.sleep(1)

    # Set labels
    p4c.set_node_label_mapping("name", style_name=style_name)

    # Node shape: gene=ELLIPSE, drug=DIAMOND
    p4c.set_node_shape_mapping(
        "node_type",
        table_column_values=["gene", "drug"],
        shapes=["ELLIPSE", "DIAMOND"],
        style_name=style_name,
    )

    # Gene nodes: color by community
    gene_nodes = [n for n, d in G_drug.nodes(data=True) if d.get("node_type") == "gene"]
    drug_nodes = [n for n, d in G_drug.nodes(data=True) if d.get("node_type") == "drug"]

    # Color gene nodes by community using bypass
    for gene in gene_nodes:
        comm = int(G_drug.nodes[gene].get("community", 0))
        color = COMMUNITY_COLORS.get(comm, COMMUNITY_COLORS_DEFAULT)
        p4c.set_node_color_bypass([gene], [color])

    # Color drug nodes by approval status
    approved_drugs = [n for n in drug_nodes if G_drug.nodes[n].get("approved") == "True"]
    unapproved_drugs = [n for n in drug_nodes if G_drug.nodes[n].get("approved") != "True"]
    if approved_drugs:
        p4c.set_node_color_bypass(approved_drugs, ["#27AE60"] * len(approved_drugs))
    if unapproved_drugs:
        p4c.set_node_color_bypass(unapproved_drugs, ["#F4D03F"] * len(unapproved_drugs))

    # Size gene nodes by PPI degree, drug nodes uniform
    for gene in gene_nodes:
        deg = G_drug.nodes[gene].get("ppi_degree", 10)
        size = max(30, min(70, 30 + deg * 0.4))
        p4c.set_node_size_bypass([gene], [size])
    if drug_nodes:
        p4c.set_node_size_bypass(drug_nodes, [30] * len(drug_nodes))

    # Gene node border by disease
    for gene in gene_nodes:
        disease = G_drug.nodes[gene].get("disease", "Both")
        border_color = DISEASE_COLORS.get(disease, "#666666")
        p4c.set_node_border_color_bypass([gene], [border_color])
    if drug_nodes:
        p4c.set_node_border_color_bypass(drug_nodes, ["#1B7A3D"] * len(drug_nodes))

    # Edge color by interaction class
    edge_class_colors = {
        "inhibitor": "#E74C3C",
        "agonist": "#3498DB",
        "modulator": "#F39C12",
        "unknown": "#95A5A6",
    }
    p4c.set_edge_color_mapping(
        "edge_class",
        table_column_values=list(edge_class_colors.keys()),
        colors=list(edge_class_colors.values()),
        mapping_type="d",
        style_name=style_name,
    )

    # Gene label font larger
    if gene_nodes:
        p4c.set_node_font_size_bypass(gene_nodes, [13] * len(gene_nodes))
    if drug_nodes:
        p4c.set_node_font_size_bypass(drug_nodes, [9] * len(drug_nodes))

    # Layout: hub genes in a ring, drugs positioned outside near their targets
    print("Applying layout...")

    # Position hub genes in a ring
    hub_list = sorted(gene_nodes)
    n_hubs = len(hub_list)
    radius = 300
    center_x, center_y = 0, 0
    hub_x = []
    hub_y = []
    hub_angles = {}
    for i, gene in enumerate(hub_list):
        angle = 2 * math.pi * i / n_hubs
        hub_angles[gene] = angle
        hub_x.append(center_x + radius * math.cos(angle))
        hub_y.append(center_y + radius * math.sin(angle))
    p4c.set_node_position_bypass(hub_list, hub_x, hub_y)
    time.sleep(1)

    # Position drugs outside the ring, grouped by their primary target gene
    # Build a map of gene -> list of drugs targeting it
    gene_to_drugs = {}
    for drug_node in drug_nodes:
        neighbors = [n for n in G_drug.neighbors(drug_node) if n in gene_nodes]
        if neighbors:
            primary = neighbors[0]
            if primary not in gene_to_drugs:
                gene_to_drugs[primary] = []
            gene_to_drugs[primary].append(drug_node)

    drug_positions_x = []
    drug_positions_y = []
    drug_names = []
    for gene, drugs in gene_to_drugs.items():
        base_angle = hub_angles.get(gene, 0)
        n_drugs = len(drugs)
        # Spread drugs in a fan around the gene's angle
        spread = min(0.3, 0.15 * n_drugs)
        for j, drug_node in enumerate(drugs):
            if n_drugs == 1:
                angle = base_angle
            else:
                angle = base_angle + spread * (j / (n_drugs - 1) - 0.5)
            outer_r = radius + 140 + (j % 3) * 50
            drug_names.append(drug_node)
            drug_positions_x.append(center_x + outer_r * math.cos(angle))
            drug_positions_y.append(center_y + outer_r * math.sin(angle))

    if drug_names:
        p4c.set_node_position_bypass(drug_names, drug_positions_x, drug_positions_y)
    time.sleep(1)

    # Export (Cytoscape 3.10+ API)
    p4c.fit_content()
    time.sleep(1)
    out_path = str(OUTPUT_DIR / "drug_hub_target_figure.png")
    print(f"Exporting to {out_path}")
    p4c.export_image(
        out_path,
        type="PNG",
        zoom=400,
        all_graphics_details=True,
        overwrite_file=True,
    )
    time.sleep(1)
    print(f"Figure 2 exported: {out_path}")
    return suid


# ---------------------------------------------------------------------------
# Phase D: Figure 3 — Community Overview
# ---------------------------------------------------------------------------


def prepare_community_network(data):
    """Prepare the full network with proper community attributes."""
    print("\n--- Figure 3: Preparing community overview ---")

    G_full = data["G_full"]
    hub_genes = data["hub_genes"]

    # Map singleton communities (0, 8-13) to nearest major community
    # by checking which major community their single PPI neighbor belongs to
    G_comm = G_full.copy()
    for node in G_comm.nodes():
        attrs = G_comm.nodes[node]
        comm = int(attrs.get("community", 0))
        if comm not in COMMUNITY_COLORS:
            # Singleton — assign to most-connected major community neighbor
            best_comm = 1
            best_count = 0
            for neighbor in G_comm.neighbors(node):
                n_comm = int(G_comm.nodes[neighbor].get("community", 0))
                if n_comm in COMMUNITY_COLORS:
                    count = sum(
                        1
                        for nn in G_comm.neighbors(node)
                        if int(G_comm.nodes[nn].get("community", 0)) == n_comm
                    )
                    if count > best_count:
                        best_count = count
                        best_comm = n_comm
            comm = best_comm

        disease = attrs.get("disease", "Both")
        is_hub = node in hub_genes
        G_comm.nodes[node]["community"] = str(comm)
        G_comm.nodes[node]["disease"] = disease
        G_comm.nodes[node]["is_hub"] = is_hub
        G_comm.nodes[node]["disease_color"] = DISEASE_COLORS.get(
            disease, DISEASE_COLORS_DEFAULT
        )
        G_comm.nodes[node]["community_color"] = COMMUNITY_COLORS.get(
            comm, COMMUNITY_COLORS_DEFAULT
        )

    print(f"Community network: {G_comm.number_of_nodes()} nodes, {G_comm.number_of_edges()} edges")
    return G_comm


def render_community_overview(G_comm, data):
    """Load community overview into Cytoscape with community-aware layout."""
    print("Loading into Cytoscape...")

    hub_genes = data["hub_genes"]
    community_df = data["community_df"]
    style_name = "CommunityOverview"

    suid = p4c.create_network_from_networkx(G_comm, title="Community Overview")
    time.sleep(3)

    # Create visual style
    defaults = {
        "NODE_SHAPE": "ELLIPSE",
        "NODE_SIZE": 12,
        "NODE_FILL_COLOR": "#CCCCCC",
        "NODE_BORDER_WIDTH": 2,
        "NODE_BORDER_PAINT": "#666666",
        "NODE_LABEL_FONT_SIZE": 9,
        "NODE_LABEL_COLOR": "#222222",
        "NODE_TRANSPARENCY": 230,
        "EDGE_WIDTH": 0.5,
        "EDGE_TRANSPARENCY": 50,
        "EDGE_STROKE_UNSELECTED_PAINT": "#AAAAAA",
        "NETWORK_BACKGROUND_PAINT": "#FFFFFF",
    }
    p4c.create_visual_style(style_name, defaults=defaults)
    p4c.set_visual_style(style_name)
    time.sleep(1)

    # Color nodes by disease
    p4c.set_node_color_mapping(
        "disease",
        table_column_values=["ASD", "SCZ", "Both"],
        colors=["#4A90D9", "#D94A4A", "#9B59B6"],
        mapping_type="d",
        default_color="#999999",
        style_name=style_name,
    )

    # Border color by community
    comm_values = [str(i) for i in range(1, 8)]
    comm_colors = [COMMUNITY_COLORS[i] for i in range(1, 8)]
    p4c.set_node_border_color_mapping(
        "community",
        table_column_values=comm_values,
        colors=comm_colors,
        mapping_type="d",
        default_color="#AAAAAA",
        style_name=style_name,
    )

    # Hub genes: larger, labeled
    all_nodes = list(G_comm.nodes())
    hub_node_list = [n for n in all_nodes if n in hub_genes]
    non_hub_list = [n for n in all_nodes if n not in hub_genes]

    if hub_node_list:
        p4c.set_node_size_bypass(hub_node_list, [35] * len(hub_node_list))
    # Non-hub: keep default size (12)

    # Labels: only hub genes
    p4c.set_node_label_mapping("name", style_name=style_name)
    time.sleep(0.5)
    if non_hub_list:
        p4c.set_node_label_bypass(non_hub_list, [""] * len(non_hub_list))
    if hub_node_list:
        p4c.set_node_font_size_bypass(hub_node_list, [14] * len(hub_node_list))

    # Community-aware layout: position nodes by community clusters
    print("Applying community-aware layout...")

    # First, compute community centroids in a circle
    community_ids = sorted(set(int(G_comm.nodes[n]["community"]) for n in G_comm.nodes()))
    n_comms = len(community_ids)
    comm_radius = 350
    comm_centroids = {}
    for i, cid in enumerate(community_ids):
        angle = 2 * math.pi * i / n_comms
        comm_centroids[cid] = (
            comm_radius * math.cos(angle),
            comm_radius * math.sin(angle),
        )

    # Group nodes by community
    comm_nodes = {}
    for node in G_comm.nodes():
        cid = int(G_comm.nodes[node]["community"])
        if cid not in comm_nodes:
            comm_nodes[cid] = []
        comm_nodes[cid].append(node)

    # Position nodes: arrange each community as a filled disk (spiral) around its centroid
    all_positioned_nodes = []
    all_x = []
    all_y = []
    for cid, nodes in comm_nodes.items():
        cx, cy = comm_centroids.get(cid, (0, 0))
        n = len(nodes)
        max_radius = max(30, math.sqrt(n) * 18)
        for j, node in enumerate(nodes):
            # Distribute in concentric rings for a filled disk look
            ring = int(math.sqrt(j))
            ring_start = ring * ring
            ring_size = 2 * ring + 1
            pos_in_ring = j - ring_start
            ring_radius = max_radius * (ring + 1) / (int(math.sqrt(n)) + 1)
            angle = 2 * math.pi * pos_in_ring / ring_size
            x = cx + ring_radius * math.cos(angle)
            y = cy + ring_radius * math.sin(angle)
            all_positioned_nodes.append(node)
            all_x.append(x)
            all_y.append(y)

    # Set all positions at once
    p4c.set_node_position_bypass(all_positioned_nodes, all_x, all_y)
    time.sleep(2)

    # Export (Cytoscape 3.10+ API)
    p4c.fit_content()
    time.sleep(1)
    out_path = str(OUTPUT_DIR / "community_overview_figure.png")
    print(f"Exporting to {out_path}")
    p4c.export_image(
        out_path,
        type="PNG",
        zoom=400,
        all_graphics_details=True,
        overwrite_file=True,
    )
    time.sleep(1)
    print(f"Figure 3 exported: {out_path}")
    return suid


# ---------------------------------------------------------------------------
# Phase E: Summary
# ---------------------------------------------------------------------------


def print_summary():
    """Print summary of generated figures."""
    print("\n" + "=" * 60)
    print("CYTOSCAPE FIGURE GENERATION COMPLETE")
    print("=" * 60)

    figures = [
        ("hub_subnetwork_figure.png", "Hub Gene Subnetwork"),
        ("drug_hub_target_figure.png", "Drug-Hub Target Map"),
        ("community_overview_figure.png", "Community Overview"),
    ]

    for fname, title in figures:
        fpath = OUTPUT_DIR / fname
        if fpath.exists():
            size_mb = fpath.stat().st_size / (1024 * 1024)
            print(f"  {title:30s} -> {fname} ({size_mb:.1f} MB)")
        else:
            print(f"  {title:30s} -> MISSING: {fname}")

    # Update kg_results_summary.json
    summary_path = OUTPUT_DIR / "kg_results_summary.json"
    if summary_path.exists():
        with open(summary_path) as f:
            summary = json.load(f)
        summary["cytoscape_figures"] = {
            "hub_subnetwork_figure": "hub_subnetwork_figure.png",
            "drug_hub_target_figure": "drug_hub_target_figure.png",
            "community_overview_figure": "community_overview_figure.png",
            "generated_with": "py4cytoscape + Cytoscape 3.10.4",
        }
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"\n  Updated: {summary_path}")

    print("\nDone.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print("=" * 60)
    print("Cytoscape Network Figure Generator")
    print("=" * 60)

    # Phase A
    verify_cytoscape()
    data = load_data()
    clear_cytoscape_session()

    # Phase B: Figure 1
    G_hub = build_hub_subnetwork(data)
    render_hub_subnetwork(G_hub, data)

    # Phase C: Figure 2
    G_drug = build_drug_hub_network(data)
    render_drug_hub_network(G_drug, data)

    # Phase D: Figure 3
    G_comm = prepare_community_network(data)
    render_community_overview(G_comm, data)

    # Phase E: Summary
    print_summary()


if __name__ == "__main__":
    main()

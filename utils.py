import os
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import OrderedDict
from copy import deepcopy
from warnings import warn
from typing import List, Optional, Dict, Union, Tuple

def load_data(
        gene:str,
        species:str="homo_sapiens",
        preprocess:bool=True,
        filepath:Optional[os.PathLike]=None,
        ) -> Tuple[pd.DataFrame]:
    """ Load data for a given gene
    
    Anticipated file structure: 
        Relative:   ./{gene}/{species}/
        Absolute:   {filepath}/{gene}/{species}/
    
    The raw codependency data is expected in the `gene` folder:
        - ./{gene}_codep.csv
        
    The following String-db output files are expected in the `species` folder:
        - ./{gene}/{species}/{gene}_string_MCL_clusters.tsv
        - ./{gene}/{species}/{gene}_network_coords_pos.tsv
        - ./{gene}/{species}/{gene}_string_interactions.tsv
    
    Args:
        gene : str
            string representation of a gene. Ensure case matches folder structure.
        species : str, default="homo_sapiens"
            Should match the above file structure
        preprocess: bool , default=True
            If True, all datasets will be preprocessed before returning.
            Not recommended if using with data outside this repository
        filepath : os.PathLike, default=None
            If specified, this filepath will be utilized rather than relative path from 
            current script. The above file structure is still expected relative to `filepath`
    
    Returns:
        data : Dict[str, pd.DataFrame]
    """
    root = os.path.join(".", f"{gene}")
    if filepath is not None:
        root = filepath

    codep = pd.read_csv(os.path.join(root, f"{gene}_codep.csv"))
    clusters = pd.read_csv(os.path.join(root, f"{species}/{gene}_string_MCL_clusters.tsv"), sep="\t")
    coords  = pd.read_csv(os.path.join(root, f"{species}/{gene}_string_network_coordinates.tsv"), sep="\t")
    interactions = pd.read_csv(os.path.join(root, f"{species}/{gene}_string_interactions.tsv"), sep="\t")
    
    if preprocess:
        codep = preprocess_codep(codep)
        clusters = preprocess_clusters(clusters)
        coords = preprocess_coords(coords)
        interactions = preprocess_interactions(interactions)

    return codep, clusters, coords, interactions


def preprocess_codep(codep:pd.DataFrame, cmap:Union[mpl.colors.Colormap, str]="viridis"):
    """ Prepreocess raw codependency data
    Not recommended as general purpose preprocessing function.

    1) subset columns
    2) fill nas
    3) calculate relative correlation strengths 
    4) Get edge color based on correlation
    5) Sort
    6) Convert `gene` column to uppercase values.

    Args:
        codep : pd.DataFrame
            The raw codependency dataframe. 
        cmap : Union[mpl.colormaps, str]
            The mpl colormap used for edge color
    """
    # Handle Union inputs
    if isinstance(cmap, str):
        cmap = mpl.colormaps[cmap]

    # Subest raw data to important columns, fill NAs:
    cols = ["gene", "correlation"]
    codep = codep[cols].fillna("na")

    # Correlation based edge color
    codep['abs_corr'] = abs(codep['correlation'])
    codep['relative_corr'] = (codep["abs_corr"] - min(codep['abs_corr'])) / (max(codep['abs_corr']) - min(codep['abs_corr']))
    codep['log_relative_corr'] = np.log(codep['relative_corr']+1.0)
    codep.sort_values(by = 'correlation', inplace=True)
    codep['edge_color'] = [(val[0], val[1], val[2]) for idx, val in enumerate(cmap(sorted(codep['correlation'])) ) ]
    codep.sort_values(by = ['relative_corr'], ascending =False, inplace = True)

    # Convert all genes to upper case 
    codep['gene'] = [gene.upper() for gene in codep['gene']]
    return codep

def preprocess_clusters(clusters:pd.DataFrame, gene_thres:int=5):
    """ Preprocess clusters dataframe.
    1) Apply gene threshold
    2) Convenience rename columns
    
    Args:
        clusters : pd.DataFrame
            The clusters dataframe. Expects a String-db string_MCL_clusters.tsv file structure
        gene_thres: int, default=5
            A cluster must have at least `gene_thres` genes to be included.
    """
    # Apply Gene count cluster threshold 
    clusters = clusters[clusters['gene count'] >= gene_thres]
    # Convenience rename
    clusters.rename({'protein name': 'protein'}, axis = 1, inplace=True)
    # Make upper:
    clusters['protein'] = [protein.upper() for protein in clusters['protein']]
    return clusters

def preprocess_coords(coords:pd.DataFrame):
    """ Preprocess coords dataframe.
    1) Convenience rename columns
    
    Args:
        coords : pd.DataFrame
            The clusters dataframe. Expects a String-db _network_coords_pos.tsv file structure
    """
    # Convenience rename
    coords.rename({"#node":"protein"}, axis = 1, inplace=True)
    # Make upper:
    coords['protein'] = [protein.upper() for protein in coords['protein']]
    return coords


def preprocess_interactions(interactions:pd.DataFrame):
    """ Preprocess interactions dataframe.
    1) Subset columns
    2) Convenience rename columns
    
    Args:
        clusters : pd.DataFrame
            The clusters dataframe. Expects a String-db _string_interactions.tsv file structure
    """
    # Subest raw data to important columns:
    cols = ["#node1", "node2", "coexpression", "combined_score"]
    interactions = interactions[cols]
    # Convenience renaming
    interactions.rename({"#node1":"protein", "node2":"n2"}, axis = 1, inplace=True)
    # Make upper:
    interactions['protein'] = [protein.upper() for protein in interactions['protein']]
    interactions['n2'] = [n2.upper() for n2 in interactions['n2']]
    return interactions


def get_within_cluster_interactions(gene:str, clusters:pd.DataFrame, interactions:pd.DataFrame)->pd.DataFrame:
    cluster_interactions = filterInput(interactions, inclusive_filters= [("protein", list(clusters['protein'].unique()) ), ("n2", list(clusters['protein'].unique()) )] )
    
    # Clkuster map
    cluster_map = {gene.upper():cluster_id for gene, cluster_id in zip(clusters['protein'], clusters['cluster number'])}
    cluster_map[f"{gene}"] = 0
    
    def _check_cluster_membership(l, r):
        return l in cluster_map.keys() and r in cluster_map.keys() and cluster_map[l] == cluster_map[r]
    
    cluster_interactions['same_cluster'] = [_check_cluster_membership(l.upper(), r.upper()) for l, r in zip(cluster_interactions['protein'], cluster_interactions['n2']) ]
    cluster_interactions = cluster_interactions[cluster_interactions['same_cluster']]
    cluster_interactions['cluster_id'] = [cluster_map[gene.upper()] for gene in cluster_interactions['protein']]
    return cluster_interactions, cluster_map


def filterInput(pdf:pd.DataFrame, inclusive_filters:List[tuple] = None, exclusive_filters: List[tuple] = None):
    """ Apply a list of tuple filters to a pandas dataframe and return copy

    Parameters
    ----------
        pdf : pd.DataFrame
            Pandas Dataframe
        inclusive_filters : List[tuple], default = None
            A list-like object containing tuple filters.
            return INCLUDES the filter values.

            Each tuple should have the structure:
                t[0] : The column of `pdf` on which to apply filter
                t[1] : The inclusion criterion
                    If the inclusion criterion is an iterable, `.isin()`
                    is applied
                    otherwise equivalence is applied (i.e. `==`)
        exclusive_filters : List[tuple], default = None
            A list-like object containing tuple filters.
            return EXCLUDES the filter values.
            Each tuple should have the structure:
                t[0] : The column of `pdf` on which to apply filter
                t[1] : The inclusion criterion
                    If the inclusion criterion is an iterable, `.isin()`
                    is applied
                    otherwise equivalence is applied (i.e. `==`)
    
    Returns
    -------
        output : pd.DataFrame
            Copy of input dataframe meeting the filtering criterion.
       
    Warns
    -----
        If an a filter column name (i.e. t[0]) is not present in `pdf` the
        current implementation just warns the user in the console and skips
        the filter. 
    
    """
    if inclusive_filters is None and exclusive_filters is None:
        warn(f"Both Inclusive and Exclusive Filters are empty... Returning unfiltered input.")
        return pdf

    output = deepcopy(pdf)

    if inclusive_filters is not None:
        for key, val in inclusive_filters:
            try:
                # Check for inclusion if iterable
                if hasattr(val, '__iter__'):
                    output = output.loc[output[f"{key}"].isin(val)]
                # o/w check equivalence
                else:
                    output = output.loc[output[f"{key}"] == val]
            except:
                warn(message = f"Uh Oh! --> {key} <-- not within the input column space! Skipping Filter...")
    
    if exclusive_filters is not None:
        for key, val in exclusive_filters:
            try:
                # Check for inclusion if iterable
                if hasattr(val, '__iter__'):
                    output = output.loc[~(output[f"{key}"].isin(val))]
                # o/w check equivalence
                else:
                    output = output.loc[~(output[f"{key}"] == val)]
            except:
                warn(message = f"Uh Oh! --> {key} <-- not within the input column space! Skipping Filter...")
            
    return output


# https://stackoverflow.com/questions/20924085/python-conversion-between-coordinates
def cart2pol(x:float, y:float):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return (rho, phi)

def pol2cart(rho:float, phi:float):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)


def fGetColorsFromPDColumn(pdf:pd.DataFrame, cmap = 'viridis'):
    """Return OrderedDict of (r, g, b) values along colormap for unique column entries
     
    Parameters
    ----------
        pdf : pd.DataFrame
            Pandas Dataframe column 
        cmap : str = 
            Name of a matplotlib.colormap 
    
    Returns
    -------
        output : collections.OrderedDict
            Sorted Dictionary of (r, g, b) linearly distributed along input colormap
       
    See Also
    --------
    >>> https://matplotlib.org/stable/users/explain/colors/colormaps.html

    """
    colors = plt.colormaps[cmap](np.linspace(0.01, 0.99, pdf.shape[0]))
    output = OrderedDict()
    for idx, val in enumerate(sorted(pdf) ):
        output[val] = (colors[idx][0], colors[idx][1], colors[idx][2])
    return output


def plot_clusters(
        node_df:pd.DataFrame, edge_df:pd.DataFrame,
        cluster_map:Dict[str,int], node_colors:Dict[str,int],
        position:Dict[str,Tuple[float,float]],
        root:Optional[str]=None,
        root_size:int=3000,
        neighbor_scaling:Optional[Union[int,float]]=2.5,
        neighbor_size:Optional[int]=None,
        clusters_to_include:Optional[List[int]]=None,
        boolUseDefaultEdgeColor:bool=True,
        boolUseDefaultNodeColor:bool=False,
        boolLabelNodes:bool=True,
        boolAddRootEdge:bool=False,
        default_node_color:str = "#00A5CF",
        default_edge_color:str = "#bfbfbf", #"#E0E0DE",
        figsize:Tuple[int,int] = (12, 12),
        edge_size_min:float=1.0,
        edge_corr_scale:float=20.0,
        boolSpringlayout:bool=False,
        springk:Optional[float]=None,
        width_col = "corr",
        
    ) -> Union[plt.Figure, plt.Axes]:
    """
    Args:
        clusters_to_include : default None
            If specified, plotted nodes will have cluster_ids in this list's values
        root_size : 
    
    """
    # Setup containers
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax = np.ravel(ax)
    G = nx.Graph()
    
    # Get scaling
    if neighbor_size is None:
        neighbor_size = root_size / neighbor_scaling

    if clusters_to_include is not None:
        genes_to_include = [g.upper() for g in edge_df[edge_df['cluster_id'].isin(clusters_to_include) ]['protein'].unique() ]
    else:
        genes_to_include = [g.upper() for g in edge_df['protein'].unique()] + [root]
    # Gene must have a position and 
    include_gene = set(genes_to_include) & set(position.keys())

    if root is not None:
        G.add_node(root, pathway = [""], size=root_size, is_neighbor=0, cluster_id=0)
    
    # Add all nodes
    for row in range(node_df.shape[0]):
        gene = node_df["gene"][row].upper()

        if gene in include_gene:
            G.add_node(
                gene,
                relative_corr = node_df["relative_corr"][row],
                size=neighbor_size, is_neighbor=1,
                cluster_id=cluster_map[gene],
                )
            if boolAddRootEdge:
                G.add_edge(root, gene, corr=node_df["correlation"][row], relcorr=node_df["relative_corr"][row])

    # Add Edges
    for n1, n2, coex, score in zip(
                edge_df['protein'],
                edge_df['n2'],
                edge_df['coexpression'],
                edge_df['combined_score']
        ):
        l, r = n1.upper(), n2.upper()
    
        # Add w/in cluster edges
        if l in genes_to_include and r in genes_to_include:
            G.add_edge(l, r, coex=coex, score=score, corr=0, relcorr=0)

    # Edge Colors 
    # Add Alpha to hex: https://gist.github.com/lopspower/03fb1cc0ac9f32ef38f4
    for u, v, edgedata in G.edges(data = True):
        if u == root or v == root:
            if edgedata['corr'] > 0:
                edgedata['edge_color'] = "#ff000066"
            else:
                edgedata['edge_color'] = "#0000ff66"
        else:
            edgedata['edge_color'] = default_edge_color #"lightgray"
    
    if boolSpringlayout:
        # Spring_layout weight
        for u, v, edgedata in G.edges(data = True):
            edgedata['weight'] = abs(edge_corr_scale*edgedata[width_col])+edge_size_min if width_col is not None else 1
        position = nx.spring_layout(G, k=springk, pos=position, weight='weight', center=position[root])

    nx.draw_networkx(
        G = G,
        ax = ax[0],
        with_labels = boolLabelNodes,
        node_size = [nodedata.get('size', neighbor_size) for node, nodedata in G.nodes(data=True)],
        width = [edge_corr_scale*edgedata[width_col]+edge_size_min for _,_ , edgedata in G.edges(data=True)] if width_col is not None else None,
        edge_color = default_edge_color if boolUseDefaultEdgeColor else [edgedata["edge_color"]  for _,_ , edgedata in G.edges(data=True)],
        node_color = default_node_color if boolUseDefaultNodeColor else [node_colors[nodedata['cluster_id']] for _, nodedata in G.nodes(data=True)],
        pos = position
        )

    return fig, ax
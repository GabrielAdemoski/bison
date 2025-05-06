import pandas as pd
import numpy as np
import os
import igraph as ig

from src.circular_plot import circular_plot
from src.parallel_plot import parallel_plot

def sequential(a,b,lst,cycle=True):
    """
    Takes in two values and a list they belong in and determines if they are sequential.
    :param a: A value in lst.
    :param b: A value in lst.
    :param lst: A list of values in sequential order.
    :param cycle: A boolean, are a and b a part of a cycle.
    :return: A boolean value where True means values are sequential and False means they are not.
    """
    seq = False
    if a in lst and b in lst:
        lenl = len(lst)
        if cycle:
            for l in range(lenl):
                if a == lst[l % lenl] and b == lst[(l + 1) % lenl]:
                    seq = True
        else:
            for l in range(lenl - 1):
                if a == lst[l] and b == lst[(l + 1)]:
                    seq = True
    return seq

def stitch(lst,cycle=True):
    """
    Retuen A list of tuples where each tuple contained 2 sequential values.
    :param lst: A list of values in order.
    :param cycle: A boolean, is lst a cycle.
    :return: A list of tuples where each tuple contained 2 sequential values.
    """
    s = []
    lenl = len(lst)
    if cycle:
        for l in range(lenl):
            s.append((lst[l % lenl], lst[(l + 1) % lenl]))
    else:
        for l in range(lenl - 1):
            s.append((lst[l], lst[(l + 1)]))
    return s

def sorted_dict(lvls, cats, cats_dict):
    """
    Sorts and pairs sequential catagory entries.
    :param lvls: A list of levels.
    :param cats: A list of catagory lists and their cycle flags.
    :param cats_dict: A dictionary of numbers and their associated catagory.
    :return seq: A dictionary of dictionaries with the first keys being levels, the second keys being catagories, the
    third keys are a consecutive pair of years. Each pair of years contains a list of csv names which are associated
    with the pair of years.
    :return lst: A list of all csv names above.
    """
    lst = []
    seq = {}
    for l in lvls:
        seq[l] = {}
        for c in range(len(cats)):
            seq[l][cats_dict[c]] = {}
            # sticth the list within cats[c] together in order.
            s = stitch(*cats[c])
            for frm, to in s:
                k = frm + '-' + to
                seq[l][cats_dict[c]][k] = []
                # Add both csv fnames to the dict and lst
                for x in (frm, to):
                    lst.append(l + "_" + frm + "-" + to + "_" + x + ".csv")
                    seq[l][cats_dict[c]][k].append(l + "_" + frm + "-" + to + "_" + x + ".csv")
    return seq, lst

def common_taxa(seq, data_dir, cats_dict):
    """
    Generate a dictionary of common taxa by level and catagory.
    :param seq: A dictionary of dictionaries with the first keys being levels, the second keys being catagories, the
    third keys are a consecutive pair of years. Each pair of years contains a list of csv names which are associated
    with the pair of years.
    :param data_dir: A path for the data directory.
    :param cats_dict: A dictionary of numbers and their associated catagory.
    :return: A dictionary where first keys are levels, second keys are catagories and values are lists of taxa names.
    """
    common = {}
    for l in lvls:
        common[l] = {}
        for cat in cats_dict.values():
            common[l][cat] = []
            for x in list(seq[l][cat].keys()):
                files = seq[l][cat][x]
                for fname in files[::-1]:
                    path = data_dir + fname
                    df = pd.read_csv(path, index_col=0, delimiter=" ")
                    common[l][cat].append(df.columns.tolist())
            common[l][cat] = list(set([j for i in common[l][cat] for j in i]))
    return common


def edge_color(weight):
    """
    Colors edges based on value.
    :param weight: A float representing an edge weight.
    :return: A string of the color name.
    """
    if weight > 0:
        return "orange"
    elif weight < 0:
        return "blue"
    else:
        return "black"


def topn(matrix, n=10):
    '''
    Selects the n highest magnitude edges per 2D numpy matrix.
    :param matrix: A 3D numpy matrix.
    :param n: An intiger for the max number of edges to keep.
    :return: A 3D numpy matrix but with at most n edges per layer.
    '''
    for z in range(matrix.shape[0]):
        array_2d = matrix[z]
        array_1d  = array_2d.flatten()
        top_n_indices = np.argpartition(np.abs(array_1d), -n)[-n:]

        # Create a mask with True for top n indices and False otherwise
        mask = np.zeros_like(array_1d, dtype=bool)
        mask[top_n_indices] = True
        result_1d = np.where(mask, array_1d, 0)

        # Reshape the 1D array back to the original 2D shape
        matrix[z] = result_1d.reshape(array_2d.shape)
    return matrix



def graph(path, figs_dir, l, c, filter=0):
    """
    Create an undirected, weighted igraph object with no singletons and colored edges.
    :param path: Path of the correlation matrix.
    :param figs_dir: A path to the figs directory.
    :param l: A sting representing the level.
    :param c: A string representing the catagory.
    :param filter: A float value [0,1] which is used to filter out values.
    :return:
    """
    if os.path.exists(path):
        # load a correlation matrix
        corr = pd.read_csv(path, delimiter=" ")
        val = corr.values
        mx = np.max(np.abs(val))
        # keep only correlation values whose magnitude are in the top 'filter' percent of edge weights.
        adj = ((np.abs(val) > mx*filter) * val).tolist()
        # create an undirected, weighted igraph object
        graph = ig.Graph.Weighted_Adjacency(adj, mode="undirected")
        graph.vs['name'] = list(corr.columns)

        # drop singletons
        to_delete = [v.index for v in graph.vs if v.degree() == 0]
        graph.delete_vertices(to_delete)

        # color the egdes of the graph
        edge_colors = [edge_color(weight) for weight in graph.es['weight']]
        fname = figs_dir + l + "_" + c + '_' + str(filter) + ".png"

        # A dictionary of visual style elements for igraph.plot
        visual_style = {}
        visual_style["vertex_shape"] = 'circle'
        # visual_style["vertex_size"] = graph.betweenness()
        visual_style["layout"] = graph.layout("kk")
        visual_style["bbox"] = (1000, 1000)
        visual_style["vertex_label_dist"] = -3
        visual_style["margin"] = 120
        visual_style["label_dist"] = 20
        visual_style["vertex_label"] = graph.vs["name"]
        visual_style["edge_color"] = edge_colors
        visual_style["edge_width"] = [5 * (abs(w) / mx) for w in graph.es['weight']]

        ig.plot(graph, target=fname, **visual_style)







# corr
data_dir = "data/corr/"
figs_dir = "figs/corr/"
os.makedirs(figs_dir, exist_ok=True)

lvls = ["lvl2", "lvl6"]
cats = [[["2018", "2019", "2020", "2021"], False], [["Fall", "Winter", "Spring", "Summer"], True]]

for l in lvls:
    for cat in cats:
        for c in cat[0]:
            data_path = data_dir + l + "_" + c + ".csv"
            graph(data_path, figs_dir, l, c, filter=0)
            if l == 'lvl6':
                graph(data_path, figs_dir, l, c, filter=.9)

data_dir = "data/delta_corr/"
figs_dir = "figs/delta/"
os.makedirs(figs_dir, exist_ok=True)

files = os.listdir(data_dir)
# sort by level and by catagory
files.sort()

lvls = ["l2", "l6"]
lvls_dict = {'l2': 'Phylum', 'l6': 'Genus'}
cats_dict = {0:"Years", 1:"Seasons"}
cats = [[["2019", "2020", "2021"], False], [["Fall", "Winter", "Spring", "Summer"], True]]

seq, lst = sorted_dict(lvls, cats, cats_dict)
common = common_taxa(seq, data_dir, cats_dict)

non_seq = []
for file in files:
    if file not in lst:
        non_seq.append(file)

change = {}
for l in lvls:
    change[l] = {}
    for cat in cats_dict.values():
        change[l][cat] = {}
        for x in list(seq[l][cat].keys()):
            files = seq[l][cat][x]
            diff = pd.DataFrame(index=sorted(common[l][cat]), columns=sorted(common[l][cat]))
            # iterate through files from last to first
            for fname, idx in zip(files[::-1], np.arange(len(files[::-1]))):
                path = data_dir + fname
                df = pd.read_csv(path, index_col=0, delimiter=" ")
                if idx == 0:
                    for col in df.columns.tolist():
                        for row in df.index.tolist():
                            diff.loc[row, col] = df[col][row]
                else:
                    for col in df.columns.tolist():
                        for row in df.index.tolist():
                            diff.loc[row, col] -= df[col][row]
            change[l][cat]["taxa"] = diff.columns.tolist()
            change[l][cat][x] = np.triu(diff.infer_objects(copy=False).fillna(0), 1)

# delta
n=10
verbose = False
for l in lvls:
    for cat in cats_dict.keys():
        cycle = cats[cat][-1]
        matrix = list(change[l][cats_dict[cat]].values())[1:]
        taxa = list(change[l][cats_dict[cat]]["taxa"])
        matrix = np.dstack(matrix)
        matrix = np.rollaxis(matrix, -1)
        fname = lvls_dict[l] + " " + cats_dict[cat]

        if cycle:
            if verbose:
                print(fname)
                for z in range(matrix.shape[0]):
                    print(np.sum(matrix[z] > 0), np.sum(matrix[z]< 0), np.sum(np.abs(matrix[z]) > 0))
                print(np.sum(matrix > 0), np.sum(matrix < 0), np.sum(np.abs(matrix) > 0))
                print(taxa)
            circular_plot(matrix, taxa, disp=False, path=figs_dir, fname=fname, title=lvls_dict[l])

            topn_matrix = topn(matrix, n=n)
            fname = fname + " Top " + str(n)
            if verbose:
                print(fname)
                for z in range(topn_matrix.shape[0]):
                    print(np.sum(topn_matrix[z] > 0), np.sum(topn_matrix[z]< 0), np.sum(np.abs(topn_matrix[z]) > 0))
                print(np.sum(topn_matrix > 0), np.sum(topn_matrix < 0), np.sum(np.abs(topn_matrix) > 0))
                print(taxa)

            circular_plot(topn_matrix, taxa, disp=False, path=figs_dir, fname=fname, title=lvls_dict[l])
        else:
            if verbose:
                print(fname)
                for z in range(matrix.shape[0]):
                    print(np.sum(matrix[z] > 0), np.sum(matrix[z] < 0), np.sum(np.abs(matrix[z]) > 0))
                print(np.sum(matrix > 0), np.sum(matrix < 0), np.sum(np.abs(matrix) > 0))
                print(taxa)
            parallel_plot(matrix, taxa, disp=False, path=figs_dir, fname=fname, title=lvls_dict[l])

        
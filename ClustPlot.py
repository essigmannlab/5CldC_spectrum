# ClustPlot
#
# B.Fedeles
#
# Functions that enable clustering and plotting of histograms for mutational data.
# Based on code by L. Kim and C. Valentine.
#
# Histogram plotting uses seaborn library    
#
# Last updated: 2020-08-16.


from collections import OrderedDict
from scipy.cluster.hierarchy import dendrogram
from fastcluster import linkage
from sklearn.metrics.pairwise import cosine_similarity

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt



def refine(mut_sig):
    # ensures OrderedDicts are in same order w.r.t keys
    # rewritten to allow for COSMIC dict input
    # cleans up dict so there's no impossibilities (ex: ACA, T>A)

    mutations = ['C>A','C>G','C>T','T>A','T>C','T>G']
    contexts = ['ACA','ACC','ACG','ACT',
                'CCA','CCC','CCG','CCT',
                'GCA','GCC','GCG','GCT',
                'TCA','TCC','TCG','TCT',
                'ATA','ATC','ATG','ATT',
                'CTA','CTC','CTG','CTT',
                'GTA','GTC','GTG','GTT',
                'TTA','TTC','TTG','TTT']
    new_sig = OrderedDict()
    if isinstance(list(mut_sig.keys())[0],str):
        for sig in list(mut_sig.keys()):
            inner_sig = OrderedDict()
            new_sig[sig] = inner_sig
            for m in mutations:
                for c in contexts:
                    if m[0] == c[1]:
                        inner_sig[(m,c)] = mut_sig[sig].setdefault((m,c),0)
        return new_sig
    for m in mutations:
        for c in contexts:
            if m[0] == c[1]:
                new_sig[(m,c)] = mut_sig.setdefault((m,c),0)
    return new_sig


def import_cosmic_sigs(in_file):
    # read COSMIC v2 sigs from the original file.
    
    spect_dict = OrderedDict()
    firstLine = True
    f = open(in_file, 'r')
    for line in f:
        line = line.strip().split('\t')
        if firstLine:
            for num in line[3:]: # for each signature, past first three columns
                sig = "signature_" + str(num[10:]) # signature number
                spect_dict[sig] = OrderedDict()
            firstLine = False
            continue
        for col in range(1,31):
            sig = "signature_" + str(col)
            spect_dict[sig][(line[0],line[1])] = float(line[col+2]) # key is (mutation, context)
    return refine(spect_dict)




def import_mutations(in_file):
    # Imports data from a text (CSV) file, with data in the order: substitution type, context, count, frequency.

    spect_dict = OrderedDict()
    firstLine = True
    f = open(in_file, 'r')
    for line in f:
        line = line.strip().split(',')[:4]
        if firstLine:
            firstLine = False
            continue
        #print(line)
        sub, con, count, freq = line
        spect_dict[(sub,con)] = float(freq)
    return refine(spect_dict)


def uhc_cluster(spec_list, ref_sig=None):
    # Unsupervised hierarchical clusters (uhc) from a collection of spectra.
    # Spec_list is a dictionary;
    # Uses linkage from fastcluster

    if ref_sig is not None:
        spectra = [list(ref_sig.values())] # so ref signature is value 0
    else:
        spectra = []
    for sig in spec_list:
        spectra.append(list(spec_list[sig].values()))
    return linkage(spectra, method='ward', metric='cosine')

def plot_uhc_heatmap(spec_list, cluster_names, isText, file_name, fmt, st_col):

    # Plot a clustermap with dendrogram and histogram heatmap.
    # Uses clustermap from seaborn
    # spec_list: a dictionary of spectra to be compared/plotted
    # cluster_names: the names of the spectra above; used for labeling plot rows/cols.
    # isText: a flag to display cossim values (True), or no display (False)
    # file_name: name of the figure file to be saved (with extension)
    # fmt: format of figure file (e.g. svg, pdf, eps)

    spectra = []
    for sig in list(spec_list.keys()):
        spectra.append(list(spec_list[sig].values()))
    linkages = uhc_cluster(spec_list)
    with plt.rc_context({'lines.linewidth': 1.25}):
        columns = cluster_names
        grid = sns.clustermap(pd.DataFrame(cosine_similarity(spectra), index=[columns], columns=[columns]),
                            method='weighted',
                            row_cluster=True, col_cluster=True,
                            row_linkage=linkages, col_linkage=linkages,
                            annot=isText, fmt='.2f', # change this line to show values in heatmap;
                            #fontsize='large',                 # .2 refers to number of decimals to show
                            cmap=sns.cubehelix_palette(start=st_col, rot=-0.1, dark=0.15, light=.55, as_cmap=True))

    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), size=14)
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), size=14, rotation=0)

    grid.ax_heatmap.get_xaxis().set_label_text(' ')
    grid.ax_heatmap.get_yaxis().set_label_text(' ')


  #  plt.show()  # change backend for this to work; use TkAgg instead of agg, for example.

    plt.savefig(file_name, format=fmt, dpi=450, bbox_inches='tight')

    plt.close(plt.gcf())  # important to close the figure once it's done...

def plot_uhc_dendrogram(spec_list, cluster_names, file_name, fmt):

    # As above but just the dendrogram.

    spectra = []
    for sig in spec_list:
        spectra.append(list(spec_list[sig].values()))
    data = spectra
    labels = cluster_names
    linkages = uhc_cluster(spec_list)
    with plt.rc_context({'lines.linewidth':0.5}):
        dendrogram(linkage(pd.DataFrame(data, [labels, ['']*len(labels)]),
                                 method='ward',
                                 metric='cosine'),
                         no_labels=False,
                         labels=labels,
                         leaf_rotation=90,
                         leaf_font_size=4,
                         link_color_func=lambda x:'k')
  #  plt.show()
    plt.savefig(file_name, format=fmt, dpi=450)
  
    plt.close(plt.gcf()) # important to close the figure once it's done...


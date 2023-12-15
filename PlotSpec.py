#!/usr/bin/env python3
#
# PlotSpec
# by Bogdan Fedeles
# v1.0 2021-05-12  
# 
# Library of functions to plot mutational spectra. Takes a .csv file with spectrum data and plots the 96bar spectrum.
# Features:
# - X-axis can show purines or pyrimidine mutations.
# - Normalize to a set of triplet contexts abundances. (provided in a separate file).
# - Display error bars (these need to be included in the input file.)
# - Title options for graph
#
# Code adapts some functions from draw_spectra.py by LK and CV.


#John's colors
# Pantone RGB colors:  AT R173/G169/B166 (#ADA9A6).  AG R150/G214/B77 (#96D64D).  AC R240/G156/B162 (#F09CA2).
# GT R82/G194/B242 (#52C2F2).  GC R36/G31/B33 (#241F21).  GA R230/G33/B36 (#E62124).  Error bars R181/G175/B175



##import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from statistics import stdev
from Bio import SeqIO
from collections import OrderedDict
from itertools import cycle, product
from matplotlib.patches import Rectangle
#from openpyxl.styles import Font, colors, PatternFill

dna_bases = ['A','C','G','T']

py_muts = ('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')
pu_muts = ('G>T', 'G>C', 'G>A', 'A>T', 'A>G', 'A>C')

colormap = {'COSMIC':['#52C3F1', '#231F20', '#E62223', '#CBC9C8', '#97D54C', '#EDBFC2'],
            'COSMIC2':['#2A60D2', '#15171F', '#BE2528', '#DBDDE3', '#1DBE06', '#F7A5D9'],
            'COSMIC3':['#52C2F2', '#241F21', '#E62124', '#ADA9A6', '#96D64D', '#F09CA2'],
            'custom':['blue', 'black', 'red', 'gray', 'green', 'pink']}

ccons =    ['ACA', 'ACC', 'ACG', 'ACT',
            'CCA', 'CCC', 'CCG', 'CCT',
            'GCA', 'GCC', 'GCG', 'GCT',
            'TCA', 'TCC', 'TCG', 'TCT']
tcons =    ['ATA', 'ATC', 'ATG', 'ATT',
            'CTA', 'CTC', 'CTG', 'CTT',
            'GTA', 'GTC', 'GTG', 'GTT',
            'TTA', 'TTC', 'TTG', 'TTT']

def fasta_to_dict(ref_file):
    """
    Import a fasta file into a dictionary.
    """
    with open(ref_file, 'r') as fi:
        return SeqIO.to_dict(SeqIO.parse(fi, 'fasta'))

def axes_onoff(ax, switch=False):
    ax.get_xaxis().set_visible(switch)
    ax.get_yaxis().set_visible(switch)
    return(ax)

def ticks_onoff(ax, switch1=False, switch2=False):
    for tic in ax.yaxis.get_major_ticks() + ax.xaxis.get_major_ticks():
        tic.tick1On = switch1
        tic.tick2On = switch2
    return(ax)

def spines_onoff(ax, switch=False):
    for spine in ['top', 'left', 'bottom', 'right']:
        ax.spines[spine].set_visible(switch)
    return ax

def init_chart(ax):
    """
    Setup a chart for spectra plotting.
    """
    ax.grid(False)
    ax.patch.set_facecolor('white')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ticks_onoff(ax, True, False)
    return ax

def rev_comp(sequence):
    """
    Returns the reverse complement of the sequence, a DNA sequence string.
    """
    
    basepair = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    revcomp = ''
    for base in sequence[::-1]:
        revcomp += basepair[base]
    return revcomp

def import_kmer_counts(file, notation='pyrimidine', k=3):
    """
    Read a jellyfish or bbmap generated file with counts of all kmers from a sequence of length k.
    Returns a dictionary of sequence counts for notation=purine/pyrimidine centric context.
    """

    allcontexts = {}
    f = open(file, 'r')
    i=0
    total=0
    for line in f:
        tri,count = line.split()
        allcontexts[tri] = int(count)
        i+=1
        total+=int(count)
   
    expected = 4**k      #4 to the 3rd power.
    found = i
    print ('Expected # of kmers :', expected)
    print ('Read from file : ', found)
    print ('Total # of occurences :', total)

    print ('Notation is :', notation)

    labels = ccons+tcons    
    if notation=='pyrimidine':
        contexts_keys = labels
    if notation=='purine':
        contexts_keys=[rev_comp(trimer) for trimer in labels]
  
    print(contexts_keys)

    contexts = {}
    for key in contexts_keys:
        contexts[key] = allcontexts[key]+allcontexts[rev_comp(key)]

    return contexts

def init_spec_dict(notation='pyrimidine'):
    """
    Initialize a spectrum OrderedDict, with 96 keys.
    Each key is a tuple (mutation, sequence context), ie: (C>A, AGA) etc.
    Notation indicates if mutations and contexts are purine/pyrimidine centric.
    """

    cmuts = list(py_muts[0:3])
    tmuts = list(py_muts[3:6])

    amuts = list(pu_muts[3:6])
    gmuts = list(pu_muts[0:3])
    
    if notation=='purine':
        #purine code
        acons = [rev_comp(con) for con in tcons]
        gcons = [rev_comp(con) for con in ccons]
        dictkeys = [(m,c) for m in gmuts for c in gcons] + [(m,c) for m in amuts for c in acons]  
    else:
        #pyrimidine code
        dictkeys = [(m,c) for m in cmuts for c in ccons] + [(m,c) for m in tmuts for c in tcons]
        
    specdict = OrderedDict.fromkeys(dictkeys, 0)
    return specdict


def read_csv_file(input_file):
    """
    Import a csv file with a spectrum, in either purine or pyrimidine notation.
    It will check for incorect entries like invalid mutation/context pair or wrong contexts.
    Any context unspecified is set to 0.
    Returns a OrderedDict spectrum with counts for each entry.
    """

    acons = [rev_comp(con) for con in tcons]
    gcons = [rev_comp(con) for con in ccons]

    specdict = init_spec_dict()  # Default is pyrimidine centric.

    with open(input_file, 'r') as fi:
        lines = fi.readlines()

        for line in lines[1:]:
            oneline = line.strip().split(',')

            # Expecting minimum 3 columns: mutation, context, count and proportion
            if len(oneline)>4:
                oneline=oneline[0:4]

            mut, con, count, prop = oneline

            if con == "Context":
                continue
            
            # Check if py or pu contexts
            if con in acons+gcons:  
                #Found a purine contexts
                conval = rev_comp(con)
                if mut in pu_muts:
                    mutval = rev_comp(mut[0])+'>'+rev_comp(mut[2])
                else:
                    print('Invalid mutation+context pair {} {}'.format(mut, con))
                    continue
            elif con in ccons+tcons:
                #Found a pyrimidine context
                conval=con
                if mut in py_muts:
                    mutval = mut
                else:
                    print('Invalid mutation+context pair {} {}'.format(mut, con))
                    continue
            else:
                print('Invalid sequence context {}'.format(con))
                continue

            #print('Found mutation {} in context {}, with value {}'.format(mutval, conval, count))
            specdict[(mutval,conval)] = specdict[(mutval,conval)] + float(count)           

    return specdict

def read_msp_file(input_file):
    """
    Import a msp file with a spectrum (which is a csv file with a header, in either purine or pyrimidine notation.
    It will check for incorect entries like invalid mutation/context pair or wrong contexts.
    Any context unspecified is set to 0.
    Returns a OrderedDict spectrum with counts for each entry.
    """

    acons = [rev_comp(con) for con in tcons]
    gcons = [rev_comp(con) for con in ccons]

    specdict = init_spec_dict()  # Default is pyrimidine centric.

    with open(input_file, 'r') as fi:
        lines = fi.readlines()

        # Additional code here for msp files. Because of header, they start on line 8
        # Later, this will be updated to read a variable length header....

        for line in lines[8:]:
            oneline = line.strip().split(',')

            # Expecting minimum 3 columns: mutation, context, count and proportion
            if len(oneline)>4:
                oneline=oneline[0:4]

            mut, con, count, norm = oneline  # count is not normalized; norm is

            if con == "Context":
                continue
            
            # Check if py or pu contexts
            if con in acons+gcons:  
                #Found a purine contexts
                conval = rev_comp(con)
                if mut in pu_muts:
                    mutval = rev_comp(mut[0])+'>'+rev_comp(mut[2])
                else:
                    print('Invalid mutation+context pair {} {}'.format(mut, con))
                    continue
            elif con in ccons+tcons:
                #Found a pyrimidine context
                conval=con
                if mut in py_muts:
                    mutval = mut
                else:
                    print('Invalid mutation+context pair {} {}'.format(mut, con))
                    continue
            else:
                print('Invalid sequence context {}'.format(con))
                continue

            #print('Found mutation {} in context {}, with value {}'.format(mutval, conval, count))
            specdict[(mutval,conval)] = specdict[(mutval,conval)] + float(norm)           

    return specdict



def unit_norm(spec):
    """
    Unit normalize a spectrum supplied as a dictionary of counts (pyrimidine style).
    Returns a new dictionary with float values.
    """

    newspec = init_spec_dict()

    if type(list(spec.values())[0])is tuple:
        # Code for avg, std type dictionary.
        #print('Found tuple')
        total = sum(val[0] for val in list(spec.values()))
        for key in newspec.keys():
            newspec[key]=(float(spec[key][0]/total), float(spec[key][1]/total))
    else:
        # Code for count type dictionary
        #print('No tuple')
        total = sum(spec.values())
        for key in newspec.keys():
            newspec[key]=float(spec[key]/total)

    return newspec

def normalize_spec(spec, contexts):
    """
    Takes a spectrum dictionary spec (py-centric) and contexts counts (py-centric).
    Normalizes to the contexts counts and then unit normalizes the result.
    Returns a unit-normalized spectrum.

    Spec dictionary contains tuple (mut, context) keys. Context dict contains only context keys.
    """

    newspec = init_spec_dict()

    if type(list(spec.values())[0])is tuple:
        # Code for avg, std type dictionary.
        for key in newspec.keys():
            newspec[key] = (float(spec[key][0]/contexts[key[1]]), float(spec[key][1]/contexts[key[1]]))
    else:
        # Code for count type dictionary.
        for key in newspec.keys():
            newspec[key] = float(spec[key])/float(contexts[key[1]])

    return(unit_norm(newspec))

def combine_csv_files(files_list, op='avg'):
    """
    # Read each file in the list of csv files.
    # Returns a new spectrum with the sum/avg of the input spectra.
    # Op - operation is 'sum' - returns sum of counts
    #                   'avg' - unitnormalizes inputs, takes avg and std. Returns tuples (avg, std).
    """
    
    dictlist=[]
    for file in files_list:
        dictlist.append(read_csv_file(file))

    #print(dictlist)
    tot = len(dictlist)
    newspec = init_spec_dict()

    if op=='sum':
        # sum of counts only
        for key in newspec.keys():
            newspec[key]=sum([dict[key] for dict in dictlist])
    else:
        #assume op is avg+std
        # Unit norm dicts
        dictlistnorm = [unit_norm(dict) for dict in dictlist]

        for key in newspec.keys():
            avg = sum([dict[key] for dict in dictlistnorm])/tot
            std = stdev([dict[key] for dict in dictlistnorm])
            newspec[key]= (avg, std)
    
    return newspec

def save_csv_file(filename, spec):
    """
    Save a spec dictionary to a csv file.
    If spec has counts values, it saves Mutation, Context, Counts columns
    If spec has tuple values (avg, std), it saves Mutation, Context, Average, Stdev columns
    """

    if type(list(spec.values())[0])is tuple:
        # Tuples saving
        with open(filename, 'w') as fo:
            fo.write('Mutation, Context, Average, Stdev, \n')
            for (mut,con),(avg,std) in spec.items():
                fo.write(','.join([mut, con, str(avg), str(std), '\n']))
    else:
        # Counts saving
        with open(filename, 'w') as fo:
            fo.write('Mutation, Context, Counts, \n')
            for (mut,con), count in spec.items():
                fo.write(','.join([mut, con, str(count), '\n']))

    return

def subtract_background(spec, bgrspec, **kwargs):
    """
    Subtract background from a spec, using bgrspec. Spec and bgrspec are spec OrderedDicts.
    Kwargs specify additional modes of subtraction.
    Default mode is 'direct', which involves direct subtraction of the counts for each context.
    Another mode is 'weighted', which requires 'ratio' to be set to the ratio of the bgrspec in spec.
    """

    # Determining the mode of subtraction

    mode = kwargs.pop('mode', 'direct')  # default is 'direct'

    if mode=='weighted':
        #Code for weighted subtraction

        print('Weighted mode')

        return spec

    if mode=='direct':
        # Code for direct subtraction

        print('Direct mode')

        newspec = init_spec_dict()
        negvalues=0

        for key in newspec.keys():
            value = spec[key]-bgrspec[key]
            print('At {}, value is {}'.format(key, value))

            negvalues+=1 if value<0 else 0
            newspec[key]=max(value, 0)
            
        print('Total negative values: {}'.format(negvalues))    
        return newspec

    else:
        print("Invalid mode selected. Input is unmodified.")
        return spec



def average2spec(spec1, spec2):
    """
    Average the counts from two spectra. Return the average spec.
    """

    newspec = init_spec_dict()

    for key in newspec.keys():

        value = (spec1[key]+spec2[key])/2
        newspec[key] = value
        print('{} {} {} avg {}'.format(key, spec1[key], spec2[key], value))

    return newspec



        




###################
#Plotting functions
###################

def colored_bins(division, ax=None, colors=None, labels=None, padding=0):
    if ax is None:
        ax = plt.gca()

    ax = axes_onoff(spines_onoff(ticks_onoff(init_chart(ax))))

    colors = cycle(['0.8']) if colors is None else cycle(colors)

    for bin in range(division):
        ax.add_patch(Rectangle(xy=(bin / division, 0),
                               width=(1 / division) - padding,
                               height=1,
                               color=next(colors)))

    if labels is not None:
        ax.get_xaxis().set_visible(True)
        xticks = [(x - 0.5) / division for x in range(1, division + 1)]
        ax.set_xticks(xticks)
        ax.set_xticklabels(labels)
    return ax


def spec_barplot(heights, ax=None, xlabels=None, y_max=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    #ax = ticks_onoff(init_chart(ax), True)
    ax = init_chart(ax)
    
    
    #ax.yaxis.grid(True, color=str(0.8), ls='-') #Gridlines

    ax.set_xlim([-0.5, len(heights)])
    if y_max is not None:
        ax.set_ylim(ax.get_ylim()[0], y_max)

    bar_width = kwargs.pop('bar_width', 0.65)

    errorbars = kwargs.pop('errorbars', None)

    colorscheme = kwargs.pop('colorscheme', 'COSMIC2')
    
    if errorbars is None:
        bars = ax.bar(x=range(len(heights)),
                      height=heights,
                      width=bar_width,
                      zorder=3)
    else:
        errkw={'capsize':3, 'ecolor':str(0.5), 'linewidth':1} # dictionary for the error bars parameters

        bars = ax.bar(x=range(len(heights)),
                      height=heights,
                      width=bar_width,
                      zorder=3, 
                      yerr=errorbars, error_kw=errkw)

    for i, color in enumerate([c for c in colormap[colorscheme] for _ in range(16)]):
        bars[i].set_color(color)

    ax.set_xticks([tick - 0.32 + bar_width / 2 for tick in range(len(heights))])
    ax.set_xticklabels(xlabels, family='monospace', rotation=90, fontsize=10)
    ax.set_ylabel(kwargs.pop('ylabel', 'Percent of Mutations'), fontsize=12)
    ax.set_title(kwargs.pop('title', None), y=0.84)
    return ax

def spec_figure(nrow, ncol, heights, xlabels=None, labels=None, y_max=None,
                 titles=None, x_inches=16, ylabel=None, errorbars=None, colorscheme='COSMIC3'):


    aspect = 4 / 11
    fig, axes = plt.subplots(nrow * 2, ncol,
                             figsize=(x_inches * ncol,
                                      x_inches * nrow * aspect),
                             gridspec_kw={'height_ratios': nrow * [28, 1],
                                          'hspace': 0.2,
                                          'wspace': 0.07})

    axes_iter = iter(enumerate(axes.flatten()))

    j=0 # plot index

    for row in range(nrow * 2):
        if row % 2 == 0:
            for i, ax in [next(axes_iter) for _ in range(ncol)]:
                if j >= len(heights) * 2:
                    ax = axes_onoff(ax)
                    continue

                if errorbars is None:
                    err2d=None
                else:
                    #print(errorbars[j])
                    err2d1 = [0 for _ in range(len(errorbars[j]))]
                    err2d2 = errorbars[j]
                    err2d =[err2d1, err2d2]

                print('Current ', j)

                spec_barplot(
                    heights[j],
                    ax=ax,
                    y_max=y_max,
                    ylabel=ylabel,
                    xlabels=None if xlabels is None else xlabels[j],
                    title=None if titles is None else titles[j],
                    errorbars=err2d,
                    colorscheme=colorscheme)
                j+=1
        else:
            for i, ax in [next(axes_iter) for _ in range(ncol)]:
                if i >= (len(heights)*2+1):
                    ax = axes_onoff(ax)
                    continue

                colored_bins(
                    division=6,
                    ax=ax,
                    colors=colormap[colorscheme],
                    labels=labels)
    return fig, axes


##files = ['Dev/testspec1.csv', 'Dev/testspec2.csv', 'Dev/testspec3.csv']
##
##new =combine_csv_files(files, 'avg')


#save_csv_file('testout3.txt', new)
#kmer = import_kmer_counts('Dev/bbmap.count.EG10c.txt')
#new2 = normalize_spec(new, kmer)

###vals = [list(new.values())]
##vals =[list(zip(*new.values()))[0]]
##stds =[list(zip(*new.values()))[1]]
##xlab = [list(zip(*new.keys()))[1]]
###print(xlab)
##mpl.rc("savefig", dpi=320)
###spectrum(vals, xlabels=xlab)
##spec_figure(1, 1, vals, xlabels=xlab, labels=py_muts, errorbars=stds)
##plt.savefig('testfig.png')


##spec1 = read_csv_file('Dev/ndmaspec1.csv')
##bgr = read_csv_file('Dev/testspec1.csv')
##subspec = subtract_background(spec1, bgr, mode='weighted')
##save_csv_file('subspec1.csv', subspec)

##spec1 = read_csv_file('Dev/testspec1.csv')
##spec2 = read_csv_file('Dev/testspec2.csv')
##
##avgspec = average2spec(spec1, spec2)
##save_csv_file('Dev/avgspec.csv', avgspec)
##
##print('spec1 has {} muts'.format(sum(spec1.values())))
##print('spec2 has {} muts'.format(sum(spec2.values())))
##print('avgspec has {} muts'.format(sum(avgspec.values())))
          





def make_figures(data, kmer_counts, sample, format, notation, proportions, ymax=None):
    # Set up local variables
    mpl.rc("savefig", dpi=320)

    image_file1 = sample + '-' + '-freq.' + format
    image_file2 = sample + '-' + '-prop.' + format
    data_file = sample + '.csv'

    total_muts = sum(data.values())
    contexts = context_dict(kmer_counts)
    
    # Divide by trinucleotide context frequencies
    count_dict = {} # contains total counts, unnormalized frequencies
    for key in data.keys():
        sub, con = key
        count_dict[key] = data[key]
        data[key] = data[key]/float(contexts[con])

    # Normalize
    total = sum(data.values())
    for key in data.keys():
        data[key] = data[key]/total

##    # Format and save to CSV
##    with open(data_file, 'w') as fo:
##        fo.write('Substitution,Context,Mutation Count,Normalized Proportion,\n')
##        for (substitution, context), counts in data.items():
##            fo.write(','.join([substitution,context,str(count_dict[(substitution,context)]),str(counts),'\n']))
    
    if ymax is not None:
        ymax = float(ymax)

    # Render plots and save
    if notation == 'pyrimidine':
        labels = sorted(set(list(zip(*data.keys()))[0]))
    if notation == 'purine':
        labels = ['G>T','G>C','G>A','A>T','A>G','A>C']
    if proportions == "proportions":
        spectrum_map(nrow=1, ncol=1,
                     heights=[list(data.values())],
                     xlabels=[list(zip(*data.keys()))[1]],
                     labels=labels,
                     titles=None,
                     ylabel='Proportion of Mutations',
                     y_max=ymax)
        plt.savefig(image_file2)

    if proportions == "frequencies":
        spectrum_map(nrow=1, ncol=1,
                     heights=[list(count_dict.values())],
                     xlabels=[list(zip(*data.keys()))[1]],
#                     labels=sorted(set(list(zip(*data.keys()))[0])),
                     labels=labels,
                     titles=None,
                     ylabel='Frequency of Total Mutations\nNot Normalized',
                     y_max=ymax)
        plt.savefig(image_file1)

def main():
    info = ("PlotSpec - plotting 96bar mutation spectra from a .csv file."
            "By Bogdan Fedeles, April 2021.")

    return

    parser = argparse.ArgumentParser(description=info)
    parser.add_argument("-i", "--input",
                        action = "store",
                        type = str,
                        dest = "input",
                        help = "File paths to input .csv file with spectrum information.",
                        required = False)
    parser.add_argument("-o", "--output_name",
                        action = "store",
                        type = str,
                        dest = "output",
                        help = "Name for the output figures. If no name specified, the input file name will be used",
                        required = False)
    parser.add_argument("-f", "--format",
                        action="store",
                        type=str,
                        dest="format",
                        default='png',
                        help=("Format of the output image file [png]"),
                        required=False)
    parser.add_argument("-y", "--ymax",
                        action="store",
                        dest="ymax",
                        default=None,
                        help="Sets the Y_MAX value in figure.",
                        required=False)
    parser.add_argument("-l", "--labels",
                        action="store",
                        type=str,
                        dest="labels",
                        default="pyrimidine",
                        help=("X-axis mutation labels [purine] or [pyrimidine]"),
                        required=False)
    parser.add_argument("-e", "--errorbars",
                        action="store",
                        dest="errorbars",
                        default=False,
                        help="Plots errorbars on spectrum, read from the .csv file. It assumes the file has an additional column for error bars. ",
                        required=False)
    parser.add_argument("-k", "--kmer_counts",
                        action="store",
                        type=str,
                        dest="kmer_counts",
                        help="The reference kmer count, in a .csv format.",
                        required=False)
    parser.add_argument("-n", "--normalize",
                        type=str,
                        dest="normalize",
                        default="yes",
                        help=("Normalization of the plot. If yes, it requires -k to specify kmer count file. If no, mutational counts are plotted."),
                        required=False)
    
    
    

    args = parser.parse_args()
    inps = args.input.split(',')
    data = extract_mutcsv(args.input, args.labels, args.errorbars) # extract mutational data from .csv file.

    print(str(sum(data.values())) + " total mutations found.")

    make_figures(data, args.kmer_counts, args.output, args.format, args.labels, args.normalize, args.ymax)

if __name__ == '__main__':
    main()


##################
# Code for testing
##################

##dicto=import_kmer_counts('bbmap.count.EG10c.txt', 'pyrimidine')
##
##print(dicto)
##
##with open('output-py.txt', 'w') as f:
##    for key in dicto.keys():
##        f.write(' '.join([key, str(dicto[key]), '\n']))
##
##file = 'EG10c3.fasta'
##
##d2 = fasta_to_dict(file)
##
##print(d2.keys())
##print(len(d2['EG10']))


####
# Things to do
# 1. mutpos_to_csv: modify draw_spectra function to parse mutpos (unique) files and generate .csv files. No plotting. 
# 2. read_csv: read spectrum data (pyrimidine-centric only). DONE!
# 3a. unit_norm: unit normalize a spec. DONE!
# 3b. normalize_spec: function to normalize data given frequencies file. DONE!
# 4. plot normalized spectrum. DONE!
# 5. combine_spectra: read multiple csv's and generate a new .csv with mean and std info. DONE!
# 6. subtract_backgr: read spec and background and do a backgr. subtraction. DONE!
 




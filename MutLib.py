#!/usr/bin/env python3

##MutLib 1.0
##-----------------
##By Bogdan I. Fedeles
##Last update: 2023-08-02
##
##Library of functions to analyze .mutpos files generated by the Loeb DCS pipeline and .mut files generated by TwinStrand pipeline.
##
##List of functions
##
##
##
## Update 2023-08-02. Add table_to_mut function.


import os
import argparse

from Bio import SeqIO
from collections import OrderedDict
from itertools import cycle, product

##Global vars
purines, pyrimidines = ('A', 'G'), ('C', 'T')
dna_bases = sorted(purines + pyrimidines)

py_muts = ('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')
pu_muts = ('G>T', 'G>C', 'G>A', 'A>T', 'A>G', 'A>C')

all_muts = py_muts+pu_muts


##Alexandrov/Stratton color scheme
sig_colors = ['#52C3F1', '#231F20', '#E62223', '#CBC9C8', '#97D54C', '#EDBFC2']


def fasta_to_dict(ref_file):

    with open(ref_file, 'r') as handle:
        return SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))


def from_mutpos(mutpos_file, ref_file, clonality=(0, 1), min_depth=100, kmer=3,
                chromosome=None, start=0, end=None, notation='pyrimidine',
                fmt='essigmann', verbose=False):
    mutations = []
    record_dict = fasta_to_dict(ref_file)

    # Open the mutpos_file and read data line by line
    with open(mutpos_file, 'r') as handle:
        for line in handle:
            # Strip newline characters and split on tabs
            line = line.strip().split('\t')
#            if isinstance(line[0],str):
#                continue

            # Unpack line and cast to proper data type
            chrom, ref = str(line[0]), str(line[1]).upper()
            position, depth = int(line[2]) - 1, int(line[3])

            if chromosome is not None and chrom != chromosome:
                continue
            if position < start:
                continue
            if end is not None and position >= end:
                continue

            if fmt == 'essigmann':
                A, C, G, T, N = map(int, line[4:9])
            elif fmt == 'wesdirect':
                # N = 0
                thisMut = {base: 0 for base in dna_bases}
                thisMut[str(line[5])] = int(line[4])
                A = thisMut['A']
                C = thisMut['C']
                G = thisMut['G']
                T = thisMut['T']
            elif fmt == 'loeb':
                T, C, G, A = map(int, line[5:9])
                # N = int(line[11])
            else:
                raise ValueError('Format must be essigmann, loeb, wesdirect')

            # If we see no observed base substitutions, skip this loop
            if sum([A, C, G, T]) == 0:
                continue

            # Read mapped depth (minus Ns) at this position must be greater
            # than min_depth, if not, skip this loop
            if min_depth > depth:
                continue

            # Get the kmer context at the given position in the given
            # chromosome as looked up in record_dict. If we encounter an
            # edge case (return of None) we'll skip this loop
            context = str(get_kmer(record_dict, chrom, position, kmer)).upper()
            print(context)
            if context is None:
                continue

            # If the base is not in our intended labeling set: let's get
            # the complement of the reference base, the reverse complement
            # of the trinucleotide context, and the mutation counts for the
            # complement of the base substitution observed
            if (
                (notation == 'purine' and ref in pyrimidines) or
                (notation == 'pyrimidine' and ref in purines)
            ):
                ref = str(reverse_complement(ref))
                context = str(reverse_complement(context))
                A, G, C, T = T, C, G, A

            for base, num_mutations in zip(dna_bases, [A, C, G, T]):
                # Clonality is defined as the frequency of any base
                # substitution at this one genomic position. For example, if
                # there are 5% T, 10% C, and 100% G in a reference position of
                # A, we can seletion a clonality filter of (0.1, 0.5) to
                # eliminate the rare T mutations and clonal G mutations.
                base_clonality = num_mutations / depth
                if not min(clonality) <= base_clonality <= max(clonality):
                    continue

                for _ in range(num_mutations):
                    mutation = Mutation(ref, base, chrom, position, context)
                    mutation.depth, mutation.clonality = depth, base_clonality
                    mutations.append(mutation)
                    print(mutation.context)

    if verbose is True:
        print('Found {} Mutations'.format(len(mutations)))
    return mutations


def get_kmer(record_dict, chromosome, position, k=3, pos='mid'):
    """
    Given a dictionary (in memory) of fasta sequences this function will
    return a kmer of length k centered about pos at a certain genomic position
    in an identified dictionary key i.e. chromosome.

    Parameters
    ----------
    record_dict : dict of Bio.Seq objects
        Dictionary where keys are chromosomes and values are Bio.Seq
        sequences.
    chromosome : str
        Key to use when accessing the record_dict.
    position : int
        Position in sequence to query.
    k : int
        Odd number for length of kmer.
    pos : 'mid' or int
        Position in kmer that the genomic position should be centered on.

    Returns
    -------
    kmer : str
        kmer of length k centered around position in chromosome
    """

    if pos == 'mid' and k % 2 != 1:
        raise ValueError("Even length DNA has no midpoint")

    if pos == 'mid':
        pos = int((k - 1) / 2)

    assert k > pos, "Cannot index past length of k"

    start = position - pos
    end = position + (k - pos)

    try:
        chromosome_length = len(str(record_dict[chromosome].seq))
    except ValueError:
        raise ValueError(
            'Chromosome {} not found in reference'.format(chromosome))

    # It is necessary to protect for the two edge cases now
    if start >= 0 and end <= chromosome_length:
        try:
            return str(record_dict[chromosome].seq[start:end])
        except ValueError:
            ValueError(
                'Position {} not found in chromosome {}'.format(
                    start, chromosome))
    else:
        return None

def rev_comp(sequence):
    """
    Returns the reverse complement of the sequence, a DNA sequence string.
    Assumes letters are all caps.
    """
    
    basepair = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    revcomp = ''
    for base in sequence[::-1]:
        revcomp += basepair[base]
    return revcomp



def extract_mutpos_contexts (mutpos_file, ref_file, outfile, mut_type, window):
    # Extract the sequence contexts around a specific type of mutation (mut_type)
    # Window specifies the +/- size around central base (ie. window=1 is trinucleotide, window=2, pentanucleotides, etc)
    # Outfile is a text output.
    # mut_type is a string, in the form C>T, etc.

    ref_dict = fasta_to_dict(ref_file)

    # Check for valid mutation type

    if mut_type not in all_muts:
        raise ValueError('Mut_type is not valid')

    #find complementary mutation
    if mut_type in py_muts:
        comp_mut = pu_muts[py_muts.index(mut_type)]
    else:
        comp_mut = py_muts[pu_muts.index(mut_type)]  
    
        
    
    # Calculate kmer length
    kmer_len = window*2+1
    
    # Open output file
    fo = open(outfile, "w")


    # Open the mutpos_file and read data line by line
    with open(mutpos_file, 'r') as handle:
        for line in handle:
            # Strip newline characters and split on tabs
            line = line.strip().split('\t')
#            if isinstance(line[0],str):
#                continue

            # Unpack line and cast to proper data type
            chrom, ref = str(line[0]), str(line[1]).upper()
            position, depth = int(line[2]) - 1, int(line[3])

##            if chromosome is not None and chrom != chromosome:
##                continue
##            if position < start:
##                continue
##            if end is not None and position >= end:
##                continue

           #format is Loeb
            T, C, G, A = map(int, line[5:9])
            # N = int(line[11])
            
            # If we see no observed base substitutions, skip this loop
            if sum([A, C, G, T]) == 0:
                continue

            if T>0:
                mut = ref+">T"
            if C>0:
                mut = ref+">C"
            if G>0:
                mut = ref+">G"
            if A>0:
                mut = ref+">A"


            context = get_kmer(ref_dict, chrom, position, kmer_len)

            if mut == mut_type:

##              print (mut, " ", context)
##                print(context)
                fo.write(context+"\n")
                continue

            if mut == comp_mut:
                revcontext = rev_comp(context)
##                print(revcontext+"--rev")
                fo.write(revcontext+"\n")


def extract_mut_contexts (mut_file, ref_file, outfile, mut_type, window):
    # Extract the sequence contexts around a specific type of mutation (mut_type)
    # Window specifies the +/- size around central base (ie. window=1 is trinucleotide, window=2, pentanucleotides, etc)
    # Outfile is a text output.
    # mut_type is a string, in the form C>T, etc.
    #
    # Expects .mut file type organization: chr, start, end, sample, var_type, ref, alt, alt_depth, depth, N, subtype, context
    # 
    # The fasta ref file is assumed to have the coordinates of the probes (ie chr1:xxx-yyy), rather than the entire chrom
    # sequence.


    ref_dict = fasta_to_dict(ref_file)

    # extracting key names and start coordinates of probes
    #
    dictkeys = list(ref_dict.keys())

##    print(dictkeys)

    probestartlist=[]
    chrlist=[]

    for key in dictkeys:
        str1, str2 = key.split(':')
        str3 = str2.split('-')[0]
        chrlist.append(str1)
        probestartlist.append(str3)

##    print(chrlist)
##    print (probestartlist)

    

    # Check for valid mutation type

    if mut_type not in all_muts:
        raise ValueError('Mut_type is not valid')

    #find complementary mutation
    if mut_type in py_muts:
        comp_mut = pu_muts[py_muts.index(mut_type)]
    else:
        comp_mut = py_muts[pu_muts.index(mut_type)]  
    
        
    
    # Calculate kmer length
    kmer_len = window*2+1
    
    # Open output file
    fo = open(outfile, "w")


    # Open the mut_file and read data line by line
    with open(mut_file, 'r') as handle:
        for line in handle:
            # Strip newline characters and split on tabs
            line = line.strip().split('\t')
#            if isinstance(line[0],str):
#                continue

            #Check for snv only
            if line[4] != 'snv':
                continue

            # Unpack line and cast to proper data type
            chrom1, ref, alt = str(line[0]), str(line[5]).upper(), str(line[6]).upper()
            pos1, depth = int(line[1]), int(line[8])
            trinuc = line[11].upper()

##            if chromosome is not None and chrom != chromosome:
##                continue
##            if position < start:
##                continue
##            if end is not None and position >= end:
##                continue

            mut = ref+'>'+alt

            # Before calling get_kmer - we need to modify chrom to match the ref keys, and offset position.

            pos = 0 

            # Special case chr1
            if chrom1=='chr1': # could be one of two keys
                if pos1 < 100000000:
                    chrom = dictkeys[0]
                    pos = pos1 - int(probestartlist[0])
                else:
                    chrom = dictkeys[1]
                    pos = pos1 - int(probestartlist[1])
                    
            # it's one of the other chromosomes. 
            elif chrom1 not in chrlist:  # Check if valid
                print('Invalid chromosome label')
                continue
            else:
                ind = chrlist.index(chrom1)
                chrom = dictkeys[ind]
                pos = pos1 - int(probestartlist[ind])
                

            if mut == mut_type:

                context = get_kmer(ref_dict, chrom, pos, kmer_len)
##              print (mut, " ", context)
##                print(context.upper()+"\t"+trinuc)
##                print('{} {} {} {} {} {} {}'.format(mut, pos1, chrom1, pos, chrom, trinuc, context))
                fo.write(context.upper()+"\n")
                continue

            if mut == comp_mut:
                context = get_kmer(ref_dict, chrom, pos, kmer_len)
                revcontext = rev_comp(context.upper())
##                print(revcontext+"--rev\t"+trinuc)
##                print('{} {} {} {} {} {} {} {}'.format(mut, pos1, chrom1, pos, chrom, trinuc, context, revcontext))
                fo.write(revcontext+"\n")

def table_to_mut(table_file, ref_file, outfile):
    # Parse a table of mutations, and select only snvs.
    # Add end pos, sample name, alt_depth(1), depth (100), N (0), subtype column for snvs ("C>T" etc).

    print('Start parsing {}'.format(table_file))

    # open outfile
    fo = open(outfile, "w")
    
    header = ['contig', 'start', 'end', 'sample', 'var_type', 'ref', 'alt', 'alt_depth', 'depth', 'N', 'subtype', 'context', 'filter']
    
    fo.write('\t'.join(header)+'\n')

    
    ref_dict = fasta_to_dict(ref_file)
    # extracting key names and start coordinates of probes
    dictkeys = list(ref_dict.keys())
    print(dictkeys)

    kmer_len = 3 # trincleotides

    #test = 100
    ind = 0
    
    with open(table_file, 'r') as handle:
        for line in handle:
            # Strip newline characters and split on tabs
            line = line.strip().split('\t')
            # expected entries: chr, pos, ref, alt, filter

            if str(line[0])=='CHROM': #header line
                continue

            chrom, ref, alt, filt = str(line[0]), str(line[2]).upper(), str(line[3]).upper(), str(line[4])
            pos = int(line[1])

            # check for non-snvs
            if len(ref)>1:
                continue
            if len(alt)>1:
                continue

            mut = ref+'>'+alt

            context = get_kmer(ref_dict, chrom, pos-1, kmer_len)  # need pos-1 to correct python index 

            #print('{} {} refbase is:{}  found:{} 3mer:{}'.format(chrom, str(pos), ref, context[1], context))

            line_out = [chrom, str(pos), str(pos+1), 'JME', 'snv', ref, alt, str(1), str(100), str(0), mut, context, filt]
            
            fo.write('\t'.join(line_out)+'\n')

            #test+=-1

##            if test==0:
##                fo.close()
##                return
            ind+=1

            if ind % 1000 == 0:
                print('Processed {} mutations'.format(str(ind)))
                

    fo.close()
    print('Parsing done!')


    #def table_to_mut (table_file, ref_file, outfile):
    # Convert a table of mutations file to a .mut file
    # Table contains: Chrom, pos, var_type, ref, alt, subtype, filter
    # Need to add: end, alt_depth (1), depth (100), N (0), context (trinucl)
    #
    # Expects .mut file type organization: chr, start, end, sample, var_type, ref, alt, alt_depth, depth, N, subtype, context
    # 
    # The fasta ref file is assumed to have the coordinates of the probes (ie chr1:xxx-yyy), rather than the entire chrom
    # sequence.


    
def check_mut_integrity(mut_file):

    # Parse a mut file and check that reference base is the same as the middle base of the trinucleotide context
    # Expect .mut file type columns: chr, start, end, sample, var_type, ref, alt, alt_depth, depth, N, subtype, context

    print('Start parsing {}'.format(mut_file))

    ind = 0
    with open(mut_file, 'r') as handle:
        for line in handle:
            # Strip newline characters and split on tabs
            line = line.strip().split('\t')
            # expected entries: chr, start, end, sample, var_type, ref, alt, alt_depth, depth, N, subtype, context

            if line[0] == 'contig': #header line
                continue
            ref = line[5]
            context = line[11]

            if ref!=context[1]:
                print ('Discrepancy found on {}, position {}. Ref:{} Context: {}'.format(line[0], line[1], ref, context))

            ind+=1

            if ind % 1000==0:
               print('Processed {} mutations'.format(str(ind)))
    print('Parsing done!')


############
### MAIN ###
############
         
            
def main():

    #test code

##    ref_gpt = "Dev/EG10_rgc_Corrected.fasta"
##    ref_twnstr = "Dev/mm10_twst-probes50.fa"
##    
##    mut_file = "Dev/8217_mgmt.mut"
##
##    out_path = "Output/"
##
##    extract_mut_contexts(mut_file, ref_twnstr, out_path+"8217_13base_G-A.txt", "G>A", 6) 
    

if __name__ == '__main__':
    main()

#!/usr/bin/env python3

# pLogoSequencesScript
#
# Generate the kmers necessary for the pLogo plotting
#
# Bogdan Fedeles, Nov 2022.
# 5ClC update - June 2023


import MutLib as ml


ref_gpt = "EG10_rgc_Corrected.fasta"
ref_twnstr = "mm10_twst-probes50.fa"

process_mut=True
process_mutpos=True

mutfilespath = "MutFiles/"
name_suff = ".1.consensus.variant-calls.mut"
out_path = "pLogoFiles/"
CT_suff="CtoT_15base_"
TC_suff="TtoC_15base_"

filenames = ['5ClC', 'Control', 'dC']
mutposfiles = ['5ClC-gptDCS', 'Control-gptDCS']


GAmut = "G>A"
AGmut = "A>G"
CTmut = "C>T"
TCmut = "T>C"


interval = 7


if process_mut:
    
    for f in filenames:
    
        mut_file = mutfilespath+f+name_suff
        CT_file = out_path+CT_suff+f+'.txt'
        TC_file = out_path+TC_suff+f+'.txt'

        print('Proceesing file '+mut_file)

        ml.extract_mut_contexts(mut_file, ref_twnstr, CT_file, CTmut, interval) 
        ml.extract_mut_contexts(mut_file, ref_twnstr, TC_file, TCmut, interval)


if process_mutpos:

    for f in mutposfiles:

        mutpos_file = mutfilespath+f+'.mutpos'
        CT_file = out_path+CT_suff+f+'.txt'
        TC_file = out_path+TC_suff+f+'.txt'

        print('Proceesing file '+mutpos_file)

        ml.extract_mutpos_contexts(mutpos_file, ref_gpt, CT_file, CTmut, interval) 
        ml.extract_mutpos_contexts(mutpos_file, ref_gpt, TC_file, TCmut, interval)



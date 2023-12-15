#!/usr/bin/env python3
#
# PlotSpecScript
#
# Analyze data with PlotSpec functions.
#
#

import PlotSpec as ps
import ClustPlot as cp
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import MutScript as ms
import os as os

print(ps.pu_muts)
print(ps.py_muts)

allmuts = ps.pu_muts+ps.py_muts

print(allmuts)


def generate_spectra(group, name, title, kmerfile):

    # Function to plot average mutational spectra.
    # Group is a list of msp files
    # Name is Figure/Save filename
    # Title is a string to be added to the plot
    # kmerfile is a file with trinucleotide proportions (for normalization)


    #### Define here paths to save files:
    datafilepath="Datafiles2/"
    picfilepath="Pics2/"
    ####

    kmer = ps.import_kmer_counts(kmerfile)

    # Printing with purine labels
    pu_dict=ps.init_spec_dict('purine')
    xlab = list(zip(*pu_dict.keys()))[1]
    
    val_list=[]
    std_list=[]
    xlab_list=[]

    print(group)

    new = ps.combine_csv_files(group, op='avg')
    ps.save_csv_file(datafilepath+name[:-5]+'.csv', new)

    #new2=new
    new2 = ps.normalize_spec(new, kmer)
    ps.save_csv_file(datafilepath+name[:-5]+'-norm.csv', new2)

    vals =list(zip(*new2.values()))[0]
    stds =list(zip(*new2.values()))[1]
    #xlab = [list(zip(*new2.keys()))[1]]

    val_list.append(vals)
    std_list.append(stds)
    xlab_list.append(xlab)

    #plotting individual plots
    mpl.rc("savefig", dpi=320)
    ps.spec_figure(1, 1, [vals], xlabels=[xlab], labels=ps.pu_muts, errorbars=[stds],
                   titles=[title], ylabel='Proportion of mutations')
    plt.savefig(picfilepath+name+'.svg')
    plt.savefig(picfilepath+name+'.png')

    return


def plot_msp(mspfile, title, kmerfile):

    # Plot one single msp file
    # The name of the figure will be same as msp file, but with svg and png extension
    # Title is a string for the plot title
    # kmerfile - file with trinucl counts for normalization


    ###Data paths
    datafilepath="Datafiles/"
    picfilepath="Pics/"

    kmer = ps.import_kmer_counts(kmerfile)

    # Printing with purine labels
    #pu_dict=ps.init_spec_dict('purine')
    py_dict=ps.init_spec_dict()
    xlab = list(zip(*py_dict.keys()))[1]

    spec = ps.read_msp_file(datafilepath+mspfile)

    #normalization
    newspec = ps.normalize_spec(spec, kmer)
    ps.save_csv_file(datafilepath+mspfile[:-5]+'-norm.csv', newspec)

    val_list=[]
    #vals =list(zip(*newspec.values()))[0]
    vals = list(newspec.values())
    val_list.append(vals)

    #plotting
    mpl.rc("savefig", dpi=320)
    ps.spec_figure(1, 1, [vals], xlabels=[xlab], labels=ps.py_muts, errorbars=None,
                   titles=[title], ylabel='Proportion of mutations')
    plt.savefig(picfilepath+mspfile[:-5]+'.svg')
    plt.savefig(picfilepath+mspfile[:-5]+'.png')

    return
    

def main():

    filespath = 'Datafiles/'
    #filespath2= 'Datafiles2/'
    picspath = 'Pics/'
    #picspath2= 'Pics2/'
    filesuff1 = 'msp.csv'

    #kmerfile = 'bbmap.count.EG10c.txt'
    kmerfile = 'twnstr-mouse-contexts.txt'

    plotspectra=True

    hist= False
    cossim = False
    cossimhist = False

    generatemsp = False

    cossimcosmic = False

    strbias = False

    
    if plotspectra:

        mspfiles=['5AC-Day2', '5ClC', 'ControlDay2', 'dC']
        ext = 'msp.csv'

        mspfiles2=[file+ext for file in mspfiles]

        for file, title in zip(mspfiles2, mspfiles):
             plot_msp(file, title, kmerfile)

        print ('done!')
             

    if hist:
    # Clustplot of heatmap'

        print('Making hist')

        datasamples = ['NDMA-WT', 'NDMA-Mgmt','Temozolomide', 'Streptozotocin', 'MNU', 'SBS11']

        datafiles = ['all_NDMA_sum_bgsub', 'all-mgmt-ndma-avg', 'Temozolomide_bgsub', 'Streptozotocin_bgsub',
                     'MNU_bgsub', 'SBS11_32']
        datafiles2= [file+'_norm.csv' for file in datafiles]

        pu_dict=ps.init_spec_dict('purine')
        xlab = list(zip(*pu_dict.keys()))[1]

        val_list=[]
        xlab_list=[xlab]
        spec_list={}

        for sam, file in zip(datasamples, datafiles2):
            spec_list[sam] = ps.read_csv_file(filespath2+file)


        ext1='png'
        ext2='svg'
        
        for st_col in [2.3]:
##        for st_col in [2.2, 2.4, 2.6, 2.8]:
            plotfile1 = picspath2+'histplot'+str(st_col)+'.'+ext1
            plotfile2 = picspath2+'histplot'+str(st_col)+'.'+ext2

            cp.plot_uhc_heatmap(spec_list, datasamples, True, plotfile1, ext1, st_col)  
            cp.plot_uhc_heatmap(spec_list, datasamples, True, plotfile2, ext2, st_col)  


    if cossim:

        # Cossim between various spectra

        wt_ctrl = ps.read_csv_file(filespath2+'wt-ctrl-all-norm.csv')
        wt_ctrl_vals = list(wt_ctrl.values())

        mgmt_ctrl = ps.read_csv_file(filespath2+'all-mgmt-ctrl-avg-norm.csv')
        mgmt_ctrl_vals = list(mgmt_ctrl.values())
        
        #for sig, file in zip(csigs,cfiles2):

##        sigspec = ps.read_msp_file(file)
##        sigspecvals = list(sigspec.values())

        cossim = cp.cosine_similarity([wt_ctrl_vals], [mgmt_ctrl_vals])
        print('{}  \t  {}'.format('wt vs mgmt ctrl', cossim[0][0])) 


    if cossimhist:

        # cossim between a collection of spectra

        datasamples = ['WT-ctrl', 'WT-males-ctrl', 'WT-fem-ctrl', 'Mgmt-ctrl', 'Mgmt-males-ctrl', 'Mgmt-fem-ctrl']

        datafiles = ['wt-ctrl-all', 'wt-ctrl-males', 'wt-ctrl-females', 'all-mgmt-ctrl-avg', 'mgmt-ctrl-males-avg',
                     'mgmt-ctrl-females-avg']

        datafiles2 = [file+'-norm.csv' for file in datafiles]

        pu_dict=ps.init_spec_dict('purine')
        xlab = list(zip(*pu_dict.keys()))[1]

        val_list=[]
        xlab_list=[xlab]
        spec_list={}

        for sam, file in zip(datasamples, datafiles2):
            spec_list[sam] = ps.read_csv_file(filespath2+file)


        ext1='png'
        ext2='svg'

        plotfile1 = picspath2+'ctrl-histplot'+str(2.3)+'.'+ext1
        plotfile2 = picspath2+'ctrl-histplot'+str(2.3)+'.'+ext2

        cp.plot_uhc_heatmap(spec_list, datasamples, True, plotfile1, ext1, 2.3)  
        cp.plot_uhc_heatmap(spec_list, datasamples, True, plotfile2, ext2, 2.3)  


    if generatemsp:

        # Generate msp files from .mut files

        mutfiles=['5AC-Day2', '5ClC', 'ControlDay2', 'dC']
        ext = '.mut'

        mutfiles2=[filespath+file+ext for file in mutfiles]

        outfiles=[filespath+file+filesuff1 for file in mutfiles]

        for mutfile, outfile in zip(mutfiles2, outfiles):

            ms.mut_to_msp(mutfile, outfile, maxratio=0.1, mindepth=1000)


        print ('done!')




    if cossimcosmic:

        # comparison with Cosmic sigs

        cosmicpath="C:/Users/bogdan/Dropbox (Personal)/BF_RESEARCH/Code/Data/Cosmic/CosmicV3.1/"

        cosmicfiles=os.listdir(cosmicpath)

        file5ClC="DataFiles/5ClCmsp.csv"

        kmer = ps.import_kmer_counts(kmerfile)


        spec5ClC=ps.read_msp_file(file5ClC)
        spec5ClCnorm=ps.normalize_spec(spec5ClC, kmer)
        
        vals5ClC=list(spec5ClCnorm.values())
        vals5ClC_CT=vals5ClC[32:48]

        print('5ClC spectrum cossim comparison')
         
        for file in cosmicfiles:

            sigspec = ps.read_msp_file(cosmicpath+file)
            sigspecvals = list(sigspec.values())
            sigspecvals_CT=sigspecvals[32:48]
    
        
            cossim = cp.cosine_similarity([vals5ClC], [sigspecvals])
            cossim2 = cp.cosine_similarity([vals5ClC_CT], [sigspecvals_CT])

            
            print('{}{}  \t  {} \t {}'.format('vs ',file, cossim[0][0], cossim2[0][0])) 



    if strbias:

        # Strand bias calculation

        mutfile='DataFIles/5ClC.mut'

        num = list(range(1,20))
        chrlist = ['chr'+str(r) for r in num]

        print('chromosome', 'total', 'type', 'strand+', 'strand-', 'bias', sep='\t')

        for c in chrlist:
            ms.mut_count(mutfile, verbose=False, selection=[c])



        return
        

if __name__=="__main__":
    main()


#!/usr/bin/env python3
#
# PlotSpecScript Cell Culture 5ClC
#
# Analyze 5ClC Cell Culture data with PlotSpec functions.
#
# Bogdan Fedeles, June 2023.

import os
import PlotSpec as ps
import ClustPlot as cp
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

print(ps.pu_muts)

def main():
    #print("Hello there!")

    #setup list of custom programs
    
    proglist = ['Read COSMIC', 'Do 5ClC plots (4)', 'Make hist (4)', 'Do 5ClCgpt plot', 'Do 5ClC plots (5)', 'Make hist (5)']

    proglist2 =[str(idx+1)+'. '+prog for idx,prog in enumerate(proglist)]
    
    prompt ='Select program:\n'+'\n'.join(proglist2)+'\n'
    a=input(prompt)

    print ('You have selected: '+a)

    choices=[str(i+1) for i in range(len(proglist))]

    prog=0

    while True:
        if a not in choices:
            print('Invalid selection. Try again!')
            print('Choices are: '+' '.join(choices))
            a = input()
        else:
            print('You have selected '+proglist[int(a)-1])
            prog=int(a)
            break

    print('Ready to run')
    print(prog)

    filespath = 'DataFiles/'
    picspath = 'Pics/'
    filesuff = '_FINAL2_unique2.csv'
    filesuff2 = '_FINAL3_unique2.csv'

##    kmerfile = 'bbmap.count.EG10c.txt'
##    kmer = ps.import_kmer_counts(kmerfile)

    sample = ['S9_control', 'NDMA1', 'NDMA2', 'Vehicle_control', 'Temozolomide',
               'Untreated_control', 'Streptozotocin', 'MNU', 'Aflatoxin_B1', 'Sterigmatocystin']


    datafiledict={}

    for i,sam in enumerate(sample):
        datafiledict[sam]=filespath+str(i+48)+filesuff2

##    print(datafiledict)


    if prog==1:
        # read a cosmic sig, then save as a regular spec file.

        #sigfile = "v3.2_SBS11_PROFILE.txt"
        #sigfile = "v3.3_SBS84_PROFILE.txt"
        sigfile = "v3.3_SBS42_PROFILE.txt"


        spec = ps.init_spec_dict()

        with open(filespath+sigfile, 'r') as f:
            lines = f.readlines()  # read all content
            topline = lines[0].strip().split(',')

        for l in lines[1:]:
            linelist = l.strip().split('\t')

            # first col is mut and context 'A[C>A]A'
            # second col is the frequency
            first=linelist[0]

            mut=first[2:5]
            con=first[0]+first[2]+first[6]

            print('mutation {} context {} freq {}'.format(mut, con, linelist[1]))
            
            spec[(mut, con)]=linelist[1]

        
        ps.save_csv_file(filespath+'SBS42_33.csv', spec)
        
        humkmerfile="hg38_ref_counts.txt"
        hkmer = ps.import_kmer_counts(humkmerfile)

        print(hkmer)

        normspec = ps.normalize_spec(spec, hkmer)
        ps.save_csv_file(filespath+'SBS42_33_norm.csv', normspec)


    if prog==2:

        # plotting normalized spectra of 5ClC, dC, SBS84, SBS42

        datasamples=['5ClCms-norm', 'dCms-norm', 'SBS84_33_norm', 'SBS42_33_norm']

        pu_dict=ps.init_spec_dict('pyrimidine')
        xlab = list(zip(*pu_dict.keys()))[1]

        val_list=[]
        xlab_list=[xlab]

        for sam in datasamples:
            spec=ps.unit_norm(ps.read_csv_file(filespath+sam+'.csv'))
    

##            if sam=='SBS11_32_norm':
##                spec=ps.unit_norm(ps.read_csv_file(filespath+sam+'.csv'))
##            else:
##                spec=ps.read_csv_file(filespath+sam+'_bgsub_norm.csv')


            vals= list(spec.values())
            val_list.append(vals)
            xlab_list.append(xlab)

            #plotting individual plots
##            mpl.rc("savefig", dpi=320)
##            ps.spec_figure(1, 1, [vals], xlabels=[xlab], labels=ps.pu_muts,
##                           titles=[sam], ylabel='Proportion of mutations')
##            plt.savefig(picspath+sam+'_bgsub_norm.png')

        #plotting composite figs
        mpl.rc("savefig", dpi=600)
        ps.spec_figure(4, 1, val_list, xlabels=xlab_list, labels=ps.pu_muts,
                       titles=datasamples, ylabel='Proportion of mutations', colorscheme='COSMIC3')

        plt.savefig(picspath+'5ClC-specs_4x1.svg')
        plt.savefig(picspath+'5ClC-specs_4x1.png')
    
    

    
    if prog==3:
        # Clustplot of heatmap'

        print('Making hist')

        datasamples=['5CldC', 'dC', 'SBS84', 'SBS42']

        datafiles = ['5ClCms-norm', 'dCms-norm', 'SBS84_33_norm', 'SBS42_33_norm']
        datafiles2= [file+'.csv' for file in datafiles]

        pu_dict=ps.init_spec_dict('pyrimidine')
        xlab = list(zip(*pu_dict.keys()))[1]

        val_list=[]
        xlab_list=[xlab]
        spec_list={}

        for sam, file in zip(datasamples, datafiles2):

            spec_list[sam] = ps.read_csv_file(filespath+file)


        ext1='png'
        ext2='svg'
        
        for st_col in [2.3]:
##        for st_col in [2.2, 2.4, 2.6, 2.8]:
            plotfile1 = picspath+'histplot'+str(st_col)+'.'+ext1
            plotfile2 = picspath+'histplot'+str(st_col)+'.'+ext2

            cp.plot_uhc_heatmap(spec_list, datasamples, True, plotfile1, ext1, st_col)  
            cp.plot_uhc_heatmap(spec_list, datasamples, True, plotfile2, ext2, st_col)  



    if prog==4:
        # import the 5ClC-gpt csv file
        # normalize, save and plot

        kmerfile = 'bbmap.count.EG10c.txt'
        kmer = ps.import_kmer_counts(kmerfile)

        py_dict=ps.init_spec_dict()
        xlab = list(zip(*py_dict.keys()))[1]

        csvfile='5ClC-gptDCS.csv'

        spec = ps.read_csv_file(filespath+csvfile)

        #normalization
        newspec = ps.normalize_spec(spec, kmer)
        ps.save_csv_file(filespath+csvfile[:-4]+'-norm.csv', newspec)

        val_list=[]
        #vals =list(zip(*newspec.values()))[0]
        vals = list(newspec.values())
        val_list.append(vals)

        #plotting
        mpl.rc("savefig", dpi=320)
        ps.spec_figure(1, 1, [vals], xlabels=[xlab], labels=ps.py_muts, errorbars=None,
                       titles=['5ClC-gptDCS'], ylabel='Proportion of mutations')
        plt.savefig(picspath+csvfile[:-4]+'.svg')
        plt.savefig(picspath+csvfile[:-4]+'.png')


    if prog==5:
        # plotting normalized spectra of 5ClC, 5ClCgpt, dC, SBS84, SBS42

        datasamples=['5ClCms-norm', '5ClC-gptDCS-norm', 'dCms-norm', 'SBS84_33_norm', 'SBS42_33_norm']

        pu_dict=ps.init_spec_dict('pyrimidine')
        xlab = list(zip(*pu_dict.keys()))[1]

        val_list=[]
        xlab_list=[xlab]

        for sam in datasamples:
            spec=ps.unit_norm(ps.read_csv_file(filespath+sam+'.csv'))

            vals= list(spec.values())
            val_list.append(vals)
            xlab_list.append(xlab)

        #plotting composite figs
        mpl.rc("savefig", dpi=600)
        ps.spec_figure(5, 1, val_list, xlabels=xlab_list, labels=ps.pu_muts,
                       titles=datasamples, ylabel='Proportion of mutations', colorscheme='COSMIC3')

        plt.savefig(picspath+'5ClC-specs_5x1.svg')
        plt.savefig(picspath+'5ClC-specs_5x1.png')
        

    if prog==6:
    
        # Clustplot of heatmap, 5 inputs'

        print('Making hist')

        datasamples=['5CldC', '5CldC-gpt', 'dC', 'SBS84', 'SBS42']

        datafiles = ['5ClCms-norm', '5ClC-gptDCS-norm', 'dCms-norm', 'SBS84_33_norm', 'SBS42_33_norm']
        datafiles2= [file+'.csv' for file in datafiles]

        pu_dict=ps.init_spec_dict('pyrimidine')
        xlab = list(zip(*pu_dict.keys()))[1]

        val_list=[]
        xlab_list=[xlab]
        spec_list={}

        for sam, file in zip(datasamples, datafiles2):

            spec_list[sam] = ps.read_csv_file(filespath+file)


        ext1='png'
        ext2='svg'
        
        for st_col in [2.3]:
##        for st_col in [2.2, 2.4, 2.6, 2.8]:
            plotfile1 = picspath+'histplot'+str(st_col)+'.'+ext1
            plotfile2 = picspath+'histplot'+str(st_col)+'.'+ext2

            cp.plot_uhc_heatmap(spec_list, datasamples, True, plotfile1, ext1, st_col)  
            cp.plot_uhc_heatmap(spec_list, datasamples, True, plotfile2, ext2, st_col)  

        

    if prog==9:

        # plotting before and after bkgr subtraction

        # Printing with purine labels
        pu_dict=ps.init_spec_dict('purine')
        xlab = list(zip(*pu_dict.keys()))[1]

        val_list=[]
        xlab_list=[xlab]
        xlab_list.append(xlab)

        for sam in datasamples:

            spec1=ps.read_csv_file(filespath+sam+'.csv')
            spec2=ps.read_csv_file(filespath+sam+'_bgsub.csv')

            vals1=list(spec1.values())
            vals2=list(spec2.values())

            val_list=[vals1, vals2]

            mpl.rc("savefig", dpi=320)
            ps.spec_figure(2, 1, val_list, xlabels=xlab_list, labels=ps.pu_muts,
                            ylabel='Mutation counts')
            plt.savefig(picspath+sam+'.png')
    
    if prog==7:
        
        # plotting normalized spectra with SBS11 and NDMA all

##        datasamples = ['NDMA_all', 'Temozolomide', 'Streptozotocin', 'MNU', 'SBS11']
        datasamples = ['Temozolomide', 'MNU', 'SBS11', 'NDMA_all', 'Streptozotocin']


##        datafiles = ['all_NDMA_sum_bgsub', 'Temozolomide_bgsub', 'Streptozotocin_bgsub',
##                     'MNU_bgsub', 'SBS11_32']

        datafiles = ['Temozolomide_bgsub', 'MNU_bgsub', 'SBS11_32', 'all_NDMA_sum_bgsub', 'Streptozotocin_bgsub']
        
        datafiles2= [file+'_norm.csv' for file in datafiles]

        pu_dict=ps.init_spec_dict('purine')
        xlab = list(zip(*pu_dict.keys()))[1]

        val_list=[]
        xlab_list=[xlab]

        for file in datafiles2:

            spec = ps.read_csv_file(filespath+file)
            vals= list(spec.values())
            val_list.append(vals)
            xlab_list.append(xlab)
        
        #plotting composite figs
        mpl.rc("savefig", dpi=320)
        ps.spec_figure(5, 1, val_list, xlabels=xlab_list, labels=ps.pu_muts,
                       titles=datasamples, ylabel='Proportion of mutations')

        plt.savefig(picspath+'Alkyl-agents_5x1_order2.svg')
        plt.savefig(picspath+'Alkyl-agents_5x1_order2.png')



    if prog==8:
        # average controls
        spec1 = ps.read_csv_file(datafiledict['Vehicle_control'])
        spec2 = ps.read_csv_file(datafiledict['Untreated_control'])

        avgspec = ps.average2spec(spec1, spec2)
        ps.save_csv_file(filespath+'avg_control3.csv', avgspec)
	


        for sam in datasamples:
            spec = ps.read_csv_file(datafiledict[sam])
            subspec = ps.subtract_background(spec, avgspec)

            print('spec has {} muts'.format(sum(spec.values())))
            print('bgr has {} muts'.format(sum(avgspec.values())))
            print('subspec has {} muts'.format(sum(subspec.values())))
            
            ps.save_csv_file(filespath+sam+'3.csv', spec)
            ps.save_csv_file(filespath+sam+'3_bgsub.csv', subspec)
        


        
    

    if prog==9:
        # Cosine similarity of NDMA spec with all COSMIC v3.1 sigs

        print('Cossim with COSMIC v3.1')

        cpath = 'C:/Users/bogdan/Dropbox (Personal)/BF_RESEARCH/Code/Data/Cosmic/CosmicV3.1/'
        cfiles = os.listdir(cpath)
        cfiles2=[cpath+file for file in cfiles]
        csigs=[file[20:-4] for file in cfiles]
        #print(csigs)

        ndmafile = filespath+'all_NDMA_sum_bgsub_norm.csv'
        #print(ndmafile)

        ndmaspec = ps.read_csv_file(ndmafile)
        ndmaspecvals = list(ndmaspec.values())
        
        for sig, file in zip(csigs,cfiles2):

            sigspec = ps.read_msp_file(file)
            sigspecvals = list(sigspec.values())

            cossim = cp.cosine_similarity([ndmaspecvals], [sigspecvals])
            print('{}  \t  {}'.format(sig, cossim[0][0]))


if __name__=="__main__":
    main()


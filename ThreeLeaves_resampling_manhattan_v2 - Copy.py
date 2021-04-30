# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 13:35:50 2021

@author: Michael Tross
"""

import matplotlib.pyplot as plt

total_reps = 100
chr_length_file = "chr_lengths.csv"
candidate_file = "Candidate_genes.csv"
nam_qtl_file = "NAM_CandidateGenes.csv"
mycolors = ['Red','Green','Blue','Purple']
mycolors = mycolors + mycolors + mycolors


result_files = ['MLM_Leaf1sigSNPCount.csv','MLM_Leaf2sigSNPCount.csv','MLM_Leaf3sigSNPCount.csv']
fig = plt.figure(figsize=(20,15)) 

plotCount = 1
for file in result_files: 
    resampling_result_file = file
    
    # Load chromosome length file
    fh = open(chr_length_file)
    chr_lengths = {}
    for x in fh:
        y = x.strip().split(',')
        chr_lengths[y[0]] = int(y[1])
    
    mychrs = list(chr_lengths)
    mychrs.sort(key=lambda a:int(a[3:]))
    additions = {}
    count = 0
    xticks = []
    for achr in mychrs:
        additions[achr] = count
        xticks.append(count + chr_lengths[achr]/2)
        count += chr_lengths[achr]
    
    
    #load GWAS results
    fh = open(resampling_result_file)
    fh.readline()
    chrsnps = {}
    maxscore = 0
    for achr in additions:
        chrsnps[achr] = {"locs":[],"scores":[]}
    for x in fh:
        y = x.strip().replace('"','').split(',')
        z = y[1].split('_')
        mychr = "Chr" + z[0][1:]
        myloc = int(z[1])
        myscore = int(y[2])/float(total_reps)
        myplot_loc = myloc + additions[mychr]
        chrsnps[mychr]["locs"].append(myplot_loc)
        chrsnps[mychr]["scores"].append(myscore)
        if myscore > maxscore:
            maxscore = myscore
    
    #load candidate genes
    fh = open(candidate_file)
    fh.readline()
    gene_names = []
    gene_pos = []
    gene_jitter = []
    myjitter = 0
    for x in fh:
        y = x.strip().split(',')
        myname = y[0]
        mychr = "Chr" + y[1].zfill(2)
        z = y[2].split("-")
        mypos = (int(z[0]) + int(z[1]))/2
        mypos_plot = mypos + additions[mychr]
    
        gene_names.append(myname)
        gene_pos.append(mypos_plot)
        gene_jitter.append(-0.04)
    
    
    #load NAM QTLS
    fh = open(nam_qtl_file)
    fh.readline()
    qtl_names = []
    qtl_xpos = []
    qtl_ypos = []
   
    for x in fh:
        y = x.strip().split(',')
        qtl_name = y[0]
        chrNum = y[1].split('_')[0].split('S')[1]
        qtl_mychr = "Chr" + chrNum.zfill(2)
        qtl_mypos = int(y[1].split('_')[1])
        qtl_mypos_plot = qtl_mypos + additions[qtl_mychr]
    
        qtl_names.append(qtl_name)
        qtl_xpos.append(qtl_mypos_plot)
        qtl_ypos.append(-0.01)
        
    
    # Add subplot
    ax = fig.add_subplot('31{}'.format(plotCount))
  
    # Plot SNPs according to significance
    for achr,acolor in zip(mychrs,mycolors[:len(mychrs)]):
        ax.scatter(chrsnps[achr]['locs'],chrsnps[achr]['scores'],c=acolor,rasterized=True)
    
    # Plot threshold level
    ax.plot([0,count],[0.1,0.1],'--',c="k")
    
    # Plot QTL positions
    ax.scatter(qtl_xpos,qtl_ypos, color='gold',marker='^',rasterized=True,s=50.0)
    
    # Plot specific candidate genes on each plot
    for i, txt in enumerate(gene_names):
        if plotCount == 1:
            if not txt == 'ZmTAC1': continue
            ax.scatter(gene_pos[i],gene_jitter[i],color="k",marker='^',rasterized=True,s=50.0)
            ax.annotate(txt, (gene_pos[i], gene_jitter[i]-0.05),ha='center',rasterized=True,fontsize=13)
            
        if plotCount == 2:
            if txt == 'Dw3' or txt == 'ZmRAVL1':
                ax.scatter(gene_pos[i],gene_jitter[i],color="k",marker='^',rasterized=True,s=50.0)
                ax.annotate(txt, (gene_pos[i], gene_jitter[i]-0.05),ha='center',rasterized=True,fontsize=13)
        
        if plotCount == 3:
            if txt == 'Dw3':
                ax.scatter(gene_pos[i],gene_jitter[i],color="k",marker='^',rasterized=True,s=50.0)
                ax.annotate(txt, (gene_pos[i], gene_jitter[i]-0.05),ha='center',rasterized=True,fontsize=13)
    
    # Format subplots    
    ax.set_ylabel("Bootstrap Support",size=15)
    ax.set_title("Leaf {} Angle".format(plotCount),size=15)
    ax.set_xticks(xticks)
    ax.set_xticklabels(mychrs)
    ax.set_xlim([0,count])
    ax.set_ylim([-.35,maxscore * 1.1])
    ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    
    plotCount+=1
   
plt.show()
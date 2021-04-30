# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 13:35:50 2021

@author: Michael Tross
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


#configuration
total_reps = 100
chr_length_file = "chr_lengths.csv"
resampling_result_file = "FarmCPU_resamp_result.csv"
gemma_result_file = 'MLM_GWAS_result.csv'
candidate_file = "Candidate_genes.csv"
mycolors = ['Red','Green','Blue','Purple']
nam_qtl_file = "NAM_QTLs.csv"
mycolors = mycolors + mycolors + mycolors



fig = plt.figure(figsize=(20,10))


############################################################################
###################### Plot  GEMMMA  Result ################################
############################################################################



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
    count += chr_lengths[achr] + 10000000 


# load result file
fh = open(gemma_result_file)
fh.readline() # skip headers
chrsnps = {}
maxscore = 0
for achr in additions: 
    chrsnps[achr] = {"locs":[],"scores":[]} 
for x in fh:
    y = x.strip().split(',')
    if y[2] == 'nan' or y[2] == '': continue 
    z = y[0].split('_') 
    mychr = "Chr" + z[0][1:] 
    myloc = int(z[1]) 
    myscore = float(y[2])
    if myscore < 0:
        myscore = 0.0
    myplot_loc = myloc + additions[mychr] 
    chrsnps[mychr]["locs"].append(myplot_loc) 
    chrsnps[mychr]["scores"].append(myscore) 
    if myscore > maxscore: 
        maxscore = myscore
            
# Plot in subpanel 1 of figure
ax = fig.add_subplot('211')


# -log10 bonferroni threshold value
bthreshold = np.negative(np.log10(0.05/78251.75))


# Plot SNPs according to significance
for achr,acolor in zip(mychrs,mycolors[:len(mychrs)]): 
    ax.scatter(chrsnps[achr]['locs'],chrsnps[achr]['scores'],c=acolor,rasterized=True) 

# Plot significance cutoff/threshold
ax.plot([0,count],[bthreshold,bthreshold],'--',c="k") 


# Format figure
ax.set_ylabel("$-log_{10}$(p)",size=15)
ax.set_xticks(xticks) 
ax.set_xticklabels(mychrs, size=15) 
ax.set_xlim([-5000000,count]) 
ax.set_ylim([-12,maxscore * 1.1]) 
ax.set_yticks([0,5,10,15]) 
ax.spines['bottom'].set_position(('data', -0.5)) 
ax.spines['left'].set_position(('data', -5000000)) 
ax.spines['top'].set_visible(False) 
ax.spines['right'].set_visible(False) 
ax.text(-0.1, 1.0, '(A)', transform=ax.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')



################################################################################
###################### Plot FarmCPU Bootstrap Result ###########################
################################################################################


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
    myjitter = -1 * int(y[3])
    gene_jitter.append(myjitter/10.0)   
    

# Add subplot
ax2 = fig.add_subplot('212')

# Plot SNPs according to significance
for achr,acolor in zip(mychrs,mycolors[:len(mychrs)]):
    ax2.scatter(chrsnps[achr]['locs'],chrsnps[achr]['scores'],c=acolor)
    
# Plot significance cutoff/threshold
ax2.plot([0,count],[0.1,0.1],'--',c="k")

# Plot locations of QTLs
ax2.scatter(gene_pos,gene_jitter,color="k",marker='^')
ax2.scatter(qtl_xpos,qtl_ypos, color='gold',marker='^',rasterized=True) 

# Plot gene names
for i, txt in enumerate(gene_names):
    ax2.annotate(txt, (gene_pos[i], gene_jitter[i]-.05),ha='center')


# Format figure
ax2.set_ylabel("Bootstrap Support",size=15)
ax2.set_xticks(xticks)
ax2.set_xticklabels(mychrs)
ax2.set_xlim([0,count])
ax2.set_ylim([-.35,maxscore * 1.1])
ax2.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
ax2.spines['bottom'].set_position(('data', 0))
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.text(-0.1, 1.0, '(B)', transform=ax2.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')


plt.tight_layout()
plt.show()
    

    
    
    
    
    
    
    
    
    
    

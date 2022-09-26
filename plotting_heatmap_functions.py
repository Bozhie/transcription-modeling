import pandas as pd
import numpy as np
import bbi
from gffutils.helpers import asinterval
from gtfparse import read_gtf
import bioframe as bf
import matplotlib.pyplot as plt
import matplotlib.colors as colors


# Uses pybbi to measure signal over set of input
def generate_signal_matrix(chip_seq_file, interval_chrom, interval_start, interval_end, windowSize, windowType, nbins):
    
    if windowType == 'extend':
        interval_start = interval_start - windowSize
        interval_end = interval_end + windowSize
        
    if windowType == 'centered':
        interval_start = (interval_start + interval_end)/2 - windowSize
        interval_end = (interval_start + interval_end)/2 + windowSize
    
    with bbi.open(chip_seq_file) as f:
        matrix = f.stackup(interval_chrom, interval_start, interval_end, bins=nbins)
        
    return matrix
    
def plot_avg_signal(DE_results_df, stackup_matrix, title, ax=None, DE_value_col='log2Fold_Change', cutoff='padj', cutoff_val=0.05, 
                    agg_key='DE_status', agg_categories=['up', 'down', 'nonsig'], color_categories=['r', 'b', 'k'], windowSize=1000, nbins=40):
    
    if ax == None:
        ax = plt.subplot()
    
    for category, color in zip(agg_categories, color_categories):
        
        cat_ix = np.where(DE_results_df[agg_key] == category)
        cat_matrix = stackup_matrix[cat_ix]
        
        ax.plot(np.nanmean(cat_matrix, axis=0), color = color, label=category)
        
    ax.set(xticks=np.arange(0, nbins+1, 10),
    xticklabels=(np.arange(0, (windowSize*2)+1, (windowSize/2))-(windowSize*2)//2),
    xlabel='Distance from boundary (bp)',
    ylabel='ChIP-Seq mean fold change over input')
    ax.set_title(title)


# Note: todo: DE_value_col='log2FoldChange', sort_by_DE=True ==> collect into sort_by_col = 'log2FoldChange'

## todo: do heatmap for DE_values in a separate function (like originally outlined) ** this can be less flexible/specific to diverging up_vs_down and maybe other color/heatmaps can have their own functions for other types of categories
def plot_binned_signal_heatmap(DE_results_df, stackup_matrix, title, ax=None, DE_value_col='log2FoldChange', sort_by_DE=True, 
                               agg_key='DE_status', agg_categories=['up', 'down'], include_category_map=True, color_categories=['r', 'b'], 
                               windowSize=1000, nbins=40):
    
    # if include_category_map --> change dimensions here
    if ax == None:
        if include_category_map:
            fig,ax=plt.subplots(1,2)
        else:
            ax = plt.subplot()
        
        
    # build the heatmap matrices 
    ordered_heatmap = []
    ordered_values = []
        
    # rearrange the matrix for plotting everything together
    for cat in agg_categories:
        
        # get first category
        cat_ix = np.where(DE_results_df[agg_key] == cat)
        sub_matrix = stackup_matrix[cat_ix]
        sub_results = DE_results_df.iloc[cat_ix]
        
        # if sort_by_DE
        # sort according to DE_value in descending order
        ordering = (-sub_results[DE_value_col]).argsort()
        ordered_heatmap.append(sub_matrix[ordering])
        ordered_values.append(sub_results.iloc[ordering][DE_value_col])

        
        
    #logFPKM = np.transpose(np.expand_dims(sig_DE['log2FoldChange'], axis=0))
    #ax.imshow(np.vstack(ordered_heatmap), cmap='gray_r', aspect='auto', vmin=0, vmax=100)
    
    if include_category_map:
        
        # collecting the set of DE values and normalizing for plotting
        change_vals = np.transpose(np.expand_dims(np.concatenate(ordered_values), axis=0))
        minval=np.min(change_vals)
        maxval=np.max(change_vals)
        
        divnorm=colors.TwoSlopeNorm(vmin=-1, vcenter=0., vmax=1)
        hotcoldmap = plt.cm.get_cmap('RdBu').reversed()
        occ = ax[0].imshow(change_vals, cmap=hotcoldmap, norm=divnorm, aspect='auto')
        
        cbar = ax[0].figure.colorbar(occ, ax=ax[0], extend='max', location='left')
        
        ax[1].imshow(np.vstack(ordered_heatmap), cmap='gray_r', aspect='auto', vmin=0, vmax=100)    
        ax[1].set(xticks=np.arange(0, nbins+1, 10),
        xticklabels=(np.arange(0, (windowSize*2)+1, (windowSize/2))-(windowSize*2)//2),
        xlabel='Distance from boundary (bp)')
        ax[1].set_title(title)
        
    else:
        
        ax.imshow(np.vstack(ordered_heatmap), cmap='gray_r', aspect='auto', vmin=0, vmax=100)
        ax.set(xticks=np.arange(0, nbins+1, 10),
        xticklabels=(np.arange(0, (windowSize*2)+1, (windowSize/2))-(windowSize*2)//2),
        xlabel='Distance from boundary (bp)')
        ax.set_title(title)
        
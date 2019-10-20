# l# last updated: July 19 2019 by Chris Mathy, chris.mathy@ucsf.edu
# code originally written by Colm Ryan
# Ryan, C. J., Roguev, A., Patrick, K., Xu, J., Jahari, H., Tong, Z., et al. (2012). Hierarchical modularity and the evolution of genetic interactomes across species. Molecular Cell, 46(5), 691–704. http://doi.org/10.1016/j.molcel.2012.05.028

import csv # 1.0
import numpy as np # 1.15.4
import pylab # 1.15.4
import pandas as pd # 0.24.1
import sys # 3.7.2
import os # 3.7.2
from scipy.interpolate import UnivariateSpline # 1.2.1
from scipy.stats import gaussian_kde
from statsmodels.graphics.gofplots import qqplot_2samples # 0.9.0
import matplotlib.pyplot as plt # 3.0.2

plt.rc('xtick',labelsize=6)
plt.rc('ytick',labelsize=6)
plt.rc('lines', markersize=2)


def get_ordered_tuple(item_a, item_b) :
    """
    Takes two items as input & returns an ordered tuple as the result. 
    Useful for interactions so that (A,B) & (B,A) are not treated seperately
    """
    if (item_a>item_b)-(item_a<item_b) <= 0:
        return (item_a, item_b)
    else :
        return (item_b, item_a)
        
def read_square_dataset_small(file_location, missing_string, separator,split = False, converter = lambda x : x, profiles = True,loop=0) :
    """
    Reads in a square data matrix from a symmetric emap file
    Returns a map of present interactions to scores, interaction profiles for each gene,
    and a list of names of genes present in the dataset.
    """
    interactions = {}
    gene_profiles = {}
    csv.register_dialect('custom', delimiter=separator)
    reader=csv.reader(open(file_location, 'r',10000000),'custom')
    gene_list = next(reader)
    gene_count = len(gene_list)
    #Removes the ' - DAMP' or ' DELETION' part
    for i in range(1,len(gene_list)) :
        if split :
            gene_list[i] = converter(gene_list[i].upper().split()[0].split('(')[0])
        else :
            gene_list[i] = converter(gene_list[i].upper())    
    for x in range(1,gene_count) :
        if gene_count % 100 == 0 :
            print(x, " of ", gene_count)
        linesx = next(reader)
        if split :      
            gene_x = converter(linesx[0].upper().split()[0].split('(')[0])
        else :
            gene_x = converter(linesx[0].upper())    
        for y in range(x+1+loop,gene_count) :
            gene_y = gene_list[y]
            if y < len(linesx) :
                interaction = linesx[y]
            else :
                interaction = missing_string
            if interaction != missing_string :
                score = float(interaction)
                xy = get_ordered_tuple(gene_x,gene_y)
                interactions[xy] = score
                if gene_x not in gene_profiles :
                    gene_profiles[gene_x] = {}
                if gene_y not in gene_profiles :
                    gene_profiles[gene_y] = {}
                if profiles :
                    gene_profiles[gene_y][gene_x] = score
                    gene_profiles[gene_x][gene_y] = score
    gene_list = gene_list[1:]
    return [interactions, gene_profiles, gene_list]

def output_delimited_text(filename,queries,arrays,correls,symm=False,converter = lambda x : x, sep="\t",missing = "",truncate = False):
    f = open(filename,'w')
    f.write("Gene")
    for i in arrays :
        f.write("%s%s" % (sep,converter(i)))
    f.write("\n")
    for i in queries :
        f.write("%s" % converter(i))
        for j in arrays :
            if symm :
                ij = get_ordered_tuple(i,j)
            else :
                ij = (i,j)
            if ij in correls :
                if truncate :
                    f.write("%s%4.4f" % (sep,correls[ij])) 
                else :
                    f.write("%s%s" % (sep,correls[ij]))
            else :
                f.write("%s%s" % (sep,missing))
        f.write("\n")
    f.close()
    return

def density_scatter_plot(intsA, intsB, outfile, xlabel, ylabel):
    """
    Plots a scatter plot showing shared interactions between two datasets
    
    input:
        intsA, intsB - dicts mapping pairs of genes to a score, format = ('str', 'str'): float
    
    returns:
        figure is printed out
        x - the list of shared scores from dataset A
        y - the list of shared scores from dataset B
        z - the point density
    """
    dfA = pd.DataFrame(data=list(intsA.items()), columns=['pair', 'scoreA'])
    dfB = pd.DataFrame(data=list(intsB.items()), columns=['pair', 'scoreB'])

    df = pd.merge(dfA, dfB, how='left', on='pair').dropna(how='any')

    x = list(df.scoreA)
    y = list(df.scoreB)

    # Calculate the point density
    xy = np.vstack([x,y])

    z = gaussian_kde(xy)(xy)

    fig, ax = plt.subplots(figsize=(2, 2))
    ax.scatter(x, y, c=z, s=2, edgecolor='')
    ax.set_xlabel(xlabel, fontname='Helvetica', fontsize=6)
    ax.set_ylabel(ylabel, fontname='Helvetica', fontsize=6)
    plt.savefig(outfile, format='png', transparent=True, bbox_inches='tight',dpi=300)

    
def main():

    outdir = '../../Supplemental_Figures/SGA_Scaling/scaler_output'

    # make output folder
    try:
        os.makedirs(outdir)
    except FileExistsError:
        pass

    # define datasets, datasetB is scaled to match datasetA
    datasetA = '../../Data/SGA_Scaling/cF3.txt'
    datasetB = '../../Data/SGA_Scaling/SGA_NxN_avg.txt'

    # read in the two datasets   
    ints, profs, genes = read_square_dataset_small(datasetA,
        "","\t",split=True,profiles = False)
    b_ints, b_profs, b_genes = read_square_dataset_small(datasetB,
        "","\t",split=True,profiles = False)

    datasetA = datasetA.split('/')[-1].split('.')[0]
    datasetB = datasetB.split('/')[-1].split('.')[0]

    avalues = []; bvalues=[]
    for i in ints : 
        if i in b_ints :
            avalues.append(ints[i])
            bvalues.append(b_ints[i])

    asorted = sorted(avalues)
    bsorted = sorted(bvalues)

    # shift datasetB so that it has the same number of negative values
    # as datasetA (makes it a little easier to scale)
    adjustment = -bsorted[len([x for x in asorted if x < 0])]
    bsorted = [x + adjustment for x in bsorted]

    # plot scatter plot showing shared interactions
    density_scatter_plot(ints,b_ints,outdir+'/unscaled_scatter.png',
                         xlabel='S-score', ylabel='SGA score')

    # record dataset information in log
    with open(outdir+'/scaler_log.txt', 'w') as f:
        f.write("cF3 EMAP has {} interactions\n".format(len(ints)))
        f.write("SGA_NxN has {} interactions\n".format(len(b_ints)))
        f.write("The sets have {} interactions in common\n".format(
                len(avalues)))
        f.write("Dataset correlation = {}\n".format(
                np.corrcoef(avalues,bvalues)[0][1]))
        f.write("Adjustment so that the SGA_NxN shared interaction "
                "set has the same number of negative values as the "
                "cF3 EMAP.\nadjustment={}\n".format(adjustment))

    ## Computing scaling values
    #essentially the data is partitioned into 100 overlapping bins
    #the mean value of bin[0] in datasetB is divided by the mean value of bin[0] from datasetA
    #this gives a scaling factor for values in the range (min(bin[0]), max(bin[0]))
    #values close to zero give unpredictable scaling factors, so they are ignored. 
    #Depending on the size of your overlap you may want to tweak the number of bins

    bins = 500
    binsize = len(avalues) / bins
    score = []; scale = []

    lower_threshold = 0.05
    upper_threshold = 0.99

    for i in np.arange(1,bins*lower_threshold) :
        start = int(i*binsize - binsize)
        end = int(i*binsize + binsize)
        score.append(np.mean(bsorted[start:end]))
        scale.append(np.mean(asorted[start:end])/np.mean(bsorted[start:end]))
    for i in np.arange(bins*upper_threshold,bins) :
        start = int(i*binsize - binsize)
        end = int(i*binsize + binsize)
        score.append(np.mean(bsorted[start:end]))
        scale.append(np.mean(asorted[start:end])/np.mean(bsorted[start:end]))

    # This function creates a curve which maps scores to scaling factors
    # the s=0.02 defines how close the curve fits your data points
    # large values give crap curves, small values may overfit your data
    # it's best to look at the resulting curve and tweak s= as appropriate
    svalue = 0.02
    s = UnivariateSpline(score, scale, s=svalue)

    #displays the scaling values(in red) and the fitted curve (in black)
    fig = plt.figure(figsize=(2, 2), dpi= 300, facecolor='w', edgecolor='k')
    plt.plot(np.arange(min(score),max(score),0.01), # changed from scatter
                [s(x) for x in np.arange(min(score),max(score),0.01)],
                color="red", linewidth=1)
    plt.scatter(score, scale, color="black", s=3)
    plt.xlim(1.1*min(score), 1.1*max(score))
    plt.ylim(0.9*min(scale), 1.1*max(scale))
    plt.ylabel('Scaling Factor', fontname='Helvetica', fontsize=6)
    plt.xlabel('SGA Score', fontname='Helvetica', fontsize=6)
    pylab.savefig(outdir+"/scaling_factor_curve.png", format='png',
                  transparent=True, bbox_inches='tight', dpi=300)

    # if the value to be scaled is larger than any value in our training set, we use
    # the scaling factor from the largest observed value
    def s_bounded(x) :
        if x<min(score) :
            x = min(score)
        elif x > max(score) :
            x=max(score)
        return s(x)

    #This function applies our scaling factor to a given value
    g= lambda x : (x + adjustment) * s_bounded(x + adjustment)

    for i in b_ints :
        b_ints[i] = float(g(b_ints[i]))

    scaled_dataset_file = "../../Data/SGA_Scaling/SGA_NxN_scaled_to_cF3.txt"
    output_delimited_text(scaled_dataset_file,b_genes,b_genes,b_ints,True)

    # save scaling info to log
    with open(outdir+'/scaler_log.txt', 'a') as f:
        f.write("Number of bins used: {}\n".format(bins))
        f.write("Lower threshold for bins: {}\n".format(
                lower_threshold))
        f.write("Upper threshold for bins: {}\n".format(
                upper_threshold))
        f.write("S value for fitting spline: {}\n".format(
                svalue))
        f.write("max_score={}\n".format(max(score)))
        f.write("min_score={}\n".format(min(score)))

    # save spline for scaling full SGA in R
    scores = np.arange(min(score),max(score),0.01)
    scales = [float(s(x)) for x in np.arange(min(score),max(score),0.01)]
    spline = pd.DataFrame(data=np.stack((scores,scales)).T,
                          columns=['score', 'scale'])
    spline.to_csv(outdir+'/spline.txt', sep='\t', index=False)

    # Plot scatter plot of shared interactions after scaling
    density_scatter_plot(ints, b_ints, outdir+'/scaled_scatter.png',
                         xlabel='S-score', ylabel='Scaled SGA score')
    
    # Make QQ Plots using the interactions before and after scaling
    avalues_after_scaling = []; bvalues_after_scaling = []
    for i in ints : 
        if i in b_ints :
            avalues_after_scaling.append(ints[i])
            bvalues_after_scaling.append(b_ints[i])

    # qqplot_2samples puts the "2nd Sample" on the x-axis
    # see documentation for statsmodels.graphics.gofplots
    fig, ax = plt.subplots(figsize=(2, 2))
    fig = qqplot_2samples(np.array(bvalues_after_scaling),
                                   np.array(avalues_after_scaling),
                                   line='r',
                                   ax=ax)
    ax.set_xlabel('S-score Quantiles',
                  fontname='Helvetica', fontsize=6)
    ax.set_ylabel('Scaled SGA score Quantiles',
                  fontname='Helvetica', fontsize=6)
    pylab.savefig(outdir+"/qq_scaled.png",
                  transparent=True, bbox_inches='tight', dpi=300)

    fig, ax = plt.subplots(figsize=(2, 2))
    fig = qqplot_2samples(np.array(bvalues),
                          np.array(avalues),
                          line='r',
                          ax=ax)
    ax.set_xlabel('S-score Quantiles', fontname='Helvetica', fontsize=6)
    ax.set_ylabel('SGA score Quantiles', fontname='Helvetica', fontsize=6)
    pylab.savefig(outdir+"/qq_unscaled.png",
                  transparent=True, bbox_inches='tight', dpi=300)

main()
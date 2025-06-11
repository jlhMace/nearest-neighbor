##########
# Code for graphing computational cost comparisons
##########


import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy




def plot_graph():
    pass

def main():
    ## CSV to DataFrame
    df_bins = pd.read_csv('data-graphing/data_spring25_bins.csv',header=0)
    df_nobins = pd.read_csv('data-graphing/data_spring25_nobins.csv',header=0)

    ## Polyfit
    coefficients = np.polyfit(df_nobins['atoms'], df_nobins['Total elapsed time'], 3)
    pl = np.poly1d(coefficients)

    ## Plot graphs
    #bins_total = sns.scatterplot(data=df_bins,x='atoms', y='Total elapsed time',label='Bin sort total time')
    #sns.regplot(data=df_bins, x='atoms', y='Total elapsed time',line_kws={"color": "purple"})
    #slope, intercept, r, p, sterr = scipy.stats.linregress(x=bins_total.get_lines()[0].get_xdata(),y=bins_total.get_lines()[0].get_ydata())
    #plt.text(12000,100, 'y = ' + str(round(intercept,3)) + ' + ' + str(round(slope,3)) + 'x')

    #less_total = sns.scatterplot(data=df_nobins,x='atoms', y='Total elapsed time',label='Without bin sort total time')
    #sns.regplot(data=df_nobins, x='atoms', y='Total elapsed time',line_kws={"color": "blue"},order=2)
    #plt.plot(df_nobins['atoms'], pl(df_nobins['atoms']), label='Quadratic Fit', color='blue')

    ## Plot bin sort times
    bins_sort = sns.scatterplot(data=df_bins, x='atoms', y='Bin sort', color='orange')
    bins_sort_reg = sns.regplot(data=df_bins, x='atoms', y='Bin sort',line_kws={"color": "purple"})
    slope, intercept, r, p, sterr = scipy.stats.linregress(x=bins_sort_reg.get_lines()[0].get_xdata(),y=bins_sort_reg.get_lines()[0].get_ydata())
    plt.text(12000,0.05, 'y = ' + str(round(intercept,3)) + ' + ' + str(round(slope,3)) + 'x')

    plt.ylabel('Total time (s)')
    plt.xlabel('Number of atoms')
    plt.title('Computational Cost of Bin Sort')
    #plt.legend()
    plt.show()
        

if __name__=='__main__':
    main()
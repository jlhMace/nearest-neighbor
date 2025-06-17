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
    df = pd.read_csv('data-graphing/orthorhombic_thin.csv',header=0)
    df_nb = pd.read_csv('data-graphing/orthorhombic_thin_nobin.csv',header=0)

    ## Plot graphs for bins
    df_one = df[df["width"]==1]
    df_two = df[df["width"]==2]
    bins_one = sns.regplot(data=df_one,x='atoms', y='Total elapsed time',label='10A thick, bins',color="blue")
    bins_two = sns.regplot(data=df_two,x='atoms', y='Total elapsed time',label='20A thick, bins',color="green")

    ## Plot graphs for no bins
    df_nb_one = df_nb[df["width"]==1]
    df_nb_two = df_nb[df["width"]==2]
    bins_nb_one = sns.regplot(data=df_nb_one,x='atoms', y='Total elapsed time',label='10A thick, no bins',color="red",order=2)
    bins_nb_two = sns.regplot(data=df_nb_two,x='atoms', y='Total elapsed time',label='20A thick, no bins',color='orange',order=2)

    ## For regression formula
    #slope, intercept, r, p, sterr = scipy.stats.linregress(x=df_one.get_lines()[0].get_xdata(),y=df_one.get_lines()[0].get_ydata())
    #plt.text(12000,100, 'y = ' + str(round(intercept,3)) + ' + ' + str(round(slope,3)) + 'x')

    plt.ylabel('Total time (s)')
    plt.xlabel('Number of atoms')
    plt.title('Computational Cost of Thin System')
    plt.legend()
    plt.savefig('data-graphing/orthorhombic_thin.png')
        

if __name__=='__main__':
    main()
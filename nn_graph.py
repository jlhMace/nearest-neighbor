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
    df = pd.read_csv('data-graphing/orthorhombic_thin_bin_2.csv',header=0)
    df_nb = pd.read_csv('data-graphing/orthorhombic_thin_nobin_2.csv',header=0)

    ## Plot graphs for bins
    df_one = df[df["width"]==1]
    df_two = df[df["width"]==2]
    bins_one = sns.regplot(data=df_one,x='atoms', y='Total elapsed time',label='10A thick, bins',color="blue")
    bins_two = sns.regplot(data=df_two,x='atoms', y='Total elapsed time',label='20A thick, bins',color="green")

    ## Plot graphs for no bins
    df_nb_one = df_nb[df_nb["width"]==1]
    df_nb_two = df_nb[df_nb["width"]==2]
    bins_nb_one = sns.regplot(data=df_nb_one,x='atoms', y='Total elapsed time',label='10A thick, no bins',color="red",order=2)
    bins_nb_two = sns.regplot(data=df_nb_two,x='atoms', y='Total elapsed time',label='20A thick, no bins',color='orange',order=2)

    ## For poly regression formula
    #x1 = df_nb[['atoms']].to_numpy().flatten()
    #y1 = df_nb[['Total elapsed time']].to_numpy().flatten()
    #eq1 = np.poly1d(np.polyfit(x1,y1, 2))
    #print(eq1)
    #polyspace = np.linspace(df[['atoms']].agg(['min']),df[['Total elapsed time']].agg(['max']))
    #plt.text(6000,6000, 'y = ' + str(round(eq1[2],3)) + ' + ' + str(round(eq1[1],3)) + 'x' + str(round(eq1[0],3)) + 'x^2')

    ## For linear regression formula
    #x2 = df[['atoms']].to_numpy().flatten()
    #y2 = df[['Total elapsed time']].to_numpy().flatten()
    #eq2 = np.polyfit(x2,y2,1)
    #plt.text(12000,50, 'y = ' + str(round(eq2[1],3)) + ' + ' + str(round(eq2[0],3)) + 'x')



    plt.ylabel('Total time (s)')
    plt.xlabel('Number of atoms')
    plt.title('Computational Cost of Thin System')
    plt.legend()
    plt.savefig('data-graphing/orthorhombic_thin_2.png')
        

if __name__=='__main__':
    main()
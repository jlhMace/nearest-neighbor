import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy


def plot_graph():
    pass

def main():
    ## Data to CSV
    df_bins = pd.DataFrame([[1350,27,0.01826039599836804,54.412598416995024,54.43223269599548],
                        [3200,64,0.04716818800079636,123.822724379992,123.8713450840005],
                        [6250,125,0.07624013899476267,238.63210353199975,238.7100271179952],
                        [10800,216,0.13078330700227525,455.56458601700433,455.6973391920037],
                        [17150,343,0.20886128500569612,707.5550146750029,707.766414194004]],
                        columns=['atoms','bins','Bin sort','Nearest neighbor','Total elapsed time'])
    csv_data = df_bins.to_csv('data.csv', index = False) 

    df_less = pd.DataFrame([[1350,49.903866401000414],
                            [3200,272.63898685600725],
                            [6250,970.8816703689954],
                            [10800,3178.949861598012],
                            [17150,7069.98548579801]],
                            columns=['atoms','Total elapsed time'])

    ## Polyfit
    coefficients = np.polyfit(df_less['atoms'], df_less['Total elapsed time'], 3)
    pl = np.poly1d(coefficients)

    ## Plot graphs
    #bins_total = sns.scatterplot(data=df_bins,x='atoms', y='Total elapsed time',label='Bin sort total time')
    #sns.regplot(data=df_bins, x='atoms', y='Total elapsed time',line_kws={"color": "purple"})
    #slope, intercept, r, p, sterr = scipy.stats.linregress(x=bins_total.get_lines()[0].get_xdata(),y=bins_total.get_lines()[0].get_ydata())
    #plt.text(12000,100, 'y = ' + str(round(intercept,3)) + ' + ' + str(round(slope,3)) + 'x')

    #less_total = sns.scatterplot(data=df_less,x='atoms', y='Total elapsed time',label='Without bin sort total time')
    #sns.regplot(data=df_less, x='atoms', y='Total elapsed time',line_kws={"color": "blue"},order=2)
    #plt.plot(df_less['atoms'], pl(df_less['atoms']), label='Quadratic Fit', color='blue')

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
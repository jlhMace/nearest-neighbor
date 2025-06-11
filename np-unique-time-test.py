import numpy as np
import random
from timer import Timer
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt



def add_data(list_length,range_length):
    timer = Timer("total")
    temp_list = []
    for i in range(list_length):
        temp_list.append(random.randint(0,range_length))
    timer.start()
    np.unique(temp_list)
    temp_time = timer.stop()
    return temp_time


def main():
    data = []
    list_length = 100
    range_length = 1000
    for i in range(list_length):
        for j in range(range_length):
            data.append([i, j, add_data(list_length, range_length)])

    df = pd.DataFrame(data,columns=['list','range','time'])
    sns.scatterplot(data=df, x='list', y='range', hue='time')
    plt.show()

if __name__=='__main__':
    main()
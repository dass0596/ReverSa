'''
Created on 31 de jul. de 2016

@author: Daniela Sanchez
'''
import matplotlib.pyplot as plt
import numpy as np

def graph(xList,yList , label, trueEvents):
    
    x = np.array(xList)
    y = np.array(yList)
    limit = float(trueEvents)
    y = y / limit 
   
    
    plt.plot(x,y, marker='o', linestyle='--', color='r')
    plt.ylim(0, 1)
    plt.xticks(x)
    plt.xlabel(label, fontsize=16)
    plt.ylabel('recovery rate', fontsize=16)
    plt.savefig('testReversa.pdf')
    #plt.show()

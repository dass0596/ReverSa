'''
Created on 31 de jul. de 2016

@author: daniela
'''

import matplotlib.pyplot as plt
import numpy as np

def graph(allResults, label, trueEvents):
    
    x = np.array(allResults.keys())
    y = np.array(allResults.values())
    limit = float(trueEvents)
    y = y / limit 
    #print y
    
    plt.figure(figsize=(8, 6), dpi=80)
    plt.plot(x,y, 'blue')
    plt.ylim(0, 1)
    plt.xticks(x)
    plt.xlabel(label, fontsize=16)
    plt.ylabel('recovery rate', fontsize=16)
    plt.savefig('testReversa.pdf')
    #plt.show()
    

'''
Created on 3 de jun. de 2016

@author: Daniela Sanchez
'''
import os
import shutil

def cleaning(queryFiles):
    currentDir = os.getcwd()
    os.remove('needle.txt')
    os.remove('needle1.txt')
    os.remove('water.txt')
    os.remove('distance.csv')
    for i in queryFiles:
        os.remove(i)
    dirPath = "%s/testing" %currentDir
    shutil.rmtree(dirPath) 
   
    pplacerDir = "%s/pplacer" %currentDir
    shutil.rmtree(pplacerDir)    
   

'''
Created on 3 de jun. de 2016

@author: Daniela Sanchez
'''
import os
import shutil

def cleaning(queryFiles):
    currentDir = os.getcwd()
    reFiles = ['needle.txt', 'needle1.txt', 'water.txt']
    
    for a in reFiles:
        if os.path.isfile(a):
            os.remove(a)

    for i in queryFiles:
        os.remove(i)
    dirPath = "%s/testing" %currentDir
    #shutil.rmtree(dirPath) 
   
    pplacerDir = "%s/pplacer" %currentDir
    #shutil.rmtree(pplacerDir)    
   

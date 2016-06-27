'''
Created on 3 de jun. de 2016

@author: Daniela Sanchez
'''
import os

def cleaning(queryFiles):
    currentDir = os.getcwd()
    os.remove('needle.txt')
    os.remove('needle1.txt')
    os.remove('water.txt')
    os.remove('distance.csv')
    for i in queryFiles:
        os.remove(i)
    dirPath = "%s/testing" %currentDir
    fileList = os.listdir(dirPath)
    for fileName in fileList:
        os.remove(dirPath+"/"+fileName)
        
    pplacerDir = "%s/pplacer" %currentDir
    files = os.listdir(pplacerDir)
    for f in files:
        os.remove(pplacerDir+"/"+f)

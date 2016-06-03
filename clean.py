'''
Created on 3 de jun. de 2016

@author: Daniela Sanchez
'''
import os

def cleaning():
    currentDir = os.getcwd()
    print currentDir
    dirPath = "%s/testing" %currentDir
    fileList = os.listdir(dirPath)
    for fileName in fileList:
        os.remove(dirPath+"/"+fileName)

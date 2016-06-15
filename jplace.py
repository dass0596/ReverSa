'''
Created on 9 de may. de 2016

@author: Daniela Sanchez
'''
import re

import json
import pandas as pd
import numpy as np
from ete3 import Tree
import statsmodels.formula.api as sm


class ParseJplace():
    
    def __init__(self, jfileComp, jfileMinus):
        self.jfileComp = jfileComp
        self.jfileMinus = jfileMinus
        self.perfil_comp = None 
        self.perfil_minus = None
        self.dir = None
        self.df= None
    
    
    def parseTree(self,tree, minus_sequence, edge_num):
        dict1 = {}
        win_id = []
        for i in edge_num:
            #add a name to nodes
            patt = re.compile(r'\)\:\d\.[0-9e\-]+\{%s\}' %i)
            obj_search = patt.search(tree)
            if obj_search is not None:
                pat_sub = re.compile(r'(\:\d\.[0-9e\-]+\{%s\})' %i)
                tree = re.sub(pat_sub, r'%s\1' %i, tree)    
            #extract the id of each branch or nodes
            patt1 = re.compile(r'(\w+)\:\d\.[0-9e\-]+\{%s\}' %i)
            obj1 = patt1.search(tree)
            dict1[i] = str(obj1.group(1))
        #Remove "{number}" from the tree 
        pattern=re.compile(r'\{\d+\}')
        new_tree = str(pattern.sub('',tree))
        nw = Tree(new_tree, format=1)
        #extract branch names for each node
        for node in nw.traverse("postorder"):
            if not node.is_leaf() and node.up: # For all internal nodes
                if node.name != '':
                    dict1[int(node.name)] = node.get_leaf_names()
                
        for y, x in dict1.items():
            if y in edge_num:
                if type(x) is list:
                    for k in x:
                        win_id.append(k)
                else:
                    win_id.append(x)
                        
        dist_fin =[]
        self.df = pd.read_csv('distance.csv', index_col = 0)       
        for a in win_id:
            p_dist = self.df.loc[minus_sequence,a]
            dist_fin.append(p_dist)  
        return dist_fin
    
    def parseMinus(self):
        dict1 = {}
        pattern = re.compile(r'_')
        pattern1 = re.compile(r'\w+_\w+_minus(.*)\.jplace')
        keys = ['n_place',  'phy_dist', 'maxdist'] #, 'pendant_len']
        self.dir = 'pplacer/'
        for x in self.jfileMinus:
            f_exp = pattern1.search(x)
            minus_win = f_exp.group(1)
            jason = '%s%s' %(self.dir, x)
            with open(jason) as data_file:
                data = json.load(data_file)
            
            for i in data['placements']: 
                for a in i['nm']:
                    id_win =  a[0]
                    id_seq =  pattern.split(id_win)[0]
            
                    if id_seq == minus_win:
                        n_place = len(i['p'])
                        values = [n_place]
                        if n_place == 1:
                            values.append(0)
                            for y in i['p']:
                                #values.append(y[0])
                                #pendant
                                #values.append(y[5])
                                #Append max distance
                                edgeVal= [y[1]]
                                distance = self.parseTree(data['tree'], minus_win, edgeVal)
                                distMax = max(distance)
                                values.append(distMax)                 
                        else:
                            #APPEND DIST
                            edge_val = [y[1] for y in i['p']]
                            dist = self.parseTree(data['tree'], minus_win, edge_val)
                            meanDist = np.mean(dist)
                            values.append(meanDist)
                            #print edge_val
                            #d_len = np.mean([y[0] for y in i['p']])
                            #p_len = np.mean([y[5] for y in i['p']]) 
                            #values.append(d_len)
                            #values.append(p_len)
                            #add max distance
                            maxDist = max(dist)
                            values.append(maxDist)
                        #print values 
                        dict1[id_win]= dict(zip(keys, values))   
                
        self.perfil_minus = pd.DataFrame.from_dict(dict1, orient = 'index')
        self.perfil_minus.to_csv('%sperfil_minus.csv' %self.dir)

    
    def parseSimple(self):
        self.parseMinus()
        dict1 = {}
        pattern = re.compile(r'_')
        keys = ['n_place', 'phy_dist', 'maxdist'] #, 'pendant_len']
        for x in self.jfileComp:
            jason = '%s%s' %(self.dir, x)
            with open(jason) as data_file:
                data = json.load(data_file)
                
            for i in data['placements']:
                n_place = len(i['p'])
                for a in i['nm']:
                    values = [n_place]
                    id_win =  a[0]
                    id_tree =  pattern.split(id_win)[0]
                    #print id_tree
                    
                    if n_place > 1:
                        edge_val = [y[1] for y in i['p']]
                        dist = self.parseTree(data['tree'], id_tree, edge_val)
                        meanDist = np.mean(dist)
                        values.append(meanDist)
                        #print edge_val
                        #d_len = np.mean([y[0] for y in i['p']])
                        #p_len = np.mean([y[5] for y in i['p']]) 
                        #values.append(d_len)
                        #values.append(p_len) 
                          
                                            
                    else:
                        values.append(0)
                        #for y in i['p']:
                            #values.append(y[0])
                            #values.append(y[5])
                            ##Add max distance
                            ##values.append(0)
                            
                    #Add max distance
                    distance= list(self.df.loc[id_tree,])
                    maxDist = max(distance)
                    values.append(maxDist)
                            
                    #print values                    
                    dict1[id_win]= dict(zip(keys, values))
            
        self.perfil_comp = pd.DataFrame.from_dict(dict1, orient = 'index')
        self.perfil_comp.to_csv('%sperfil_completo.csv' %self.dir)

    
    def correlation(self):
        self.parseSimple()
        dict_cor = {}
        for i in self.perfil_comp.index.tolist():
            comp = list(self.perfil_comp.loc[i,])
            minus = list(self.perfil_minus.loc[i,])
            df = pd.DataFrame({"A":minus, "B":comp})
            result = sm.ols(formula="A ~ B", data=df).fit()
            inter = result.params.Intercept
            dict_cor[i] = inter
        corr_coef = pd.DataFrame.from_dict(dict_cor, orient = 'index')
        corr_coef.columns=['intercept']
        corr_coef.to_csv('%sintercepts.csv' %self.dir)  
        return corr_coef
        
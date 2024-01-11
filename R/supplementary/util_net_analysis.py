# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 15:35:43 2023

@author: arnau
"""

import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np



def bipartite_network_constr(dataframe):
    # Network construction 
    net = nx.from_pandas_edgelist(dataframe,source='geneName', target='teName', edge_attr='coef', create_using=None, edge_key=None)
    # Set the bipartition attributes.
    ZNFnodes = list(dataframe.loc[:, 'geneName'].unique())
    TEnodes = list(dataframe.loc[:, 'teName'].unique())
    ZNFdict = dict.fromkeys(ZNFnodes,1)
    TEdict = dict.fromkeys(TEnodes,0)
    bipartitedict = dict(ZNFdict, **TEdict)
    nx.set_node_attributes(net, bipartitedict, name = 'bipartite' )
    #check
    print('Is the network connected ?')
    print(nx.is_connected(net))
    print('Is the network directed ?')
    print(nx.is_directed(net))
    print('Is the network bipartite ?')
    print(nx.is_bipartite(net))
    
    return net, ZNFnodes, TEnodes

def nodes_metrics(net, ZNFnodes, TEnodes):
    top_nodes = ZNFnodes
    # Dataframe creation.
    degreeZNF, degreeTE = bipartite.degrees(net, top_nodes)
    strengthZNF, strengthTE = bipartite.degrees(net, top_nodes, weight = 'coef')
    degreeZNF = pd.DataFrame(degreeZNF, columns = ['name', 'degree'])
    degreeZNF.set_index('name', inplace = True)
    degreeTE = pd.DataFrame(degreeTE,  columns = ['name', 'degree'])
    degreeTE.set_index('name', inplace = True)
    strengthZNF = pd.DataFrame(strengthZNF,  columns = ['name', 'strength'])
    strengthZNF.set_index('name', inplace = True)
    strengthTE = pd.DataFrame(strengthTE ,columns = ['name', 'strength'])
    strengthTE.set_index('name', inplace = True)
    ZNFdf = degreeZNF.merge(strengthZNF, left_index =True, right_index=True)
    ZNFdf.loc[:, "deg_centr"] = ZNFdf.loc[:, "degree"]/len(TEnodes)
    ZNFdf.loc[:, "strength_centr"] = ZNFdf.loc[:, "strength"]/len(TEnodes)
    ZNFdf.loc[:, "type"] = 'ZNF'
    TEdf = degreeTE.merge(strengthTE, left_index =True, right_index=True)
    TEdf.loc[:, "deg_centr"] = TEdf.loc[:, "degree"]/len(ZNFnodes)
    TEdf.loc[:, "strength_centr"] = TEdf.loc[:, "strength"]/len(ZNFnodes)
    TEdf.loc[:, "type"] = 'TE'
    node_metrics = pd.concat([ZNFdf,TEdf])
    node_metrics = node_metrics.merge(pd.Series(bipartite.betweenness_centrality(net, top_nodes), name = 'betweeness'), left_index=True, right_index=True)
    node_metrics = node_metrics[['degree', 'strength', 'deg_centr', 'strength_centr', 'betweeness','type']]
    # Partition attribution checking.
    first_znf_name = node_metrics.loc[node_metrics['type']=='ZNF'].iloc[1].name
    if first_znf_name in TEnodes:
        node_metrics['type'] = np.where(node_metrics['type']=='ZNF', 'TE', 'ZNF')

    return node_metrics

def violin_met(node_met):
    node_met["all"] = ""
    my_pal = {"ZNF": "purple", "TE" : "green"}
    for i in node_met.columns.to_list()[:-2] : 
        sns.violinplot(data=node_met, x = 'all', y = i, hue ='type', split = True,  palette = my_pal)
        plt.show()

def hub_list(node_met, threshold):
    hub_lists = {}
    for i in node_met.columns.to_list()[:-2] : 
        asc = False
        if i in ['strength', 'strength_centr']:
            asc = True
        hub_df = node_met.sort_values(i, axis=0, ascending = asc)
        n_lin = round(len(hub_df)*threshold)
        hub_df = hub_df.iloc[:n_lin]
        hlist = list(hub_df.index)
        hub_lists[i] = hlist
    return hub_lists

def hub_dict(node_met, threshold, save = False, save_name = '' ):
    hub_listed = hub_list(node_met, threshold)
    hub_dicted = {k: node_met.loc[v] for k, v in hub_listed.items()}
    # Saving.
    if save==True:
        with pd.ExcelWriter('../results/'+ save_name +'.xlsx') as writer:
            for df_name, df in hub_dicted.items():
                df.to_excel(writer, sheet_name=df_name)
    return hub_dicted
  
    
def assortativity(net, dataframe):
    znf_age = dataframe.loc[:,['geneName', 'kznf_age']]
    znf_age.set_index('geneName', inplace = True)
    znf_age = znf_age.to_dict()
    result = {}
    for key,value in znf_age.items():
        if value not in result.values():
            result[key] = value
    znf_age = result['kznf_age']
    te_age = dataframe.loc[:,['teName', 'te_age']]
    te_age.set_index('teName', inplace = True)
    te_age = te_age.to_dict()
    result = {}
    for key,value in te_age.items():
        if value not in result.values():
            result[key] = value
    te_age = result['te_age']
    agedict = dict(znf_age, **te_age)
    nx.set_node_attributes(net, agedict, name = 'age' )
    assort = nx.attribute_assortativity_coefficient(net, 'age', nodes=None)
    return assort

# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 16:11:39 2023

@author: arnau
"""
import pandas as pd
import networkx as nx
import itertools
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import scipy


def node_overl(df1,df2,df1_name = '', df2_name = ''):
    znf_1 = list(df1.loc[:, 'geneName'].unique())
    te_1 = list(df1.loc[:,  'teName'].unique())
    znf_2 = list(df2.loc[:, 'geneName'].unique())
    te_2 = list(df2.loc[:,  'teName'].unique())
    # Node intersections.
    common_znf = list(set(znf_1).intersection(znf_2))
    common_te = list(set(te_1).intersection(te_2)) 
    # Node specifics to 1 network.
    znf_only1 = list(set(znf_1).difference(znf_2))
    znf_only2 = list(set(znf_2).difference(znf_1))
    te_only1 = list(set(te_1).difference(te_2))
    te_only2 = list(set(te_2).difference(te_1))
    total_znf = list(itertools.chain(znf_only1, znf_only2, common_znf))
    total_te = list(itertools.chain(te_only1 , te_only2 , common_te))
    # DIctionnary construction.
    node_dict = {'znf_' + df1_name : znf_1 ,
                 'te_'+ df1_name : te_1,
                 'znf_' + df2_name: znf_2,
                 'te_' + df2_name : te_2,
                 'common_znf' : common_znf,
                 'common_te' : common_te,
                 'znf_only' + df1_name : znf_only1,
                 'znf_only' + df2_name : znf_only2 ,
                 'te_only' + df1_name : te_only1,
                 'te_only' + df2_name : te_only2,
                 'total_znf' : total_znf,
                 'total_te' : total_te }
    # Venn diagrams.
    venn2(subsets = (len(node_dict['znf_only' + df1_name]),  len(node_dict['znf_only' + df2_name]), len(node_dict['common_znf'])), set_labels = (df1_name + ' ZNF', df2_name + 'ZNF'))
    plt.show()
    venn2(subsets = (len(node_dict['te_only' + df1_name]),  len(node_dict['te_only' + df2_name]), len(node_dict['common_te'])), set_labels = (df1_name + ' te', df2_name + ' te'))
    plt.show()
    # Prints.
    print('\n')
    print('The ' + df1_name + ' network size is :')
    print(len(node_dict['znf_' + df1_name]) + len(node_dict['te_' + df1_name]))
    print('ZNF : ' + str(len(node_dict['znf_' + df1_name])) + ' and TE : ' + str(len(node_dict['te_' + df1_name])))
    
    print('\n')
    print('The ' + df2_name + ' network size is :')
    print(len(node_dict['znf_' + df2_name]) + len(node_dict['te_' + df2_name]))
    print('ZNF : ' + str(len(node_dict['znf_' + df2_name])) + ' and TE : ' + str(len(node_dict['te_' + df2_name])))
    
    print('\n')
    total_overlap = (len(node_dict['common_znf']) + len(node_dict['common_te'])) / ((len(node_dict['total_znf']) + len(node_dict['total_te'])))
    print('The Total overlap is : ' + str(total_overlap))
    znf_overlap = len(node_dict['common_znf']) / len(node_dict['total_znf'])
    print('The ZNF overlap is : ' + str(znf_overlap))
    te_overlap = len(node_dict['common_te']) / len(node_dict['total_te'])
    print('The TE overlap is : ' + str(te_overlap))
    
    return node_dict, total_overlap, znf_overlap, te_overlap

def overlap_node_df(df_list, df_names):
    default_df = pd.DataFrame(index = df_names, columns = df_names)
    result_dict = {'Total overlap' : default_df ,
                   'ZNF overlap' : default_df.copy(deep= True) ,
                    'TE overlap' : default_df.copy(deep= True) 
        }
    
    n = 1
    for a in result_dict.keys(): 
        for i in range(len(df_names)):
            for j in range(len(df_names)):
                result_dict[a].loc[df_names[i], df_names[j]] = node_overl(df_list[i], df_list[j], '', '' )[n]
        n += 1
    return result_dict

def dendo_funct(simil_datframe):
    dist_mat = 1 - simil_datframe
    dist_mat = dist_mat.to_numpy()
    dist_mat = scipy.spatial.distance.squareform(dist_mat)
    link = scipy.cluster.hierarchy.linkage(dist_mat, "single")
    scipy.cluster.hierarchy.dendrogram(link, labels = simil_datframe.index)
    plt.show()
    return None

def edge_overl(net1,net2):
    inters_net = nx.intersection(net1, net2)
    inters_size = inters_net.number_of_edges()
    compos_net = nx.compose(net1, net2)
    compos_size =  compos_net.number_of_edges()
    overlap = inters_size/compos_size
    return overlap


def overlap_edge_df(net_list, net_names):
    results_df = pd.DataFrame(index = net_names, columns = net_names)
    for i in range(len(net_list)):
        for j in range(len(net_list)):
            results_df.loc[net_names[i], net_names[j]] = edge_overl(net_list[i], net_list[j])
    return results_df

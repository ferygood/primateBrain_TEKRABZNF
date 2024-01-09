# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 15:28:30 2023

@author: Arnaud Maupas
"""
import pandas as pd
import matplotlib.pyplot as plt
import util_net_analysis as util
import util_comp
import scipy

# Loading.
human_df = pd.read_csv('../datasets/human_c1.csv')
bonobo_df = pd.read_csv('../datasets/bonobo_c1.csv') 
chimp_df = pd.read_csv('../datasets/chimpanzee_c1.csv') 
macaque_df = pd.read_csv('../datasets/macaque_c1.csv') 

# Bipartite network construction.
hum_net = util.bipartite_network_constr(human_df)[0]
bonobo_net = util.bipartite_network_constr(bonobo_df)[0]
chimp_net = util.bipartite_network_constr(chimp_df)[0]
macaque_net = util.bipartite_network_constr(macaque_df)[0]


# Overlaping nodes between 2 dataframes.
util_comp.node_overl(bonobo_df, chimp_df, 'macaque', 'bonobo' )

# Overlaping nodes for the whole set.
df_dict = util_comp.overlap_node_df([human_df, bonobo_df, chimp_df, macaque_df],['human', 'bonobo', 'chimp', 'macaque'] )
with pd.ExcelWriter('../results/node_overl.xlsx') as writer:
    for df_name, df in df_dict.items():
        df.to_excel(writer, sheet_name=df_name)

# Dendograms representations(see also util_comp.dendo_funct, sometimes an unidentified memory problem make it crash
#If this is the case, run one line after the other and print the dist _mat in between).
for a in df_dict.keys():
    dist_mat = 1 - df_dict[a]
    dist_mat = dist_mat.to_numpy()
    dist_mat = scipy.spatial.distance.squareform(dist_mat)
    link = scipy.cluster.hierarchy.linkage(dist_mat, "single")
    scipy.cluster.hierarchy.dendrogram(link, labels = df_dict[a].index)
    plt.title(a)
    plt.show()

#util_comp.dendo_funct(df_dict[a])

# Overlaping edges between two networks.
util_comp.edge_overl(hum_net,bonobo_net)
# Overlaping edges for the whole set.
edge_over_df = util_comp.overlap_edge_df([hum_net, bonobo_net, chimp_net, macaque_net],['human', 'bonobo', 'chimp', 'macaque'] )
edge_over_df.to_excel('../results/edge_over.xlsx')

# Dendogram ploting (same potential problem as previously mentioned).
dist_mat = 1 - edge_over_df
dist_mat = dist_mat.to_numpy()
dist_mat = scipy.spatial.distance.squareform(dist_mat)
link = scipy.cluster.hierarchy.linkage(dist_mat, "single")
scipy.cluster.hierarchy.dendrogram(link, labels = ['human', 'bonobo', 'chimp', 'macaque'])
plt.title(a)
plt.show()

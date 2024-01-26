from soma import aims 
import numpy as np 
import matplotlib.pyplot as plt 
from skimage.morphology import skeletonize_3d
import networkx as nx
import pandas as pd 
from scipy.spatial.distance import euclidean

scheme = pd.read_csv('../../SourceData/atlas_labelling_scheme.csv')

subj_id = 'CC00063AN06'
session_id = 15102

#subj_id = 'CC00065XX08'
#session_id = 18600

#for hemi, h in zip(['left', 'right'], ['L', 'R']):
for hemi, h in zip(['left'], ['L']):

    iLabel = aims.read('../../DerivedData/subjects/sub-{}/ses-{}/sub-{}_ses-{}_{}_texture_cortex_voronoi_majority.gii'.format(subj_id, session_id,subj_id, session_id,   hemi))
    iMesh = aims.read('/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-{}/ses-{}/anat/sub-{}_ses-{}_T2w_{}white_bv_transformed.gii'.format(subj_id, session_id,subj_id, session_id,  h))
    iCing = aims.read('../../DerivedData/subjects/sub-{}/ses-{}/sub-{}_ses-{}_hemi-{}_pathKmeans_128.label.gii'.format(subj_id, session_id,subj_id, session_id,  hemi))

    segm = np.array(iLabel)[:]
    vert, poly = np.array(iMesh.vertex()) ,np.array(iMesh.polygon())
    cing = np.array(iCing)[:]
    
    ## this should be done after removing cingulate area
    cing_l = np.unique(cing)[-1]
    segm[cing == cing_l] = 0

    ## this should be done after removing cingulate area
    res = {}

    for label in np.unique(segm)[1:]: 
        res[label] = []
        arg = np.argwhere(segm[segm==label]).ravel()
        partial = poly[arg]
    
   
        edges = [(edge[0], edge[1]) for edge in partial]
        edges.extend([(edge[1], edge[2]) for edge in partial])
        edges.extend([(edge[2], edge[0]) for edge in partial])
    
        G = nx.Graph()
        G.add_nodes_from(arg)
        G.add_edges_from(edges)
        #print(label, nx.number_connected_components(G))
        for component in list(nx.connected_components(G)):
            res[label].append(len(component)*100./len(arg))
            #print(len(component)*100./len(arg))
    
    fig, ax = plt.subplots(figsize=(14,6))

    for i, label in enumerate(np.unique(segm)[1:]):
        x = np.empty_like(np.array(res[label]))
        x[:] = i
    
        ax.scatter(x, res[label], c='gray', alpha=0.5)
    ax.set_xticks(range(len(np.unique(segm)[1:])))

    labels = [item.get_text() for item in ax.get_xticklabels()]


    names = [scheme[scheme.Label == label]['Abbreviation\xc2\xa0'].values[0] for label in np.unique(segm)[1:]]
    names = [ word.encode('ascii', 'replace').split('?')[0] for word in names]

    labels[:] = names
    ax.set_xticklabels(labels, rotation=90)
    #ax.set_xticklabels(np.unique(segm)[1:])
    #plt.scatter([res[label] for label in np.unique(segm)[1:]], )
    ax.set_title(hemi)
    ax.set_ylabel('% of vertices in fragment', fontsize=16)
    plt.show()



### Testing reassign 
hemi, h = 'left', 'L'
iLabel = aims.read('../../DerivedData/subjects/sub-{}/ses-{}/sub-{}_ses-{}_{}_texture_cortex_voronoi_majority.gii'.format(subj_id, session_id,subj_id, session_id,   hemi))
iMesh = aims.read('/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-{}/ses-{}/anat/sub-{}_ses-{}_T2w_{}white_bv_transformed.gii'.format(subj_id, session_id,subj_id, session_id,  h))
iCing = aims.read('../../DerivedData/subjects/sub-{}/ses-{}/sub-{}_ses-{}_hemi-{}_pathKmeans_128.label.gii'.format(subj_id, session_id,subj_id, session_id,  hemi))

segm = np.array(iLabel)[0][:]
vert, poly = np.array(iMesh.vertex()) ,np.array(iMesh.polygon())
cing = np.array(iCing)[0][:]

## this should be done after removing cingulate area
cing_l = np.unique(cing)[-1]
segm[cing == cing_l] = 0

## first, for every polygon create edges: 
edges = [(edge[0], edge[1]) for edge in poly]
edges.extend([(edge[1], edge[2]) for edge in poly])
edges.extend([(edge[2], edge[0]) for edge in poly])

## for every edge compute weight based on euclidean distance 
wg_edges = [(edge[0], edge[1], np.round(euclidean(vert[edge[0]], vert[edge[1]]),3)) for edge in edges]

G = nx.Graph()
nodes = [(i, {'label' : segm[i]}) for i in range(len(vert))]
G.add_nodes_from(nodes)
G.add_weighted_edges_from(wg_edges)


## test - should not be there 
l = 81
to_relabel = np.argwhere(segm == l).ravel()

for label in to_relabel:
    
    length = nx.single_source_dijkstra_path_length(G, label)
    
    ar = [[key, dji] for key, dji in zip(length.keys(), length.values())]
    ar.sort(key = lambda x: int(x[1]))

    for point in ar:
        if G.node[point[0]]['label'] != int(l):
            
            inNum = G.node[point[0]]['label']
            segm[label] = inNum
            print('Found new label: ', inNum)

            break 
    print(point)



### requires bv bash environment for brainvisa tools

import os, shutil
import numpy as np 
from soma import aims
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
from sklearn import cluster, mixture, neighbors
from skimage.measure import mesh_surface_area
from scipy.spatial.distance import pdist, squareform, cdist
from brainvisa.cortical_surface.surface_tools import PDE_tools as pdeTls
from scipy.spatial.distance import pdist, squareform

from sklearn.cluster.k_means_ import  _labels_inertia
from sklearn.metrics.pairwise import pairwise_distances_argmin_min
from scipy.optimize import minimize
from scipy.spatial.distance import braycurtis, euclidean, mahalanobis
import networkx as nx 

## ========================= RANDOM Parcellation =================================== ##

iTplFolder='/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/SourceData/DHCP_surf_templates/' 

## surface templates 
L_wm_template = os.path.join( iTplFolder, 'week-40_hemi-left_space-dhcpSym_dens-32k_wm.surf.gii')
R_wm_template = os.path.join( iTplFolder, 'week-40_hemi-right_space-dhcpSym_dens-32k_wm.surf.gii')

## cingulate manual parcellation (to exclude from the mesh)
L_cing_template =  os.path.join( iTplFolder, 'left_Cing_manual_parcel.gii')
R_cing_template = os.path.join( iTplFolder, 'right_Cing_manual_parcel.gii')

## ========================== General Functions =================================== ##

def get_non_cingulate_idx(file_path):
    
    text_cing = aims.read(file_path)
    atext_cing = text_cing[0].arraydata()

    return np.argwhere(atext_cing != 100)



def write_texture(mesh, labels, clustering, hemi, n_clusters, tmp_dir):
    tex = aims.TimeTexture('FLOAT')
    for i in range(mesh.size()):
        vert = np.asarray(mesh.vertex(i))
        tex[i].assign(np.zeros((len(vert),), dtype=np.float32))
        t = np.asarray(tex[i])
        t[:] = labels[:]
        
    if not os.path.exists(tmp_dir): 
        os.mkdir(tmp_dir)
        
    aims.write(tex, os.path.join(tmp_dir, 
                                 '{}_template_random_parc_{}_{}_clusters.tex.gii'.\
                                 format(hemi, clustering, n_clusters)))
    print('{} texture written'.format(clustering))

## ========================= K-means, reassign based on path + weight by number ===== ##
def to_minimize(weights, initial_distance):
    
    post_distance = initial_distance * weights
    labels, mindist = _arg_min_reduce(post_distance)
    
    return _error_function(labels) 

def _arg_min_reduce(dist):
    indices = dist.argmin(axis=1)
    values = dist[np.arange(dist.shape[0]), indices]
    return indices, values


def _error_function(labels):
    
    error = np.full(len(np.unique(labels)), -1, np.int32)
    expected = np.ceil(len(labels)/len(np.unique(labels)))
    
    error = np.asarray([(len(labels[labels == label]) - expected)**2 for label in np.unique(labels)])
    ssq =  np.sqrt(np.sum(error))
    print(ssq)
    return ssq

## get best centroids using K-means
def get_best_vertex_to_centroids_idx(vertices, distance_metric= 'mahalanobis'):
    k_means = cluster.MiniBatchKMeans( 
            n_clusters=n_clusters, init='k-means++', max_iter=200, 
            batch_size=200, verbose=False, compute_labels=True, 
            max_no_improvement=None, n_init=10, reassignment_ratio=0.5, random_state = 42)
    k_means.fit(vertices)
    
    centroids = k_means.cluster_centers_ 
    dist = cdist(vertices, centroids, metric= distance_metric)
    new_centroid_idx = np.argmin(dist, axis=0)
    return new_centroid_idx

def calculate_dijkstra_path_matrix(vertices, faces, new_centroid_idx):
    G = nx.Graph()

    ## edges 
    edges1 = [ (faces[i,0], faces[i,1], euclidean(vertices[faces[i,0]], vertices[faces[i,1]]) ) for i in range(len(faces))]
    edges2 = [ (faces[i,1], faces[i,2], euclidean(vertices[faces[i,1]], vertices[faces[i,2]]) ) for i in range(len(faces))]
    edges3 = [ (faces[i,2], faces[i,0], euclidean(vertices[faces[i,2]], vertices[faces[i,0]]) ) for i in range(len(faces))]

    edges = edges1 + edges2 + edges3

    G.add_weighted_edges_from(edges)
    
    # get paths lengths
    djikstra = np.empty((len(vertices), len(new_centroid_idx)))
    djikstra[:] = np.nan
    
    for idx in range(len(new_centroid_idx)):
        ctr = new_centroid_idx[idx] 
        #print(idx+1)
        length = nx.single_source_dijkstra_path_length(G, ctr, weight='weight')
    
        for vertex in range(len(vertices)):
            djikstra[vertex, idx] = length[vertex]
        
    return djikstra


## ======================== Run & save ========================================== ##

out_dir='/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/random_parcellation/'
 
#clusters = [64,128,256,512]
### this is for the replication in case reviewers ask for the same as neonate atlas number + for networks paper
clusters = [26,36]

for n_clusters in clusters:
    print('***{} clusters***'.format(n_clusters))  
    
    for hemi in ['left', 'right']:
        print('***{}***'.format(hemi))
    
        if hemi == 'left':
            mesh = aims.read(L_wm_template)
            bcg_idx = get_non_cingulate_idx(L_cing_template)   
            
        else:
            mesh = aims.read(R_wm_template)
            bcg_idx = get_non_cingulate_idx(L_cing_template) 
    
        vertices, faces = np.asarray(mesh.vertex()), np.asarray(mesh.polygon())
    
        ###
        samples = vertices[bcg_idx][:,0]
    
        print('getting centroids')
        new_centroid_idx = get_best_vertex_to_centroids_idx(vertices=samples, 
                                                        distance_metric= 'mahalanobis')
        #new_centroid_idx, init_weights = get_initial_weightset_best_vertex_to_centroids_idx(vertices=samples, 
        #                                                    n_clusters=n_clusters, 
        #                distance_metric= 'mahalanobis', init_weights=True)
        
        previous_centroid_idx = bcg_idx[new_centroid_idx].ravel()
    
        print('getting paths')
        distance_matrix = calculate_dijkstra_path_matrix(vertices=vertices, 
                                                         faces=faces, 
                                                     new_centroid_idx=previous_centroid_idx)
    
        
    
        new_dist = distance_matrix[bcg_idx][:,0]
        
    
        ### optimizing
        init_weights= np.ones(n_clusters, dtype=np.int16)
        
        fun = lambda x: to_minimize(x, new_dist)
        res = minimize(fun,  init_weights,method='Nelder-Mead', bounds=(0,2))
        final_dist = new_dist  * res.x
        l, mins = _arg_min_reduce(final_dist)
               
        out_labels = np.full(len(vertices),n_clusters+1)
        np.put(out_labels,bcg_idx, l)
        
        
        write_texture(mesh=mesh, 
                      labels=out_labels, 
                      clustering='pathKmeans', 
                      hemi=hemi, 
                      n_clusters=n_clusters, 
                      tmp_dir = out_dir)



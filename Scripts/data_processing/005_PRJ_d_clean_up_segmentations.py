# requires bv 
from soma import aims 
import numpy as np 
#import matplotlib.pyplot as plt 
#from skimage.morphology import skeletonize_3d
import networkx as nx
import pandas as pd 
from scipy.spatial.distance import euclidean
import datetime, os, sys
from networkx.algorithms.distance_measures import center

### segmentation
R_labels = [67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99,101,103,105,107,109,111,113,117,121]
L_labels = [68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,118,122]


def get_main_graph(subj_id, session_id, hemi):
    
    if hemi == 'left':
        h = 'L'
    else:
        h = 'R'
        
    # inDATA 
    iLabel = aims.read('../../DerivedData/subjects/sub-{}/ses-{}/sub-{}_ses-{}_{}_texture_cortex_voronoi_majority.gii'.format(subj_id, session_id,subj_id, session_id,   hemi))
    iMesh = aims.read('/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-{}/ses-{}/anat/sub-{}_ses-{}_T2w_{}white_bv_transformed.gii'.format(subj_id, session_id,subj_id, session_id,  h))
    iCing = aims.read('../../DerivedData/subjects/sub-{}/ses-{}/sub-{}_ses-{}_hemi-{}_pathKmeans_128.label.gii'.format(subj_id, session_id,subj_id, session_id,  hemi))

    
    segm = np.array(iLabel)[0][:].ravel()
    vert, poly = np.array(iMesh.vertex()) ,np.array(iMesh.polygon())
    cing = np.array(iCing)[0][:].ravel()

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
    
    return G, segm


def correct_hemi(G, segm, to_remove):
    
    '''params:
    G: wm graph 
    segm: segmentation array
    to_remove: labels to look for and relabel
    
    To speed up, this does work per 'blob' rather than vertex-by-vertex which might lead to some small errors based on 
    distance but significantly speeds up the runtime
    
    '''
    
    ## find blobs
    for label in to_remove:
        
        # to create subgraph with only those labels
        to_keep = [node for node in G.nodes() if G.node[node]['label'] == label]
        H = G.subgraph(to_keep)
        
        # get connected elements and list of blobs to relabel
        connected = list(nx.connected_components(H))
	#S = [center(H.subgraph(c).copy())[0] for c in nx.connected_components(H)]

        #for rpr, blob in zip(S,connected):
	for blob in connected:
            # get random vertex from the blob to find the closes new label based on path 
            rpr = np.random.choice(list(blob))
            print(rpr)

            length = nx.single_source_dijkstra_path_length(G, rpr)
            ar = [[key, dji] for key, dji in zip(length.keys(), length.values())]
            ar.sort(key = lambda x: int(x[1]))
            
            for point in ar:
        
                if G.node[point[0]]['label'] not in to_remove :
                    inNum = G.node[point[0]]['label']
                    #print('Blob new label:', inNum)
                    break
            
            ## relabel whole blob 
            print('Size of blob: ', len(list(blob)), 'previous label:', label, 'new label:', inNum)
            segm[list(blob)[:]] = inNum
            
    ## return new segmentation and re-labelled graph 
    ### need to relabel graph with new segmentation!!! 
    new = dict(zip(range(len(segm)), segm))
    nx.set_node_attributes(G, 'label', new)
    
    print(np.unique(segm))
        
    return G, segm 



def relabel_blobs(G, segm, ths=20):
    
    labels = np.unique(segm)[1:]
    for label in labels:
        # to create subgraph
        to_keep = [node for node in G.nodes() if G.node[node]['label'] == label]
        H = G.subgraph(to_keep)
        
        # get connected elements and list of blobs to relabel
        connected = list(nx.connected_components(H))
        sizes = [(len(sub) * 100.)/len(H.nodes()) for sub in connected]
	#S = [center(H.subgraph(c).copy())[0] for c in nx.connected_components(H)]

        remove_blobs = []
	#centre_nodes = []

        for i, size in enumerate(sizes):
            if size < ths:
                remove_blobs.append(connected[i])
		#centre_nodes.append(S[i])
                
        ## get new label per blob to speed the stuff up
        #for rpr, blob in zip(centre_nodes,remove_blobs):
	for blob in remove_blobs:
            rpr = np.random.choice(list(blob))
    	    print(rpr)
            length = nx.single_source_dijkstra_path_length(G, rpr)
            ar = [[key, dji] for key, dji in zip(length.keys(), length.values())]
            ar.sort(key = lambda x: int(x[1]))
    
            for point in ar:
        
                if G.node[point[0]]['label'] != G.node[rpr]['label'] :
                    inNum = G.node[point[0]]['label']
                    #print('Blob new label:', inNum)
                    break
            
            ## relabel whole blob 
            print('Size of blob: ', len(list(blob)), 'previous label:', label, 'new label:', inNum)
            segm[list(blob)[:]] = inNum
    
    ### need to relabel graph with new segmentation!!! 
    new = dict(zip(range(len(segm)), segm))
    nx.set_node_attributes(G, 'label', new)
    
    return G, segm 


def get_error(d):
    error = 0 
    for key in d.keys():
        for el in d[key]:
            error = error + (100-el)
    return error


def run_correction(subj_id, session_id, hemi):

	G, in_segm = get_main_graph(subj_id=subj_id, session_id=session_id, hemi=hemi)

	if hemi == 'left':
		newG, new_segm = correct_hemi(G=G.copy(), segm=in_segm.copy(), to_remove=R_labels)
	else:
		newG, new_segm = correct_hemi(G=G.copy(), segm=in_segm.copy(), to_remove=L_labels)

	## iterate
	outG, out_segm = relabel_blobs(G=newG.copy(), segm=new_segm.copy(), ths=15)

	for i in range(5):
    		outG, out_segm = relabel_blobs(G=outG.copy(), segm=out_segm.copy(), ths=15)
	## to do the final cleaning of cingulate area, don't do this, might be dangerous
	#outG, out_segm = relabel_blobs(G=outG.copy(), segm=out_segm.copy(), ths=50)

	if hemi == 'left':
		h='L'
	else: 
		h='R'

	mesh = aims.read('/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-{}/ses-{}/anat/sub-{}_ses-{}_T2w_{}white_bv_transformed.gii'.format(subj_id, session_id,subj_id, session_id,  h))
	vert, poly = np.array(mesh.vertex()) ,np.array(mesh.polygon())
	tex = aims.TimeTexture('FLOAT')
	for i in range(mesh.size()):
    		tex[i].assign(np.zeros((len(vert),), dtype=np.float32))
    		t = np.asarray(tex[i])
    		t[:] = out_segm[:]

	outf = '../../DerivedData/subjects/sub-{}/ses-{}/sub-{}_ses-{}_neoAtlas_corrected_{}.gii'.format(subj_id, session_id, subj_id, session_id, hemi)
	aims.write(tex, outf)
	if os.path.isfile(outf):
		print('CORRECTED {} hemi {}'.format(subj_id, hemi))
	else: 
		print('BAD {} hemi {}'.format(subj_id, hemi))

if __name__ == "__main__":
	print('START at: {}'.format(datetime.datetime.now()))

	scheme = pd.read_csv('../../SourceData/atlas_labelling_scheme.csv')

	### read in the file containing the list of subject ID and session IDs
	if len(sys.argv) < 2:
		print("You must provide subject file!")
		sys.exit()
	else:
		subject_file = sys.argv[1]

	subjects = pd.read_csv(subject_file, names=['subject_id', 'session_id', 'template'])
	# remove the bad metrics
	bad = pd.read_csv('../../DerivedData/failed_metric_QC.csv', names=['subject_id'])
	subjects = subjects[~ subjects.subject_id.isin(bad.subject_id.values)]
	print('number of subjects:', len(subjects))

	for i, row in subjects.iterrows():  

		subj_id = row.subject_id
		session_id = row.session_id

        	print('Running subject: {}'.format(subj_id))
		
		### left 
		print('Running LEFT:...')

		if os.path.isfile('../../DerivedData/subjects/sub-{}/ses-{}/sub-{}_ses-{}_left_texture_cortex_voronoi_majority.gii'.format(
								subj_id,session_id, subj_id, session_id	)):
	
			out_segm = run_correction(subj_id=subj_id, session_id=session_id, hemi='left')
		else: 
			with open('missing_cortex_segm_texture', 'a') as the_file:
    				the_file.write('{},{},left\n'.format(subj_id, session_id))
		
		### right 
		print('Running RIGHT:...')
		if os.path.isfile('../../DerivedData/subjects/sub-{}/ses-{}/sub-{}_ses-{}_right_texture_cortex_voronoi_majority.gii'.format(
								subj_id,session_id, subj_id, session_id	)):
	
			out_segm = run_correction(subj_id=subj_id, session_id=session_id, hemi='right')
		else: 
			with open('missing_cortex_segm_texture', 'a') as the_file:
    				the_file.write('{},{},right\n'.format(subj_id, session_id))
	print('END at: {}'.format(datetime.datetime.now()))

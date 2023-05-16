import os
from utils import *
from os.path import join
from pprint import pprint

import numpy as np
from scipy.stats import circmean, circstd

import requests
import xml.etree.ElementTree as ET

#### Paths/to/stuff ####
cwd = os.getcwd()
path_to_save = cwd



#### User-defined functions ####
def get_model(nmldb_id,save=False,path_to_save=None):
    '''
        NeuroML-DB API query for model
        TODO: Check if file is store already somewhere and load that instead of calling API.
        ----
        PARAMETER:
            - nmldb_id : (str) Contains NML ID
            - save : (bool)
            - path_to_save : (str)
        OUTPUT:
            - mldb_model_response : (str) Contains all model specifications
    '''


    nmldb_xml_url = 'https://neuroml-db.org/render_xml_file?modelID='

    model_xml_url = nmldb_xml_url + nmldb_id

    xml_name = nmldb_id+'.xml'

    filename = join(path_to_save,xml_name)
    
    model_xml_response = requests.get(model_xml_url)
    print(nmldb_id,'model requested and received...')

    # If file exists, load and return
    if os.path.isfile(filename):
           
        
        print('Returning model file as string')
        return model_xml_response.text
    
    else:
        
        # save file
        if save:

            if path_to_save:
                filename = join(path_to_save,xml_name)
            else:
                filename = xml_name

            with open(filename, 'wb') as file:
                file.write(model_xml_response.content)
                
        print('Returning model file as string')
        return model_xml_response.text

                


#### Define morphology domain of interest ####
spacing = 10

# NOTE: can set bounds to np.inf if desired
# NOTE: Works for cylindrical bounds only

# bound in Y
if spacing == 10:
    y_lb, y_ub = -155, 155
elif spacing == 20:
    y_lb, y_ub = -310,310
else:
    pass


# bound in XZ-plane
rad_ub = 50


#### Define models of interest ####
model_ids = [
    # L5_TTPC1
    'NMLCL000687',
    'NMLCL000688',
    'NMLCL000689',
    'NMLCL000690',
    'NMLCL000691',

    # L5_TTPC2
    'NMLCL000695',
    'NMLCL000692', 
    'NMLCL000693', 
    'NMLCL000694', 
    'NMLCL000696', 

    # L5_UTPC
    'NMLCL000698',
    'NMLCL000699', 
    'NMLCL000701',
    'NMLCL000700',
    'NMLCL000697',

    # L5_STPC
    'NMLCL000684',
    'NMLCL000685',
    'NMLCL000682',
    'NMLCL000686',
    'NMLCL000683', 

    # L6_TPC_L4
    'NMLCL000672',
    'NMLCL000676',
    'NMLCL000673',
    'NMLCL000674', 
    'NMLCL000675',

    # L6_TPC_L1
    'NMLCL000669',
    'NMLCL000670',
    'NMLCL000667',
    'NMLCL000671',
    'NMLCL000668', 

    # L6_IPC
    'NMLCL000662',
    'NMLCL000666',
    'NMLCL000663',
    'NMLCL000664',
    'NMLCL000665',

    # L6_UTPC
    'NMLCL000680',
    'NMLCL000677',
    'NMLCL000681',
    'NMLCL000678',
    'NMLCL000679',

    # L6_BPC
    'NMLCL000659',
    'NMLCL000660', 
    'NMLCL000657',
    'NMLCL000661',
    'NMLCL000658',

    # L4_SS
    'NMLCL000652',
    'NMLCL000653',
    'NMLCL000654',
    'NMLCL000655',
    'NMLCL000656',

    # L4_SP
    'NMLCL000647',
    'NMLCL000648',
    'NMLCL000649',
    'NMLCL000650',
    'NMLCL000651',

    # L23_PC
    'NMLCL000637',
    'NMLCL000638',
    'NMLCL000639',
    'NMLCL000640',
    'NMLCL000641',

    # L23_NBC
    'NMLCL000494',
    'NMLCL000492',
    'NMLCL000496',
    'NMLCL000495', 
    'NMLCL000493', 

    # L23_MC
    'NMLCL000487',
    'NMLCL000491', 
    'NMLCL000488', 
    'NMLCL000489', 
    'NMLCL000490', 

    # L23_BTC
    'NMLCL000469',
    'NMLCL000487', 
    'NMLCL000491',
    'NMLCL000488',
    'NMLCL000490',

    # L23_DBC
    'NMLCL000480',
    'NMLCL000479',
    'NMLCL000477',
    'NMLCL000481', 
    'NMLCL000478',

    # L23_BP
    'NMLCL000462',
    'NMLCL000466',
    'NMLCL000463', 
    'NMLCL000464', 
    'NMLCL000465', 

    # L23_LBC
    'NMLCL000484',
    'NMLCL000485', 
    'NMLCL000482',
    'NMLCL000486',
    'NMLCL000483',

    # L23_SBC
    'NMLCL000505',
    'NMLCL000502',
    'NMLCL000506',
    'NMLCL000503',
    'NMLCL000504',

    # L23_ChC
    'NMLCL000473',
    'NMLCL000474', 
    'NMLCL000476',
    'NMLCL000472',
    'NMLCL000475', 

    # L23_NGC
    'NMLCL000501',
    'NMLCL000499',
    'NMLCL000498',
    'NMLCL000497',
    'NMLCL000500'
]

cell_labels = ['L5_TTPC1_%s'%(i+1) for i in range(5)]+  \
              ['L5_TTPC2_%s'%(i+1) for i in range(5)]+   \
              ['L5_UTPC_%s'%(i+1) for i in range(5)]+  \
              ['L5_STPC_%s'%(i+1) for i in range(5)]+  \
              ['L6_TPC_L4_%s'%(i+1) for i in range(5)]+\
              ['L6_TPC_L1_%s'%(i+1) for i in range(5)]+\
              ['L6_IPC_%s'%(i+1) for i in range(5)]+ \
              ['L6_UTPC_%s'%(i+1) for i in range(5)]+ \
              ['L6_BPC_%s'%(i+1) for i in range(5)]+ \
              ['L4_SS_%s'%(i+1) for i in range(5)]+ \
              ['L4_SP_%s'%(i+1) for i in range(5)]+ \
              ['L23_PC_%s'%(i+1) for i in range(5)]+ \
              ['L23_NBC_%s'%(i+1) for i in range(5)]+ \
              ['L23_MC_%s'%(i+1) for i in range(5)]+ \
              ['L23_BTC_%s'%(i+1) for i in range(5)]+ \
              ['L23_DBC_%s'%(i+1) for i in range(5)]+ \
              ['L23_BP_%s'%(i+1) for i in range(5)]+ \
              ['L23_LBC_%s'%(i+1) for i in range(5)]+ \
              ['L23_SBC_%s'%(i+1) for i in range(5)]+ \
              ['L23_ChC_%s'%(i+1) for i in range(5)]+ \
              ['L23_NGC_%s'%(i+1) for i in range(5)]



for model_i, (NMLDB_ID,pop) in enumerate(zip(model_ids,cell_labels)):
    print('Computing morphometrics for %s'%NMLDB_ID)
    
    # get model file
    try:
        neuron_model = get_model(NMLDB_ID,path_to_save=path_to_save,save=True)
        
        # get directly from neuron_model string
        root = ET.fromstring(neuron_model)
          
    except requests.exceptions.SSLError:
    
        model_file = join(path_to_save,NMLDB_ID+'.xml')
        model_tree = ET.parse(model_file)
        root = model_tree.getroot()
    
    
    
    # cell model root and morphology elements
    cell = root[-1]
    morpho = cell[1]
    
    
    # high-level grouping of compartment domains
    segmentGroups = {}
    for i, group in enumerate(morpho.findall('{http://www.neuroml.org/schema/neuroml2}segmentGroup')):

        segmentGroups.update({i:group})
        
    # has compartment domain data
    comp_domains = {}
    for key, val in segmentGroups.items():

        if val.attrib['id'] in ['somatic','basal','axonal','apical']:

            comp_domains.update({val.attrib['id']:val})

    
    #### Find all somatic segments ####
    soma_group = comp_domains['somatic'][0].attrib['segmentGroup']
    soma_segments = {}
    for segment in morpho.findall('{http://www.neuroml.org/schema/neuroml2}segment'):
        if soma_group in segment.attrib['name']:
            soma_segments.update({segment.attrib['id']:segment.attrib['name']})
            
    ###################       
    #### Feature 1 ####
    ###################
    print('...soma stuff')
    these_segments = soma_segments
#     soma_length = 0
    soma_parents = {} # used elsewhere
    
    for segment in morpho.findall('{http://www.neuroml.org/schema/neuroml2}segment'):
        seg_id = segment.attrib['id']
        if seg_id in these_segments.keys():

            # this part only applies to soma compartments
            if 'Seg0' in segment.attrib['name']: 
                proximal = segment[0]
                distal = segment[1]

                # most medial soma segment
                soma_y_prox = np.float(proximal.attrib['y'])
                soma_y_dist = np.float(distal.attrib['y'])


                # get direction of prox->dist
                if soma_y_prox<soma_y_dist:
                    lb = soma_y_prox
                    ub = soma_y_dist # this needs to be updated
                    soma_dir = 1 # prox = more medial, dist = more lateral
                    soma_root = lb # seg set to origin in simulation
                else:
                    lb = soma_y_dist # this needs to be updated
                    ub = soma_y_prox 
                    soma_dir = -1 # prox = more lateral, dist = more medial
                    soma_root = ub # seg set to origin in simulation
            else:

                # find parent
                parent_id = segment[0].attrib['segment']

                proximal = soma_parents[parent_id]
                distal = segment[1]

                # update until find most lateral soma segment
                if soma_dir==1 and np.float(distal.attrib['y'])>ub:
                    ub = np.float(distal.attrib['y'])
                if soma_dir==-1 and np.float(distal.attrib['y'])<lb:
                    lb = np.float(distal.attrib['y'])

            soma_parents.update({seg_id:distal})

    #         soma_length += compute_segment_length(proximal,distal)

    # compute soma center shift
    soma_center = lb + (ub-lb)/2 # this is used to define upper/lower groups of neurites
    soma_length = ub-lb
    soma_lb = lb
    soma_ub = ub

#     for segment in morpho.findall('{http://www.neuroml.org/schema/neuroml2}segment'):
#         seg_id = segment.attrib['id']
#         if seg_id in these_segments.keys():

#             # this part only applies to soma compartments, other you have to find parent
#             if 'Seg0' in segment.attrib['name']:
#                 proximal = segment[0]
#                 distal = segment[1]

#                 # most medial soma segment
#                 soma_y_prox = proximal.attrib(['y'])
#                 soma_y_dist = distal.attrib(['y'])

#                 # get direction of prox->dist
#                 if soma_y_prox<soma_y_dist:
#                     lb = soma_y_prox
#                     ub = soma_y_dist # this needs to be updated
#                     soma_dir = 1 # prox = more medial, dist = more lateral
#                     soma_root = lb # seg set to origin in simulation
#                 else:
#                     lb = soma_y_dist # this needs to be updated
#                     ub = soma_y_prox 
#                     soma_dir = -1 # prox = more lateral, dist = more medial
#                     soma_root = ub # seg set to origin in simulation
#             else:

#                 # find parent
#                 parent_id = segment[0].attrib['segment']

#                 proximal = soma_parents[parent_id]
#                 distal = segment[1]

#                 # update until find most lateral soma segment
#                 if soma_dir==1 and distal.attrib(['y'])>ub:
#                     ub = distal.attrib(['y'])
#                 if soma_dir==-1 and disal.attrib(['y'])<lb:
#                     lb = distal.attrib(['y'])

#             soma_parents.update({seg_id:distal})

#             soma_length += compute_segment_length(proximal,distal)

#     # compute soma center shift
#     soma_center = lb + (ub-lb)/2 # this is used to define upper/lower groups of neurites

    #### Find all stems ####
    print('...stem stuff')
    stems = {}

    for segment in morpho.findall('{http://www.neuroml.org/schema/neuroml2}segment'):
        seg_id = segment.attrib['id']
        seg_name = segment.attrib['name']

        if not seg_name in soma_segments.values():

            parent_id = segment[0].attrib['segment']

            if parent_id in soma_parents.keys():
                stems.update({seg_id:seg_name})
                
    # compute stems cross-sectional areas
    stem_names = []
    stem_diams = []
    
    # group CSAs and initial angles
#     apic_csas = []
    upper_axon_csas = []
    lower_axon_csas = []
    basal_csas = []
    upper_csas = []
    lower_csas = []
    upper_basal_csas = []
    lower_basal_csas = []

    for segment in morpho.findall('{http://www.neuroml.org/schema/neuroml2}segment'):
        seg_name = segment.attrib['name']

        if seg_name in stems.values():

            proximal = segment[1] # always connected to soma
            distal = segment[2]
            diameter = float(proximal.attrib['diameter'])
            distal_y = float(distal.attrib['y'])
            proximal_y = float(proximal.attrib['y'])
            
            stem_names.append(seg_name)
            stem_diams.append(diameter)
            
            # separate upper and lower basal dendrites
            if proximal_y>=soma_center and not 'axon' in seg_name: # soma_center used to be 0 and compared to distal_y (deprecated)
                
                # only basal dendrites
                if 'dend' in seg_name:
                    upper_basal_csas.append(compute_stem_csa(diameter))
                
                # used to compute the remainder for non-basal segments
                upper_csas.append(compute_stem_csa(diameter))
                 
            elif proximal_y<soma_center and not 'axon' in seg_name: # soma_center used to be 0 and compared to distal_y (deprecated)
                
                # only basal dendrites
                if 'dend' in seg_name:
                    lower_basal_csas.append(compute_stem_csa(diameter))
                
                lower_csas.append(compute_stem_csa(diameter))
                
                
            else:
                pass
            
            
            if 'axon' in seg_name and proximal_y>=soma_center: # soma_center used to be 0 and compared to distal_y (deprecated)
                upper_axon_csas.append(compute_stem_csa(diameter))
            elif 'axon' in seg_name and proximal_y<soma_center: # soma_center used to be 0 and compared to distal_y (deprecated)
                lower_axon_csas.append(compute_stem_csa(diameter))
            else:
                basal_csas.append(compute_stem_csa(diameter))
            
    
    stem_diams = np.array(stem_diams)
    stem_csas = compute_stem_csa(stem_diams)
    
    ###################
    #### Feature 2 ####
    ###################
    num_stems = len(stem_csas)
    
    for stem_id,name in stems.items():
    
        prox_list, dist_list = get_segment_coord_lists({stem_id:name},morpho,stems)

        if 'axon' in name:
            
            ####################
            #### Feature 3 #####
            ####################
            axon_init_angle = compute_average_orientation(prox_list,dist_list)
            axon_init_rel_angle = compute_average_relative_orientation(prox_list,dist_list)
            axon_length = compute_total_length(prox_list,dist_list)
            axon_yi = float(prox_list[0].attrib['y'])
            axon_yf = float(dist_list[-1].attrib['y'])
        else:
            continue
            

            
    ######################        
    #### Features 4-6 ####
    ######################
    total_csa = np.sum(stem_csas)
    if len(upper_axon_csas)>0:
        upper_axon_csa = np.sum(upper_axon_csas)
    else:
        upper_axon_csa = 0.0
        
    if len(lower_axon_csas)>0:
        lower_axon_csa = np.sum(lower_axon_csas)
    else:
        lower_axon_csa = 0
        
    basal_csa = np.sum(basal_csas)
    
    upper_basal_csa = np.sum(upper_basal_csas)
    lower_basal_csa = np.sum(lower_basal_csas)
    
    upper_csa = np.sum(upper_csas)
    lower_csa = np.sum(lower_csas)
    
    # this is really just apical, but depends since a neuron has 0, 1, or 2 apical-like dendrites
    upper_non_basal_csa = upper_csa - upper_basal_csa
    lower_non_basal_csa = lower_csa - lower_basal_csa
    
    
    #### Stem specific stuff ####
    stem_groups = [s.split('_')[1]+'_'+s.split('_')[2] for s in stems.values()]


#     soma_y_ub, soma_y_lb = get_domain_extent(comp_domains,morpho,stems,
#                                             which_domain='somatic')
    axon_y_ub, axon_y_lb = get_domain_extent(comp_domains,morpho,stems,
                                            which_domain='axonal')
        
    axonal_locs, axonal_diams = get_domain_locs(comp_domains,morpho,stems,which_domain='axonal')    
    
    
    #### Find all basal stuff ####
    print('...basal stuff')
    basal_groups = get_domain_groups('basal',comp_domains)
    basal_init_groups = [group for group in basal_groups if group in stem_groups]
    
    basal_y_ub, basal_y_lb = get_domain_extent(comp_domains,morpho,stems,
                                            which_domain='basal')
    
    # only the initial compartments (multiple segments) of the stem
    basal_root_metrics = get_stem_root_metrics(comp_domains,morpho,stems,
                                      which_stem_domain='basal',
                                      soma_center=soma_center)
    
    basal_yis, basal_yfs, basal_xzis, basal_xzfs, [basal_upper_angle, basal_lower_angle, basal_upper_rel_angle, basal_lower_rel_angle, basal_upper_length, basal_lower_length, basal_upper_dist, basal_lower_dist, basal_upper_y_dist, basal_lower_y_dist] = basal_root_metrics
    
    
    
    # only the terminal segments of the domains
    basal_term_metrics = get_stem_terminal_metrics(comp_domains,morpho,stems,
                                                   which_stem_domain='basal',
                                                   soma_center=soma_center)
    
    
    basal_term_lower_yfs,basal_term_upper_yfs,basal_term_lower_xzfs,basal_term_upper_xzfs, [basal_term_lower_angle,basal_term_upper_angle,
                                                                                            basal_term_upper_rel_angle, basal_term_lower_rel_angle,
                                                                                            basal_term_lower_y_dist,basal_term_upper_y_dist] = basal_term_metrics
    num_basal_terms = len(basal_term_lower_yfs+basal_term_upper_yfs)
    
    
    basal_locs, basal_diams = get_domain_locs(comp_domains,morpho,stems,which_domain='basal')
        
        
        
    #### Find all apical stuff ####
    print('...apical stuff')
    apical_groups = get_domain_groups('apical',comp_domains)
    
    if len(apical_groups): # non-zero (has apical dends)
        apical_init_groups = [group for group in apical_groups if group in stem_groups]

        apic_y_ub, apic_y_lb = get_domain_extent(comp_domains,morpho,stems,
                                                which_domain='apical')

        # only the initial comparment (multiple segments) of the stem
        apical_root_metrics = get_stem_root_metrics(comp_domains,morpho,stems,
                                           which_stem_domain='apical',
                                           soma_center=soma_center)
        
        # only terminal segments of the domain
        apical_term_metrics = get_stem_terminal_metrics(comp_domains,morpho,stems,
                                                   which_stem_domain='apical',
                                                   soma_center=soma_center)
        
        num_apic_terms = len(apical_term_metrics[0]+apical_term_metrics[1])
        
        apical_locs, apical_diams = get_domain_locs(comp_domains,morpho,stems,which_domain='apical')
        
    else:
        apical_root_metrics = [np.nan],[np.nan],[np.nan],[np.nan],[np.nan for _ in range(10)]
        apical_term_metrics = [[np.nan],[np.nan],[np.nan],[np.nan],[np.nan for _ in range(6)]]
        num_apic_terms = 0
        apical_locs, apical_diams = [np.nan],[np.nan]
        
        
    apic_yis, apic_yfs, apic_xzis, apic_xzfs, [apic_upper_angle, apic_lower_angle,apic_upper_rel_angle,apic_lower_rel_angle,apic_upper_length, apic_lower_length, apic_upper_dist, apic_lower_dist, apic_upper_y_dist, apic_lower_y_dist] = apical_root_metrics
    
    apic_term_lower_yfs,apic_term_upper_yfs,apic_term_lower_xzfs,apic_term_upper_xzfs, [apic_term_lower_angle,apic_term_upper_angle,
                                                                                        apic_term_lower_rel_angle,apic_term_upper_rel_angle,
                                                                                        apic_term_lower_y_dist,apic_term_upper_y_dist] = apical_term_metrics
    
    
    
    # define feature vector
    features = [cell_labels[model_i],soma_root,soma_center,soma_length,soma_lb,soma_ub,
                num_stems,axon_init_angle,axon_init_rel_angle,axon_length,axon_yi,axon_yf,
                total_csa,upper_axon_csa,lower_axon_csa, upper_basal_csa, upper_non_basal_csa,lower_basal_csa,lower_non_basal_csa,
                axon_y_ub,axon_y_lb,basal_y_ub,basal_y_lb,apic_y_ub,apic_y_lb,
                axonal_locs,axonal_diams,basal_locs,basal_diams,apical_locs,apical_diams,
                
                # basal stem + term features
                num_basal_terms,basal_yis, basal_yfs, basal_xzis, basal_xzfs, 
                basal_upper_angle, basal_lower_angle, basal_upper_rel_angle,basal_lower_rel_angle,
                basal_upper_length, basal_lower_length, 
                basal_upper_dist, basal_lower_dist, basal_upper_y_dist, basal_lower_y_dist,
                basal_term_lower_yfs,basal_term_upper_yfs,basal_term_lower_xzfs,basal_term_upper_xzfs,
                basal_term_lower_angle,basal_term_upper_angle,basal_upper_rel_angle,basal_lower_rel_angle,
                basal_term_lower_y_dist,basal_term_upper_y_dist, 
                
                # apical + non-basal stem/term features
                num_apic_terms,apic_yis, apic_yfs, apic_xzis, apic_xzfs, 
                apic_upper_angle, apic_lower_angle, apic_upper_rel_angle,apic_lower_rel_angle,
                apic_upper_length, apic_lower_length, 
                apic_upper_dist, apic_lower_dist, apic_upper_y_dist, apic_lower_y_dist,
                apic_term_lower_yfs,apic_term_upper_yfs,apic_term_lower_xzfs,apic_term_upper_xzfs,
                apic_term_lower_angle,apic_term_upper_angle,apic_term_upper_rel_angle,apic_term_lower_rel_angle,
                apic_term_lower_y_dist,apic_term_upper_y_dist]
                
    
    print('...computed %s local morphometrics'%len(features))

    
    if model_i==0:
        columns = ['cell_name','soma_root','soma_center','soma_length','soma_y_lb','soma_y_ub',
                   'num_stems','axon_init_angle','axon_init_rel_angle','axon_length','axon_yi','axon_yf',
                   'total_csa','upper_axon_csa','lower_axon_csa','upper_basal_csa','upper_non_basal_csa','lower_basal_csa','lower_non_basal_csa',
                    'axon_y_ub','axon_y_lb','basal_y_ub','basal_y_lb','apic_y_ub','apic_y_lb',
                   'axon_locs','axon_diams','basal_locs','basal_diams','apical_locs','apical_diams',
                   
                    'num_basal_terms','basal_yis', 'basal_yfs', 'basal_xzis', 'basal_xzfs', 
                    'basal_upper_angle', 'basal_lower_angle','basal_upper_rel_angle','basal_lower_rel_angle',
                   'basal_upper_length', 'basal_lower_length', 
                    'basal_upper_dist','basal_lower_dist', 'basal_upper_y_dist', 'basal_lower_y_dist',
                    'basal_term_lower_yfs','basal_term_upper_yfs','basal_term_lower_xzfs','basal_term_upper_xzfs',
                   'basal_term_lower_angle','basal_term_upper_angle','basal_term_upper_rel_angle','basal_term_lower_rel_angle',
                   'basal_term_lower_y_dist','basal_term_upper_y_dist',
                   
                    'num_apic_terms','apic_yis', 'apic_yfs', 'apic_xzis', 'apic_xzfs', 
                    'apic_upper_angle','apic_lower_angle',  'apic_upper_rel_angle','apic_lower_rel_angle',
                   'apic_upper_length', 'apic_lower_length', 
                    'apic_upper_dist','apic_lower_dist', 'apic_upper_y_dist', 'apic_lower_y_dist',
                    'apic_term_lower_yfs','apic_term_upper_yfs','apic_term_lower_xzfs','apic_term_upper_xzfs',
                   'apic_term_lower_angle','apic_term_upper_angle','apic_term_upper_rel_angle','apic_term_lower_rel_angle',
                   'apic_term_lower_y_dist','apic_term_upper_y_dist']

        morphometrics_df = pd.DataFrame(columns=columns)
        
        for col, feat in zip(columns,features):
            morphometrics_df[col] = [feat]
    
    else:
        temp_df = pd.DataFrame(columns=columns)
        
        for col, feat in zip(columns,features):
            temp_df[col] = [feat]
            
        join_frames = [morphometrics_df,temp_df]
        morphometrics_df = pd.concat(join_frames,ignore_index=True)

        

# save
filename = 'final_model_morphometrics.pkl'
f = os.path.join(path_to_save,filename)
morphometrics_df.to_pickle(f)

        

            

    

    

    

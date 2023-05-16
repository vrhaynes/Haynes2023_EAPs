import numpy as np
from scipy.stats import variation, circmean, circstd

def find_radius_of_sphere(SA):
    return np.sqrt(SA)/(2.*np.sqrt(np.pi))

def is_in_interval(pt,lb,ub):
    if pt<=ub and pt>=lb:
        return True
    return False

def is_not_in_cylinder(x,y,z,SA):
    '''
        Needs to be mapped to equivalent cylinder to conform to NEURON implementations.
        Note models may have multiple soma sections which would complicate all of this.
        How did Justas compute surface area for those models?
    '''
    pass



def is_not_in_sphere(x,y,z,SA):
    
    r_bound = find_radius_of_sphere(SA)
    
    if np.sqrt(x**2 + y**2 + z**2)>r_bound: return True
    
    else: return False

def compute_segment_length(prox,dist):
    
    prox_x, prox_y, prox_z = float(prox.attrib['x']),float(prox.attrib['y']),float(prox.attrib['z'])
    dist_x, dist_y, dist_z = float(dist.attrib['x']),float(dist.attrib['y']),float(dist.attrib['z'])
    
    return np.sqrt((dist_x-prox_x)**2 + (dist_y-prox_y)**2 + (dist_z-prox_z)**2)    

def compute_stem_csa(diam):
    return np.pi*np.square(diam/2)

def get_domain_groups(domain,comp_domains):
    
    domain_groups = []
    for seg_group in comp_domains[domain]:
        try:
            domain_groups.append(seg_group.attrib['segmentGroup'])
        except KeyError:
            continue # there can be one extra for path length
        
    return domain_groups

def get_domain_segments(morpho,domain_groups,bounds=None):
    domain_segments  = {}
    
    for segment in morpho.findall('{http://www.neuroml.org/schema/neuroml2}segment'):
        # leverage standardized segment naming
        seg_name = segment.attrib['name'].split('_')
        this_group = seg_name[1]+'_'+seg_name[2]
        
        if this_group in domain_groups:
            
            if not bounds is None:
                [y_lb,y_ub] = bounds[0]
                rad_ub = bounds[1]
                
                # always the last field
                distal = segment[-1]
                x,y,z  = float(distal.attrib['x']),float(distal.attrib['y']),float(distal.attrib['z'])
                
                # project onto xz-plane
                xz = np.sqrt(x**2+z**2)
                
                # cylinder interval
                if is_in_interval(y,y_lb,y_ub) and xz<=rad_ub:
                    domain_segments.update({segment.attrib['id']:segment.attrib['name']})
                else:
                    continue  
                    
            # only do if bounds not included
            else:
                domain_segments.update({segment.attrib['id']:segment.attrib['name']})
                
    return domain_segments


def get_segment_coord_lists(domain_segments,morpho,stems,verbose=False):

    segment_parents = {}
    no_parent = {}

    # for compute domain-level props
    proximal_list = []
    distal_list = []

    for segment in morpho.findall('{http://www.neuroml.org/schema/neuroml2}segment'):
        seg_id = segment.attrib['id']
        seg_name = segment.attrib['name']


        if seg_name in domain_segments.values():

            parent_id = segment[0].attrib['segment']

            # if not stem of soma, look at parent for proximal
            if seg_name in stems.values():

                # parent is segment[0]
                proximal = segment[1]
                distal = segment[2]
                
                segment_parents.update({seg_id:distal})
            else:
                # no proximal so distal is segment[1]
                distal = segment[1]
                segment_parents.update({seg_id:distal})
                
                # sometimes parent segment is out of bounds so not included
                try:
                    proximal = segment_parents[parent_id]
            
                except KeyError:
                    no_parent.update({seg_id:{}})
                    no_parent[seg_id].update({'Name':seg_name})
                    no_parent[seg_id].update({'Parent ID':parent_id})
                    
                    # don't append anything to list
                    continue
                    
            proximal_list.append(proximal)
            distal_list.append(distal)
    if verbose:
        print('Removed %s segments from list...'%(len(domain_segments) - len(proximal_list)))
    return proximal_list, distal_list

def compute_total_length(all_prox,all_dist):
    
    prox_xs = np.zeros(len(all_prox))
    prox_ys = np.zeros(len(all_prox))
    prox_zs = np.zeros(len(all_prox))
    
    dist_xs = np.zeros(len(all_dist))
    dist_ys = np.zeros(len(all_dist))
    dist_zs = np.zeros(len(all_dist))
    
    for i, prox in enumerate(all_prox):
        prox_xs[i], prox_ys[i], prox_zs[i] = float(prox.attrib['x']),float(prox.attrib['y']),float(prox.attrib['z'])

    for i, dist in enumerate(all_dist):
        dist_xs[i], dist_ys[i], dist_zs[i] = float(dist.attrib['x']),float(dist.attrib['y']),float(dist.attrib['z'])


    total_length = 0
    for i in range(len(prox_xs)):
        total_length += np.sqrt((dist_xs[i]-prox_xs[i])**2 + (dist_ys[i]-prox_ys[i])**2 + (dist_zs[i]-prox_zs[i])**2)
        
    return total_length

        

def compute_average_orientation(all_prox,all_dist,rads=False):
    
    dist_xs = np.zeros(len(all_dist))
    dist_ys = np.zeros(len(all_dist))
    dist_zs = np.zeros(len(all_dist))
    
    for i, dist in enumerate(all_dist):
        dist_xs[i], dist_ys[i], dist_zs[i] = float(dist.attrib['x']),float(dist.attrib['y']),float(dist.attrib['z'])
        
        
    # project all segments onto XZ-plane
    dist_xzs = np.sqrt(np.square(dist_xs)**2 + np.square(dist_zs)**2)
    
    # compute signed angle in radians
    all_angles = np.arctan2(dist_ys,dist_xzs)
    
    # circular average
    mean = circmean(all_angles)
    
    if not rads:
        mean = np.degrees(mean)

    return mean

def compute_dispersion_orientation(all_prox,all_dist,rads=False):
    
    dist_xs = np.zeros(len(all_dist))
    dist_ys = np.zeros(len(all_dist))
    dist_zs = np.zeros(len(all_dist))

    for i, dist in enumerate(all_dist):
        dist_xs[i], dist_ys[i], dist_zs[i] = float(dist.attrib['x']),float(dist.attrib['y']),float(dist.attrib['z'])
        
        
    # project all segments onto XZ-plane
    dist_xzs = np.sqrt(np.square(dist_xs)**2 + np.square(dist_zs)**2)
    
    # compute signed angle in degrees
    all_angles = np.arctan2(dist_ys,dist_xzs)
    disp = circstd(all_angles)
    
    if not rads:
        disp = np.degrees(disp)
    
    return disp


def compute_average_distance(all_dist):
    
    dist_xs = np.zeros(len(all_dist))
    dist_ys = np.zeros(len(all_dist))
    dist_zs = np.zeros(len(all_dist))

    for i, dist in enumerate(all_dist):
        dist_xs[i], dist_ys[i], dist_zs[i] = float(dist.attrib['x']),float(dist.attrib['y']),float(dist.attrib['z'])
        
        
    # get location
    xyz_locations = np.sqrt(np.square(dist_xs)**2 + np.square(dist_ys)**2 + np.square(dist_zs)**2)
    
    return np.mean(xyz_locations)

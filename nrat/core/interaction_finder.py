from structure_framework import *
import geometry_engine as geo
from ..io import file_manager as fm
from ..data import reference_tools as ref

interaction_distance=15
'''
Maximum Distance for interactions
'''

def atms_centroid(atoms=list[atom]) -> list[float]:
    points=[]
    weights=[]
    for atom in atoms:
        points.append(atom.pos)
    return centroid(points)
        

def centroid(points:list[list[float]]) -> list[float]:
    return np.mean(points,axis=0)

centroid_distance_threshold:float = 6.0
plane_angle_threshold:float = 30.0
centroid_plane_distance_threshold:float = 0.75



def pi_stacking(ring1:list[atom], ring2:list[atom]):
    '''
    Determines if two aromatic rings could be in a pi-stacking interaction
    '''
    #print(ring1)
    #Centroids of each of the rings
    r1_centroid=atms_centroid(ring1)
    r2_centroid=atms_centroid(ring2)

    # Finds the distance between the centroids. This is a quick and inexpensive way to rule out rings that are obviously not interacting
    centroid_distance=geo.euclidean_distance(list(r1_centroid),list(r2_centroid)) # Distance between centroids
    #print(centroid_distance)
    if centroid_distance < centroid_distance_threshold:
  
        # Plane setup
        r1_points=[atom.pos for atom in ring1]
        r1_plane=geo.fit_plane(r1_points)
        r1_norm=r1_plane[:-1]
        
        r2_points=[atom.pos for atom in ring2]
        r2_plane=geo.fit_plane(r2_points)
        r2_norm=r2_plane[:-1]

        # Plane angle calculations
        planes_dot=np.dot(r1_norm,r2_norm)
        angle_rad=np.arccos(np.clip(planes_dot,-1.0,1.0))
        angle_deg=np.degrees(angle_rad)
        
        if angle_deg>90.0:
            angle_deg=180-angle_deg
        # Precalculates centroid-plane distances
        c1_p2=geo.point_to_plane(r1_centroid,r2_plane)
        c2_p1=geo.point_to_plane(r2_centroid,r1_plane)
        #print(angle_deg)
        if angle_deg < plane_angle_threshold: # Face-Face interaction
            # Checks to make sure that both c-p distances obey the threshold distance
            if all([distance > centroid_plane_distance_threshold for distance in (c1_p2,c2_p1)]):
                stack_type = 'face-face'
                return[(ring1,ring2),stack_type]
            else:
                return False
        else: # Edge-Face interaction
            # Checks to make sure that only one d-p distance is over the threshold. (set will contain both True and False values)
            #print(c1_p2,c2_p1)
            if len(set([distance > centroid_plane_distance_threshold for distance in (c1_p2,c2_p1)])) > 1:
                stack_type = 'edge-face'
                return[(ring1,ring2),stack_type]
            #elif(c1_p2 > centroid_plane_distance_threshold) and (c2_p1 > centroid_plane_distance_threshold): #Exeption for certain things. More what we care about is if one is significantly larger than the other.
            #    larger=max(c1_p2,c2_p1)
            #    smaller=min(c1_p2,c2_p1)
            #    if larger * 0.75 > smaller:
            #        stack_type = 'edge-face'
            #        return[(ring1,ring2),stack_type]
            else:
                return False
    else:
        return False

centroid_cation_distance_threshold:float = 5.0
cation_angle_threshold:float = None

def pi_cation(ring:list[atom], cation:atom):

    ring_centroid=atms_centroid(ring)
    centroid_cation_distance=geo.euclidean_distance(ring_centroid,cation.pos)

    if centroid_cation_distance < centroid_cation_distance_threshold:
        ring_pts = [atom.pos for atom in ring]
        ring_plane=geo.fit_plane(ring_pts)
        ring_norm=ring_plane[:-1]

        centroid_cation_vector=np.array(cation.pos)-np.array(ring_centroid)
        
        c_c_magnitude=np.linalg.norm(centroid_cation_vector)
        ring_norm_magnitude=np.linalg.norm(ring_norm)
        cation_norm_angle=np.degrees(np.arccos(np.dot(ring_norm,centroid_cation_vector)/(c_c_magnitude * ring_norm_magnitude)))

        if cation_norm_angle < cation_angle_threshold:
            return[(ring, cation),cation_norm_angle]
        else:
            return False

hydrogen_bond_distance_threshold:float=None
hydrogen_bond_angle_tolerance:float=None


def h_bonding(res: residue, ligand: ligand):
    '''
    Find hydrogen bonds between residue and ligand.
    '''
    interactions = []
    for donor_atom_name in ref.hBondDonors.get(ref.AminoacidDict[res.type], []):
        if donor_atom_name not in res.atoms:
            continue
        donor = res.atoms[donor_atom_name]
        carbon = res.atoms.get(ref.refCarbons[ref.AminoacidDict[res.type]].get(donor_atom_name))
        if not carbon:
            continue
        for latom in ligand.atomlist:
            if latom.element not in ['O', 'N', 'F']:
                continue
            d_to_a = geo.euclidean_distance(donor.pos, latom.pos)
            if d_to_a > 3.5:
                continue
            angle = geo.vector_angles(carbon.pos, donor.pos, donor.pos, latom.pos)
            ideal = ref.aproxAngles[ref.AminoacidDict[res.type]][donor_atom_name]
            if abs(angle - ideal) < 30:
                interactions.append(((donor, latom), angle))
    return interactions if interactions else False
            
            


def all_interactions(res:residue,ligand:ligand):
    interactions={}
    interactions['pi-stack']=[]
    interactions['pi-cation']=[tuple(rring,cation)]
    interactions['hbond'] = h_bonding(res, ligand) or []
    for rring in res.rings():
        for lring in ligand.rings():
            test=pi_stacking(rring,lring)
            if test:
                interactions['pi-stack'].append(test)
    #    for cation in ligand.cations():
    #        test=pi_cation(rring,cation)
    #        if test:
    #            interactions['pi-cation'].append(test)
    for cation in res.cations():
        for lring in ligand.rings():
            test=pi_cation(lring,cation)
            if test:
                interactions['pi-cation'].append(test)
    
    #interactions['hbond']=h_bonding(res,ligand)

    return interactions if interactions['pi-stack']!=[] else False
    
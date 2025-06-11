import structure_framework as sf
import networkx as nx
import numpy as np
import math

def atom_distance(atom1:sf.atom, atom2:sf.atom) -> float:
    '''
    **Finds the distance between two atom objects.**

    Uses a basic 3D Pythagorean theorem calculation:\n

    `\sqrt{(x_{2}-x_{1})^{2}+(y_{2}-y_{1})^{2}+(z_{2}-z_{1})^{2}}`
    
    '''
    xc1=atom1.x
    yc1=atom1.y
    zc1=atom1.z
    xc2=atom2.x
    yc2=atom2.y
    zc2=atom2.z

    distance=np.sqrt(((round(xc2-xc1,5))**2)+((round(yc2-yc1,5))**2)+((round(zc2-zc1,5))**2))
    return float(distance)


def res_distance(res1:sf.residue,res2:sf.residue, res1Atom='CA', res2Atom='CA') -> float:
    '''
    **Finds the distance between two residue objects.**

    Uses the CA coordinates by default, use res1/res2atom arguments to change this behavior.
    Passes desired atoms to `atom_distance`
    '''

    resAtm1=res1.atoms[res1Atom]
    resAtm2=res2.atoms[res2Atom]

    return(atom_distance(resAtm1,resAtm2))


def detect_aromatic_rings(object:sf.ligand):
    '''
    **Detects aromatic rings in ligands**
    '''
    G = nx.Graph()
    for atom in object.atomlist:
        G.add_node(atom.name)
    for i, atom1 in enumerate(object.atomlist):
        for j, atom2 in enumerate(object.atomlist):
            if i < j and atom_distance(atom1, atom2) <= 1.6:
                G.add_edge(atom1.name, atom2.name)

    cycles=nx.cycle_basis(G)
    rings=[]
    for cycle in cycles:
        coords=[[object.atoms[atom_name].x, object.atoms[atom_name].y, object.atoms[atom_name].z] for atom_name in cycle]
        if planar(coords) and len(cycle) in [5,6]:
            rings.append([object.atoms[atom_name] for atom_name in cycle])
    return rings



def planar(points:list[list[float]], tolerance=0.1) -> bool:
    '''
    **Checks to see if a set of points is coplanar**

    Default tolerance is 0.1 units, this can be changed with the `tolerance` parameter

    A set of fewer than three points will always return False, and a set of three points will always return True.
    '''

    if len(points) < 3:
        return False
    if len(points) == 3:
        return True
    
    coords=np.array(points)

    v1=coords[1]-coords[0]
    v2=coords[2]-coords[1]

    normal=np.cross(v1,v2)
    normal=normal/np.linalg.norm(normal)

    distances=np.dot(coords-coords[0], normal)

    return np.all(np.abs(distances) < tolerance)

def centroid(points:list[list[float]]) -> list[float]:
    return np.mean(points,axis=0)


def heronArea (a:float,b:float,c:float):
    s=(a+b+c)/2
    area=np.sqrt(s*(s-a)*(s-b)*(s-c))
    return area

def euclidean_distance(point1:list,point2:list):
    point1 = np.array(point1)
    point2 = np.array(point2)
    distance = np.linalg.norm(point1-point2)
    return float(distance)

def vector_angles(p1:list,p2:list,q1:list,q2:list) -> float:
    p_vector=np.array(p2)-np.array(p1)
    q_vector=np.array(q2)-np.array(q1)

    magnitudes=np.linalg.norm(p_vector) * np.linalg.norm(q_vector)
    dp=np.dot(p_vector,q_vector)
    
    cos_theta=np.clip(dp / magnitudes, -1.0, 1.0)

    theta=np.arccos(cos_theta)
    return(np.degrees(theta))

def fit_plane(points):
    pts_centroid=centroid(points)
    centered_points=points-pts_centroid

    _, _, vh = np.linalg.svd(centered_points)

    normal=vh[-1]
    A,B,C = normal
    D=-(np.dot(normal,pts_centroid))
    return A,B,C,D

def point_to_plane(point,plane) -> float:
    A, B, C, D = plane
    x, y, z = point
    numerator = abs(A * x + B * y + C * z + D)
    denominator = math.sqrt(A**2 + B**2 + C**2)
    #print(distance)
    return numerator/denominator
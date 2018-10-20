import numpy as np
import scipy
import scipy.sparse
import scipy.linalg
import pandas as pd
from scipy.spatial.distance import cdist
from collections import defaultdict

class Autocorrelation_run(object):
    """ Class which runs and holds the results of a 3D autocorrelation
    analysis of a protein.
    
    Parameters
    ----------
    protein : Protein
      Protein object to be analysed
    angles:
      Specifies whether angles are included in the stiffness matrix
    dihedrals:
      Specifies whether dihedrals are included in the stiffness matrix 

    """

    def __init__(self, resprotein):
        """
       
        """

        self.resprotein = resprotein
        self.interaction_results = []
        self.angle_results = []
        self.residue_neighbour_autocorrelation = {}
        self.residue_average_neighbour_autocorrelation = []

    
    def calculate_autocorrelation(self, include_angles=False):
        """
        bonds_out: calculate autocorrelation for bonds
        angles_out: calculate autocorrelation for angles
        """

        nresidues = len(self.resprotein.residues) 
        ninteraction = len(self.resprotein.interactions)
       

        #two-centre (bond) adjacency matrix
        A_interactions, G_interactions = generate_interaction_matrices(self.resprotein)
        K_interactions = (A_interactions.T).dot(G_interactions).dot(A_interactions)

        #three-centre (bond) adjacency matrix
        A_angles, G_angles = generate_angle_matrices(self.resprotein)
        K_angles = (A_angles.T).dot(G_angles).dot(A_angles)

    
        print('Calculating pseudoinverse of the stiffness matrix...')
        K_pinv = np.linalg.pinv(K_interactions.toarray())
        print('..done.')



        #'trick' for getting the diagonal elements only of the autocorrelation matrix
        autocorrelation = G_interactions.dot(np.multiply(((A_interactions.todense()).dot(K_pinv)), A_interactions.todense()).sum(1))
        autocorrelation = [element[0] for element in autocorrelation.tolist()]

        interaction_ids = [interaction.id for interaction in self.resprotein.interactions]

        interaction_names = [make_interaction_name(interaction) for interaction in self.resprotein.interactions]

        res1_pdbnum = [interaction.residue1.PDBnum for interaction in self.resprotein.interactions]
        res2_pdbnum = [interaction.residue2.PDBnum for interaction in self.resprotein.interactions]

        res1_chain = [interaction.residue1.chain for interaction in self.resprotein.interactions]
        res2_chain = [interaction.residue2.chain for interaction in self.resprotein.interactions]

        interaction_results_frame = pd.DataFrame(
            {
            'interaction_name' : interaction_names,
            'atom1_name': res1_pdbnum,
            'atom2_name': res2_pdbnum,
            'res1_chain' : res1_chain,
            'res2_chain' : res2_chain, 
            'autocorrelation' : autocorrelation, 
            }, 
            index=interaction_ids, 
        )
        self.interaction_results = interaction_results_frame[['atom1_name','atom2_name', 'res1_chain', 'res2_chain', 'interaction_name', 'autocorrelation']]
    

    def calculate_angle_autocorrelation(self):

        nresidues = len(self.resprotein.residues) 
        nangles = len(self.resprotein.angles)

        #two-centre (bond) adjacency matrix
        A_interactions, G_interactions = generate_interaction_matrices(self.resprotein)
        K_interactions = (A_interactions.T).dot(G_interactions).dot(A_interactions)

        #three-centre (bond) adjacency matrix
        A_angles, G_angles = generate_angle_matrices(self.resprotein)
        K_angles = (A_angles.T).dot(G_angles).dot(A_angles)

        print('Calculating pseudoinverse of the stiffness matrix...')
        K_pinv = np.linalg.pinv((K_interactions + K_angles).toarray())
        print('..done.')


        #'trick' for getting the diagonal elements only of the autocorrelation matrix
        autocorrelation = G_angles.dot(np.multiply(((A_angles.todense()).dot(K_pinv)), A_angles.todense()).sum(1))
        autocorrelation = [element[0] for element in autocorrelation.tolist()]

        angle_ids = [angle.id for angle in self.resprotein.angles]
        angle_names = [make_angle_name(angle) for angle in self.resprotein.angles]

        res1_pdbnum = [angle.residue1.PDBnum for angle in self.resprotein.angles]
        res2_pdbnum = [angle.residue2.PDBnum for angle in self.resprotein.angles]
        res3_pdbnum = [angle.residue3.PDBnum for angle in self.resprotein.angles]
        

        res1_chain = [angle.residue1.chain for angle in self.resprotein.angles]
        res2_chain = [angle.residue2.chain for angle in self.resprotein.angles]
        res3_chain = [angle.residue3.chain for angle in self.resprotein.angles]        

        angle_results_frame = pd.DataFrame(
            {
            'angle_name' : angle_names,
            'atom1_name': res1_pdbnum,
            'atom2_name': res2_pdbnum,
            'atom3_name': res3_pdbnum,            
            'res1_chain' : res1_chain,
            'res2_chain' : res2_chain,
            'res3_chain' : res3_chain,              
            'autocorrelation' : autocorrelation,     
            }, 
            index=angle_ids, 
        )
        self.angle_results = angle_results_frame[['atom1_name','atom2_name', 'atom3_name', 'res1_chain', 'res2_chain', 'res3_chain', 'angle_name', 'autocorrelation']]



    
    #not sure about these at the moment, bit arbitrary
    def calculate_residue_neighbour_autocorrelation(self, cutoff=10):

        interaction_xyz = np.array([interaction.xyz for interaction in self.resprotein.interactions])  
        res_xyz = np.array([residue.xyz for residue in self.resprotein.residues])
        dist_mat = cdist(interaction_xyz, res_xyz)

        temp = np.where(dist_mat < cutoff)
        temp2 = list(zip(temp[1],temp[0]))

        neighbour_dict = defaultdict(list)
        for k, v in temp2:
            neighbour_dict[k].append(v)

        neighbour_autocorr_dict = defaultdict(list)
        for i in range(len(self.resprotein.residues)):
            neighbour_autocorr = [self.interaction_results['autocorrelation'][j] for j in neighbour_dict[i]]
            neighbour_autocorr_dict[i] = neighbour_autocorr

        self.residue_neighbour_autocorrelation = neighbour_autocorr_dict


    def calculate_residue_average_neighbour_autocorrelation(self):
        residue_average_results = []
        for i,j in enumerate(self.resprotein.residues):
            average_autocorr = np.mean([j for j in self.residue_neighbour_autocorrelation[i]])
            residue_average_results.append(average_autocorr)

        residue_average_results_frame = pd.DataFrame(
            {
            'res_name': [residue.res_name for residue in self.resprotein.residues],
            'res_num': [residue.res_num for residue in self.resprotein.residues],
            'res_chain' : [residue.chain for residue in self.resprotein.residues], 
            'autocorr_average' : residue_average_results, 
            }, 
            index=self.resprotein.residues.id())

    
        self.residue_average_neighbour_autocorrelation = residue_average_results_frame


def generate_interaction_matrices(resprotein):
    """
    Returns the (sparse) m by n incidence matrix A for the two-centre
    interaction (i.e. the bonds) and the m by m diagonal matrix of
    force constants.
    """
    nresidues = len(resprotein.residues)
    ninteractions = len(resprotein.interactions)

    A = np.zeros([ninteractions, 3*nresidues])
    force_constants = np.zeros(ninteractions)
    for interaction in resprotein.interactions:
        
        res1_id = interaction.residue1.id
        res2_id = interaction.residue2.id

        res1_xyz = interaction.residue1.xyz
        res2_xyz = interaction.residue2.xyz

        interaction_length = np.linalg.norm(res1_xyz - res2_xyz)

        row = A[interaction.id]
        row[[3*res1_id, (3*res1_id)+1, (3*res1_id)+2]] = (res1_xyz - res2_xyz)/interaction_length
        row[[3*res2_id, (3*res2_id)+1, (3*res2_id)+2]] = (res2_xyz - res1_xyz)/interaction_length

        force_constant = interaction.force_constant 
        force_constants[interaction.id] = force_constant

    A = scipy.sparse.csr_matrix(A)
    G = scipy.sparse.diags(force_constants) 

    return (A, G)




def make_interaction_name(interaction):
    """ Returns a string giving the name of the two residues in 
    a given interaction.
    
    Parameters
    ----------
    interaction : :class:'Interaction'
      interaction object for which string is to be returned
    
    Returns
    -------
    result : string
      Name of the interaction in the form 'Res1 Name : Res2 Name'
    
    """
    return (interaction.residue1.res_name +
            interaction.residue1.res_num + ' ' + 
            interaction.residue1.chain + ' ' + 
            interaction.residue1.name + ' : ' + 
            interaction.residue2.res_name +
            interaction.residue2.res_num + ' ' + 
            interaction.residue2.chain + ' ' + 
            interaction.residue2.name)


def generate_angle_matrices(resprotein):
    """
    Returns the (sparse) m by n incidence matrix for the three-centre
    interaction (i.e. the angles) and the m by m diagonal matrix of
    force constants.
    """

    #double check maths for this to be safe (particularly signs)

    nresidues = len(resprotein.residues)
    nangles = len(resprotein.angles)

    #A = np.zeros([nangles, 3*natoms])
    A = scipy.sparse.lil_matrix((nangles, 3*nresidues))    

    force_constants = np.zeros(nangles)
    for angle in resprotein.angles:

        residue1_id = angle.residue1.id
        residue2_id = angle.residue2.id
        residue3_id = angle.residue3.id

        residue1_xyz = angle.residue1.xyz
        residue2_xyz = angle.residue2.xyz
        residue3_xyz = angle.residue3.xyz

        three_centre_length = np.linalg.norm(residue1_xyz - residue3_xyz)

        #row = A[angle.id]
        A[angle.id ,[3*residue1_id, (3*residue1_id)+1, (3*residue1_id)+2]] = (residue2_xyz - residue3_xyz)/three_centre_length
        A[angle.id ,[3*residue2_id, (3*residue2_id)+1, (3*residue2_id)+2]] = -((residue2_xyz - residue1_xyz) + (residue2_xyz - residue3_xyz))/three_centre_length
        A[angle.id ,[3*residue3_id, (3*residue3_id)+1, (3*residue3_id)+2]] = (residue2_xyz - residue1_xyz)/three_centre_length

        force_constant = angle.force_constant
        force_constants[angle.id] = force_constant
    
    A = scipy.sparse.csr_matrix(A)
    G = scipy.sparse.diags(force_constants)

    return (A,G)



def make_angle_name(angle):
    """ Returns a string giving the name of the two atoms in 
    a given bond.
    
    Parameters
    ----------
    bond : :class:'Bond'
        bond object for which string is to be returned
    
    Returns
    -------
    result : string
        Name of the bond in the form 'Atom1 Name : Atom 2 Name'
    
    """

    return (angle.residue1.res_name +
            angle.residue1.res_num + ' ' + 
            angle.residue1.chain + ' ' + 
            angle.residue1.name + ' : ' + 
            angle.residue2.res_name +
            angle.residue2.res_num + ' ' + 
            angle.residue2.chain + ' ' + 
            angle.residue2.name + ' : ' + 
            angle.residue3.res_name +
            angle.residue3.res_num + ' ' + 
            angle.residue3.chain + ' ' + 
            angle.residue3.name)  
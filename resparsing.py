import numpy as np
import residuenetwork.residues
from itertools import combinations

"""Don't really need to strip/add anything here as only parsing the C-alpha atoms anyway. """

#####################################
#                                   #
# The following functions deal with #
#      parsing of the PDB file      #
#                                   #
#####################################

class PDBParser(object):
    """Class for loading information on atoms and residues 
    from specified PDB file. 

    Parameters
    -----------
    pdb_filename : str 
      PDB ID or file name
    """

    def __init__(self, pdb_filename):
        self.pdb = pdb_filename
        self.test = None

        
    def parse(self, resprotein, ligands, bio=False, model=None, chain='all'):

        """ Takes a proteingraph3 Protein object and loads it with data 
        from the given PDB file.

        Parameters
        ----------
        bio : bool
          Indicates whether the file is a biological pdb file.
          If True, multiple models in the same file will be combined.
        model : str
          Model to be loaded.  Should be specified if there
          are multiple models in the same file (e.g. NMR structures)
        chain : tuple of str
          Tuple of chains to load.
        ligands : tuple of ligands in PDB file by chain and ID in a list of tuples  
        """

        #ligands can be a list of ligands chosen by the user.  Put all ligand atoms into a dict
        #where the keys are the (chain, ID) tuples then sort out after residue info has been loaded into a ResidueList

        with open(self.pdb) as fin:
            lines = fin.readlines()
            id_ = 0
            residues = residuenetwork.residues.ResidueList()
            ligand_dict = {key:[] for key in ligands}

            for line in lines:
                #load in residues
                if line[0:4] == 'ATOM' and line[13:15] == 'CA' and (line[21], line[22:27].strip()) not in ligands: #check this, CA could appear in preamble...
                    name = line[12:16].strip()
                    PDBnum = int(line[6:11])
                    chain = line[21]
                    res_num = line[22:27].strip()
                    res_name = line[17:20].strip() + line[26].strip()
                    
                    try:
                        bfactor = float(line[60:66])
                    except ValueError:
                        bfactor = None

                    coordinates = np.array([float(line[30:38]),
                                            float(line[38:46]),
                                            float(line[46:54])])

                    residues.append(residuenetwork.residues.Residue(id_, name, PDBnum, res_name, chain, res_num, coordinates, bfactor))
                    id_ = id_ + 1

                #ligand atoms need to be loaded individually then turned into a single bead
                elif line.startswith('ATOM') or line.startswith('HETATM'):
                    if (line[21], line[22:27].strip()) in ligands:
                        name = line[12:16].strip()
                        PDBnum = int(line[6:11])
                        chain = line[21]
                        res_num = line[22:27].strip()
                        res_name = line[17:20].strip() + line[26].strip()

                        coordinates = np.array([float(line[30:38]),
                                                float(line[38:46]),
                                                float(line[46:54])])

                        ligand_dict[(chain, res_num)].append((coordinates, res_name))
            

            for ligand in ligands:
                chain = ligand[0]
                res_num = ligand[1]

                res_name = list(set([entry[1] for entry in ligand_dict[ligand]]))[0]
                coordinates = [entry[0] for entry in ligand_dict[ligand]]

                centre_of_mass = np.array(coordinates).sum(axis=0)/len(coordinates)

                residues.append(residuenetwork.residues.Residue(id_, name, PDBnum, res_name, chain, res_num, centre_of_mass, None))
                
                id_ = id_ + 1




                



        resprotein.residues = residues


#####################################
#                                   #
# The following functions deal with #
# generation of the protein network #
# using distances                   #
#                                   #
#####################################

class Residue_Network_Generator(object):

    def __init__(self):
        self.interactions = residuenetwork.residues.InteractionList()
        self.res_list = residuenetwork.residues.ResidueList()
        self.pdb = ''

    def generate_cutoff_network(self, resprotein, cutoff_distance = 10):

        interactions = residuenetwork.residues.InteractionList()
        connected_residues = {i:[] for i,_ in enumerate(resprotein.residues)}
        id_ = 0
        for i,_ in enumerate(resprotein.residues):
            for j in range(i):
                distance = np.linalg.norm(resprotein.residues[i].xyz - resprotein.residues[j].xyz)
                if distance < cutoff_distance:
                    interaction_xyz = (resprotein.residues[i].xyz + resprotein.residues[j].xyz)/2
                    interactions.append(
                    residuenetwork.residues.Interaction(id_, resprotein.residues[i], resprotein.residues[j], interaction_xyz, 1))
                    connected_residues[i].append(j)
                    connected_residues[j].append(i)                    
                    id_ = id_ + 1

        resprotein.interactions = interactions
        for i,res in enumerate(resprotein.residues):
            res.connected_residues = connected_residues[i]     

    def generate_angle_interactions(self, resprotein):
        
        angles = residuenetwork.residues.InteractionList()
        id_ = 0
        for res in resprotein.residues:
            angle_pairs = list(combinations(res.connected_residues,2))
            for i in angle_pairs:
                angles.append(
                    residuenetwork.residues.Angle_interaction(id_, resprotein.residues[i[0]], res, resprotein.residues[i[1]], 0.1))
                id_ = id_ + 1
        
        resprotein.angles = angles

        




    def generate_expo_network(self, resprotein, cutoff_distance = 7, precision=1e-6):

        interactions = residuenetwork.residues.InteractionList()
        id_ = 0
        for i,_ in enumerate(resprotein.residues):
            for j in range(i):
                distance = np.linalg.norm(resprotein.residues[i].xyz - resprotein.residues[j].xyz)
                forceconstant = np.exp(-distance**2/cutoff_distance**2)
                if forceconstant > precision:
                    interaction_xyz = (resprotein.residues[i].xyz + resprotein.residues[j].xyz)/2
                    interactions.append(
                    residuenetwork.residues.Interaction(id_, resprotein.residues[i], resprotein.residues[j], interaction_xyz, forceconstant))
                    id_ = id_ + 1

        resprotein.interactions = interactions            



        

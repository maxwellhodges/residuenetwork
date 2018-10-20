"""Generate network that includes only the C-alpha atoms of each residue"""

import os
import numpy as np
import warnings

class ResProtein(object):
    """ Stores information about residues.
    """

    def __init__(self, **kwargs):
        self.name = []
        self.pdb_id = []
        self.interactions = []
        # dict for mapping chains to atom ids
        self.chains = []
        # dict for mapping residues to atom ids
        self.residues = []
        self.angles = []

    def __repr__(self):
        return "<ResProtein object containing {0} residues>".format(len(self.residues))


class Residue(object):

    def __init__(self, id_, name, PDBnum, res_name, chain, res_num, coordinates, bfactor, connected_residues=[]):
        self.id = id_
        self.res_name = res_name
        self.chain = chain
        self.res_num = res_num
        self.bfactor = bfactor
        #coords of C-alpha
        self.xyz = coordinates
        #PDBNum of alpha-carbon
        self.PDBnum = PDBnum
        self.name = name
        self.connected_residues = connected_residues


    def __repr__(self):
        output = ('<Residue ID: {0}.  {1} {2} in chain {3}>'
                  .format(self.id, 
                          self.res_name, 
                          self.res_num,
                          self.chain,
                          ))
        return output

class ResidueList(object):

    def __init__(self):
        self.residues = []

    def append(self, residue):
        self.residues.append(residue)

    def id(self):
        return [residue.id for residue in self.residues]
    
    def __getitem__(self, idx):
        return self.residues[idx]

    def __getitem__(self, idx):
        if isinstance(idx, int):
            return self.residues[idx] 
        elif isinstance(idx, list):
            return ResidueList([self.residues[i] for i in idx])
        elif isinstance(idx, slice):
            return ResidueList(self.residues[idx])
        
    def __len__(self):
        return len(self.residues)

    def __repr__(self):
        return "<ResidueList containing %s residues>" % len(self.residues)

    def index(self, residue):
        return self.residues.index(residue)

    def __contains__(self, residue):
        return residue in self.residues

    def __add__(self, x):
        if isinstance(x, ResidueList):
            return ResidueList(self.residues + x.residues)
        elif isinstance(x, Residue):
            new_residues = list(self.residues).append(x)
            return ResidueList(new_residues)

    def __iter__(self):
        return iter(self.residues)



class Interaction(object):
    """Create an object representing an interaction

    Parameters
    ----------
    id : int
      unique id number for the interaction
    residue1 : :proteingraph3:Residue
      first atom forming the interaction
    residue2 : :proteingraph3:Residue 
      second atom forming the interaction
    force_constant : float 
      force constant of the interaction (in J/A^2)
    interaction_type : list of str
      Types of interaction contributing to the interaction
    """

    def __init__(self, id_, residue1, residue2, interaction_xyz, force_constant):
        self.id = id_
        self.residue1 = residue1
        self.residue2 = residue2
        self.xyz = interaction_xyz
        self.force_constant = force_constant


    def __repr__(self):
        output = ('<Interaction ID: {0}.  Interaction between residues {1} and {2}>'
                  .format(self.id, self.residue1, self.residue2))
        return output


class InteractionList(object):
    """ Create a list of interactions

    Parameters
    ----------
    interactions : list of interactions
      Interactions to initialise list with.  If None list will initially be empty.
    """

    def __init__(self, interactions=None):
        if interactions:
            self.interactions = interactions
        else:
            self.interactions = []

    def append(self, interactions):
        self.interactions.append(interactions)

    def extend(self, interactions):
        self.interactions.extend(interactions)

    def id(self):
        return [interaction.id for interaction in self.interactions]

    def return_by_id(self, idx):
        """ Return a list of interactions whose id matches those
        in idx.
        
        Parameters
        ----------
        idx : list of ints 
          List of ids for the interactions you would like returned
        """
        if isinstance(idx, int):
            interaction_out = [interaction for interaction in self.interactions if interaction.id == idx]
            if len(interaction_out) > 0:
                IOError('More than one interaction in the list has the same ID')
            if len(interaction_out) == 0:
                IOError('No interaction found with that id')
            else: 
                return interaction_out[0]

        elif isinstance(idx, list):
            return InteractionList([interaction for interaction in self.interactions if interaction.id in idx])
        else:
            IOError('Input should be a single integer or a list of integers')
           
    def __getitem__(self, idx):
        if isinstance(idx, int):
            return self.interactions[idx] 
        elif isinstance(idx, list):
            return InteractionList([self.interactions[i] for i in idx])
        elif isinstance(idx, slice):
            return InteractionList(self.interactions[idx])
    
    def __len__(self):
        return len(self.interactions)

    def __repr__(self):
        return "<InteractionList containing %s interactions>" % len(self.interactions)

    def index(self, interaction):
        return self.interactions.index(interaction)

    def __contains__(self, interaction):
        return interaction in self.interactions

    def __add__(self, x):
        if isinstance(x, InteractionList):
            return InteractionList(self.interactions + x.interactions)
        elif isinstance(x, Interaction):
            new_interactions = list(self.interactions).append(x)
            return InteractionList(new_interactions)

    def __iter__(self):
        return iter(self.interactions)


    
class Angle_interaction(object):
    """Create an object representing an interaction

    Parameters
    ----------
    id : int
      unique id number for the interaction
    residue1 : :proteingraph3:Residue
      first atom forming the interaction
    residue2 : :proteingraph3:Residue 
      second atom forming the interaction
    force_constant : float 
      force constant of the interaction (in J/A^2)
    interaction_type : list of str
      Types of interaction contributing to the interaction
    """

    def __init__(self, id_, residue1, residue2, residue3, force_constant):
        self.id = id_
        self.residue1 = residue1
        self.residue2 = residue2
        self.residue3 = residue3        
        self.force_constant = force_constant


    def __repr__(self):
        output = ('<Interaction ID: {0}.  Interaction between residues {1}, {2} and {3}>'
                  .format(self.id, self.residue1, self.residue2, self.residue3))
        return output

    

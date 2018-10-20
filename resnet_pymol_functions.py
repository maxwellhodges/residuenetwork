import colorsys,sys,re
from pymol import cmd
import csv

def load_interactions(filename, prefix='interaction', mol=None):

    f = open(filename, 'r')
    reader = csv.reader(f)
    header = reader.next()
    
    # if header != ['bond_id', 'atom1_id', 'atom2_id',
    #               'bond_name', 'weight', 'distance',
    #               'propensity', 'qs', 'qs_test_set']:
    #     raise ValueError('header is not correct format')
    
    for row in reader:
        interaction_id, res1_id, res2_id, res1_chain, res2_chain = row[0], row[1], row[2], row[3], row[4]
    
        cmd.distance(prefix + interaction_id,
                    'id ' + res1_id + ' and chain ' + res1_chain,
                    'id ' + res2_id + ' and chain ' + res2_chain)
    cmd.hide('labels')
    cmd.set('dash_gap',0)

cmd.extend('load_interactions', load_interactions)


def color_interactions(filename, color_by='autocorrelation', name='interaction', 
                power=1, gradient="bwr", scale=None):

    power = float(power)    

    interactions = []

    f = open(filename, "r")
    reader = csv.reader(f)

    header = reader.next()
    
    column = [i for i, col in enumerate(header)
              if col == color_by]
    if len(column) > 1:
        raise ValueError('More than one {0} column'.format(color_by))
    else:
        column = column[0]

    for line in reader:
        interactions.append([line[0], float(line[column])])
    f.close()
    
    num_interactions = len(interactions)
        
    scores = [interaction[1] for interaction in interactions]
    if scale:
        score_min = float(scale[0])
        score_max = float(scale[1])
    else:
        score_max = max(scores)
        score_min = min(scores)

    for interaction in interactions:
        id_ = interaction[0]
        score = interaction[1]

        scaled_score = (score - score_min)/(score_max - score_min)
    
        curve1 = round((2**power)*(scaled_score**power), 5)
        curve2 = round((2**power)*((1-scaled_score)**power), 5)
        
        rgb = [min(1, curve1), 
               min(curve1, curve2),
               min(curve2, 1)]
    
        cmd.set_color('b_color' + name + id_, rgb)
        cmd.color('b_color' + name + id_, name + id_)


def show_interactions(filename, threshold, show_by='autocorrelation', name='interaction', 
               mol=None, hide = 'no', residues='no'):

    threshold = float(threshold)

    if hide == 'yes':
        cmd.hide('dashes')
        cmd.hide('sticks')

    f =  open(filename, "r")
    reader = csv.reader(f)
    
    header = reader.next()
    
    interaction_list = []
    for line in reader:
        interaction_list.append(
            [line[0], line[1], line[2],
             float(line[6])]
        )
    f.close()

    for interaction in interaction_list:
        if interaction[3] >= threshold:
            cmd.show('dashes', name + str(interaction[0]))


def load_angles(filename, prefix='angle', mol=None):

    f = open(filename, 'r')
    reader = csv.reader(f)
    header = reader.next()
    
    # if header != ['bond_id', 'atom1_id', 'atom2_id',
    #               'bond_name', 'weight', 'distance',
    #               'propensity', 'qs', 'qs_test_set']:
    #     raise ValueError('header is not correct format')
    
    for row in reader:
        interaction_id, res1_id, res2_id, res3_id, res1_chain, res2_chain, res3_chain = row[0], row[1], row[2], row[3], row[4], row[5], row[6]
    
        cmd.distance(prefix + interaction_id,
                    'id ' + res1_id + ' and chain ' + res1_chain,
                    'id ' + res3_id + ' and chain ' + res3_chain)
    cmd.hide('labels')
    cmd.set('dash_gap',0)


def show_angles(filename, threshold, show_by='autocorrelation', name='angle', 
               mol=None, hide = 'no', residues='no'):

    threshold = float(threshold)

    if hide == 'yes':
        cmd.hide('dashes')
        cmd.hide('sticks')

    f =  open(filename, "r")
    reader = csv.reader(f)
    
    header = reader.next()
    
    angle_list = []
    for line in reader:
        angle_list.append(
            [line[0], line[1], line[3],
             float(line[8])]
        )
    f.close()

    for angle in angle_list:
        if angle[3] >= threshold:
            cmd.show('dashes', name + str(angle[0])



cmd.extend('load_interactions', load_interactions)
cmd.extend('show_interactions',show_interactions)
cmd.extend('color_interactions',color_interactions)
cmd.extend('load_angles', load_angles)
cmd.extend('show_angles', show_angles)

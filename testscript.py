#import pandas as pd
import elasticnetwork as en

protein = en.molecules.Protein()

parser = en.parsing.PDBParser('3orz_mm1.pdb')
parser.parse(protein, strip={'res_name': ['HOH']})

ggenerator = pg.parsing.FIRST_Network_Generator()
ggenerator.generate_elastic_network(pdk1, angles=False, dihedrals=False)


source_residues = [('1','A')] 
propensity_calcs = pg.bondbond.BondBond_run(pdk1, source_residues)
propensity_calcs.calculate_3D_propensities()
(propensity_calcs.bond_results).to_csv('pdk1_bondbond_bonds.csv')
#propensity_calcs.calculate_bond_quantile_scores()


pdk1_autocorr = pg.autocorrelation3D.Autocorrelation3D_run(pdk1, angles=False, dihedrals=False)
pdk1_autocorr.calculate_autocorrelation()
(pdk1_autocorr.bond_results).to_csv('pdk1_autocorr_bonds.csv')

infrig_run = pg.infrig.Infrig_run(pdk1)
infrig_run.calculate_infrig(angles=False, dihedrals=False)
(infrig_run.results).to_csv('pdk1_infrig_bonds')



"""Residue network """

#import pandas as pd
import elasticnetwork as en

resprotein = en.residues.ResProtein()

ligands = [('A','1'),('A','360')]

parser = en.resparsing.PDBParser('test_pdb.pdb')
parser.parse(resprotein, ligands)

ggenerator = en.resparsing.Residue_Network_Generator()
ggenerator.generate_cutoff_network(resprotein)

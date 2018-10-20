import residuenetwork as rn

"""with ligands """

resprotein = rn.residues.ResProtein()

#A 360 - BI4 (active)  A 1 - 2A2 (allosteric)
ligands = [('A','1'),('A','360')]

parser = rn.resparsing.PDBParser('3orz_mm1.pdb')
parser.parse(resprotein, ligands)

ggenerator = rn.resparsing.Residue_Network_Generator()
ggenerator.generate_cutoff_network(resprotein)

autocorrelation_run = rn.resautocorrelation.Autocorrelation_run(resprotein)
autocorrelation_run.calculate_autocorrelation()

(autocorrelation_run.interaction_results).to_csv('autocorr_testing_3orz_w_ligand.csv')


""" without ligands """

resprotein = rn.residues.ResProtein()
ligands = []

parser = rn.resparsing.PDBParser('3orz_mm1.pdb')
parser.parse(resprotein, ligands)

ggenerator = rn.resparsing.Residue_Network_Generator()
ggenerator.generate_cutoff_network(resprotein)

autocorrelation_run = rn.resautocorrelation.Autocorrelation_run(resprotein)
autocorrelation_run.calculate_autocorrelation()

(autocorrelation_run.interaction_results).to_csv('autocorr_testing_3orz.csv')
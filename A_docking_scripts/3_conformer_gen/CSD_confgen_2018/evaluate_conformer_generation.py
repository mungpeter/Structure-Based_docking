"""
    Evaluate a test set of structures, starting from 3D coordinates
    This script expects two input folders containing matching sets (by filename) of
    expected 3D geometries and input 3D geometries.

    The input geometry is used to generate conformers and then the conformers are compared
    to the expected geometries.

    By default the code generates results coinsistent with the definitions used in this publication:

      Freely Available Conformer Generation Methods: How Good Are They?

      Jean-Paul Ebejer,Garrett M. Morris,and Charlotte M. Dean
      J.Chem. Inf. Model. 2012, 52(5), pp 1146-1158
      http://pubs.acs.org/doi/abs/10.1021/ci2004658

   The script expects 2 folders as input, each containing sets of mol2 files. One set is
   the set of 3D input geometries, generated however you choose, the other is the set of
   expected geometries.

   The code will pair these files up and for all pairs run the conformer generator and write results to
   a CSV file to the output stream specified.

   The results contain a field expressing the RMSD of the first and the best conformer, rankings where a conformer of
   a specified threshold occured first and a count of the number of rotatable bonds (consistent with the definition used in
   the paper above, note.)

   The results can be pasted into an associated xlsx file (see results_template.xlsx) to generate various visual representations of
   the results in excel.

"""

from __future__ import division, absolute_import, print_function

from ccdc import conformer as conformer_module
from ccdc import io
from ccdc.descriptors import MolecularDescriptors

import string
import sys
import argparse
import os
import csv
import time

try:
    from ccdc_internal.molecule import Molecule
except:
    # 1D->3D will not be available
    from ccdc.molecule import Molecule

from ccdc.io import MoleculeReader

from ccdc.io import MoleculeWriter

from ccdc.utilities import Logger

from rotatable_bond_calculations import RotatableBondCalculator

try:
    import etkdg_generator
    from rdkit.Chem import AllChem
    etkdg_available = True
except:
    etkdg_available = False


global global_logger
global_logger = Logger()

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class CGEvaluation:
    """
       Class to run conformer generation evaluations.

       Takes a list of numerical thresholds as a parameter and will find the ranks of the first conformer
       that falls within that threshold
    """
    def __init__(self, rotatable_bond_source, thresholds = (0.5,1.0,1.5,2.0)):

        self.rotatable_bond_calculator = RotatableBondCalculator(rotatable_bond_source)
        self.thresholds = thresholds

    def _molecule(self,filename,no_assign_bond_types):
        """  Implementation detail - Get the first molecule from a file
        """
        m = MoleculeReader(filename)
        mol = m[0]
        mol.normalise_labels()
        if not no_assign_bond_types:
            mol.assign_bond_types()
        return mol

    def _create_input_map(self, input_folder, expected_folder):
        """  Implementation detail - pair up files from the input folder and the expected folder
        """
        def create_input_map_of_mol2_files(folder, files):
            themap = {}
            for file in files:
                basename = os.path.basename(file)
                parts = os.path.splitext(basename)
                if parts[1] == ".mol2" or parts[1] == ".mol":
                    themap[parts[0]] = os.path.abspath(os.path.join(folder,file))
            return themap

        files_in_input    = create_input_map_of_mol2_files(input_folder,os.listdir(input_folder))
        files_in_expected = create_input_map_of_mol2_files(expected_folder,os.listdir(expected_folder))

        input_expected_map = {}
        for key in sorted(files_in_input.keys()):
            if files_in_expected.has_key(key):
                input_expected_map[key] = (files_in_input[key],files_in_expected[key])

        return input_expected_map

    def rdkit_rmsd(self, conformer, expected):
        try:
            rd_conformer = etkdg_generator.ccdc_to_rdkit(conformer)
            rd_expected = etkdg_generator.ccdc_to_rdkit(expected)
            rmsd = AllChem.GetBestRMS(rd_expected, rd_conformer)
            aligned_conf = etkdg_generator.rdkit_conformer_to_ccdc(rd_conformer, -1)
            return (aligned_conf, rmsd)
        except:
            return MolecularDescriptors.overlay_rmsd_and_rmsd_tanimoto(conformer,expected,rotate_torsions=False)

    def evaluate_molecule(self, molid, input_molecule, expected_molecule, results_folder, cg_settings, save_best, save_all, generator, rmsd_method):
        """
           Runs the conformer generator on the input molecule and then reports statistics
           for the resultant conformers as compared to the expected molecule.

           returns a tuple containing overlay results for the first conformer generated,
        """

        #global_logger.set_ccdc_log_level(1)

        # Assess Input RMSD

        input_overlay                = MolecularDescriptors.overlay_rmsd_and_rmsd_tanimoto(input_molecule, expected_molecule, rotate_torsions=True)
        input_overlay_with_inversion = MolecularDescriptors.overlay_rmsd_and_rmsd_tanimoto(input_molecule, expected_molecule, invert = True, rotate_torsions=True)

        without_inversion_input_rmsd = input_overlay[1]
        best_input_rmsd = min(without_inversion_input_rmsd,input_overlay_with_inversion[1])

        settings = conformer_module.ConformerSettings()
        if cg_settings is not None:
            settings = cg_settings
        print('Generating', molid)
        start = time.time()
        if generator == 'etkdg':
            if not etkdg_available:
                raise runtime_error('etkdg is not available')
            conformer_generator = etkdg_generator.EtkdgGenerator(settings)
        else:
            conformer_generator = conformer_module.ConformerGenerator(settings)
        conformers = conformer_generator.generate(input_molecule)
        end = time.time()
        gentime = end - start
        print(molid, 'Took', gentime)

        #global_logger.set_ccdc_log_level(0)

        conformer_number = 0
        rank_of_best = 0
        best_overlay_result  = None
        best_overlay_molecule = None
        first_overlay_result = None
        rank_of_first_with_rmsd_less_than = {}
        for threshold in self.thresholds:
            rank_of_first_with_rmsd_less_than[threshold] = -1

        ensemble_size = len(conformers)

        for conformer in conformers:
            conformer_number = conformer_number + 1

            if rmsd_method == "rdkit":
                overlay = self.rdkit_rmsd(conformer.molecule,expected_molecule)
            else:
                overlay = MolecularDescriptors.overlay_rmsd_and_rmsd_tanimoto(conformer.molecule,expected_molecule,rotate_torsions=False)

            if conformer_number == 1:
                first_overlay_result = overlay

            if best_overlay_result is None or overlay[1] < best_overlay_result[1]:
                best_overlay_result = overlay
                rank_of_best = conformer_number - 1
                best_overlay_molecule = conformer.molecule

            for key in rank_of_first_with_rmsd_less_than.keys():
                if rank_of_first_with_rmsd_less_than[key] == -1:
                    if overlay[1] < key:
                        rank_of_first_with_rmsd_less_than[key] = conformer_number

        if results_folder is not None:
          if save_best:
            filename = os.path.join(results_folder,"%s_best.mol2" % molid)
            with MoleculeWriter(filename) as mol_writer:
                mol_writer.write(conformers[rank_of_best].molecule)

          if save_all:
            filename = os.path.join(results_folder,"%s_all.mol2" % molid)
            with MoleculeWriter(filename) as mol_writer:
                for conf in conformers:
                    mol_writer.write(conf.molecule)

        input_with_best = MolecularDescriptors.overlay_rmsd_and_rmsd_tanimoto(best_overlay_molecule,input_molecule,rotate_torsions=False)
        probably_input = input_with_best[1] < 0.001

        degrees_of_freedom = conformers.n_flexible_rings_in_molecule + conformers.n_rotamers_in_molecule

        return {
            'first_overlay_result': first_overlay_result,
            'best_overlay_result': best_overlay_result,
            'rank_of_first_with_rmsd_less_than': rank_of_first_with_rmsd_less_than,
            'ensemble_size': ensemble_size,
            'gentime': gentime,
            'probably_input': probably_input,
            'without_inversion_input_rmsd': without_inversion_input_rmsd,
            'best_input_rmsd': best_input_rmsd,
            'rank_of_best': rank_of_best+1,
            'degrees_of_freedom': degrees_of_freedom
            }


    def evaluate_library(self, input_folder, expected_folder, results_folder, results_stream, cg_settings, save_best, save_all, add_hydrogens, generator, no_assign_bond_types, slice_start, slice_end, rmsd_method):
        """
           Runs the conformer generator on a pair of folders each containing matching pairs of file names.
           Consolidates the individual results and writes them out to the results stream.
        """

        summary_stats = {}
        first = True
        best_rmsds = []

        input_expected_pairs = self._create_input_map(input_folder,expected_folder)

        keys = sorted(input_expected_pairs.keys())
        if slice_end == -1:
            slice_end = len(keys)
        keys = keys[slice_start:slice_end]
        for key in keys:
            try:
                print('running', key)
                input    = self._molecule(input_expected_pairs[key][0], no_assign_bond_types)
                expected = self._molecule(input_expected_pairs[key][1], no_assign_bond_types)
                if add_hydrogens:
                    input.add_hydrogens()
                    expected.add_hydrogens()
                data = self.evaluate_molecule(key,input,expected,results_folder,cg_settings,save_best,save_all, generator, rmsd_method)

                if first == True:
                    rankheader = ""
                    for rankkey in sorted(data['rank_of_first_with_rmsd_less_than'].keys()):
                        rankheader = rankheader + ",rank of first conformer with RMSD < %s" % str(rankkey)

                    results_stream.write('id,rmsd of first,rmsd of best,number of rotatable bonds (rdkit strict)%s,ensemble size,generation time,probably input,as read input rmsd, as read input rmsd with inversion allowed, rank of best, degrees of freedom\n' % rankheader)
                    first = False

                rankdata = ""
                for rankkey in sorted(data['rank_of_first_with_rmsd_less_than'].keys()):
                    rankdata = rankdata +",%s" % str(data['rank_of_first_with_rmsd_less_than'][rankkey])
                first_rmsd = data['first_overlay_result'][1]
                best_rmsd  = data['best_overlay_result'][1]
                ensemble_size = data['ensemble_size']
                gentime = data['gentime']
                probably_input = data['probably_input']
                as_read_input_rmsd = data['without_inversion_input_rmsd']
                best_input_rmsd = data['best_input_rmsd']
                best_rank = data['rank_of_best']
                degrees_of_freedom = data['degrees_of_freedom']
                nrot = self.rotatable_bond_calculator(input)

                best_rmsds.append(best_rmsd)
                line = string.join( [ str(x) for x in [key,first_rmsd,best_rmsd,nrot] ] , "," ) + \
                rankdata + "," + str(ensemble_size) + "," + \
                str(gentime) + "," + str(probably_input) + "," + str(as_read_input_rmsd)+ "," + str(best_input_rmsd)+ \
                "," + str(best_rank) + "," + str(degrees_of_freedom) + "\n"
                results_stream.write(line)
                sys.stdout.flush()
                sys.stderr.flush()
            except:
                print("exception, skipping", key)


    def generate_ccdc_start_points(self, id_smilescode_pairs, expected_folder, output_folder, results_stream, no_assign_bond_types):
        """
            Some simple code to use our 1D->3D methods to generate start points
        """
        if output_folder is not None:
          if not os.path.exists(output_folder):
            os.mkdir(output_folder)

        first = True
        for pair in id_smilescode_pairs:
            molid = pair[0]
            smiles = pair[1]
            mol = Molecule.from_smiles(smiles, molid, True)

            if output_folder is not None:
                output_file = os.path.join(output_folder, '%s.mol2' % molid)
                with MoleculeWriter(output_file) as mol_writer:
                    mol_writer.write(mol)

            if expected_folder is not None:
                expected_mol_filename = os.path.join(expected_folder,"%s.mol2" % molid)
                expected_molecule = self._molecule(expected_mol_filename, no_assign_bond_types)
                overlay = MolecularDescriptors.overlay_rmsd_and_rmsd_tanimoto(mol,expected_molecule,rotate_torsions=True)
                if first:
                    results_stream.write('id,rotated torsion rmsd,rotated torsion rmsd tanimoto\n')
                    first = False

                results_stream.write(string.join([molid,str(overlay[1]),str(overlay[2])],',') + '\n')

    def generate_ccdc_start_points_from_expected(self, expected_folder, output_folder, results_stream, add_hydrogens, no_assign_bond_types):
        """
            Some simple code to use our 1D->3D methods to generate start points
        """
        if output_folder is not None:
          if not os.path.exists(output_folder):
            os.mkdir(output_folder)

        first = True
        for f in os.listdir(expected_folder):
           if os.path.splitext(f)[1] == ".mol2":
                expected_mol_filename = os.path.join(expected_folder,f)
                output_mol_filename = os.path.join(output_folder,os.path.basename(f))

                input_molecule = self._molecule(expected_mol_filename, no_assign_bond_types)
                if add_hydrogens:
                    input_molecule.add_hydrogens()
                molid = input_molecule.identifier

                mol = Molecule.from_molecule(input_molecule)

                with MoleculeWriter(output_mol_filename) as mol_writer:
                    mol_writer.write(mol)

                expected_molecule = self._molecule(expected_mol_filename, no_assign_bond_types)
                overlay = MolecularDescriptors.overlay_rmsd_and_rmsd_tanimoto(mol,expected_molecule,rotate_torsions=True)
                if first:
                    results_stream.write('id,rotated torsion rmsd,rotated torsion rmsd tanimoto\n')
                    first = False

                results_stream.write(string.join([molid,str(overlay[1]),str(overlay[2])],',') + '\n')

    def read_pairs(self,csv_file):
        pairs = []
        with open(csv_file,"r") as csv_stream:
           reader = csv.reader(csv_stream)
           for row in reader:
              pairs.append( (row[0],row[1]) )
        return pairs

    def generate_ccdc_start_points_from_csv_file(self, csv_file, expected_folder, output_folder, results_stream, no_assign_bond_types):
        pairs = self.read_pairs(csv_file)
        self.generate_ccdc_start_points(pairs, expected_folder, output_folder, results_stream, no_assign_bond_types)


def parse_arguments(arguments):
    def program_mode(args):

        if args.input_folder is not None and args.expected_folder is not None:
            return "conformers"
        elif args.smiles_file is not None:
            return "1Dto3D"
        elif args.expected_folder is not None and args.expected_folder != args.results_folder:
            return "Mol2GenMol"

        raise InputError("Unable to establish mode of action for the script")

    parser = argparse.ArgumentParser( description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter )

    parser.add_argument( '-i', '--input_folder', default = None, help='folder containing mol2 files of 3D input files for conformer generation', required=False )
    parser.add_argument( '-e', '--expected_folder', default = None, help='folder containing mol2 files of 3D geometries actually observed', required=False )

    parser.add_argument( '-s', '--smiles_file', default = None, help='For 1D->3D transformation - a CSV file in the format <identifier>,<smiles code>', required=False )

    parser.add_argument( '-o', '--output_file', default = None, help='output CSV file (default writes to stdout) to write any comparison results to', required=False )

    results_folder_help = \
""" folder to write the molecule results to. If this folder if specified, when generating conformers,
    the best conformer will be written to this folder.

    When performing 1D->3D transformation the final 3D coordinates for the output will be written to this file.

    Default is '.'
"""

    parser.add_argument( '-r', '--results_folder',  default = '.', help=results_folder_help,required=False )

    parser.add_argument( '-sb', '--save_best',  action='store_true',help="whether to save the best structure when doing conformer generation", required=False )
    parser.add_argument( '-sa', '--save_all',  action='store_true',help="whether to save all structures generated when doing conformer generation", required=False )

    parser.add_argument( '-n', '--number_of_conformers', type=int, default = -1, help='Number of conformers to generate for each molecule (-1 means the default) as specified by the version of the API in use', required=False )

    parser.add_argument( '-rbf', '--rotatable_bond_count_source',  default = 'rdkit_strict',
    help='either a csv file containing counts of rotatable bonds and associated molecule identifiers for each entry or one of pybel,rdkit_strict,api or moe',
    required=False )
    parser.add_argument('--add_hydrogens', action="store_true", default=False)
    parser.add_argument('--generator', default='ccdc', help='the generator to use, either "ccdc" or "etkdg"', required=False)
    parser.add_argument('--no_assign_bond_types', action="store_true", default=False)

    parser.add_argument( '-ss', '--slice_start', type=int, default = 0, help='start of subset slice', required=False )
    parser.add_argument( '-se', '--slice_end', type=int, default = -1, help='end of subset slice', required=False )
    parser.add_argument('--rmsd', default='ccdc', help='the rmsd calculator to use, either "ccdc" or "rdkit"', required=False)

    args = parser.parse_args(arguments)

    try:
       mode = program_mode(args)
    except InputError as e:
       parser.error(e.message)

    return args,mode


def main(arguments):

    args,mode = parse_arguments(arguments)

    out = sys.stdout
    if args.output_file is not None:
        out = open(args.output_file,'w')

    evaluator = CGEvaluation(rotatable_bond_source=args.rotatable_bond_count_source)
    print(mode)
    if mode == "conformers":
        settings = conformer_module.ConformerSettings()
        if args.number_of_conformers != -1:
            settings.max_conformers = args.number_of_conformers
        evaluator.evaluate_library(args.input_folder, args.expected_folder, args.results_folder, out, settings, args.save_best, args.save_all, args.add_hydrogens, args.generator, args.no_assign_bond_types, args.slice_start, args.slice_end, args.rmsd)
    elif mode == "Mol2GenMol":
        evaluator.generate_ccdc_start_points_from_expected(args.expected_folder, args.results_folder, out, args.add_hydrogens, args.no_assign_bond_types)
    else:
        evaluator.generate_ccdc_start_points_from_csv_file(args.smiles_file, args.expected_folder, args.results_folder, out, args.no_assign_bond_types)

if __name__ == "__main__":
    main(sys.argv[1:])

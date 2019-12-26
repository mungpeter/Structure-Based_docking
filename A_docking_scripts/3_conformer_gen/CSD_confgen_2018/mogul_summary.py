from __future__ import print_function

import sys
import os
import glob
from ccdc import io
from ccdc import conformer

def mogul_summary(filename, engine):
    mr = io.MoleculeReader(filename)
    mol = mr[0]
    gmol = engine.analyse_molecule(mol)
    name = os.path.splitext(os.path.basename(filename))[0]
    out = [name]
    for geom in (gmol.analysed_angles, gmol.analysed_bonds, gmol.analysed_rings, gmol.analysed_torsions):
        out.append(len(geom))
        out.append(len([x for x in geom if x.unusual]))
    print(",".join([str(x) for x in out]))


def main():
    if len(sys.argv) != 2:
        print("mogul_summary <dir_of_molecules>")
        exit(1)
    
    print("name,n_angles,n_unusual_angles,n_bonds,n_unusual_bonds,n_rings,n_unusual_rings,n_torsions,n_unusual_torsions")

    engine = conformer.GeometryAnalyser()

    files = glob.glob(os.path.join(sys.argv[1], "*.mol2"))
    for filename in files:
        mogul_summary(filename, engine)

if __name__=="__main__":
    main()

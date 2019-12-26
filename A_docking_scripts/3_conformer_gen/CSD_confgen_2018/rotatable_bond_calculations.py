from ccdc.search import SMARTSSubstructure
from ccdc.search import SubstructureSearch
from ccdc import molecule
import csv
import os

class RotatableBondCalculator:
  def __init__(self, mode):
      self.mode = mode
      if mode not in ['pybel','rdkit_strict','api','moe']:
        if os.path.exists(mode):
           self.mode = "file"
           self.data = {}
           with open(mode,'r') as rbfile:
              reader = csv.DictReader(rbfile)
              for row in reader:
                 self.data[row['id']] = row['Number of Rotatable Bonds']
        else:
          raise ValueError("Mode must be one of 'pybel','rdkit_strict','api','moe' or an existing csv file")

  def pybel_is_rotatable(self,bond):
      """
          This is a function that establishes whether a bond is rotatable according to the conventions
          used in Pybel. It goes beyond the standard CCDC definition of a rotatable bond to say that terminal groups such as
          methyls or OH groups are not rotatable. Further it does not regard single bonds bound to an sp carbon
          as rotatable.

          It reproduces the distribution seen in this paper:

                Freely Available Conformer Generation Methods: How Good Are They?
                Jean-Paul Ebejer, Garrett M. Morris, and Charlotte M. Dean
                http://pubs.acs.org/doi/abs/10.1021/ci2004658

          Which uses Pybel for its work.
      """
      def is_terminal_h_group(atom,other_atom):
          for x in atom.neighbours:
              if x.index != other_atom.index:
                  if x.atomic_number != 1:
                      return False
          return True

      def atom_is_bound_to_triple_bond(atom):
          for bond in atom.bonds:
              if bond.bond_type == 3:
                  return False

      if bond.is_rotatable:
          if is_terminal_h_group(bond.atoms[0],bond.atoms[1]):
              return False
          if is_terminal_h_group(bond.atoms[1],bond.atoms[0]):
              return False
          if atom_is_bound_to_triple_bond(bond.atoms[0]):
              return False
          if atom_is_bound_to_triple_bond(bond.atoms[1]):
              return False
          return True

      return False

  def ester_bond_count(self,molecule):
      ss = SMARTSSubstructure("[#6][O,S]!@C(=[O,S])[#6]")
      substructure_search = SubstructureSearch()
      substructure_search.add_substructure(ss)

      hits = substructure_search.search(molecule)

      subset = {}
      for hit in hits:
         hit_atoms = hit.match_atoms()
         if not subset.has_key(hit_atoms[1].label + hit_atoms[2].label):
             subset[hit_atoms[1].label + hit_atoms[2].label] = 1

      return len(subset.keys())


  def amide_bond_count(self,molecule):
      ss = SMARTSSubstructure("[!H][NX3]!@C(=[O,S])")
      substructure_search = SubstructureSearch()
      substructure_search.add_substructure(ss)

      hits = substructure_search.search(molecule)
      subset = {}
      for hit in hits:
         hit_atoms = hit.match_atoms()
         if not subset.has_key(hit_atoms[1].label + hit_atoms[2].label):
             subset[hit_atoms[1].label + hit_atoms[2].label] = 1

      return len(subset.keys())

  def moe_is_rotatable(self,bond):
      if bond.is_cyclic:
         return False

      # Assumes all H atoms are present: strict definition should
      # be all Hatoms to fill valence.

      h_atom1 = len([a for atom in bond.atoms[0].neighbours if atom.atomic_number == 1])
      d_atom1  = len([a for atom in bond.atoms[0].neighbours if atom.atomic_number != 1])
      h_atom2 = len([a for atom in bond.atoms[1].neighbours if atom.atomic_number == 1])
      d_atom2  = len([a for atom in bond.atoms[1].neighbours if atom.atomic_number != 1])
      if h_atom1 + d_atom1 < 2:
          return False
      if h_atom2 + d_atom2 < 2:
          return False

      return True

  def __call__(self,molecule):
      """
          return the count of rotatble bonds according to the definition expressed in the function pybel_is_rotatable
      """
      if self.mode == "pybel":
          return len([b for b in molecule.bonds if self.pybel_is_rotatable(b) ])
      elif self.mode == "rdkit_strict":
          return len([b for b in molecule.bonds if self.pybel_is_rotatable(b)]) - self.ester_bond_count(molecule) - self.amide_bond_count(molecule)
      elif self.mode == "moe":
          return len([b for b in molecule.bonds if self.moe_is_rotatable(b) ])
      elif self.mode == "api":
          return len([b for b in molecule.bonds if b.is_rotatable])
      elif self.mode == "file":
          return self.data[molecule.identifier]

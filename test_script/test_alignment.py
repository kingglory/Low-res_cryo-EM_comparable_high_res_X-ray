from __future__ import absolute_import, division, print_function
import time
import mmtbx
import iotbx.pdb
import iotbx.phil
import iotbx.pdb.fetch
from libtbx.utils import null_out
from hydrogen_bond_restraints import prepare_hydrogen_restraints
from hydrogen_bond_restraints import as_pymol
from libtbx import easy_run
from phenix.programs import homology
import iotbx.bioinformatics.pdb_info

def align_tow_chain(hierarchy_A,hierarchy_B):
  for chain_A in hierarchy_A.chains():
    for chain_B in hierarchy_B.chains():
      se = chain_A.as_padded_sequence()
      sx = chain_B.as_padded_sequence()
      #print ("\n",se,"\n",sx)
      align = mmtbx.alignment.align(seq_a = se,
                                    seq_b = sx)
      alignment = align.extract_alignment()
      #print (dir(alignment))
      (i_seqs,j_seqs) = alignment.exact_match_selections()
      match_X_E = {}
      for (i,j) in  zip(i_seqs,j_seqs):
        print (i,j)
        match_X_E[j+1]=i+1
      print (match_X_E)

if __name__ == '__main__':
  pdb_inp_B = iotbx.pdb.input("B.pdb")
  pdb_inp_A = iotbx.pdb.input("A.pdb")
  h_A = homology.get_hierarchy(pdb_inp_A)
  h_B = homology.get_hierarchy(pdb_inp_B)
  model = mmtbx.model.manager(
    model_input=pdb_inp_B,
    process_input=True,
    log=null_out())
  align_tow_chain(h_A,h_B)
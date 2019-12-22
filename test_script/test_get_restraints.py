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


def test_0():
  '''
  easy_run.call("phenix.pdbtools {0} remove_alt_confs=True".format("A.pdb"))
  easy_run.call("phenix.pdbtools {0} remove_alt_confs=True".format("B.pdb"))
  pdb_inp_B = iotbx.pdb.input("B.pdb_modified.pdb")
  pdb_inp_A = iotbx.pdb.input("A.pdb_modified.pdb")
  '''
  pdb_inp_B = iotbx.pdb.input("B.pdb")
  pdb_inp_A = iotbx.pdb.input("A.pdb")
  h_A = homology.get_hierarchy(pdb_inp_A)
  h_B = homology.get_hierarchy(pdb_inp_B)
  model = mmtbx.model.manager(
    model_input=pdb_inp_A,
    process_input=True,
    log=null_out())
  as_pymol(model, prefix="A")
  prepare_hydrogen_restraints(hierarchy=h_B,hierarchy_pair=h_A,pdb_id="B")


if __name__ == '__main__':
  start = time.time()
  test_0()
  end = time.time()
  time_cost = (end - start)
  print("it cost % seconds" % time_cost)




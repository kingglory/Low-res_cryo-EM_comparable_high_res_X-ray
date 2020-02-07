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

  pdb_inp_B = iotbx.pdb.input("B.pdb")
  pdb_inp_A = iotbx.pdb.input("A.pdb")
  h_A = homology.get_hierarchy(pdb_inp_A)
  h_B = homology.get_hierarchy(pdb_inp_B)
  model = mmtbx.model.manager(
    model_input=pdb_inp_B,
    process_input=True,
    log=null_out())
  #as_pymol(model, prefix="B0")
  prepare_hydrogen_restraints(hierarchy=h_A,hierarchy_pair=h_B,pdb_id="B0")


def test_1():

  pdb_inp_B = iotbx.pdb.input("B1.pdb")
  pdb_inp_A = iotbx.pdb.input("A1.pdb")
  h_A = homology.get_hierarchy(pdb_inp_A)
  h_B = homology.get_hierarchy(pdb_inp_B)
  model = mmtbx.model.manager(
    model_input=pdb_inp_B,
    process_input=True,
    log=null_out())
  #as_pymol(model, prefix="B1")
  prepare_hydrogen_restraints(hierarchy=h_A,hierarchy_pair=h_B,pdb_id="B1")


def test_2():

  pdb_inp_B = iotbx.pdb.input("B2.pdb")
  pdb_inp_A = iotbx.pdb.input("A2.pdb")
  h_A = homology.get_hierarchy(pdb_inp_A)
  h_B = homology.get_hierarchy(pdb_inp_B)
  model = mmtbx.model.manager(
    model_input=pdb_inp_B,
    process_input=True,
    log=null_out())
  #as_pymol(model, prefix="B2")
  prepare_hydrogen_restraints(hierarchy=h_A,hierarchy_pair=h_B,pdb_id="B2")

if __name__ == '__main__':
  start = time.time()
  test_0()
  test_1()
  test_2()
  end = time.time()
  time_cost = (end - start)
  print("it cost % seconds" % time_cost)




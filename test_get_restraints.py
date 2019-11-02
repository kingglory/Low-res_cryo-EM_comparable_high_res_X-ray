from __future__ import absolute_import, division, print_function
import os
import time
import mmtbx
import iotbx.pdb
import iotbx.phil
import mmtbx.model
import mmtbx.alignment
import iotbx.pdb.fetch
from hydrogen_bond_restraints import prepare_hydrogen_restraints
from mmtbx.nci import hbond
from libtbx import easy_run
from libtbx import easy_pickle
from libtbx.utils import Sorry
from libtbx.utils import null_out
from phenix.programs import homology
import iotbx.bioinformatics.pdb_info
from libtbx.test_utils import approx_equal
from mmtbx.utils import run_reduce_with_timeout


def test_0():
  easy_run.call("phenix.pdbtools {0} remove_alt_confs=True".format("A.pdb"))
  easy_run.call("phenix.pdbtools {0} remove_alt_confs=True".format("B.pdb"))
  pdb_inp_B = iotbx.pdb.input("B.pdb_modified.pdb")
  pdb_inp_A = iotbx.pdb.input("A.pdb_modified.pdb")
  h_A = homology.get_hierarchy(pdb_inp_A)
  h_B = homology.get_hierarchy(pdb_inp_B)
  prepare_hydrogen_restraints(hierarchy=h_A,hierarchy_pair=h_B,pdb_id="A")


if __name__ == '__main__':
  start = time.time()
  test_0()
  end = time.time()
  time_cost = (end - start)
  print("it cost % seconds" % time_cost)




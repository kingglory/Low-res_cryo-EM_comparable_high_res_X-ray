from __future__ import absolute_import, division, print_function
import os
import time
import mmtbx
import logging
import traceback
import iotbx.pdb
import iotbx.phil
import mmtbx.model
import mmtbx.alignment
import iotbx.pdb.fetch
from mmtbx.nci import hbond
from libtbx import easy_run
from libtbx import easy_pickle
from libtbx.utils import Sorry
from libtbx.utils import null_out
from phenix.programs import homology
import iotbx.bioinformatics.pdb_info
from libtbx.easy_mp import pool_map
from libtbx.test_utils import approx_equal
from mmtbx.utils import run_reduce_with_timeout

aa_codes = [
"resname ALA",
"resname ARG",
"resname ASN",
"resname ASP",
"resname CYS",
"resname GLN",
"resname GLU",
"resname GLY",
"resname HIS",
"resname ILE",
"resname LEU",
"resname LYS",
"resname MET",
"resname MSE",
"resname PHE",
"resname PRO",
"resname SER",
"resname THR",
"resname TRP",
"resname TYR",
"resname VAL"
]



def get_hierarchy(pdb_inp):
  hierarchy = pdb_inp.construct_hierarchy()
  asc = hierarchy.atom_selection_cache()
  ss = " or ".join(aa_codes)
  ss = "(%s) and not (resname UNX or resname UNK or resname UNL)"%ss
  sel = asc.selection(ss)
  hierarchy = hierarchy.select(sel)
  xray_structure = hierarchy.extract_xray_structure()
  xray_structure.convert_to_isotropic()
  hierarchy.adopt_xray_structure(xray_structure)
  if(hierarchy.atoms().size()==0):
    return None
  if(len(hierarchy.models())>1):
    return None
  return hierarchy

def run(file_name):
  pdb_code = os.path.basename(file_name)[3:7]
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  data_type = pdb_inp.get_experiment_type()
  resolution = pdb_inp.resolution()
  hierarchy = get_hierarchy(pdb_inp)
  list_p=None
  if hierarchy is not None:
    rr = run_reduce_with_timeout(
      stdin_lines=hierarchy.as_pdb_string().splitlines(),
      file_name=None,  # "model.pdb.gz",
      parameters="-oh -his -flip -keep -allalt -pen9999 -",
      override_auto_timeout_with=None)
    pdb_inp = iotbx.pdb.input(source_info=None, lines=rr.stdout_lines)
    hierarchy = pdb_inp.construct_hierarchy()
    model = mmtbx.model.manager(
          model_input = None,
          pdb_hierarchy = hierarchy,
          crystal_symmetry = pdb_inp.crystal_symmetry(),
          process_input = True,
          log = null_out())
    result = mmtbx.nci.hbond.find(model=model).result
    for r in result:
      list_p = [r.d_AD,r.d_HA,r.a_YAH,r.a_DHA,data_type,resolution,pdb_code]
  return list_p


def foo(f):
  list_p =None
  pdb_inp = iotbx.pdb.input(file_name=f)
  try:
    list_p = run(f)
  except Exception as e :
    errorFile = open('log.txt','a')
    errorFile.write(f + '\n')
    errorFile.write(traceback.format_exc())
    errorFile.close()
  return list_p



if __name__ == '__main__':
  logging.basicConfig(filename='log.txt',level=logging.DEBUG,
            format='%(asctime)s - %(levelname)s - %(message)s')
  path = '/home/wensong/PDB/pdb/'
  of = open("".join([path,"INDEX"]),"r")
  files = ["".join([path,f]).strip() for f in of.readlines()]
  for f in files:
    if not os.path.exists(f):
      files.remove(f)
  of.close()
  i=0
  list_pdbcode_result = pool_map(
    processes=10,
    func=foo,
    iterable=files)
  easy_pickle.dump("hydrogen_bond_statistics.pkl",list_pdbcode_result)





    
 

















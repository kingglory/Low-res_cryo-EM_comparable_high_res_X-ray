from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from phenix.programs import homology
from mmtbx.nci import hbond
import time
import mmtbx
import mmtbx.alignment
from libtbx import easy_run
import iotbx.pdb
import iotbx.pdb.fetch
import iotbx.phil
import iotbx.bioinformatics.pdb_info
from mmtbx.utils import run_reduce_with_timeout

def add_hydrogen_for_hierarchy(hierarchy_old):
  rr = run_reduce_with_timeout(
    stdin_lines=hierarchy_old.as_pdb_string().splitlines(),
    file_name=None,  # "model.pdb.gz",
    parameters="-oh -his -flip -keep -allalt -pen9999 -",
    override_auto_timeout_with=None)
  pdb_inp = iotbx.pdb.input(source_info=None, lines=rr.stdout_lines)
  hierarchy_new = pdb_inp.construct_hierarchy()
  return hierarchy_new

def prepare_hydrogen_restraints(pdb_file):
  params = homology.get_default_params()
  params.num_of_best_pdb = 1
  pdb_code_ref = pdb_file[0:4]
  res = homology.file_perfect_pair(pdb_file, params)
  phil_result= {}
  for r in res :
    data_type_A = False
    #print (r.chain_ref,r.match[0].chain_id,r.match[0])
    pdb_code_ref = pdb_file[0:4]
    pdb_code = r.match[0].pdb_code.lower()
    chain_id = r.match[0].chain_id
    match_info = pdb_code_ref+"_"+r.chain_ref+"_"+pdb_code+"_"+chain_id
    #print (match_info)
    #easy_run.call("phenix.fetch_pdb {0}".format(pdb_code))
    # there are two conformers in chian B 233 5kdo.pdb,
    # so cann't keep only one conformer situation in hierarchy
    easy_run.call("phenix.pdbtools {0} remove_alt_confs=True".format(pdb_code+".pdb"))
    pdb_inp = iotbx.pdb.input(pdb_code+".pdb_modified.pdb")
    data_type_X =pdb_inp.get_experiment_type()
    if data_type_X == False:continue
    he = iotbx.pdb.input(file_name=pdb_file).construct_hierarchy()
    hx = iotbx.pdb.input(file_name=pdb_code + ".pdb_modified.pdb").construct_hierarchy()
    sel_e = he.atom_selection_cache().selection("protein and chain %s"%(r.chain_ref))
    he = he.select(sel_e)
    sel_x = hx.atom_selection_cache().selection("protein and chain %s" % (chain_id))
    hx = hx.select(sel_x)
    chain_E = (he.only_chain())
    chain_X = (hx.only_chain())
    hx_h = add_hydrogen_for_hierarchy(hx)
    phil_obj = align_tow_chain(chain_E,chain_X,hx_h.as_pdb_string())
    phil_result[match_info] = phil_obj
  return phil_result

def dump_phil_file(phil_result,for_phenix_refine=True):
  top = """refinement{
    geometry_restraints.edits{
      %s
    }
  }
      """
  for match_info, phil_obj in phil_result.items():
    top_str = top % phil_obj
    file_name = match_info + ".eff"
    with open(file_name, 'w') as fileobject:
      fileobject.write(top_str)

def align_tow_chain(chain_E,chain_X,str_chain_X,
              distance_ideal=None, sigma_dist=0.1,angle_ideal = None,
              sigma_angle=2, use_actual=True):
  se = chain_E.as_padded_sequence()
  sx = chain_X.as_padded_sequence()
  align = mmtbx.alignment.align(
    seq_a = se,
    seq_b = sx)
  alignment = align.extract_alignment()
  chain_E_id = chain_E.id
  chain_X_id = chain_X.id
  (i_seqs,j_seqs) = alignment.exact_match_selections()
  match_result = {}
  for (i,j) in  zip(i_seqs,j_seqs):
    match_result[j+1]=i+1
  reses_E = chain_E.residues()
  reses_X = chain_X.residues()
  #print (reses_X[0].resseq,reses_E[0].resseq)
  pdb_inp = iotbx.pdb.input(source_info=None, lines=str_chain_X)
  model = mmtbx.model.manager(
    model_input=pdb_inp,
    process_input=True,
    log=null_out())
  result = mmtbx.nci.hbond.find(model=model).get_result()
  f = "chain %s and resseq %s and name %s"
  top = "a"
  for r in result:
    print (type(r.atom_A.resseq),type(match_result[6]))
    h = f % (chain_E_id, match_result[int(r.atom_H.resseq.strip())], r.atom_H.name)
    a = f % (chain_E_id, match_result[int(r.atom_A.resseq.strip())], r.atom_A.name)
    d = f % (chain_E_id, match_result[int(r.atom_D.resseq.strip())], r.atom_D.name)
    print (h,a,d)
    if (not use_actual):
      if (r.d_HA < 2.5):
        dt = 2.05
      else:
        dt = 2.8
      if (r.a_DHA < 130):
        at = 115
      else:
        at = 160
    else:
      dt = r.d_HA
      at = r.a_DHA
    dis = """bond {
            atom_selection_1 = %s
            atom_selection_2 = %s
            symmetry_operation = %s
            distance_ideal = %f
            sigma = 0.05
            }
    """ % (h, a, str(r.symop), dt)
    if (str(r.symop) != "x,y,z"): continue
    ang = """angle {
            atom_selection_1 = %s
            atom_selection_2 = %s
            atom_selection_3 = %s
            angle_ideal = %f
            sigma = 5
            }
    """ % (a, h, d, at)
    base = dis + ang
    top = top + base
  top_final = top[1:]
  dump_phil_file(top_final)


    #return phil_obj
if __name__ == '__main__':
    start = time.time()
    prepare_hydrogen_restraints('6d9h.pdb')
    #dump_phil_file('6d9h.pdb')
    end = time.time()
    time_cost = (end - start)
    print ("it cost % seconds" % time_cost)

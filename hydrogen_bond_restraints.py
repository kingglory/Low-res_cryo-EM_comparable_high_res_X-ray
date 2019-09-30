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

def prepare_hydrogen_restraints(pdb_file,file_name="hbond.eff"):
  params = homology.get_default_params()
  params.num_of_best_pdb = 1
  pdb_code_ref = pdb_file[0:4]
  res = homology.file_perfect_pair(pdb_file, params)
  for r in res :
    data_type_A = False
    #print (r.chain_ref,r.match[0].chain_id,r.match[0])
    pdb_code_ref = pdb_file[0:4]
    pdb_code = r.match[0].pdb_code
    chain_id = r.match[0].chain_id
    # easy_run.call("phenix.fetch_pdb {0}".format(pdb_code))
    pdb_inp = iotbx.pdb.input(pdb_code+".pdb")
    data_type_A =pdb_inp.get_experiment_type()
    if data_type_A == False:continue
    he = iotbx.pdb.input(file_name=pdb_file).construct_hierarchy()
    hx = iotbx.pdb.input(file_name=pdb_code + ".pdb").construct_hierarchy()
    sel_e = he.atom_selection_cache().selection("protein and chain %s"%(r.chain_ref))
    he = he.select(sel_e)
    sel_x = hx.atom_selection_cache().selection("protein and chain %s" % (chain_id))
    hx = hx.select(sel_x)
    chain_A = (he.only_chain())
    chain_B = (hx.only_chain())
    hx_h = add_hydrogen_for_hierarchy(hx)
    #align_tow_chain(chain_A, chain_B,hx_h.as_pdb_string())
    phil_obj = align_tow_chain(chain_A,chain_B,hx_h.as_pdb_string())
    '''TO DO:
       replace chian A atoms informations by B atoms informations (res_id,chian_id)
    '''
    with open(file_name, "w") as of:
      print("geometry_restraints.edits {", file=of)
      

def align_tow_chain(chain_A,chain_B,str_chain,
              distance_ideal=None, sigma_dist=0.1,angle_ideal = None,
              sigma_angle=2, use_actual=True):
  se = chain_A.as_padded_sequence()
  sx = chain_B.as_padded_sequence()
  align = mmtbx.alignment.align(
    seq_a = se,
    seq_b = sx)
  alignment = align.extract_alignment()
  matches = alignment.identity_matches()
  a_seq = alignment.a
  b_seq = alignment.b
  ai_seq = alignment.i_seqs_a
  bi_seq = alignment.i_seqs_b
  '''
  print (dir(alignment))
  print(dir(matches))
  print (a_seq)
  print (matches)
  print (b_seq)
  print (ai_seq)
  print (bi_seq)
  '''
  match_id=[]
  i=1
  for m in matches:
    if m == "|":
        match_id.append(i)
    i = i + 1
  #print (match_id)
  pdb_inp = iotbx.pdb.input(source_info=None, lines=str_chain)
  model = mmtbx.model.manager(
    model_input=pdb_inp,
    process_input=True,
    log=null_out())
  result = mmtbx.nci.hbond.find(model=model).get_result()
  f = "chain %s and resseq %s and name %s"
  top = "a"
  for r in result:
    h = f % (r.atom_H.chain, r.atom_H.resseq, r.atom_H.name)
    a = f % (r.atom_A.chain, r.atom_A.resseq, r.atom_A.name)
    d = f % (r.atom_D.chain, r.atom_D.resseq, r.atom_D.name)
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
  top_final.replace(chain_A.id,chain_B.id)
  print (top_final)
  return top_final


   #return phil_obj
if __name__ == '__main__':
    start = time.time()
    prepare_hydrogen_restraints('6d9h.pdb')
    end = time.time()
    time_cost = (end - start)
    print ("it cost % seconds" % time_cost)

from __future__ import absolute_import, division, print_function
import os
import time
import mmtbx
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

def add_hydrogen_for_hierarchy(hierarchy_old):
  xray_structure = hierarchy_old.extract_xray_structure()
  xray_structure.convert_to_isotropic()
  hierarchy_old.adopt_xray_structure(xray_structure)
  rr = run_reduce_with_timeout(
    stdin_lines=hierarchy_old.as_pdb_string().splitlines(),
    file_name=None,  # "model.pdb.gz",
    parameters="-oh -his -flip -keep -allalt -pen9999 -",
    override_auto_timeout_with=None)
  pdb_inp = iotbx.pdb.input(source_info=None, lines=rr.stdout_lines)
  hierarchy_new = pdb_inp.construct_hierarchy()
  return hierarchy_new

def prepare_hydrogen_restraints(hierarchy,hierarchy_pair=None,pdb_id=None,
                                fetch_pdb=True,dump_phil_files=True,SS_percent=True):
  params = homology.get_default_params()
  params.num_of_best_pdb = 1
  if hierarchy_pair == None:
    #print (dir(hierarchy))
    he_h = add_hydrogen_for_hierarchy(hierarchy)
    atoms = he_h.atoms()
    atom_total = atoms.size()
    phils = "a"
    num_sum = 0
    for chain in hierarchy.chains():
      if not chain.is_protein(): continue
      sequence = chain.as_padded_sequence()
      res = homology.perfect_pair(sequence, params)
      phil_result= {}
      if res is None:continue
      for r in res :
        pdb_code = r.pdb_code.lower()
        chain_id = r.chain_id.lower()
        #match_info = pdb_id+"_"+chain.id.lower()+"_"+pdb_code+"_"+chain_id
        if fetch_pdb == True:
          easy_run.call("phenix.fetch_pdb {0}".format(pdb_code))
          pdb_inp = iotbx.pdb.input(pdb_code + ".pdb")
          data_type_X = pdb_inp.get_experiment_type()
          data_resolution = pdb_inp.resolution()
          if data_resolution > 2.0:continue
          if data_type_X != "X-RAY DIFFRACTION": continue
          # there are two conformers in chian B 233 5kdo.pdb,
          # so cann't keep only one conformer situation in hierarchy
          easy_run.call("phenix.pdbtools {0} remove_alt_confs=True".format(pdb_code + ".pdb"))
          pdb_inp = iotbx.pdb.input(pdb_code + ".pdb_modified.pdb")
        else:
          pdb_file_local = "/home/wangwensong/PDB/pdb/%s"%pdb_code[1:3] + "/pdb%s"%pdb_code + ".ent.gz"
          #print (pdb_file_local)
          #if not os.path.exists(pdb_file_local): continue
          pdb_inp = iotbx.pdb.input(pdb_file_local)
          data_type_X = pdb_inp.get_experiment_type()
          data_resolution = pdb_inp.resolution()
          if data_resolution > 2.0:continue
          #print (data_resolution,"8"*50)
          #if data_type_X != "X-RAY DIFFRACTION": continue
          easy_run.call("phenix.pdbtools {0} remove_alt_confs=True".format(pdb_file_local))
          pdb_inp = iotbx.pdb.input("pdb"+pdb_code + ".ent.gz_modified.pdb")
        hx = pdb_inp.construct_hierarchy()
        ss = " or ".join(aa_codes)
        ss = "(%s) and not (resname UNX or resname UNK or resname UNL)" % ss
        sel_x = hx.atom_selection_cache().selection("%s and chain %s" % (ss,chain_id))
        hx = hx.select(sel_x)
        #print (len(hx.chains))
        hx.write_pdb_file(file_name="test.pdb")
        chain_E = chain
        for chain in hx.chains():
          chain_X = chain
        hx_h = add_hydrogen_for_hierarchy(hx)
        (phil_obj, num_sub) = align_tow_chain(chain_E,chain_X,hx_h.as_pdb_string())
        if phil_obj is None: continue
        #if num_sub == 0: continue
        num_sum = num_sum + num_sub
        phils = phils + phil_obj
        easy_run.call("rm -f {0}".format("*"+pdb_code + "*"))
        #phil_result[match_info] = phil_obj
    if dump_phil_files==True:
      phils=phils[1:]
      dump_phil_file(phils,pdb_id)
    if SS_percent==True:
      p = (float('%.4f'%(float(num_sum) /atom_total)))*100
  if hierarchy_pair is not None:
    he_h = add_hydrogen_for_hierarchy(hierarchy)
    atoms = he_h.atoms()
    atom_total = atoms.size()
    for chain in hierarchy.chains():
      if not chain.is_protein(): continue
      chain_E = chain
      ss = " or ".join(aa_codes)
      ss = "(%s) and not (resname UNX or resname UNK or resname UNL)" % ss
      sel_x = hierarchy_pair .atom_selection_cache().selection(ss)
      hx = hierarchy_pair .select(sel_x)
      chain_X = (hx.only_chain())
      hx_h = add_hydrogen_for_hierarchy(hx)
      (phil_obj, num_sub) = align_tow_chain(chain_E, chain_X, hx_h.as_pdb_string())
      if phil_obj is None: continue
      if num_sub == 0: continue
      dump_phil_file(phil_obj, pdb_id="A")
      p = (float('%.4f' % (float(num_sub) / atom_total))) * 100
      print("\nRestraint file is B.eff; ")
      print("Percentage of B pdb atoms that have a match in A pdb file is :{0}".format(p))
  easy_run.call("rm -f {0}".format("*modified*"))
  if os.path.exists("hbond.eff"):
    easy_run.call("rm -r {0}".format("hbond.eff"))
  if os.path.exists("myprotein.fasta"):
    easy_run.call("rm -r {0}".format("myprotein.fasta"))
  return (p)

def dump_phil_file(phils,pdb_id,for_phenix_refine=True):
  top = """refinement{
  geometry_restraints.edits{
    %s
      }
    }
        """
  file_name = pdb_id + ".eff"
  top_str = top % phils
  with open(file_name, 'w') as fileobject:
    fileobject.write(top_str)

def as_pymol(model,prefix=None):
  pdb_file_name = "%s.pdb"%prefix
  with open(pdb_file_name, "w") as of:
    print(model.model_as_pdb(), file=of)
  with open("%s.pml"%prefix, "w") as of:
    print("load", "/".join([os.getcwd(), pdb_file_name]), file=of)
    for r in mmtbx.nci.hbond.find(model=model).result:
      if(str(r.symop) != "x,y,z"): continue
      ai = r.atom_H
      aj = r.atom_A
      one = "chain %s and resi %s and name %s and alt '%s'"%(
        ai.chain, ai.resseq, ai.name, ai.altloc)
      two = "chain %s and resi %s and name %s and alt '%s'"%(
        aj.chain, aj.resseq, aj.name, aj.altloc)
      print("dist %s, %s"%(one, two), file=of)
'''
def dump_phil_files(phil_result,for_phenix_refine=True):
  top = """refinement{
    geometry_restraints.edits{
    }
  }
      """
  for match_info, phil_obj in phil_result.items():
    phils = phil_obj+phils
    top_str = top % phil_obj
    file_name = match_info + ".eff"
    with open(file_name, 'w') as fileobject:
      fileobject.write(top_str)
'''
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
  match_X_E = {}
  for (i,j) in  zip(i_seqs,j_seqs):
    match_X_E[j+1]=i+1
  reses_E = chain_E.residues()
  reses_X = chain_X.residues()
  pdb_inp = iotbx.pdb.input(source_info=None, lines=str_chain_X)
  model = mmtbx.model.manager(
    model_input=pdb_inp,
    process_input=True,
    log=null_out())
  #print ((mmtbx.nci.hbond.find(model=model).result))
  result = mmtbx.nci.hbond.find(model=model).result
  num_sub = 0
  for r in reses_E:
    if (int(r.resseq.strip())-1) in i_seqs:
      num_sub = num_sub + r.atoms_size()
  f = "chain %s and resseq %s and name %s"
  top = "a"
  for r in result:
    num_h = int(r.atom_H.resseq.strip())
    num_a = int(r.atom_A.resseq.strip())
    num_d = int(r.atom_D.resseq.strip())
    if num_h not in match_X_E.keys(): continue
    if num_a not in match_X_E.keys(): continue
    if num_d not in match_X_E.keys(): continue
    h = f % (chain_E_id, match_X_E[int(r.atom_H.resseq.strip())], r.atom_H.name)
    a = f % (chain_E_id, match_X_E[int(r.atom_A.resseq.strip())], r.atom_A.name)
    d = f % (chain_E_id, match_X_E[int(r.atom_D.resseq.strip())], r.atom_D.name)
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
  phil_obj = top[1:]
  #print (phil_obj)
  return (phil_obj,num_sub)

if __name__ == '__main__':
    ### test past

    start = time.time()
    #test_files = ["new_try_pdbs.pdb","5k7l.pdb","5vms.pdb","5vm7.pdb","6gdg.pdb","5kmg.pdb","5owx.pdb"]
    #files = ["5m6s.pdb","5nbz.pdb","6h3l.pdb","6h3n.pdb","6n2p.pdb","6o1k.pdb","6o1l.pdb","6o1m.pdb"]
    #files =["6h3l.pdb"]
    files = ["5kmg.pdb","5ni1.pdb","5njg.pdb","6e7c.pdb","6gk3.pdb",
             "6n7w.pdb","6b0l.pdb","6b0i.pdb","6cvn.pdb","5vm7.pdb"]
    dic_per = {}
    for pdb_file in files:
      print (pdb_file,"*@"*60)
      pdb_inp = iotbx.pdb.input(pdb_file)
      hierarchy = homology.get_hierarchy(pdb_inp)
      p = prepare_hydrogen_restraints(hierarchy,pdb_id=pdb_file[:4])
      dic_per[pdb_file[:4]] = p
    easy_pickle.dump("ss-percent.pkl", dic_per)
    end = time.time()
    time_cost = (end - start)
    print ("it cost % seconds" % time_cost)
    percent = easy_pickle.load("ss-percent.pkl")
    for key,value in percent.items():
      print (key,value)


    '''
    start = time.time()
    dic_per = {}
    files = easy_pickle.load("cryo-EM_3_6_filename.pkl")
    i=0
    for pdb_file in files:
      i = i+1
      print(i,"   ",pdb_file, "*" * 50)
      pdb_id = os.path.basename(pdb_file)[3:7]
      pdb_inp = iotbx.pdb.input(pdb_file)
      hierarchy = homology.get_hierarchy(pdb_inp)
      p=False      
      try:
        p = prepare_hydrogen_restraints(hierarchy, pdb_id=pdb_id)
        if p:
          dic_per[pdb_file[:4]] = p
      except Exception,e:
        print (e)
    easy_pickle.dump("ss-percent_3_6.pkl", dic_per)
    percent = easy_pickle.load("ss-percent_3_6.pkl")
    for key, value in percent.items():
      print(key, value)
    end = time.time()
    time_cost = (end - start)
    print("it cost % seconds" % time_cost)
    '''





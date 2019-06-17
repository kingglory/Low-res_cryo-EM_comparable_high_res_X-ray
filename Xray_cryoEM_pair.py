from __future__ import division
import os 
import mmtbx.model
import sys
import iotbx.pdb
from libtbx import group_args
from libtbx import easy_pickle

from iotbx.phil import process_command_line_with_files
import iotbx.bioinformatics.pdb_info
from iotbx.bioinformatics import local_blast
from iotbx.bioinformatics.structure import summarize_blast_output
from iotbx.pdb.fetch import fetch
from libtbx.easy_mp import pool_map
master_params_str = """\
model_name = None
  .type = path
  .multiple = False
  .help = Model file name
high_res = 1.5
  .type = float
  .help = High resolution pdb limit
low_res = 3.5
  .type = float
  .help = Low resolution pdb limit
identity = 100
  .type = int
  .help = sequence identity
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def get_inputs(args, log, master_params):
  cmdline = process_command_line_with_files(
    args         = args,
    master_phil  = master_params
  )
  params = cmdline.work.extract()
  return params

def get_perfect_pair(model,file_path, chain_ids=None):
  pdb_info = iotbx.bioinformatics.pdb_info.pdb_info_local()
  pdb_id = os.path.basename(file_path).strip("pdb")[:4].upper()
  h = model.get_hierarchy()
  results =  []
  for chain in h.only_model().chains():
    if not chain.is_protein():
      continue
    if chain_ids is not None and chain.id not in chain_ids():
      continue
    sequence = chain.as_padded_sequence()
    l_blast = local_blast.pdbaa(seq=sequence)
    blast_xml_result = l_blast.run()
    blast_summary = summarize_blast_output("\n".join(blast_xml_result))
    pdb_ids_to_study = {}
    for hit in blast_summary:
      #print (hit.all_ids, hit.chain_id, hit.evalue,hit.hit_num,hit.identity,hit.length,hit.pdb_id)
      #hit.show(out=sys.stdout)
      if hit.identity < 100:
        continue
      pdb_ids_to_study[hit.pdb_id] = hit.chain_id 
    #print pdb_ids_to_study 
    info_list = pdb_info.get_info_list(pdb_ids_to_study.keys())
    info_list.sort(key=lambda tup: tup[1])
    #print info_list
    if info_list :
      best_pdb_id = info_list[0][0]
      best_pdb_chain = pdb_ids_to_study[info_list[0][0]]
      #print "Best pdb:", info_list[0], "chain:", pdb_ids_to_study[info_list[0][0]]
      if info_list[0][1] < params.high_res:
        result = group_args(
          low_res_pdb  = pdb_id+'_'+chain.id ,
          high_res_pdb = best_pdb_id +'_'+ best_pdb_chain)
        if result is not None:
          results.append(result)         
  return results

def run(file_path):
   result = None
   if os.path.isfile(file_path):
      try:
        model = mmtbx.model.manager(model_input=iotbx.pdb.input(file_path))
        rs = get_perfect_pair(model,file_path)
        if rs:
          for r in rs:
            if r is not None:
              result=r
      except Exception, e:
        print "%s: %s"%(tup[0], e)        
   return result

def filte_result(tup):
  #tup[0] = pdb_id+'_'+chain.id
  #tup[1] = best_pdb_id+'_'+best_pdb_chain
  work_dir="/home/wensong/pdb_bob/"
  tup1=None 
  if tup is not None:
    if tup[0] is not None:
      if tup[1] is not None:
        file_name_1 = "pdb"+tup[0][:4].lower()+".ent.gz"
        file_path_1 = work_dir + file_name_1
        file_name_2 = "pdb"+tup[1][:4].lower()+".ent.gz"
        file_path_2 = work_dir + file_name_2
        pdb_inp_1 = iotbx.pdb.input(file_path_1)
        pdb_inp_2 = iotbx.pdb.input(file_path_2)
        res_1 = pdb_inp_1.resolution()
        res_2 = pdb_inp_2.resolution()
        data_type_1 = pdb_inp_1.get_experiment_type()
        data_type_2 = pdb_inp_2.get_experiment_type()
        if "ELECTRON MICROSCOPY" in data_type_1:
          if 6.0>= res_1 >=3.0:
            if "X-RAY DIFFRACTION" in data_type_2:
              if 0 <= res_2 <= 2.0:
                tup1 = (tup[0][:4],tup[0][4:],res_1,tup[1][:4],tup[1][4:],res_2)
        if "ELECTRON MICROSCOPY" in data_type_2:
          if 6.0>=res_2>=3.0:
            if "X-RAY DIFFRACTION" in data_type_1:
              if 0<=res_1<=2.0:
                 tup1 = (tup[1][:4],tup[1][4:],res_2,tup[0][:4],tup[0][4:],res_1)
  return tup1 
if __name__ == '__main__':
    files=[]
    work_dir = "/home/wensong/pdb_bob/" 
    pdb_info = easy_pickle.load(file_name='pdb_dict.pickle')
    #print pdb_info.items(),len(pdb_info)
    for tup in pdb_info.items():
      if tup is None:continue
      if tup[0] is None:continue 
      if tup[1] is None:continue
      file_name = "pdb"+tup[0].lower()+".ent.gz"
      file_path = os.path.join(work_dir,tup[0].lower()[1:3],file_name)
     # run(file_path)
      files.append(file_path)
   
    pair_list = pool_map(
      processes=60,
      func = run,
      iterable=files)
    
    X_RAY_CYRO_EM_PAIR= pool_map(
      processes=60,
      func = filte_result,
      iterable = pair_list)
    easy_pickle.dump("X_RAY_CYRO_EM_PAIR.pkl",X_RAY_CYRO_EM_PAIR)



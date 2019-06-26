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

def get_perfect_pair(res,data_type,model,file_path, chain_ids=None):
  pdb_info = iotbx.bioinformatics.pdb_info.pdb_info_local()
  pdb_id = os.path.basename(file_path).strip("pdb")[:4].lower()
  pdb_inp = iotbx.pdb.input(file_path)
  h = model.get_hierarchy()
  low_res = res
  tup = None
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
      print file_path,(chain.id)
      hit.show(out=sys.stdout)
      if hit.identity < 95: continue
      pdb_ids_to_study[hit.pdb_id] = hit.chain_id
      print pdb_ids_to_study,"*"*39
      pdb_pair_id = hit.pdb_id.lower()
      pdb_pair_chain_id = hit.chain_id
      work_dir = "/home/wensong/pdb_bob/"
      pdb_path = work_dir +str(pdb_pair_id[1:3])+"/pdb"+ pdb_pair_id + ".ent.gz"
      pdb_inp = iotbx.pdb.input(pdb_path)
      high_res = pdb_inp.resolution()
      pdb_pair_data_type = pdb_inp.get_experiment_type()
      if(0<=high_res<=2.0)and("X-RAY DIFFRACTION" in  pdb_pair_data_type):
        tup = (pdb_id,chain.id,low_res,data_type,
              pdb_pair_id,pdb_pair_chain_id,high_res,pdb_pair_data_type,hit.identity)
        if tup is not None:
          print tup,"tup"*40         
  return tup

def run(file_path):
   result = None
   if os.path.isfile(file_path):
     pdb_inp = iotbx.pdb.input(file_path)
     data_type = pdb_inp.get_experiment_type()
     res = pdb_inp.resolution()
     if (data_type == "ELECTRON MICROSCOPY") and (3.0 <= res <=6.0):
       try:
         model = mmtbx.model.manager(model_input=iotbx.pdb.input(file_path))
         tup = get_perfect_pair(res,data_type,model,file_path)
         if tup:
           result=tup
       except Exception, e:
         print "%s: %s"%(file_path, e)        
   return result
 
if __name__ == '__main__':
    files=[]
    work_dir = "/home/wensong/pdb_bob/" 
    pdb_info = easy_pickle.load(file_name='pdb_dict.pickle')
    for key,value in pdb_info.items():
      file_name = "pdb"+key.lower()+".ent.gz"
      file_path = os.path.join(work_dir,key.lower()[1:3],file_name)
      if os.path.exists(file_path):
        files.append(file_path)
 
    pair_list=[]
    for f in files:
      try:
        result = run(f)
        if result:
          pair_list.append(result)
      except Exception,e:
        print e
    if pair_list is None: print "None pairList"*20
    if pair_list is not None:
      for pair in pair_list:
        if pair is None:
          pair_list.remove(pair)
    if pair_list:
      easy_pickle.dump("X_RAY_CYRO_EM_PAIR_98_new.pkl",pair_list)
    else:print "why None?end"*15
  

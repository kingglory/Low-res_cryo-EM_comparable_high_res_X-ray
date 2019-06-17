from __future__ import division
import os
import sys
from libtbx import easy_pickle
import iotbx.pdb
from libtbx import group_args
from libtbx.easy_mp import pool_map 



def get_experimental_pdb_info(file_name):
  pdb_info = None
  prefix = os.path.basename(file_name)[3:7]
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  r_free = pdb_inp.get_r_rfree_sigma().r_free
  r_work = pdb_inp.get_r_rfree_sigma().r_work
  resolution = pdb_inp.resolution()
  data_type = pdb_inp.get_experiment_type()
  if 0< resolution <= 6.0:
    if "X-RAY DIFFRACTION" or "ELECTRON MICROSCOPY" in data_type:
      pdb_info = group_args(
        data_type  = data_type,
        resolution = resolution,
        r_free     = r_free,
        r_work     = r_work)
    crystal_symmetry = pdb_inp.crystal_symmetry()
    if crystal_symmetry.space_group():
      pdb_info.space_group = crystal_symmetry.space_group().type().lookup_symbol()
  return (prefix,pdb_info)

if __name__ == '__main__':
    #if 0:
    rdict = {}
    path = "/home/wensong/pdb_bob/"
    #dpath = "/home/wangwensong/PDB/structure_factors/"    
    of = open("".join([path,"INDEX"]),"r")
    files = ["".join([path,f]).strip() for f in of.readlines()]
    for f in files:
      if not os.path.exists(f):
        files.remove(f)
    of.close()
    #of = open("".join([dpath,"INDEX"]),"r")
    #dfiles = [os.path.basename("".join([path,f]).strip())[1:5] for f in of.readlines()]
    #of.close()
    #pdb_code = os.path.basename(f)[3:7]
    # if(pdb_code in dfiles):
    #   have_experimental_data = True
    tup_list=pool_map(
      processes=60,
      func=get_experimental_pdb_info,
      iterable=files)
    try:
      for tup in tup_list: 
        if tup is not None:
          if tup[0] or tup[1] is not None:
            rdict[tup[0]] = tup[1]
    except Exception,e:
        print "Failed",str(e)
    easy_pickle.dump(file_name='pdb_dict.pickle', obj=rdict)
    #else:
    #info=get_experimental_pdb_info("/home/wangwensong/PDB/pdb/us/pdb1us0.ent.gz")
    # print info
    #get_experimental_pdb_info("/home/pdb/pdb/us/pdb1us0.ent.gz")

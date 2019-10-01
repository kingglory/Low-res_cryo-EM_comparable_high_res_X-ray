from __future__ import division
import os
import mmtbx.model
import iotbx.pdb
from libtbx import group_args
from libtbx import easy_pickle






















if __name__ == '__main__':
  print "start"
  path="/home/wangwensong/PDB/pdb/"
  of = open("".join([path,"INDEX"]),"r")
  files = ["".join([path,f]).strip() for f in of.readlines()]
  of.close()
  Low_res_cryo_EM=[]
  High_res_X_ray=[]
  for f in files:
    pdb_inp = iotbx.pdb.input(file_name=f)
    resolution = pdb_inp.resolution()
    data_type = pdb_inp.get_experiment_type()
    if resolution is None:continue
    if data_type is None:continue

    if resolution <= 2.0:
      if data_type == "X-RAY DIFFRACTION":
        High_res_X_ray.append(f)
        #print f, resolution, data_type, "\n"
    if 3.0<= resolution <= 6.0:
      if data_type == "ELECTRON MICROSCOPY":
        Low_res_cryo_EM.append(f)
        print f, resolution, data_type, "\n"
  easy_pickle("High-res-X-ray.pkl",High_res_X_ray)
  easy_pickle("Lowr-res-cryo-EM.pkl",Low_res_cryo_EM)
  print len(Low_res_cryo_EM),len(High_res_X_ray)



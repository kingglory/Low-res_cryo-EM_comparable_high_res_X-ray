from libtbx.easy_mp import pool_map
from libtbx import easy_pickle


i=0
pdb_pair = easy_pickle.load("X_RAY_CYRO_EM_PAIR_98_new.pkl")
for pair in pdb_pair:
  if pair is not None:
    i=i+1
    print pair
    print i

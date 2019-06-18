from libtbx.easy_mp import pool_map
from libtbx import easy_pickle



pdb_pair = easy_pickle.load("X_RAY_CYRO_EM_PAIR.pkl")
for pair in pdb_pair:
  if pair is not None:
    print pair

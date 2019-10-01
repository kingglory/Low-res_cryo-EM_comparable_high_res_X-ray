from libtbx import easy_pickle


i=0
pdb_pair = easy_pickle.load("X_RAY_CYRO_EM_PAIR_100_whole.pkl")
for pair in pdb_pair:
  if pair is not None:
    i=i+1
    print pair
    f = open("pair_100_new.txt","a")
    f.write(str(pair)+"\n")
    f.close()
    print i

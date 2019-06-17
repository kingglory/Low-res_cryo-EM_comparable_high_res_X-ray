from libtbx import easy_pickle
Low_res_cryo_EM=easy_pickle.load("Low-res-cryo-EM.pkl")
High_res_X_ray=easy_pickle.load("High-res-X-ray.pkl")
for tup in High_res_X_ray:
  f = tup[0]
  #print tup
for tup in Low_res_cryo_EM:
  #print tup
  f = tup[0]
print (len(High_res_X_ray),len(Low_res_cryo_EM))

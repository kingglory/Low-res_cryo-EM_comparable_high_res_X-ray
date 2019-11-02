from libtbx import easy_pickle

percent = easy_pickle.load("ss-percent_3_6.pkl")
print (len(percent.items()))
i=0
for key, value in percent.items():
  i=i+1
  print(i,key, value)

from libtbx import easy_pickle
fo = easy_pickle.load("hydrogen_bond_statistics.pkl")
list_0_9=[]
list_1_2=[]
list_3_5=[]
for i in fo:
  if i is not None:
    if 0<i[5]<0.9:
      list_0_9.append(i)
    if 0<i[5]<=1.2:
      list_1_2.append(i)
    if i[5]>=3.5:
      list_3_5.append(i)
print (len(list_0_9),len(list_1_2),len(list_3_5))
easy_pickle.dump("res_0_9.pkl",list_0_9)
easy_pickle.dump("res_1_2.pkl",list_1_2)
easy_pickle.dump("res_3_5.pkl",list_3_5)

      

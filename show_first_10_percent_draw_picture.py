from libtbx import easy_pickle
from operator import itemgetter
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
whole_num_of_cryo_EM = 4040
def exercise_0():
  f1= easy_pickle.load("ss-percent_3_6_high_res.pkl")
  f_1 = sorted(f1.items(),key=itemgetter(1),reverse=True)
  labels = ["matching","unmatching"]
  X = [len(f_1),whole_num_of_cryo_EM-len(f_1)]
  fig = plt.figure()
  plt.pie(X,labels=labels,autopct="%1.2f%%")
  plt.title("3_6_high_percent")
  plt.show()
  plt.savefig("PieChart_3_6_high.jpg")
  print (f_1)
def exercise_1():
  f2= easy_pickle.load("ss-percent_all_high_res.pkl")
  f_2 = sorted(f2.items(),key=itemgetter(1),reverse=True)
  labels = ["matching","unmatching"]
  X = [len(f_2),whole_num_of_cryo_EM-len(f_2)]
  fig = plt.figure()
  plt.pie(X,labels=labels,autopct="%1.2f%%")
  plt.title("all_high_percent")
  plt.show()
  plt.savefig("PieChart_all_high.jpg")
  print (f_2)
f3= easy_pickle.load("ss-percent_3_6.pkl")
print (len(f3))



f4= easy_pickle.load("ss-percent.pkl")
print (len(f4))

exercise_0()
#exercise_1()
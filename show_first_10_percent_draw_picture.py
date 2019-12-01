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
  list=[]
  for key,value in f1.items():
    list.append(value)
  plt.hist(list,bins=12,color="c",density=False,range=(0,100))
  plt.title("3_6_high_hist")
  plt.show()
  list_0_20=[]
  list_20_40=[]
  list_40_60=[]
  list_60_100=[]
  for i in list:
    if 0<=i<20:
      list_0_20.append(i)
    if 20<=i<40:
      list_20_40.append(i)
    if 40<=i<60:
      list_40_60.append(i)
    if 60<=i<=100:
      list_60_100.append(i)
  labels = ["0~20","20~40","40~60","60~100"]
  Y = [len(list_0_20),len(list_20_40),len(list_40_60),
      len(list_60_100)]
  print (Y)
  fig = plt.figure()
  plt.pie(Y, labels=labels, autopct="%1.2f%%")
  plt.title("3_6_matching_percent")
  plt.show()
def exercise_1():
  f2= easy_pickle.load("ss-percent_all_high_res.pkl")
  f_2 = sorted(f2.items(),key=itemgetter(1),reverse=True)
  labels = ["matching","unmatching"]
  X = [len(f_2),whole_num_of_cryo_EM-len(f_2)]
  fig = plt.figure()
  plt.pie(X,labels=labels,autopct="%1.2f%%")
  plt.title("all_high_percent")
  plt.show()
def exercise_2():
  f3= easy_pickle.load("ss-percent_3_6.pkl")
  print (len(f3))


def exercise_3():
  f4= easy_pickle.load("ss-percent_all.pkl")
  f_4 = sorted(f4.items(), key=itemgetter(1), reverse=True)
  labels = ["matching", "unmatching"]
  X = [len(f_4), whole_num_of_cryo_EM - len(f_4)]
  fig = plt.figure()
  plt.pie(X, labels=labels, autopct="%1.2f%%")
  plt.title("all_pie")
  plt.show()
  list = []
  for key, value in f4.items():
    list.append(value)
  plt.hist(list, bins=12, color="c", density=False, range=(0, 100))
  plt.title("all_hist")
  plt.show()
  list_0_20 = []
  list_20_40 = []
  list_40_60 = []
  list_60_100 = []
  for i in list:
    if 0 <= i < 20:
      list_0_20.append(i)
    if 20 <= i < 40:
      list_20_40.append(i)
    if 40 <= i < 60:
      list_40_60.append(i)
    if 60 <= i <= 100:
      list_60_100.append(i)
  labels = ["0~20", "20~40", "40~60", "60~100"]
  Y = [len(list_0_20), len(list_20_40), len(list_40_60),
       len(list_60_100)]
  print (Y)
  fig = plt.figure()
  plt.pie(Y, labels=labels, autopct="%1.2f%%")
  plt.title("all_matching_percent")
  plt.show()

exercise_3()
#exercise_1()
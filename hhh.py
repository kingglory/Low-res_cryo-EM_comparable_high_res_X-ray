import math
def distant(n):
  if n==1:
    z=0
  else:
    m=int(math.sqrt(n))
    h=int(n-(math.pow(m,2)))
    x=h%(m+1)
    z= abs(x-((m+1)/2))+((m+1)/2)
  print (z)

distant(100000)
distant(2345678)
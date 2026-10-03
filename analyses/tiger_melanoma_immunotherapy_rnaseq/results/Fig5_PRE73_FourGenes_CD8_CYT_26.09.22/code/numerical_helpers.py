import math
import numpy as np

def rank(x):
 a=np.asarray(x);order=np.argsort(a,kind='stable');r=np.empty(len(a),float);i=0
 while i<len(a):
  j=i+1
  while j<len(a) and a[order[j]]==a[order[i]]:j+=1
  r[order[i:j]]=(i+1+j)/2;i=j
 return r

def betacf(a,b,x):
 qab=a+b;qap=a+1;qam=a-1;c=1.;d=1.-qab*x/qap
 if abs(d)<1e-300:d=1e-300
 d=1/d;h=d
 for m in range(1,401):
  m2=2*m;aa=m*(b-m)*x/((qam+m2)*(a+m2));d=1+aa*d
  if abs(d)<1e-300:d=1e-300
  c=1+aa/c
  if abs(c)<1e-300:c=1e-300
  d=1/d;h*=d*c;aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2));d=1+aa*d
  if abs(d)<1e-300:d=1e-300
  c=1+aa/c
  if abs(c)<1e-300:c=1e-300
  d=1/d;delta=d*c;h*=delta
  if abs(delta-1)<3e-14:return h
 raise AssertionError('beta failed convergence')

def ibeta(a,b,x):
 if x==0:return 0.
 if x==1:return 1.
 bt=math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log1p(-x))
 return bt*betacf(a,b,x)/a if x<(a+1)/(a+b+2) else 1-bt*betacf(b,a,1-x)/b

from numpy.linalg import inv,svd,solve,eig
from numpy import real,imag,dot
from numpy import pi
from pandas import read_csv
import matplotlib.pyplot as plt
from spsolver import spsolver
import time

t0=time.time()
Udata=read_csv('../DataSets/GFUdata.csv',header=None,index_col=0)
print(time.time()-t0)
Udata=Udata.values

rk=380;

U0=Udata[:,:rk]
U1=Udata[:,1:(rk+1)]

t0=time.time()
A0=spsolver(U0,U1,rk,1e-4)
print(time.time()-t0)

plt.spy(A0)
plt.show()

t0=time.time()
l0,w0=eig(A0)
print(time.time()-t0)

plt.plot(real(l0),imag(l0),'r.')
plt.show()

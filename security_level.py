from estimator.lwe_parameters import LWEParameters as LWEest
from estimator import *
import math
etas=2
etae=2
sig=2*math.sqrt(etae/2)
Q=130561 #2**17
Q=2**21+2**19+2**18 #=2^26
Q=2**27 #=2^32
Q=2**30 #=2^38
Q=2**34 #=2^38

Q=2**44+2**43
N=256*10
params = LWEest(
     n=N,
     q=Q,#1280,256 bit
     Xs=ND.CenteredBinomial(2*etas), # s
     Xe=ND.CenteredBinomial(etae), # 2e
     m=N,
     tag="example"
)
print(f"Q={Q},N={N}")
r=LWE.estimate(params)

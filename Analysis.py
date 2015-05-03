#!/usr/local/bin/python

import numpy
from hmm import *
from config import *

if __name__ == "__main__":
    # State 0 -> No Aneuploidy
    # State 1 -> Aneuploidy
    V = [str(x) for x in range(1,10)]
    A = numpy.array([[.999,.001],[.999,.001]])
    B = numpy.array([
        [.1, .1, .1, .1, .2, .1, .1, .1, .1],
        [.1, .1, .1, .1, .2, .1, .1, .1, .1]])
    Pi = numpy.array([.99, .01])

    del11m = HMM(2, A=A, B=B, V=V, Pi=Pi)
    del11p = HMM(2, A=A, B=B, V=V, Pi=Pi)
    del13m = HMM(2, A=A, B=B, V=V, Pi=Pi)
    del13p = HMM(2, A=A, B=B, V=V, Pi=Pi)
    completem = HMM(2, A=A, B=B, V=V, Pi=Pi)
    completep = HMM(2, A=A, B=B, V=V, Pi=Pi)
    longdm = HMM(2, A=A, B=B, V=V, Pi=Pi)
    longdp = HMM(2, A=A, B=B, V=V, Pi=Pi)
    dup11m = HMM(2, A=A, B=B, V=V, Pi=Pi)
    dup11p = HMM(2, A=A, B=B, V=V, Pi=Pi)
    standard = HMM(2, A=A, B=B, V=V, Pi=Pi)

    obs_l = []
    for i in range(10):
        with open(OUTPUTPATH + 'newtestmm' + str(i)) as f:
            obs_l.append(f.readline())
    
    numpy.set_printoptions(precision=2)

    print(forward(del11m, obs_l[0]))
    del11m = baum_welch(del11m, obs_l, epochs=2, scaling = False, updatePi = True, verbose=True)
    print(forward(del11m, obs_l[0]))

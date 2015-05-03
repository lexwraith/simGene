#!/usr/local/bin/python

import numpy
from hmm import *
from config import *
from pprint import pprint

def trainHMM(hmm, filepaths):
    obs_l = []
    for elem in filepaths:
        with open(elem) as f:
            obs_l.append(f.readline())

    return baum_welch(hmm, obs_l, epochs = 5, updatePi = True, verbose = True)

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
    for i in range(5):
        with open(OUTPUTPATH + 'newtestmm' + str(i)) as f:
            obs_l.append(f.readline())
    
    numpy.set_printoptions(precision=2)

    pre = map(lambda x: forward(del11m, x)[0], obs_l)
    pprint(pre)
    del11m = baum_welch(del11m, obs_l, epochs=2, updatePi = True, verbose=True)
    post = map(lambda x: forward(del11m, x)[0], obs_l)
    pprint(post)

    pprint(map(lambda x, y: x - y, pre, post))


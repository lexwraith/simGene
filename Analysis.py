#!/usr/local/bin/python

import numpy
import itertools
from hmm import *
from config import *
from pprint import pprint

def genPaths(filebase, gender="m"):
    return [filebase + s + str(n) for s in gender for n in range(5)]

def checkHMM(hmm, filepaths):
    obs_l = []
    for elem in filepaths:
        with open(OUTPUTPATH + elem) as f:
            obs_l.append(f.readline())

    obs_l = [forward(hmm,i)[0] for i in obs_l]
    avg = sum(obs_l)/len(obs_l)
    return avg

def trainHMM(hmm, filepaths):
    obs_l = []
    for elem in filepaths:
        with open(OUTPUTPATH + elem) as f:
            obs_l.append(f.readline())

    return baum_welch(hmm, obs_l, epochs = 2, updatePi = True, verbose = True)

if __name__ == "__main__":
    numpy.set_printoptions(precision=2)

    # State 0 -> No Aneuploidy
    # State 1 -> Aneuploidy
    V = [str(x) for x in range(1,10)]
    A = numpy.array([[.999,.001],[.999,.001]])
    B = numpy.array([
        [.1, .1, .1, .1, .2, .1, .1, .1, .1],
        [.1, .1, .1, .1, .2, .1, .1, .1, .1]])
    Pi = numpy.array([.99, .01])

    del_11_m = genPaths('22q11del')
    del_11_p = genPaths('22q11del', 'p')
    del_13_m = genPaths('22q13del')
    del_13_p = genPaths('22q13del', 'p')
    dup_11_m = genPaths('22q11dup')
    dup_11_p = genPaths('22q11dup', 'p')
    longd_m = genPaths('longd')
    longd_p = genPaths('longd', 'p')
    complete_m = genPaths('complete')
    complete_p = genPaths('complete', 'p')
    none = genPaths('none', 'p')
    

    ALL = list(itertools.chain(del_11_m, del_11_p, del_13_m, del_13_p, dup_11_m, dup_11_p, longd_m, longd_p, complete_m, complete_p, none))

    del11m = HMM(2, A=A, B=B, V=V, Pi=Pi)
    trainHMM(del11m, del_11_m)

    del11p = HMM(2, A=A, B=B, V=V, Pi=Pi)
    trainHMM(del11p, del_11_p)

    del13m = HMM(2, A=A, B=B, V=V, Pi=Pi)
    trainHMM(del13m, del_13_m)

    del13p = HMM(2, A=A, B=B, V=V, Pi=Pi)
    trainHMM(del13p, del_13_p)
    
    """
    completem = HMM(2, A=A, B=B, V=V, Pi=Pi)
    trainHMM(completem, complete_m)

    completep = HMM(2, A=A, B=B, V=V, Pi=Pi)
    trainHMM(completep, complete_p)

    longdm = HMM(2, A=A, B=B, V=V, Pi=Pi)
    trainHMM(longdm, longd_m)

    longdp = HMM(2, A=A, B=B, V=V, Pi=Pi)
    trainHMM(longdp, longd_p)
    """
    
    dup11m = HMM(2, A=A, B=B, V=V, Pi=Pi)
    trainHMM(dup11m, dup_11_m)

    dup11p = HMM(2, A=A, B=B, V=V, Pi=Pi)
    trainHMM(dup11p, dup_11_p)

    standard = HMM(2, A=A, B=B, V=V, Pi=Pi)
    trainHMM(standard, none)

    data = [del_11_m, del_11_p, del_13_m, del_13_p, dup_11_m, dup_11_p] 
    names = ['DEL-11-m', 'DEL-11-p', 'DEL-13-m', 'DEL-13-p', 'DUP-11-m', 'DUP-11-p']
    hmms = [del11m, del11p, del13m, del13p, dup11m, dup11p]
    colnames = ",".join(['DEL-11-m', 'DEL-11-p', 'DEL-13-m', 'DEL-13-p', 'DUP-11-m', 'DUP-11-p'])
    colnames = "," + colnames + "\n"
    print(colnames)
    with open("heatmap.csv", 'w+') as f:    
        f.write(colnames)
        for i in range(len(hmms)):
            print(i)
            f.write(names[i] + ",")
            for j in range(len(data)):
                g = checkHMM(hmms[i], data[j])
                f.write(str(g) + ",")
            f.write("\n")


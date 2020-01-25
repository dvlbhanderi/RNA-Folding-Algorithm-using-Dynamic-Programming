

import numpy
import json
import ast

 
DEBUG_FILL = False          
DEBUG_TRACEBACK = False

class Structure :

    def __init__(self):
        self.pair = []
        self.unpair = []

    def merge(self, s):
        for i in s.pair :
            self.pair.append(i)
        for i in s.unpair :
            self.unpair.append(i)

    def __str__(self):
        return "pair\n------\n" +\
               str(self.pair) + "\n" +\
               "unpair\n-------\n" +\
               str(self.unpair)

def alphabetLookup(c) :

    if (c=='A') : return 0
    if (c=='C') : return 1
    if (c=='U') : return 2
    if (c=='G') : return 3
    print "error - don't know this char: " + c


def cost(i,j,F,sgm,seq,comp):

    
    h = F[i,j-1]
    v = F[i+1,j]
    d = F[i+1,j-1] + sgm[alphabetLookup(seq[i]), alphabetLookup(seq[j])]
    
    return comp(h,v,d)

def fill(F, sgm, seq, comp):
    
    n = len(seq)

    for k in range(1,n):
        for i in range(0,n-k):
            j = k+i

            c = cost(i,j,F, sgm, seq, comp)
            bifurcate = []
            for h in range(i,j):
                bifurcate.append(F[i,h] + F[h+1,j])
                
            if (DEBUG_FILL) :
                print "doing " + str(i) + "," + str(j)
                print "cost is " + str(c)
                bfr = comp(bifurcate)
                msg = "bifurcate is: " + str(bfr)
                if comp(comp(bifurcate),c) == bfr and bfr != c : msg += "*\n"
                else : msg += "\n"
                print msg

            F[i,j] = comp(comp(bifurcate), c)


sgm = numpy.zeros([4,4], float)
sgm.fill(0)
sgm[alphabetLookup('A'), alphabetLookup('U')] = 1
sgm[alphabetLookup('U'), alphabetLookup('A')] = 1
sgm[alphabetLookup('G'), alphabetLookup('C')] = 1
sgm[alphabetLookup('C'), alphabetLookup('G')] = 1


def traceback(F, sgm, seq, comp, initi, initj) :

    i = initi
    j = initj

    struct = Structure()
    
    while (j != i-1):
        if F[i,j] == F[i,j-1] :
            # unpair
            if DEBUG_TRACEBACK : print "the " + str(j) + "th char is not paired"
            struct.unpair.append(j)
            j -= 1
        elif F[i,j] == F[i+1,j] :
            # unpair
            if DEBUG_TRACEBACK : print "the " + str(i) + "th char is not paired"
            struct.unpair.append(i)
            i += 1
        elif F[i,j] == F[i+1,j-1] + sgm[alphabetLookup(seq[i]), alphabetLookup(seq[j])] :
            # pair
            if DEBUG_TRACEBACK : print "the " + str(i) + "th char and the " + str(j) + "th chars are paired"
            struct.pair.append((i,j))
            i += 1
            j -= 1
        else :
            best = None
            for h in range(i,j):
                cr = F[i,h] + F[h+1,j]
                if best == None or comp(cr,best) == cr :
                    best = cr
                    bi = i
                    bj = j
                    bh = h
                    if DEBUG_TRACEBACK : msg = "combined F[" + str(i) + "," + str(h) + "] and F[" + str(h+1) + "," + str(j) + "]"
            if DEBUG_TRACEBACK : print "there was a bifurcation: " + msg
            s1 = traceback(F, sgm, seq, comp, bi, bh)
            s2 = traceback(F, sgm, seq, comp, bh+1, bj)
            struct.merge(s1)
            struct.merge(s2)
            break

    return struct


#seq = "GAAGUUGGGCCG"
#seq = "GGGAAAUCC"
#seq = "AUGUCAGCGUU"
with open("textfile.txt", "r") as myfile:
    seq = myfile.readline().strip()   
n = len(seq)



F = numpy.zeros([n,n],float)

fill(F, sgm, seq, max)

s = traceback(F, sgm, seq, max, 0, n-1)

res = str(s.pair)
result = ast.literal_eval(res)

print "max count of pairs "
print len(result)

def print_struct(seq, result):
    brace = ["." for _ in range(len(seq))]
    for tup in result:
        brace[min(tup)] = "{"
        brace[max(tup)] = "}"
    return "".join(brace)

d = print_struct(seq, result)
print seq
print d

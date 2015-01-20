#!/usr/bin/python
import sys
from random import randint

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)



def BinaryTree(n, dist, iter):
    nodes = [[] for i in xrange(n)]
    for i in range(n):
        nodes[i].append(i)
    while len(nodes) > 1:
        n0 = randint(0, len(nodes)-1)
        n1 = randint(0, len(nodes)-2)
        if n1 >= n0:
            n1+=1
 #       print "n0=" + str(n0) + " n1=" + str(n1)
 #       dist[len(nodes[n0])+len(nodes[n1])-1] += float(len(nodes[n0])*len(nodes[n1]))/(iter*n)
        dist[len(nodes[n0])+len(nodes[n1])-1] += len(nodes[n0])*len(nodes[n1])
        nodes[n0] = nodes[n0] + nodes[n1]
        del nodes[n1]
        print nodes
    return


def BinaryTree2(n, dist, iter):
    nodes = [[] for i in xrange(n)]
    nodes[0] = [1,1]
    nodes[1] = [1,1]
    for i in range(2, n):
        nodes[i] = [0,1]
    while len(nodes) > 1:
        n0 = randint(0, len(nodes)-1)
        n1 = randint(0, len(nodes)-2)
        if n1 >= n0:
            n1+=1
 #       print str(n0) + " " + str(n1)
  #      print nodes
        if nodes[n0][0] == 1 and nodes[n1][0] == 1:
            dist[nodes[n0][1]+nodes[n1][1]-1] += 1
            return
        else:
            nodes[n0][0] = nodes[n0][0] + nodes[n1][0]
            nodes[n0][1] = nodes[n0][1] + nodes[n1][1]
        del nodes[n1]
    return




print "Binary trees simulator"
#n=10
try:
    n = num(sys.argv[1])
except IndexError:
    print "Please set the tree size"
    sys.exit()

try:
    iter = num(sys.argv[2])
except IndexError:
    print "The number of iterations was set to " + str(n)
    iter = n

dist = [0.0]*n

for i in range(0,iter):
    BinaryTree2(n, dist, iter)

print "Distances"
sum = 0.0
for i in range(1,n):
    print float(dist[i])/iter
    sum += float(dist[i])/iter
print sum
#    print str (i+1) + " " + str(int(dist[i]))

#int dist
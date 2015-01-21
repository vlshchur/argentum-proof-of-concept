#!/usr/bin/python
import sys
import re
import os


def RWrite3DListInFile(arr, dir, fname):
    n = len(arr)
    for i in range(n):
        if i % 200 == 0:
            dn = dir + str(i/200)
            try:
                os.makedirs(dn)
            except:
                print "Impossible to create directory " + dn
                sys.exit()
        fn = dn + "/" + fname + "_" + str(i).zfill(3) + ".csv"
        f = open(fn, 'w')
        colNames = "\"Size\",\"Frequency\",\"Haplotype1\",\"Haplotype2\"\n"
        for j in range(n):
            if i == j:
                arr[i][j][0] = 1
            else:
                arr[i][j][0] = sum(arr[i][j])
                if arr[i][j][0] == 0:
                    arr[i][j][0] = 1
        for j in range(n):
            for k in range(1,n):
                fr = float(arr[i][j][k])/arr[i][j][0]
                datStr = str(k+1) + "," + str(fr) + "," +str(i) + "," + str(j) + "\n"
                f.write(datStr)
    f.close()




def InitializeLists(arr, n):
    for i in range(n):
        den1 = []
        for j in range(n):
            den = [0]*n
            den1.append(den)
        arr.append(den1)


try:
    fname = sys.argv[1]
except IndexError:
    print "Please specify the filename"
    sys.exit()


n = -1
counter = 0
density = []

with open(fname) as f:
    for line in f:
        line = line.rstrip()
        s = list(line)
        if len(s) == 0:
            continue
        if n == -1:
            n = len(s)
            print str(n) + " haplotypes."
            InitializeLists(density, n)
        elif len(s) != n & len(s) > 0:
            break
        fr = s.count('1')
        if fr + s.count('0') != n:
            break
        counter+=1
        for i in range(n):
            if s[i] == '1':
                for j in range(i+1, n):
                    if s[j] == '1':
                        density[i][j][fr-1] += 1
                        density[j][i][fr-1] += 1

RWrite3DListInFile(density, "frequency", "freq")


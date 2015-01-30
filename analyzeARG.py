#!/usr/bin/python
import sys
import re
import os

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)



def MinBranchesAverage(string, density, n):
    s = re.compile('[0-9]+')
    string = s.sub('1', string)
    for i in range(n-1):
        s = re.compile('\(([0-9]+),([0-9]+)\)')
        m = re.search(s, string)
        density[int(m.group(1))+int(m.group(2)) - 1] += int(m.group(1))*int(m.group(2))
        su = '\('+m.group(1)+','+m.group(2)+'\)'
        s = re.compile(su)
        string = s.sub(str(int(m.group(1))+int(m.group(2))), string, 1)






def MinBranches(string, density, n):
    string = re.sub(r"(\d+)", r"[\1]", string)
    for i in range(n-1):
        s = re.compile("\(\[([0-9,]+)\],\[([0-9,]+)\]\)")
        m = re.search(s, string)
        g1 = m.group(1).split(",")
        g2 = m.group(2).split(",")
        ms = len(g1) + len(g2) -1
        for h1 in g1:
            for h2 in g2:
                density[int(h1)][int(h2)][ms] += 1
                density[int(h2)][int(h1)][ms] += 1
        r = "["+m.group(1)+','+m.group(2)+']'
        string = re.sub(s, r, string, 1)




def WriteListInFile(arr, fname, header = ""):
    if header == "":
        header = fname
    f = open(fname, 'w')
    f.write(header+"\n")
    for el in arr:
        f.write(str(float(el)/ll)+'\n')
    f.close()


def Write3DListInFile(arr, dir, fname, mode):
    n = len(arr)
    if mode != "csv":
        mode = "tab"
        print "You chose tabular format for the output."
    else:
        print "You chose csv format for the output."
    for i in range(n):
        if i % 200 == 0:
            dn = dir + str(i/200)
            try:
                os.makedirs(dn)
            except:
                print "Impossible to create directory " + dn
                sys.exit()
        fn = dn + "/" + fname + "_" + str(i).zfill(3)
        if mode == "tab":
            fn = fn + ".txt"
        else:
            fn = fn + ".csv"
        f = open(fn, 'w')
        if mode == "csv":
            colName = "\"Size\",\"Frequency\",\"Haplotype1\",\"Haplotype2\"\n"
        else:
            colName = "Size\tFrequency\tHaplotype1\tHaplotype2\"\n"
        f.write(colName)
        for j in range(n):
            if i == j:
                arr[i][j][0] = 1
            else:
                arr[i][j][0] = sum(arr[i][j])
                if arr[i][j][0] == 0:
                    arr[i][j][0] = 1
        for j in range(n):
            for k in range(1,n):
                f.write(str(k+1))
                if mode == "csv":
                    f.write(",")
                else:
                    f.write("\t")
                f.write(str(float(arr[i][j][k])/arr[i][j][0]))
                if mode == "csv":
                    f.write(",")
                else:
                    f.write("\t")
                f.write(str(i))
                if mode == "csv":
                    f.write(",")
                else:
                    f.write("\t")
                f.write(str(j))
                f.write("\n")
        f.close()





def InitializeLists(mode, arr, n):
    if mode == "-a":
        for i in range(n):
            arr.append(0)
    if mode == "-i":
        for i in range(n):
            den1 = []
            for j in range(n):
                den = [0]*n
                den1.append(den)
            arr.append(den1)



print "Clustering method: read and analyze trees in Newick format"
#n=10
try:
    mode = sys.argv[1]
except IndexError:
    print "Set mode: -a for the average over the population, -i for one-to-one, -h for help"
    sys.exit()

if mode != "-a" and mode != "-i" and mode != "-h":
    print "Set mode: -a for the average over the population, -i for one-to-one, -h for help"
    sys.exit()

if mode == "-h":
    print "First parameter is the mode. Use -a for the average over the population, -i for one-to-one, -h for help."
    print "Second parameter specifies the filename."
    print "-ll [integer]: set the line limit (the number of trees to be processed)."
    print "-sl [integer]: set the sampling parameter."
    print "-of [tab, csv]: set the output format (\'csv\' by default)."
    print "-lc [integer]: set the line counter."
    sys.exit()

try:
    fname = sys.argv[2]
except IndexError:
    print "Please specify the filename"
    sys.exit()

ll = -1
sl = 1
of = "csv"
lc = 500
for i in range(3, len(sys.argv)):
    if sys.argv[i] == "-ll":
        i += 1
        try:
            ll = num(sys.argv[i])
        except IndexError:
            print "Cannot set line limit, the whole data will be processed."
    elif sys.argv[i] == "-sl":
        i += 1
        try:
            sl = num(sys.argv[i])
        except IndexError:
            print "Cannot set sampling parameter, the whole data will be processed."
    elif sys.argv[i] == "-lc":
        i += 1
        try:
            lc = num(sys.argv[i])
        except IndexError:
            print "Cannot set sampling parameter, the the default value \'500\' will be used."
    elif sys.argv[i] == "-of":
        i += 1
        try:
            of = sys.argv[i]
        except IndexError:
            print "Cannot set output format, \'csv\' format will be used."
        if of != "csv" and of != "tab":
            print "Cannot set output format, \'csv\' format will be used."

f = open(fname, 'r')

n = -1
counter = 0
density = []
densityAverage = []

for line in f:
    if line.startswith("NEWICK_TREE:"):
        counter += 1
        if (counter-1)%sl != 0:
            continue
        s = re.compile('NEWICK_TREE:\t\[[0-9]+\]|\:[0-9]+\.[0-9e\-]+')
        line = s.sub('', line)
        if n == -1:
            n = line.count("(") + 1
            print str(n) + " leaves."
            if mode == "-a":
                InitializeLists(mode, densityAverage, n)
            if mode == "-i":
                InitializeLists(mode, density, n)
        else:
            if n-1 != line.count("("):
                print "Wrong number of leaves. Aborted."
                sys.exit()
#        line = '(((0,1),(2,3)),(4,5))'
#        print line[0:100]
        if mode == "-a":
            MinBranchesAverage(line, densityAverage, n)
        if mode == "-i":
            MinBranches(line, density, n)
        if counter%lc == 0:
            print str(counter) + ' lines analyzed'
        if counter == ll*sl:
            break
f.close()
if mode == "-a":
    WriteListInFile(densityAverage, "densityAverage.txt", header = "MACS_correct")
if mode == "-i":
    Write3DListInFile(density, "data", "densityOneToOne", of)

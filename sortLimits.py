from sys import argv
from operator import itemgetter
import os
#------------------------------------------------------
def readdata(document):
    file = open(document,"r")
    lines = file.readlines()
    file.close()
    lists = []
    for line in lines:
        p = line.split()
        L = []
        for i in range(100):
            try:
                L.append(p[i])
            except IndexError:
                break
        lists.append(L)
    return lists
#------------------------------------------------------
data=readdata(str(argv[1]))
limits=[]
for i in range(len(data)):
    if (i%2 != 0):
        cut_combo=[]
        cut_combo.append(data[i-1])
        cut_combo.append(data[i])
        span=float(data[i][-1]) - float(data[i][-2])
        cut_combo.append(span)
        limits.append(cut_combo)
limits.sort(key=itemgetter(-1))
for i in range(len(limits)):
    print limits[i]
path=os.path.dirname(str(argv[1]))
print path
f=open(str(path)+"/sortedLimits", "w")
for x in limits:
    f.write(str(x[0][0]) + "\n")
    entry=""
    for y in x[1]:
        entry += str(y)+" "
    f.write(entry + "\n")
    f.write("Span = "+str(x[2]) + "\n")


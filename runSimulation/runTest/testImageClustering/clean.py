#!/usr/bin/python2.7
import random

fpos = open("paper720_30cam_pos_origin.txt","r");
fposOut = open("paper720_30cam_pos.txt", "w");

i = 0;
posList = [ ];
tmpList = []
for line in fpos:
    if i == 0:
        fposOut.write(line)
    else:
        x = float(line.split()[0])
        y = float(line.split()[1])
        tmpList.append(x)
        tmpList.append(y)
        posList.append(tmpList);
        
    i  = i + 1
    tmpList = []

for j in range(len(posList)):
    for k in range(len(posList)):
        if j != k and posList[k] == posList[j]:
            posList[k][0] = posList[k][0] + random.uniform(0,0.01)
for j in range(len(posList)):
    tmpLine = str(posList[j][0])+"\t"+str(posList[j][1])+"\n"
    fposOut.write(tmpLine)

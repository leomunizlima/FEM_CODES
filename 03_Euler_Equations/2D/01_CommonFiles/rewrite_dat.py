#importing a library
import sys
import math

file_name = sys.argv[1]
fp = open(file_name,"r")
lines = fp.readlines()
fp.close()
i = 0
nnodes = int(lines[0].split('\n')[0])
nelem = int(lines[nnodes+1].split('\n')[0])
x = []
y = []
b1 = []
b2 = []
b3 = []
b4 = []
e1 = []
e2 = []
e3 = []
remove = [2,5,6,7,8,9,10,11,12]

i = 0
for l in lines:   
    if i>0 and i<=nnodes and i-1 not in remove:
        x.append(l.split('\t')[0])
        y.append(l.split('\t')[1])
        b1.append(l.split('\t')[2])         
        b2.append(l.split('\t')[3])         
        b3.append(l.split('\t')[4])         
        b4.append(l.split('\t')[5])
    i=i+1        

i = 0

for l in lines:
    if i>nnodes+1: 
        v1 = int(l.split('\t')[0])
        v2 = int(l.split('\t')[1])
        v3 = int(l.split('\t')[2])
        e1.append(v1)
        e2.append(v2)
        e3.append(v3)
    i = i + 1

for r in range(len(remove)):
    for i in range(len(e1)):
        if e1[i]>=remove[r]:
            e1[i] = e1[i]-1
        if e2[i]>=remove[r]:
            e2[i] = e2[i]-1
        if e3[i]>=remove[r]:
            e3[i] = e3[i]-1
    

file_name = file_name.split('/')[-1]
fp = open("new"+file_name,"w")
print(nnodes-len(remove),file=fp)
for i in range(nnodes-len(remove)):
    print(x[i],'\t',y[i],'\t',b1[i],'\t',b2[i],'\t',b3[i],'\t',b4[i],file=fp)

print(nelem,file=fp)
for i in range(nelem):
    print(e1[i],'\t',e2[i],'\t',e3[i],'\t',"-1",file=fp)

fp.close()
    








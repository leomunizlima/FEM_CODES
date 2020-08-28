#importing a library
import sys
import math

file_name = sys.argv[1]
fp = open(file_name,"r")
lines = fp.readlines()
fp.close()
nnodes = int(lines[0])
Elements = lines[nnodes+2:]
nelem = int(lines[nnodes+1])
LostNodes = []
TOTAL_LINES = nnodes + nelem + 2
for j in range(nnodes):
    achei = 0
    for e in Elements:
        if int(e.split('\t')[0])==j or int(e.split('\t')[1])==j or int(e.split('\t')[2])==j:
            achei = 1
            break
    if (achei==0):
        print ("Nao achei n=",j)
        LostNodes.append(j)
for i in LostNodes:
    lines[i+1] =lines[i+1].split('\t')[0]+'\t'+lines[i+1].split('\t')[1]+'\t'+"0\t0\t0\t0\n"

print(LostNodes)

file_name_new = "new"+file_name.split('/')[-1]
fp = open(file_name_new,"w")
for i in range(TOTAL_LINES):
    print(lines[i][:-1],file=fp)
fp.close()


    




    



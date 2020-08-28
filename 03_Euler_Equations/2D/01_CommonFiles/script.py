#importing a library
import sys
import math

file_name = sys.argv[1]
fp = open(file_name,"r")
lines = fp.readlines()
fp.close()
t = []
r1 = []
r2 = []
r3 = []
r4 = []

r1a = []
r2a = []
r3a = []
r4a = []

for x in lines:
    t.append(float(x.split('\t')[0]))
    r1.append(-1*math.log(10,float(x.split('\t')[1])))
    r2.append(-1*math.log(10,float(x.split('\t')[2])))
    r3.append(-1*math.log(10,float(x.split('\t')[3])))
    r4.append(-1*math.log(10,float(x.split('\t')[4])))
    r1a.append(float(x.split('\t')[1]))
    r2a.append(float(x.split('\t')[2]))
    r3a.append(float(x.split('\t')[3]))
    r4a.append(float(x.split('\t')[4]))

file_name_out="log"+file_name
f = open(file_name_out, 'w')
for i in range(len(t)):
    print(t[i],"\t",r1[i],"\t",r2[i],"\t",r3[i],"\t",r4[i],file=f)
f.close()

for i in range(len(t)):
    print(t[i],"\t",r1a[i],"\t",r2a[i],"\t",r3a[i],"\t",r4a[i])




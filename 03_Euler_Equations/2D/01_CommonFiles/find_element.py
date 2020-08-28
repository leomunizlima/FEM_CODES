#importing a library
import sys
import math

file_name = sys.argv[1]
n = int(sys.argv[2])
print(n)
fp = open(file_name,"r")
lines = fp.readlines()
fp.close()
i = 0
achei = 0

for x in lines:
    if int(x.split('\t')[0])==n or int(x.split('\t')[1])==n or int(x.split('\t')[2])==n:
        print("Achei em", i)
        achei = 1
    i = i + 1

if (achei==0):
    print ("Nao achei!")




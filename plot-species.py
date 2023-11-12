import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
np.set_printoptions(threshold=sys.maxsize)

#file = open('tmp.out', 'r')
file = open('species.out', 'r')
filelines = file.readlines()

species = []
uniquespecies = []
count1 = 0
for line in filelines:
	count1 += 1
	print('Line Number =', count1)
	splitted = line.split()
	#print(splitted[0])
	#print(line)
	count2 = 0
	for column in splitted:
		#print('Column',count2,splitted[count2])
		if count2 > 3 and splitted[0] == '#':
			#species.append(splitted[count2])
			if splitted[count2] not in uniquespecies:
				uniquespecies.append(splitted[count2])
		count2 += 1
print(uniquespecies)
print(len(uniquespecies))

file.seek(0)
filelines = file.readlines()

count1 = 0
molecule = []
numberofmolecules = np.zeros([int(len(filelines)/2),len(uniquespecies)])
molecule = np.zeros(1000,int)

for line in filelines:
	print('Line Number DataFrame', count1)
	#print(line)
	splitted = line.split()
	count2 = 0
	for column in splitted:
		if count2 > 3 and splitted[0] == '#':
			count3 = 0
			for i in uniquespecies:
				#print(i,splitted[count2])
				if splitted[count2] == i:
					molecule[count2]=count3
					#print('OK',molecule[count2],count2,count3)
				count3 += 1
		if count2 > 2 and splitted[0] != '#':
			numberofmolecules[int(count1/2),molecule[count2+1]]=splitted[count2]
			#print(count1,count2+1,molecule[count2+1],splitted[count2])
		count2 += 1
	count1 += 1
#print(numberofmolecules)

avemol = np.mean(numberofmolecules, axis=0)
maxmol = np.amax(numberofmolecules, axis=0)
print(uniquespecies)
print(avemol)
print(maxmol)
print(len(numberofmolecules[:,0]))
print(int(len(filelines)/2))
print(numberofmolecules.shape)

j = 101
k = 0
print(int(len(filelines)/200))
submolnum = np.zeros([(int(len(filelines)/200))+1,len(uniquespecies)])
for i in range(len(numberofmolecules[:,0])):
	if j == 101:
		submolnum[k,:] = numberofmolecules[i,:]
		print(k,j,i)
		j = 1
		k += 1
	j += 1
print(len(submolnum[:,0]))

fontP = FontProperties()
fontP.set_size('xx-small')
#linestyles = ['k-','k--','k-.','k:','b-','b--','b-.','b:','r-','r--','r-.','r:','c-','c--','c-.','c:','g-','g--','g-.','g:','m-','m--','m-.','m:','y-','y--','y-.','y:','gray-','gray--','gray-.','gray:','lime-','lime--','lime-.','lime:','violet-','violet--','violet-.','violet:','pink-','pink--','pink-.','pink:',]
#linestyles = ['k-','k:','b-','b:','r-','r:','c-','c:','g-','g:','m-','m:','y-','y:'] 
#linedash = ['solid','dashed','dashdot','dotted','solid','dashed','dashdot','dotted','solid','dashed','dashdot','dotted','solid','dashed','dashdot','dotted','solid','dashed','dashdot','dotted','solid','dashed','dashdot','dotted','solid','dashed','dashdot','dotted','solid','dashed','dashdot','dotted','solid','dashed','dashdot','dotted','solid','dashed','dashdot','dotted']
#linecolor = ['black','black','black','black','blue','blue','blue','blue','green','green','green','green','red','red','red','red','cyan','cyan','cyan','cyan','magenta','magenta','magenta','magenta','gray','gray','gray','gray','lime','lime','lime','lime','violet','violet','violet','violet','pink','pink','pink','pink'] 
#marker = ['.','v','<','>','1','2','3','4','s','p','P','*','+','x','.','v','<','>','1','2','3','4','s','p','P','*','+','x','.','v','<','>','1','2','3','4','s','p','P','*','+','x']
#markercolor = ['black','black','black','black','black','black','black','black','black','black','black','black','black','black','blue','blue','blue','blue','blue','blue','blue','blue','blue','blue','blue','blue','blue','blue','red','red','red','red','red','red','red','red','red','red','red','red','red','red']

x = np.linspace(0,int(len(filelines)/2),int(len(filelines)/2))
l = -1
m=0
np.savetxt('species.dat',submolnum)
#dataout = open('species-header.dat',"w")
#dataout.write(uniquespecies)
#dataout.close()

import csv 
with open('species.csv', 'w') as f: 
    write = csv.writer(f) 
    write.writerow(uniquespecies)
    write.writerows(submolnum) 

#for molecule in uniquespecies:
#	l += 1
#	if maxmol[l] > 1: #avemol[l] > 1.5E-3:
#		#plt.plot(x[::10],numberofmolecules[:,l],linestyle=linedash[m],color=linecolor[m],label=molecule)
#		plt.plot(x[::100],submolnum[:,l],linestyle=linedash[m],color=linecolor[m],label=molecule)
#		#plt.scatter(x,numberofmolecules[:,l],color=markercolor[m],label=molecule,marker=marker[m],markevery=10)
#		m += 1
#plt.legend(bbox_to_anchor=(1.10, 0.5),loc="center right", prop=fontP)
#plt.xlabel("Step")
#plt.ylabel("# of molecules")
##plt.xlim([0,256])
##plt.ylim([0,1.2])
#plt.savefig('species.png', dpi=1200)
#


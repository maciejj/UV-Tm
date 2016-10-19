#!/usr/bin/env python

########################################
#Adjustable parameters:
#Temperature in C for energy analysis:
Tenergy = 20

#Reaction type
#monomolecular=1; bimolecular(selfcomplementary)=2; bimolecular(non-self-complementary)=3
reactionType=1

#concentration in M; needed only when reactionType is not equal to 1
concentration=0.000002

#Temperature range in C for analyze. Zero values mean, it analyses whole curves
Tlow=0
Tmax=0

#Range for the lnK fitting. Usually 0.2 - 0.8 is OK 
#If not try for example 0.3 to 0.7
lnKlow=0.2
lnKhigh=0.8
#######################################

from scipy import constants as const
import numpy as np
import os
from string import *
import sys #pozwoli na uzywanie parametrow z lini polecen jako nazw plikow in i out
#import sip
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

Tenergy = Tenergy+273.15
def line_fitting(Temper, Abs):
	T1 = np.vstack([Temper, np.ones(len(Temper))]).T
	a, b = np.linalg.lstsq(T1, Abs)[0]
	return a,b

def hills_fit(T, b):
	return ((T**b)/(melting_temp**b + T**b))


def enthalpy_fit(x, H, S):
	H = H*const.calorie*1000
	S = S*const.calorie*1000
	return -H/(x*const.R) + S/const.R
	
def run_for_one_set(Te,Ab,name):
	zadowolony = "n"
	while (zadowolony!="y"):
		name = name.replace(' ', '-')
		print "Working on "+name
		plt.clf()
		T = np.array(Te)
		#T = T #+ 273.15
		A = np.array(Ab)
		#constraing T for analysis
		if Tlow or Tmax:
			A = A[(T<Tmax)&(T>Tlow)]
			T = T[(T<Tmax)&(T>Tlow)]
		plt.grid(True)
		plt.plot(T,A)
		# Save the figure in a separate file
		# plt.show()
		plt.savefig('AB.png')
		print ("Give 4 points separated by space to fit lines ")
		print ("T starts in " + str(T[0]) + ", ends in "+ str(T[-1]))
		os.system("display AB.png")
		s = raw_input("")
		punkty = map(float, s.split())
		#szukanie indeksow dla zadanych puktow
		indexy = []
		for i in range(len(punkty)):
			indexy.append(np.nonzero(T>=punkty[i])[0][0])
		Tparts=np.split(T,indexy)
		Aparts=np.split(A,indexy)
		if punkty[0] != punkty[1]:
			dolna_prosta = line_fitting(Tparts[1],Aparts[1])
		else:
			dolna_prosta = [0,A[indexy[0]]]
		
		gorna_prosta = line_fitting(Tparts[3],Aparts[3])
		dol = T*dolna_prosta[0] + dolna_prosta[1]
		gora = T*gorna_prosta[0] + gorna_prosta[1]
		f = (gora - A) / (gora - dol)
		Tm = T[np.nonzero(f<=0.5)[0][0]]
		#print Tm	
		#creating global melting_temp to access it in the Hills fit function
		global melting_temp
		melting_temp= Tm
		f2 = 1-f
		HillsParams = curve_fit(hills_fit, T, f2)
		#print HillsParams 
		# preparing f in range 0.2 > f < 0.8
		fEner = f[(lnKlow < f)&(f < lnKhigh)]
		TEner = T[(lnKlow < f)&(f < lnKhigh)]
		TEner = TEner+273.15
		if reactionType == 1:
			KEner = fEner/(1-fEner)
		elif reactionType == 2: 
			KEner = fEner/2*((1-fEner)**2)*concentration
		elif reactionType == 3:
			KEner = 2*fEner/((1-fEner)**2)*concentration
		lnK = np.log(KEner)
		#print lnK
		energy_params = curve_fit(enthalpy_fit, TEner, lnK)
		#print energy_params
		#print energy_params[0][0], energy_params[0][1]*Tenergy, energy_params[0][0]-energy_params[0][1]*Tenergy
		
		plt.clf()
		plt.plot(T, A, 'o', markersize=2)
		plt.plot(T, dol, 'r')
		plt.plot(T, gora, 'r')
		#Tu wstawic nazwe pliku z f
		plt.savefig('WB.png')
		np.savetxt(sys.argv[1]+'_'+name+'_f.out', np.c_[T,f], fmt='%10.5f')
		np.savetxt(sys.argv[1]+'_'+name+'_fit.out', np.c_[T,dol,gora,A], fmt='%10.5f')
		print("Check fitting")		
		os.system("display WB.png")	
		plt.clf()
		plt.plot(1/TEner, lnK, 'o', markersize=2)
		lnKlin=(1/TEner)*energy_params[0][0]*-1000*const.calorie/const.R + energy_params[0][1]*1000*const.calorie/const.R
		np.savetxt(sys.argv[1]+'_'+name+'_lnK_fit.out', np.c_[(1/TEner),lnK,lnKlin], fmt='%10.8f')
		plt.plot(1/TEner, lnKlin)
		plt.savefig('lnK.png')
		print("Is fiting OK? y/n")
		os.system("display lnK.png")	
	
		zadowolony = raw_input("")
	return Tm, name, HillsParams, energy_params
	
infile = open(sys.argv[1], "r")
per_row = []
for line in infile:
	#print line
	line = line.replace(',', '.')
	per_row.append(line.split('\t'))
per_column = zip(*per_row)
infile.close()

#Find a number of 'Cycle row'
cycle_line=per_column[0].index('Cycle')
data_line=cycle_line+7

#Find wavelengths
indexes = [i for i,x in enumerate(per_row[cycle_line]) if x == 'Cycle']
serie=[]
for i in indexes:
	serie.append(per_row[cycle_line][i+1].split())
dl_fali = list(set(zip(*serie)[1]))

#Wybor ktore krzywe analizowac
print 'Found ' + str(len(serie)) + ' series of measurements containing ' + str(len(dl_fali)) + ' wavlengths' + str(dl_fali)

if len(dl_fali) == 1:
	choice = raw_input('Press a number to choose what to analyze: \n ALL:    1 \n heating only: 2 \n')
elif len(dl_fali) == 2:
	choice = raw_input('Press a number to choose what to analyze: \n ALL:    1 \n '+str(dl_fali[0])+' heating: 4 \n'+str(dl_fali)+' heating: 2 \n')

#choice = 1

#zmiana wartosci na floaty i pozbycie sie pustych elementow
per_column2 = []
for i in range(len(serie*4)):
	per_column[i] = per_column[i][data_line:]
	for j in range(len(per_column[i])):
		per_column2.append([])
		try:
			per_column2[i].append(float(per_column[i][j]))
			#print per_column[i][j]
		except ValueError:
			pass
per_column = per_column2

#run
step = 4*int(choice)
Results = []
wartTm = []
wartHill = []
wartH = []
wartTdS = []
wartG = []
for d in range(1, len(per_row[cycle_line]), step):
#for d in range(1, 17, step):
	Results.append(run_for_one_set(per_column[d], per_column[d+1], per_row[cycle_line][d]))
	wartTm.append(Results[-1][0])
	wartHill.append(Results[-1][2][0][0])
	wartH.append(Results[-1][3][0][0])
	wartTdS.append(Results[-1][3][0][1]*Tenergy)
	wartG.append(Results[-1][3][0][0] - Results[-1][3][0][1]*Tenergy)

	with open(sys.argv[1].replace('txt', 'Tms'), "a") as myfile:
		myfile.write(str(Results[-1][1])+"  "+str(Results[-1][0])+"\n")
	with open(sys.argv[1].replace('txt', 'Hills'), "a") as myfile:
		myfile.write(str(Results[-1][1])+"  "+str(Results[-1][2][0][0])+"  ("+str(Results[-1][2][1][0][0])+") \n")
	with open(sys.argv[1].replace('txt', 'Energy'), "a") as myfile:
		myfile.write(str(Results[-1][1])+"  "+str(Results[-1][3][0][0])+"  ("+str(Results[-1][3][1][0][0])+")  "+str(Results[-1][3][0][1]*Tenergy)+"  ("+str(Results[-1][3][1][1][1]*Tenergy)+")  "+str((Results[-1][3][0][0] - Results[-1][3][0][1]*Tenergy))+"\n")

with open(sys.argv[1].replace('txt', 'Tms'), "a") as myfile:
	myfile.write(str(np.mean(wartTm))+"  "+"average \n" + str(np.sqrt(np.var(wartTm)))+"  "+"stDev\n")
with open(sys.argv[1].replace('txt', 'Hills'), "a") as myfile:
	myfile.write(str(np.mean(wartHill))+"  "+"average \n" + str(np.sqrt(np.var(wartHill)))+"  "+"stDev\n")
with open(sys.argv[1].replace('txt', 'Energy'), "a") as myfile:
	myfile.write("H  "+str(np.mean(wartH))+"  "+"average \n" + str(np.sqrt(np.var(wartH)))+"  "+"stDev\n")
	myfile.write("TdS  "+str(np.mean(wartTdS))+"  "+"average \n" + str(np.sqrt(np.var(wartTdS)))+"  "+"stDev\n")
	myfile.write("G  "+str(np.mean(wartG))+"  "+"average \n" + str(np.sqrt(np.var(wartG)))+"  "+"stDev\n")


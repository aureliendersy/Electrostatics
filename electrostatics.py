#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  electrostatics.py
#  
#  Copyright 2019 Aurélien <aurelien@Aureliens-mac.local>
#  
#  
#  GOAL: Distribute charges on the lattice such that addition 
# of 1 electron gives a charge of -1 on the lattice site considered 
# and +1/L^2 for all sites. We then wish to look for the configurations 
# that minimize the potential energy at a 
# given number of added electrons. We only add electrons no positrons
import numpy as np


					


class Lattice:
	
	def __init__(self,fls):
		self.fls=fls		#fluid lattice spacing
		self.sls=fls*fls	#SU(N) lattice spacing
		self.lattice=[[0]*fls for i in range(fls)]
		self.energy=0
		print('Initialized the lattice')
	
	
	def display_lattice(self):
		for lines in self.lattice:
			for col in lines:
				if col>0: #Only display the positive charges
					print('|'+str(int(round(col))),end='')
				elif col==0:
					print('|'+str(0),end='')
				else:		
					print('|'+'~',end='')
	 
			print('|')
		print('\n')    	
		
	
	def add_charge(self,Q,x,y):
		for i in range(0,self.fls):
			for j in range(0,self.fls):
				if i==x and j==(self.fls-y-1):
					self.lattice[j][i]+=Q-np.sign(Q)*(1/self.sls)
				else:
					self.lattice[j][i]+=-np.sign(Q)*(1/self.sls)

		
		print('Added a charge centered at ('+str(x)+','+str(y)+')')
	
	def charge_at_site(self,x,y):
		print('Charge at ('+str(x)+','+str(y)+') is '+str(self.lattice[self.fls-y-1][x]))
	
	def charge_energy(self):
		etemp=0
		for i in range(0,self.fls):
			for j in range(0,self.fls):
				etemp+=np.power(abs(self.lattice[i][j]),2)
		return etemp*fund_energy(self)
		
	def update_energy(self):
		energy_temp=0
		for i in range(0,self.fls):
			for j in range(0,self.fls):
				for i2 in range(0,self.fls):
					for j2 in range(0,self.fls):
						if [i,j]!= [i2,j2]:
							energy_temp -=2*np.pi*self.lattice[i][j]*self.lattice[i2][j2]*np.log(norm_reshaped2(i,j,i2,j2,self)/self.sls)
		
		self.energy= (0.5)*energy_temp+ self.charge_energy()
		print('Energy updated')
	
	def show_energy(self):
		self.update_energy()
		print('Energy of the system is '+str(self.energy))	
							
def umklapp_coord(x,lat): #coord go from 0 to fls-1
	if x<0 or x>=lat.fls:
		print('Invalid coordinate for this lattice')
	else:
		return x-(lat.fls-1)/2
	
def norm_reshaped(x,y,lat):
	return np.linalg.norm([umklapp_coord(x,lat),umklapp_coord(y,lat)])

def mod_lattice(x,lat):
	if x>(lat.fls-1)/2:
		return x-lat.fls
	elif x<-(lat.fls-1)/2:
		return x+lat.fls
	else:
		return x
	
	
def norm_reshaped2(x1,y1,x2,y2,lat):
	return np.linalg.norm([mod_lattice(umklapp_coord(x1,lat)-umklapp_coord(x2,lat),lat),mod_lattice(umklapp_coord(y1,lat)-umklapp_coord(y2,lat),lat)])		

def fund_energy(lat):
	etemp=0
	for i in range(0,lat.fls):
		for j in range(0,lat.fls):
			if [umklapp_coord(i,lat),umklapp_coord(j,lat)]!=[0,0]:
				etemp+=1/(norm_reshaped(i,j,lat)*norm_reshaped(i,j,lat))
	return etemp	
			
def main():
	lattice1=Lattice(11)
	lattice1.add_charge(1,1,0)
	lattice1.add_charge(1,0,1)
	lattice1.display_lattice()
	lattice1.charge_at_site(1,0)
	lattice1.show_energy()
	print(lattice1.charge_energy())
	print('Efund= '+ str(fund_energy(lattice1)))
	return 0
    
    

if __name__ == '__main__':
    import sys
    sys.exit(main())

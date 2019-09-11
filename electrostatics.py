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
# given number of added electrons
import numpy as np


					


class Lattice:
	
	def __init__(self,fls):
		self.fls=fls		#fluid lattice spacing
		self.sls=fls*fls	#SU(N) lattice spacing
		self.lattice=[[0]*fls for i in range(fls)]
		
		print('Initialized the lattice')
	
	
	def display_lattice(self):
		for lines in self.lattice:
			for col in lines:
				if abs(col)>=1-self.fls/self.sls:
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
				if i==x and j==y:
					self.lattice[j][i]+=Q-np.sign(Q)*(1/self.sls)
				else:
					self.lattice[j][i]+=-np.sign(Q)*(1/self.sls)

		
		print('Added a charge centered at ('+str(x)+','+str(y)+')')
	 	
def main():
	lattice1=Lattice(11)
	lattice1.display_lattice()
	lattice1.add_charge(1,1,0)
	lattice1.add_charge(1,0,1)
	lattice1.display_lattice()
	return 0
    
    

if __name__ == '__main__':
    import sys
    sys.exit(main())

##
# BASED ON: ctml_writer.py
#
# Cantera .cti input file processor
#
# The functions and classes in this module process Cantera .cti input
# files and run a simulation. It can be imported as a module, or used
# as a script.
#
# script usage:
# 
# python cti2prime.py infile.cti
# 
# This will read in the cti file (default if not specified = 'chem.cti') 
# and produce a prime formatted thermo polynomial xml file for each species 
#
# Richard West 2009

from string import Template

try:
	import quantities as pq
except ImportError:
	print "ERROR - you need the Quantities package which can be found at http://packages.python.org/quantities/"
	raise
	
pq.UnitQuantity('kilocalories', pq.cal*1e3, symbol='kcal')
pq.UnitQuantity('kilojoules', pq.J*1e3, symbol='kJ')
pq.UnitQuantity('kilomoles', pq.mol*1e3, symbol='kmol')


# import all the classes from ctml_writer
#reload(ctml_writer)
from ctml_writer import *
import ctml_writer
# then we can modify the ones we want to do more than just write ctml


# define some more quantities units
_uConc=pq.Quantity(1,ctml_writer._umol) / pq.Quantity(1,ctml_writer._ulen)**3
_uTime=pq.Quantity(1,ctml_writer._utime)

# get the ctml_writer lists into the local namespace
_species=ctml_writer._species
_reactions=ctml_writer._reactions
_speciesnames=ctml_writer._speciesnames


class outsideValidRangeError(Exception):
	"""Was not within valid range for expression"""
	pass

# Classes should be named with CamelCase but need to stick with Cantera .ctml convention, hence lowercase:
class species(ctml_writer.species):
	"""a species. Inherit from the ctml_writer.species and add some methods."""
	def getCorrectNasaPoly(self,T):
		"""Returns whichever nasa polynomial is appropriate for a given temperature."""
		for nasapoly in self._thermo:
			if nasapoly.ValidTemperature(T): return nasapoly
		return None
		
	def getLowerNasaPoly(self):
		"""Get the lower (in terms of valid temperature range) of two Nasa polynomials"""
		if self._thermo[0]._t[0]<self._thermo[1]._t[0]:
			return self._thermo[0]
		else:
			return self._thermo[1]
	def getUpperNasaPoly(self):
		"""Get the upper (in terms of valid temperature range) of two Nasa polynomials"""
		if self._thermo[0]._t[0]<self._thermo[1]._t[0]:
			return self._thermo[1]
		else:
			return self._thermo[0]	
	def getLowerT(self):
		"""Get the lowest T covered by two Nasa polynomials"""
		return self.getLowerNasaPoly()._t[0]
	def getMiddleT(self):
		"""Get the changeover (middle) T for two Nasa polynomials"""
		assert self.getLowerNasaPoly()._t[1] == self.getUpperNasaPoly()._t[0], "Nasa polynomials don't meet at a single temperature"
		return self.getLowerNasaPoly()._t[1]
	def getUpperT(self):
		"""Get the highest T covered by two Nasa polynomials"""
		return self.getUpperNasaPoly()._t[1]
				
	def getEnthalpy(self,T):
		"""getEnthalpy(T) returns  Enthalpy at temperature T"""
		return self.getCorrectNasaPoly(T).getEnthalpy(T)
	
	def getThermo(self,T):
		"""Gets (HeatCapacityOverR, EnthalpyOverRT, EntropyOverR) at a given T"""
		thermo=self.getCorrectNasaPoly(T).getThermo(T)
		return thermo
		
	

class NASA(ctml_writer.NASA):
	"""NASA polynomial representation of thermo. Inherit from ctml_writer.NASA and add some methods."""
	def ValidTemperature(self,temperature):
		"""True if this polynomial is valid at a given temperature, else False"""
		if temperature<self._t[0]: return False
		if temperature>self._t[1]: return False
		return True

	def getThermo(self,T):
		"""getThermo(T) returns (HeatCapacityOverR, EnthalpyOverRT, EntropyOverR) for a given temperature T
		    Raises outsideValidRangeError exception if T is not within range of polynomial"""
		if not self.ValidTemperature(T): raise outsideValidRangeError
		if self._pref > 0.0: raise Exception("not sure what to do with customised standard state pressure")
		# http://www.me.berkeley.edu/gri-mech/data/nasa_plnm.html
		# Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
		# H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
		# S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
		c=self._coeffs
		HeatCapacityOverR=0.0
		EnthalpyOverRT=0.0
		EntropyOverR=0.0
		for i in range(5):
			HeatCapacityOverR+= c[i] * T**i
			EnthalpyOverRT+= c[i] * T**i / (i+1)
		EnthalpyOverRT+=c[5]/T
		EntropyOverR = c[0]*math.log(T) + c[1]*T + c[2]*T*T/2 + c[3]*T*T*T/3 + c[4]*T*T*T*T/4 + c[6]
		return(HeatCapacityOverR,EnthalpyOverRT,EntropyOverR)
		
	def getEnthalpy(self,T):
		"""getEnthalpy(T) returns Enthalpy at temperature T"""
		(HeatCapacityOverR,EnthalpyOverRT,EntropyOverR)=self.getThermo(T)
		Enthalpy=EnthalpyOverRT*pq.constants.R*T * pq.K  ## really T should have units already. oops! it's being passed around as a dimensionless float everywhere else in the code, so we just multiply by 1 Kelvin here. 
		Enthalpy.units='J/mol'
		return Enthalpy
		


if __name__ == "__main__":
    
	#reload(ctml_writer)
	
	# Load (run) the cantera input (.cti) file
	import sys, os
	if len(sys.argv)>1:
		file = sys.argv[1]
	else:
		print "using default file, chem.cti"
		file = 'chem.cti' # default input file
	base = os.path.basename(file)
	root, ext = os.path.splitext(base)
	dataset(root)
	if not _species:
		execfile(file)
	# we have now got the chem.cti file loaded

	# load the template
	XMLtemplate=Template( open('thp00000001.template.xml').read() )

	# iterate through the species
	for s in _species:
		# create a dictionary of values to go into the template
		d={
			'speciesname': s._name,
			'lowerT': s.getLowerT(),
			'middleT': s.getMiddleT(),
			'upperT': s.getUpperT()
		 }
		
		if s.getLowerT()<298:
			H= s.getEnthalpy(298)
			H.units='J/mol' # make sure units are right
			d['dfH298inJmol'] = float(H) # strip units
		elif s.getLowerT()<=300:
			H= s.getEnthalpy(300)
			H.units='J/mol' # make sure units are right
			d['dfH298inJmol'] = str(float(H)) + '<!-- polynomials invalid at 298K so this was evaluated at 300K -->' # strip units
		else:
			print "polynomials not valid at 298K so can't find standard enthlapy of formation"
			d['dfH298inJmol'] = '<!-- polynomials invalid at 298K so Hf unknown -->'
		
		for i,coeff in enumerate(s.getLowerNasaPoly()._coeffs):
			d['a%d'%(i+1)] = coeff
		for i,coeff in enumerate(s.getUpperNasaPoly()._coeffs):
			d['b%d'%(i+1)] = coeff
		
		# create output XML by substituting the dictionary d into the template
		output = XMLtemplate.substitute(d)
		print output
		#save it to a file
		outfile = open('thp00000001.%s.xml'%s._name,'w')
		outfile.write(output)
		outfile.close()



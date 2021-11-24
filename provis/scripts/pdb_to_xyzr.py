#!/usr/bin/python
# Python script based on msms pdb_to_xyzr  bash nawk script by Mike Pique
# File written by: Kristof Czirjak
#
# Extract x,y,z coordinates from PDB file and
# generate radius of for atom from atmtypenumbers.
# Hydrogens are presumed to be missing ("united atom" approach) unless -h given
# later: options for alternate radius and pattern files
# --- Mike Pique, The Scripps Research Institute
#
# input: bool: 1 to show solvent 0 to emit (HETATM)
# input: pdb file as argument or stdin
# output: new xyzr file to stdout
#
# Options:
#
# NOT SURE HOW THIS PART WORKS ####, from original file
# -h: use explicit hydrogen radii instead of default united atom radii
#
# examples:
#   python pdb_to_xyzr 1 crambin.pdb > crambin.xyzr
#	python pdb_to_xyzr 0 crambin.pdb > crambin.xyzr
#	
#

import sys
from linecache import getline
import re


pdb_file = str(sys.argv[2])
# Can easily modify to write to file directly
# out_file = "2fd7_new.xyzr"
solvent = int(sys.argv[1])
####
# if str(sys.argv[1]) != "-h":
# 	solvent = int(sys.argv[1])
# else:
# 	solvent = int(sys.argv[2])

numparams = len(sys.argv)
# set which field to use for radii:
h_select=5
####
# don't quite understand how this whole thing works
# implementation works, but commented out for now
# if numparams >= 2: #if number of parameters >= 1
# 	if str(sys.argv[1]) == "-h": 
# 		h_select=4
# 	# shift


# read radius table and patterns from supplied file
npats=0
atm_file = "/data/atmtypenumbers"
import os
new_path = os.getcwd() + atm_file

explicit_rad = [None] * 100
united_rad = [None] * 100
respat = [None] * 1000
atmpat = [None] * 1000
atmnum = [None] * 1000
# open file
with open(new_path) as atm_f:
	# for (line number, line) in file
	tempstr = ""
	for nr, line in enumerate(atm_f):
		# if line only has newline charactor or is comment
		if len(line) < 2 or line[0] == "#": 
			#skip
			continue
		# create string list from line
		line_list = line.split()
		# print(line_list)
		# if line starts with radius
		# set radius explicitly
		if line_list[0] == "radius":
			n = int(line_list[1]) #- 1
			#default case
			explicit_rad[n] = float(line_list[3])
			if len(line_list) <= 4 or line_list[4][0] == "#":
				united_rad[n] = float(explicit_rad[n])
			else:
				united_rad[n] = float(line_list[4])
			continue
		# residue name extracted
		respat[npats] = line_list[0]
		# no idea why conversions are needed/done
		if respat[npats] == "*":
			respat[npats] = ".*"
		# respat[npats] = "^" + respat[npats] + "$"
		tempstr = str(line_list[1]) #str("^" + line_list[1] + "$")
		tempstr = re.sub("_", " ", tempstr)
		atmpat[npats] = tempstr
		# atom number extracted
		atmnum[npats] = int(line_list[2])
		# print(atmnum[npats] , " and ", explicit_rad)
		# if current atom number not in atmnum list, fill with dummy
		if not (atmnum[npats] in atmnum):
			# the key has no radius --- complain and fake one
			# print("pdb_to_xyzr: error in library file" ,atm_file,
			# 	"entry ",line_list[0],line_list[1],line_list[2], "has no corresponding radius value" \
			# 	, "cat 1>&2") # write to stderr
			explicit_rad[atmnum[npats]] = 0.01
			united_rad[atmnum[npats]] = 0.01
		npats += 1
		# the above part runs the same as on linux (at least while testing)
		
with open(pdb_file) as pdb_f:
	for nr, line in enumerate(pdb_f):
		line_list = line.split()
		# we only care about "ATOM" and (potentially) "HETATM" lines
		if line_list[0]=="ATOM" or line_list[0]=="atom" or ((line_list[0]=="HETATM" or line_list[0]=="hetatm") and solvent):
			# save and store coordinates and name
			x = line_list[6]
			y = line_list[7]
			z = line_list[8]
			resname=line_list[3]

			# try to figure out the atom type
			aname = line_list[2]
			# special handling needed for hydrogens in PDB files: they start with
			# digits not the letter "H"
			A = line_list[2]
			B = re.compile("[1-9HhDd]+")
			m = B.search(A)
			if m is not None:
				aname="H"
			# However, some bogus PDP files have the H in column 13 so we allow
			# those too, which means we will treat as Hydrogen helium and hafnium 
			# but we protect HG ... ... mp
			A = line_list[2]
			B = re.compile("[Hh][^Gg]")
			m = B.search(A)
			if m is not None:
				aname="H"

			resnum=line_list[5]

			# trim any blanks
			resname = re.sub(" ", "", resname)
			aname = re.sub(" ", "", aname)
			# check for each atom and residue
			for pat in range(0, npats):
				# that we have a record for given atom in given residue
				m = re.search(atmpat[pat],aname)
				n = re.search(respat[pat],resname)
				if m is not None and n is not None: 
					break 
				pat+=1
			# with open(out_file, "a") as out_f:
			if pat==npats:
				# Not found
				# print("pdb_to_xyzr: error, file",pdb_file,"line",nr,"residue",
				# resnum,
				# "atom pattern",resname, aname,"was not found in ", pdb_file,
				# "cat 1>&2") # write to stderr
				newline = str(x) + " " + str(y) + " " + str(z) + " 0.01" #+ "\n"
				print(newline) #############
				# out_f.write(newline)
				pass
			else:
				temp = 0
				if h_select == 5:
					temp = united_rad[atmnum[pat]]
				else:
					temp = explicit_rad[atmnum[pat]]
				newline = str(x) + " " + str(y) + " " + str(z) + " " + str(temp) #+ "\n"
				print(newline) ###############
				# out_f.write(newline)

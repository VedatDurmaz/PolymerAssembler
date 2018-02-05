#!/usr/bin/env python

#########################################################################################################
#	 author: Vedat Durmaz
#							
#    polymer.py assembles a polymer from monomeric units defined in the dictionary @monomers and given as
#    pdb files and using a markov model for the frequency of their occurencies. the polymer is 
#    represented and visualized as a graph as well as saved as a pdb file with atom types according to
#    the amber force field which needs to modified accordingly: residuetypes.dat, aminoacids.hdb, 
#    aminoacids.rtp, specbonds.dat. one initial (root) and one ending (leaf) unit type must be specified.
#
#	command line argument of polymer.py:
#		<numUnits> = integer
#		<shuffleSites> = boolean: add units systematically (False) or randomly (True)
#		<forceLinearChain> = boolean: build globular (False) or linear (True) polymer
#         <pdbFilename> = string: name for pdb file and graphical output
#		<graphFilename> = string: name for graphical output (formats: eps, png, pdf)
#        
#         example:
#         "polymer.py 50 False False polymer.pdb polymer.eps"
#         "polymer.py 88 y f polymer.pdb polymer.eps"
#					
#	important global settings below (among others):
#		%monomers: dictionary of units,bindig sites (in/out), pdb file names
#         $fn<unit1>, $fn<unit2>, ...: path/name of monomers' pdb files
#         $mm<unit1>, $mm<unit2>, ...: transition probabilities from unit to others
#         $bondLen_cCOc, $angle_COC, $angle_OCC, $dihedAngle: geometry settings for new bonds
#		
#    additionally required packages:
#        python-pygraphviz
#        Biopython
#
#########################################################################################################

import sys,os
import numpy as np
import random as rd
import pygraphviz as pgv
import yaml
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
import lib.general_settings as genset





############## --------------------- subroutines ------------------------------

def checkArgs(argv):

	#argsList = ['<nUnits>', '<shuffleSites>', '<linear>', '<file.pdb>', '<file' + '/'.join(listFigure) + '>']
	argsList = ['<nUnits>', '<file.pdb>', '<file' + '/'.join(genset.listFigure) + '>']
	
	if not len(argv) == len(argsList)+1:
		sys.exit(" ERROR: wrong number of arguments. usage:\n " \
		     +sys.argv[0].split('/')[-1] + " " + ' '.join(argsList))
	elif not argv[2][-4:] == ".pdb":
		sys.exit(" ERROR: wrong file extension for pdb output. usage:\n " \
		     +sys.argv[0].split('/')[-1] + " " + ' '.join(argsList))
	elif argv[3][-4:] not in genset.listFigure:
		sys.exit(" ERROR: wrong file extension for graphical output. usage:\n " \
		     +sys.argv[0].split('/')[-1] + " " + ' '.join(argsList))
	else:
		try:
			argv = argv[:1] + map(int,[argv[1]]) + argv[2:]
		except ValueError:
			sys.exit(" ERROR: incorrect integer  \"" + str(argv[1]) + "\" ")
		else:
			return argv


### read pdb file using Bio.PDB module
def readFileBioPDB(struct,fn):
	p = PDBParser()
	try:
		with open(fn,'r') as infile:
			struct_in  = p.get_structure(struct, infile)
	except Exception, e:
		print " "+str(e)
		sys.exit(1)
	return struct_in


### write pdb file using Bio.PDB module
def writeFileBioPDB(struct,fn):
	io = PDBIO()
	io.set_structure(struct)
	try:	
		io.save(fn)
	except Exception, e:
		print " "+str(e)
		sys.exit(1)
	else:	
		testvar = 1


### read file as ASCII
def readFile(fn):
	try:
		data = [line.rstrip('\n') for line in open(fn,'r')]	
	except Exception, e:
		print " "+str(e)
		sys.exit(1)
	return data


### write file as ASCII
def writeFile(fn,data,mod='w'):
	try:
		with open(fn,mod) as outfile:
			outfile.write(data)
	except Exception, e:
		print " "+str(e)
		sys.exit(1)
	else:
		outfile.close()


### get markov model for polymer structure
def markovModel(usr_cfg): 
	
	MM = np.ones(len(usr_cfg['transmatrix'][usr_cfg['polymer']['type']]))
	for elem in usr_cfg['transmatrix'][usr_cfg['polymer']['type']]:
		MM = np.vstack([MM,elem])
	MM = np.delete(MM,(0),axis=0)
	
	for i in range(0,len(MM[:,1])-1):
		if np.sum(MM[i,:]) > 0:
			MM[i,:] = MM[i,:]/np.sum(MM[i,:])
		
	return MM


### build polymer with much more bifurcations
def polymerize(size,usr_cfg,pol_cfg):
	
	### get settings
	pol_type = usr_cfg['polymer']['type']
	shuffleSites = usr_cfg['polymer']['shuffleSites']
	forceLinearChain = usr_cfg['polymer']['forceLinearChain']
	code_unit = pol_cfg[pol_type]['code_unit']
	
	### get row-wisely cumulated HMM matrix
	MM = np.cumsum(markovModel(usr_cfg),1)
	
	### shuffle list of root monomer's binding sites twice or not!
	init_unit = pol_cfg[pol_type]['init_unit']
	srcList = pol_cfg[pol_type]['monomers'][init_unit]['out'][:]
	
	if shuffleSites in genset.listTrue:
		rd.shuffle(srcList)
		rd.shuffle(srcList)
	
	### create stack of added monomers' lists
	dictGraph = {1:{"type": init_unit, "par": 0, "site": "", "child": set([])}}
	numMonomers = 1
	
	### append all binding sites of current unit to list
	stackRoot = []
	for src in srcList:   ## TODO:  BEFORE using unshuffled: pol_cfg[pol_type]['monomers'][init_unit]['out']: 
		stackRoot.append([numMonomers,init_unit,src])
		
	
	while any([len(stackRoot)>0, numMonomers<size]):

		if shuffleSites in genset.listTrue:
			rd.shuffle(stackRoot)
			rd.shuffle(srcList)

		if forceLinearChain in genset.listTrue:
			bondPar, bondUnit, bondSrc = stackRoot.pop()
		else:
			#print stackRoot
			bondPar, bondUnit, bondSrc = stackRoot.pop(0)
		
		### next random item 
		if numMonomers + len(stackRoot) + 1 >= size:
			childUnit = pol_cfg[pol_type]['end_unit']
		else:
			rand = rd.random()
			idxUnit = code_unit.index(bondUnit)
			
			### apply markov model
			for i in range(len(MM[idxUnit])-1):
				if rand <= MM[idxUnit][i+1]:
					childUnit = code_unit[i+1]
					break

		### append new monomer to polymer graph
		numMonomers += 1
		dictGraph[numMonomers] = {"type":childUnit,"par":bondPar,"site":bondSrc,"child": set([])}
		
		### add child in parent's child list
		dictGraph[bondPar]["child"].add(numMonomers)
		
		### update stack
		srcList = pol_cfg[pol_type]['monomers'][childUnit]['out'][:]
		
		### adding all binding sites of new monomer to list
		if shuffleSites in genset.listTrue:
			rd.shuffle(srcList)
			rd.shuffle(srcList)
		for src in pol_cfg[pol_type]['monomers'][childUnit]['out']:
			stackRoot.append([numMonomers,childUnit,src])

	if len(stackRoot) != 0:
		print " !! ATTENTION: " + str(len(stackRoot)) + " free valencies found!"
	elif numMonomers != size:
		print " !! ATTENTION: " + str(numMonomers) + " units added instead of " + str(size) + "!"
	
	return dictGraph


### labels combined with '-' used for "resID-atom" and "elem-elem"
def joinLabel2(lab1,lab2): return str(lab1) +genset.labelSep + str(lab2)
def splitLabel(lab): return str(lab).split(genset.labelSep)


### build rotation matrix for axis-angle notation
def setAxAng(axis,angle):
	
	x = axis[0]
	y = axis[1]
	z = axis[2]
	
	cos  = np.cos(angle)
	sin  = np.sin(angle)
	cos1 = 1 - cos
	
	row1 = np.array( (cos+(x**2)*cos1, x*y*cos1-z*sin, x*z*cos1+y*sin) )
	row2 = np.array( (y*x*cos1+z*sin, cos+(y**2)*cos1, y*z*cos1-x*sin) )	
	row3 = np.array( (z*x*cos1-y*sin, z*y*cos1+x*sin, cos+(z**2)*cos1) )	

	return np.array([row1,row2,row3])


### find all atoms connected to $start using depth first search 
def dfs(graph, start, done=None):
	if done is None:
		done = set()
	done.add(start)
	#for next in graph[start] - done:	
	for next in graph[start]["child"] - done: #subtraction of set elements!!
		dfs(graph, next, done)
	return done


### get cross product of (v,w); normalized if desired (,,True)
def getCrossProd(v,w,norm=False):
	
	cross = np.cross(v,w)
	if norm:
		return cross/np.linalg.norm(cross)
	else:
		return cross


### get inner product of (v,w); normalized if desired (,,True)
def getInnerProd(v,w,norm=False):
	
	dot = np.dot(v,w)
	if norm:
		return dot/np.linalg.norm(v)/np.linalg.norm(w)
	else:
		return dot
		
		
def get_angle(a, b):
	w = getCrossProd(a,b)
	wlen  = np.linalg.norm(w)
	s = getInnerProd(a,b)
	return( np.arctan2(wlen,s) )


### compute dihedral angle
def get_dihedral(q1,q2,q3,q4):
	r_ij = q2-q1
	r_kj = q2-q3
	r_kl = q4-q3

	m = getCrossProd(r_ij, r_kj)
	n = getCrossProd(r_kj, r_kl)
	phi = get_angle(m,n)
	ipr = getInnerProd(r_ij,n)
	sign= np.sign(ipr)
	if sign == 0:
		sign = 1
	return( np.array(sign*phi) )


### rotate dihedral of res atoms to angle $phi
def rotateDihed(center,R,res):

	for a in res:
		a.transform(genset.unity3, -center)
		a.set_coord(np.dot(R,a.get_coord()))
		a.transform(genset.unity3, center)
	return res


### bend residue $res around $axis through atom $q2 by $angle
def bendAngle(center,R,res):
	
	for a in res:						
		a.transform(genset.unity3, -center)
		a.set_coord(np.dot(R,a.get_coord()))
		a.transform(genset.unity3, center)
	return res


### build next unit
def buildMonomer(res,resPar,a1,a2,a3,a4,a5,a6,pol_geom):

	bondLen_cCOc = pol_geom['bondLen_cCOc']
	angle_COC = pol_geom['angle_COC']
	angle_OCC = pol_geom['angle_OCC']
	dihedAngle = pol_geom['dihedAngle']
	
	q1 = resPar[a1].get_coord()
	q2 = resPar[a2].get_coord()
	q3 = resPar[a3].get_coord()
	
	### translate to parent's link
	bond = bondLen_cCOc * (q2-q1) / np.linalg.norm(q2-q1)
	q4 = np.array(res[a4].get_coord())
	for a in res:
		a.transform(genset.unity3, q3+bond-q4)			
	
	### compute first dihedral angle (related to new bond) and rotate to trans
	q4 = res[a4].get_coord()
	phi = get_dihedral(q1,q2,q3,q4)
	axis = (q3-q2)/np.linalg.norm(q3-q2)
	R = setAxAng(axis,phi-dihedAngle*np.pi/180)
	res = rotateDihed(q3,R,res)
	
	### compute second dihedral angle and rotate to trans
	q4 = res[a4].get_coord()
	q5 = res[a5].get_coord()
	phi = get_dihedral(q2,q3,q4,q5)	
	axis = (q4-q3)/np.linalg.norm(q4-q3)
	R = setAxAng(axis,phi-dihedAngle*np.pi/180)
	res = rotateDihed(q4,R,res)
	
	### compute third dihedral angle and rotate to trans
	q4 = res[a4].get_coord()
	q5 = res[a5].get_coord()
	q6 = res[a6].get_coord()
	phi = get_dihedral(q3,q4,q5,q6)	
	axis = (q5-q4)/np.linalg.norm(q5-q4)
	R = setAxAng(axis,phi-dihedAngle*np.pi/180)
	res = rotateDihed(q5,R,res)
	
	### compute first bond angle and bend angle around axis $getCrossProd(r_ij, r_kj)
	q4 = res[a4].get_coord()
	theta = get_angle(q2-q3, q4-q3)
	axis = getCrossProd(q2-q3, q4-q3, True)
	R = setAxAng(axis,angle_COC*np.pi/180-theta)
	res = bendAngle(q3,R,res)
	
	### compute second bond angle and bend angle around axis $plane2
	q4 = res[a4].get_coord()
	q5 = res[a5].get_coord()
	theta = get_angle(q3-q4, q5-q4)
	axis = getCrossProd(q3-q4, q5-q4, True)
	R = setAxAng(axis,angle_OCC*np.pi/180-theta)
	res = bendAngle(q4,R,res)


### draw graph: import pygraphviz as pgv
#http://pygraphviz.github.io/documentation/latest/tutorial.html
def drawGraphviz(graph,outfile,code_colors,dirG,labelG):
	
	G = pgv.AGraph(strict=True,directed=dirG,label=labelG,landscape='False',ranksep='0.5')
	#G.graph_attr.update(landscape='True')
	
	for node in sorted(graph.keys()):
		G.add_node(node, label='',color=code_colors[graph[node]["type"]],style='filled')
		#G.add_node(node, label=str(graph[node]["type"]),color=codeColors[graph[node]["type"]],style='filled')
		par = graph[node]["par"]		
		if par > 0:
			weight = 1
			if G.get_node(par).attr['label'] == graph[node]["type"]:
				weight = 10
			site = graph[node]["site"].split('-')[2]
			#G.add_edge(par,node,label=site,color=codeColors[site])
			G.add_edge(par,node,weight=weight,label='',color=code_colors[site])	
	G.layout()
	#G.draw(outfile,prog='dot')
	G.draw(outfile)


	
############## --------------------- main ------------------------------

if __name__ == '__main__':
	
	### check and assign argument list
	argv = checkArgs(sys.argv)
	maxMonomers = argv[1]
	fn_out = argv[2]
	fn_out_viz = argv[3]
	scriptPath = os.path.abspath(os.path.dirname(sys.argv[0]))

	### read user config file
	with open(os.path.join(scriptPath,genset.usr_cfg_fn), 'r') as usr_file:
		usr_cfg = yaml.load(usr_file)
	
	### read user config file
	with open(os.path.join(scriptPath,genset.pol_cfg_fn), 'r') as pol_file:
		pol_cfg = yaml.load(pol_file)
	
	pol_type = usr_cfg['polymer']['type']
	pol_geom = pol_cfg[pol_type]['geometry']
	
	fn_tmp = "temp.pdb"
	if os.path.isfile(fn_out):
		os.remove(fn_out) 

	### compute graph of monomers 
	print "\n building graph of polymer topology ... "
	graph = polymerize(maxMonomers,usr_cfg,pol_cfg)
	
	### assemble polymer according to dfs-list of monomers as defined in @monomers
	print " computing polymer coordinates ... "
	count = 0
	polymer = {}
	for i in dfs(graph, 1):
		count += 1
		
		### read current pdb unit, store in polymer{} and set residue number
		unit = graph[i]["type"]
		polymer[i] = readFileBioPDB(unit, os.path.join(scriptPath,pol_cfg[pol_type]['pdb_fn'][unit]))
		res = polymer[i].get_residues().next()
		res.id = (" ",count," ")

		### initial unit (1) doesn't need to be translated/rotated
		if count > 1:
			
			### load parental monomer
			par = graph[i]["par"]
			resPar = polymer[par].get_residues().next()

			### start building with atoms involved in 
			a1,a2,a3 = splitLabel(graph[i]["site"])
			#a4,a5,a6 = splitLabel(monomers[unit]["tgt"][0])
			
			a4,a5,a6 = splitLabel(pol_cfg[pol_type]['monomers'][unit]['in'][0])
			
			buildMonomer(res,resPar,a1,a2,a3,a4,a5,a6,pol_geom)
			
			#print(a1,a2,a3,a4,a5,a6)
			
			
			
			
		### append current monomer to pdb file
		writeFileBioPDB(polymer[i],fn_tmp)
		lines = file(fn_tmp, 'r').readlines()
		if any([lines[-1] == "TER\n", lines[-1] == "END\n"]):
			del lines[-1]
		if any([lines[-1] == "TER\n", lines[-1] == "END\n"]):
			del lines[-1]
		file(fn_out, 'a').writelines(lines)
		#print " monomer " + str(i)	+ " (" + unit + ") built"
	
	
	drawGraphviz(graph,fn_out_viz, pol_cfg[pol_type]['code_colors'],usr_cfg['graph']['dirG'],usr_cfg['graph']['labelG'])
			
	if os.path.exists(fn_tmp):
		os.remove(fn_tmp)
	
	print " successfully added " + str(count) + " residues ... adieu!"
			

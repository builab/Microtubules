#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Using RANSAC to fit line/slightly curve line on centered helical coordinate
# Image option will miss the 1st micrograph

"""
Created on Sat Jun  6 17:35:42 2020

@author: kbui2

Using zhefan's RANSAC code to fit
lost of stuff to fix
"""

import os, sys, argparse, os.path, glob, math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from utils.parsers import *
from utils.filament_fit import *

def learnstarheader(infile):
	"""Learn which column contains which information from an already open starfile"""
	infile.seek(0) # Go to the beginning of the starfile
	doneheader = False
	doneprelabels = False
	headerlabels = []
	while not doneprelabels:
		line=infile.readline()
		if line.startswith('loop_'):
			doneprelabels = True # read until 'loop_'
	while not doneheader:
		line=infile.readline()
		if not line.startswith('_'): # read all lines the start with '_'
			doneheader = True
		else:
			headerlabels += [line] 
	infile.seek(0) # return to beginning of starfile before return
	return headerlabels

def writestarheader(outfile,headerlabels):		  
	"""With an already opened starfile write a header"""
	outfile.write('\ndata_\n\nloop_\n')
	for label in headerlabels:
		outfile.write(label)

def readstarline(infile):
	"""Read a record (line) from an already open starfile and return XXX"""
	line=infile.readline()
	records = line.split()
	return records

def writestarline(outfile,records):
	"""Write a record (line) to an already open starfile"""
	for item in records:
		outfile.write(item+'  ')
	outfile.write('\n')

def starcol_exact_label(starlabels, label):
	"""New function to do exact match of relion label such as _rlnImageCol"""
	for i, s in enumerate(starlabels):
		record=s.split()
		if label == record[0]:
			return i
	return -1

def writestarblock(outfile,recordblock):
	"""Write a record (line) to an already open starfile"""
	for record in recordblock:
		writestarline(outfile, record)

def interpol_helix(helicalrecord, binfactor, spacing, helicalid):
	"""Interpolate to the periodicity of dist"""
	#print(helicalid)
	#print(*helicalrecord)
	tiltlist = np.array([])
	spacingpixel = spacing*binfactor/(float(helicalrecord[0][detpixelsizecol])*10000/float(helicalrecord[0][magcol]))

	for i in range(len(helicalrecord)):
		tiltlist = np.append(tiltlist, [float(helicalrecord[i][tiltcol])])
		shiftX = float(helicalrecord[i][originxcol])*binfactor
		shiftY = float(helicalrecord[i][originycol])*binfactor
		if i == 0:
			origin = np.array([float(helicalrecord[i][coordxcol]) + shiftX, float(helicalrecord[i][coordycol]) + shiftY])
		else:
			origin = np.vstack((origin, [float(helicalrecord[i][coordxcol]) + shiftX, float(helicalrecord[i][coordycol]) + shiftY]))



	try:
		poly_o = ransac_fit.polyfit(origin, 2, 1, disable_linear=False, directory_mode=False)
		arclength_o = ransac_fit.arclength(poly_o)
		x = ransac_fit.spacing(arclength_o, spacingpixel)

		y = poly_o["model"].predict(x)
		x = [(item[0]) for item in x]

	     # swap the x and y predicted if they were intially swapped by polyfit
		if (poly_o["swapped"] == True):
			tmp = x
			x = y
			y = tmp
	except ValueError:
		print ('No problem')	
		return
		
	psi = calculatepsi(x, y)
	medtilt = np.median(tiltlist)

	partandstack=helicalrecord[0][imagecol].split('@')

	
	fittedhelicalrecord = []
	for npart in range(len(x)):
			#print(npart)
			singledataline = helicalrecord[0].copy()
			singledataline[coordxcol] = "{:.6f}".format(x[npart]) 
			singledataline[coordycol] = "{:.6f}".format(y[npart]) 
			# New tracklength, image & psi
			singledataline[helicaltracklengthcol] = "{:.6f}".format(spacingpixel*npart) 
			singledataline[psicol] = "{:.6f}".format(psi[npart]) 
			singledataline[psiflipratiocol] = "{:.6f}".format(0.5) 
			singledataline[psipriorcol] = "{:.6f}".format(psi[npart]) 
			singledataline[tiltpriorcol] = "{:.6f}".format(medtilt) 
			singledataline[tiltcol] = "{:.6f}".format(medtilt) 
			singledataline[originxcol] = "{:.6f}".format(0) 
			singledataline[originycol] = "{:.6f}".format(0) 
			singledataline[imagecol]=str(npart + 1).zfill(6)+'@'+partandstack[1]
			fittedhelicalrecord.append(singledataline)
	#print(*fittedhelicalrecord, "\n")
	return fittedhelicalrecord
		

def calculatepsi(x, y):
	"""Calculate the psi from the tangent"""
	xd = np.diff(x)
	yd = np.diff(y)
	psirad = np.arctan2(yd,xd)*-1
	psi = np.array(psirad)*180/math.pi
	psi = np.hstack([psi, [psi.max()]])
	return psi
	#print(psi)

	
if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	print('WARNING: Only compatible with Relion 3.0 star file')
	
	parser = argparse.ArgumentParser(description='Plot coordinate of star file')
	parser.add_argument('--istar', help='Input particle star file',required=True)
	parser.add_argument('--ostar', help='Output particle star file',required=True)
	parser.add_argument('--spacing', help='Distance in pixel',required=True)
	parser.add_argument('--ibin', help='Bin in current star file',required=True,default=5.079)
	parser.add_argument('--minpart', help='Minimum number of particles for fitting',required=False,default=5)
	parser.add_argument('--im', help='Directory for output fitted image',required=False,default="")


	#parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	
	
	instar = open(args.istar, 'r')
	outstar= open(args.ostar, 'w')
	binfactor = float(args.ibin)
	spacing = float(args.spacing)
	minpart = float(args.minpart)
	
		
	
	starlabels = learnstarheader(instar)
	coordxcol = starcol_exact_label(starlabels, '_rlnCoordinateX')
	coordycol = starcol_exact_label(starlabels, '_rlnCoordinateY')
	originxcol = starcol_exact_label(starlabels, '_rlnOriginX')
	originycol = starcol_exact_label(starlabels, '_rlnOriginY')
	microcol = starcol_exact_label(starlabels, '_rlnMicrographName')
	imagecol = starcol_exact_label(starlabels, '_rlnImageName')
	helicalidcol = starcol_exact_label(starlabels, '_rlnHelicalTubeID')
	helicaltracklengthcol = starcol_exact_label(starlabels, '_rlnHelicalTrackLength')
	psicol = starcol_exact_label(starlabels, '_rlnAnglePsi')
	psiflipratiocol = starcol_exact_label(starlabels, '_rlnAnglePsiFlipRatio')
	psipriorcol = starcol_exact_label(starlabels, '_rlnAnglePsiPrior')
	tiltpriorcol = starcol_exact_label(starlabels, '_rlnAngleTiltPrior')
	tiltcol = starcol_exact_label(starlabels, '_rlnAngleTilt')
	rotcol = starcol_exact_label(starlabels, '_rlnAngleRot')
	dfucol = starcol_exact_label(starlabels, '_rlnDefocusU')
	dfvcol = starcol_exact_label(starlabels, '_rlnDefocusV')
	dfacol = starcol_exact_label(starlabels, '_rlnDefocusAngle')
	magcol = starcol_exact_label(starlabels, '_rlnMagnification')
	detpixelsizecol = starcol_exact_label(starlabels, '_rlnDetectorPixelSize')

	writestarheader(outstar, starlabels)
	
	fig1 = plt.figure(1)
	
	helicalid = 0
	microlist ={}
	prevhelicalid = 0
	micronum = 1
	helicalrecord = []
	for line in instar:
		record = line.split()
		if len(record)==len(starlabels): # if line looks valid
			microname=record[microcol]
			microname = os.path.basename(microname)
			# Create a dictionary, if microname not exist as key, insert
			# Parsing each helicaid into a block of record
			if microlist.get(microname):
				if prevhelicalid != record[helicalidcol]:
					if len(helicalrecord) > minpart:
						fittedhelicalrecord = interpol_helix(helicalrecord, binfactor, spacing, helicalid)
						writestarblock(outstar, fittedhelicalrecord)
						# image option
						ax1 = None
						if args.im != "":
							ax1 = fig1.gca()
							x = [float(row[coordxcol]) for row in fittedhelicalrecord]
							y = [float(row[coordycol]) for row in fittedhelicalrecord]
							ax1.plot(x, y, color="purple")
							orix = [float(row[0]) for row in helicalrecord]
							oriy = [float(row[1]) for row in helicalrecord]
							ax1.scatter(orix, oriy)
					
					helicalid += 1
					prevhelicalid = record[helicalidcol]
					helicalrecord = []
					helicalrecord.append(record)
				else:
					helicalrecord.append(record)

			else:
				microlist[microname] = micronum
				print(" Processing Micrograph ", micronum)
				if micronum > 1:
					if len(helicalrecord) > minpart:
						fittedhelicalrecord = interpol_helix(helicalrecord, binfactor, spacing, helicalid)
						writestarblock(outstar, fittedhelicalrecord)
						if args.im != "":
							ax1 = fig1.gca()
							x = [float(row[coordxcol]) for row in fittedhelicalrecord]
							y = [float(row[coordycol]) for row in fittedhelicalrecord]
							ax1.plot(x, y, color="purple")
							orix = [float(row[0]) for row in helicalrecord]
							oriy = [float(row[1]) for row in helicalrecord]
							ax1.scatter(orix, oriy)
							ax1.set_title(microname)
							print('Writing fitted plot for micrograph ', micronum)
							fig1.savefig(args.im + "/" + microname.replace(".mrc", ".png"))
							ax1 = None
					# Write out
					helicalrecord = []
					
				micronum += 1
				helicalid += 1
				prevhelicalid = record[helicalidcol]
				helicalrecord.append(record)

			#if helicalid > 5:
			#	break

			
	instar.close()
	outstar.close()

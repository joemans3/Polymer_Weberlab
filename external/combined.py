import numpy as np
import sys
import os
from glob import glob
from os.path import dirname, basename
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d as ax3d
import matplotlib.animation
import matplotlib.colors as colors
import math
from pathlib import Path


from cmpaircorrelation import *
from scipy.optimize import curve_fit
import shutil
import pdb


ID_name = []
ID_stack = []
gt_temp = []
gt_ms = []
gt_msc = []
gt_rh = []

def load_data(path):
	ID_name.append(basename(path))
	dirs = glob("{0}/*/".format(path))
	dirs2 = []
	for k in dirs:
		dirs2.append(k[:-1])
	dirs = dirs2

	ID_stack.append([basename(entry) for entry in dirs])

	g_temp = []
	g_ms = []
	g_msc = []
	g_rh = []

	for i in range(len(dirs)):

		temper_l = []
		msd_s_l = []
		msd_s_cm_l = []
		rh_l = []



		dirs_T = glob("{0}/*/".format(dirs[i]))

		for j in range(len(dirs_T)):
			dirs_T[j] = dirs_T[j][:-1]
			temper_l.append(basename(dirs_T[j]))
			msd_s_l.append(np.loadtxt(dirs_T[j]+"/s_msd_s_meana.txt"))
			msd_s_cm_l.append(np.loadtxt(dirs_T[j]+"/s_msd_s_meana_cm.txt"))
			rh_l.append(np.loadtxt(dirs_T[j]+"/s_radius_holder.txt"))


		g_temp.append(temper_l)
		g_ms.append(msd_s_l)
		g_msc.append(msd_s_cm_l)
		g_rh.append(rh_l)
	gt_temp.append(g_temp)
	gt_ms.append(g_ms)
	gt_msc.append(g_msc)
	gt_rh.append(g_rh)

	return 


def msd_line(D_A,length):
	return D_A[0]*np.array(list(range(length)))**D_A[1]

def runit():
	for i in range(len(ID_name)):
		for j in range(len(ID_stack[i])):
			for k in range(len(gt_temp[i][j])):
				plt.plot(list(range(1,len(gt_ms[i][j][k])+1)),gt_ms[i][j][k],label = "{0}".format(ID_stack[i][j]+"_"+gt_temp[i][j][k]))
	plt.plot(msd_line([0.07,1.0],300),label = "a=1")
	plt.plot(msd_line([0.05,0.45],300),label = "a=0.6")
	plt.xscale("log")
	plt.yscale("log")
	plt.legend()
	plt.show()
	return 

def runit_dt():
	for i in range(len(ID_name)):
		for j in range(len(ID_stack[i])):
			ar =np.array(gt_ms[i][j])[:,0]
			art = np.array(gt_temp[i][j],dtype = float)
			plt.plot(art,ar,'.',label = "{0}".format(ID_stack[i][j]))
	plt.xscale("log")
	plt.yscale("log")
	plt.legend()
	plt.show()
	return 
#*int(ID_stack[i][j][2:])
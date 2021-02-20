from __future__ import print_function
import numpy as np
import rmsd
import sys, os, glob

def rotation_matrix(sigma):
    radians = sigma * np.pi / 180.0
    r11 = np.cos(radians)
    r12 = -np.sin(radians)
    r21 = np.sin(radians)
    r22 = np.cos(radians)
    return np.array([[r11, r12], [r21, r22]])
def getxyzs(filename):
    xyz = []
    xyz_noH = []
    lines = open(filename,"r").readlines()
    for i in range(2,len(lines)):
        vector = np.array([float(lines[i].split()[1]),
                           float(lines[i].split()[2]),
                           float(lines[i].split()[3])])
        if lines[i].split()[0] != "H":
            xyz_noH.append(vector)
        xyz.append(vector)
    return np.array(xyz), np.array(xyz_noH), float(lines[1].split(";")[0])

molnames = []
molxyzs = []
molxyzs_noH = []
molenergies = []
for file in glob.glob("./*xyz"):
    A, AnoH, ener = getxyzs(file)
    A -= rmsd.centroid(A)
    AnoH -= rmsd.centroid(AnoH)
    molnames.append(file)
    molxyzs.append(A)
    molxyzs_noH.append(AnoH)
    molenergies.append(ener)

loglines = []
for i in range(len(molxyzs)):
    for j in range(i+1,len(molxyzs)):
        U_full = rmsd.kabsch(molxyzs[i], molxyzs[j])
        Temp_full = np.dot(molxyzs[i], U_full)
        U_noH = rmsd.kabsch(molxyzs_noH[i], molxyzs_noH[j])
        Temp_noH = np.dot(molxyzs_noH[i], U_noH)
        loglines.append("%s,%s,%f,%f,%f" % (molnames[i], molnames[j],
                                         (molenergies[i]-molenergies[j])*627.5094740,
                                         rmsd.rmsd(Temp_full, molxyzs[j]),
                                         rmsd.rmsd(Temp_noH, molxyzs_noH[j])))


outfile = open("dataset.csv", "w")
outfile.write("\n".join(loglines))
outfile.close()
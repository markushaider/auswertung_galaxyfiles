#!/usr/bin/python

#clusterNumbers = ['600', '601', '602', '603', '604', '605', '607', '608', '609', '610']
clusterNumbers = ['602']

lpath='/RAID3/markus/clusterdata/galaxy_files/auswertung/'
spath='/RAID3/markus/clusterdata/galaxy_files/auswertung/'

#Welche Videos
video=('diskmetal_histo','distance_histo','mcold_histo',\
       'pos3d','pos_x','pos_y','pos_z','rdisk_histo',\
       'stf_histo','vel_x_histo','vel_y_histo','vel_z_histo','vel_histo')

for name in video:
	for cNumber in clusterNumbers:
		from subprocess import call
		call(["mencoder",
		      "mf://"+lpath+cNumber+"/plots/"+name+"*.png",
		      "-mf",
		      "type=png:w=800:h=600:fps=15",
		      "-ovc",
		      "lavc",
		      "-lavcopts",
		      "vcodec=mpeg4",
		      "-oac",
		      "copy",
		      "-o",
		      spath+name+".avi"]

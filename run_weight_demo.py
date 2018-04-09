#__________________________________#
#      GLaSS                       #
#      Peter Taylor                #
#      Mullard Space Laboratory    #
#      University College London   #
#      2018                        #
#__________________________________#



import sys
import os
import GLaSS
import numpy as np


#__________________________________#
# setup and run cosmosis pipeline  #
#__________________________________#



cosmosis_path = os.environ['COSMOSIS_SRC_DIR']
sys.path.append(cosmosis_path)
#you can also hard code the path to cosmosis
#in manually if above doesn't work
#sys.path.append('/Users/peter/codes/cosmosis/')
from cosmosis.runtime.config import Inifile
from cosmosis.runtime.pipeline import LikelihoodPipeline
cwd = os.getcwd()
file = cwd + '/demo-spherical-bessel.ini'
ini = Inifile(file)
pipeline = LikelihoodPipeline(ini)
pipeline.quiet = True
pipeline.debug = False
pipeline.timing = True


data = pipeline.run_parameters([])




#__________________________________#
# get back the data                #
#__________________________________#

k_samps = data["3d_weak_lensing_output", "k_samps"]
r_samps = data["3d_weak_lensing_output", "r_samps"]
l_input_array = data["3d_weak_lensing_output", "l_input_array"]
resolution = 200




#__________________________________#
# make the weights                 #
#__________________________________#


# use spherical-bessel for example #
# note the cosmology is not the
# same as for the fiducial bessel
# weight so will not be exactly the
# same as default spherical-bessel
from scipy.special import spherical_jn


print 'computing and saving the weights'
for li in l_input_array:
	weight = np.zeros((resolution, resolution))
	for i in range(resolution):
		for j in range(resolution):
			weight[i,j] = spherical_jn(li, k_samps[j] * r_samps[i])
	np.savetxt('weights/weight_%s.txt' %li, weight)



#__________________________________#
# run cosmosis with weighted 
# lensing
#__________________________________#

file = cwd + '/demo-weighted-lensing.ini'
ini = Inifile(file)
pipeline = LikelihoodPipeline(ini)
pipeline.quiet = True
pipeline.debug = False
pipeline.timing = True


data = pipeline.run_parameters([])




#__________________________________#
# get back the data                #
#__________________________________#

c_l_storage_array = data["3d_weak_lensing_output", "c_l_storage_array"]
shot_noise_storage_array = data["3d_weak_lensing_output", "shot_noise_storage_array"]
k_samps = data["3d_weak_lensing_output", "k_samps"]
z_samps = data["3d_weak_lensing_output", "z_samps"]
l_input_array = data["3d_weak_lensing_output", "l_input_array"]








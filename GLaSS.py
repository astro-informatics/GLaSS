#__________________________________#
#      GLaSS                       #
#      Peter Taylor                #
#      Mullard Space Laboratory    #
#      University College London   #
#      2018                        #
#__________________________________#




import numpy as np
from cosmosis.datablock import names, option_section
import lensing_calculation
import os



#__________________________________#
# This will just run once in       
# cosmosis.  Load in run-mode
# options + calculate Bessel 
# data
#
#__________________________________#

def setup(options):
    
    print "You are running GLaSS!"




    #__________________________________#
    # custom weights on or off   
    # tomographic should be false 
    # to use these
    #__________________________________#

 
    try:
        weighted_lensing = options.get_bool(option_section, "weighted_lensing")
    except:
        weighted_lensing = False
    if weighted_lensing == True:
        print "You are using weighted lensing"
    else:
        pass



    #__________________________________#
    # curved or flat cosmology         #
    #__________________________________#


    try:
        curved_lens_kernel = options.get_bool(option_section, "curved_lens_kernel")
        if curved_lens_kernel == True:
            print "You are using the lensing kernel to account for curved cosmology"
        else:
            curved_lens_kernel = False
    except:
        curved_lens_kernel = False
    try: 
        omega_k = options.get_double(option_section, "omega_k")
    except:
        omega_k = 0.
    try:
        hyper_spherical = options.get_bool(option_section, "hyper_spherical")
    except:
        hyper_spherical = False
    if hyper_spherical == True:
        print "You are using hyper_spherical bessel funcitons"
    else:
        pass
    



    #__________________________________#
    # l-modes to sample                #
    # l_max must be larger than        # 
    #largest element of l_input_array  #
    #__________________________________#


    try:
        l_input_array = options.get_int_array_1d(option_section, "l_input_array")
    except:
        l_input_array = 100
        print "CAUTION: You have not specified which l-modes to sample"
    try:
        l_max = options.get_int(option_section, "l_max")
    except:
        l_max = 1000
        print "CAUTION: You have not specified a value for l_max"
    


    #__________________________________#
    # computation grid res             #
    #__________________________________#

    try:
        resolution = options.get_int(option_section, "resolution")
    except:
        resolution = 500
        print "CAUTION: You have not specified a value for resolution"
    resolution_sub_samps = resolution 

    


    #__________________________________#
    # bessel grid resolution
    # probably not a good idea to touch 
    # these
    #__________________________________#

    try:
        resolution_bes = options.get_int(option_section, "resolution_bes")
    except:
        resolution_bes = 30000
    try:
        resolution_bes_lim = options.get_int(option_section, "resolution_bes_lim")
    except:
        resolution_bes_lim = 30000
    try:
        tolerance = options.get(option_section, "tolerance")
    except:
        tolerance = 1e-10

    


    #__________________________________#
    # photometric redshift error
    # parameters
    # 
    #__________________________________#


    try:
        c_cal = options.get_double(option_section, "c_cal")
    except:
        c_cal = 1.
        print "CAUTION: You have not specified a value for c_cal"
    try:
        z_bias = options.get_double(option_section, "z_bias")
    except:
        z_bias = 0.
        print "CAUTION: You have not specified a value for z_bias"
    try:
        A_sigma = options.get_double(option_section, "A_sigma")
    except:
        A_sigma = 0.05
        print "CAUTION: You have not specified a value for A_sigma"





    #__________________________________#
    # n(z) parameters.
    # If you want to use custom n(z)
    # see the comment above 
    # redshift_density_function function
    # in lensing_calculation file
    #__________________________________#

    try:
        a1 = options.get_double(option_section, "a1")
    except:
        a1 = 1.5
    try:
        c1 = options.get_double(option_section, "c1")
    except:
        c1 = 0.2
    try:
        b1 = options.get_double(option_section, "d1")
    except:
        b1 = 0.32
    try:
        d1 = options.get_double(option_section, "d1")
    except:
        d1 = 0.46
    




    #__________________________________#
    # computation grid params
    # only things I would touch
    # are z_max and kh_max (increase only)
    # effect of increasing kh_max is 
    # small
    #__________________________________#
    
    try:
        z_max = options.get_double(option_section, "z_max")
    except:
        z_max = 3.
    try:
        z_min = options.get_double(options_section, "z_min")
    except:
        z_min = 3.5 * 10. ** (-6.)
    try:
        kh_min = options.get_double(option_section, "kh_min")
    except:
        kh_min = 10. ** (-3.)
    try:
        kh_max = options.get_double(option_section, "kh_max")
    except:
        kh_max = 10.






    #__________________________________#
    # limber on off. Even when limber
    # = False, Limber will still turn on 
    # at an l-mode of limber begin
    #__________________________________#


    try:
        limber = options.get_bool(option_section, "limber")
    except:
        limber = True
        print "CAUTION: You have not specified a value for limber"
    try:
        limber_begin = options.get_int(option_section, "limber_begin")
    except:
        limber_begin = 100
        print "CAUTION: You have not specified a value for limber_begin"






    #__________________________________#
    # tomography on off
    # + number of bins 
    # + binning strategy
    #__________________________________#

    try:
        tomographic = options.get_bool(option_section, "tomographic")
    except:
        tomographic = False
        print "CAUTION: You have not specified a value for tomographic"
    if tomographic == True:
        try: 
            equal_z_tomo_bins = options.get_bool(option_section, "equal_z_tomo_bins")
        except:
            equal_z_tomo_bins = False
            print "CAUTION: You are using equal number of gals per bin"
    else:
        equal_z_tomo_bins = False
    print equal_z_tomo_bins
    try:
        num_tomo_bins = options.get_int(option_section, "num_tomo_bins")
    except:
        num_tomo_bins = 10
        print "CAUTION: You have not specified a value for num_tomo_bins"
    





    

    #__________________________________#
    # survey specific parameters
    # 
    #__________________________________#

    try:
        rms_ellipticity = options.get_double(option_section, "rms_ellipticity")
    except:
        rms_ellipticity = 0.3
        print "CAUTION: You have not specified a value for rms_ellipticity"
    try:
        number_of_galaxies = options.get_double(option_section, "number_of_galaxies")
    except:
        number_of_galaxies = 1.6e9
        print "CAUTION: You have not specified a value for number_of_galaxies"
    





    #__________________________________#
    # specify where you want to get
    # the power spectrum from.
    # To read from a custom file,
    # just specify a directory. You
    # will need files: z.txt, a.txt, 
    # d_m.tx, p_k.txt, k_h.txt
    # These should be in the same format
    # as in the cosmosis demo1 output
    #__________________________________#


    try:
        linear = options.get_bool(option_section, "linear")
    except:
        linear = True
        print "CAUTION: You are using the linear power spectrum"
    try:
        read_from_file = options.get_bool(option_section, "read_from_file")
    except:
        read_from_file = False
        print "Power Spectra will be read from the pipeline."
        print "To read from file set read_from_file = True"
    if read_from_file == True:
        try:
            file_directory = options.get(option_section, "file_directory")
            print "You are reading data from files in:"
            print file_directory
        except:
            file_directory = 0
            print "File directory has not been specified"
    else:
        file_directory = 0
    try:
        use_default_cosmology = options.get_bool(option_section, "use_default_cosmology")
    except:
        use_default_cosmology = False
    if use_default_cosmology == True:
        print "CAUTION: Using default cosmology"
    else:
        pass
    






    #__________________________________#
    # load the weights
    # if you are doing weighted 
    # lensing
    #__________________________________#

    if weighted_lensing == False:
        weight_list = 0.
    else:
        print "you are loading the weights"
        weight_list = [ i for i in range((l_max + 1))]
        for li in l_input_array:
            weight_list[li] = np.loadtxt('weights/weight_%s.txt' %li)






    #__________________________________#
    # pre-calculate the bessel data
    # just once. don't calculate
    # the bessel data at all
    # if we are using limber
    # without spherical bessel 
    # lensing
    #__________________________________#

    if limber == False:
        bessel_data, arguments = lensing_calculation.compute_bessel_data(resolution, resolution_sub_samps, resolution_bes, z_min, z_max, kh_min, kh_max, resolution_bes_lim, tolerance, l_max, limber_begin)
    else:
        if tomographic == True or weighted_lensing == True:
            bessel_data, arguments = np.zeros((resolution_bes, l_max + 1)), np.zeros(resolution_bes)
        else:
            bessel_data, arguments = lensing_calculation.compute_bessel_data(resolution, resolution_sub_samps, resolution_bes, z_min, z_max, kh_min, kh_max, resolution_bes_lim, tolerance, l_max, limber_begin)






    #__________________________________#
    # load default files to calculate
    # the exterior bessel funcitons
    # using a look up table 
    #__________________________________#

    
    cwd = os.getcwd()
    z_default_file = np.loadtxt(cwd + "/z.txt")
    d_m_default_file = np.loadtxt(cwd + "/d_m.txt")





    #__________________________________#
    # return all the variables
    # to pass them to cosmosis
    # exectute
    #__________________________________#



    return z_default_file, d_m_default_file, curved_lens_kernel, omega_k, hyper_spherical,  weighted_lensing, weight_list, l_input_array,  resolution, resolution_sub_samps, resolution_bes, resolution_bes_lim, tolerance, l_max, c_cal, z_bias, A_sigma, a1, c1, b1, d1, z_max, z_min, kh_min, kh_max, limber,  tomographic, equal_z_tomo_bins, num_tomo_bins, linear, limber_begin, bessel_data, arguments,  rms_ellipticity, number_of_galaxies,  read_from_file, file_directory, use_default_cosmology









#__________________________________#
# executed for each new set of
# cosmological parameters
# 
#__________________________________#


def execute(block, config):


    #__________________________________#
    # get the variable from config
    #__________________________________#
 

    z_default_file, d_m_default_file, curved_lens_kernel, omega_k, hyper_spherical,weighted_lensing, weight_list, l_input_array,  resolution, resolution_sub_samps, resolution_bes, resolution_bes_lim, tolerance, l_max, c_cal, z_bias, A_sigma, a1, c1, b1, d1, z_max, z_min, kh_min, kh_max, limber,  tomographic, equal_z_tomo_bins, num_tomo_bins, linear, limber_begin, bessel_data, arguments,  rms_ellipticity, number_of_galaxies,  read_from_file, file_directory, use_default_cosmology = config
    



    #__________________________________#
    # physical constants
    #__________________________________#
    omega_m = block["cosmological_parameters","omega_m"]
    hubble_constant = 100. * block["cosmological_parameters","h0"]
    #km/s
    light_speed = 299792.




    #__________________________________#
    # get and format power spectrum
    # and other cosmological arrays
    # for custom input arrays
    # use the format provided with
    # fiducial_cosmology subdirectory
    #__________________________________#    
    if use_default_cosmology == False:
        if read_from_file == False:
            if linear == True:
                a = block["distances", "a"]
                d_m = block["distances", "d_m"]
                k_h = block["matter_power_lin", "k_h"]
                p_k = block["matter_power_lin", "p_k"]
                z = block["matter_power_lin", "z"]
                # new version of cosmosis
                # arrays are reversed in distance
                # section of data block. This case division
                # handles both versions of Cosmosis
                if np.all(np.diff(d_m) > 0) == True and np.all(np.diff(a) > 0) == False:
                    d_m, a = d_m[::-1], a[::-1]
                else:
                    pass
            elif linear == False:
                a = block["distances", "a"]
                d_m = block["distances", "d_m"]
                k_h = block["matter_power_nl", "k_h"]
                p_k = block["matter_power_nl", "p_k"]
                z = block["matter_power_nl", "z"]
                # new version of cosmosis
                # arrays are reversed in distance
                # section of data block. This case division
                # handles both versions of cosmosis
                if np.all(np.diff(d_m) > 0) == True and np.all(np.diff(a) > 0) == False:
                    d_m, a = d_m[::-1], a[::-1]
                else:
                    pass
            else:
                print "linear must be either true or false"
        else:
            a = np.loadtxt(file_directory + "/a.txt")
            d_m = np.loadtxt(file_directory + "/d_m.txt")
            k_h = np.loadtxt(file_directory +"/k_h.txt")
            p_k = np.loadtxt(file_directory + "/p_k.txt")
            z = np.loadtxt(file_directory + "/z.txt")
            z = z[::-1]
    else:
        file_directory = "fiducial_cosmology"
        a = np.loadtxt(file_directory + "/a.txt")
        d_m = np.loadtxt(file_directory + "/d_m.txt")
        k_h = np.loadtxt(file_directory +"/k_h.txt")
        p_k = np.loadtxt(file_directory + "/p_k.txt")
        z = np.loadtxt(file_directory + "/z.txt")
        z = z[::-1]
    p_k = p_k[::-1,:]
    d_m, a, p_k, z  = d_m[::-1], a[::-1], p_k[::-1,:], z[::-1]







    #__________________________________#
    # calculate the Cls
    #__________________________________#

    
    c_l_storage_array, shot_noise_storage_array, z_sub_samps, k_sub_samps, p_k_samps, z_samps, k_samps, r_samps, a_samps = lensing_calculation.main(z_default_file, d_m_default_file, curved_lens_kernel, omega_k, hyper_spherical,weighted_lensing, weight_list, l_input_array,  resolution, resolution_sub_samps, resolution_bes, resolution_bes_lim, tolerance, l_max, c_cal, z_bias, A_sigma, a1, c1, b1, d1, z_max, z_min, kh_min, kh_max, limber,  tomographic, equal_z_tomo_bins, num_tomo_bins, linear, limber_begin, bessel_data, arguments,  rms_ellipticity, number_of_galaxies,  read_from_file, file_directory, use_default_cosmology, omega_m, hubble_constant, light_speed, d_m, k_h,  a, p_k, z)

        


    #__________________________________#
    # save the output to the datablock
    #__________________________________#

 

    block["3d_weak_lensing_output", "c_l_storage_array"]= c_l_storage_array
    block["3d_weak_lensing_output", "shot_noise_storage_array"]= shot_noise_storage_array
    block["3d_weak_lensing_output", "z_samps"]= z_samps
    block["3d_weak_lensing_output", "r_samps"]= r_samps
    # Caution these are not in units of h
    block["3d_weak_lensing_output", "k_samps"]= k_samps
    block["3d_weak_lensing_output", "p_k_samps"]= p_k_samps
    block["3d_weak_lensing_output", "l_input_array"]= l_input_array


    
    return 0






#__________________________________#
# no clean up
# to be done
#__________________________________#

def cleanup(config):
    pass



#__________________________________#
#      GLaSS                       #
#      Peter Taylor                #
#      Mullard Space Laboratory    #
#      University College London   #
#      2018                        #
#__________________________________#



import numpy as np
import scipy
from ctypesGsl import sf_bessel_jl_array
from scipy.special import spherical_jn
import math
from scipy import integrate
from scipy.integrate import dblquad
import os
import scipy.interpolate
import itertools
from scipy import interpolate






###################_______________________################
#
#
#      Section 1: Pre-compute Bessel Data in l and x  
#
#
#
###################_______________________################





#__________________________________#
# get extent of k and r 
# these are used to put
# limits on x when calculating
# j_l(x)
#_________________________________#


def create_fiducial_grid_for_bessel(resolution, resolution_sub_samps, z_min, z_max, kh_min, kh_max):
    cwd = os.getcwd()
    z = np.loadtxt(cwd + "/z.txt")
    d_m = np.loadtxt(cwd + "/d_m.txt")
    #reverse the arrays to do interpolation
    z, d_m = z[::-1], d_m[::-1]
    z_samps = np.linspace(z_min, z_max, num = resolution)
    r_samps   = np.interp(z_samps, z, d_m)
    k_samps = np.logspace(np.log10(kh_min), np.log10(kh_max), resolution)
    h_o = 0.67
    k_samps = k_samps * h_o
    return  r_samps, k_samps



#__________________________________#
# use millers algorithem
# so we need to calculate
# max l for each x to 
# avoid underflow
#_________________________________#


def find_max_l(resolution_bes, resolution_bes_lim, resolution, r_samps, k_samps, tolerance, l_max, limber_begin):
    x_min = np.min(k_samps) * np.min(r_samps)
    x_max = np.max(k_samps) * np.max(r_samps)
    x = np.logspace(np.log10(x_min), np.log10(x_max), num = resolution_bes_lim)
    number_x_samples = np.shape(x)[0]
    l_end = np.zeros(number_x_samples)
    for i in range(number_x_samples):
        if abs(spherical_jn(l_max, x[i])) > tolerance:
            l_end[i] = l_max
        else:
            li = 0
            while li < l_max:
                if abs(spherical_jn(li, x[i])) < tolerance and abs(spherical_jn(li+1, x[i])) < tolerance:
                    l_end[i] = li
                    li += l_max
                elif li == l_max-1  and abs(spherical_jn(li+1, x[i])) < tolerance:
                    l_end[i] = li
                    li += l_max
                else:
                    li += 1
    arguments = np.logspace(math.log10(x_min), math.log10(x_max), num = resolution_bes)
    l_end = np.interp(arguments, x, l_end)
    return l_end.astype(np.int64), arguments





#__________________________________#
# compute bessel data
# usin ultra fast ctypes GSL
# implementation 
# and store in an array
#  in (l,x)
# 
#_________________________________#



def create_bessel_data(l_end, arguments, l_max, resolution):
    bessel_data = np.zeros((np.shape(arguments)[0], l_max + 1))
    zero_array = np.zeros((l_max + 1))
    for i in range(np.shape(arguments)[0]):
        fixed_x_bessel_sub_array = sf_bessel_jl_array(l_end[i], arguments[i])
        fixed_x_bessel_array = zero_array.copy()
        fixed_x_bessel_array[:len(fixed_x_bessel_sub_array)] += fixed_x_bessel_sub_array
        bessel_data[i][:] = fixed_x_bessel_array
    return bessel_data




#__________________________________#
# this is the function  
# that is called by GLaSS.py to
#  produce all the bessel data.
# calls all the above functions in 
# this section in order.
#_________________________________#



def compute_bessel_data(resolution, resolution_sub_samps, resolution_bes, z_min, z_max, kh_min, kh_max, resolution_bes_lim, tolerance, l_max, limber_begin):
    r_samps, k_samps = create_fiducial_grid_for_bessel(resolution, resolution_sub_samps, z_min, z_max, kh_min, kh_max) 
    l_end, arguments = find_max_l(resolution_bes, resolution_bes_lim, resolution, r_samps, k_samps, tolerance, l_max, limber_begin)
    bessel_data = create_bessel_data(l_end, arguments, l_max, resolution)
    return bessel_data, arguments






















































###################_______________________################
#
#
#      Section 2: Functions needed to compute C_l  
#
#
#
###################_______________________################





#__________________________________#
# lay down a computation grid
#  in z and k.
# sub samps is all legacy
# so z_samps and z_sub_samps
# are now interchangeable
#_________________________________#

def create_grid(resolution, resolution_sub_samps, z_min, z_max, kh_min, kh_max, z, d_m, hubble_constant) :
    #reverse the z array, but not the d_m array. Legacy from old cosmosis
    z = z[::-1]
    z_samps = np.linspace(z_min, z_max, num = resolution)
    r_samps   = np.interp(z_samps, z, d_m)
    k_samps = np.logspace(np.log10(kh_min), np.log10(kh_max), resolution)
    h_o = hubble_constant / 100.
    k_samps = k_samps * h_o
    spacing = int( resolution / resolution_sub_samps)
    z_sub_samps = z_samps[::spacing]
    k_sub_samps = k_samps[::spacing]
    return z_samps, r_samps, k_samps, spacing, z_sub_samps, k_sub_samps






#__________________________________#
# For the exterior Bessel Function
# in 3D cosmic shear
# the mapping from z to r should
# not change with cosmology
# so we compute default r values
# on grid
#_________________________________#

def create_grid_exterior(z_samps, z_default_file, d_m_default_file):
    z, d_m = z_default_file[::-1], d_m_default_file[::-1]
    r_samps_exterior   = np.interp(z_samps, z, d_m)
    return r_samps_exterior






#__________________________________#
# Interpolate a and p_k onto
# the calculation grid
#_________________________________#

def interpolation(a, d_m, k_h, p_k, z, z_samps, r_samps, k_samps,  resolution, hubble_constant):
    h_o = hubble_constant / 100.
    k_h = k_h * h_o
    z = z[::-1]
    f = interpolate.RectBivariateSpline(z, k_h, p_k)
    p_k_samps = f(z_samps, k_samps)
    a_samps = np.interp(r_samps, d_m, a)
    return p_k_samps, a_samps




#__________________________________#
# create look up table that 
# maps bessel data stored in x
# to k and r
#_________________________________#


def create_bessel_look_up_table(arguments, resolution, r_samps, k_samps):
    grid = r_samps[:,None] * k_samps[None,:]
    grid = grid.ravel()
    look_up_table = find_closest(arguments, grid)
    look_up_table = look_up_table.reshape((resolution,resolution))
    return look_up_table.astype(np.int64)



#__________________________________#
# same as the above but with 
# the exterior bessel function
# which does not change with 
# cosmology
#_________________________________#

def create_bessel_look_up_table_exterior(arguments, resolution, r_samps_exterior, k_samps):
    grid = r_samps_exterior[:,None] * k_samps[None,:]
    grid = grid.ravel()
    look_up_table = find_closest(arguments, grid)
    look_up_table = look_up_table.reshape((resolution,resolution))
    return look_up_table.astype(np.int64)



#__________________________________#
# A function that is called
#  in the two functions above
# to make the lookup table
#_________________________________#


def find_closest(A, target):
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = idx-1
    right = idx
    idx -= target - left < right - target
    return idx




#__________________________________#
# Guassian photometric uncertainty
#________________________________#

def photometric_error(z_samps, c_cal, z_bias, A_sigma, z_min, z_max):
    sigma_z = A_sigma * ( 1 + z_samps)
    prob_z_prime_z_p_matrix = 1. / ( 2. * math.pi) ** 0.5 / sigma_z[:,None] * np.exp( - (z_samps[None,:] - c_cal * z_samps[:,None] + z_bias) ** 2. / (2. * sigma_z[:,None] ** 2.) )
    return prob_z_prime_z_p_matrix







#__________________________________#
# currently using analytic n(z)
# but there are plans to use 
# numerical n(z) shortly.
# In the meantime, if you want to
# use custom n(z) just do a 
# find for variable a1 and
# replace default n(z) everywhere 
# where that variable appears in
# this file.
#________________________________#

def redshift_density_function(z_samps, a1, c1, b1, d1, resolution, z_min, z_max, num_tomo_bins, tomographic, equal_z_tomo_bins):
    if tomographic == False:
        redshift_density_pdf = a1 / c1 * np.exp( - (z_samps - 0.7) ** 2. / (b1 ** 2.) ) + np.exp( - (z_samps - 1.2) ** 2. / (d1 ** 2.) )
        #not used in the non-tomographic case
        break_indices = np.array([0])
        noise_diagonal = np.array([0])
    elif tomographic == True:
        redshift_density_pdf = a1 / c1 * np.exp( - (z_samps - 0.7) ** 2. / (b1 ** 2.) ) + np.exp( - (z_samps - 1.2) ** 2. / (d1 ** 2.) )
        #needed for tomo binning
        pdf = lambda z: a1 / c1 * np.exp( - (z - 0.7) ** 2. / (b1 ** 2.) ) + np.exp( - (z - 1.2) ** 2. / (d1 ** 2.) )
        norm = integrate.quad(pdf, z_min, z_max)[0]
        norm_pdf = lambda z: (a1 / c1 * np.exp( - (z - 0.7) ** 2. / (b1 ** 2.) ) + np.exp( - (z - 1.2) ** 2. / (d1 ** 2.) )) / norm
        integral = lambda z: integrate.quad(norm_pdf, z_min, z)[0]
        #create the hat functions( see kithing limits of cosmic shear)
        if equal_z_tomo_bins == False:
            test_indices = np.zeros(1920)
            test_integrals = np.zeros(1920)
            for i in range(1920):
                test_indices[i] = np.abs(z_samps - 0.0015625 * i).argmin()
                test_indices = test_indices.astype(np.int64)
                test_integrals[i] = integral(z_samps[test_indices[i]])
            break_indices = np.zeros(num_tomo_bins)
            for i in range(num_tomo_bins):
                value = ((i+1.)/(num_tomo_bins))
                break_indices[i] = test_indices[np.abs(test_integrals - value).argmin()]
            #break_indices = np.append(break_indices,resolution-1)
            break_indices = np.insert(break_indices,0,0)
        else:
            break_indices = np.linspace(0, resolution - 1, num_tomo_bins + 1)
        #print break_indices
        noise_diagonal = np.zeros(num_tomo_bins)
        for i in range(np.shape(break_indices)[0] - 1):
            #print [break_indices[i]], [break_indices[i+1]]
            noise_diagonal[i] = ((integrate.quad(norm_pdf, z_samps[int(break_indices[i])], z_samps[int(break_indices[i+1])])[0]))
            #print  0.3 ** 2. / 2. / 3.14 ** 2.
    else:
        print "tomographic variable must be either True or False"
    return redshift_density_pdf, break_indices.astype(np.int64), noise_diagonal








#__________________________________#
# normalises analytic n(z)p(z|z_p)
# in the signal.
#________________________________#


def normalise_signal(z_min, z_max, A_sigma, c_cal, z_bias, a1, c1, b1, d1):
    signal_norm = dblquad(lambda z_p, z_prime: 1. / ( 2. * math.pi) ** 0.5 / (A_sigma * ( 1. + z_p)) * np.exp( - (z_prime - c_cal * z_p + z_bias) ** 2. / (2. * (A_sigma * ( 1 + z_p)) ** 2.) ) * (a1 / c1 * np.exp( - (z_p - 0.7) ** 2. / (b1 ** 2.) ) + np.exp( - (z_p - 1.2) ** 2. / (d1 ** 2.) )), z_min, z_max, lambda z_prime: z_min, lambda z_prime: z_max)[0]
    return signal_norm






#__________________________________#
# normalises n(z) in the noise.
# 
#________________________________#


def normalise_noise(z_min, z_max, a1, c1, b1, d1):
    pdf = lambda z: a1 / c1 * np.exp( - (z - 0.7) ** 2. / (b1 ** 2.) ) + np.exp( - (z - 1.2) ** 2. / (d1 ** 2.) )
    noise_norm = integrate.quad(pdf, z_min, z_max)[0]
    return noise_norm






#__________________________________#
#  Create the top hats
#  for tomographic lensing
#_________________________________#


def create_top_hats_tomo(resolution, spacing, break_indices, tomographic, equal_z_tomo_bins, num_tomo_bins):
    if tomographic == True:
        top_hat_dict = {}
        for i in range(num_tomo_bins):
            top_hat_array = np.zeros(resolution)
            top_hat_array[break_indices[i]:break_indices[i+1]] = 1
            top_hat_2d = top_hat_array[:,None] * top_hat_array[None,:]
            top_hat_dict[i] = top_hat_2d
        top_hat = sum(top_hat_dict.itervalues())
    elif tomographic == False:
        top_hat = None
    else:
        "tomographic must be true of false."
    return top_hat








#__________________________________#
# Compute hyperspherical Bessel Functions
# For non-flat cosmologies only
# First need to install Class 
# and kindly ask Thomas Tram for
# the hyperspherical Bessel Wrappers
# There are lots of rescalings in this
# function by c and H etc to be consistent
# with Tram paper
# computation of Hyperspherical Bessel Functions
#_________________________________#


def compute_hyper_spherical(l_input_array, r_samps, k_samps, omega_k, hubble_constant, light_speed, l_max):
    from classy import Class
    cosmo = Class()
    hyper_spherical_k_r_l = np.zeros( ((np.shape(k_samps)[0]), np.shape(r_samps)[0], np.shape(l_input_array)[0]))
    Ksign = np.sign(- omega_k)
    r_samps = r_samps / light_speed
    ##renormalisation##
    #see http://cds.cern.ch/record/354727/files/9805173.pdf page 3# 
    beta_samps = np.sqrt(k_samps ** 2. * light_speed ** 2. + hubble_constant ** 2. * (-omega_k))
    print 'computing hyper-spherical bessel functions'
    if Ksign == -1:
        for k in range(np.shape(k_samps)[0]):
            hyper_spherical_k_r_l[k,:,:] = cosmo.hyperspherical_recurrence(l_input_array.max(), beta_samps[k], Ksign, r_samps)[:,l_input_array]
    elif Ksign == 1:
        for k in range(np.shape(k_samps)[0]):
            l_max_allowed = int(beta_samps[k] - 1) 
            if l_max_allowed < l_max:
                beta_samp = int(round(beta_samps[k]))
                hyper_spherical_k_r_l[k,:,:] = np.hstack((cosmo.hyperspherical_recurrence(l_max_allowed, beta_samp, Ksign, r_samps), np.zeros(( np.shape(r_samps)[0], l_input_array.max() - l_max_allowed + 1) )))[:,l_input_array]
            else:
                hyper_spherical_k_r_l[k,:,:] = cosmo.hyperspherical_recurrence(l_input_array.max(), beta_samps[k], Ksign, r_samps)[:,l_input_array]
    else:
        print "error K should either be 1 or -1"
    return hyper_spherical_k_r_l




#__________________________________#
#  Use the look up table to 
#  get j_l(kr) from j_l(x)
#  for some fixed x
#_________________________________#


def create_bessel_grid(bessel_data, look_up_table, li):
    bessel_grid = bessel_data[look_up_table, li]
    return bessel_grid




#__________________________________#
#  same as above but just for
#  the exterior Bessel funciton
#  in spherical-Bessel lensing
#  this does not change with 
#_________________________________#


def create_bessel_grid_exterior(bessel_data, look_up_table_exterior, li):
    bessel_grid = bessel_data[look_up_table_exterior, li]
    return bessel_grid






#__________________________________#
#  calculate the U-matrix
#_________________________________#
def create_u_matrix(bessel_grid, a_samps, p_k_samps, r_samps, k_samps, resolution, limber, li, limber_begin, curved_lens_kernel, omega_k, light_speed, hubble_constant, hyper_spherical):
    if limber == False:
        if li < limber_begin:
            u_matrix = np.zeros((resolution, resolution))
            lensing_kernel = lensing_kernel_func(r_samps[:,None], r_samps[None,:], curved_lens_kernel, omega_k, light_speed, hubble_constant)
            lensing_kernel[lensing_kernel < 0.] = 0.
            deltas = r_samps - np.append(r_samps[0], r_samps[:-1])
            r_r_prime_matrix = np.array([deltas]*resolution) * lensing_kernel / np.array([a_samps]*resolution)
            r_prime_k_matrix = bessel_grid * (p_k_samps ** 0.5 )
            u_matrix = np.dot(r_r_prime_matrix, r_prime_k_matrix)
        else:
            nu_samps = (li + 0.5)/ k_samps
            nearest_r_samps = find_closest(r_samps, nu_samps)
            lensing_kernel = lensing_kernel_func(r_samps[:,None], nu_samps[None,:], curved_lens_kernel, omega_k, light_speed, hubble_constant)
            if curved_lens_kernel == True:
                if omega_k >= 0.:
                    pass
                # curvature only: K = -1 numerical errors need to be cleaned up
                else:
                    lensing_kernel_test = lensing_kernel_func(r_samps[:,None], nu_samps[None,:], curved_lens_kernel, 0, light_speed, hubble_constant)
                    lensing_kernel[list(np.where(lensing_kernel_test < 0.))] = 0.
            else:
                pass
            lensing_kernel[lensing_kernel < 0.] = 0.
            denominator = k_samps[None,:] * a_samps[nearest_r_samps]
            numerical_factor = math.sqrt(math.pi / (2 * (li + 0.5)))
            k_samples = np.linspace(0,resolution-1, resolution).astype(np.int64)
            power_spec = p_k_samps[nearest_r_samps,k_samples] ** 0.5
            u_matrix = lensing_kernel * numerical_factor * power_spec[None,:] / denominator
            if hyper_spherical == False:
                pass
            else:
                beta_samps = np.sqrt(k_samps ** 2. * light_speed ** 2. + hubble_constant ** 2. * (-omega_k))
                k_hat = -np.sign(omega_k)
                prefactor = np.zeros((resolution, resolution))
                for ri in range(resolution):
                    for ki in range(resolution):
                        #only valid when l < \beta for minus case U should be exactly 0 
                        #in this case but due to numerics this needs to be inforced#
                        if li < beta_samps[ki]:
                            prefactor[ri,ki] = (1.- k_hat * li ** 2. / beta_samps[ki] ** 2.) ** (-0.25)
                        else:
                            pass 
                u_matrix = prefactor * u_matrix 
    elif limber == True:
        nu_samps = (li + 0.5)/ k_samps
        nearest_r_samps = find_closest(r_samps, nu_samps)
        lensing_kernel = lensing_kernel_func(r_samps[:,None], nu_samps[None,:], curved_lens_kernel, omega_k, light_speed, hubble_constant)
        lensing_kernel[lensing_kernel < 0.] = 0.
        denominator = k_samps[None,:] * a_samps[nearest_r_samps]
        numerical_factor = math.sqrt(math.pi / (2 * (li + 0.5)))
        k_samples = np.linspace(0,resolution-1, resolution).astype(np.int64)
        power_spec = p_k_samps[nearest_r_samps,k_samples] ** 0.5
        u_matrix = lensing_kernel * numerical_factor * power_spec[None,:] / denominator
        if hyper_spherical == False:
            pass
        else:
            beta_samps = np.sqrt(k_samps ** 2. * light_speed ** 2. + hubble_constant ** 2. * (-omega_k))
            k_hat = -np.sign(omega_k)
            prefactor = np.zeros((resolution, resolution))
            for ri in range(resolution):
                for ki in range(resolution):
                        if li < beta_samps[ki]:
                            prefactor[ri,ki] = (1.- k_hat * li ** 2. / beta_samps[ki] ** 2.) ** (-0.25)
                        else:
                            pass 
            u_matrix = prefactor * u_matrix 
    else:
        "print mode limber must be true or false"
    return u_matrix





#__________________________________#
#  get lensing kernel for 
#  U-matrix
#_________________________________#


def lensing_kernel_func(r, r_prime, curved_lens_kernel, omega_k, light_speed, hubble_constant):
    if curved_lens_kernel == False:
        return (r - r_prime) /r / r_prime
    else:
        if omega_k == 0.:
            return (r - r_prime) /r / r_prime
        elif omega_k > 0.:
            K = - (hubble_constant / light_speed) ** 2. * omega_k
            return co_moving_angular_dist(K, (r - r_prime)) / co_moving_angular_dist(K, r) / co_moving_angular_dist(K, r_prime)
        else:
            K = - (hubble_constant / light_speed) ** 2. * omega_k
            return co_moving_angular_dist(K, (r - r_prime)) / co_moving_angular_dist(K, r) / co_moving_angular_dist(K, r_prime)




#__________________________________#
# get the correct co-moving dist.
# for lensing kernel when run
# in non-flat cosmo mode
#_________________________________#

def co_moving_angular_dist(K, r):
    if K > 0.:
        return K ** (-0.5) * np.sin(K ** 0.5 * r)
    else: 
        return (-K) ** (-0.5) * np.sinh((-K) ** 0.5 * r)





#__________________________________#
# Integrate over z_prime
# intermediate matrix calculation
# just like the G-matrix
# call the ouput of this the H-matrix
#_________________________________#
def integrate_z_prime(u_matrix, z_samps, prob_z_prime_z_p_matrix, l_max, resolution):
    deltas = z_samps - np.append(z_samps[0], z_samps[:-1])
    z_p_z_prime_matrix = np.array([deltas]*resolution) * prob_z_prime_z_p_matrix.transpose()
    H_matrix = np.zeros((resolution, resolution))
    H_matrix = np.dot( z_p_z_prime_matrix, u_matrix)
    return H_matrix






#__________________________________#
# Calculate the G_matrix
#__________________________________#

def create_G_matrix(bessel_grid_exterior, spacing, redshift_density_pdf, H_matrix, z_samps, k_samps, resolution, tomographic, equal_z_tomo_bins, top_hat, resolution_sub_samps, weighted_lensing, weight_list,  li,  curved_lens_kernel, omega_k, hyper_spherical, r_samps):
    if tomographic == False:
        deltas = z_samps - np.append(z_samps[0], z_samps[:-1])
        z_p_array = deltas * redshift_density_pdf
        G_matrix = np.zeros((resolution_sub_samps, resolution))
        if weighted_lensing == False:
            bessel_grid_sub_samps = bessel_grid_exterior[:,::spacing]
            k1_z_p_matrix = np.array([z_p_array]* resolution_sub_samps) * bessel_grid_sub_samps.transpose()
            G_matrix = np.dot(k1_z_p_matrix, H_matrix)
        else:
            weight = weight_list[li]
            k1_z_p_matrix = np.array([z_p_array]* resolution_sub_samps) * weight.transpose()
            G_matrix = np.dot(k1_z_p_matrix, H_matrix)
    elif tomographic == True:
        deltas = z_samps - np.append(z_samps[0], z_samps[:-1])
        z_p_array = deltas * redshift_density_pdf
        G_matrix = np.zeros((resolution, resolution))
        z_z_p_matrix = np.array([z_p_array]*resolution) * top_hat
        G_matrix = np.dot(z_z_p_matrix, H_matrix)
    else:
        print "Error: tomographic must be either true or false" 
    return G_matrix








#__________________________________#
# Calculate the C_l
#__________________________________#
def c_l_function(G_matrix, resolution_sub_samps, k_samps,  resolution, li, tomographic, equal_z_tomo_bins, break_indices, num_tomo_bins, r_samps, hubble_constant):
    if tomographic == False:
        deltas = k_samps - np.append(k_samps[0], k_samps[:-1])
        c_l = np.zeros((resolution_sub_samps, resolution_sub_samps))
        measure = np.array([deltas / k_samps ** 2.]*resolution_sub_samps)
        k1_k_prime_matrix = measure * G_matrix
        c_l =  np.dot(k1_k_prime_matrix, G_matrix.transpose())
    elif tomographic == True:
        deltas = k_samps - np.append(k_samps[0], k_samps[:-1])
        c_l_duplicate_rows = np.zeros((resolution, resolution))
        measure = np.array([deltas / k_samps ** 2]*resolution)
        z1_k_prime_matrix = measure * G_matrix
        c_l_duplicate_rows = np.dot(z1_k_prime_matrix, G_matrix.transpose())
        c_l = np.zeros((num_tomo_bins, num_tomo_bins))
        #get rid of duplicate rows or columns
        for i in range(num_tomo_bins):
            for j in range(num_tomo_bins):
                c_l[i,j] = c_l_duplicate_rows[break_indices[i],break_indices[j]]
    else:
        print "Error: tomographic must be either true or false"
    return c_l










#__________________________________#
# compute the shot noise
#__________________________________#
def shot_noise_fnc(resolution, rms_ellipticity, bessel_grid_exterior, redshift_density_pdf, z_samps, k_samps, tomographic, equal_z_tomo_bins, resolution_sub_samps, spacing, num_tomo_bins, break_indices, r_samps, hubble_constant, li, weighted_lensing, weight_list, noise_diagonal):
    if tomographic == False:
        amplitude_shot_noise = rms_ellipticity  ** 2. / 2. / math.pi ** 2.
        deltas = z_samps - np.append(z_samps[0], z_samps[:-1])
        if weighted_lensing == False:
            k_z_matrix = np.array([deltas * redshift_density_pdf] * resolution_sub_samps) * bessel_grid_exterior[:,::spacing].transpose()
            z_k_prime_matrix = bessel_grid_exterior[:,::spacing]
        else:
            weight = weight_list[li]
            k_z_matrix = np.array([deltas * redshift_density_pdf] * resolution_sub_samps)  *  weight[::spacing,:] 
            z_k_prime_matrix = weight[::spacing,:].transpose()
        shot_noise_l =  amplitude_shot_noise * np.dot(k_z_matrix, z_k_prime_matrix)
    else:
        amplitude_shot_noise = rms_ellipticity ** 2. / 2. / math.pi ** 2.
        diagonal = noise_diagonal * amplitude_shot_noise
        shot_noise_l = np.diag(diagonal)
    return shot_noise_l

































































###################_______________________################
#
#
#      Section 3: Big Main Function Called 
#               by GLaSS.py
#       Runs all the Functions in Section 2
#       in correct order
#
###################_______________________################




def main(z_default_file, d_m_default_file, curved_lens_kernel, omega_k, hyper_spherical,weighted_lensing, weight_list, l_input_array,  resolution, resolution_sub_samps, resolution_bes, resolution_bes_lim, tolerance, l_max, c_cal, z_bias, A_sigma, a1, c1, b1, d1, z_max, z_min, kh_min, kh_max, limber,  tomographic, equal_z_tomo_bins, num_tomo_bins, linear, limber_begin, bessel_data, arguments,  rms_ellipticity, number_of_galaxies,  read_from_file, file_directory, use_default_cosmology, omega_m, hubble_constant, light_speed, d_m, k_h, a, p_k, z):


    #__________________________________#
    #     Do everything that we 
    #     don't want to loop over
    #_________________________________#


    h0 = hubble_constant / 100.
    z_samps, r_samps, k_samps, spacing, z_sub_samps, k_sub_samps = create_grid(resolution, resolution_sub_samps, z_min, z_max, kh_min, kh_max, z, d_m, hubble_constant)
    r_samps_exterior = create_grid_exterior(z_samps, z_default_file, d_m_default_file)
    p_k_samps, a_samps = interpolation(a, d_m, k_h, p_k, z, z_samps, r_samps, k_samps,  resolution, hubble_constant)
    look_up_table = create_bessel_look_up_table(arguments, resolution, r_samps, k_samps)
    look_up_table_exterior = create_bessel_look_up_table_exterior(arguments, resolution, r_samps_exterior, k_samps)
    prob_z_prime_z_p_matrix = photometric_error(z_samps, c_cal, z_bias, A_sigma, z_min, z_max)
    redshift_density_pdf, break_indices, noise_diagonal = redshift_density_function(z_samps, a1, c1, b1, d1, resolution,  z_min, z_max, num_tomo_bins, tomographic, equal_z_tomo_bins) 
    signal_norm = normalise_signal(z_min, z_max, A_sigma, c_cal, z_bias, a1, c1, b1, d1)
    noise_norm = normalise_noise(z_min, z_max, a1, c1, b1, d1)
    top_hat = create_top_hats_tomo(resolution, spacing, break_indices, tomographic, equal_z_tomo_bins, num_tomo_bins)
    if hyper_spherical == False:
        pass
    else:
        hyper_spherical_k_r_l = compute_hyper_spherical(l_input_array, r_samps, k_samps, omega_k, hubble_constant, light_speed, l_max)
    #create a big storage array that we dump all the c_l_data into
    if tomographic == False:
        c_l_storage_array = np.zeros(( resolution_sub_samps, resolution_sub_samps * np.shape(l_input_array)[0] ))
        shot_noise_storage_array = np.zeros(( resolution_sub_samps, resolution_sub_samps * np.shape(l_input_array)[0] ))
    else:
        c_l_storage_array = np.zeros(( num_tomo_bins, num_tomo_bins * np.shape(l_input_array)[0] ))
        shot_noise_storage_array = np.zeros(( num_tomo_bins, num_tomo_bins * np.shape(l_input_array)[0] ))
    

    #__________________________________#
    #       now loop over l-modes
    #_________________________________#
    for li in range(l_max + 1):
        if li in l_input_array:
            

            #__________________________________#
            #  calculate nested integrals
            #_________________________________#
            index = np.where(l_input_array == li)[0][0]
            bessel_grid = create_bessel_grid(bessel_data, look_up_table, li)
            bessel_grid_exterior = create_bessel_grid_exterior(bessel_data, look_up_table_exterior, li)
            if hyper_spherical == False:
                pass
            else:
                bessel_grid = hyper_spherical_k_r_l[:,:,index].T
            u_matrix = create_u_matrix(bessel_grid, a_samps, p_k_samps, r_samps, k_samps, resolution, limber, li, limber_begin, curved_lens_kernel, omega_k, light_speed, hubble_constant, hyper_spherical)
            H_matrix = integrate_z_prime(u_matrix, z_samps, prob_z_prime_z_p_matrix, l_max, resolution)
            G_matrix = create_G_matrix(bessel_grid_exterior, spacing, redshift_density_pdf, H_matrix, z_samps, k_samps, resolution, tomographic, equal_z_tomo_bins, top_hat, resolution_sub_samps, weighted_lensing, weight_list,  li,  curved_lens_kernel, omega_k, hyper_spherical, r_samps)


            #__________________________________#
            #  get all units correct etc
            #  for the pre-factors
            #_________________________________#
            amplitude = 9. * omega_m ** 2. * hubble_constant ** 4. / math.pi ** 4. / light_speed ** 4. * 1. / 16.
            amplitude = amplitude * 1. / h0 ** 3.  
            l_prefactor = (li + 2) * (li + 1) * li * (li - 1)



            #__________________________________#
            #  Finally get the Cls and 
            #  put all the ouput into a big
            #  2D array to pass through 
            #  cosmosis
            #_________________________________#
            if tomographic == False:
                    c_l =  1./ signal_norm ** 2. * number_of_galaxies**2. * l_prefactor * amplitude   * c_l_function(G_matrix, resolution_sub_samps, k_samps,  resolution, li, tomographic, equal_z_tomo_bins, break_indices, num_tomo_bins, r_samps, hubble_constant)
                    c_l_storage_array[:,index * resolution_sub_samps : (index + 1) * resolution_sub_samps] =  c_l
                    n_l = 1./ noise_norm * number_of_galaxies * shot_noise_fnc(resolution, rms_ellipticity, bessel_grid_exterior, redshift_density_pdf, z_samps, k_samps, tomographic, equal_z_tomo_bins, resolution_sub_samps, spacing, num_tomo_bins, break_indices, r_samps, hubble_constant, li, weighted_lensing, weight_list, noise_diagonal)
                    shot_noise_storage_array[:,index * resolution_sub_samps : (index + 1) * resolution_sub_samps] = n_l
            else:
                    c_l_storage_array[:,index * num_tomo_bins : (index + 1) * num_tomo_bins]=   1./ signal_norm ** 2. * number_of_galaxies **2. * l_prefactor * amplitude   * c_l_function(G_matrix, resolution_sub_samps, k_samps,  resolution, li, tomographic, equal_z_tomo_bins, break_indices, num_tomo_bins, r_samps, hubble_constant)
                    shot_noise_storage_array[:,index * num_tomo_bins : (index + 1) * num_tomo_bins] = number_of_galaxies * shot_noise_fnc(resolution, rms_ellipticity, bessel_grid_exterior, redshift_density_pdf, z_samps, k_samps, tomographic, equal_z_tomo_bins, resolution_sub_samps, spacing, num_tomo_bins, break_indices, r_samps, hubble_constant, li, weighted_lensing, weight_list, noise_diagonal)



    #__________________________________#
    #  save c_l data in a 2d array
    #_________________________________#
    return c_l_storage_array, shot_noise_storage_array, z_sub_samps, k_sub_samps, p_k_samps, z_samps, k_samps, r_samps, a_samps





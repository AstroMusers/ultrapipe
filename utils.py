from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.mast import Observations
import numpy as np
import pandas as pd
from dustmaps.bayestar import BayestarWebQuery
import math
import csv
from expecto import get_spectrum
from scipy.integrate import simpson

radial_pos = [0,0.42,0.84,1.27,1.69,2.11,2.53,2.95,3.38,3.8,4.22,4.64,5.06,5.49,5.91,6.33,6.75,7.17,7.6,8.02,8.44,8.86,9.28,9.71,10.13]
source_T = [9600,9040,8750,8310,7920,29200,23000,17600,15200,12300,11400,7350,7050,6700,6550,6300,6050,5800,5660,5440,5240,4960,4800,4600,4400,4000,3750,3700,3600,3500,3400,3200,3100,54000,37800]
def find_nearest(array, value):
    """A utility function that returns the index of the value in an array that is closest to the value provided.

    Args:
        array (np.array or similar): Array to find the nearest value.
        value (float or int): Value to find the closest value in the array

    Returns:
        index (int): The index of the element in the array that is closest to the provided value.
    """    ''''''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def TESS_MAG_to_Jy(TESS_Mag):
    """Converts from TESS magnitude to spectral flux density in the TESS band (600-1000 nm).

    Args:
        TESS_Mag (float): The TESS magnitude.

    Returns:
        Fv (float): The spectral flux density in Jy in the TESS band.
    """    ''''''
    return (2416)*10**(-0.4*TESS_Mag)

def Jy_to_AB_MAG(Jy):
    """Converts from the spectral flux density in a certain passband to the conventional AB magnitude.

    Args:
        Jy (float): Spectral flux density

    Returns:
        m_ab (float): Returns the AB magnitude based on the spectral flux density provided.
    """    ''''''
    return -2.5*np.log10(Jy) + 8.9

def SatMag(source_position,source_temp, SatMag_dir):
    """A function to be used with the 'SatMag.dat' file for ULTRASAT. Returns the saturation AB ULTRASAT magnitude
    of the provided source position (radially from the center of the detector) given a source temperature.

    Args:
        source_position (float): Radial position from the center of the detector.
        source_temp (float, int): Effective temperature of the source (K)

    Returns:
        Sat Mag (float): The saturation AB magnitude in the ULTRASAT band
    """    
    Saturation_Mag_Array = list(csv.reader(open(SatMag_dir)))
    Saturation_Mag_Array = [[float(y) for y in x] for x in Saturation_Mag_Array]
    position = find_nearest(radial_pos, source_position)
    temperature = find_nearest(source_T,source_temp)
    return Saturation_Mag_Array[temperature][position]

def LimMag(source_position,source_temp, LimMag_dir):
    """A function to be used with the 'LimMag.dat' file for ULTRASAT. Returns the limiting ULTRASAT AB magnitude
    of the provided source position (radially from the center of the detector) given a source temperature.

    Args:
        source_position (float): Radial position from the center of the detector.
        source_temp (float, int): Effective temperature of the source (K)

    Returns:
        Lim Mag (float): The limiting AB magnitude in the ULTRASAT band
    """   
    Limiting_Mag_Array = list(csv.reader(open(LimMag_dir)))
    Limiting_Mag_Array = [[float(y) for y in x] for x in Limiting_Mag_Array]
    position = find_nearest(radial_pos, source_position)
    temperature = find_nearest(source_T,source_temp)
    return Limiting_Mag_Array[temperature][position]

def cost_function(CVZ_RA, CVZ_DEC, CVZ_radius, targets_RAs, targets_DECs, targets_AB_mag, targets_temps, reddening=None,target_weight_list=None, data_weight_list=None, extinction_coeff = 5.487):
    """A function that provides a cost function on a particular CVZ FOV based on a target's RA, DEC, ULTRASAT AB mag, and effective temperature. Also accounts for dust extinction.

    Args:
        CVZ_RA (float): The RA where the proposed CVZ is centered.
        CVZ_DEC (float): The DEC where the proposed CVZ is centered.
        CVZ_radius (float): The radius of the proposed CVZ (assumed 7 degrees).
        targets_RAs (array of floats): The list of RAs of the proposed targets.
        targets_DECs (array of floats): The list of DECs of the proposed targets.
        targets_AB_mag (array of floats): The list of AB magnitude in the ULTRASAT band of the proposed targets.
        targets_temps (array of floats): The effective temperatures of the proposed targets
        reddening (float, optional): The dust extinction suffered by the target
        target_weight_list (array of floats, optional): A list corresponding to the weights given to a given target.
        data_weight_list (array of floats, optional): A list corresponding to the weights given to the associated data of a given target.
        extinction_coeff (float, optional): The coefficient relating exinction in the Pan-STARRS passband with a different passband.
        By default, the passband of Hubble's WFC3 F275W (~240-300 nm) is assumed, as well as an intersellar Rv value of 3.1.
        Refer to Table 6 in Shlafly & Finkbeiner 2011 https://iopscience.iop.org/article/10.1088/0004-637X/737/2/103#apj398709t6
        for alternative passbands
    Returns:
        Weight (float): Gives the aggregate weight of the proposed CVZ based on the passed targets.
    """    
    CVZ_coord = SkyCoord(CVZ_RA, CVZ_DEC, unit='deg', frame='icrs')
    CVZ_RA = float(CVZ_coord.ra.degree)
    CVZ_DEC = float(CVZ_coord.dec.degree)
    target_count = 0
    for index in range(len(targets_RAs)):
        target_RA = targets_RAs[index] - CVZ_RA
        target_DEC = targets_DECs[index] - CVZ_DEC
        if target_RA < 0:
            target_RA += 360
        if target_RA > 360:
            target_RA -= 360
        if CVZ_DEC == 90:
            if target_DEC > (90-CVZ_radius):
                dist = 0
        if CVZ_DEC == -90:
            if target_DEC < (-90+CVZ_radius):
                dist = 0
        dist = np.sqrt(target_RA**2 + target_DEC**2)
        if dist < CVZ_radius:
            sat_mag = SatMag(dist, targets_temps[index])
            lim_mag = LimMag(dist, targets_temps[index])
            
            red = reddening[index]
            if math.isnan(red) == True:
                red = 0
            target_weight = target_weight_list[index]
            data_weight = data_weight_list[index]
            AB_mag = targets_AB_mag[index]
            total_AB_Mag = AB_mag + extinction_coeff*red
            if total_AB_Mag > 18 and AB_mag < 18:
                print(AB_mag, total_AB_Mag)
            if (total_AB_Mag < 13) and (total_AB_Mag > sat_mag):
                target_count += target_weight*1*((4-1)*(data_weight - 1)/(2-1) + 1)
            elif (total_AB_Mag < 18) and (total_AB_Mag > sat_mag):
                target_count += target_weight*0.5*((4-1)*(data_weight - 1)/(2-1) + 1)
            elif (total_AB_Mag < 21) and (total_AB_Mag > sat_mag):
                target_count += target_weight*0.25*((4-1)*(data_weight - 1)/(2-1) + 1)
            # targets.append((target_name[index], targets_RAs[index], targets_DECs[index], targets_temps[index], total_AB_Mag, TESS_Mag[index], CVZ_Name))
    if target_count == 0:
        return -1
    return target_count

def cost_function_map(targets_RAs, targets_DECs, target_AB_MAG,target_distance, target_tmp, target_weight_list, data_weight_list, radius=7, output_dir=None, resolution=50, dust_map_version = 'bayestar2019'):
    """Produces a full-sky map of the cost function returned by the list of targets and the specified CVZ FOV. The function iterates a circular FOV with a provided radius
    iteratively across the sky. It then checks the lists of targets provided, and returns the weight described in the function 'cost_function'.
    It automatically corrects for magnitude increases due interstellar extinction in the ULTRASAT band.

    Args:
        targets_RAs (array of floats): An array of RAs of the targets.
        targets_DECs (array of floats): An array of DECs of the targets.
        target_AB_MAG (array of floats): An array of AB magntiudes in the ULTRASAT band of the targets.
        target_distance (array of floats): An array of target distances, in parsecs.
        target_tmp (array of floats): An array of target effective temperatures, in Kelvin.
        target_weight_list (array of floats): An array of assigned weights to each target.
        data_weight_list (array of floats): An array of assigned weights to each target for available data (can be passed as an array of 0s, if desired).
        radius (int, optional): The radius, in degrees, of the proposed FOV to iterate over the sky. Defaults to 7.
        output_dir (string, optional): The output directory of the cost function map. Defaults to None.
        resolution (int, optional): The N resolution of the resulting NxN map.  Defaults to 50.
        dust_map_version (str, optional): The dust map version to use to compute reddening. Options are bayestar2015, bayestar2017, and bayestar2019. Defaults to 'bayestar2019'.
    """    
    bayestar = BayestarWebQuery(version=dust_map_version)
    coords = SkyCoord(np.array(targets_RAs)*u.deg, np.array(targets_DECs)*u.deg, distance=np.array(target_distance)*u.pc,  frame='icrs')
    reddening = bayestar(coords, mode='median')    
    CF_map = []
    RA_span = np.linspace(0, 360, num=resolution)
    DEC_span = np.linspace(-90, 90, num=resolution)
    for i in range(len(DEC_span)):
        values_DEC = []
        for j in range(len(RA_span)):
            a = cost_function(RA_span[j], DEC_span[i], radius, targets_RAs, targets_DECs, target_AB_MAG, target_tmp, reddening, target_weight_list, data_weight_list)
            values_DEC.append(a)
        CF_map.append(values_DEC)
        print(i, a)
    if output_dir != None:
        np.savetxt(output_dir, CF_map, delimiter=',')

def data_coverage_map(output_dir, radius=7, resolution=50, project='TESS'):
    """Summary
    Generates a 2D Cartesian grid of available mission data based on the resolution and radius specified.
    Args:
        output_dir (string, _type_): The directory to save the data coverage map.
        radius (int, optional): The radius FOV (in degrees) to which to sum over the file counts.
        resolution (int, optional): The resolution of the 2D grid (e.g., 50 will output a 50x50 grid, evenly spaced in RA/DEC)
    """    
    data_map = []
    RA_span = np.linspace(0, 360, num=resolution)
    DEC_span = np.linspace(-90, 90, num=resolution)
    for i in range(len(DEC_span)):
        values_DEC = []
        for j in range(len(RA_span)):
            c = SkyCoord(ra=RA_span[j]*u.degree, dec=DEC_span[i]*u.degree, frame='icrs')
            count =  Observations.query_criteria_count(coordinates=c, radius = radius, project=project)
            values_DEC.append(count)
        data_map.append(values_DEC)
    if output_dir != None:
        np.savetxt(output_dir, data_map, delimiter=',')
        
def dust_map(resolution = 50, dist=50, extinction_coeff = 5.487, output_dir=None):
    """Produces an NxN dust map of interestellar extinction at a given distance.

    Args:
        resolution (int, optional): The size N of the NxN map produced. Defaults to 50.
        dist (int, optional): The distance assumed for the dust map, in parsecs. Defaults to 50.
        extinction_coeff (float, optional): The coefficient relating exinction in the Pan-STARRS passband with a different passband.
        By default, the passband of Hubble's WFC3 F275W (~240-300 nm) is assumed, as well as an intersellar Rv value of 3.1.
        Refer to Table 6 in Shlafly & Finkbeiner 2011 https://iopscience.iop.org/article/10.1088/0004-637X/737/2/103#apj398709t6
        for alternative passbands
        output_dir (string, optional): The file directory to which to save the dust map to. Defaults to None.
    """    
    RA_span = np.linspace(0, 360, num=resolution)
    DEC_span = np.linspace(-90, 90, num=resolution)
    tot_DEC = []
    tot_RA = []
    for DECs in DEC_span:
        for RAs in RA_span:
            tot_RA.append(RAs)
            tot_DEC.append(DECs)
            if DECs < 0:
                print(DECs)
    tot_DEC = np.array(tot_DEC)
    tot_RA = np.array(tot_RA)
    print(tot_DEC.shape, tot_RA.shape)
    dist_array = np.ones(len(tot_RA))*dist
    bayestar = BayestarWebQuery(version='bayestar2019')
    coords = SkyCoord(tot_RA*u.deg, tot_DEC*u.deg, distance=dist_array*u.pc,  frame='icrs')
    reddening = bayestar(coords, mode='median')      
    print(reddening)
    dust_map = []
    dust = []
    count = 0
    for extinction in reddening:
        if count == len(RA_span):
            dust = []
            dust_map.append(dust)
            count = 0
        if math.isnan(extinction) == True:
            extinction = 0
        dust.append(extinction*extinction_coeff)
        count += 1
    dust_map.append(dust)
    if output_dir != None:
        np.savetxt(output_dir, dust_map, delimiter=',') 

def TESS_ULTRASAT_color_transform(Teff, logg):
    """Returns the ratio of the integrated PHOENIX spectrum over the ULTRASAT and TESS bandpasses, respectively, of stars with the provided Teff and log(g).
    REQUIRES AN INTERNET CONNECTION.
    Args:
        Teff (array of float): Array of stellar effective temperatures, in Kelvin
        logg (array of float): Array of stellar log(g).
    Returns:
        Ratio (array of float): The ratio (ULTRASAT/TESS) of of the integrated stellar spectra.
    """    
    ratio = []    
    for index, value in enumerate(Teff):
        try: 
            spectrum = get_spectrum(
                T_eff=value, log_g=logg[index], cache=True
            )
            wavelength = np.array(spectrum.wavelength)
            flux = np.array(spectrum.flux)
            ULTRA_lower = find_nearest(wavelength, 2300)
            ULTRA_upper = find_nearest(wavelength, 2900)
            TESS_lower = find_nearest(wavelength, 6000)
            TESS_upper = find_nearest(wavelength, 10000)
            ratio.append(simpson(flux[ULTRA_lower:ULTRA_upper], wavelength[ULTRA_lower:ULTRA_upper])/simpson(flux[TESS_lower:TESS_upper], wavelength[TESS_lower:TESS_upper]))
        except:
            ratio.append(np.nan)
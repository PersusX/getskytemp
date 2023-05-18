#!/usr/bin/env python

# this sopftware is to get the galaxy sky temperature
# the data file is haslam408_dsds_Remazeilles2014.fits, Remazeilles et al. 2014
# Contact PersusXie@outlook.com if in doubt


import sys
import math
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm


def read_healpix_map(filename):
    map_data, header = hp.read_map(filename, h=True, verbose=True)
    nside = hp.get_nside(map_data)
    hdu = fits.open(filename)
    ordering = hdu[1].header['ORDERING']
    coord_sys = hdu[1].header['COORDSYS']
    hdu.close()
    
    return map_data, nside, ordering, coord_sys



def get_temp(map_data, nside, ordering, gl, gb, frequencies):
    npix = hp.nside2npix(nside)
    M_PI = math.pi
    for freq in frequencies:
        theta =  (-gb) * (np.pi / 180.0) + M_PI/2.
        phi = gl * (np.pi / 180.0) #* (np.pi / 180.0)
        i_val = hp.ang2pix(nside, theta, phi)
        scaled = (map_data[i_val] - 2.725) * (freq / 408.0)**(-2.6) + 2.725
        print(f"gl = {gl} gb = {gb} temperature = {np.round(map_data[i_val],3)} K at 408 MHz or {np.round(scaled,3)} K at {freq} MHz")
        
        return scaled

def usage():

    print("<","--"*10,">\n")
    print("-h","   the help information")
    print("e.g.: python getskytemp_v1.0.py haslam408_dsds_Remazeilles2014_ns2048.fits 10 12 1250")
    print("gl: 10 deg, gb:12 deg, obs freq: 1250 MHz")
    print("eg: python getskytemp_v1.0.py haslam408_dsds_Remazeilles2014_ns2048.fits -plot")
    print("gl: 10 deg, gb:12 deg, obs freq: 1250 MHz")
    print("\n<","--"*10,">\n")

def main():

    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
        
    filename = sys.argv[1]
    map_data, nside, ordering, coord_sys = read_healpix_map(filename)
    
    
    
    
    if "-plot" in sys.argv:
    
        gl = np.arange(-180,180.25,0.25)
        gb = np.arange(-90,90.25,0.25)
        temp = np.zeros((len(gl),len(gb)))
        
        for i in range(len(gb)):
            igb = gb[i]
            for j in range(len(gl)):
                igl = gl[j]
                temp[j][i] = get_temp(map_data, nside, ordering,igl,igb,[408])

        temp = np.array(temp)
        plt.figure(figsize=(12,9),dpi=100)
        GL, GB = np.meshgrid(gl, gb)
        plt.pcolor(GL, GB, temp.T, cmap='jet', norm=LogNorm(vmin=temp.min(), vmax=temp.max()/50))

        plt.show()
        sys.exit(1)
        
    
    
    if len(sys.argv) <= 2:
        usage()
        sys.exit(1)
    
    
    else:
        
        try:
            gl = float(sys.argv[2])
            gb = float(sys.argv[3])
            freq = float(sys.argv[4])
            get_temp(map_data, nside, ordering,gl,gb,[freq])
        
        
        except:
            usage()
            sys.exit(1)
        

if __name__ == '__main__':
    main()





#TEMP = np.reshape(temp,(len(gl),len(gb)))


# -*- coding: utf-8 -*-
"""
Infill CMIP5 atmospheric fields over land

Created on Mon Dec  7 08:39:36 2015

@author: jdha
"""
import numpy as np
from scipy.interpolate import Rbf as Rbf
from netCDF4 import Dataset
import warnings

def infill_Rbf(atm_file,var_name,msk_file):

    """ The land values in variable var_name are removed and infilled with 
        nearest neighbour averaging. The input file is then updated."""

    # set pointers to the data files
    try:
        nc_fid = Dataset(atm_file,mode='a')
    except KeyError:
        print "\t\tWARNING: %s cannot be opened" % atm_file
        
    try:
        nc_msk = Dataset(msk_file,mode='r')
    except KeyError:
        print "\t\tWARNING: %s cannot be opened" % msk_file
    
    # get variables
    
    nc_fid_var = nc_fid.variables
    nc_msk_var = nc_msk.variables

    # get dimensions
    
    nc_fid_dim = nc_fid.dimensions
    nc_msk_dim = nc_msk.dimensions 
    
    # check number of dims

    nt = len(nc_fid_dim["time"])
    nx = len(nc_fid_dim["lon"])
    ny = len(nc_fid_dim["lat"])
    
    nx_m = len(nc_msk_dim["lon"])
    ny_m = len(nc_msk_dim["lat"])
    
    # check dims match
    
    if (nx_m != nx) or (ny_m != ny) :
        print "dimensions do not match"
        exit()
    
    # produce mesh to work from
    
    x,y = np.meshgrid(np.arange(nx),np.arange(ny))
    
    # extract the lsm

    lsm = nc_msk_var["sftlf"][:,:]
    
    for t in range(nt):
        var = nc_fid_var[var_name][t,:,:]
        f = Rbf(x[lsm==0], y[lsm==0], var[lsm==0], function='linear')
        var[lsm!=0] = f(x[lsm!=0],y[lsm!=0])
        nc_fid_var[var_name][t,:,:] = var
  
    nc_fid.close()
    nc_msk.close()
    
def infill_flood(atm_file,var_name,msk_file):

    """ The land values in variable var_name are removed and infilled with 
        nearest neighbour averaging. The input file is then updated."""

    # set pointers to the data files
    try:
        nc_fid = Dataset(atm_file,mode='a')
    except KeyError:
        print "\t\tWARNING: %s cannot be opened" % atm_file
        
    try:
        nc_msk = Dataset(msk_file,mode='r')
    except KeyError:
        print "\t\tWARNING: %s cannot be opened" % msk_file
    
    # get variables
    
    nc_fid_var = nc_fid.variables
    nc_msk_var = nc_msk.variables

    # get dimensions
    
    nc_fid_dim = nc_fid.dimensions
    nc_msk_dim = nc_msk.dimensions 
    
    # check number of dims

    nt = len(nc_fid_dim["time"])
    nx = len(nc_fid_dim["lon"])
    ny = len(nc_fid_dim["lat"])
    
    nx_m = len(nc_msk_dim["lon"])
    ny_m = len(nc_msk_dim["lat"])
    
    # check dims match
    
    if (nx_m != nx) or (ny_m != ny) :
        print "dimensions do not match"
        exit()
    
    # produce mesh to work from
    
    x,y = np.meshgrid(np.arange(nx),np.arange(ny))
    
    # extract the lsm

    lsm = nc_msk_var["sftlf"][:,:]
###
    for t in range(nt):
        ind = np.where(lsm != 0)
        var = nc_fid_var[var_name][t,:,:]
        var[lsm != 0] = np.nan 
        
        while len(ind[0])>1 :

            i = ind[1] 
            j = ind[0]
        
            # create indices of neighbouring points

            
            jp1 = np.clip(j+1, 0, ny-1)
            ip1 = np.clip(i+1, 0, nx-1)
            im1 = np.clip(i-1, 0, nx-1) 
            jm1 = np.clip(j-1, 0, ny-1) 

            # create 1D indices as quicker than looping over i and j

            ind_e = ( nx * j ) + im1
            ind_w = ( nx * j ) + ip1
            ind_s = ( nx * jm1 ) + i
            ind_n = ( nx * jp1 ) + i
            ind_o = ( nx * j ) + i

            # replace missing data

            # I expect to see RuntimeWarnings in this block

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
    
                var.ravel()[ind_o] = np.nanmean(np.array([var.ravel()[ind_e], 
                                                          var.ravel()[ind_w],
                                                          var.ravel()[ind_s],
                                                          var.ravel()[ind_n]]), axis=0)

            # find new indices for next loop
    
            s_ind = len(ind[0])

            ind   = np.where(np.isnan(var))

            if len(ind[0]) == s_ind :
                
                print "no change"
                exit()
                
        tmp_pad = np.pad(var,1,'constant',constant_values=(np.nan,np.nan)) 

        tmp = np.nanmean(
                         np.concatenate( (
                                    tmp_pad[:-2,:-2,None],
                                    tmp_pad[1:-1,:-2,None],
                                    tmp_pad[2:,:-2,None],
                                    tmp_pad[:-2,1:-1,None],
                                    tmp_pad[1:-1,1:-1,None],
                                    tmp_pad[2:,1:-1,None],
                                    tmp_pad[:-2,2:,None],  
                                    tmp_pad[1:-1,2:,None],  
                                    tmp_pad[2:,2:,None] ), axis=2 ), axis=2) ;
      
        var[lsm != 0] = tmp[lsm != 0] 
        nc_fid_var[var_name][t,:,:] = var
        
    nc_fid.close()
    nc_msk.close()

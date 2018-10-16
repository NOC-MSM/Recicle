'''
Copyright 2017 James Harle 

Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do 
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.
'''

# Import required python modules
import os
import re
import ftplib
import netrc
import numpy as np
import scipy as sp
from netCDF4 import Dataset, netcdftime
from calendar import isleap
from nco import Nco
nco = Nco()
from nco.custom import Limit

def cmip_get(ddir, var_list, model_list, experiment):
    '''Simple routine to grab cmip data from the BADC archive
    '''

    # if directory doesn't exist make it
    if not os.path.isdir(ddir):
        os.mkdir(ddir)

    # change the local directory to where you want to put the data
    os.chdir(ddir)

    # get uid and pwd from file
    
    HOST = 'ftp.ceda.ac.uk'
    credentials = netrc.netrc()
    uid, acc, pwd = credentials.authenticators( HOST )

    # login to FTP
    f = ftplib.FTP(HOST, uid, pwd)

    # loop through variables
    
    for var in var_list:
        print var
        for model in model_list:
            print 'processing model: '+model
            try:
                f.cwd("/badc/cmip5/data/cmip5/output1/%s/%s/mon/ocean/Omon/r1i1p1/latest/%s/" % (model,experiment,var))
                fname = f.nlst()
                print fname
                pattern='.+.r1i1p1_.+.nc'
                for file in fname:
                    if re.match(pattern, file) :
                        f.retrbinary("RETR %s" % file, open(file,"wb").write)
            except:
                pass
    # close FTP connection
    f.close()

def process_data(src_dir,dst_dir,var):
    '''function to subsample a region of interest over a period of time and 
       output as a new single file'''
       
    # Change to source directory 
    
    os.chdir(src_dir)
    
    # Gather model metadata
    
    model_dict = get_model_info()
    
    # Loop of models to find available files are process
    
    for model in model_dict.keys():
    #for model in ['GISS-E2-R-CC',]:
        model_present = False
        fname = os.listdir('./')
        pattern = var + '.+.' + model + '_' +'.+.r1i1p1.+.nc'
        yr_st, mn_st, yr_en, mn_en = 9999, 0, 0, 12
        yr_st0, mn_st0, yr_en0, mn_en0 = 1960, 1, 2099, 12
        
        opt = [
                Limit(model_dict[model]['lon'], model_dict[model]['imin'], 
                                                model_dict[model]['imax'], 1),
                Limit(model_dict[model]['lat'], model_dict[model]['jmin'], 
                                                model_dict[model]['jmax'], 1)
              ]
        
        for file in fname:
            if re.match(pattern, file) :
                print file
                model_present = True
                yr = int(file[-16:-12])
                yr_st = min(yr_st, yr)
                if yr_st == yr:
                    mn_st = int(file[-12:-10])
                yr = int(file[-9:-5])
                yr_en = max(yr_en, yr)
                if yr_en == yr:
                    mn_en = int(file[-5:-3])
                    
                nco.ncks(input=file, output=dst_dir+'tmp_'+file, options=opt)    
                    
        if ( yr_st < yr_en and model_present ):
             
            # need to write in here to extract time based on netcdf time
            # would make the following more generic
            
            if ( yr_st >= yr_st0 ):
                mn_st0 = max(mn_st0,mn_st)
            if ( yr_en <= yr_en0 ):
                mn_en0 = min(mn_en0,mn_en)
            yr_st0 = max(yr_st0,yr_st)
            yr_en0 = min(yr_en0,yr_en)
            
            file_out = '{:s}_{:s}_rcp85_{:04d}{:02d}_{:04d}{:02d}.nc'.format(
                                         var,model,yr_st0,mn_st0,yr_en0,mn_en0)
            
            ind_st = (yr_st0 * 12 + mn_st0) - (yr_st * 12 + mn_st )
            ind_en = (yr_en0 * 12 + mn_en0) - (yr_st * 12 + mn_st )
            
            nco.ncrcat(input=dst_dir+'tmp_*.nc',output=dst_dir+'tmp.nc')
            opt = [
                    Limit('time', ind_st, ind_en, 1)
                  ]
    
            nco.ncks(input = dst_dir+'tmp.nc', output=dst_dir+file_out, 
                                                                   options=opt)
            
        for file in os.listdir(dst_dir):
            if re.match('tmp_.+.nc', file) :
                os.remove(dst_dir+file)
            if re.match('tmp.nc', file) :
                os.remove(dst_dir+file)
                
def process_atmos_data(src_dir,dst_dir,var):
    '''function to subsample a region of interest over a period of time and 
       output as a new single file'''
       
    # Change to source directory 
    
    os.chdir(src_dir)
    
    # Gather model metadata
    
    model_dict = get_model_info()
    
    # Loop of models to find available files are process
    
    for model in model_dict.keys():
    #for model in ['GISS-E2-R-CC',]:
        model_present = False
        fname = sorted(os.listdir('./'))
        pattern = var + '.+.' + model + '_' +'.+.r1i1p1.+.nc'
        yr_st, mn_st, yr_en, mn_en = 9999, 0, 0, 12
        yr_st0, mn_st0, yr_en0, mn_en0 = 1960, 1, 2099, 12
        dy_st0, dy_en0 = 1, 31 
        opt = [
                Limit("lat",  20.,80.),
                Limit("lon", 260.,60.)
              ]
        
        for file in fname:
            if re.match(pattern, file) :
                print file
                # gather meta data to work out the calendar used
                meta = nco.ncks(input=file,options=['-M', '-m'])
                cal  = re.search(r'(.?calendar.*?)\n',meta).group(1).split()[-1]
                model_present = True
                yr = int(file[-20:-16])
                yr_st = min(yr_st, yr)
                if yr_st == yr:
                    mn_st = int(file[-16:-14])
                yr = int(file[-11:-7])
                yr_en = max(yr_en, yr)
                if yr_en == yr:
                    mn_en = int(file[-7:-5])
                    dy_en = int(file[-5:-3])
                    
                nco.ncks(input=file, output=dst_dir+'tmp_'+file, options=opt)    
                    
        if ( yr_st < yr_en and model_present ):
             
            # need to write in here to extract time based on netcdf time
            # would make the following more generic
            
            if ( yr_st >= yr_st0 ):
                mn_st0 = max(mn_st0,mn_st)
            if ( yr_en <= yr_en0 ):
                mn_en0 = min(mn_en0,mn_en)
                dy_en0 = min(dy_en0,dy_en)
            yr_st0 = max(yr_st0,yr_st)
            yr_en0 = min(yr_en0,yr_en)
            
            file_out = '{:s}_{:s}_rcp85_{:04d}{:02d}{:02d}_{:04d}{:02d}{:02d}.nc'.format(
                                         var,model,yr_st0,mn_st0,dy_st0,yr_en0,mn_en0,dy_en0)
            
            tmp_cal = netcdftime.utime('days since {:04d}-{:02d}-{:02d}'.format(yr_st,mn_st,1),cal)
            
            ind_st = int(tmp_cal.date2num(netcdftime.datetime(yr_st0,mn_st0,1)))
            ind_en = int(tmp_cal.date2num(netcdftime.datetime(yr_en0,mn_en0,dy_en0)))                        
                        
            print ind_st, ind_en
            
            for yr in np.arange(yr_st0,yr_en0+1):
                
                ind_st = int(tmp_cal.date2num(netcdftime.datetime(yr,mn_st0,1)))
                ind_en = int(tmp_cal.date2num(netcdftime.datetime(yr,mn_en0,dy_en0))) 
                
                file_out = '{:s}_{:s}_rcp85_y{:04d}.nc'.format(
                                         var,model,yr)
            #nco.ncrcat(input=dst_dir+'tmp_*.nc',output=dst_dir+'tmp.nc')
                opt = [
                        Limit('time', ind_st, ind_en, 1)
                      ]
    
            #nco.ncks(input = dst_dir+'tmp.nc', output=dst_dir+file_out, 
            #                                                       options=opt)
                nco.ncrcat(input = dst_dir+'tmp_*.nc', output=dst_dir+file_out, 
                                                                   options=opt)
            
        for file in os.listdir(dst_dir):
            if re.match('tmp_.+.nc', file) :
                os.remove(dst_dir+file)
            if re.match('tmp.nc', file) :
                os.remove(dst_dir+file)
                
def time_stretch(fname,vname,year):
    
    rootgrp = Dataset(fname, "a", format="NETCDF4")
    
    nt = netcdftime.utime('days since 1960-01-01','gregorian')
    
    vid   = rootgrp.variables[vname]
    tid   = rootgrp.variables['time']
    
    if isleap(year):
        y_len = 366.
    else:
        y_len = 365.
    
    v     = vid[:]
    t_old = np.arange(0.,y_len,(y_len-1.)/(len(tid[:])-1.))
    t_new = np.arange(0,y_len)
    
    f = sp.interpolate.interp1d(t_old,v,axis=0)
    v_new = f(t_new)
    vid[:,:,:] = v_new
    
    os = nt.date2num(netcdftime.datetime(year,01,01,12))
    
    tid[:]     = t_new + os 
    tid.setncattr('calendar','gregorian')
    tid.setncattr('units','days since 1960-01-01')
    
    rootgrp.close()
    
              
def mesh_gen(tfile,ufile,vfile,file_out,model,zero_merid=True):

    # would attempt this using NCO but the ncrename doesn't appear to work with 
    # netcdf4 files correctly
    
    # example:
    # for m in model_dict.keys():
    #   os.chdir(m)
    #   cp.mesh_gen('thetao_'+m+'_rcp85_196001_209912.nc','uo_'+m+
    #               '_rcp85_196001_209912.nc','vo_'+m+
    #               '_rcp85_196001_209912.nc','fx_'+m+'.nc',m)
    #   os.chdir('../')
    
    
    model_dict = get_model_info()
       
    # create output file
 
    ncfile = Dataset(file_out, 'w', format="NETCDF4")
 
    # open first input file
    
    rootgrp = Dataset(tfile, "r", format="NETCDF4")
    
    lat = rootgrp.variables['lat'][:]
    lon = rootgrp.variables['lon'][:]
    lev = rootgrp.variables['lev'][:]
    
    if zero_merid:
        lon = np.where(lon>=180., lon-360., lon)
        lon = np.where(lon<-180., lon+360., lon)
    else:
        lon = np.where(lon>=360., lon-360., lon)
        lon = np.where(lon<   0., lon+360., lon)
    
    if lat.ndim == 2:
        ny, nx = lat.shape
    else:
        ny = len(lat)
        nx = len(lon)
        
    nz = len(lev)
       
    # create the x, y and z dimensions in the outout file

    ncfile.createDimension('x',nx)
    ncfile.createDimension('y',ny)
    ncfile.createDimension('z',nz)

    # create variables for output
    
    glam = ncfile.createVariable('glamt','f8',('y','x'))
    gphi = ncfile.createVariable('gphit','f8',('y','x'))
    gdep = ncfile.createVariable('gdept','f8',('z'))
    
    if model_dict[model]['dim'] == 1:

        glam[:], gphi[:] = np.meshgrid(lon,lat)
        
    else:
        
        gphi[:] = lat
        glam[:] = lon
        
    gdep[:] = lev
    
    # close input file
    
    rootgrp.close()

    # repeat for ufile and vfile
    
    files = {'u':ufile, 'v':vfile}
    
    for grd, fname in files.items():
        
        rootgrp = Dataset(fname, "r", format="NETCDF4")
  
        lat = rootgrp.variables['lat'][:]
        lon = rootgrp.variables['lon'][:]
        lev = rootgrp.variables['lev'][:]
        
        if zero_merid:
            lon = np.where(lon>=180., lon-360., lon)
            lon = np.where(lon<-180., lon+360., lon)
        else:
            lon = np.where(lon>=360., lon-360., lon)
            lon = np.where(lon<   0., lon+360., lon)
        
        glam = ncfile.createVariable('glam'+grd,'f8',('y','x'))
        gphi = ncfile.createVariable('gphi'+grd,'f8',('y','x'))
        gdep = ncfile.createVariable('gdep'+grd,'f8',('z'))
        
        if model_dict[model]['dim'] == 1:
            
            glam[:], gphi[:] = np.meshgrid(lon,lat)
            
        else:
                
            gphi[:] = lat
            glam[:] = lon
                
        gdep[:] = lev
        
        # close input file
        
        rootgrp.close()
        
    # add masking arrays
    
    files = {'u':ufile, 'v':vfile, 't':tfile}
    
    for grd, fname in files.items():
        
        rootgrp = Dataset(fname, "r", format="NETCDF4")
     
        if grd=='t':
            vname = 'thetao'
        else:
            vname = grd+'o'
            
        var = rootgrp.variables[vname][0,:,:,:]
        mask = ncfile.createVariable(grd+'mask','i4',('z','y','x'))
        mask[:] = (var.data<10000.) * 1
        
    # close files
    
    ncfile.close()

def get_model_info():

    model_dict = {
              'bcc-csm1-1-m':    {'lat': 'rlat', 'lon': 'rlon', 
                                  'imin':  200, 'imax':  339,
                                  'jmin':  161, 'jmax':  227,
                                  'dim' : 2, 'grd' : 'b', 
                                  'b_i1':   94, 'b_j1':   35,
                                  'b_i2':   88, 'b_j2':   38,
                                  'bcon': False,
                                  'bu_i': None, 'bu_j': None},
              "bcc-csm1-1":      {'lat': 'rlat', 'lon': 'rlon', 
                                  'imin': 200, 'imax': 339,
                                  'jmin': 161, 'jmax': 227,
                                  'dim' : 2},
              "BNU-ESM":         {'lat': 'lat', 'lon': 'lon', 
                                  'imin': 280, 'imax': 59,
                                  'jmin': 129, 'jmax': 189,
                                  'dim' : 1},
              "CanESM2":         {'lat': 'lat', 'lon': 'lon', 
                                  'imin': 199, 'imax':  42,
                                  'jmin': 117, 'jmax': 181,
                                  'dim' : 1},
              "CMCC-CESM":       {'lat': 'j', 'lon': 'i', 
                                  'imin': 102, 'imax': 170,
                                  'jmin':  92, 'jmax': 142,
                                  'dim' : 2},
              "CNRM-CM5":        {'lat': 'j', 'lon': 'i', 
                                  'imin': 208, 'imax': 347,
                                  'jmin': 182, 'jmax': 285,
                                  'dim' : 2},
              "inmcm4":          {'lat': 'rlat', 'lon': 'rlon', 
                                  'imin': 280, 'imax': 43,
                                  'jmin': 199, 'jmax': 325,
                                  'dim' : 2},
              "IPSL-CM5A-LR":    {'lat': 'j', 'lon': 'i', 
                                  'imin': 102, 'imax': 170,
                                  'jmin':  92, 'jmax': 142,
                                  'dim' : 2},
              "IPSL-CM5A-MR":    {'lat': 'j', 'lon': 'i', 
                                  'imin': 102, 'imax': 170,
                                  'jmin':  92, 'jmax': 142,
                                  'dim' : 2},
              "IPSL-CM5B-LR":    {'lat': 'j', 'lon': 'i', 
                                  'imin': 102, 'imax': 170,
                                  'jmin':  92, 'jmax': 142,
                                  'dim' : 2},
              "MIROC-ESM-CHEM":  {'lat': 'lat', 'lon': 'lon', 
                                  'imin': 199, 'imax':  42,
                                  'jmin': 126, 'jmax': 184,
                                  'dim' : 1},
              "MIROC-ESM":       {'lat': 'lat', 'lon': 'lon', 
                                  'imin': 199, 'imax':  42,
                                  'jmin': 126, 'jmax': 184,
                                  'dim' : 1},
              "HadGEM2-CC":      {'lat': 'lat', 'lon': 'lon', 
                                  'imin': 280, 'imax':  60,
                                  'jmin': 145, 'jmax': 204,
                                  'dim' : 1},
              "HadGEM2-ES":      {'lat': 'lat', 'lon': 'lon', 
                                  'imin': 280, 'imax':  60,
                                  'jmin': 145, 'jmax': 204,
                                  'dim' : 1},
              "MPI-ESM-LR":      {'lat': 'j', 'lon': 'i', 
                                  'imin':   0, 'imax': 255,
                                  'jmin':   0, 'jmax': 117,
                                  'dim' : 2},
              "MPI-ESM-MR":      {'lat': 'j', 'lon': 'i', 
                                  'imin': 439, 'imax': 748,
                                  'jmin':  10, 'jmax': 155,
                                  'dim' : 2},                          
              "GISS-E2-H-CC":    {'lat': 'lat', 'lon': 'lon',  
                                  'imin': 280, 'imax': 59,
                                  'jmin': 110, 'jmax': 169,
                                  'dim' : 1},
              "GISS-E2-R-CC":    {'lat': 'lat', 'lon': 'lon', 
                                  'imin': 224, 'imax':  47,
                                  'jmin': 110, 'jmax': 169,
                                  'dim' : 1},
              "NorESM1-ME":      {'lat': 'j', 'lon': 'i',  
                                  'imin':   0, 'imax': 319,
                                  'jmin': 255, 'jmax': 383,
                                  'dim' : 2},
              "GFDL-ESM2G":      {'lat': 'rlat', 'lon': 'rlon',  
                                  'imin': 200, 'imax': 339,
                                  'jmin': 127, 'jmax': 201,
                                  'dim' : 2},
              "GFDL-ESM2M":      {'lat': 'rlat', 'lon': 'rlon', 
                                  'imin': 200, 'imax': 339,
                                  'jmin': 129, 'jmax': 195,
                                  'dim' : 2},
              "CESM1-BGC":       {'lat': 'j', 'lon': 'i', 
                                  'imin':   0, 'imax': 319,
                                  'jmin': 255, 'jmax': 383,
                                  'dim' : 2},
             }
              
    return model_dict 

def get_model_list():
    
    model_list = ["BCC/bcc-csm1-1-m",
                  "BCC/bcc-csm1-1",
                  "BNU/BNU-ESM",
                  "CCCma/CanESM2",
                  "CMCC/CMCC-CESM",
                  "CNRM-CERFACS/CNRM-CM5",
                  "INM/inmcm4",
                  "IPSL/IPSL-CM5A-LR",
                  "IPSL/IPSL-CM5A-MR",
                  "IPSL/IPSL-CM5B-LR",
                  "MIROC/MIROC-ESM-CHEM",
                  "MIROC/MIROC-ESM",
                  "MOHC/HadGEM2-CC",
                  "MOHC/HadGEM2-ES",
                  "MPI-M/MPI-ESM-LR",
                  "MPI-M/MPI-ESM-MR",
                  "NASA-GISS/GISS-E2-H-CC",
                  "NASA-GISS/GISS-E2-R-CC",
                  "NCC/NorESM1-ME",
                  "NOAA-GFDL/GFDL-ESM2G",
                  "NOAA-GFDL/GFDL-ESM2M",
                  "NSF-DOE-NCAR/CESM1-BGC",
                 ]
    return model_list
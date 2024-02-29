#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import argparse


def parse_args():
   p = argparse.ArgumentParser(description='Flatten a lat-lon to 1D')
   p.add_argument('-i','--input',type=str,help='input file',default=None)
   p.add_argument('-o','--output',type=str,help='output file',default=None)
   return vars(p.parse_args())

                #------------------
                # Opening the file
                #------------------
comm_args    = parse_args()
Input_file   = comm_args['input']
Output_prefix  = comm_args['output']

nc_input =  Dataset(Input_file,mode='r')
cube_size = len(nc_input.dimensions['Xdim'])
nx  = cube_size*2
nxp = (cube_size*2)+1
ny  = nx
nyp = nxp


for n in range(6):
   output_file=Output_prefix+".tile"+str(n+1)+".nc"
   ncfid=Dataset(output_file,mode='w',format='NETCDF4')
   lon = nc_input.variables['lons'][:]
   lat = nc_input.variables['lats'][:]
   lon_c = nc_input.variables['corner_lons'][:]
   lat_c = nc_input.variables['corner_lats'][:]

   xdim=ncfid.createDimension('nx',nx)
   ydim=ncfid.createDimension('ny',ny)
   xdimc=ncfid.createDimension('nxp',nxp)
   ydimc=ncfid.createDimension('nyp',nyp)

   x = ncfid.createVariable('x','f8',('nyp','nxp'))
   setattr(x,'standard_name','geographic_longitude')
   setattr(x,'units','degree_east')
   y = ncfid.createVariable('y','f8',('nyp','nxp'))
   setattr(x,'standard_name','geographic_latitude')
   setattr(x,'units','degree_north')
   x[0:nx+1:2,0:nx+1:2]  = lon_c[n,:,:]
   y[0:nx+1:2,0:nx+1:2]  = lat_c[n,:,:]
   x[1:nx:2,1:nx:2]  = lon[n,:,:]
   y[1:nx:2,1:nx:2]  = lat[n,:,:]
   ncfid.close()



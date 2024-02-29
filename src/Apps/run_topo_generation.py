#!/usr/bin/env python

import os
import sys
import argparse
import subprocess as sp

# Regular cubed-sphere descriptors
_ScriptFiles_  = "/discover/nobackup/bmauer/Cubed_Sphere_Grids_new_stretched"
# Stretched cubed-sphere descriptors
#_ScriptFiles_  = "/discover/nobackup/bmauer/Cubed_Sphere_Grids_stretched"
_TopoPackage_  = "/discover/swdev/bmauer/packages/Topo-NCAR_Topo_1_1"
_c3000file_    = "output_3000.gmted_fixedanarticasuperior.nc"

def _add_zero_str(tres):

   im=int(tres)
   jm=im*6
   if im < 100:
      ims = "00"+str(im)
   elif im < 1000:
      ims = "0"+str(im)
   elif im < 10000:
      ims=str(im)

   if jm < 100:
      jms = "000"+str(jm)
   elif jm < 1000:
      jms = "00"+str(jm)
   elif jm < 10000:
      jms="0"+str(jm)
   elif jm < 100000:
      jms=str(jm)
   return ims+"x"+jms

def _bin_to_cube(topo):
   
   nl_file="&binparams\n  raw_latlon_data_file='"+topo+"'\n  output_file='"+_c3000file_+"'\n  ncube=3000\n/"
   f = open("bin_to_cube.nl",'w')
   f.write(nl_file)
   f.close()

   os.system("ln -fs "+_TopoPackage_+"/bin_to_cube/landm_coslat.nc .")
   os.system(_TopoPackage_+"/bin_to_cube/bin_to_cube")

def _cube_to_target(tres,smooth_file):

   im = int(tres)
   jm = int(tres)*6
   nl_file="&topoparams\n"
   grid_descriptor="/PE"+str(im)+"x"+str(jm)+"-CF.nc4"
   nl_file=nl_file+"  grid_descriptor_fname = '"+_ScriptFiles_+grid_descriptor+"'\n"
   nl_file=nl_file+"  intermediate_cubed_sphere_fname = '"+_c3000file_+"'\n"

   if smooth_file == "/dev/null":
      output_fname = "gmted_to_c"+str(im)+".nc"
      nl_file=nl_file+"  output_fname = 'gmted_to_c"+str(im)+".nc'\n"
      nl_file=nl_file+"  externally_smoothed_topo_file  = '"+smooth_file+"'\n"
      nl_file=nl_file+"  ladjustmean = .false.\n"
      nl_file=nl_file+"  lsmooth_terr = .false.\n"

   else:
      output_fname = "gmted_to_c"+str(im)+"_smoothed.nc"
      nl_file=nl_file+"  output_fname = 'gmted_to_c"+str(im)+"_smoothed.nc'\n"
      nl_file=nl_file+"  externally_smoothed_topo_file  = '"+smooth_file+"'\n"
      nl_file=nl_file+"  ladjustmean = .true.\n"
      nl_file=nl_file+"  lsmooth_terr = .true.\n"

   nl_file=nl_file+"  lexternal_smooth_terr = .true.\n"
   nl_file=nl_file+"  lzero_out_ocean_point_phis = .false.\n"
   nl_file=nl_file+"  lsmooth_on_cubed_sphere = .false.\n"
   nl_file=nl_file+"  ncube_sph_smooth_coarse = 20\n"
   nl_file=nl_file+"  ncube_sph_smooth_fine = 1\n"
   nl_file=nl_file+"  lfind_ridges = .false.\n"
   nl_file=nl_file+"  nwindow_halfwidth = 14\n"
   nl_file=nl_file+"  nridge_subsample = 14\n/\n"
   f = open("cube_to_target.nl","w")
   f.write(nl_file)
   f.close()
   os.system(_TopoPackage_+"/cube_to_target/cube_to_target")
   return output_fname

def _convert_to_gmao(tres,output_fname):
   print "bma converting ",output_fname
   os.system(_TopoPackage_+"/convert_to_GMAO/convert_to_output_gmao.exe -i "+output_fname+" --im "+tres)

def _do_smoothing(tres):

   smooth_bin = "/gpfsm/dswdev/bmauer/packages/UFS_UTILS/install/bin"
   grid_descrip_dir = "/discover/nobackup/bmauer/Cubed_Sphere_Grids_new_stretched_geosformat"
   mosaic_def_dir = "/discover/nobackup/bmauer/tmp/ufs_stretched_grids/my_grids"
  
   im = int(tres)
   jm = im*6
   # tile orography
   inputf = "gmted_DYN_ave_"+str(im)+"x"+str(jm)+".data"
   output_prefix = "oro.C"+tres
   run_cmd = smooth_bin+"/mosaic_topo.py" + " -i "+inputf+" -o "+output_prefix+" -c "+tres
   print run_cmd
   sp.call(run_cmd,shell=True)

   # tile grid
   inputf = grid_descrip_dir+"/c"+tres+"_coords.nc4"
   output_prefix = "C"+tres+"_grid"
   run_cmd = smooth_bin+"/mosaic_grid.py" + " -i "+inputf+" -o "+output_prefix
   sp.call(run_cmd,shell=True)

   # run smoother
   nl_file="&filter_topo_nml\n"
   nl_file=nl_file+"  grid_file = '"+mosaic_def_dir+"/C"+tres+"/C"+tres+"_mosaic.nc'\n"
   nl_file=nl_file+"  topo_file = 'oro.C"+tres+"'\n"
   nl_file=nl_file+"  mask_field = 'land_frac'\n"
   nl_file=nl_file+"  regional = .false.\n"
   nl_file=nl_file+"  stretch_fac = 1.0\n"
   nl_file=nl_file+"  res = "+tres+"\n"
   nl_file=nl_file+"  /"
   f = open("input.nml","w")
   f.write(nl_file)
   f.close()

   run_cmd = smooth_bin+"/filter_topo"
   sp.call(run_cmd,shell=True)
   foutput_format=_add_zero_str(tres)

   # convert back to form we know
   input_prefix = "oro.C"+tres
   #output_file = "topo_DYN_ave_"+foutput_format+".bin"
   output_file = "smooth_ncarformat_c"+tres+".nc"
   run_cmd = smooth_bin+"/combineMosaic.py" + " -i "+input_prefix+" -o "+output_file
   sp.call(run_cmd,shell=True)
  
   return output_file

def _convert_to_ncar(tres,inputf):

   outputf = "smooth_ncarformat_c"+tres+".nc"
   print inputf
   print outputf
   os.system(_TopoPackage_+"/convert_to_GMAO/convert_bin_to_netcdf.exe -i "+inputf+" --im "+tres+" --ncar "+outputf)
   return outputf

def _copy_final_geosoutput(outputdir,tres):
   im=int(tres)
   jm=im*6
   res_str=str(im)+"x"+str(jm)
   os.system("mv gmted*"+res_str+"* "+outputdir+"/")


if __name__ == "__main__":

   parser = argparse.ArgumentParser(description="generate topograph")

   parser.add_argument("--res",dest="tres",help="output cube resolution")
   parser.add_argument("--hres_topo",dest="hres_topo",help="high resoution topography")
   parser.add_argument("--workdir",dest="workdir",help="working directory")
   parser.add_argument("--outputdir",dest="outputdir",help="output directory")

   options = parser.parse_args()
   tres      = options.tres
   hres_topo = options.hres_topo
   workdir   = options.workdir
   outputdir = options.outputdir

   os.chdir(workdir)
   _bin_to_cube(hres_topo)

   output_fname=""
   ncar_file=""

   print "producing unsmoothed topography\n"
   output_fname = _cube_to_target(tres,"/dev/null")
   print "converting unsmoothed topography\n"
   _convert_to_gmao(tres,output_fname)

   print "smoothing topography\n"
   #output_fname = _do_smoothing(tres)
   ncar_file = _do_smoothing(tres)

   #ncar_file = _convert_to_ncar(tres,output_fname)

   if not os.path.exists(outputdir+"/unsmoothed"):
      os.mkdir(outputdir+"/unsmoothed")
   _copy_final_geosoutput(outputdir+"/unsmoothed",tres)

   print "producing smoothed topography\n"
   output_fname = _cube_to_target(tres,ncar_file)
   print "converting smoothed topography\n"
   _convert_to_gmao(tres,output_fname)
   if not os.path.exists(outputdir+"/smoothed"):
      os.mkdir(outputdir+"/smoothed")
   _copy_final_geosoutput(outputdir+"/smoothed",tres)



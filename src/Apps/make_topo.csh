#!/bin/csh -f
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --constraint=cas
#SBATCH --ntasks-per-node=45
#SBATCH --job-name=make_topo
#SBATCH --account=s1873
#SBATCH --mail-type=ALL

#setenv topo_package /home/bmauer/swdev/models/topo_workflow/GEOS_topography_workflow/install/bin
setenv topo_package /home/bmauer/swdev/packages/latest_ufs_utils/UFS_UTILS
#source /gpfsm/dswdev/bmauer/packages/UFS_UTILS/load_mod_discover.csh

set modinit = /usr/share/modules/init/csh
source $modinit
source ${topo_package}/load_spack.sh


foreach n (24)

   set output_dir = /discover/nobackup/bmauer/tmp/ufs_grid_test
   set   work_dir = /discover/nobackup/bmauer/tmp/ufs_grid_test
   if (! -e $output_dir) then
      mkdir $output_dir
   endif
   if (! -e $work_dir) then
      mkdir $work_dir
   endif
   cd $work_dir

   ${topo_package}/driver_scripts/driver_grid.discover.sh
   exit
   cd $work_dir
   set config_file = GenScrip.rc
   @ jm = $n * 6
   echo $n
   echo $jm
   set scripfile = PE${n}x${jm}-CF.nc4
   echo $scripfile
   cat << _EOF_ > ${config_file}
CUBE_DIM: $n
output_scrip: ${scripfile}
output_geos: geos_cube.nc4
_EOF_
   mpirun -np 6 ${topo_package}/install/bin/generate_scrip_cube.x
   rm GenScrip.rc

   exit
   cd $output_dir
   echo "running c$n"
   ${topo_package}/run_topo_generation.py --res $n --hres_topo /discover/nobackup/bmauer/gmted_topo/gmted_fix_superior/gmted_fixed_anartica_superior_caspian.nc4 --outputdir $output_dir --topo_install $topo_package --workdir $work_dir

end


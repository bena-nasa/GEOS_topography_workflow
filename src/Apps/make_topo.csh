#!/bin/csh -f
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --constraint=cas
#SBATCH --ntasks-per-node=45
#SBATCH --job-name=make_topo
#SBATCH --account=s1873
#SBATCH --mail-type=ALL

setenv topo_package /home/bmauer/swdev/models/topo_workflow/GEOS_topography_workflow/install/bin
#source /gpfsm/dswdev/bmauer/packages/UFS_UTILS/load_mod_discover.csh
source ${topo_package}/g5_modules

foreach n (180)

   set output_dir = /discover/nobackup/bmauer/tmp/topo_work
   if (! -e $output_dir) then
      mkdir $output_dir
   endif
   cd $output_dir
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

   mpirun -np 6 ${topo_package}/generate_scrip_cube.x

   echo "running c$n"
   ${topo_package}/run_topo_generation.py --res $n --hres_topo /discover/nobackup/bmauer/gmted_topo/gmted_fix_superior/gmted_fixed_anartica_superior_caspian.nc4 --outputdir $output_dir --topo_install $topo_package

end


#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --constraint=cas
#SBATCH --ntasks-per-node=45
#SBATCH --job-name=make_topo
#SBATCH --account=s1873
#SBATCH --mail-type=ALL

#setenv topo_package /home/bmauer/swdev/models/topo_workflow/GEOS_topography_workflow/install/bin
export topo_package=/home/bmauer/swdev/packages/latest_ufs_utils/UFS_UTILS
export topo_bin=/home/bmauer/swdev/packages/latest_ufs_utils/UFS_UTILS/install/bin
#source /gpfsm/dswdev/bmauer/packages/UFS_UTILS/load_mod_discover.csh

export modinit=/usr/share/modules/init/bash
source $modinit
source ${topo_package}/load_spack.sh
module load python/3.10.13
module load py-pyyaml
module load py-numpy
module load py-scipy
module load py-netcdf4


#for n in 1120;
for n in 90;
do

   export output_dir=/discover/nobackup/bmauer/tmp/ufs_grid_test_c90
   export   work_dir=/discover/nobackup/bmauer/tmp/ufs_grid_test_c90
   if [[ ! -e $output_dir ]]; then
      mkdir $output_dir
   fi 
   if [[ ! -e $work_dir ]]; then
      mkdir $work_dir
   fi

   #cd $work_dir
   #export res=$n
   #export other things like schmidt, directories to use
   #${topo_package}/driver_scripts/driver_grid.discover.sh

   cd $work_dir

   config_file=GenScrip.rc
   let jm=$n*6
   echo $n
   echo $jm
   scripfile=PE${n}x${jm}-CF.nc4
   echo $scripfile
   echo $config_file
   cat << _EOF_ > ${config_file}
CUBE_DIM: $n
output_scrip: ${scripfile}
output_geos: c${n}_coords.nc4
_EOF_
   mpirun -np 6 ${topo_package}/install/bin/generate_scrip_cube.x
   rm GenScrip.rc

   cd $output_dir
   echo "running c$n"
   ${topo_package}/install/bin/run_topo_generation.py --res $n --hres_topo /discover/nobackup/bmauer/gmted_topo/gmted_fix_superior/gmted_fixed_anartica_superior_caspian.nc4 --outputdir $output_dir --topo_install $topo_bin --workdir $work_dir

done


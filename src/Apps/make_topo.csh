#!/bin/csh -f
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --constraint=cas
#SBATCH --ntasks-per-node=45
#SBATCH --job-name=make_topo
#SBATCH --account=s1873
#SBATCH --mail-type=ALL


source /gpfsm/dswdev/bmauer/packages/UFS_UTILS/load_mod_discover.csh

foreach n (270 540 1080 2160)

   set output_dir = /discover/nobackup/bmauer/gmted_topo/NCAR_TOPO_GMTED_UFS_SMOOTHING_STRETCH/c$n
   if (! -e $output_dir) then
      mkdir $output_dir
   endif
   echo "running c$n"
   ./run_stretched_topography_ufssmoothing.py --res $n --hres_topo /discover/nobackup/bmauer/gmted_topo/gmted_fix_superior/gmted_fixed_anartica_superior_caspian.nc4 --workdir /discover/nobackup/bmauer/tmp/topo_working_dir_ufs_smooth --outputdir $output_dir

end


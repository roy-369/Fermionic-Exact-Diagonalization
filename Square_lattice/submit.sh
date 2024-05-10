#!/bin/bash

# Nx 2 - 4 ,Nxstep = 1
# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10
source activate ~/anaconda3/envs/pyenv

for Nx in $(seq 3 1 4)
do
	mkdir -p Nx_$Nx
	cd Nx_$Nx
       
        for Ny in $(seq 2 1 4)
        do
               mkdir -p Ny_$Ny
               cd Ny_$Ny

               for U in $(seq 10 1 10)
	       do
	        	mkdir -p U_$U
	        	cd U_$U
          
                        for t1 in $(seq 0 0.5 1)
                        do 
                   
                        	mkdir -p t1_$t1
                                cd t1_$t1          

mydir=$(pwd)
rawjob="$(cat <<EOF
#!/bin/sh
#SBATCH --job-name=HF-Nx$Nx-Ny$Ny
#SBATCH -N 1
#SBATCH --partition=trivedi
#SBATCH --ntasks=8
#SBATCH -t 168:00:00 			
#SBATCH --mem=168G
#SBATCH --mail-type=FAIL			
#SBATCH --mail-user=roy.369@osu.edu				
hostname
#SBATCH --no-requeue
source activate ~/anaconda3/envs/pyenv
cd \$SLURM_SUBMIT_DIR
time
#source ~/.bashrc
python /home/roy.369/Research/Mott_Insulators/2_dimension/Square_lattice/Post_March_Meeting/FHM_2d_spectrum_OBC_cluster.py $Nx $Ny $U $t1 2000 &> output.out
time
EOF
)"	
         		echo "$rawjob" &> job.bat
	        	sbatch job.bat

		        	cd ..
	        	done	
	        	cd ..
                 done
                 cd ..
        done
        cd ..
	echo "Nx_$Nx, Ny_$Ny, U_$U, t1_$t1"
done





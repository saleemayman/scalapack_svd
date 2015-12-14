#!/bin/bash
#@ job_name = pos_ag
# job_type = parallel
#@ job_type = MPICH
#@ output = job.out
#@ error = job.err
#@ class = test
#@ node = 1
#@ total_tasks = 1
#@ energy_policy_tag = di29vud2_energy_tag
#@ minimize_time_to_solution = yes
#@ island_count = 1
#@ wall_clock_limit = 00:30:00
#@ queue
source /etc/profile.d/modules.sh
module load scalapack
module load valgrind
#module load mpi.intel

#mpiexec -n 1 valgrind --leak-check=full ./bin/testRead ../../data/mymatrix.dat , 1 1 9 9 2
mpiexec -n 4 ./bin/testRead ../../data/mymatrix.dat , 2 2 9 9 2
echo "The job has finished!"

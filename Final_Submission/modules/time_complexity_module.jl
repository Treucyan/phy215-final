module time_complexity
using Distributions
include("serialized_dla_modules.jl")
include("parallelized_DLA.jl")

#Function for calculating the Run Time of Serial DLA
"""
# Description
Calculates and saves the run time for the DLA algorithm given a number of particles for three trials. 
Returns a Nx4 array where the first column contains the the number of particles N and the succeeding columns are the elapsed time to run the DLA algorithm for the given N.

## Args
    particle_number (Int64): number of particles used to build the cluster
    maximum_radius (Float64): maximum radius the death circle and the birth 
        circle can grow into.
    sticking_prob (Float64): probability of a particle to stick on the cluster.
    N_points (Int64): number of data points desired

## Returns
    serial_time (Matrix{Float64}): 1 x 4 array with first column containing particle_number and second to fourth column containing elapsed time to run DLA algorithm
"""

function run_time_serial(particle_number::Int64, maximum_radius::Float64, sticking_prob::Float64)
    serial_time = zeros(1,4) #Memory Allocation
    for i in 1:size(serial_time,1)
        serial_time[i,1] += particle_number
        for j in 2:4
            time_cluster = @elapsed Random_walker.serialized_dla(particle_number, maximum_radius, sticking_prob)
            serial_time[i,j] += time_cluster
        end
    end
    return serial_time

end


#Function for calculating the Run Time of Serial DLA
"""
# Description
Calculates and saves the run time for the DLA algorithm given a number of particles for three trials. 
Returns a Nx4 array where the first column contains the the number of particles N and the succeeding columns are the elapsed time to run the DLA algorithm for the given N.

## Args
    desired_cluster_mass (Int64): desired number of particles in the cluster
    spawn_density (Float64): Must be between 0 and 1, exclusive (usually set to 0.8).
Prevents there being more particles than choices of where to place them on the birth circle, especially towards the start of the simulation.
    num_steps (Int64): Set to a large number >= 10_000. Number of steps that each particle can walk when generated in parallel.
    max_radius (Int64): Set to about 150-300. Maximum radius of the birth and death circles.
    spawn_cap (Int64): Caps the maximum number of particles that can be generated at one time.
    deposit_prob (Float64): Must be between 0 and 1, inclusive. The probability that a particle attaches to the cluster when there is contact between the two.

## Returns
    parallel_time (Matrix{Float64}): 1 x 4 array with first column containing particle_number and second to fourth column containing elapsed time to run DLA algorithm
"""

function run_time_parallel(desired_cluster_mass::Int64, spawn_density::Float64, num_steps::Int64, max_radius::Int64, spawn_cap::Int64, deposit_prob::Int64)
    parallel_time = zeros(1,4) #Memory Allocation
    for i in 1:size(parallel_time,1)
        parallel_time[i,1] += desired_cluster_mass
        for j in 2:4
            time_cluster = @elapsed DLA_parallel.run_parallel_DLA(
                                        desired_cluster_mass, spawn_density,
                                        num_steps, max_radius, spawn_cap, 
                                        deposit_prob)
            parallel_time[i,j] += time_cluster
        end
    end
    return parallel_time

end

#Function for calculating the Run Time Mean and Standard Deviation
"""
# Description
Calculates the mean and standard deviation of the run time given an array of number of particles and the corresponding three trials of run time.
Returns a N_pointsx3 array where the first column contains the the number of particles N, second column contains mean run time, and third column contains standard deviation of run time.

## Args
    serial_time (Matrix{Float64}): array of particle number and run time of three trials.

## Returns
   time_stats (Matrix{Float64}): size(run_time,1) x 3 array with first column containing particle_number, and second column contains mean run time, and third column contains standard deviation of run time.
"""

function run_time_stats(run_time::Matrix{Float64})
    time_stats = zeros(size(run_time,1),3)
    for i in 1:size(time_stats,1)
        time_stats[i,1] = run_time[i,1]
        time_stats[i,2] = mean(run_time[i,2:4])
        time_stats[i,3] = std(run_time[i,2:4])
    end
    return time_stats

end

end

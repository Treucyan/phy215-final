module time_complexity
using Distributions
include("serialized_dla_modules.jl")

#Function for calculating the Run Time
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
    serial_time (Matrix{Float64}): N_points x 4 array with first column containing particle_number and second to fourth column containing elapsed time to run DLA algorithm
"""

function run_time(particle_number::Int64, maximum_radius::Float64, sticking_prob::Float64,N_points::Int64)
    particle_number_array = [k for k in range(0,particle_number,N_points+1)]
    serial_time = zeros(N_points,4) #Memory Allocation
    for i in 1:size(serial_time,1)
        L_cluster = Int(particle_number_array[i+1])
        serial_time[i,1] += L_cluster
        for j in 2:4
            time_cluster = @elapsed Random_walker.serialized_dla(L_cluster, maximum_radius, sticking_prob)
            serial_time[i,j] += time_cluster
        end
    end
    return serial_time

end

#Function for calculating the Run Time Mean and Standard Deviation
"""
# Description
Calculates the mean and standard deviation of the run time given an array of number of particles and the corresponding three trials of run time.
Returns a N_pointsx3 array where the first column contains the the number of particles N, second column contains mean run time, and third column contains standard deviation of run time.

## Args
    serial_time (Matrix{Float64}): array of particle number and run time of three trials.

## Returns
    serial_time_stats (Matrix{Float64}): size(serial_time,1) x 3 array with first column containing particle_number, and second column contains mean run time, and third column contains standard deviation of run time.
"""

function run_time_stats(run_time::Matrix{Float64})
    serial_time_stats = zeros(size(run_time,1),3)
    for i in 1:size(serial_time_stats,1)
        serial_time_stats[i,1] = run_time[i,1]
        serial_time_stats[i,2] = mean(run_time[i,2:4])
        serial_time_stats[i,3] = std(run_time[i,2:4])
    end
    return serial_time_stats

end

end

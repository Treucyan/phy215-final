module time_complexity_calculation
using Distributions
using DelimitedFiles
include("serialized_dla_modules.jl")
include("parallelized_DLA.jl")
include("time_complexity_module.jl")

#Function for saving the Run Time of Serial DLA
"""
# Description
Saves the run time for the DLA algorithm given a number of particles for three trials. 
Returns a 1x4 array where the first column contains the the number of particles N and the succeeding columns are the elapsed time to run the DLA algorithm for the given N.
"""
function serial_time_save(init,final,step)
    init_particle_number = init
    final_particle_number = final
    steps = step
    for i = init_particle_number:steps:final_particle_number
        particle_number = i
        serial_time = time_complexity.run_time_serial(particle_number,100.0,1.0)
        writedlm("raw_data/serial_time particle_number = $particle_number.txt", serial_time)
    end 
end

#Function for saving the Run Time of Parallel DLA
"""
# Description
Saves the run time for the DLA algorithm given a number of particles for three trials. 
Returns a 1x4 array where the first column contains the the number of particles N and the succeeding columns are the elapsed time to run the DLA algorithm for the given N.
"""
function parallel_time_save(init,final,step)
    init_particle_number = init
    final_particle_number = final
    steps = step
    for i = init_particle_number:steps:final_particle_number
        particle_number = i
        parallel_time = time_complexity.run_time_parallel(particle_number,0.8,10000,150,300,1)
        writedlm("raw_data/parallel_time particle_number = $particle_number.txt", parallel_time)
    end
end


end

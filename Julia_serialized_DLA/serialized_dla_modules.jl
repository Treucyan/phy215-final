module Random_walker
using Distributions


#function for creating a random walker
"""
# Description
Initializes a random walker in a birth ring with an inner radius of inner_walk_radius
and outer radius of outer_walk_radius.

## Args
    birth_radius (Float64): radial distance where the walker is initialized.

## Returns
    [floor(x_position), floor(y_position)] (Vector{Float64}): initial postion of the walker.
        within the square lattice 
"""
function initialize_randomwalker(birth_radius::Float64)
    anglular_position = rand(Uniform(-pi, pi))
    x_position = birth_radius * cos(anglular_position)
    y_position = birth_radius * sin(anglular_position)

    return [floor(x_position), floor(y_position)]

end


#function for updating position of random walker
"""
# Description
Updates the position of the walker for a single step, whether it will go up, down, left, or right.
    If the walker exits the walking circle of radius `death_radius`, the walker is re-initialize_randomwalker
    back to the birth ring with outer radius of `death_radius` and inner radius of `birth_radius`.

## Args
    walker (Vector{Float64}): current position of the walker
    death_radius (Float64): radius of the death circle where the walker is allowed to walk. It is
        also the outer radius of the birth ring.
    birth_radius (Float64): radial distance where the walker is initialized.

## Returns
    new_position (Vector{Float64}): updated position of the walker

"""
function walker_update_position(walker::Vector{Float64}, death_radius::Float64, birth_radius::Float64)

    #updating position 
    x_or_y_update = rand(0:1)
    if x_or_y_update == 0
        new_position = walker + [rand([-1.0, 1.0]), 0.0]
    else
        new_position = walker + [0.0, rand([-1.0, 1.0])]
    end
    
    if new_position[1]^2 + new_position[2]^2 > death_radius^2
        new_position = initialize_randomwalker(birth_radius)
    end

    return new_position 
    
end



#sample function for generating a random walk
"""
# Description
Creates an array containing a single random walk trajectory of a single particle.

## Args
    step_number (Int64): number of steps a particle does in the random walk.
    death_radius (Float64): radius of the death circle where the walker is allowed to walk. It is
        also the outer radius of the birth ring.
    birth_radius (Float64): inner radius of the birth ring

## Returns
    walker_trajectory (Matrix{Int64}): 2 x step_number matrix containing the x and y position 
        of the particle throughout the random walk.
"""
function random_walk_generator(step_number::Int64, death_radius::Float64, birth_radius::Float64)

    walker_trajectory = zeros(Int64, (2, step_number))
    walker_position = initialize_randomwalker(birth_radius)
    for i in 1:step_number
        walker_trajectory[:, i] = walker_position
        walker_position = walker_update_position(walker_position, death_radius, birth_radius)
    end
    return walker_trajectory
end




function serialized_dla(particle_number::Int64, maximum_radius::Float64)

    #initializing constants
    cluster_aggregate = zeros(Float64, (2, particle_number + 1))


    #variables to be dynamically updated
    cluster_particle_number = 1
    death_radius = 1.0
    birth_radius = 0.0

    #creating the cluster
    for particle in 1:particle_number
        walker_position = initialize_randomwalker(birth_radius)
        far_from_cluster = true
        
        

        #random walk of a single particle
        while far_from_cluster

            #checking whether the system entered 
            walker_cluster_vector = (cluster_aggregate .- walker_position).^2
            distance_from_cluster = zeros(Float64, cluster_particle_number)

            for i in 1:cluster_particle_number
                distance_from_cluster[i] = sum(walker_cluster_vector[:, i])
            end

            if abs(minimum(distance_from_cluster) - 1) <= 1e-6
                cluster_aggregate[:, cluster_particle_number] = walker_position
                cluster_particle_number += 1
                far_from_cluster = false
            end
            
            walker_position = walker_update_position(walker_position, death_radius, birth_radius)
        end

        #calculating cluster particle distances 
        cluster_particle_distance_array = zeros(Float64, cluster_particle_number)
        for i in 1 : cluster_particle_number
            cluster_particle_distance_array[i] = sum(cluster_aggregate[:, i].^2)
        end

        
        if death_radius < maximum_radius
            birth_radius = maximum(cluster_particle_distance_array).^0.5
            death_radius = 2.0 * birth_radius
        else
            if birth_radius <= maximum_radius
                birth_radius = maximum(cluster_particle_distance_array).^0.5
                death_radius = maximum_radius
            else
                birth_radius = maximum_radius
                death_radius = maximum_radius
            end
        end

        print("\r Walker # $particle done!")
    end

    return cluster_aggregate 
end


end #end of Random_walker module

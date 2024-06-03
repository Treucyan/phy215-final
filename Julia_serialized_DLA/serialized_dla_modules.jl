module Random_walker
using Distributions


#function for creating a random walker
"""
# Description
Initializes a random walker in a birth ring with an inner radius of inner_walk_radius
and outer radius of outer_walk_radius.

## Args
    outer_walk_radius (Float64): outer radius of the birth ring.
    inner_walk_radius (Float64): inner radius of the birth ring. 

## Returns
    [floor(x_position), floor(y_position)] (Vector{Float64}): initial postion of the walker.
        within the square lattice 
"""
function initialize_randomwalker(outer_birth_radius::Float64, inner_birth_radius::Float64)
    radial_position = rand(Uniform(inner_birth_radius, outer_birth_radius))
    anglular_position = rand(Uniform(-pi, pi))
    x_position = radial_position * cos(anglular_position)
    y_position = radial_position * sin(anglular_position)

    return [floor(x_position), floor(y_position)]

end


#function for updating position of random walker
"""
# Description
Updates the position of the walker for a single step, whether it will go up, down, left, or right.
    If the walker exits the death circle of radius `death_radius`, the walker is re-initialize_randomwalker
    back to the birth ring with outer radius of `outer_birth_radius` and inner radiues of `inner_birth_radius`.

## Args
    walker (Vector{Float64}): current position of the walker
    death_radius (Float64): radius of the death circle where the walker is allowed to walk
    outer_birth_radius (Float64): outer radius of the birth ring
    inner_birth_radius (Float64): inner radius of the birth ring

## Returns
    new_position (Vector{Float64}): updated position of the walker

"""
function walker_update_position(walker::Vector{Float64}, death_radius::Float64, outer_birth_radius::Float64,
    inner_birth_radius::Float64)

    #updating position 
    x_or_y_update = rand(0:1)
    if x_or_y_update == 0
        new_position = walker + [rand([-1.0, 1.0]), 0.0]
    else
        new_position = walker + [0.0, rand([-1.0, 1.0])]
    end
    
    if new_position[1]^2 + new_position[2]^2 > death_radius^2
        new_position = initialize_randomwalker(outer_birth_radius, inner_birth_radius)
    end

    return new_position 
    
end



#sample function for generating a random walk
"""
# Description
Creates an array containing a single random walk trajectory of a single particle.

## Args
    step_number (Int64): number of steps a particle does in the random walk.
    death_radius (Float64): radius of the death circle where the walker is allowed to walk
    outer_birth_radius (Float64): outer radius of the birth ring
    inner_birth_radius (Float64): inner radius of the birth ring

## Returns
    walker_trajectory (Matrix{Int64}): 2 x step_number matrix containing the x and y position 
        of the particle throughout the random walk.
"""
function random_walk_generator(step_number::Int64, death_radius::Float64, outer_birth_radius::Float64, 
    inner_birth_radius::Float64)

    walker_trajectory = zeros(Int64, (2, step_number))
    walker_position = initialize_randomwalker(outer_birth_radius, inner_birth_radius)
    for i in 1:step_number
        walker_trajectory[:, i] = walker_position
        walker_position = walker_update_position(walker_position, death_radius, outer_birth_radius, inner_birth_radius)
    end
    return walker_trajectory
end




function serialized_dla(particle_number::Int64, starting_death_radius::Float64, starting_outer_birth_radius::Float64)

    #initializing constants
    cluster_aggregate = zeros(Float64, (2, particle_number + 1))


    #variables to be dynamically updated
    cluster_particle_number = 1
    death_radius = starting_death_radius
    outer_birth_radius = starting_outer_birth_radius
    inner_birth_radius = 0.0
    

    #creating the cluster
    for particle in 1:particle_number
        walker_position = initialize_randomwalker(outer_birth_radius, inner_birth_radius)
        far_from_cluster = true

        #random walk of a single particle
        while far_from_cluster
            for i in 1 : cluster_particle_number
                distance_from_cluster_particle = sum((cluster_aggregate[:, i] - walker_position).^2)
                if abs(distance_from_cluster_particle - 1) < 1e-6
                    cluster_particle_number += 1

                    cluster_aggregate[:, cluster_particle_number] = walker_position
                    far_from_cluster = false
                    break
                end
            end 
            walker_position = walker_update_position(walker_position, death_radius, outer_birth_radius, inner_birth_radius)
        end


        #updating birth and death birth_radius
        cluster_particle_distance_array = zeros(Float64, cluster_particle_number)
        for i in 1 : cluster_particle_number
            cluster_particle_distance_array[i] = sum(cluster_aggregate[:, i].^2)
        end

        inner_birth_radius = maximum(cluster_particle_distance_array).^0.5
        if inner_birth_radius^2 > (0.5 * outer_birth_radius)^2
            outer_birth_radius += 0.5 * inner_birth_radius
        end

        if outer_birth_radius^2 > (0.75 * death_radius)^2
            death_radius += 10
        end

        print("\r Walker # $particle done!")
    end

    return cluster_aggregate 
end


end #end of Random_walker module


# using ..Random_walker

# const outer_birth_radius = 10.0
# const inner_birth_radius = 5.0
# const death_radius = 20.0
# const step_number = 3

# initial_position = Random_walker.initialize_randomwalker(outer_birth_radius, inner_birth_radius)
# updated_position = Random_walker.walker_update_position(initial_position, death_radius, outer_birth_radius, inner_birth_radius)
# random_walk = Random_walker.random_walk_generator(step_number, death_radius, outer_birth_radius, inner_birth_radius)

# println(initial_position, typeof(initial_position) )
# println(updated_position, typeof(updated_position))
# println(typeof(random_walk), "", random_walk)
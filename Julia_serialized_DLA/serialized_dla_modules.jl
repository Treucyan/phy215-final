module Random_walker



#function for creating a random walker
function initialize_randomwalker(birth_radius::Int64)
    return rand(-1 * birth_radius : birth_radius, 2)
end


#function for updating position of random walker
function walker_update_position(walker::Array, death_radius::Int64, birth_radius::Int64)
    #updating position 
    x_or_y_update = rand(0:1)
    if x_or_y_update == 0
        new_position = walker + [rand([-1, 1]), 0]
    else
        new_position = walker + [0, rand([-1, 1])]
    end
    
    if new_position[1]^2 + new_position[2]^2 > death_radius^2
        new_position = initialize_randomwalker(birth_radius)
    end

    return new_position 
    
end

#sample function for generating a random walk
function random_walk_generator(step_number::Int64, death_radius::Int64, birth_radius::Int64)
    walker_trajectory = zeros(Int64, (2, step_number))
    walker_position = initialize_randomwalker(birth_radius)
    for i in 1:step_number
        walker_trajectory[:, i] = walker_position
        walker_position = walker_update_position(walker_position, death_radius, birth_radius)
    end
    return walker_trajectory
end




function serialized_dla(starting_birth_radius::Int64, starting_death_radius::Int64, particle_number::Int64)
    cluster_aggregate = zeros(Int64, (2, particle_number + 1))


    #variables to be dynamically updated
    cluster_particle_number = 1
    birth_radius = starting_birth_radius
    death_radius = starting_death_radius


    #creating the cluster
    for particle in 1:particle_number
        walker_position = initialize_randomwalker(birth_radius)
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
            walker_position = walker_update_position(walker_position, death_radius, birth_radius)
        end


        #updating birth and death birth_radius
        cluster_particle_distance_array = zeros(Float64, cluster_particle_number)
        for i in 1 : cluster_particle_number
            cluster_particle_distance_array[i] = sum(cluster_aggregate[:, i].^2)
        end

        if maximum(cluster_particle_distance_array) > birth_radius^2
            birth_radius += 1
            if birth_radius^2 > 0.5 * death_radius^2
                death_radius += 10
            end
        end

        print("\r Walker # $particle done!")
    end

    return cluster_aggregate 
end


end #end of Random_walker module


# using ..Random_walker

# const starting_birth_radius = 5
# const starting_death_radiua = 10
# const particle_number = 3

# println(Random_walker.serialized_dla(starting_birth_radius, starting_death_radiua, particle_number))
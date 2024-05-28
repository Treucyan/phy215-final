#creating object for a random walker

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


function random_walk_generator(step_number::Int64, death_radius::Int64, birth_radius::Int64)
    walker_trajectory = zeros(Int64, (2, step_number))
    walker_position = initialize_randomwalker(birth_radius)
    for i in 1:step_number
        walker_trajectory[:, i] = walker_position
        walker_position = walker_update_position(walker_position, death_radius, birth_radius)
    end
    return walker_trajectory
end


end #end of Random_walker module



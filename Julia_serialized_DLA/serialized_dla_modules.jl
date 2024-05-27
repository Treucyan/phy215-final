#creating object for a random walker

module Random_walker

#function for creating a random walker
function initialize_randomwalker(birth_radius::Int64)
    return rand(-1 * birth_radius : birth_radius, 2)
end


#function for updating position of random walker
function walker_update_position(walker::Array, death_radius::Int64)
    #updating position 
    x_or_y_update = rand(0:1)
    if x_or_y_update == 0
        new_position = walker + [rand([-1, 1]), 0]
    else
        new_position = walker + [0, rand([-1, 1])]
    end
    
    if new_position[1]^2 + new_position[2]^2 > death_radius^2
        new_position = initialize_randomwalker(death_radius)
    end

    return new_position 
    
end


function random_walk_generator(step_number::Int64)
    print(0)
end


end #end of Random_walker module






#main file
using ..Random_walker

const death_radius = 10
const birth_radius = 2

first_walker = Random_walker.initialize_randomwalker(birth_radius)
walker_first_step = Random_walker.walker_update_position(first_walker, death_radius)

println("function is working")
println(first_walker, " next step : ", walker_first_step)


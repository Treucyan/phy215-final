#creating object for a random walker
struct RandomWalker
    x_position::Int64
    y_position::Int64
end


function initialize_randomwalker(death_radius::Int64)
    return RandomWalker(rand(0:death_radius), rand(0:death_radius))
end

first_walker = initialize_randomwalker(10)

println("function is working")
println("x position: ", first_walker.x_position, ";  y position: ",first_walker.y_position);
# ----------------------------------------------
#  Utils
# ----------------------------------------------
module Utils
using Memoize

"""
    get_circle_points(r)

Implements Bresenham's midpoint circle algorithm to obtain the grid points that best approximate a circle of radius `r`.

# Source
Linus Arver (2021). https://funloop.org/post/2021-03-15-bresenham-circle-drawing-algorithm.html
"""
@memoize function get_circle_points(r)
    points = []
    sizehint!(points, Int(floor(2*pi*r)))  # upper bound
    
    x = 0
    y = -r
    F_M = 1 - r
    d_e = 3
    d_ne = -(r << 1) + 5

    push!(points, [x, y])
    push!(points, [-x, y])
    push!(points, [x, -y])
    push!(points, [-x, -y])
    
    push!(points, [y, x])
    push!(points, [-y, x])
    push!(points, [y, -x])
    push!(points, [-y, -x])

    while x < -y
        if F_M <= 0
            F_M += d_e
        else
            F_M += d_ne
            d_ne += 2
            y += 1
        end
        d_e += 2
        d_ne += 2
        x += 1
        
        push!(points, [x, y])
        push!(points, [-x, y])
        push!(points, [x, -y])
        push!(points, [-x, -y])
    
        push!(points, [y, x])
        push!(points, [-y, x])
        push!(points, [y, -x])
        push!(points, [-y, -x])
    end
    return unique(points)
end

end  # module end


# ----------------------------------------------
#  DLA_parallel
# ----------------------------------------------
module DLA_parallel
using LinearAlgebra
using ..Utils


"""
    generate_trajectories(num_walkers, num_steps, birth_radius)

Returns a matrix of vectors (with dimension num_steps+1 by num_walkers), where each vector represents a site visited by a randomly walking particle. 
Each column of the matrix comprises the trajectory of a particular particle which starts on a random grid point on a circle of radius `birth_radius`.
"""
function generate_trajectories(num_walkers, num_steps, birth_radius)
    directions = [[1, 0], [-1, 0], [0, 1], [0, -1]]  # constrains the movement of the particle to lattice sites

    trajectories = rand(directions, (num_steps+1, num_walkers))
    trajectories[1, :] = rand(Utils.get_circle_points(birth_radius), num_walkers)
    
    return cumsum(trajectories, dims=1)
end

"""
    is_adjacent_to_cluster(step, cluster)

Checks whether the position of randomly walking particle `step` is adjacent to any of the points in the cluster.
"""
function is_adjacent_to_cluster(step, cluster)
    directions = [[1,0], [-1,0], [0,1], [0,-1]]
    for deposit in cluster
        for direction in directions
            if step == (deposit + direction)
                return true
            end
        end
    end
    return false
end


"""
    step_parallel_DLA(trajectories, cluster)

Executes one step of the parallelized DLA algorithm where each particle may interfere.
"""
function step_parallel_DLA(trajectories, cluster)
    for trajectory in eachcol(trajectories)
        for step in trajectory
            if is_adjacent_to_cluster(step, cluster)
                push!(cluster, step)
                break
            end
            # TODO: Skip trajectory if it passes some death circle
            # though, this might need re-calcuating cluster_radius at each attachment event
        end
    end
    return cluster
end


"""
"""
function parallelized_DLA(num_walkers, num_steps, desired_cluster_size, birth_radius)

    # Initialize cluster and place seed at the origin
    cluster = []
    sizehint!(cluster, desired_cluster_size+1)
    push!(cluster, [0,0])

    while length(cluster) < desired_cluster_size
        cluster_radius = maximum(map(x -> norm(x), cluster))
        trajectories = generate_trajectories(num_walkers, num_steps, birth_radius)
        cluster = step_parallel_DLA(trajectories, cluster)

        if birth_radius < cluster_radius
            birth_radius = Int(floor(cluster_radius + 3))  # arbitrary choice of 3
        end
        print("\r The current size of the cluster is $(length(cluster)). ")
    end
    return cluster
end


end  # module end
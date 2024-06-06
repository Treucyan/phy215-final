# ----------------------------------------------
#  Utils
# ----------------------------------------------
module Utils
using Memoize

"""
    get_circle_points(r)

Implements Bresenham's midpoint circle algorithm to obtain the grid points that best approximate a circle of radius `r`.

# Source
Rosetta code (GNU Free Document License 1.3) [Original code translated from Python to Julia]. https://rosettacode.org/wiki/Bitmap/Midpoint_circle_algorithm#Python
"""
@memoize function get_circle_points(r)
    points = []
    sizehint!(points, Int(floor(5.66*r)))  # constant determined by linear regression
    
    x = 0
    y = r
    f= 1 - r
    ddf_x = 1
    ddf_y = -2 * r

    push!(points, [x, y])
    push!(points, [-x, y])
    push!(points, [x, -y])
    push!(points, [-x, -y])
    
    push!(points, [y, x])
    push!(points, [-y, x])
    push!(points, [y, -x])
    push!(points, [-y, -x])

    while x < y
        if f >= 0
            y -= 1
            ddf_y += 2
            f += ddf_y
        end
        x += 1
        ddf_x += 2
        f += ddf_x
        
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
    walker_cluster_distances = norm.(cluster .- Ref(step))
    if abs(minimum(walker_cluster_distances) - 1) <= 1e-6
        return true
    end
    return false
end

"""
    run_parallel_DLA(desired_cluster_mass, spawn_density,  num_steps, max_radius, spawn_cap, deposit_prob)

Args:
- `desired_cluster_mass::Int`: Desired number of particles in the cluster
- `spawn_density::Float`: Must be between 0 and 1, exclusive (usually set to 0.8).
Prevents there being more particles than choices of where to place them on the birth circle, especially towards the start of the simulation.
- `num_steps::Int`: Set to a large number >= 10_000. Number of steps that each particle can walk when generated in parallel.
- `max_radius::Int`: Set to about 150-300. Maximum radius of the birth and death circles.
- `spawn_cap::Int`: Caps the maximum number of particles that can be generated at one time.
- `deposit_prob::Float`: Must be between 0 and 1, inclusive. The probability that a particle attaches to the cluster when there is contact between the two.
"""
function run_parallel_DLA(desired_cluster_mass, spawn_density, num_steps, max_radius, spawn_cap = 300, deposit_prob = 1)

    # Initialize cluster and place a seed at the origin
    cluster = []
    sizehint!(cluster, desired_cluster_mass+1)
    push!(cluster, [0,0])

    # Initialize dynamic variables
    birth_radius = 1
    death_radius = 3
    num_walkers = 3

    while length(cluster) < desired_cluster_mass

        if num_walkers < spawn_cap
            num_walkers = Int(floor(spawn_density * 5.661 * birth_radius))  # 5.66 determined by linear regression of the radii versus
                                                                            #  the number of grid points on a circle of that radius
        else
            num_steps = spawn_cap
        end
        
        # Generate random walks
        trajectories = generate_trajectories(num_walkers, num_steps, birth_radius)

        # Determine fate of each walker
        for trajectory in eachcol(trajectories)
            for step in trajectory

                # Skip walkers beyond a death circle
                if norm(step) > death_radius
                    break
                end

                # Get only the cluster points closest to the walker
                closest_cluster_points = filter(
                    x -> norm(x - step) < 2,  # 2 is arbitrary 
                    cluster
                )

                # If there are no nearby cluster points, skip to next step
                if length(closest_cluster_points) == 0
                    continue
                end

                if is_adjacent_to_cluster(step, closest_cluster_points)
                    if rand() < deposit_prob
                        push!(cluster, step)
                        break
                    end
                end

            end
            
            # End when cluster achieves desired mass
            if length(cluster) >= desired_cluster_mass
                break
            end

            print("\r Mass = $(length(cluster)), r_B = $(birth_radius), r_D = $(death_radius) ")
        end

        # Update birth and death radii
        cluster_radius = maximum(norm.(cluster))
        if death_radius < max_radius
            birth_radius = Int(ceil(cluster_radius + 1))
            death_radius = 2 * birth_radius
        else
            if birth_radius <= max_radius
                birth_radius = Int(ceil(cluster_radius + 1))
            else
                birth_radius = max_radius
            end
            death_radius = max_radius
        end
    end
    return cluster
end


end  # module end
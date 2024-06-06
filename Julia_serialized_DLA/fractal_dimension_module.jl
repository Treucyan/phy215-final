module fractal_dimension


#Function for calculating the Fractal Dimension
"""
# Description
Calculates the number of particles in a fractal per radius r = 5:5:30. 
Returns a 6x2 array where the first column contains the radii and the second columns contains the masses.

## Args
    cluster_aggregate (Matrix{Float64}): The input fractal whose dimension is to be calculated.

## Returns
    [radiusArray mass_per_radius] (Matrix{Float64}): Mass per radius
"""
function mass_per_radius(cluster_aggregate::Matrix{Float64})
    distances = sqrt.(sum(cluster_aggregate.^2,dims = 1))
    radiusArray = 5:5:30
    mass_per_radius = zeros(Float64, 6)
        for k in 1:6
            criterion = x -> x < radiusArray[k]
            mass = count(criterion, distances)
            mass_per_radius[k] = mass
        end
    return [radiusArray mass_per_radius]
end

end

# Fractals from Diffusion Limited Aggregation
## Contibutors:
- GUINTO, Mikael Gian
- JARA, Roy Jr.
- MACATANGAY, David Marick
- OIDEM, James Marwin ``JM”


## Required Packages
The functions in the jl files (`serialized_dla_modules.jl`,`parallelized_DLA.jl`,`time_complexity_module.jl`,`fractal_dimension_module.jl`,) and notebooks (`serialized_dla_notebook.ipynb`, `parallelized_DLA_notebook.ipynb`, and `serialized_dla_fractal_dimensions.ipynb`, and `215_Final_Project...`) rely on the following packages:

- `DelimitedFiles` for accessing the raw data
- `GLM`, `DataFrames`, and `Polynomials` for best fit line
- `Plots` for plotting results
- `LinearAlgebra` for norm calculation in the parallel DLA 
- `Distribution package` for initialization ng random walker in the serial DLA
- `BenchmarkTools` for benchmarking functions

Make sure the following packages are accessible to avoid any errors in the code.


## Main Files Description
#### `serialized_dla_modules.jl` 

This jl file contains all the functions to simulate serial DLA. It has one module:

  - The `Random_walker` module contains the function for creating a random walker, updating position of random walker, generating a random walk, calculating distance of walker from the cluster, calculating cluster particle distance from the origin, and the serial DLA.

#### `parallelized_DLA.jl`
This jl file contains all the functions to simulate parallel DLA. It has one module `Random_walker` with 6 functions:

  - `initialize_randomwalker` for creating a random walker 
  - `walker_update_position` for updating position of random walker
  - `random_walker_generator` for generating a random walk
  - `walker_distance_from_cluster` for calculating distance of walker from the cluster 
  - `cluster_distance_from_origin` for calculating cluster particle distance from the origin 
  - `serialized_dla` for implementation of serial DLA.


#### `time_complexity_module.jl`
This jl file contains all the functions to simulate serial DLA. It has one module `time_complexity` with 3 functions:

  - `run_time_serial` for obtaining elapsed time of runtime of serial DLA
  - `run_time_parallel` for obtaining elapsed time of runtime of parallel DLA
  - `run_time_stats` for obtaining mean and standard deviation of elapsed time for three runs of serial / parallel DLA


#### `fractal_dimension_module.jl` 
This jl file contains one function to calculate the Fractal Dimension of the DLA. It has one module `fractal_dimension` with the function:

  - `mass_per_radius` the function for calculating the Fractal Dimension of the DLA.


#### `215_Final_Project_Guinto_Jara_Macatangay_Oidem.ipynb` 
This is a compilation of all our benchmarking and simulation results, all organized according to the key results in the proposal. This serves as the main notebook for reporting our progress in the final project.


#### `serialized_dla_notebook.jl` 
This is a sandbox notebook, where we explored our functions for serial DLA in greater detail. It also contains the main code for generating the clusters in the `215_Final_Project...` notebook.


#### `parallelized_DLA_notebook.ipynb` 
This is a sandbox notebook, where we explored our functions for parallel DLA in greater detail. It also contains the main code for generating the clusters in the `215_Final_Project...` notebook.

#### `serialized_dla_fractal_dimensions.ipynb`
This is a sandbox notebook, where we explored our function for obtaining fractal dimension of DLA. It contains the log-log plot for different sticking probabilities also found in the `215_Final_Project...` notebook.

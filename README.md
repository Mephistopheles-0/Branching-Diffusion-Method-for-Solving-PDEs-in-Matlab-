# Branching Diffusion Method for PDEs
This repository contains Matlab code for implementing the branching diffusion method for solving partial differential equations (PDEs). The code uses Monte Carlo simulation and the branching process to approximate the solution of PDEs.

## Usage
To use the code, first run the Branching_Matlab function, which sets the parameters for the model and algorithm and calls the MC_BM function to perform the branching method. The BM_Eval function evaluates the branching process for a given set of parameters and the BP function calculates the branching path.

## Parameters
The parameters for the model and algorithm can be adjusted in the Branching_Matlab function. These include:

- T: final time
- t0: initial time
- x0: initial position
- d: dimension
- m: number of particles
- mu: mean of the diffusion process
- sigma: covariance matrix of the diffusion process
- a: vector of weights
- g: function to be evaluated

## Output
The Branching_Matlab function prints the terminal condition, branching method result, estimated standard deviation, estimated L2 approximation error, and elapsed runtime.

## Credits
This code was developed by [Your Name] and is based on the research paper "Branching diffusion processes for solving PDEs" by [Author Name].

## License
This code is released under the MIT License.

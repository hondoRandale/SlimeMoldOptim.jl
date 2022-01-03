
<a id='SimpleFWA.jl'></a>

<a id='SimpleFWA.jl-1'></a>

# SimpleFWA.jl


[![TagBot](https://github.com/hondoRandale/SlimeMoldOptim.jl/actions/workflows/TagBot.yml/badge.svg)](https://github.com/hondoRandale/SlimeMoldOptim.jl/actions/workflows/TagBot.yml)


[![documentation](https://github.com/hondoRandale/SlimeMoldOptim.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/hondoRandale/SlimeMoldOptim.jl/actions/workflows/documentation.yml)


<a id='Introduction'></a>

<a id='Introduction-1'></a>

## Introduction


This solver is an implementation of the slime mould algorithm. It belongs to the family of swarm-intelligence solvers. The method doesn't require first or second order derivatives.    The original paper on the algorithm is called: " Slime mould algorithm: A new method for stochastic optmimization "


___


<a id='calling-convention'></a>

<a id='calling-convention-1'></a>

## calling convention


each Objective function passed to SlimeMoldOptim has to comply with the following    simple parameter convention f( x; kwargs ) where f is the objective    function to be minimized. This convention ensures SlimeMoldOptim can be used with    time-series-problems, classification-problems, regression-problems.    Univariate as well as multivariate target sets are admissible.


___


<a id='performance-considerations'></a>

<a id='performance-considerations-1'></a>

## performance considerations


SlimeMoldOptim uses several cores, to make best use of your sytem, set Threads.nthreads() to the max. number.


___


<a id='example'></a>

<a id='example-1'></a>

## example


```julia
   using SlimeMoldOptim
   using Test
   Easom(x;kwargs) = -cos( x[1] ) * cos( x[2] ) *
                     exp( -( (x[1]-π)^2 + (x[2]-π)^2 ) )
   lower    = [ -10.0f0, -10.0f0 ];
   upper    = [ 10.0f0, 10.0f0 ];
   sma( objFunction ) = slimeMoldAlgo( 30, 1000, ();
                                        lower       = lower,
                                        upper       = upper,
                                        objFunction = objFunction,
                                        maxiter     = 1000,
                                        ϵ_conv      = 1f-20,
                                        z           = 0.03f0 )                            
   solutionFWA = sma( Easom );
   @test isapprox( solutionFWA.x_b[1], π; atol=0.01 )
   @test isapprox( solutionFWA.x_b[2], π; atol=0.01 )                             
```


___


<a id='function-reference'></a>

<a id='function-reference-1'></a>

## function reference

<a id='SlimeMoldOptim.slimeMoldAlgo' href='#SlimeMoldOptim.slimeMoldAlgo'>#</a>
**`SlimeMoldOptim.slimeMoldAlgo`** &mdash; *Function*.



```julia
     slimeMoldAlgo( nMolds::Int,
                    n_vc::Int,
                    kwargs...;
                    lower::Vector{Float32},
                    upper::Vector{Float32},
                    objFunction::Function,
                    maxiter::Int,
                    ϵ_conv::Float32,
                    z::Float32=0.03f0 )


  minimize objective function objFunction, the solution space is limited by lower and upper bound.
  - The optimization algorithm utilized is the slime mold algorithm.
  - The nMolds parameter governs the number of molds being modeled.
  - n_vc is the number of steps learning parameter vb needs to reach 0.0.
  - maxiter denotes the overall number of algo iteraions, it has to be >= n_vc.
  the number of sparks per firework, in remains constant foreach firework.
  - ϵ_A is the smoothing parameter controlling the variance of amplitudes computed foreach fw.
  - C_a ist the upscaling parameter for explosion amplitudes.
  - C_r is the downscaling parameter for explosion amplitudes.
  - maxiter is the maximum number iteraions.
  - ϵ_conv denotes the convergence parameter.
```


<a target='_blank' href='https://github.com/hondoRandale/SlimeMoldOptim.jl/blob/7ad970665f085ed32e501655d4011db12e9d57a1/src/SlimeMoldOptim.jl#L111-L136' class='documenter-source'>source</a><br>


SlimeMolds struct


| Parameter  | Description                         | Type              |
|:---------- |:----------------------------------- |:----------------- |
| S          | each column denotes a mold position | Matrix{Float32}   |
| fitnessS   | fitness of each mold                | Vector{Float32}   |
| x*b*global | best found solution                 | Vector{ Float32 } |
| minFGlobal | best found function value           | Float32           |
| x*w*global | worst found solution                | Vector{ Float32 } |
| maxFGlobal | worst found function value          | Float32           |
| iter       | number of iterations executed       | Int               |


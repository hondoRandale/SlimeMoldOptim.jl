module SlimeMoldOptim

using  InvertedIndices
using  Random
using  Statistics
export slimeMoldAlgo, SlimeMolds

  struct SlimeMolds
    S::Matrix{Float32};
    fitnessS::Vector{Float32};
    x_b_global::Vector{Float32};
    minFGlobal::Float32;
    x_w_global::Vector{Float32};
    maxFGlobal::Float32;
    iter::Int;
  end

  function drawMoldsPositionsBoundary!( ;rng::RandomDevice,
                                        X::Matrix{Float32},
                                        rand_mat::Matrix{Float32},
                                        upper::Vector{Float32},
                                        lower::Vector{Float32} )
    @assert length( upper ) == length( lower )
    @assert all( size( X )  == size( rand_mat ) )

    nFireworks = size( X, 2 );
    d          = length( lower );
    rand_mat  .= Float32.( rand( rng, d, nFireworks ) );
    for iter ∈ 1:1:nFireworks
      @inbounds x_sp     = view( X,        :, iter );
      @inbounds rand_vec = view( rand_mat, :, iter );
      for k ∈ 1:1:d
        @inbounds lu      = lower[k]::Float32;
        @inbounds x_sp[k] = ( rand_vec[k] * ( upper[k] - lu ) + lu )::Float32;
      end
    end
  end

  function smellIndex!( fitnessS::Vector{Float32},
                        W::Vector{Float32},
                        rng::RandomDevice )

    nMolds = length( fitnessS )::Int;
    τ      = median( fitnessS )::Float32;
    bF     = minimum( fitnessS )::Float32;
    wF     = maximum( fitnessS )::Float32;
    denom  = ( bF - wF )::Float32;
    if denom == 0.0f0
      denom = 1f-20;
    end
    rr     = rand( rng, Float32, nMolds );
    @simd for ii ∈ 1:1:nMolds
      if fitnessS[ii] <= τ
        @inbounds W[ii] = 1.0f0 + rr[ii] * log( ( ( bF - fitnessS[ii] ) / denom ) + 1.0f0 );
      else
        @inbounds W[ii] = 1.0f0 - rr[ii] * log( ( ( bF - fitnessS[ii] ) / denom ) + 1.0f0 );
      end
    end
  end

  function updateMoldPositions!( S::Matrix{Float32},
                                 vb::Float32,
                                 vc::Float32,
                                 W::Vector{Float32},
                                 z::Float32,
                                 rng::RandomDevice,
                                 lower::Vector{Float32},
                                 upper::Vector{Float32},
                                 x_b_global::Vector{Float32},
                                 minFGlobal::Float32,
                                 fitnessS::Vector{Float32} )
    nMolds = size( S, 2 );
    d      = size( S, 1 );
    rr     = rand( rng, Float32, nMolds );
    Threads.@threads for ii ∈ 1:1:nMolds
      r = rr[ii];
      if r < z
        ## reinit slime mold
        @inbounds S[:,ii] .= ( rand( rng, Float32, d ) .* ( upper .- lower ) ) .+ lower;
      else
        p = tanh( abs( fitnessS[ii] - minFGlobal  ) );
        if r < p
          ## execute exploitation step
          idx                = collect( 1:1:nMolds )[Not( ii )][rand( 1:d-1, 1 )][1];
          @inbounds S[:,ii] .= x_b_global .+ vb * ( W[idx] * ( S[:,idx] .- S[:,ii] ) )
        else
          ## execute exploration step
          @inbounds S[:,ii] .= vc * S[:,ii];
        end
        ## check if we are out of bounds, if true, reinit
        if any( S[:,ii] .< lower ) || any( S[:,ii] .> upper )
          @inbounds S[:,ii] .= ( rand( rng, Float32, d ) .* ( upper .- lower ) ) .+ lower;
        end
      end
    end
  end

  function evaluateCandidatesParallel( ;kwargs::Tuple,
                                       f::Function,
                                       X::Matrix{Float32},
                                       fitness::Vector{Float32} )
    @assert size( X, 2 ) == size( fitness, 1 )
    n = size( X, 2 )::Int;
    @assert n > 0
    Threads.@threads for ii ∈ 1:1:n
      x            = view( X, :, ii );
      fitness[ii]  = f( x; kwargs=kwargs );
    end
  end

  """


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

        the number of sparks per firework, in remains constant foreach firework. ϵ_A is the smoothing parameter
        controlling the variance of amplitudes computed foreach fw. C_a ist the upscaling parameter for explosion
        amplitudes. C_r is the downscaling parameter for explosion amplitudes. XPrimary is the feature set of the
        primary algorithm to be tuned. yPrimary is the target set of the primary algorithm. maxiter is the maximum
        number iteraions.  ϵ_conv denotes the convergence parameter.
  """
  function slimeMoldAlgo( nMolds::Int,
                          n_vc::Int,
                          kwargs...;
                          lower::Vector{Float32},
                          upper::Vector{Float32},
                          objFunction::Function,
                          maxiter::Int,
                          ϵ_conv::Float32,
                          z::Float32=0.03f0 )
    @assert nMolds >  0
    @assert n_vc   >  0
    @assert n_vc   >= maxiter
    @assert length( lower ) == length( upper )
    @assert all( lower .< upper )
    @assert maxiter >  0
    @assert ϵ_conv  >  0.0f0
    @assert z       >= 0.0f0
    @assert z       <  1.0f0

    d = length( lower );
    @assert d > 0
    rng = RandomDevice();

    δ        = 1.0f0 / Float32( n_vc - 1 );
    vc_full  = collect( 1.0f0:-δ:0.0f0 );
    S        = Matrix{Float32}( undef, d, nMolds );
    rand_mat = Matrix{Float32}( undef, d, nMolds );
    fitnessS = Vector{Float32}( undef, nMolds );
    W        = Vector{Float32}( undef, nMolds );
    x_b      = Vector{Float32}( undef, d );
    p        = -1.0f0;
    drawMoldsPositionsBoundary!( ;rng=rng, X=S, rand_mat=rand_mat, upper=upper, lower=lower );
    maxFGlobal = Float32( objFunction( S[:,1]; kwargs=() ) );
    minFGlobal = Float32( objFunction( S[:,1]; kwargs=() ) );
    x_b_global = zeros( Float32, d );
    x_w_global = zeros( Float32, d );
    iter       = 1;
    while iter <= maxiter
      evaluateCandidatesParallel(;kwargs   = kwargs,
                                 f         = objFunction,
                                 X         = S,
                                 fitness   = fitnessS );
      smellIndex!( fitnessS, W, rng );
      if maximum( fitnessS ) > maxFGlobal
        @inbounds maxFGlobal  = fitnessS[argmax( fitnessS )]::Float32;
        @inbounds x_w_global .= S[:,argmax( fitnessS )];
      end
      if minimum( fitnessS ) < minFGlobal
        @inbounds minFGlobal  = fitnessS[argmin( fitnessS )]::Float32;
        @inbounds x_b_global .= S[:,argmin( fitnessS )];
      end
      a             = Float32( atanh( -( iter / maxiter ) + 1.0f0 ) );
      vb            = ( ( rand( rng, Float32, 1 )[1] * 2.0f0 * a ) - a )::Float32;
      @inbounds vc  = vc_full[iter]::Float32;
      updateMoldPositions!( S, vb, vc, W, z, rng, lower, upper, x_b_global, minFGlobal, fitnessS );
      iter += 1;
    end
    return SlimeMolds( S, fitnessS, x_b_global, minFGlobal, x_w_global, maxFGlobal, iter - 1 )
  end

end

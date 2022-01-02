# MIT License

#Copyright (c) 2021 hondoRandale <jules.rasetaharison@tutanota.com> and contributors

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

include( "benchmarkOptimProblems.jl" )

using SlimeMoldOptim
using Test

@testset "SlimeMoldOptim.jl" begin
    # 2-D objective functions
    @test Beale_2D( [3.0f0,0.5f0]; kwargs=() )                     ≈ 0.0f0
    @test typeof( Beale_2D( [3.0f0,0.5f0]; kwargs=() ) )          == Float32

    @test Himmelblau_2D( [ 3.0f0, 2.0f0 ]; kwargs=() )             ≈ 0.0f0
    @test typeof( Himmelblau_2D( [ 3.0f0, 2.0f0 ]; kwargs=() ) )  == Float32
    @test Himmelblau_2D( [ -2.805118f0,  3.131312f0 ]; kwargs=() ) ≈ 0.0f0 atol=1f-3
    @test Himmelblau_2D( [ -3.779310f0, -3.283186f0 ]; kwargs=() ) ≈ 0.0f0 atol=1f-10
    @test Himmelblau_2D( [ 3.584428f0, -1.848126f0 ]; kwargs=() )  ≈ 0.0f0 atol=1f-5

    @test Levi_no13_2D( [ 1.0f0, 1.0f0 ]; kwargs=() )                      ≈ 0.0f0 atol=1f-15
    @test typeof( Levi_no13_2D( [ 1.0f0, 1.0f0 ]; kwargs=() ) )           == Float32

    @test Ackley_2D( [ 0.0f0, 0.0f0 ]; kwargs=() )                         ≈ 0.0f0
    @test typeof( Ackley_2D( [ 0.0f0, 0.0f0 ]; kwargs=() ) )              == Float32

    @test Booth_2D(  [ 1.0f0, 3.0f0 ]; kwargs=() )                         ≈ 0.0f0
    @test typeof( Booth_2D(  [ 1.0f0, 3.0f0 ]; kwargs=() ) )              == Float32

    @test Rosenbrock_2D( [ 1.0f0, 1.0f0 ]; kwargs=() )                     ≈ 0.0f0
    @test typeof( Rosenbrock_2D( [ 1.0f0, 1.0f0 ]; kwargs=() ) )          == Float32

    @test Matyas_2D( [ 0.0f0, 0.0f0 ]; kwargs=() )                         ≈ 0.0f0
    @test typeof( Matyas_2D( [ 0.0f0, 0.0f0 ]; kwargs=() ) )              == Float32

    @test Three_hump_camel_2D( [ 0.0f0, 0.0f0 ]; kwargs=() )               ≈ 0.0f0
    @test typeof( Three_hump_camel_2D( [ 0.0f0, 0.0f0 ]; kwargs=() ) )    == Float32

    @test Easom_2D( [ Float32(π), Float32(π) ]; kwargs=() )                ≈ -1.0f0
    @test typeof( Easom_2D( [ Float32(π), Float32(π) ]; kwargs=() ) )     == Float32

    @test Sphere_2D( [0.0f0, 0.0f0]; kwargs=() )                           ≈ 0.0f0
    @test typeof( Sphere_2D( [0.0f0, 0.0f0]; kwargs=() ) )                == Float32

    @test Schaffer_function_no2_2D( [0.0f0, 0.0f0]; kwargs=() )            ≈ 0.0f0
    @test typeof( Schaffer_function_no2_2D( [0.0f0, 0.0f0]; kwargs=() ) ) == Float32

    @test Bukin_function_no6_2D( [-10.0f0, 1.0f0 ]; kwargs=() )            ≈ 0.0f0
    @test typeof( Bukin_function_no6_2D( [-10.0f0, 1.0f0 ]; kwargs=() ) ) == Float32

    # N-D objective functions
    @test Michalewicz_ND( [ 2.20f0, 1.57f0 ]; kwargs=() )                          ≈ -1.8081267f0
    @test typeof( Michalewicz_ND( [ 2.20f0, 1.57f0 ]; kwargs=() ) )               == Float32

    @test Zakharov_ND( [ 0.0f0, 0.0f0 ]; kwargs=() )                               ≈ 0.0f0
    @test typeof( Zakharov_ND( [ 0.0f0, 0.0f0 ]; kwargs=() ) )                    == Float32

    @test AxisParallelHyperEll_ND( [ 0.0f0, 0.0f0 ]; kwargs=() )                   ≈ 0.0f0
    @test typeof( AxisParallelHyperEll_ND( [ 0.0f0, 0.0f0 ]; kwargs=() ) )        == Float32

    @test Sphere_ND( [ 0.0f0, 0.0f0 ]; kwargs=() )                                 ≈ 0.0f0
    @test typeof( Sphere_ND( [ 0.0f0, 0.0f0 ]; kwargs=() ) )                      == Float32

    @test Rastrigin_ND( [ 0.0f0, 0.0f0 ]; kwargs=() )                              ≈ 0.0f0
    @test typeof( Rastrigin_ND( [ 0.0f0, 0.0f0 ]; kwargs=() ) )                   == Float32

    @test Styblinski_Tang_ND( [ -2.903534f0, -2.903534f0 ]; kwargs=() )            ≈ -78.33233f0
    @test typeof( Styblinski_Tang_ND( [ -2.903534f0, -2.903534f0 ]; kwargs=() ) ) == Float32

    @test Exponential_ND( [ 0.0f0, 0.0f0 ]; kwargs=() )                            ≈ -1.0f0
    @test typeof( Exponential_ND( [ 0.0f0, 0.0f0 ]; kwargs=() ) )                 == Float32

    @test Griewank_ND( [ 0.0f0, 0.0f0 ]; kwargs=() )                               ≈ 0.0f0
    @test typeof( Griewank_ND( [ 0.0f0, 0.0f0 ]; kwargs=() ) )                    == Float32

    lower    = [ -10.0f0, -10.0f0 ];
    upper    = [ 10.0f0, 10.0f0 ];

    sma( objFunction ) = slimeMoldAlgo( 30, 1000, ();
                                        lower       = lower,
                                        upper       = upper,
                                        objFunction = objFunction,
                                        maxiter     = 1000,
                                        ϵ_conv      = 1f-20,
                                        z           = 0.03f0 )


   @test sum( [ sma( fObj ).minFGlobal for fObj ∈ obj_functions_2D ] ) < 0.0f0
   @test sum( [ sma( fObj ).minFGlobal for fObj ∈ obj_functions_Nd ] ) < 750.0f0
end

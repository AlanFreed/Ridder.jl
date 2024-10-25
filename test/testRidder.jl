#=
Created on Wed 12 Jun 2024
Updated on Fri 25 Oct 2024

This test program verifies the root finding algorithm of Ridder implemented in
Ridder.jl. This is done by finding a root for the cosine function.
=#

#------------------------------------------------------------------------------

module testRidder

using
    Ridder

export
    run

#=
-------------------------------------------------------------------------------
=#

function run()
    function fn(x::Float64)::Float64
        return cos(x)
    end # fn

    xL = 1.0
    xR = 2.0

    root = findRoot(xL, xR, fn)

    error = abs(root - π/2)
    println("The error in finding the root for cos(x) at π/2 was ", error)
end # run

end # testRidder

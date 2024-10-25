#=
Created on Tue 11 Jun 2024
Updated on Thr 24 Oct 2024
Translated from Python (code dated 09/24/2019).
=#

"""
# Ridder

Applys Ridder's method to find a root of a function, i.e., find a *x* such that
```julia
f(x) = 0.
```

## Function

To find a root to a default error tolerance of 1.0E-9, call:
```julia
root = findRoot(xleft, xright, f)
```
To find a root to a specified error tolerance, call:
```julia
root = findRoot(xleft, xright, f, tol)
```
where

    xleft   is the lower boundary of the search window
    xright  is the upper boundary of the search window
    f       is the function whose root is being sought
    tol     is the error tolerance in input variable x

**Note**: *xleft* and *xright* must be of opposite sign to ensure a root exits between them.
"""
module Ridder

export
    findRoot

#=
-------------------------------------------------------------------------------
=#

N = 30  # Maximum number of iterations allowed.

function findRoot(xleft::Float64, 
                  xright::Float64, 
                  f::Function, 
                  tol::Float64=1.0E-9)::Float64

    # analyze the original interval
    xl = xleft
    fl = convert(Float64, f(xl))
    if abs(fl) < tol
        return xl
    end
    xr = xright
    fr = convert(Float64, f(xr))
    if abs(fr) < tol
        return xr
    end
    if sign(fl) == sign(fr)
        error("Error: root is not bracketed.")
    end
    
    # Perform Ridder's algorithm for root finding.
    n    = 1
    ξₙ₋₁ = 0.5(xl + xr)
    while n < N
        xm = 0.5(xl + xr)    # x at the midpoint
        fm = convert(Float64, f(xm))
        s = sqrt(fm^2 - fl*fr)
        if s == 0.0
            error("Ridder's method became singular.")
        end
        dx = (xm - xl) * fm / s
        if (fl - fr) < 0.0
            dx = -dx
        end
        ξₙ = xm + dx        # an estimate for the root
        fξ = convert(Float64, f(ξₙ))
        # Test for convergence.
        if abs(fξ) < tol || abs(ξₙ - ξₙ₋₁) < tol*max(abs(ξₙ), 1.0)
            return ξₙ
        end
        
        # Re-bracket the root as tightly as possible.
        ξₙ₋₁ = ξₙ
        if sign(fm) == sign(fξ)
            if sign(fl) ≠ sign(fξ)
                xr = ξₙ
                fr = fξ
            else
                xl = ξₙ
                fl = fξ
            end
        else
            xl = xm
            fl = fm
            xr = ξₙ
            fr = fξ
        end
        
        # Increment the counter.
        n += 1
    end

    println("Warning: iterations exceeded a preset maximum in findRoot.")
    return ξ
end # findRoot

end # Ridder


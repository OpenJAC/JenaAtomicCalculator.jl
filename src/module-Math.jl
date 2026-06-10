
"""
`module JAC.Math`  
... a submodel of JAC that contains various (adopted) procedures from different libraries and sources.
"""
module Math


using GSL
using JenaAtomicCalculator, ..Radial


"""
`Math.finiteDifferenceWeights(position, npoints; order=1)`  
    ... computes the weights for a n-point finite-difference approximation. position=0 corresponds to the centered stencil, 
        while positive and negative values refer to the non-symmetric stencils, respecitvely. The order specifies for which 
        derivative the coefficients are computed, the default ist the first derivative.
"""
function finiteDifferenceWeights(position, npoints; order=1)

    if !isodd(npoints)
        error("Number of points must be odd")
    end

    x = collect(-Int64((npoints-1)/2):Int64((npoints-1)/2))
    c = Array{Float64}(undef, order+1, npoints)

    c1 = 1
    c4 = x[1] - position
    c[1, 1] = 1

    for i = 2:npoints
        mn = min(i, order+1)
        c2 = 1
        c5 = c4
        c4 = x[i] - position
        for j = 1:i-1
            c3 = x[i] - x[j]
            c2 = c2 * c3
            if j == i - 1
                c[2:mn, i] = c1 / c2 * ( collect(1:mn -1 ) .* c[1:mn-1, i-1] - c5 * c[2:mn, i-1])
                c[1, i] = -c1 * c5 * c[1, i-1]/c2
            end
            c[2:mn, j] = 1 / c3 * (c4 * c[2:mn, j] - collect(1:mn-1) .* c[1:mn-1, j])
            c[1, j] = c4 * c[1, j] / c3
        end
        c1 = c2
    end

    return( c )
end


"""
`Math.derivative(input::Array{Float64}, i::Int64)`  
    ... computes the derivative of a given input array at position i; a Float64 is returned. An equally spaced grid is 
        assumed, the result must be multiplied by the inverse lattice spacing h manually. If an exponential grid is 
        used, the transformation must be applied too.
"""
function derivative(input::Array{Float64}, i::Int64) :: Float64
    
    n = JenaAtomicCalculator.Defaults.FINITE_DIFFERENCE_NPOINTS
    result ::Float64 = 0.
    
    if i <= n
    
        for j = 1:2*n + 1
            result += JenaAtomicCalculator.Defaults.weights[i, j] * input[j]
        end
        return result
    end
    
    if size(input, 1) - i < n
    
        for j = 1:2*n + 1
            result += JenaAtomicCalculator.Defaults.weights[2*n + 1 + i - size(input, 1), j] * input[size(input, 1) - 2*n - 1 + j]
        end
    
        return result
    end
    
    for j = i - n:i+n
        result += JenaAtomicCalculatorDefaults.weights[n + 1, j - i + n + 1] * input[j]
    end
    
    return result
end


"""
`Math.integrateFitTransform(f::Function, mtp::Int64, grid::Radial.Grid)`  
    ... to integrate a (point-wise) given function over mtp grid points on an exponential grid. The transformation of the 
        exponential grid is automatically performed such that no modifications to the integrand are necessary; a value::Float64 
        is returned.
"""
function integrateFitTransform(f::Function, mtp::Int64, grid::Radial.Grid)
    n      = size(JenaAtomicCalculator.Defaults.newtonCotesCoefficients, 1)
    newtonCotesCoefficientsOnGrid = copy(JenaAtomicCalculator.Defaults.newtonCotesCoefficients) * grid.h
    i0     = 2
    gamma  = 0.
    result = 0.

    while abs(f(i0)) == 0 || f(i0)/f(i0+1) < 0
        i0 += 1
    end
    gamma = log(f(i0)/f(i0+1)) / log(grid.r[i0]/grid.r[i0+1])
    result = 1/(gamma + 1) * grid.r[i0] * f(i0)

    result += f(i0) * grid.rp[i0] * newtonCotesCoefficientsOnGrid[n]
    for i = i0+1:mtp
        result += newtonCotesCoefficientsOnGrid[(i - i0) % (n - 1) + 1] * f(i) * grid.rp[i]
    end

    if (mtp - i0) % (n - 1) + 1 == 1
        result -= newtonCotesCoefficientsOnGrid[n] * f(mtp) * grid.rp[mtp]
    end

    return( result )
end


"""
`Math.integrateFit(f::Function, mtp::Int64, grid::Radial.Grid)`  
    ... to integrate a (point-wise) given function over mtp grid points on an exponential grid. This function assumes that the 
        ingegrand was already transformed to the exponential grid by multiplying it with the derivative of the transformation 
        function; a value::Float64 is returned.
"""
function integrateFit(f::Function, mtp::Int64, grid::Radial.Grid)
    n      = size(JenaAtomicCalculator.Defaults.newtonCotesCoefficients, 1)
    newtonCotesCoefficientsOnGrid = copy(JenaAtomicCalculator.Defaults.newtonCotesCoefficients) * grid.h

    i0     = 2
    gamma  = 0.
    result = 0.

    while abs(f(i0)) == 0 || f(i0)/f(i0+1) < 0
        i0 += 1
    end
    gamma = log(f(i0)/f(i0+1) * grid.rp[i0+1] / grid.rp[i0]) / log(grid.r[i0]/grid.r[i0+1])
    result = 1/(gamma + 1) * grid.r[i0] * f(i0) / grid.rp[i0]

    result += f(i0) * newtonCotesCoefficientsOnGrid[n]
    for i = i0+1:mtp
        result += newtonCotesCoefficientsOnGrid[(i - i0) % (n - 1) + 1] * f(i)
    end

    if (mtp - i0) % (n - 1) + 1 == 1
        result -= newtonCotesCoefficientsOnGrid[n] * f(mtp)
    end

    return result
end



"""
`Math.integrateTransform(f::Function, n0::Int64, mtp::Int64, grid::Radial.Grid)`  
    ... to integrate a (point-wise) given function over mtp grid points ...; a value::Float64 is returned.
"""
function integrateTransform(f::Function, n0::Int64, mtp::Int64, grid::Radial.Grid)
    newtonCotesCoefficientsOnGrid             = Array{Float64}([2 * 7, 32, 12, 32, 7]) * 2 / 45  # * grid.h
    ## newtonCotesCoefficientsOnGrid            = Array{Float64}([2 * 41, 216, 27, 272, 27, 216, 41]) * 1 / 140
    ## newtonCotesCoefficientsOnGrid             = Array{Float64}([2 * 989, 5888, -928, 10496, -4540, 10496, -928, 5888, 989]) * 4 / 14175
    newtonCotesCoefficientsOnGrid *= grid.h
    n      = size(newtonCotesCoefficientsOnGrid, 1)

    i0     = n0
    gamma  = 0.
    result = 0.

    result += f(i0) * grid.rp[i0] * newtonCotesCoefficientsOnGrid[n]
    for i = i0+1:mtp
        result += newtonCotesCoefficientsOnGrid[(i - i0) % (n - 1) + 1] * f(i) * grid.rp[i]
    end

    if (mtp - i0) % (n - 1) + 1 == 1
        result -= newtonCotesCoefficientsOnGrid[n] * f(mtp) * grid.rp[mtp]
    end

    return( result )
end


"""
`Math.polylogExp(x::Float64, s::Int64)`  ... to evaluate the polylog(s, -exp(x)) function.
"""
function polylogExp(x::Float64, s::Int64)
    return -sf_fermi_dirac_int(s-1, x)
end

end # module


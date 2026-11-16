#
# probeDifferentiability.jl -- priority item 25: is JAC accidentally differentiable, and if not, by what?
#
#   Run from the JAC root:   julia --project=. tools/probeDifferentiability.jl
#
# A dual number x + eps*x', with eps^2 = 0, makes ordinary arithmetic reproduce the chain rule: a function
# written generically in its argument's type returns the derivative for free when that argument is seeded
# with one. DFTK.jl became differentiable largely by being type-generic throughout; what breaks a port of
# that idea is mundane -- one `::Float64` argument annotation, one `::Float64` struct field, one call into
# compiled BLAS/LAPACK code -- and where it breaks matters more than whether it breaks, since a mundane
# annotation is cheap to relax and a BLAS call is not.
#
# So this probe does not ask "does it run" but "does the number that comes back match a derivative this
# script can verify independently" -- and where it fails, it names the exact call that stopped it, rather
# than guessing.
#
# ONE PATH, as the item scoped it: the one-electron hydrogenic spectrum in a local (point-nucleus) Dirac
# potential, `Basics.computeDiracEnergy`. Its own closed-form Sommerfeld expression is differentiated BY
# HAND below as the oracle -- not by finite differences, which is only as good as its own step size, and
# not by the same dual machinery under test, which would be circular.
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

# ===== 1. A minimal dual number -- just enough arithmetic to carry a derivative through this one formula

struct Dual{T<:Real} <: Real
    val::T
    eps::T
end

Base.promote_rule(::Type{Dual{T}}, ::Type{S}) where {T<:Real,S<:Real} = Dual{promote_type(T,S)}
Base.convert(::Type{Dual{T}}, x::Real) where T = Dual(T(x), zero(T))
Base.convert(::Type{Dual{T}}, x::Dual) where T = Dual(T(x.val), T(x.eps))

Base.:+(a::Dual, b::Dual) = Dual(a.val + b.val, a.eps + b.eps)
Base.:-(a::Dual, b::Dual) = Dual(a.val - b.val, a.eps - b.eps)
Base.:-(a::Dual)          = Dual(-a.val, -a.eps)
Base.:*(a::Dual, b::Dual) = Dual(a.val * b.val, a.eps*b.val + a.val*b.eps)
Base.:/(a::Dual, b::Dual) = Dual(a.val / b.val, (a.eps*b.val - a.val*b.eps) / b.val^2)
Base.sqrt(a::Dual)        = (s = sqrt(a.val); Dual(s, a.eps / (2s)))
Base.literal_pow(::typeof(^), a::Dual, ::Val{2}) = Dual(a.val^2, 2*a.val*a.eps)
Base.isless(a::Dual, b::Real) = isless(a.val, b)
Base.isless(a::Real, b::Dual) = isless(a, b.val)
Base.isless(a::Dual, b::Dual) = isless(a.val, b.val)
Base.:<(a::Dual, b::Real)  = a.val < b
Base.:<(a::Real, b::Dual)  = a < b.val
Base.:<(a::Dual, b::Dual)  = a.val < b.val
Base.:<=(a::Dual, b::Real) = a.val <= b
Base.:<=(a::Real, b::Dual) = a <= b.val
Base.:<=(a::Dual, b::Dual) = a.val <= b.val

# Deliberately NO `convert(::Type{Float64}, ::Dual)`: a lossy conversion back to plain Float64 is exactly
# the silent failure this probe exists to catch, so a struct field typed `::Float64` must raise, not drop
# the eps slot quietly.

# ===== 2. Self-test of the dual type alone, against ordinary calculus, before trusting it on anything else

function selfTest()
    x  = Dual(1.7, 1.0)
    y  = sqrt(x^2 + 3*x)
    dy = (2*1.7 + 3) / (2*sqrt(1.7^2 + 3*1.7))

    @printf("self-test   d/dx sqrt(x^2+3x) at x=1.7:   dual = %.15f   calculus = %.15f   diff = %.2e\n",
            y.eps, dy, abs(y.eps - dy))
end

# ===== 3. The SAME closed form `Basics.computeDiracEnergy` uses, reproduced here with Z::Real instead of
#          Z::Float64 -- the one signature change that would let the real function accept a dual, nothing
#          else touched. NOT a src/ edit; this file never writes to src/.

function diracEnergyGeneric(sh, Z::Real; mass::Real=1.0)
    jPlusHalf = (Basics.subshell_2j(sh) + 1) / 2;   nr = sh.n - jPlusHalf
    alpha = Defaults.getDefaults("alpha")
    wa = sqrt(jPlusHalf^2 - Z^2 * alpha^2)
    wa = sqrt(1.0 + Z^2 * alpha^2 / (nr + wa)^2)
    wa = mass * Defaults.getDefaults("speed of light: c")^2 * (1/wa - 1.0)

    return( wa )
end

# ===== 4. The independent oracle: dE/dZ of the SAME closed form, differentiated BY HAND -- owes nothing
#          to the dual machinery above and nothing to a finite-difference step size

function diracEnergyDerivativeAnalytic(sh, Z::Real; mass::Real=1.0)
    jPlusHalf = (Basics.subshell_2j(sh) + 1) / 2;   nr = sh.n - jPlusHalf
    alpha = Defaults.getDefaults("alpha");          c1 = alpha^2
    wa1   = sqrt(jPlusHalf^2 - Z^2*c1)
    D     = nr + wa1
    wa2   = sqrt(1.0 + Z^2*c1/D^2)
    dwa1  = -Z*c1/wa1
    dg    = 2*Z*c1*(D - Z*dwa1) / D^3
    de    = -mass * Defaults.getDefaults("speed of light: c")^2 * dg / (2*wa2^3)

    return( de )
end

# ===== 5. Run it

println("JAC differentiability probe (priority item 25) -- one path, three questions.\n")
selfTest()
println()

sh = Subshell("1s_1/2");   Z = 26.0

println("Q1: does the REAL `Basics.computeDiracEnergy`, called exactly as src/ calls it, accept a dual?")
try
    e = Basics.computeDiracEnergy(sh, Dual(Z, 1.0))
    @printf("    SUCCEEDED, dE/dZ = %.10f\n", e.eps)
catch err
    println("    FAILED at dispatch -- its signature is `Z::Float64`, a concrete annotation a dual cannot match:")
    println("    ", sprint(showerror, err))
end
println()

println("Q2: if that ONE annotation were relaxed to `Z::Real` and nothing else touched, would the SAME")
println("    formula differentiate correctly?")
generic  = diracEnergyGeneric(sh, Dual(Z, 1.0))
analytic = diracEnergyDerivativeAnalytic(sh, Z)
@printf("    dual-seeded generic copy:        dE/dZ = %.10f\n", generic.eps)
@printf("    independent hand-derived oracle:  dE/dZ = %.10f\n", analytic)
@printf("    agreement: %.2e  --  %s\n", abs(generic.eps - analytic),
        abs(generic.eps - analytic) < 1e-9 ? "YES, to 9 digits" : "NO")
println()

println("Q3: one level out -- `Nuclear.Model(Z::Real)` is ALREADY written generically. Does a dual survive")
println("    its own body, which is JAC code this script does not control or relax?")
try
    nm = Nuclear.Model(Dual(Z, 1.0))
    @printf("    SUCCEEDED, which was not expected -- nm.Z = %s\n", string(nm.Z))
catch err
    println("    FAILED inside the function body, despite the ::Real signature accepting the dual:")
    println("    ", sprint(showerror, err))
end

#
# probe-fineStructureSpectator.jl   --   the reproducer for priority item 27.
#
# A CONFIGURATION INTERACTION THAT SHARES AN OPEN SPECTATOR SUBSHELL WITH j > 1/2 INVERTS THE FINE STRUCTURE OF
# THE TERM IT ACTS ON.  Adding ONE configuration to a {3s^2, 3s3p, 3p^2} CI either leaves the Mg-like 3s3p 3P*
# term alone or turns it inside out, and which of the two happens is decided by the SPECTATOR:
#
#   * a single replacement whose spectator is 3s  (j = 1/2 only)          -- correct
#   * a DOUBLE replacement, sharing no open subshell                      -- correct
#   * a single replacement whose spectator is 3p  (j = 1/2 and 3/2)       -- INVERTED
#
# THE TEST NEEDS NO EXPERIMENTAL DATA, which is what makes it a proof rather than a comparison.  This interaction
# is diagonal in L and S, so in the LS limit it shifts every J of the term by the SAME amount and cannot change
# the interval at all; the Lande interval rule then fixes 1-2 / 0-1 = 2 for a normal 3P.  Z = 12 is run because
# intermediate coupling there is far too weak to move that ratio -- Mg I's real intervals are 20.1 and 40.7 cm^-1.
# Delta-l is NOT the trigger and must not be reported as one: 3s4p (dl = 0) is correct, 3p4s (dl = 0) is inverted.
#
# The second table follows the spurious shift along the sequence, to show it is systematic rather than a low-Z
# accident: it grows roughly linearly with Z and at Z = 26 takes Fe XV from +6 % agreement with NIST to the
# wrong sign.  Parity is the built-in control -- 3s3d and 3d^2 are EVEN, cannot touch the odd 3s3p, and must
# move nothing.
#
using JenaAtomicCalculator, Printf

c(s)  = Configuration("1s^2 2s^2 2p^6 " * s)
base  = [c("3s^2"), c("3s 3p"), c("3p^2")]

"""
`intervals(refs::Array{Configuration,1}, Z::Float64)`
    ... runs an AL mean field plus CI on the given reference configurations and returns the two `3s3p 3P*`
        fine-structure intervals, `(E(J=1)-E(J=0), E(J=2)-E(J=1))` in cm^-1, as a Tuple{Float64,Float64}.
        Levels are identified by their DOMINANT configuration and not by energy order, since a correlation
        level can lie between the ones wanted.
"""
function intervals(refs::Array{Configuration,1}, Z::Float64)
    grid = Basics.recommendedGrid(refs, Nuclear.Model(Z); rbox = Z < 16 ? 25.0 : 10.0)
    set  = AsfSettings(AsfSettings(); scField=Basics.ALField(), eeInteraction=CoulombInteraction(),
                                      eeInteractionCI=CoulombInteraction(), gridStopper=false)
    mp   = SelfConsistent.performSCF(refs, Nuclear.Model(Z), grid, set; printout=false)
    bs   = mp.levels[1].basis
    cof  = [ Basics.extractConfiguration(Basics.FromBasis(), bs, cc)  for cc in bs.csfs ]
    i2   = findall(x -> x == c("3s^2"), cof);    ip = findall(x -> x == c("3s 3p"), cof)
    wt(l, ix) = sum( l.mc[r]^2  for r in ix; init=0.0 )
    sy(J, p)  = [ l  for l in mp.levels  if l.J == AngularJ64(J) && l.parity == p ]
    pk(J, p, ix) = ( ls = sy(J, p);  isempty(ls) ? nothing : ls[argmax([wt(l, ix) for l in ls])] )
    g  = pk(0, Basics.plus, i2);   e0 = pk(0, Basics.minus, ip);   e2 = pk(2, Basics.minus, ip)
    l1 = sort( [ l for l in sy(1, Basics.minus) if wt(l, ip) >= 0.20 ], by = l -> l.energy )
    (isnothing(g) || isnothing(e0) || isnothing(e2) || isempty(l1))  &&  return( (NaN, NaN) )
    k(x) = Defaults.convertUnits("energy: from atomic to Kayser", x - g.energy)

    return( (k(l1[1].energy) - k(e0.energy), k(e2.energy) - k(l1[1].energy)) )
end

sets = [ ("(nothing added)",              base,                          "-",  "correct"),
         ("+ 3s4p   3p->4p, spectator 3s", vcat(base, [c("3s 4p")]),     "0",  "correct"),
         ("+ 3s4f   3p->4f, spectator 3s", vcat(base, [c("3s 4f")]),     "2",  "correct"),
         ("+ 3s5p   3p->5p, spectator 3s", vcat(base, [c("3s 5p")]),     "1",  "correct"),
         ("+ 4s4p   BOTH replaced",        vcat(base, [c("4s 4p")]),     "-",  "correct"),
         ("+ 3s3d   even, cannot couple",  vcat(base, [c("3s 3d")]),     "-",  "control"),
         ("+ 3d^2   even, cannot couple",  vcat(base, [c("3d^2")]),      "-",  "control"),
         ("+ 3p3d   3s->3d, spectator 3p", vcat(base, [c("3p 3d")]),     "2",  "INVERTED"),
         ("+ 3p4d   3s->4d, spectator 3p", vcat(base, [c("3p 4d")]),     "2",  "INVERTED"),
         ("+ 3p4s   3s->4s, spectator 3p", vcat(base, [c("3p 4s")]),     "0",  "INVERTED"),
         ("+ 3p5s   3s->5s, spectator 3p", vcat(base, [c("3p 5s")]),     "1",  "INVERTED") ]

println("\n", "="^96)
println("Z = 12, near-pure LS coupling: the 3P ratio 1-2 / 0-1 MUST stay near 2 for every row")
println("="^96)
@printf("%-32s %4s %11s %11s %8s   %s\n", "added configuration", "dl", "0-1", "1-2", "ratio", "expected")
for  (tag, rf, dl, verdict)  in  sets
    (a, b) = intervals(rf, 12.0)
    @printf("%-32s %4s %11.1f %11.1f %8.2f   %s\n", tag, dl, a, b, b/a, verdict);   flush(stdout)
end

println("\n", "="^96)
println("The spurious shift along the Mg-like sequence: complete n=3 complex against {3s^2, 3s3p, 3p^2}")
println("="^96)
full = vcat(base, [c("3s 3d"), c("3p 3d"), c("3d^2")])
@printf("%5s  %11s %11s   %11s %11s   %12s\n", "Z", "base 0-1", "base 1-2", "complex 0-1", "complex 1-2", "shift in 0-1")
for  Z  in  [12.0, 14.0, 18.0, 22.0, 26.0]
    (a1, a2) = intervals(base, Z);    (b1, b2) = intervals(full, Z)
    @printf("%5.0f  %11.1f %11.1f   %11.1f %11.1f   %+12.1f\n", Z, a1, a2, b1, b2, b1 - a1);   flush(stdout)
end
println("\nNIST for scale:  Mg I  3P*_0-1 = 20.1, 1-2 = 40.7 cm^-1;   Fe XV  0-1 = 5818, 1-2 = 14160 cm^-1.")

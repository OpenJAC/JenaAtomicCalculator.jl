
println("Dv) Two-photon IONIZATION of a one-electron system: generalized cross sections, the intermediate sum,")
println("    and the resonance guard.")

using Printf

# WRITTEN 25-Aug-2026, first implementation of module-MultiPhotonIonization.jl.  What stood there before computed
# nothing; the one file with content produced partial amplitudes and stopped short of any cross section.  All of it
# was discarded and written again.
#
# WHERE THIS PROCESS BELONGS, since three modules could plausibly claim it.  MultiPhotonTransition is exclusively
# BOUND-BOUND and says so itself: "As soon as the final state carries a free electron the process belongs to
# MultiPhotonIonization."  PhotonScattering is for processes with a photon on BOTH sides.  Here two photons go in
# and an electron comes out, which is neither.
#
# THE ESCAPING ELECTRON CHANGES TWO THINGS relative to two-photon absorption between bound levels:
#   + the final state is a continuum partial wave, so the final channels are distinct states and add
#     INCOHERENTLY, while the intermediate symmetries that reach one such channel still add coherently;
#   + the electron energy is NOT free.  At one photon energy, conservation fixes eps = 2*omega - I_P, so there is
#     no energy-sharing scan of the kind a two-photon EMISSION line needs; omega is the variable to scan.
#
# UNITS, WHICH ARE NOT AN AREA.  A two-photon process needs two photons at once, so its rate goes as the SQUARE of
# the photon flux, W = sigma^(2) F^2 with F in photons cm^-2 s^-1.  Dimensionally sigma^(2) therefore carries
# cm^4 s, and is quoted in GM (1 GM = 1e-50 cm^4 s) after Maria Goeppert-Mayer, who predicted two-photon
# absorption in her 1931 thesis.  One atomic unit, a0^4 t_au, comes out at 1.8968 GM -- the tables print that line
# themselves, together with the single-beam convention, which is pinned as in MultiPhotonTransition by requiring
# the monochromatic result to agree with the bichromatic one as omega_1 -> omega_2.


if  true
    # Last visit:      25-Aug-2026
    # Last successful:
    #
    # NO DATE, DELIBERATELY, and the reason is in the branch itself: nothing here is checked against anything
    #   outside this module.  The magnitudes are plausible and the two gauges agree, but "it ran and the number
    #   looks about right" is precisely what Rule 7 keeps a date away from.  Branches b and c carry checks that
    #   could actually fail; this one does not, and it stays blank until an absolute reference is found for it.
    #
    # Branch a: THE REFERENCE COMPUTATION -- two-photon ionization of hydrogen from 1s, in the clean window between
    #   the two-photon threshold at omega = 0.25 a.u. and the 1s -> 2p resonance at 0.375 a.u.  Both gauges are
    #   computed, and their agreement is the only internal check this module has.
    #
    # REPORT (25-Aug-2026): with 34 bound intermediates (n <= 8) and a 90 a.u. box, sigma^(2) in GM:
    #        omega      eps(e-)     linear Coul   linear Bab   ratio     unpolarized Coul / Bab
    #        0.30       0.100       6.52215e-01   7.46943e-01   0.873     5.64360e-01  6.39188e-01
    #        0.32       0.140       6.45269e-01   6.17444e-01   1.045     5.57974e-01  5.28415e-01
    #        0.35       0.200       1.30345e+00   1.27714e+00   1.021     1.12613e+00  1.10981e+00
    #   Magnitudes of order 1 GM are right for an atomic two-photon cross section, and the two gauges agree to
    #   2-13%.  THE ABSOLUTE SCALE IS NOT VERIFIED against any published value, so these numbers are not to be
    #   quoted; what is dated here is that the machinery runs, that the gauges are close, and that the shape and
    #   size are not absurd.
    #
    #   THE BOX FOR THE EJECTED ELECTRON IS THE DOMINANT PRACTICAL CONTROL, and it is easy to get wrong in a way
    #   that looks like physics.  Near threshold the electron is slow and its wavelength long: at omega = 0.26 the
    #   excess energy is 0.02 a.u., a de Broglie wavelength of 31 a.u., and Continuum refuses outright with
    #   "enlarge box-size" -- which is the good case.  The bad case is a box that is merely too small: at
    #   rbox = 60 a.u. the same computation ran without complaint and returned gauge ratios as far from unity as
    #   0.2.  Branch b shows what that did to a convergence study.
    setDefaults("method: continuum, Galerkin");   setDefaults("method: normalization, pure sine")
    setDefaults("print summary: open", "zzz-MultiPhotonIonization-Dv-reference.sum")

    Z = 1.0
    grid = Radial.Grid(Radial.Grid(false); rnt=2.0e-6, h=5.0e-2, hp=1.0e-2, rbox=90.0)
    nm   = Nuclear.Model(Z, Nuclear.PointNucleus(), 0., 0., AngularJ64(0), 0., 0., 0.)
    prim = Bsplines.generatePrimitives(grid)
    shs  = Subshell[]
    for n = 1:8
        push!(shs, Subshell(n,-1))
        if n >= 2   push!(shs, Subshell(n,1));   push!(shs, Subshell(n,-2))   end
        if n >= 3   push!(shs, Subshell(n,2));   push!(shs, Subshell(n,-3))   end
    end
    orbs = Bsplines.generateOrbitalsHydrogenic(shs, nm, prim; printout=false)
    pot  = Nuclear.nuclearPotential(nm, grid)
    MultiPhotonIonization.computeLines2pOneElectron(Subshell("1s_1/2"), [0.30, 0.32, 0.35], [E1], orbs, pot)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      25-Aug-2026
    # Last successful: 25-Aug-2026
    #
    # Branch b: DOES THE INTERMEDIATE SUM CONVERGE?  A second-order amplitude sums over the WHOLE spectrum, bound
    #   and free.  Only the bound part is available here, so the question is how much that costs -- and the gauge
    #   ratio answers it without any external reference, since length and velocity forms agree exactly only for a
    #   COMPLETE intermediate set.
    #
    # REPORT (25-Aug-2026), hydrogen 1s at omega = 0.30 a.u., 90 a.u. box, sigma^(2) linear in ATOMIC units:
    #        n <= 4    14 orbitals   2.91666e-01   3.67928e-01   ratio 0.7927
    #        n <= 6    24 orbitals   3.18452e-01   3.96701e-01   ratio 0.8028
    #        n <= 8    34 orbitals   3.43852e-01   3.93793e-01   ratio 0.8732
    #        n <= 10   44 orbitals   3.85381e-01   4.24246e-01   ratio 0.9084
    #   THE SUM CONVERGES, and slowly: the ratio climbs monotonically towards 1 but is still 9% away with 44
    #   orbitals.  That residual is the missing CONTINUUM part of the intermediate sum, which no number of bound
    #   states can supply, and it is the reason a Green's-function expansion rather than an orbital list is the
    #   natural route here.  A one-electron system is where such an expansion can actually be constructed.
    #
    #   THIS BRANCH FIRST GAVE THE OPPOSITE ANSWER, and the mistake is worth keeping.  Run in a 60 a.u. box, the
    #   ratios came out 0.775, 0.553, 0.396, 0.198 -- moving monotonically AWAY from unity, which would have said
    #   that adding bound states makes matters worse and that the whole bound route is unsound.  That was not
    #   physics: the box was too small for the ejected electron, and the corrupted continuum orbital grew worse as
    #   the intermediate sum grew larger.  The two effects were indistinguishable in the answer.  Only repeating
    #   the study in a box wide enough for the outgoing wave separated them.
    setDefaults("method: continuum, Galerkin");   setDefaults("method: normalization, pure sine")
    setDefaults("print summary: open", "zzz-MultiPhotonIonization-Dv-convergence.sum")

    Z = 1.0
    grid = Radial.Grid(Radial.Grid(false); rnt=2.0e-6, h=5.0e-2, hp=1.0e-2, rbox=90.0)
    nm   = Nuclear.Model(Z, Nuclear.PointNucleus(), 0., 0., AngularJ64(0), 0., 0., 0.)
    prim = Bsplines.generatePrimitives(grid);    pot = Nuclear.nuclearPotential(nm, grid)
    println("\n  convergence of the bound intermediate sum, H 1s, omega = 0.30 a.u. [atomic units]\n")
    @printf("     %-8s %-10s %14s %14s %10s\n", "n <=", "orbitals", "linear Coul", "linear Bab", "ratio")
    for  nMax  in  [4, 6, 8, 10]
        shs = Subshell[]
        for n = 1:nMax
            push!(shs, Subshell(n,-1))
            if n >= 2   push!(shs, Subshell(n,1));   push!(shs, Subshell(n,-2))   end
            if n >= 3   push!(shs, Subshell(n,2));   push!(shs, Subshell(n,-3))   end
        end
        orbs = Bsplines.generateOrbitalsHydrogenic(shs, nm, prim; printout=false)
        lines = redirect_stdout(devnull) do
                    MultiPhotonIonization.computeLines2pOneElectron(Subshell("1s_1/2"), [0.30], [E1], orbs, pot)
                end
        c = lines[1].csLinear.Coulomb;    b = lines[1].csLinear.Babushkin
        @printf("     %-8d %-10d %14.5e %14.5e %10.4f\n", nMax, length(shs), c, b, c/b)
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      25-Aug-2026
    # Last successful: 25-Aug-2026
    #
    # Branch c: THE RESONANCE GUARD, which refuses rather than skips.  If a photon energy brings a real
    #   intermediate state onto the energy shell, the non-resonant second-order expression does not describe the
    #   process at all: the amplitude diverges, and what actually happens is a sequential two-step excitation
    #   through a real state with its own lifetime.  A bound-bound module with hundreds of intermediates may
    #   reasonably drop the offending term; here it would remove most of the amplitude and return a small number
    #   instead of no number.  So this module names the state and stops.
    #
    # REPORT (25-Aug-2026): asking for hydrogen at omega = 0.375 a.u., which is the 1s -> 2p resonance, the guard
    #   trips and names 2p_1/2.  Below the resonance, at omega = 0.35, the same computation runs normally.
    #
    #   THE TOLERANCE HAD TO BE MADE RELATIVE, and the first version was useless.  With an ABSOLUTE tolerance of
    #   1e-6 a.u. the guard did NOT trip at omega = 0.375: the computed 1s -> 2p transition sits at 0.3750046, so
    #   the denominator was 4.6e-6 -- larger than 1e-6, and therefore "safe" -- while inflating that single term
    #   by a factor of 2e5.  A guard that lets the textbook resonance of the textbook atom through is not a guard.
    #   It is now relative to the photon energy, refusing at |denominator| < tolerance * omega with a default
    #   tolerance of 1e-3, about 0.01 eV at optical energies.  WHERE "NEAR RESONANCE" BEGINS IS PHYSICS AND NOT
    #   ARITHMETIC, so the number is a parameter of the computation rather than a constant of the code.
    setDefaults("method: continuum, Galerkin");   setDefaults("method: normalization, pure sine")
    setDefaults("print summary: open", "zzz-MultiPhotonIonization-Dv-resonance.sum")

    Z = 1.0
    grid = Radial.Grid(Radial.Grid(false); rnt=2.0e-6, h=5.0e-2, hp=1.0e-2, rbox=90.0)
    nm   = Nuclear.Model(Z, Nuclear.PointNucleus(), 0., 0., AngularJ64(0), 0., 0., 0.)
    prim = Bsplines.generatePrimitives(grid);    pot = Nuclear.nuclearPotential(nm, grid)
    shs  = [Subshell("1s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), Subshell("3p_1/2"), Subshell("3p_3/2")]
    orbs = Bsplines.generateOrbitalsHydrogenic(shs, nm, prim; printout=false)
    println("\n  E(1s) = ", orbs[Subshell("1s_1/2")].energy, "     E(2p) = ", orbs[Subshell("2p_1/2")].energy)
    println("  the 1s -> 2p resonance sits at omega = ",
            orbs[Subshell("2p_1/2")].energy - orbs[Subshell("1s_1/2")].energy, " a.u.\n")
    for  omega  in  [0.35, 0.375]
        println("  ---- omega = $omega a.u. ----")
        try     redirect_stdout(devnull) do
                    MultiPhotonIonization.computeLines2pOneElectron(Subshell("1s_1/2"), [omega], [E1], orbs, pot)
                end
                println("    computed normally.")
        catch e println("    REFUSED: ", first(split(sprint(showerror, e), "\n"))[1:min(end,120)], " ...")
        end
    end
    setDefaults("print summary: close", "")
    #
end

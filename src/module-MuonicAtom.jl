
"""
`module  JenaAtomicCalculator.MuonicAtom`  
    ... a submodel of JAC that computes the bound states of a MUON in an atom, and the muonic X-ray energies that
        follow from them.

        THE PHYSICAL PICTURE, and it is what makes this module small.  A muon is an electron with 206.77 times the
        mass: same charge, same spin, the same Dirac equation.  It is treated here as a SINGLE PARTICLE moving in
        the field of the nucleus plus, optionally, the mean field of the electrons.  It never enters an amplitude
        together with the electrons, and there is no muon-electron exchange: the two systems see each other only
        through their average charge distributions.  That is an approximation, and a good one, because the muon
        orbits some two hundred times closer to the nucleus than any electron.

        WHY THE NUCLEUS MATTERS MORE THAN ANYTHING ELSE HERE.  The muon's innermost orbit has a radius of roughly
        256/Z fm, so for a heavy atom it lies INSIDE the nuclear charge distribution:

              Z = 6  (C)     muon 1s ~ 43 fm     nucleus ~ 2.8 fm     outside
              Z = 26 (Fe)    muon 1s ~ 9.8 fm    nucleus ~ 4.6 fm     outside
              Z = 82 (Pb)    muon 1s ~ 3.1 fm    nucleus ~ 7.1 fm     INSIDE

        The finite size of the nucleus is therefore not a small correction but the dominant effect above Z ~ 40:
        measured here, it shifts the muon 1s level by 0.4% in carbon and by 50% in lead.  That is precisely why
        muonic X-rays are a classical way of measuring nuclear charge radii, and it is the reason this module takes
        a Nuclear.Model seriously rather than assuming a point charge.

        WHAT IS NOT INCLUDED, and the first two matter at the per-cent level for heavy elements:
        + VACUUM POLARISATION, the dominant QED correction in muonic atoms, which increases the binding.  Its
          absence is visible: the 2p -> 1s energies computed here come out about 1% BELOW the measured ones for
          heavy elements, which is the right size and the right sign for this omission.
        + NUCLEAR POLARISATION and the finer nuclear-structure corrections.
        + transition RATES.  Only transition ENERGIES are formed; the muonic X-ray intensities need a
          one-particle multipole operator that is not written yet.
        + the muon CASCADE.  A muon is captured into a high orbit and cascades down; that belongs in a cascade
          scheme of its own and is deliberately not attempted here.
"""
module MuonicAtom

using  Printf, ..Basics, ..Bsplines, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings


"""
`MuonicAtom.recommendedGrid(Z::Float64, nMax::Int64)`  
    ... returns a Radial.Grid suited to the MUON's orbitals of an atom with nuclear charge Z, up to principal
        quantum number nMax.  A grid::Radial.Grid is returned.

        This exists because Rule 12 -- the box must be matched to the orbitals -- bites twice as hard for a muon.
        A muon orbit is ~207 times smaller than the corresponding electron orbit, so an electron grid is both far
        too coarse near the nucleus and enormously too wide; the B-spline basis then spends its splines on empty
        space and misrepresents exactly the innermost orbital one is after.  The outer turning point of a
        hydrogenic (n,l) orbit is r+ = (n^2/Z)(1 + sqrt(1 - l(l+1)/n^2)) / m, and the box is taken as 2.5 r+ for
        the l = n-1 member of the highest shell requested.
"""
function recommendedGrid(Z::Float64, nMax::Int64)
    mass = Defaults.getDefaults("mass: muon")
    rPlus = (nMax^2 / Z) * (1.0 + sqrt(max(0., 1.0 - (nMax-1)*nMax / nMax^2))) / mass
    return( Radial.Grid(Radial.Grid(false); rnt = 1.0e-9, h = 4.0e-2, hp = 0., rbox = 2.5*rPlus) )
end


"""
`MuonicAtom.computeOrbitals(subshells::Array{Subshell,1}, nm::Nuclear.Model, grid::Radial.Grid;
                            electronPotential::Union{Nothing,Radial.Potential}=nothing, printout::Bool=false)`  
    ... solves the one-particle Dirac equation for a MUON in the potential of the nucleus nm, optionally with the
        mean field of the electrons added.  A dictionary orbitals::Dict{Subshell,Orbital} is returned.

        `electronPotential` is the electrons' contribution in the Zr convention of Radial.Potential and is simply
        ADDED to the nuclear one; passing `nothing` leaves the muon in the bare nuclear field.  The screening it
        provides is small -- the muon sits far inside the electron cloud and sees almost the whole nuclear charge
        -- but it is what makes the muon levels depend on the charge state of the ion.
"""
function computeOrbitals(subshells::Array{Subshell,1}, nm::Nuclear.Model, grid::Radial.Grid;
                         electronPotential::Union{Nothing,Radial.Potential}=nothing, printout::Bool=false)
    mass = Defaults.getDefaults("mass: muon")
    pot  = Nuclear.nuclearPotential(nm, grid)
    if  !isnothing(electronPotential)
        if  length(electronPotential.Zr) != length(pot.Zr)
            error("MuonicAtom.computeOrbitals(): the electron potential is given on a different grid than the " *
                  "nuclear one; both must use the muon grid.")
        end
        pot = Radial.Potential(pot.name * " + electron mean field", pot.Zr + electronPotential.Zr, grid)
    end
    primitives = Bsplines.generatePrimitives(grid)

    return( Bsplines.generateOrbitals(subshells, pot, nm, primitives; printout=printout, mass=mass) )
end


"""
`MuonicAtom.screeningPotential(orbital::Orbital, grid::Radial.Grid)`  
    ... forms the potential that the bound MUON presents to the electrons, in the Zr convention of
        Radial.Potential.  A potential::Radial.Potential is returned, to be ADDED to the nuclear one when the
        electronic structure of a muonic atom is computed.

        The classical statement about muonic atoms falls out of this: because the muon orbits some two hundred
        times closer than any electron, the electrons cannot resolve it from the nucleus and see one unit of
        charge glued on.  Zr tends to 1 at electron distances, so a muonic atom of charge Z behaves chemically
        like an ordinary atom of charge Z-1 -- and that limit is the natural check on this function.
"""
function screeningPotential(orbital::Orbital, grid::Radial.Grid)
    nr = min(length(orbital.P), length(grid.r))
    D  = [ orbital.P[i]^2 + orbital.Q[i]^2  for i = 1:nr ]      # radial density of the muon
    Zr = zeros( length(grid.r) )
    ## Zr(r) = INT_0^r D dr'  +  r INT_r^inf (D/r') dr'   -- the standard Hartree screening of one particle.
    ## The integrals use the grid's OWN Gauss-Legendre weights, grid.wr, and not a hand-rolled trapezoidal rule:
    ## on the exponential mesh a trapezoid left the large-r limit at 1.000374 instead of 1, i.e. it broke the one
    ## exact statement this function has to satisfy.  A quadrature error is indistinguishable from a physics error
    ## in the answer, and only the exact limit tells them apart.
    tail = zeros(nr)
    outer = 0.
    for  i = nr:-1:2
        outer     = outer + D[i]/grid.r[i] * grid.wr[i]
        tail[i-1] = outer
    end
    inner = 0.
    for  i = 2:nr
        inner  = inner + D[i] * grid.wr[i]
        Zr[i]  = inner + grid.r[i] * tail[i]
    end
    for  i = nr+1:length(grid.r)    Zr[i] = Zr[nr]    end

    return( Radial.Potential("muon screening", Zr, grid) )
end


"""
`MuonicAtom.displayLevels(stream::IO, orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model)`  
    ... displays the muon binding energies of the given orbitals, beside the point-nucleus values they would have
        if the nuclear charge were concentrated at the origin.  Nothing is returned.
"""
function displayLevels(stream::IO, orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model)
    mass = Defaults.getDefaults("mass: muon");    nx = 96
    shs  = sort(collect(keys(orbitals)), by = sh -> (sh.n, sh.kappa))
    println(stream, " ")
    println(stream, "  Muon binding energies for Z = $(nm.Z), nuclear model $(nm.model):")
    println(stream, " ")
    println(stream, "  ", "-"^nx)
    println(stream, "     subshell        energy [keV]      point nucleus [keV]     finite-size shift")
    println(stream, "  ", "-"^nx)
    for  sh  in  shs
        en = Defaults.convertUnits("energy: from atomic to eV", orbitals[sh].energy) * 1.0e-3
        ep = Defaults.convertUnits("energy: from atomic to eV",
                                   Basics.computeDiracEnergy(sh, nm.Z; mass=mass)) * 1.0e-3
        println(stream, "     " * rpad(string(sh), 12) * @sprintf("%16.4f", en) * @sprintf("%20.4f", ep) *
                        @sprintf("%18.2f", 100*(en-ep)/abs(ep)) * " %")
    end
    println(stream, "  ", "-"^nx)
    println(stream, "    The finite-size shift is the whole of the physics here for a heavy nucleus: the muon 1s")
    println(stream, "    orbit has a radius of about 256/Z fm and so lies INSIDE the nuclear charge above Z ~ 40.")

    return( nothing )
end


"""
`MuonicAtom.displayTransitions(stream::IO, orbitals::Dict{Subshell,Orbital})`  
    ... displays the muonic X-ray energies, i.e. the differences of the muon binding energies, for every pair of
        the given subshells with a lower final level.  Nothing is returned.

        ENERGIES ONLY.  No rate or intensity is formed, so nothing here says which of these lines is strong; the
        multipole selection rules are not applied either, and a pair that E1 cannot connect is still listed.
"""
function displayTransitions(stream::IO, orbitals::Dict{Subshell,Orbital})
    nx  = 96
    shs = sort(collect(keys(orbitals)), by = sh -> orbitals[sh].energy)
    println(stream, " ")
    println(stream, "  Muonic X-ray energies (differences of the levels above):")
    println(stream, " ")
    println(stream, "  ", "-"^nx)
    println(stream, "     transition                 energy [keV]        energy [MeV]")
    println(stream, "  ", "-"^nx)
    for  i = 1:length(shs),  f = 1:length(shs)
        if  orbitals[shs[i]].energy <= orbitals[shs[f]].energy    continue    end
        en = Defaults.convertUnits("energy: from atomic to eV", orbitals[shs[i]].energy - orbitals[shs[f]].energy) * 1.0e-3
        println(stream, "     " * rpad(string(shs[i]) * " --> " * string(shs[f]), 24) *
                        @sprintf("%16.4f", en) * @sprintf("%20.6f", en*1.0e-3))
    end
    println(stream, "  ", "-"^nx)
    println(stream, "    Energies only: no rates, no intensities, and the multipole selection rules are NOT")
    println(stream, "    applied, so a pair that no E1 photon can connect is still listed above.")

    return( nothing )
end

end # module

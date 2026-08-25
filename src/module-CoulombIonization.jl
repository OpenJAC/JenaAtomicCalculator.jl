
"""
`module  JenaAtomicCalculator.CoulombIonization`  
    ... a submodel of JAC that computes Coulomb-IONIZATION cross sections and the alignment of the residual ion, for
        the ionization of target or projectile electrons by fast ion impact.

        THE THEORY IS THAT OF CoulombExcitation, and this module follows it deliberately -- the same amplitude, the
        same momentum-transfer integration, the same alignment parameters -- with one physical difference and its
        consequences:

            excitation      A^(q+) + |i(N)>   -->   A^(q+) + |f(N)>
            ionization      A^(q+) + |i(N)>   -->   A^(q+) + |f(N-1)> + e^-

        THE ELECTRON LEAVES, AND THAT CHANGES THREE THINGS.

        + The final state is not a level but a level TIMES a free electron.  It is built here by attaching a
          continuum orbital to the residual ion and coupling the two to a total symmetry, so that both sides of the
          matrix element carry N electrons and the many-electron machinery of CoulombExcitation applies unchanged.
        + THE NORMALISATION.  A continuum orbital is normalised per unit ENERGY, so what this module returns is a
          cross section DIFFERENTIAL in the energy of the ejected electron, d(sigma)/d(epsilon), at each requested
          electron energy -- NOT a total cross section.  A total needs an integration over epsilon, which is left to
          the caller because the sensible energy range depends on the ion and the collision.  Every table says so.
        + The minimum momentum transfer grows with the electron energy: q0 = (E_f - E_i + epsilon)/(beta c), since
          the projectile must supply the binding energy AND the electron's kinetic energy.

        THE ALIGNMENT PARAMETERS ARE THE SAME EXPRESSIONS, and that is worth stating because it is not quite
        obvious.  A_2 and A_4 describe the residual ION, so they are formed from the cross section resolved by the
        ion's magnetic quantum number M_f -- exactly as in CoulombExcitation.  What is new is that the computed
        amplitudes are labelled by the TOTAL symmetry of (ion + electron), so reaching sigma(M_f) requires that
        combination to be decoupled again:

            K(M_f, m)  =  SUM_Jt  <J_f M_f, j m | J_t M_t>  K_(J_t M_t)         with  M_t = M_f + m

        a COHERENT sum over the total symmetries.  The total cross section does not need this -- the
        Clebsch-Gordan transformation is unitary, so summing |K|^2 over (J_t, M_t) and over (M_f, m) give the same
        number -- which is why the decoupling is done only when the alignment is actually requested.

        APPROXIMATIONS, all inherited from CoulombExcitation unless stated:
        + first-order perturbation theory in the projectile-electron interaction; no higher-order or
          coupled-channel effects, so this is a fast-collision treatment.
        + all cross sections are computed for a projectile of unit charge, Z_p = 1; multiply by Z_p^2 otherwise.
        + the projectile moves on a straight line and is not deflected.
        + the ejected electron is described in the STATIC field of the residual ion, with no exchange between it
          and the bound electrons beyond what the mean field provides.
        + the partial-wave sum over the ejected electron is truncated at `lValues`, and a truncation that is too
          severe does not merely lower a cross section -- it can distort its shape; see example-Dp.jl branch b.
"""
module CoulombIonization


using  Printf, ..AngularMomentum, ..Basics, ..Continuum, ..CoulombExcitation, ..Defaults, ..ManyElectron, ..Nuclear,
       ..Radial, ..RadialIntegrals, ..SpinAngular, ..TableStrings

## The reduced matrix elements and the projectile velocity are TAKEN FROM CoulombExcitation rather than copied.
## The theory is the same one -- Surzhykov et al., Phys. Rev. A 77, 042722 (2008), Eq. (8) -- and the corrections
## made there, in particular the -i on every magnetic term that RATIP carries and the paper does not show, are
## corrections to this module too.  A copy would have had to be corrected twice.

"""
`struct  CoulombIonization.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing Coulomb-ionization amplitudes and the
        energy-differential cross sections that follow from them.

    + ionEnergies             ::Array{Float64,1}    ... List of projectile ion energies [MeV/u].
    + electronEnergies        ::Array{Float64,1}    ... List of energies of the EJECTED electron, in the currently
                                                         selected energy unit.
    + calcAlignment           ::Bool                ... True, if the alignment of the residual ion is to be computed.
    + printBefore             ::Bool                ... True, if all energies and lines are printed before evaluation.
    + lineSelection           ::LineSelection       ... Specifies the selected levels, if any.
    + zerosGL                 ::Int64               ... Number of Gauss-Legendre zeros in the momentum-transfer integral.
    + lValues                 ::Array{Int64,1}      ... Orbital angular momenta of the ejected electron to be included.
"""
struct Settings  <:  AbstractProcessSettings 
    ionEnergies               ::Array{Float64,1} 
    electronEnergies          ::Array{Float64,1}
    calcAlignment             ::Bool 
    printBefore               ::Bool 
    lineSelection             ::LineSelection
    zerosGL                   ::Int64
    lValues                   ::Array{Int64,1}
end 


"""
`CoulombIonization.Settings()`  ... constructor for the default values of Coulomb-ionization line computations.
"""
function Settings()
    Settings(Float64[], Float64[], false, false, LineSelection(), 0, Int64[])
end


"""
`CoulombIonization.Settings(set::CoulombIonization.Settings;`

        ionEnergies=..,         electronEnergies=..,        calcAlignment=..,       printBefore=..,
        lineSelection=..,       zerosGL=..,                 lValues=..)
                    
    ... constructor for modifying the given CoulombIonization.Settings by 'overwriting' the previously selected
        parameters.  A settings::CoulombIonization.Settings is returned.
"""
function Settings(set::CoulombIonization.Settings;    
    ionEnergies::Union{Nothing,Array{Float64,1}}=nothing,           electronEnergies::Union{Nothing,Array{Float64,1}}=nothing,
    calcAlignment::Union{Nothing,Bool}=nothing,                     printBefore::Union{Nothing,Bool}=nothing,
    lineSelection::Union{Nothing,LineSelection}=nothing,            zerosGL::Union{Nothing,Int64}=nothing,
    lValues::Union{Nothing,Array{Int64,1}}=nothing)  
    
    if  isnothing(ionEnergies)         ionEnergiesx      = set.ionEnergies       else  ionEnergiesx      = ionEnergies       end 
    if  isnothing(electronEnergies)    electronEnergiesx = set.electronEnergies  else  electronEnergiesx = electronEnergies  end 
    if  isnothing(calcAlignment)       calcAlignmentx    = set.calcAlignment     else  calcAlignmentx    = calcAlignment     end 
    if  isnothing(printBefore)         printBeforex      = set.printBefore       else  printBeforex      = printBefore       end 
    if  isnothing(lineSelection)       lineSelectionx    = set.lineSelection     else  lineSelectionx    = lineSelection     end 
    if  isnothing(zerosGL)             zerosGLx          = set.zerosGL           else  zerosGLx          = zerosGL           end 
    if  isnothing(lValues)             lValuesx          = set.lValues           else  lValuesx          = lValues           end 
    
    Settings( ionEnergiesx, electronEnergiesx, calcAlignmentx, printBeforex, lineSelectionx, zerosGLx, lValuesx)
end


# `Base.show(io::IO, settings::CoulombIonization.Settings)`  ... prepares a proper printout of the settings.
function Base.show(io::IO, settings::CoulombIonization.Settings) 
    println(io, "ionEnergies:              $(settings.ionEnergies)  ")
    println(io, "electronEnergies:         $(settings.electronEnergies)  ")
    println(io, "calcAlignment:            $(settings.calcAlignment)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "zerosGL:                  $(settings.zerosGL)  ")
    println(io, "lValues:                  $(settings.lValues)  ")
end


"""
`struct  CoulombIonization.Channel`  
    ... defines a type for one q-, kappa- and symmetry-dependent Coulomb-ionization amplitude.

    + kappa          ::Int64                ... partial wave of the ejected electron.
    + m              ::AngularM64           ... magnetic quantum number of the ejected electron.
    + symmetry       ::LevelSymmetry        ... total symmetry J_t of (residual ion + ejected electron).
    + phase          ::Float64              ... phase of the continuum partial wave.
    + q              ::Float64              ... momentum transfer of this amplitude.
    + w              ::Float64              ... Gauss-Legendre weight of this q in the integration.
    + amplitude      ::Complex{Float64}     ... K^(Coulion) amplitude of this channel.
"""
struct  Channel
    kappa            ::Int64
    m                ::AngularM64
    symmetry         ::LevelSymmetry
    phase            ::Float64
    q                ::Float64  
    w                ::Float64 
    amplitude        ::Complex{Float64} 
end


# `Base.show(io::IO, channel::CoulombIonization.Channel)`  ... prepares a proper printout of the channel.
function Base.show(io::IO, channel::CoulombIonization.Channel) 
    println(io, "kappa:         $(channel.kappa)  ")
    println(io, "m:             $(channel.m)  ")
    println(io, "symmetry:      $(channel.symmetry)  ")
    println(io, "phase:         $(channel.phase)  ")
    println(io, "q:             $(channel.q)  ")
    println(io, "w:             $(channel.w)  ")
    println(io, "amplitude:     $(channel.amplitude)  ")
end


"""
`struct  CoulombIonization.MagneticLine`  
    ... defines a type for one pair (M_i, M_f) of magnetic quantum numbers of the initial level and the RESIDUAL
        ION, together with the differential cross section that belongs to it and its channels.

    + Mi             ::AngularM64           ... magnetic quantum number of the initial level.
    + Mf             ::AngularM64           ... magnetic quantum number of the residual ion.
    + partialCs      ::Float64              ... d(sigma)/d(epsilon) for this pair, summed over the ejected
                                                 electron's kappa and m.
    + channels       ::Array{CoulombIonization.Channel,1}  ... channels of this magnetic line.
"""
struct  MagneticLine
    Mi               ::AngularM64
    Mf               ::AngularM64
    partialCs        ::Float64
    channels         ::Array{CoulombIonization.Channel,1}
end


# `Base.show(io::IO, mLine::CoulombIonization.MagneticLine)`  ... prepares a proper printout of the magnetic line.
function Base.show(io::IO, mLine::CoulombIonization.MagneticLine) 
    println(io, "Mi:            $(mLine.Mi)  ")
    println(io, "Mf:            $(mLine.Mf)  ")
    println(io, "partialCs:     $(mLine.partialCs)  ")
    println(io, "channels:      $(mLine.channels)  ")
end


"""
`struct  CoulombIonization.Line`  
    ... defines a type for a Coulomb-ionization line, i.e. one initial level, one level of the residual ion, one
        projectile energy and one energy of the ejected electron.

    + initialLevel   ::Level                  ... initial level, with N electrons.
    + finalLevel     ::Level                  ... level of the residual ion, with N-1 electrons.
    + ionEnergy      ::Float64                ... projectile energy [MeV/u].
    + electronEnergy ::Float64                ... energy of the ejected electron [a.u.].
    + q0             ::Float64                ... minimum momentum transfer for this transition and energy.
    + dCs            ::Float64                ... d(sigma)/d(epsilon), summed over all magnetic sublevels.
    + alignmentA2    ::Float64                ... alignment A_2 of the residual ion.
    + alignmentA4    ::Float64                ... alignment A_4 of the residual ion.
    + mLines         ::Array{MagneticLine,1}  ... list of magnetic lines.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    ionEnergy        ::Float64
    electronEnergy   ::Float64
    q0               ::Float64
    dCs              ::Float64
    alignmentA2      ::Float64
    alignmentA4      ::Float64
    mLines           ::Array{CoulombIonization.MagneticLine,1}
end


# `Base.show(io::IO, line::CoulombIonization.Line)`  ... prepares a proper printout of the line.
function Base.show(io::IO, line::CoulombIonization.Line) 
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "ionEnergy:         $(line.ionEnergy)  ")
    println(io, "electronEnergy:    $(line.electronEnergy)  ")
    println(io, "q0:                $(line.q0)  ")
    println(io, "dCs:               $(line.dCs)  ")
    println(io, "alignmentA2:       $(line.alignmentA2)  ")
    println(io, "alignmentA4:       $(line.alignmentA4)  ")
    println(io, "mLines:            $(line.mLines)  ")
end


"""
`CoulombIonization.computeAmplitude(q::Float64, q0::Float64, beta::Float64, Mi::AngularM64, Mt::AngularM64,
                                    cLevel::Level, iLevel::Level, grid::Radial.Grid)`  
    ... computes the Coulomb-ionization amplitude K^(Coulion)(vec q; alpha_i J_i M_i --> (alpha_f J_f, eps kappa) J_t M_t)
        for one momentum transfer q.  An amplitude::ComplexF64 is returned.

        This is CoulombExcitation.computeAmplitude with the final level replaced by the COUPLED state of the
        residual ion and the ejected electron; the angular algebra is untouched, which is the whole point of
        building the final state that way.
"""
function  computeAmplitude(q::Float64, q0::Float64, beta::Float64, Mix::Float64, Mtx::Float64,
                           cLevel::Level, iLevel::Level, grid::Radial.Grid, key::Int64,
                           cache::Dict{Tuple{Symbol,Int64,Int64,Int64},ComplexF64})
    Ji  = iLevel.J;                   Jt  = cLevel.J;             amplitude = ComplexF64(0.)
    Jix = AngularMomentum.oneJ(Ji);   Jtx = AngularMomentum.oneJ(Jt)

    if  abs(Mix) > Jix  ||  abs(Mtx) > Jtx    return( amplitude )    end
    for  t  in  AngularMomentum.j_values(Ji, Jt)
        tx   = AngularMomentum.oneJ(t)
        ## the transferred M cannot exceed the rank t; the coefficient is zero there, but the Wigner package
        ## RAISES on an invalid (j,m) pair rather than returning zero, so it must not be reached at all
        if  abs(Mtx - Mix) > tx    continue    end
        wa   = AngularMomentum.ClebschGordan(Jix, Mix, tx, Mtx-Mix, Jtx, Mtx) / sqrt(2*Jtx + 1.0)
        if  wa == 0.   continue   end
        Mval = Mtx - Mix;    wb = ComplexF64(0.)
        for  L  in  AngularMomentum.j_values(t, AngularJ64(1))
            Lx = AngularMomentum.oneJ(L);   Lint = Int64(Lx);   Mint = Int64(Mval)
            if  abs(Mint) > Lint   continue   end
            wc = ComplexF64(0.)
            ## THE REDUCED MATRIX ELEMENTS DO NOT DEPEND ON ANY MAGNETIC QUANTUM NUMBER, so they are cached on
            ## (state, rank, L).  Without this they would be recomputed for every (M_i, M_f, m) combination --
            ## a many-electron matrix element each time, and the dominant cost of the whole module.
            if  t == L
                kY = (:Y, key, Lint, Lint)
                if  !haskey(cache, kY)   cache[kY] = CoulombExcitation.computeKjYme(cLevel, L, iLevel, q, grid)   end
                wc = wc + cache[kY]
            end
            kT = (:T, key, Int64(tx), Lint)
            if  !haskey(cache, kT)   cache[kT] = CoulombExcitation.computeKjTme(cLevel, t, L, iLevel, q, grid)   end
            wc = wc + im * beta * AngularMomentum.ClebschGordan(Lx, Mval, 1., 0., tx, Mval) * cache[kT]
            wc = im^Lx * conj( AngularMomentum.sphericalYlm(Lint, Mint, acos(q0/q), 0.) ) * wc
            wb = wb + wc
        end
        amplitude = amplitude + wa * wb
    end

    return( amplitude )
end


"""
`CoulombIonization.generateContinuumLevels(line::CoulombIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                           settings::CoulombIonization.Settings)`  
    ... builds, once per line, the coupled states of the residual ion and the ejected electron.  A tuple
        (iLevel::Level, states::Array{Tuple{Int64,Float64,LevelSymmetry,Level},1}) is returned, whose entries are
        (kappa, phase, total symmetry, coupled level) and whose iLevel carries the placeholder subshell that makes
        the two sides of the matrix element share one subshell list.

        This is the expensive part of a computation and it is done ONCE, outside the loops over magnetic quantum
        numbers and momentum transfers, which reuse the same continuum orbitals.
"""
function  generateContinuumLevels(line::CoulombIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                  settings::CoulombIonization.Settings)
    contSettings = Continuum.Settings(false, grid.NoPoints-50)
    nuclearPot   = Nuclear.nuclearPotential(nm, grid)
    ## THE BASIS MUST BE REDUCED TO THE LEVEL'S OWN SYMMETRY FIRST.  A multiplet's basis carries the CSFs of every
    ## symmetry of its configuration -- 1s2s holds both J=0 and J=1 -- and Basics.generateLevelWithExtraElectron
    ## requires that every CSF share the level's own J and parity, raising "Improper symmetry of CSF." otherwise.
    ## The He-like case hides this, since its final level has a single CSF; the first many-level case does not.
    subshells    = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
    iLev         = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshells)
    fLev         = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshells)
    states       = Tuple{Int64,Float64,LevelSymmetry,Level,Level}[]
    kappas       = Int64[]
    for  l  in  settings.lValues
        push!(kappas, -(l+1));    if  l > 0    push!(kappas, l)    end
    end
    for  kappa  in  unique(kappas)
        cSubshell = Subshell(101, kappa)
        ## Each kappa needs its OWN placeholder on the initial side, since the placeholder IS the subshell
        ## (101, kappa); carrying one from a previous kappa would leave the two sides with different subshell
        ## lists, which is exactly what the placeholder exists to prevent.
        newiLevel = Basics.generateLevelWithExtraSubshell(cSubshell, iLev)
        cOrbital, phase = Continuum.generateOrbitalForLevel(line.electronEnergy, cSubshell, fLev, nm, grid,
                                                            contSettings; nuclearPot=nuclearPot)
        ## the coupled state has the ion's parity times (-1)^l of the ejected electron
        cParity = iseven( Basics.subshell_l(cSubshell) ) ? line.finalLevel.parity :
                                                           Basics.invertParity(line.finalLevel.parity)
        for  Jt  in  AngularMomentum.j_values(line.finalLevel.J, AngularMomentum.kappa_j(kappa))
            symt      = LevelSymmetry(Jt, cParity)
            newcLevel = Basics.generateLevelWithExtraElectron(cOrbital, symt, fLev)
            push!( states, (kappa, phase, symt, newcLevel, newiLevel) )
        end
    end

    return( (iLev, states) )
end



"""
`CoulombIonization.computeAmplitudesProperties(line::CoulombIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                               settings::CoulombIonization.Settings; printout::Bool=true)`  
    ... computes all amplitudes, the energy-differential cross section and, if requested, the alignment of the
        residual ion for the given line.  A newLine::CoulombIonization.Line is returned.
"""
function  computeAmplitudesProperties(line::CoulombIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                      settings::CoulombIonization.Settings; printout::Bool=true)
    beta   = CoulombExcitation.betaProjectile(line.ionEnergy)
    alpha  = Defaults.getDefaults("alpha")
    (iLev, states) = CoulombIonization.generateContinuumLevels(line, nm, grid, settings)
    if  length(states) == 0
        return( CoulombIonization.Line(line.initialLevel, line.finalLevel, line.ionEnergy, line.electronEnergy,
                                        line.q0, 0., 0., 0., CoulombIonization.MagneticLine[]) )
    end
    gl     = Radial.GridGL(Radial.GridGaussLegendreFinite(), line.q0, 10*line.q0, settings.zerosGL)
    qs     = gl.t;    ws = gl.wt;    nq = length(qs)
    cache  = Dict{Tuple{Symbol,Int64,Int64,Int64},ComplexF64}()
    Ji     = line.initialLevel.J;    Jf = line.finalLevel.J
    Jfx    = AngularMomentum.oneJ(Jf)
    kappas = unique( st[1] for st in states )
    #
    newmLines = CoulombIonization.MagneticLine[]
    for  Mi  in  AngularMomentum.m_values(Ji)
        Mix = AngularMomentum.oneM(Mi)
        for  Mf  in  AngularMomentum.m_values(Jf)
            Mfx = AngularMomentum.oneM(Mf);    chs = CoulombIonization.Channel[];    wSum = 0.
            for  kappa  in  kappas
                j   = AngularMomentum.kappa_j(kappa);    jx = AngularMomentum.oneJ(j)
                for  m  in  AngularMomentum.m_values(j)
                    mx  = AngularMomentum.oneM(m);    Mtx = Mfx + mx
                    for  iq = 1:nq
                        ## COHERENT over the total symmetries J_t, incoherent over everything else: the (ion +
                        ## electron) state is decoupled here so that the ION's magnetic quantum number M_f is
                        ## resolved, which is what the alignment is about.
                        coh = ComplexF64(0.)
                        for  (si, st)  in  enumerate(states)
                            if  st[1] != kappa   continue   end
                            Jtx = AngularMomentum.oneJ(st[3].J)
                            ## M_t = M_f + m is only a magnetic quantum number of those total symmetries with
                            ## J_t >= |M_t|; the others contribute nothing and must not be asked for
                            if  abs(Mtx) > Jtx    continue    end
                            cg  = AngularMomentum.ClebschGordan(Jfx, Mfx, jx, mx, Jtx, Mtx)
                            if  cg == 0.   continue   end
                            amp = CoulombIonization.computeAmplitude(qs[iq], line.q0, beta, Mix, Mtx, st[4], st[5],
                                                                     grid, (si-1)*nq + iq, cache)
                            coh = coh + cg * amp
                            push!( chs, CoulombIonization.Channel(kappa, m, st[3], st[2], qs[iq], ws[iq], amp) )
                        end
                        if  coh == ComplexF64(0.)   continue   end
                        wa   = qs[iq]^2 - line.q0^2 * beta^2
                        wSum = wSum + qs[iq] * ws[iq] / (wa^2) * abs2(coh)
                    end
                end
            end
            partialCs = wSum * 2pi * (8pi * alpha / beta)^2 / (Basics.twice(Ji) + 1)
            push!( newmLines, CoulombIonization.MagneticLine(Mi, Mf, partialCs, chs) )
        end
    end
    #
    dCs = 0.;   for mLine in newmLines   dCs = dCs + mLine.partialCs   end
    alignmentA2 = 0.;   alignmentA4 = 0.
    ## The SAME expressions as CoulombExcitation, Eqs. (9)-(10) of Surzhykov et al. (2008): the alignment belongs
    ## to the residual ion, so it is formed from the cross section resolved by that ion's M_f.
    if  settings.calcAlignment  &&  dCs > 0.
        MfList  = AngularMomentum.m_values(Jf)
        sigmaMf = Dict{AngularM64, Float64}( Mf => 0. for Mf in MfList )
        for  mLine in newmLines   sigmaMf[mLine.Mf] = sigmaMf[mLine.Mf] + mLine.partialCs   end
        wa2 = 0.;   wa4 = 0.
        for  Mf in MfList
            Mfx   = AngularMomentum.oneM(Mf);    phase = (-1)^Int64( round(Jfx - Mfx) )
            wa2   = wa2 + phase * AngularMomentum.ClebschGordan(Jfx, Mfx, Jfx, -Mfx, 2., 0.) * sigmaMf[Mf]
            wa4   = wa4 + phase * AngularMomentum.ClebschGordan(Jfx, Mfx, Jfx, -Mfx, 4., 0.) * sigmaMf[Mf]
        end
        alignmentA2 = sqrt(2*Jfx + 1) * wa2 / dCs
        alignmentA4 = sqrt(2*Jfx + 1) * wa4 / dCs
    end

    return( CoulombIonization.Line(line.initialLevel, line.finalLevel, line.ionEnergy, line.electronEnergy,
                                    line.q0, dCs, alignmentA2, alignmentA4, newmLines) )
end


"""
`CoulombIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                                  settings::CoulombIonization.Settings)`  
    ... determines the list of Coulomb-ionization lines for the given multiplets, projectile energies and
        ejected-electron energies.  An Array{CoulombIonization.Line,1} is returned, with all properties still zero.
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::CoulombIonization.Settings)
    lines = CoulombIonization.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  !Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)   continue   end
            for  ionEnergy  in  settings.ionEnergies
                beta = CoulombExcitation.betaProjectile(ionEnergy)
                for  en  in  settings.electronEnergies
                    electronEnergy = Defaults.convertUnits("energy: to atomic", en)
                    ## THE PROJECTILE MUST SUPPLY THE BINDING ENERGY AND THE ELECTRON'S KINETIC ENERGY, so the
                    ## minimum momentum transfer grows with the electron energy -- unlike the excitation case,
                    ## where it is fixed by the transition alone.
                    deltaE = fLevel.energy - iLevel.energy + electronEnergy
                    if  deltaE <= 0.   continue   end
                    q0 = deltaE / (beta * Defaults.getDefaults("speed of light: c"))
                    push!( lines, CoulombIonization.Line(iLevel, fLevel, ionEnergy, electronEnergy, q0, 0., 0., 0.,
                                                          CoulombIonization.MagneticLine[]) )
                end
            end
        end
    end

    return( lines )
end


"""
`CoulombIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                grid::Radial.Grid, settings::CoulombIonization.Settings; output=true)`  
    ... computes the Coulomb-ionization amplitudes and all requested properties.  An
        Array{CoulombIonization.Line,1} is returned if output is true, and nothing otherwise.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                       grid::Radial.Grid, settings::CoulombIonization.Settings; output=true)
    println("")
    printstyled("CoulombIonization.computeLines(): The computation of Coulomb ionization cross sections starts now ... \n",
                color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------ \n",
                color=:light_green)
    println("")
    sa = "\n* Coulomb-ionization cross sections are computed on the same footing as Coulomb excitation, with the " *
         "following assumptions: \n" *
         "\n    + Projectile energies are given in [MeV/u] and converted into a relative velocity beta = v/c. " *
         "\n    + Cross sections are computed for a projectile of unit charge Z_p = 1; multiply by Z_p^2 otherwise. " *
         "\n    + WHAT IS RETURNED IS d(sigma)/d(epsilon), DIFFERENTIAL in the energy of the ejected electron, " *
         "because a continuum orbital is normalised per unit energy.  A total cross section needs an integration " *
         "over epsilon and is NOT formed here. " *
         "\n    + The ejected electron is described in the static field of the residual ion, with partial waves " *
         "limited to lValues. \n"
    println(sa)

    lines = CoulombIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    if  settings.printBefore    CoulombIonization.displayLines(stdout, lines)    end
    newLines = CoulombIonization.Line[]
    for  line in lines
        push!( newLines, CoulombIonization.computeAmplitudesProperties(line, nm, grid, settings) )
    end
    CoulombIonization.displayCrossSections(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    CoulombIonization.displayCrossSections(iostream, newLines, settings)   end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`CoulombIonization.displayLines(stream::IO, lines::Array{CoulombIonization.Line,1})`  
    ... lists the selected lines before their evaluation.  Nothing is returned.
"""
function  displayLines(stream::IO, lines::Array{CoulombIonization.Line,1})
    nx = 104
    println(stream, " ")
    println(stream, "  Selected Coulomb-ionization lines:")
    println(stream, " ")
    println(stream, "  ", "-"^nx)
    println(stream, "       i-level-f        i--J^P--f        E_p [MeV/u]     eps(e-) [eV]        q0 [a.u.]")
    println(stream, "  ", "-"^nx)
    for  line in lines
        sa = "     " * TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index) * "   " *
             TableStrings.symmetries_if(LevelSymmetry(line.initialLevel.J, line.initialLevel.parity),
                                        LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity))
        println(stream, sa * @sprintf("%14.3f", line.ionEnergy) *
                @sprintf("%17.3f", Defaults.convertUnits("energy: from atomic to eV", line.electronEnergy)) *
                @sprintf("%18.5e", line.q0))
    end
    println(stream, "  ", "-"^nx)

    return( nothing )
end


"""
`CoulombIonization.displayCrossSections(stream::IO, lines::Array{CoulombIonization.Line,1},
                                        settings::CoulombIonization.Settings)`  
    ... displays the energy-differential cross sections and, if requested, the alignment of the residual ion.
        Nothing is returned.
"""
function  displayCrossSections(stream::IO, lines::Array{CoulombIonization.Line,1},
                               settings::CoulombIonization.Settings)
    nx = 118
    println(stream, " ")
    println(stream, "  Coulomb-ionization cross sections, DIFFERENTIAL in the energy of the ejected electron:")
    println(stream, " ")
    println(stream, "  ", "-"^nx)
    sa = "       i-level-f        i--J^P--f        E_p [MeV/u]    eps(e-) [eV]      q0 [a.u.]     " *
         "d(sigma)/d(eps) [b/eV]"
    if  settings.calcAlignment    sa = sa * "        A_2         A_4"    end
    println(stream, sa)
    println(stream, "  ", "-"^nx)
    for  line in lines
        ## the cross section is in atomic units per atomic unit of energy; report it in barn per eV, which is what
        ## an ejected-electron spectrum is usually plotted in
        dcs = Defaults.convertUnits("cross section: from atomic to barn", line.dCs) /
              Defaults.convertUnits("energy: from atomic to eV", 1.0)
        sa  = "     " * TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index) * "   " *
              TableStrings.symmetries_if(LevelSymmetry(line.initialLevel.J, line.initialLevel.parity),
                                         LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)) *
              @sprintf("%14.3f", line.ionEnergy) *
              @sprintf("%16.3f", Defaults.convertUnits("energy: from atomic to eV", line.electronEnergy)) *
              @sprintf("%16.5e", line.q0) * @sprintf("%18.5e", dcs)
        if  settings.calcAlignment
            sa = sa * @sprintf("%13.5f", line.alignmentA2) * @sprintf("%12.5f", line.alignmentA4)
        end
        println(stream, sa)
    end
    println(stream, "  ", "-"^nx)
    println(stream, "    d(sigma)/d(eps) is DIFFERENTIAL in the ejected-electron energy and is NOT a total cross")
    println(stream, "    section; integrate it over eps to obtain one.  Computed for a projectile of unit charge,")
    println(stream, "    Z_p = 1; multiply by Z_p^2 for any other projectile.")
    if  settings.calcAlignment
        println(stream, "    A_2 and A_4 describe the alignment of the RESIDUAL ION and are formed from exactly the")
        println(stream, "    expressions used for Coulomb excitation; reaching them requires the (ion + electron)")
        println(stream, "    state to be decoupled again, which is done coherently over the total symmetries.")
    end

    return( nothing )
end


end # module

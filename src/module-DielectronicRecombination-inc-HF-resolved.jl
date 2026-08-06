
## Hyperfine-resolved dielectronic recombination -- a strict analogue of the fine-structure route in
## module-DielectronicRecombination-inc-FS-resolved.jl.
##
## REWRITTEN COMPLETELY on 05-Aug-2026. The previous version had never run: computeHyperfineAmplitudes,
## displayHyperfineResults and displayHyperfineRateCoefficients were each called but never defined (the
## amplitude function was spelled computecomputeHyperfineAmplitudes), displayHyperfinePassages was handed an
## undefined `passages`, the non-distributed path passed a Multiplet where an Hfs.HfMultiplet was required, and
## it carried the threadid()/nthreads() crash fixed elsewhere. It also constructed a DielectronicRecombination
## .Passage, which was the last reason that retired type was still alive. See git history for the old file.
##
## THE IDEA. Both operators of DR are purely ELECTRONIC -- the Coulomb interaction for the capture step, the
## multipole field for the photon step -- so neither touches the nuclear spin. The nucleus is a spectator, and a
## hyperfine amplitude is therefore the corresponding ELECTRONIC amplitude times a recoupling coefficient. No
## electronic amplitude is ever recomputed here: this file runs the fine-structure route to obtain them, builds
## the hyperfine multiplets, and recouples. That is what keeps the hyperfine route affordable, and it is also
## what makes it a genuine analogue rather than a second implementation that can drift out of step.
##
## WHERE THE PHYSICS ACTUALLY IS. The INITIAL ion carries the full nuclear moments: its hyperfine splitting is
## observable, the capture resonances split by F_i, and their strengths carry (2F_i+1) in place of (2J_i+1). The
## intermediate and final levels split too, but negligibly, which is imitated by mu = Q = 0. That is not an
## approximation coded around: with mu = Q = 0 the hyperfine Hamiltonian vanishes identically, so
## Hfs.computeHyperfineRepresentation returns pure |(I J) F> states with mc = [.. 1.0 ..] of its own accord.
## Their F is then pure bookkeeping and is summed over before anything is displayed.


"""
`struct  DielectronicRecombination.HfCaptureLine`
    ... defines a type for a hyperfine-resolved dielectronic capture line, i.e. the step i + e- --> m between two
        hyperfine levels, together with the two total widths of the intermediate hyperfine level. It is the exact
        pendant of DielectronicRecombination.CaptureLine, with Hfs.HfLevel in place of Level.

    + initialLevel      ::Hfs.HfLevel   ... initial hyperfine level
    + intermediateLevel ::Hfs.HfLevel   ... intermediate (resonant) hyperfine level
    + electronEnergy    ::Float64       ... energy of the captured electron
    + captureRate       ::Float64       ... A_a(m --> i) for THIS channel only
    + totalAugerRate    ::Float64       ... Gamma_a(m), summed over all initial hyperfine levels
    + totalPhotonRate   ::EmProperty    ... Gamma_r(m), summed over all final hyperfine levels
    + resonanceStrength ::EmProperty    ... S(i,m)

        Unlike CaptureLine this type carries NO channel list. The channels live on the electronic capture line
        from which this one was recoupled; keeping a second copy here would duplicate state that cannot be kept
        consistent, which is the failure mode this whole re-structuring exists to remove.
"""
struct  HfCaptureLine
    initialLevel        ::Hfs.HfLevel
    intermediateLevel   ::Hfs.HfLevel
    electronEnergy      ::Float64
    captureRate         ::Float64
    totalAugerRate      ::Float64
    totalPhotonRate     ::EmProperty
    resonanceStrength   ::EmProperty
end


"""
`DielectronicRecombination.HfCaptureLine()`  ... constructor for an 'empty' HfCaptureLine.
"""
function HfCaptureLine()
    em = EmProperty(0., 0.)
    HfCaptureLine(Hfs.HfLevel(), Hfs.HfLevel(), 0., 0., 0., em, em)
end


# `Base.show(io::IO, line::DielectronicRecombination.HfCaptureLine)`  ... prepares a proper printout.
function Base.show(io::IO, line::DielectronicRecombination.HfCaptureLine)
    println(io, "initialLevel (F):           $(line.initialLevel.F)  ")
    println(io, "intermediateLevel (F):      $(line.intermediateLevel.F)  ")
    println(io, "electronEnergy:             $(line.electronEnergy)  ")
    println(io, "captureRate:                $(line.captureRate)  ")
    println(io, "totalAugerRate:             $(line.totalAugerRate)  ")
    println(io, "totalPhotonRate:            $(line.totalPhotonRate)  ")
    println(io, "resonanceStrength:          $(line.resonanceStrength)  ")
end


"""
`struct  DielectronicRecombination.HfPhotonLine`
    ... defines a type for one hyperfine-resolved radiative stabilization m --> f + hv; the pendant of
        DielectronicRecombination.PhotonLine.

    + intermediateLevel ::Hfs.HfLevel   ... intermediate (resonant) hyperfine level
    + finalLevel        ::Hfs.HfLevel   ... final hyperfine level
    + photonEnergy      ::Float64       ... energy of the emitted photon
    + photonRate        ::EmProperty    ... A_r(m,f)
"""
struct  HfPhotonLine
    intermediateLevel   ::Hfs.HfLevel
    finalLevel          ::Hfs.HfLevel
    photonEnergy        ::Float64
    photonRate          ::EmProperty
end


"""
`DielectronicRecombination.HfPhotonLine()`  ... constructor for an 'empty' HfPhotonLine.
"""
function HfPhotonLine()
    HfPhotonLine(Hfs.HfLevel(), Hfs.HfLevel(), 0., EmProperty(0., 0.))
end


# `Base.show(io::IO, line::DielectronicRecombination.HfPhotonLine)`  ... prepares a proper printout.
function Base.show(io::IO, line::DielectronicRecombination.HfPhotonLine)
    println(io, "intermediateLevel (F):      $(line.intermediateLevel.F)  ")
    println(io, "finalLevel (F):             $(line.finalLevel.F)  ")
    println(io, "photonEnergy:               $(line.photonEnergy)  ")
    println(io, "photonRate:                 $(line.photonRate)  ")
end


"""
`DielectronicRecombination.hfPhotonRecoupling(spinI::AngularJ64, Ja::AngularJ64, Fa::AngularJ64,
                            Jb::AngularJ64, Fb::AngularJ64, L::Int64)`
    ... to return the factor that converts an ELECTRONIC reduced matrix element of rank L into the corresponding
        hyperfine one, for an operator that acts on the electrons alone:

            <(I J_b) F_b || T^L || (I J_a) F_a>
                = (-1)^(I+J_b+F_a+L) sqrt((2F_a+1)(2F_b+1)) {J_b F_b I; F_a J_a L} <J_b || T^L || J_a>

        A value::Float64 is returned. Derived from the standard result for a tensor acting on one part of a
        coupled system (Edmonds 7.1.7), with the nucleus as the untouched spectator.

        VERIFIED (work/diag-recoupling.jl, 05-Aug-2026) over 454 combinations of I, J_a, J_b, F_a and L = 1, 2:
        sum over F_b of the squared factor equals (2F_a+1)/(2J_a+1) to better than 1e-15. That identity is
        exactly what leaves a radiative width unchanged by hyperfine coupling, since the rate carries a
        compensating 1/(2F_a+1): Gamma_r(F_a) = Gamma_r(J_a). For I = 0 the factor reduces to 1.
"""
function  hfPhotonRecoupling(spinI::AngularJ64, Ja::AngularJ64, Fa::AngularJ64,
                             Jb::AngularJ64, Fb::AngularJ64, L::Int64)
    ## MOVED to Hfs on 06-Aug-2026 and kept here only as a name: the identical factor is needed by the
    ## hyperfine-INDUCED transitions too, and two copies of one formula is how they begin to diverge.
    return( Hfs.recouplingElectronicOperator(spinI, Ja, Fa, Jb, Fb, L) )
end


"""
`DielectronicRecombination.hfCaptureRecoupling(spinI::AngularJ64, Ji::AngularJ64, Fi::AngularJ64,
                            je::AngularJ64, Jm::AngularJ64, Fm::AngularJ64)`
    ... to return the coefficient that connects the two ways of coupling nucleus, initial ion and free electron:

            <((I J_i) F_i, j_e) F_m | (I, (J_i j_e) J_m) F_m>
                = (-1)^(I+J_i+j_e+F_m) sqrt((2F_i+1)(2J_m+1)) {I J_i F_i; j_e F_m J_m}

        A value::Float64 is returned. This is NOT an operator matrix element but a recoupling of three angular
        momenta (Edmonds 6.1.5): the physical state has the free electron coupled to the hyperfine level F_i,
        whereas the electronic capture amplitude is computed with it coupled to the electronic J_i to give J_m.

        VERIFIED (work/diag-recoupling.jl, 05-Aug-2026) over 515 combinations: at FIXED F_m, the sum over J_m of
        the squared coefficient is 1 to better than 1e-15, as it must be, since the two coupling schemes are
        complete bases of the same space and the transformation between them is unitary. Note that summing over
        F_m as well is NOT an identity -- that counts independent final states. For I = 0 the factor reduces to 1.
"""
function  hfCaptureRecoupling(spinI::AngularJ64, Ji::AngularJ64, Fi::AngularJ64,
                              je::AngularJ64, Jm::AngularJ64, Fm::AngularJ64)
    wa = AngularMomentum.phaseFactor([spinI, +1, Ji, +1, je, +1, Fm])
    wb = sqrt( (Basics.twice(Fi) + 1) * (Basics.twice(Jm) + 1) )
    wc = AngularMomentum.Wigner_6j(spinI, Ji, Fi, je, Fm, Jm)
    return( wa * wb * wc )
end


"""
`DielectronicRecombination.electronicComponents(hfLevel::Hfs.HfLevel)`
    ... to return the list of (mc, electronic Level) pairs that make up the given hyperfine level, dropping the
        components whose weight is numerically zero. An Array{Tuple{Float64,Level},1} is returned.

        With mu = Q = 0 -- the intermediate and final levels of this route -- every hyperfine level is a pure
        |(I J) F> and this returns exactly one component with mc = 1.
"""
function  electronicComponents(hfLevel::Hfs.HfLevel)
    comps = Tuple{Float64,Level}[]
    for  k in eachindex(hfLevel.mc)
        if  abs(hfLevel.mc[k]) > 1.0e-10    push!(comps, (hfLevel.mc[k], hfLevel.hfBasisVectors[k].levelJ))    end
    end
    return( comps )
end


"""
`DielectronicRecombination.determineHfCaptureLines(intermediateHfMultiplet::Hfs.HfMultiplet,
                            initialHfMultiplet::Hfs.HfMultiplet, settings::DielectronicRecombination.Settings)`
    ... to determine the skeletons of all hyperfine capture lines (i,m) that are energetically allowed, i.e. for
        which the captured electron has a positive energy. An Array{HfCaptureLine,1} is returned.
        The pendant of DielectronicRecombination.determineCaptureLines.
"""
function  determineHfCaptureLines(intermediateHfMultiplet::Hfs.HfMultiplet, initialHfMultiplet::Hfs.HfMultiplet,
                                  settings::DielectronicRecombination.Settings)
    lines = DielectronicRecombination.HfCaptureLine[];    em = EmProperty(0., 0.)
    eShift = Defaults.convertUnits("energy: to atomic", settings.electronEnergyShift)
    for  iLevel in initialHfMultiplet.hfLevels
        for  mLevel in intermediateHfMultiplet.hfLevels
            energy = mLevel.energy - iLevel.energy + eShift
            if  energy <= 0.    continue    end
            push!( lines, DielectronicRecombination.HfCaptureLine(iLevel, mLevel, energy, 0., 0., em, em) )
        end
    end
    return( lines )
end


"""
`DielectronicRecombination.determineHfPhotonLines(finalHfMultiplet::Hfs.HfMultiplet,
                            intermediateHfMultiplet::Hfs.HfMultiplet, settings::DielectronicRecombination.Settings)`
    ... to determine the skeletons of all hyperfine radiative stabilization lines (m,f) with a positive photon
        energy. An Array{HfPhotonLine,1} is returned. The pendant of determinePhotonLines.
"""
function  determineHfPhotonLines(finalHfMultiplet::Hfs.HfMultiplet, intermediateHfMultiplet::Hfs.HfMultiplet,
                                 settings::DielectronicRecombination.Settings)
    lines = DielectronicRecombination.HfPhotonLine[]
    pShift = Defaults.convertUnits("energy: to atomic", settings.photonEnergyShift)
    minEn  = Defaults.convertUnits("energy: to atomic", settings.mimimumPhotonEnergy)
    for  mLevel in intermediateHfMultiplet.hfLevels
        for  fLevel in finalHfMultiplet.hfLevels
            energy = mLevel.energy - fLevel.energy + pShift
            if  energy <= 0.  ||  energy < minEn    continue    end
            push!( lines, DielectronicRecombination.HfPhotonLine(mLevel, fLevel, energy, EmProperty(0., 0.)) )
        end
    end
    return( lines )
end


"""
`DielectronicRecombination.computeHfCaptureAmplitudes(hfLine::DielectronicRecombination.HfCaptureLine,
                            eCaptureLines::Dict{Tuple{Int64,Int64},DielectronicRecombination.CaptureLine},
                            spinI::AngularJ64)`
    ... to obtain the capture rate of one hyperfine capture line by RECOUPLING the electronic capture amplitudes
        that the fine-structure route has already computed. A new HfCaptureLine is returned.

        For each partial wave kappa separately -- different kappa are distinct, incoherent channels -- the
        components of the two hyperfine levels are summed COHERENTLY:

            A(kappa) = sum_{p,q} mc_i[p] mc_m[q] * <((I J_i^p) F_i, j_e) F_m | (I,(J_i^p j_e) J_m^q) F_m>
                                                 * A_electronic(i^p, kappa --> m^q)
            captureRate = 2 pi * sum_kappa |A(kappa)|^2

        The 2 pi and the sum over kappa are exactly as in the fine-structure computeCaptureAmplitudes, so that the
        two routes coincide term by term when I = 0.
"""
function  computeHfCaptureAmplitudes(hfLine::DielectronicRecombination.HfCaptureLine,
                                     eCaptureLines::Dict{Tuple{Int64,Int64},DielectronicRecombination.CaptureLine},
                                     spinI::AngularJ64)
    iComps = DielectronicRecombination.electronicComponents(hfLine.initialLevel)
    mComps = DielectronicRecombination.electronicComponents(hfLine.intermediateLevel)
    Fi     = hfLine.initialLevel.F;      Fm = hfLine.intermediateLevel.F
    ## Collect the partial waves that occur in any of the contributing electronic lines
    kappas = Int64[]
    for  (mci, iLev) in iComps,  (mcm, mLev) in mComps
        eLine = get(eCaptureLines, (iLev.index, mLev.index), nothing)
        if  eLine === nothing    continue    end
        for  ch in eLine.captureChannels    if  !(ch.kappa in kappas)    push!(kappas, ch.kappa)    end    end
    end
    #
    rateA = 0.
    for  kappa in kappas
        je  = Basics.subshell_j( Subshell(101, kappa) )
        amp = ComplexF64(0.)
        for  (mci, iLev) in iComps,  (mcm, mLev) in mComps
            eLine = get(eCaptureLines, (iLev.index, mLev.index), nothing)
            if  eLine === nothing    continue    end
            wa = DielectronicRecombination.hfCaptureRecoupling(spinI, iLev.J, Fi, je, mLev.J, Fm)
            if  wa == 0.    continue    end
            for  ch in eLine.captureChannels
                if  ch.kappa != kappa    continue    end
                amp = amp + mci * mcm * wa * ch.amplitude
            end
        end
        rateA = rateA + abs(amp)^2
    end
    #
    return( DielectronicRecombination.HfCaptureLine(hfLine.initialLevel, hfLine.intermediateLevel,
                                                    hfLine.electronEnergy, 2pi * rateA, 0., EmProperty(0., 0.),
                                                    EmProperty(0., 0.)) )
end


"""
`DielectronicRecombination.computeHfPhotonAmplitudes(hfLine::DielectronicRecombination.HfPhotonLine,
                            ePhotonLines::Dict{Tuple{Int64,Int64},DielectronicRecombination.PhotonLine},
                            spinI::AngularJ64)`
    ... to obtain the radiative rate of one hyperfine photon line by recoupling the electronic photon amplitudes
        already computed by the fine-structure route. A new HfPhotonLine is returned.

        Each (multipole, gauge) is an incoherent channel; within one, the components of the two hyperfine levels
        are summed coherently with the rank-L recoupling factor. The Einstein prefactor is that of the
        fine-structure route with (2J_m+1) replaced by (2F_m+1), the statistical weight of the EMITTING level:

            Gamma_r contribution = 8 pi alpha omega / (2F_m + 1) * |A|^2

        Together with the identity sum_{F_f} |recoupling|^2 = (2F_m+1)/(2J_m+1) this gives Gamma_r(F_m) =
        Gamma_r(J_m), i.e. hyperfine coupling leaves the radiative width of a level untouched -- which is the
        first thing to check if this route is ever in doubt.
"""
function  computeHfPhotonAmplitudes(hfLine::DielectronicRecombination.HfPhotonLine,
                                    ePhotonLines::Dict{Tuple{Int64,Int64},DielectronicRecombination.PhotonLine},
                                    spinI::AngularJ64)
    mComps = DielectronicRecombination.electronicComponents(hfLine.intermediateLevel)
    fComps = DielectronicRecombination.electronicComponents(hfLine.finalLevel)
    Fm     = hfLine.intermediateLevel.F;    Ff = hfLine.finalLevel.F
    ## Collect the (multipole, gauge) channels that occur
    mpGauges = Tuple{EmMultipole,EmGauge}[]
    for  (mcm, mLev) in mComps,  (mcf, fLev) in fComps
        eLine = get(ePhotonLines, (mLev.index, fLev.index), nothing)
        if  eLine === nothing    continue    end
        for  ch in eLine.photonChannels
            if  !((ch.multipole, ch.gauge) in mpGauges)    push!(mpGauges, (ch.multipole, ch.gauge))    end
        end
    end
    #
    rateC = 0.;    rateB = 0.
    for  (mp, gauge) in mpGauges
        amp = ComplexF64(0.)
        for  (mcm, mLev) in mComps,  (mcf, fLev) in fComps
            eLine = get(ePhotonLines, (mLev.index, fLev.index), nothing)
            if  eLine === nothing    continue    end
            wa = DielectronicRecombination.hfPhotonRecoupling(spinI, mLev.J, Fm, fLev.J, Ff, mp.L)
            if  wa == 0.    continue    end
            for  ch in eLine.photonChannels
                if  ch.multipole != mp  ||  ch.gauge != gauge    continue    end
                amp = amp + mcm * mcf * wa * ch.amplitude
            end
        end
        if       gauge == Basics.Coulomb     rateC = rateC + abs(amp)^2
        elseif   gauge == Basics.Babushkin   rateB = rateB + abs(amp)^2
        elseif   gauge == Basics.Magnetic    rateB = rateB + abs(amp)^2;   rateC = rateC + abs(amp)^2
        end
    end
    #
    wa = 8.0pi * Defaults.getDefaults("alpha") * hfLine.photonEnergy / (Basics.twice(Fm) + 1)
    #
    return( DielectronicRecombination.HfPhotonLine(hfLine.intermediateLevel, hfLine.finalLevel,
                                                  hfLine.photonEnergy, EmProperty(wa * rateC, wa * rateB)) )
end


"""
`DielectronicRecombination.setHfTotalRates(hfCaptureLines::Array{HfCaptureLine,1},
                            hfPhotonLines::Array{HfPhotonLine,1})`
    ... to accumulate, for every intermediate HYPERFINE level m, the two total widths

            Gamma_a(m) = sum_i A_a(m --> i)          Gamma_r(m) = sum_f A_r(m --> f)

        and to store them, with the resulting resonance strength, on each hyperfine capture line. A new
        Array{HfCaptureLine,1} is returned. The single aggregation point of this route, exactly as
        setTotalRates is for the fine-structure one, and for the same reason: it is where the capture and the
        photon side couple, and it is the step both old implementations got wrong.

        Keyed on (F, energy) rather than a level index, because Hfs.HfLevel carries no index; two hyperfine
        levels of equal F and equal energy are the same level.
"""
function  setHfTotalRates(hfCaptureLines::Array{DielectronicRecombination.HfCaptureLine,1},
                          hfPhotonLines::Array{DielectronicRecombination.HfPhotonLine,1})
    key(lev) = (Basics.twice(lev.F), round(lev.energy, digits=10))
    totalAuger  = Dict{Tuple{Int64,Float64},Float64}()
    totalPhoton = Dict{Tuple{Int64,Float64},EmProperty}()
    for  cLine in hfCaptureLines
        k = key(cLine.intermediateLevel)
        totalAuger[k]  = get(totalAuger,  k, 0.) + cLine.captureRate
    end
    for  pLine in hfPhotonLines
        k = key(pLine.intermediateLevel)
        totalPhoton[k] = get(totalPhoton, k, EmProperty(0., 0.)) + pLine.photonRate
    end
    #
    newLines = DielectronicRecombination.HfCaptureLine[]
    for  cLine in hfCaptureLines
        k      = key(cLine.intermediateLevel)
        gammaA = get(totalAuger,  k, 0.)
        gammaR = get(totalPhoton, k, EmProperty(0., 0.))
        ## C(i,m) = pi^2/k^2 * A_a(m --> i) * (2F_m+1)/(2F_i+1);  the statistical weights are those of the
        ## HYPERFINE levels, which is the one place where the hyperfine resolution enters the strength directly.
        wavenb = Defaults.convertUnits("kinetic energy to wave number: atomic units", cLine.electronEnergy)
        factor = pi*pi / (wavenb*wavenb) * cLine.captureRate *
                 ((Basics.twice(cLine.intermediateLevel.F) + 1) / (Basics.twice(cLine.initialLevel.F) + 1))
        totC   = gammaA + gammaR.Coulomb;      totB = gammaA + gammaR.Babushkin
        sC     = totC == 0.  ?  0.  :  factor * gammaR.Coulomb   / totC
        sB     = totB == 0.  ?  0.  :  factor * gammaR.Babushkin / totB
        push!( newLines, DielectronicRecombination.HfCaptureLine(cLine.initialLevel, cLine.intermediateLevel,
                                cLine.electronEnergy, cLine.captureRate, gammaA, gammaR, EmProperty(sC, sB)) )
    end
    return( newLines )
end


"""
`DielectronicRecombination.computeHfCaptureLines(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet,
                            initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                            settings::DielectronicRecombination.Settings; output::Bool=true)`
    ... to compute hyperfine-resolved dielectronic recombination. The driver of this file and the analogue of
        computeCaptureLines. A tuple (hfCaptureLines, hfPhotonLines) is returned if output = true.

        The three ELECTRONIC multiplets are passed in exactly as for the fine-structure route; the hyperfine
        multiplets are built here. The initial ion keeps the full nuclear moments of nm, whereas the intermediate
        and final ions are given mu = Q = 0, so that their hyperfine Hamiltonian vanishes and their levels come
        out as pure |(I J) F> degenerate with their electronic parents. Their F is retained because it carries the
        angular factors, and summed over on display.
"""
function  computeHfCaptureLines(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                                nm::Nuclear.Model, grid::Radial.Grid, settings::DielectronicRecombination.Settings;
                                output::Bool=true)
    println("")
    printstyled("DielectronicRecombination.computeHfCaptureLines(): The computation of hyperfine-resolved DR starts now ... \n",
                color=:light_green)
    printstyled("----------------------------------------------------------------------------------------------------------- \n",
                color=:light_green)
    println("")
    if  length(settings.corrections) > 0
        error("\n\nDielectronicRecombination.computeHfCaptureLines():  STOP -- the hyperfine-resolved route does not \n" *
              "implement the high-n corrections; they were requested as\n\n    $(settings.corrections)\n\n"             *
              ">>> Drop them, or use the fine-structure route, where they are implemented and verified.\n")
    end
    if  nm.spinI == AngularJ64(0)
        @warn("computeHfCaptureLines(): the nuclear spin is 0, so every hyperfine level coincides with its " *
              "electronic parent and this route can only reproduce the fine-structure one.")
    end
    #
    ## (1) THE ELECTRONIC SIDE, through the fine-structure machinery. Nothing is recomputed afterwards.
    DielectronicRecombination.checkConsistentMultiplets(finalMultiplet, intermediateMultiplet, initialMultiplet)
    DielectronicRecombination.checkOrbitalRepresentation(finalMultiplet, intermediateMultiplet, initialMultiplet)
    eCaptureLines = DielectronicRecombination.determineCaptureLines(finalMultiplet, intermediateMultiplet,
                                                                    initialMultiplet, settings)
    ePhotonLines  = DielectronicRecombination.determinePhotonLines( finalMultiplet, intermediateMultiplet,
                                                                    initialMultiplet, settings)
    maxEnergy = 0.;   for cLine in eCaptureLines   maxEnergy = max(maxEnergy, cLine.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    #
    newEPhotonLines = Vector{DielectronicRecombination.PhotonLine}(undef, length(ePhotonLines))
    @threads for  p in eachindex(ePhotonLines)
        newEPhotonLines[p] = DielectronicRecombination.computePhotonAmplitudes(ePhotonLines[p], grid)
    end
    newECaptureLines = Vector{DielectronicRecombination.CaptureLine}(undef, length(eCaptureLines))
    @threads for  c in eachindex(eCaptureLines)
        newECaptureLines[c] = DielectronicRecombination.computeCaptureAmplitudes(eCaptureLines[c], nm, grid,
                                                                                 nrContinuum, settings)
    end
    ## The electronic lines are aggregated too, so that their resonance strengths exist and can serve as the
    ## reference of the F-sum rule below. Without this they carry zero strength and the sum rule compares
    ## against nothing -- which is exactly what happened on the first run of this route.
    empTreatment     = DielectronicRecombination.determineEmpiricalTreatment(finalMultiplet, intermediateMultiplet,
                                                                             nm, initialMultiplet, settings)
    newECaptureLines = DielectronicRecombination.setTotalRates(newECaptureLines, newEPhotonLines, empTreatment)
    #
    eCapDict = Dict{Tuple{Int64,Int64},DielectronicRecombination.CaptureLine}()
    for  cLine in newECaptureLines   eCapDict[(cLine.initialLevel.index, cLine.intermediateLevel.index)] = cLine   end
    ePhoDict = Dict{Tuple{Int64,Int64},DielectronicRecombination.PhotonLine}()
    for  pLine in newEPhotonLines    ePhoDict[(pLine.intermediateLevel.index, pLine.finalLevel.index)]   = pLine   end
    println(">>> $(length(newECaptureLines)) electronic capture lines and $(length(newEPhotonLines)) electronic " *
            "photon lines computed; these are now recoupled.")
    #
    ## (2) THE HYPERFINE MULTIPLETS. Full moments for the initial ion, mu = Q = 0 for the other two.
    nmZero  = Nuclear.Model(nm; mu=0., Q=0.)
    iHfMult = Hfs.computeHyperfineRepresentation(Hfs.defineHyperfineBasis(initialMultiplet, nm; printout=false), grid)
    mHfMult = Hfs.computeHyperfineRepresentation(Hfs.defineHyperfineBasis(intermediateMultiplet, nmZero; printout=false), grid)
    fHfMult = Hfs.computeHyperfineRepresentation(Hfs.defineHyperfineBasis(finalMultiplet, nmZero; printout=false), grid)
    println(">>> hyperfine levels:  initial $(length(iHfMult.hfLevels)),  intermediate $(length(mHfMult.hfLevels)), " *
            " final $(length(fHfMult.hfLevels));   nuclear spin I = $(nm.spinI)")
    #
    ## (3) RECOUPLE
    hfCaptureLines = DielectronicRecombination.determineHfCaptureLines(mHfMult, iHfMult, settings)
    hfPhotonLines  = DielectronicRecombination.determineHfPhotonLines(fHfMult, mHfMult, settings)
    newHfPhoton    = Vector{DielectronicRecombination.HfPhotonLine}(undef, length(hfPhotonLines))
    @threads for  p in eachindex(hfPhotonLines)
        newHfPhoton[p] = DielectronicRecombination.computeHfPhotonAmplitudes(hfPhotonLines[p], ePhoDict, nm.spinI)
    end
    newHfCapture   = Vector{DielectronicRecombination.HfCaptureLine}(undef, length(hfCaptureLines))
    @threads for  c in eachindex(hfCaptureLines)
        newHfCapture[c] = DielectronicRecombination.computeHfCaptureAmplitudes(hfCaptureLines[c], eCapDict, nm.spinI)
    end
    newHfCapture = DielectronicRecombination.setHfTotalRates(newHfCapture, newHfPhoton)
    #
    ## (4) DISPLAY, resolved in F_i and summed over F_m and F_f
    DielectronicRecombination.displayHfResults(stdout, newHfCapture, newECaptureLines, nm)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    DielectronicRecombination.displayHfResults(iostream, newHfCapture, newECaptureLines, nm)   end
    #
    if  output    return( (newHfCapture, newHfPhoton) )
    else          return( nothing )
    end
end


"""
`DielectronicRecombination.displayHfResults(stream::IO, hfCaptureLines::Array{HfCaptureLine,1},
                            eCaptureLines::Array{DielectronicRecombination.CaptureLine,1}, nm::Nuclear.Model)`
    ... to list the hyperfine-resolved DR resonance strengths, resolved in F_i and SUMMED over F_m. Nothing is
        returned.

        Why summed over F_m (and, inside Gamma_r, over F_f): the intermediate and final levels were given
        mu = Q = 0, so their hyperfine sublevels are exactly degenerate and their F is a bookkeeping label with
        no observable consequence. Displaying it would suggest a structure that has been set to zero and would
        bury the one splitting that IS observable, that of the initial ion.

        THE SUM RULE is evaluated and printed here rather than left to a separate script, because it is the only
        check that this route is right and it costs nothing once the numbers are in hand. With every hyperfine
        level degenerate with its electronic parent,

            sum_{F_i} (2F_i+1)/((2I+1)(2J_i+1)) * sum_{F_m} S(F_i --> F_m)   ==   sum_m S(J_i --> J_m)

        i.e. the statistically averaged hyperfine strength must reproduce the fine-structure one EXACTLY -- it is
        an identity among recoupling coefficients, not an approximation. A deviation localises immediately to a
        wrong phase, a missing sqrt(2F+1) or a mis-ordered 6-j.
"""
function  displayHfResults(stream::IO, hfCaptureLines::Array{DielectronicRecombination.HfCaptureLine,1},
                           eCaptureLines::Array{DielectronicRecombination.CaptureLine,1}, nm::Nuclear.Model)
    ## Aggregate over F_m, keeping F_i and the capture energy
    byFi = Dict{Int64,EmProperty}();    enOf = Dict{Int64,Float64}();   nOf = Dict{Int64,Int64}()
    for  cLine in hfCaptureLines
        twoFi = Basics.twice(cLine.initialLevel.F)
        byFi[twoFi] = get(byFi, twoFi, EmProperty(0., 0.)) + cLine.resonanceStrength
        enOf[twoFi] = get(enOf, twoFi, cLine.electronEnergy)
        nOf[twoFi]  = get(nOf, twoFi, 0) + 1
    end
    twoFis = sort(collect(keys(byFi)))
    #
    nx = 104
    println(stream, " ")
    println(stream, "  Hyperfine-resolved DR resonance strengths, summed over F_m and F_f:")
    println(stream, " ")
    println(stream, "  nuclear spin I = $(nm.spinI),  mu = $(nm.mu),  Q = $(nm.Q)")
    println(stream, "  the intermediate and final levels carry mu = Q = 0, so their F is summed over")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "F_i"; na=4);        sb = sb * TableStrings.center(10, "        "; na=4)
    sa = sa * TableStrings.center(10, "lines"; na=4);      sb = sb * TableStrings.center(10, "        "; na=4)
    sa = sa * TableStrings.center(16, "E_e  " * TableStrings.inUnits("energy"); na=4)
    sb = sb * TableStrings.center(16, "                "; na=4)
    sa = sa * TableStrings.center(18, "shift from lowest F_i"; na=2)
    sb = sb * TableStrings.center(18, "[meV]"; na=2)
    sa = sa * TableStrings.center(30, "Cou -- strength -- Bab"; na=4)
    sb = sb * TableStrings.center(30, TableStrings.inUnits("strength"); na=4)
    println(stream, sa);   println(stream, sb);   println(stream, "  ", TableStrings.hLine(nx))
    #
    twoI = Basics.twice(nm.spinI)
    ## The capture energy is E_m - E_i(F_i), so the resonances are shifted by MINUS the initial hyperfine
    ## splitting. Printing that shift separately, in meV, is the only way to see it: it sits many orders of
    ## magnitude below the resonance energy itself and cannot survive the %.4e of the energy column.
    enRef = minimum([enOf[t] for t in twoFis])
    for  twoFi in twoFis
        shift = Defaults.convertUnits("energy: from atomic", enOf[twoFi] - enRef) * 1000.0
        sa  = "  " * TableStrings.center(10, string(AngularJ64(twoFi//2)); na=4)
        sa  = sa * TableStrings.center(10, string(nOf[twoFi]); na=4)
        sa  = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", enOf[twoFi]))            * "     "
        sa  = sa * @sprintf("%12.6f", shift)                                                             * "        "
        sa  = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", byFi[twoFi].Coulomb))  * "   "
        sa  = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", byFi[twoFi].Babushkin))* "       "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    if  length(twoFis) > 1
        println(stream, ">>> The spread of the capture energies over F_i IS the hyperfine splitting of the INITIAL ion,")
        println(stream, "    with reversed sign: a more strongly bound F_i needs a correspondingly larger electron energy.")
    end
    #
    ## THE SUM RULE.  The statistical weight of one hyperfine level within its electronic parent is
    ## (2F_i+1)/((2I+1)(2J_i+1)); summing the hyperfine strengths with that weight must return the electronic sum.
    hfSumC = 0.;   hfSumB = 0.
    for  cLine in hfCaptureLines
        twoFi = Basics.twice(cLine.initialLevel.F)
        ## the electronic parent of this hyperfine level, i.e. its dominant component
        comps = DielectronicRecombination.electronicComponents(cLine.initialLevel)
        twoJi = length(comps) > 0 ? Basics.twice(comps[argmax([abs(c[1]) for c in comps])][2].J) : 0
        w     = (twoFi + 1) / ((twoI + 1) * (twoJi + 1))
        hfSumC = hfSumC + w * cLine.resonanceStrength.Coulomb
        hfSumB = hfSumB + w * cLine.resonanceStrength.Babushkin
    end
    eSumC = 0.;   eSumB = 0.
    for  cLine in eCaptureLines
        eSumC = eSumC + cLine.resonanceStrength.Coulomb;    eSumB = eSumB + cLine.resonanceStrength.Babushkin
    end
    println(stream, " ")
    println(stream, "  F-SUM RULE -- the statistically averaged hyperfine strength must reproduce the fine-structure one:")
    println(stream, " ")
    @printf(stream, "    sum_F (2F_i+1)/((2I+1)(2J_i+1)) * S(hyperfine)   Coulomb %.8e   Babushkin %.8e\n", hfSumC, hfSumB)
    @printf(stream, "    sum   S(fine structure)                          Coulomb %.8e   Babushkin %.8e\n", eSumC, eSumB)
    if  eSumC != 0.  &&  eSumB != 0.
        @printf(stream, "    ratio                                            Coulomb %.10f   Babushkin %.10f\n",
                hfSumC/eSumC, hfSumB/eSumB)
        dev = max(abs(hfSumC/eSumC - 1.), abs(hfSumB/eSumB - 1.))
        if  dev < 1.0e-8
            println(stream, "    >>> the sum rule holds to $(round(dev, sigdigits=2)); the recoupling is consistent.")
        else
            println(stream, "    >>> WARNING: the sum rule is violated by $(round(100*dev, digits=4)) %. This is an identity")
            println(stream, "        among recoupling coefficients, so any deviation is a defect -- look for a wrong phase,")
            println(stream, "        a missing sqrt(2F+1), or a mis-ordered 6-j, NOT for a physical explanation.")
        end
    end
    println(stream, " ")
    #
    return( nothing )
end

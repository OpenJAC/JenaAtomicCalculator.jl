
"""
`module  JAC.CrystalFieldEmission`
    ... a submodel of JAC that contains all methods for computing crystal-field-resolved multipole
        transitions ("lines") between the Stark sublevels of an initial and a final multiplet, both
        split by the same external point-charge lattice (CrystalField.Lattice).

        This is the TRANSITION/PROCESS side of the crystal-field capability, in contrast to
        CrystalField.jl's LEVEL/PROPERTY side (CrystalField.Outcome/CfMultiplet describe a single
        multiplet's own Stark splitting). JAC draws this line consistently everywhere else in the
        codebase (Basics.AbstractPropertySettings -> Outcome, e.g. Hfs, IsotopeShift, LandeZeeman,
        CrystalField itself; versus Basics.AbstractProcessSettings -> Line, e.g. PhotoEmission) --
        no existing module mixes both, so this capability lives in its own module rather than being
        added to CrystalField.jl.

        CrystalFieldEmission.Settings/Channel/Line mirror PhotoEmission.Settings/Channel/Line
        directly, with CrystalField.CfLevel in place of ManyElectron.Level, and are invoked exactly
        like any other JAC process: initialConfigs/finalConfigs are declared together in an
        Atomic.Computation, and Basics.perform(...) calls
        CrystalFieldEmission.computeLines(finalMultiplet, initialMultiplet, grid, settings)
        automatically once computation.processSettings isa CrystalFieldEmission.Settings.

        The reduced (Wigner-Eckart) one-electron/CSF machinery is entirely reused from
        PhotoEmission.amplitude(...); this module only adds the M_J-resolved expansion over
        CrystalField.CfBasisVector's (mirroring CrystalField.computeInteractionMatrix's own 3-j/phase
        pattern, but between two INDEPENDENT bases -- initial- and final-multiplet-derived -- giving
        a genuinely non-Hermitian transition matrix, not the Hermitian crystal-field matrix) and the
        final mc-sandwich.
"""
module CrystalFieldEmission


using Printf, ..AngularMomentum, ..Basics, ..CrystalField, ..Defaults, ..ManyElectron, ..PhotoEmission, ..Radial, ..TableStrings


"""
`struct  CrystalFieldEmission.Settings  <:  Basics.AbstractProcessSettings`
    ... defines the settings for computing crystal-field-resolved multipole transitions between an
        initial and a final multiplet.

    + multipoles       ::Array{EmMultipole,1}                    ... multipoles to be included.
    + gauges           ::Array{UseGauge,1}                       ... gauges to be included.
    + lattice          ::CrystalField.Lattice                    ... the external point-charge lattice (same for both sides).
    + model            ::CrystalField.AbstractCrystalFieldModel  ... the crystal-field model to be applied.
    + maxRank          ::Int64                                   ... maximum tensor rank for the crystal-field splitting itself.
    + includeJmixing   ::Bool                                    ... whether to J-mix within each side's own splitting.
    + printBefore      ::Bool                                    ... True if a list of selected lines is printed before computing.
    + minEnergy        ::Float64                                 ... minimum transition energy (Hartree) to be included.
    + maxEnergy        ::Float64                                 ... maximum transition energy (Hartree) to be included.
"""
struct  Settings  <:  Basics.AbstractProcessSettings
    multipoles       ::Array{EmMultipole,1}
    gauges           ::Array{UseGauge,1}
    lattice          ::CrystalField.Lattice
    model            ::CrystalField.AbstractCrystalFieldModel
    maxRank          ::Int64
    includeJmixing   ::Bool
    printBefore      ::Bool
    minEnergy        ::Float64
    maxEnergy        ::Float64
end


"""
`CrystalFieldEmission.Settings()`  ... constructor for the default values of crystal-field-resolved transitions.
"""
function Settings()
    Settings(EmMultipole[E1], UseGauge[Basics.UseCoulomb], CrystalField.Lattice(), CrystalField.PointChargeModel(),
              6, false, false, 0., Inf64)
end


"""
`CrystalFieldEmission.Settings(set::CrystalFieldEmission.Settings;`

        multipoles=.., gauges=.., lattice=.., model=.., maxRank=.., includeJmixing=.., printBefore=..,
        minEnergy=.., maxEnergy=..)

    ... keyword copy-constructor for re-defining selected values of a settings::CrystalFieldEmission.Settings.
"""
function Settings(set::CrystalFieldEmission.Settings;
        multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,          gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
        lattice::Union{Nothing,CrystalField.Lattice}=nothing,             model::Union{Nothing,CrystalField.AbstractCrystalFieldModel}=nothing,
        maxRank::Union{Nothing,Int64}=nothing,                            includeJmixing::Union{Nothing,Bool}=nothing,
        printBefore::Union{Nothing,Bool}=nothing,
        minEnergy::Union{Nothing,Float64}=nothing,                        maxEnergy::Union{Nothing,Float64}=nothing)
    if  isnothing(multipoles)        multipolesx     = set.multipoles     else   multipolesx     = multipoles     end
    if  isnothing(gauges)            gaugesx         = set.gauges         else   gaugesx         = gauges         end
    if  isnothing(lattice)           latticex        = set.lattice        else   latticex        = lattice        end
    if  isnothing(model)             modelx          = set.model          else   modelx          = model          end
    if  isnothing(maxRank)           maxRankx        = set.maxRank        else   maxRankx        = maxRank        end
    if  isnothing(includeJmixing)    includeJmixingx = set.includeJmixing else   includeJmixingx = includeJmixing end
    if  isnothing(printBefore)       printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  isnothing(minEnergy)         minEnergyx      = set.minEnergy      else   minEnergyx      = minEnergy      end
    if  isnothing(maxEnergy)         maxEnergyx      = set.maxEnergy      else   maxEnergyx      = maxEnergy      end

    Settings( multipolesx, gaugesx, latticex, modelx, maxRankx, includeJmixingx, printBeforex, minEnergyx, maxEnergyx )
end


# `Base.show(io::IO, settings::CrystalFieldEmission.Settings)`  ... prepares a proper printout of settings::CrystalFieldEmission.Settings.
function Base.show(io::IO, settings::CrystalFieldEmission.Settings)
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "lattice:                  $(settings.lattice)  ")
    println(io, "model:                    $(settings.model)  ")
    println(io, "maxRank:                  $(settings.maxRank)  ")
    println(io, "includeJmixing:           $(settings.includeJmixing)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "minEnergy:                $(settings.minEnergy)  ")
    println(io, "maxEnergy:                $(settings.maxEnergy)  ")
end


"""
`struct  CrystalFieldEmission.Channel`
    ... a single multipole/gauge transition channel between two Stark sublevels.  It mirrors what
        PhotoEmission.Channel was before that module was retired; this one is CrystalFieldEmission's own and
        is unaffected.

    + multipole   ::EmMultipole   ... multipole of the photon emission/absorption.
    + gauge       ::EmGauge       ... gauge used for the radiation field.
    + amplitude   ::ComplexF64    ... crystal-field-resolved amplitude of this channel.
"""
struct  Channel
    multipole    ::EmMultipole
    gauge        ::EmGauge
    amplitude    ::ComplexF64
end


"""
`struct  CrystalFieldEmission.Line`
    ... a transition between one Stark sublevel of an initial CfMultiplet and one Stark sublevel of a
        final CfMultiplet, mirroring PhotoEmission.Line with CrystalField.CfLevel in place of
        ManyElectron.Level.

    + initialCfLevel   ::CrystalField.CfLevel               ... initial Stark sublevel.
    + finalCfLevel     ::CrystalField.CfLevel               ... final Stark sublevel.
    + omega            ::Float64                            ... transition energy (Hartree).
    + channels         ::Array{CrystalFieldEmission.Channel,1}  ... list of multipole/gauge channels.
"""
struct  Line
    initialCfLevel     ::CrystalField.CfLevel
    finalCfLevel       ::CrystalField.CfLevel
    omega              ::Float64
    channels           ::Array{Channel,1}
end


#################################################################################################################################
#################################################################################################################################


"""
`CrystalFieldEmission.computeLineAmplitude(mp::EmMultipole, gauge::EmGauge, finalCfLevel::CrystalField.CfLevel, initialCfLevel::CrystalField.CfLevel, grid::Radial.Grid)`
    ... computes the crystal-field-resolved multipole transition amplitude between one final and one
        initial Stark sublevel, following the Wigner-Eckart expansion
        <Jf Mf|O^mp_q|Ji Mi> = (-1)^(Jf-Mf) (Jf mp.L Ji; -Mf q Mi) <Jf||O^mp||Ji>, where q = Mf-Mi is
        fixed by the 3-j projection selection rule and the reduced matrix element <Jf||O^mp||Ji> is
        PhotoEmission.amplitude(Emission(),...) between the two CfBasisVectors' parent levels (cached
        per distinct parent-level pair, using the PARENT levels' own energy difference as the
        emission-operator's omega -- an excellent approximation since crystal-field-induced shifts
        are always tiny compared to a genuine parent-level transition energy). The resulting
        (final CfBasisVector) x (initial CfBasisVector) matrix is sandwiched between the two levels'
        mixing vectors: adjoint(finalCfLevel.mc) * matrix * initialCfLevel.mc (adjoint, not
        transpose, since CfLevel.mc is complex). A value::ComplexF64 is returned.
"""
function computeLineAmplitude(mp::EmMultipole, gauge::EmGauge, finalCfLevel::CrystalField.CfLevel, initialCfLevel::CrystalField.CfLevel,
                               grid::Radial.Grid)
    finalBasis = finalCfLevel.cfBasis;   initialBasis = initialCfLevel.cfBasis
    nf = length(finalBasis);   ni = length(initialBasis)
    matrix = zeros(ComplexF64, nf, ni)
    redamp = Dict{Tuple{Int64,Int64},ComplexF64}()
    for  bf = 1:nf
        for  bi = 1:ni
            fp = finalBasis[bf].parentLevel;   ip = initialBasis[bi].parentLevel
            key = (fp.index, ip.index)
            if  !haskey(redamp, key)
                omega = ip.energy - fp.energy
                if  omega <= 0.   redamp[key] = ComplexF64(0)
                else              redamp[key] = PhotoEmission.amplitude(PhotoEmission.Emission(), mp, gauge, omega, fp, ip, grid)
                end
            end
            amp = redamp[key]
            if  amp == ComplexF64(0)   continue    end
            Jf = fp.J;   Ji = ip.J;   Mf = finalBasis[bf].M;   Mi = initialBasis[bi].M
            q  = Int64( (Basics.twice(Mf) - Basics.twice(Mi)) ÷ 2 )
            if  abs(q) > mp.L   continue    end
            threej = AngularMomentum.Wigner_3j(Jf, AngularJ64(mp.L), Ji, AngularM64(-Mf.num//Mf.den), AngularM64(q//1), Mi)
            if  threej == 0.   continue    end
            phase = (-1.0)^Int64( (Basics.twice(Jf) - Basics.twice(Mf)) ÷ 2 )
            matrix[bf,bi] = phase * threej * amp
        end
    end
    return  adjoint(finalCfLevel.mc) * matrix * initialCfLevel.mc
end


"""
`CrystalFieldEmission.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, settings::CrystalFieldEmission.Settings)`
    ... the standard top-level driver of the module, called automatically by
        Basics.perform(::Atomic.Computation) once computation.processSettings isa
        CrystalFieldEmission.Settings (mirroring PhotoEmission.computeLines's role exactly --
        final-multiplet-first argument order, matching PhotoEmission throughout). Splits both
        `finalMultiplet` and `initialMultiplet` by the same settings.lattice via
        CrystalField.computeOutcomes(...), then forms the Cartesian product of all
        (initial CfLevel, final CfLevel) pairs across those outcomes, skipping pairs whose energy
        difference falls outside [settings.minEnergy, settings.maxEnergy]. For each surviving pair,
        one CrystalFieldEmission.Channel per requested (multipole,gauge) combination is computed via
        CrystalFieldEmission.computeLineAmplitude(...); as in
        PhotoEmission.determineChannels, only one Magnetic-gauge channel is kept per multipole
        regardless of how many gauges were requested. An Array{CrystalFieldEmission.Line,1} is returned.
"""
function computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, settings::Settings)
    # All levels of each multiplet are included by default (mirroring PhotoEmission.determineLines,
    # which iterates over the full initialMultiplet.levels x finalMultiplet.levels product unless a
    # LineSelection restricts it); CrystalField.computeOutcomes on its own would otherwise default to
    # just the lowest level of each multiplet, which is wrong here.
    finalSelection   = LevelSelection(true, indices=[ lev.index for lev in finalMultiplet.levels   ])
    initialSelection = LevelSelection(true, indices=[ lev.index for lev in initialMultiplet.levels ])
    cfSettingsFinal   = CrystalField.Settings(CrystalField.Settings(); lattice=settings.lattice, model=settings.model,
                                               maxRank=settings.maxRank, includeJmixing=settings.includeJmixing,
                                               printBefore=settings.printBefore, levelSelection=finalSelection)
    cfSettingsInitial = CrystalField.Settings(cfSettingsFinal; levelSelection=initialSelection)
    finalOutcomes   = CrystalField.computeOutcomes(finalMultiplet,   settings.lattice, grid, cfSettingsFinal)
    initialOutcomes = CrystalField.computeOutcomes(initialMultiplet, settings.lattice, grid, cfSettingsInitial)
    #
    if  settings.printBefore
        println("CrystalFieldEmission.computeLines(): $(length(finalOutcomes)) final and " *
                 "$(length(initialOutcomes)) initial crystal-field outcome(s) selected.")
    end
    #
    lines = Line[]
    for  fOutcome in finalOutcomes,  iOutcome in initialOutcomes
        for  fCfLevel in fOutcome.cfMultiplet.cfLevels,  iCfLevel in iOutcome.cfMultiplet.cfLevels
            omega = iCfLevel.energy - fCfLevel.energy
            if  omega <= settings.minEnergy  ||  omega > settings.maxEnergy   continue    end
            #
            channels = Channel[]
            for  mp in settings.multipoles
                hasMagnetic = false
                for  useGauge in settings.gauges
                    if      string(mp)[1] == 'E'  &&  useGauge == Basics.UseCoulomb      gauge = Basics.Coulomb
                    elseif  string(mp)[1] == 'E'  &&  useGauge == Basics.UseBabushkin    gauge = Basics.Babushkin
                    elseif  string(mp)[1] == 'M'  &&  !hasMagnetic                       gauge = Basics.Magnetic;   hasMagnetic = true
                    else    continue
                    end
                    amp = computeLineAmplitude(mp, gauge, fCfLevel, iCfLevel, grid)
                    push!(channels, Channel(mp, gauge, amp))
                end
            end
            if  length(channels) == 0   continue    end
            push!(lines, Line(iCfLevel, fCfLevel, omega, channels))
        end
    end
    #
    return  lines
end


"""
`CrystalFieldEmission.displayLines(stream::IO, lines::Array{CrystalFieldEmission.Line,1})`
    ... prints the crystal-field-resolved transition lines to `stream`, mirroring
        PhotoEmission.displayLines/displayRates: one row per Line, one column per requested
        (multipole,gauge) channel showing |amplitude|^2. Rows are labeled by sublevel energy rather
        than by level index, since a CfLevel has no single index (unlike a plain Level). Nothing is
        returned.

    + stream   ::IO                                    ... the output stream.
    + lines    ::Array{CrystalFieldEmission.Line,1}     ... the computed crystal-field-resolved lines.
"""
function displayLines(stream::IO, lines::Array{Line,1})
    println(stream, "\n  Crystal-field-resolved transition lines:\n")
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(28, "initial -- final energy [Hartree]"; na=2);   sb = sb * TableStrings.hBlank(30)
    sa = sa * TableStrings.center(16, "omega [Hartree]"; na=4);                     sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.flushleft(40, "Channels: |amplitude|^2"; na=4);          sb = sb * TableStrings.hBlank(44)
    println(stream, sa);   println(stream, sb)
    for  line in lines
        sa  = "  "
        sa  = sa * @sprintf("%14.8f -- %14.8f", line.initialCfLevel.energy, line.finalCfLevel.energy) * "    "
        sa  = sa * @sprintf("%.6e", line.omega) * "    "
        for  channel in line.channels
            sa = sa * "$(channel.multipole)($(channel.gauge)): $(abs(channel.amplitude)^2)   "
        end
        println(stream, sa)
    end
    println(stream, "")
    return  nothing
end


end # module CrystalFieldEmission

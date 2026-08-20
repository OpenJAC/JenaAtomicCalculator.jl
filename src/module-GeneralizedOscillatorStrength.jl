
"""
`module  JAC.GeneralizedOscillatorStrength`
... a submodel of JAC that contains all methods for computing generalized oscillator strengths (GOS) f_n(K) between the levels of
    an initial- and a final-state multiplet.  The GOS is the central quantity of the Bethe theory of fast charged-particle
    collisions, cf. M. Inokuti, Rev. Mod. Phys. 43 (1971) 297; it is the inelastic counterpart of the atomic form factor and it
    reduces to the ordinary (optical) oscillator strength in the limit of vanishing momentum transfer.

    With the plane wave expanded into multipoles, exp(i K.r) = Sum_L i^L (2L+1) j_L(Kr) C^L_0(rhat), and after averaging over the
    magnetic substates of the initial level and summing over those of the final level, the GOS is given (in atomic units) by

        f_n(K) = 2 dE / K^2 * 1/(2J_i+1) * Sum_L (2L+1) |<alpha_f J_f || Sum_j j_L(K r_j) C^L(j) || alpha_i J_i>|^2 .

    The operator j_L(Kr) C^L carries the parity (-1)^L, so that only those L contribute which obey both the triangular rule
    |J_f - J_i| <= L <= J_f + J_i and (-1)^L = pi_i pi_f.  The sum over L is therefore FINITE and EXACT: unlike a multipole
    expansion of the radiation field, it needs no truncation and no convergence parameter.

    This version treats bound-bound (discrete) transitions only.  The density df(K,E)/dE of the generalized oscillator strength
    in the continuum -- and hence the complete Bethe surface, the Bethe sum rule, the first-Born integrated cross sections and
    the stopping power -- is deliberately left out; the data structures of this module carry over to those quantities unchanged.
"""
module GeneralizedOscillatorStrength


using  Printf, ..AngularMomentum, ..Basics, ..BiOrthogonal, ..Defaults, ..ManyElectron, ..Radial, ..RadialIntegrals,
       ..SpinAngular, ..TableStrings


"""
`struct  GeneralizedOscillatorStrength.Settings  <:  AbstractProcessSettings`
    ... defines a type for the details and parameters of computing generalized oscillator strengths.

    + qValues                  ::Array{Float64,1}   ... List of momentum transfers K [a.u.] for which the GOS is computed.
    + calcOpticalLimit         ::Bool               ... True, if the optical limit f_n (K --> 0) is evaluated in addition; this
                                                        value must agree with the (length-form) oscillator strength of the same
                                                        transition and is the natural calibration of the whole computation.
    + calcBiorthogonal         ::Bool               ... True, if the initial- and final-state multiplets are first brought into a
                                                        bi-orthogonal representation (`BiOrthogonal.computeTransformation`) before
                                                        the amplitudes are evaluated, and false if they are used as they are.
    + printBefore              ::Bool               ... True, if all selected lines are printed before their evaluation.
    + lineSelection            ::LineSelection      ... Specifies the selected level pairs, if any.
"""
struct Settings  <:  AbstractProcessSettings
    qValues                    ::Array{Float64,1}
    calcOpticalLimit           ::Bool
    calcBiorthogonal           ::Bool
    printBefore                ::Bool
    lineSelection              ::LineSelection
end


"""
`GeneralizedOscillatorStrength.Settings()`
    ... constructor for an `empty` instance of GeneralizedOscillatorStrength.Settings; a settings::Settings is returned.
"""
function Settings()
    Settings( Float64[0.01, 0.1, 0.5, 1.0, 2.0, 5.0], true, false, false, LineSelection() )
end


"""
`GeneralizedOscillatorStrength.Settings(set::GeneralizedOscillatorStrength.Settings;`

        qValues=..,             calcOpticalLimit=..,        calcBiorthogonal=..,        printBefore=..,
        lineSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::GeneralizedOscillatorStrength.Settings;
        a settings::Settings is returned.
"""
function Settings(set::GeneralizedOscillatorStrength.Settings;
    qValues::Union{Nothing,Array{Float64,1}}=nothing,               calcOpticalLimit::Union{Nothing,Bool}=nothing,
    calcBiorthogonal::Union{Nothing,Bool}=nothing,                  printBefore::Union{Nothing,Bool}=nothing,
    lineSelection::Union{Nothing,LineSelection}=nothing)

    if  isnothing(qValues)            qValuesx          = set.qValues          else   qValuesx          = qValues          end
    if  isnothing(calcOpticalLimit)   calcOpticalLimitx = set.calcOpticalLimit else   calcOpticalLimitx = calcOpticalLimit end
    if  isnothing(calcBiorthogonal)   calcBiorthogonalx = set.calcBiorthogonal else   calcBiorthogonalx = calcBiorthogonal end
    if  isnothing(printBefore)        printBeforex      = set.printBefore      else   printBeforex      = printBefore      end
    if  isnothing(lineSelection)      lineSelectionx    = set.lineSelection    else   lineSelectionx    = lineSelection    end

    Settings( qValuesx, calcOpticalLimitx, calcBiorthogonalx, printBeforex, lineSelectionx )
end


# `Base.show(io::IO, settings::GeneralizedOscillatorStrength.Settings)`  ... prepares a proper printout of settings.
function Base.show(io::IO, settings::GeneralizedOscillatorStrength.Settings)
    println(io, "qValues:                  $(settings.qValues)  ")
    println(io, "calcOpticalLimit:         $(settings.calcOpticalLimit)  ")
    println(io, "calcBiorthogonal:         $(settings.calcBiorthogonal)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
end


"""
`struct  GeneralizedOscillatorStrength.Channel`
    ... defines a type for a single (momentum transfer, multipole rank) component of a generalized oscillator strength.

    + q                        ::Float64      ... Momentum transfer K [a.u.] of this component.
    + L                        ::Int64        ... Rank of the tensor operator j_L(Kr) C^L.
    + amplitude                ::ComplexF64   ... Reduced matrix element <alpha_f J_f || Sum_j j_L(K r_j) C^L(j) || alpha_i J_i>.
"""
struct  Channel
    q                          ::Float64
    L                          ::Int64
    amplitude                  ::ComplexF64
end


"""
`GeneralizedOscillatorStrength.Channel()`
    ... constructor for an `empty` instance of GeneralizedOscillatorStrength.Channel; a channel::Channel is returned.
"""
function Channel()
    Channel( 0., 0, ComplexF64(0.) )
end


# `Base.show(io::IO, channel::GeneralizedOscillatorStrength.Channel)`  ... prepares a proper printout of channel.
function Base.show(io::IO, channel::GeneralizedOscillatorStrength.Channel)
    println(io, "q:                        $(channel.q)  ")
    println(io, "L:                        $(channel.L)  ")
    println(io, "amplitude:                $(channel.amplitude)  ")
end


"""
`struct  GeneralizedOscillatorStrength.Line`
    ... defines a type for a generalized oscillator strength, i.e. for one pair of an initial and a final level, together with the
        GOS values on the whole grid of momentum transfers.

    + initialLevel             ::Level                ... initial-(state) level
    + finalLevel               ::Level                ... final-(state) level
    + deltaEnergy              ::Float64              ... Transition energy E_f - E_i [a.u.] of this line.
    + qValues                  ::Array{Float64,1}     ... Momentum transfers K [a.u.] at which the GOS has been evaluated.
    + gosValues                ::Array{Float64,1}     ... Generalized oscillator strengths f_n(K), parallel to qValues.
    + opticalLimit             ::Float64              ... Optical oscillator strength f_n, i.e. the limit of f_n(K) for K --> 0.
    + channels                 ::Array{GeneralizedOscillatorStrength.Channel,1}   ... List of all (q, L) components of this line.
"""
struct  Line
    initialLevel               ::Level
    finalLevel                 ::Level
    deltaEnergy                ::Float64
    qValues                    ::Array{Float64,1}
    gosValues                  ::Array{Float64,1}
    opticalLimit               ::Float64
    channels                   ::Array{GeneralizedOscillatorStrength.Channel,1}
end


"""
`GeneralizedOscillatorStrength.Line()`
    ... constructor for an `empty` instance of GeneralizedOscillatorStrength.Line; a line::Line is returned.
"""
function Line()
    Line( Level(), Level(), 0., Float64[], Float64[], 0., GeneralizedOscillatorStrength.Channel[] )
end


"""
`GeneralizedOscillatorStrength.Line(initialLevel::Level, finalLevel::Level, deltaEnergy::Float64, qValues::Array{Float64,1})`
    ... constructor for a generalized-oscillator-strength line between a specified initial and final level, and for which no
        amplitude has yet been evaluated; a line::Line is returned.
"""
function Line(initialLevel::Level, finalLevel::Level, deltaEnergy::Float64, qValues::Array{Float64,1})
    Line( initialLevel, finalLevel, deltaEnergy, qValues, zeros(length(qValues)), 0., GeneralizedOscillatorStrength.Channel[] )
end


# `Base.show(io::IO, line::GeneralizedOscillatorStrength.Line)`  ... prepares a proper printout of line.
function Base.show(io::IO, line::GeneralizedOscillatorStrength.Line)
    println(io, "initialLevel:             $(line.initialLevel)  ")
    println(io, "finalLevel:               $(line.finalLevel)  ")
    println(io, "deltaEnergy:              $(line.deltaEnergy)  ")
    println(io, "qValues:                  $(line.qValues)  ")
    println(io, "gosValues:                $(line.gosValues)  ")
    println(io, "opticalLimit:             $(line.opticalLimit)  ")
    println(io, "channels:                 $(line.channels)  ")
end


"""
`GeneralizedOscillatorStrength.amplitude(L::Int64, q::Float64, finalLevel::Level, initialLevel::Level, grid::Radial.Grid;`
                                         `display::Bool=false)`
    ... to compute the (many-electron) reduced matrix element <alpha_f J_f || Sum_j j_L(q r_j) C^L(j) || alpha_i J_i> of the L-th
        multipole component of the plane wave exp(i K.r), for the given final and initial level.  The single-electron kernel is
        the Grant reduced matrix element <kappa_a || C^L || kappa_b> -- which is the same for the large and the small components,
        both requiring l_a + l_b + L to be even -- times the radial integral int dr j_L(qr) [P_a P_b + Q_a Q_b].  The operator
        carries the parity (-1)^L, so that the matrix element vanishes unless (-1)^L = pi_i pi_f.  A value::ComplexF64 is returned.
"""
function amplitude(L::Int64, q::Float64, finalLevel::Level, initialLevel::Level, grid::Radial.Grid; display::Bool=false)
    # The operator j_L(qr) C^L has parity (-1)^L: for even L the two levels must share the same parity, for odd L they must
    # differ.  A blind "parity must differ" guard would be correct for L = 1 only and would let a same-parity pair through.
    if       iseven(L)  &&  finalLevel.parity != initialLevel.parity     return( ComplexF64(0.) )
    elseif   isodd(L)   &&  finalLevel.parity == initialLevel.parity     return( ComplexF64(0.) )
    end
    # Bring both levels into a common subshell basis, if they do not already share one
    if  initialLevel.basis.subshells == finalLevel.basis.subshells
        iLevel = initialLevel;   fLevel = finalLevel
    else
        subshells = Basics.merge(initialLevel.basis.subshells, finalLevel.basis.subshells)
        iLevel    = Level(initialLevel, subshells);     fLevel = Level(finalLevel, subshells)
    end
    #
    nf = length(fLevel.basis.csfs);     ni = length(iLevel.basis.csfs)
    if  display   printstyled("Compute the GOS matrix of dimension $nf x $ni for L = $L, q = $q  ... ", color=:light_green)   end
    matrix = zeros(Float64, nf, ni)
    #
    for  r = 1:nf
        for  s = 1:ni
            if  fLevel.mc[r] == 0.  ||  iLevel.mc[s] == 0.    continue    end

            subshellList = iLevel.basis.subshells
            opa = SpinAngular.OneParticleOperator(L, Basics.plus, true)
            wa  = SpinAngular.computeCoefficients(opa, fLevel.basis.csfs[r], iLevel.basis.csfs[s], subshellList)
            #
            for  coeff in wa
                orba = fLevel.basis.orbitals[coeff.a];      orbb = iLevel.basis.orbitals[coeff.b]
                ja   = Basics.subshell_2j(orba.subshell)
                tamp = AngularMomentum.CL_reduced_me(orba.subshell, L, orbb.subshell) / sqrt( ja + 1.0 ) *
                       RadialIntegrals.GrantJL(L, q, orba, orbb, grid)
                matrix[r,s] = matrix[r,s] + coeff.T * tamp
            end
        end
    end
    if  display   printstyled("done. \n", color=:light_green)   end
    # The spin-angular coefficients of SpinAngular do NOT come in a single normalization: computeCoefficientsScalar (rank 0)
    # divides by sqrt(2J+1) and therefore yields the ORDINARY matrix element, as it must for the one-body Hamiltonian that is
    # its main client, whereas computeCoefficientsNonScalar (rank > 0) divides by sqrt(2L+1) and yields the reduced one up to
    # a factor sqrt(2J_f+1).  Both were calibrated against the exact one-electron reduced matrix element
    # <kappa_f||C^L||kappa_i> * int j_L(qr)[P_f P_i + Q_f Q_i] dr of hydrogen, for J_f = 1/2, 3/2 and 5/2 and for L = 0 ... 4:
    # the ratio is sqrt(2J_f+1) for every L >= 1 and (2J_f+1) for L = 0, with no residual L dependence.
    if  L == 0   wNorm = Basics.twice(fLevel.J) + 1.0
    else         wNorm = sqrt( Basics.twice(fLevel.J) + 1.0 )
    end
    amplitude = ComplexF64( wNorm * transpose(fLevel.mc) * matrix * iLevel.mc )
    #
    if  display
        sa = @sprintf("%.5e", amplitude.re) * "  " * @sprintf("%.5e", amplitude.im)
        println("    < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || j_$L($q r) C^($L) ||" *
                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
    end

    return( amplitude )
end


"""
`GeneralizedOscillatorStrength.computeAmplitudesProperties(line::GeneralizedOscillatorStrength.Line, grid::Radial.Grid,`
                                                          `settings::GeneralizedOscillatorStrength.Settings; printout::Bool=true)`
    ... to compute all reduced matrix elements and the generalized oscillator strengths of the given line, for every momentum
        transfer of the line and every contributing multipole rank L.  A newLine::GeneralizedOscillatorStrength.Line is returned
        for which the channels, the GOS values and the optical limit have now been evaluated.
"""
function computeAmplitudesProperties(line::GeneralizedOscillatorStrength.Line, grid::Radial.Grid,
                                     settings::GeneralizedOscillatorStrength.Settings; printout::Bool=true)
    Ls          = GeneralizedOscillatorStrength.determineRanks(line.finalLevel, line.initialLevel)
    newChannels = GeneralizedOscillatorStrength.Channel[]
    gosValues   = Float64[]
    #
    for  q  in  line.qValues
        qChannels = GeneralizedOscillatorStrength.Channel[]
        for  L  in  Ls
            amp = GeneralizedOscillatorStrength.amplitude(L, q, line.finalLevel, line.initialLevel, grid; display=printout)
            push!( qChannels, GeneralizedOscillatorStrength.Channel(q, L, amp) )
        end
        append!( newChannels, qChannels )
        push!( gosValues, GeneralizedOscillatorStrength.computeGos(q, line.deltaEnergy, line.initialLevel, qChannels) )
    end
    #
    # The optical limit is evaluated at a momentum transfer small enough that j_L(qr) is indistinguishable from its leading term
    # (qr)^L/(2L+1)!! over the whole orbital range, so that only L = 1 survives and f_n(K) has reached its K --> 0 plateau.
    opticalLimit = 0.
    if  settings.calcOpticalLimit
        qSmall    = 1.0e-3
        qChannels = GeneralizedOscillatorStrength.Channel[]
        for  L  in  Ls
            amp = GeneralizedOscillatorStrength.amplitude(L, qSmall, line.finalLevel, line.initialLevel, grid; display=false)
            push!( qChannels, GeneralizedOscillatorStrength.Channel(qSmall, L, amp) )
        end
        opticalLimit = GeneralizedOscillatorStrength.computeGos(qSmall, line.deltaEnergy, line.initialLevel, qChannels)
    end
    newLine = GeneralizedOscillatorStrength.Line( line.initialLevel, line.finalLevel, line.deltaEnergy, line.qValues,
                                                 gosValues, opticalLimit, newChannels )

    return( newLine )
end


"""
`GeneralizedOscillatorStrength.computeGos(q::Float64, deltaEnergy::Float64, initialLevel::Level,`
                                          `channels::Array{GeneralizedOscillatorStrength.Channel,1})`
    ... to combine the reduced matrix elements of all multipole ranks L, evaluated at one and the same momentum transfer q, into
        the generalized oscillator strength f_n(q) = 2 dE/q^2 * 1/(2J_i+1) * Sum_L (2L+1) |<f|| j_L(qr) C^L ||i>|^2, with all
        quantities in atomic units.  A value::Float64 is returned.
"""
function computeGos(q::Float64, deltaEnergy::Float64, initialLevel::Level,
                    channels::Array{GeneralizedOscillatorStrength.Channel,1})
    wa = 0.
    for  ch  in  channels   wa = wa + (2*ch.L + 1) * abs2(ch.amplitude)    end
    gos = 2 * deltaEnergy / (q*q) / (Basics.twice(initialLevel.J) + 1) * wa

    return( gos )
end


"""
`GeneralizedOscillatorStrength.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid,`
                                            `settings::GeneralizedOscillatorStrength.Settings; output=true)`
    ... to compute the generalized oscillator strengths for all selected pairs of levels of the given initial- and final-state
        multiplet, and as specified by the given settings.  A list of lines::Array{GeneralizedOscillatorStrength.Line,1} is
        returned, and the results are printed in neat tables to screen and to the summary file.
"""
function computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid,
                      settings::GeneralizedOscillatorStrength.Settings; output=true)
    if  settings.calcBiorthogonal
        initialMultiplet, finalMultiplet = BiOrthogonal.computeTransformation(initialMultiplet, finalMultiplet, grid)
    end
    println("")
    printstyled("GeneralizedOscillatorStrength.computeLines(): The computation of the GOS starts now ... \n", color=:light_green)
    printstyled("----------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = GeneralizedOscillatorStrength.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    GeneralizedOscillatorStrength.displayLines(stdout, lines)    end
    # Calculate all amplitudes and requested properties
    newLines = GeneralizedOscillatorStrength.Line[]
    for  line in lines
        newLine = GeneralizedOscillatorStrength.computeAmplitudesProperties(line, grid, settings; printout=false)
        push!( newLines, newLine)
    end
    # Print all results to screen
    GeneralizedOscillatorStrength.displayResults(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    GeneralizedOscillatorStrength.displayResults(iostream, newLines, settings)   end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`GeneralizedOscillatorStrength.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet,`
                                              `settings::GeneralizedOscillatorStrength.Settings)`
    ... to determine a list of lines for the computation of the generalized oscillator strengths between the levels of the given
        initial- and final-state multiplet.  A pair of levels is kept only if it is selected by the line selection and if at
        least one multipole rank L contributes.  An Array{GeneralizedOscillatorStrength.Line,1} is returned.
"""
function determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::GeneralizedOscillatorStrength.Settings)
    lines = GeneralizedOscillatorStrength.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                deltaEnergy = fLevel.energy - iLevel.energy
                if  deltaEnergy == 0.                                                          continue   end
                if  length( GeneralizedOscillatorStrength.determineRanks(fLevel, iLevel) ) == 0  continue   end
                push!( lines, GeneralizedOscillatorStrength.Line(iLevel, fLevel, deltaEnergy, settings.qValues) )
            end
        end
    end

    return( lines )
end


"""
`GeneralizedOscillatorStrength.determineRanks(finalLevel::Level, initialLevel::Level)`
    ... to determine all multipole ranks L that contribute to the generalized oscillator strength between the given initial and
        final level.  Two rules apply and, together, make the list finite and complete: the triangular rule
        |J_f - J_i| <= L <= J_f + J_i, and the parity of the operator j_L(Kr) C^L, which requires (-1)^L = pi_i pi_f.  No
        truncation and no convergence parameter is involved.  An Array{Int64,1} is returned.

        The monopole term L = 0 is INCLUDED and is not a trivial one: j_0(Kr) = sin(Kr)/(Kr) is a genuine operator that connects
        two different levels of the same J and parity, and it is what carries, for instance, the whole 1s --> 2s generalized
        oscillator strength of atomic hydrogen.  Only in the limit K --> 0 does j_0(Kr) approach unity, whereupon the matrix
        element collapses to the overlap of two orthogonal states -- which is precisely why f_n(K --> 0) vanishes for a
        monopole-only transition.
"""
function determineRanks(finalLevel::Level, initialLevel::Level)
    Ls    = Int64[]
    Ji2   = Basics.twice(initialLevel.J);    Jf2 = Basics.twice(finalLevel.J)
    Lmin  = div( abs(Jf2 - Ji2), 2 );        Lmax = div( Jf2 + Ji2, 2 )
    same  = (finalLevel.parity == initialLevel.parity)
    for  L = Lmin:Lmax
        if  iseven(L) == same   push!( Ls, L )    end
    end

    return( Ls )
end


"""
`GeneralizedOscillatorStrength.displayLines(stream::IO, lines::Array{GeneralizedOscillatorStrength.Line,1})`
    ... to display a list of the level pairs that have been selected for the computations, together with their transition energy
        and the multipole ranks L that will contribute.  A neat table is printed but nothing is returned otherwise.
"""
function displayLines(stream::IO, lines::Array{GeneralizedOscillatorStrength.Line,1})
    nx = 108
    println(stream, " ")
    println(stream, "  Selected generalized-oscillator-strength lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Energy"; na=4)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(20, "Multipole ranks L"; na=2);              sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(20, "No. of q-values"; na=2);                sb = sb * TableStrings.hBlank(22)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  line in lines
        sa   = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                      fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.deltaEnergy)) * "    "
        Ls = GeneralizedOscillatorStrength.determineRanks(line.finalLevel, line.initialLevel)
        sa = sa * TableStrings.flushleft(20, string(Ls); na=2)
        sa = sa * TableStrings.flushleft(20, string(length(line.qValues)); na=2)
        println(stream,  sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`GeneralizedOscillatorStrength.displayResults(stream::IO, lines::Array{GeneralizedOscillatorStrength.Line,1},`
                                              `settings::GeneralizedOscillatorStrength.Settings)`
    ... to display the generalized oscillator strengths of all lines, one block per line: the momentum transfer both as K and as
        the variable (Ka_0)^2 in which the Bethe theory is usually plotted, the total f_n(K), and the separate contribution of
        each multipole rank L.  The optical limit is appended if it has been computed.  A neat table is printed but nothing is
        returned otherwise.
"""
function displayResults(stream::IO, lines::Array{GeneralizedOscillatorStrength.Line,1},
                        settings::GeneralizedOscillatorStrength.Settings)
    nx = 116
    println(stream, " ")
    println(stream, "  Generalized oscillator strengths f_n(K)  [f_n(K --> 0) is the optical oscillator strength]:")
    println(stream, " ")
    #
    for  line  in  lines
        isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        Ls   = GeneralizedOscillatorStrength.determineRanks(line.finalLevel, line.initialLevel)
        sa   = "  Line " * TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index) * "   " *
               TableStrings.symmetries_if(isym, fsym) * "    dE = " *
               @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.deltaEnergy)) * " " *
               TableStrings.inUnits("energy") * "    ranks L = " * string(Ls)
        println(stream, " ");    println(stream, sa)
        println(stream, "  ", TableStrings.hLine(nx))
        sb = "     " * TableStrings.center(16, "K [a.u.]"; na=2) * TableStrings.center(16, "(K a_0)^2"; na=2) *
             TableStrings.center(16, "f_n(K)"; na=2)
        for  L in Ls   sb = sb * TableStrings.center(16, "f_n(K), L=$L"; na=2)   end
        println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
        #
        for  (i, q)  in  enumerate(line.qValues)
            sc = "     " * @sprintf("%16.6e", q) * "  " * @sprintf("%16.6e", q*q) * "  " *
                 @sprintf("%16.6e", line.gosValues[i]) * "  "
            for  L in Ls
                qChannels = filter( ch -> ch.q == q  &&  ch.L == L, line.channels )
                gosL      = GeneralizedOscillatorStrength.computeGos(q, line.deltaEnergy, line.initialLevel, qChannels)
                sc        = sc * @sprintf("%16.6e", gosL) * "  "
            end
            println(stream, sc)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        if  settings.calcOpticalLimit
            println(stream, "     optical limit  f_n (K --> 0)  = " * @sprintf("%.6e", line.opticalLimit))
        end
    end
    println(stream, " ")

    return( nothing )
end

end # module

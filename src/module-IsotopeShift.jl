
"""
`module  JAC.IsotopeShift`  
... a submodel of JAC that contains all methods for computing isotope-dependent properties for some level(s).
"""
module IsotopeShift


using Printf, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals,
                ..SpinAngular, ..TableStrings

"""
`struct  IsotopeShift.Settings  <:  AbstractPropertySettings`  
    ... defines a type for the details and parameters of computing isotope-shift M and F parameters.

    + calcNMS                  ::Bool             ... True if mass-shift parameters M_nmn need to be calculated, and false otherwise.
    + calcSMS                  ::Bool             ... True if mass-shift parameters M_sms need to be calculated, and false otherwise.
    + calcF                    ::Bool             ... True if the field-shift parameter need to be calculated, and false otherwise.
    + calcBoson                ::Bool             ... True if the boson-field shift parameter need to be calculated, and false otherwise.
    + printBefore              ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
    + bosonMass                ::Float64          ... mass of the scalar boson [e_electron].
    + levelSelection           ::LevelSelection   ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings 
    calcNMS                    ::Bool
    calcSMS                    ::Bool
    calcF                      ::Bool
    calcBoson                  ::Bool
    printBefore                ::Bool
    bosonMass                  ::Float64
    levelSelection             ::LevelSelection
end 


"""
`IsotopeShift.Settings(; calcNMS::Bool=true,` calcSMS::Bool=false, calcF::Bool=false, calcBoson::Bool=false, 
                            printBefore::Bool=true, bosonMass::Float64=0., levelSelection::LevelSelection=LevelSelection()) 
    ... keyword constructor to overwrite selected value of isoshift computations.
"""
function Settings(; calcNMS::Bool=true, calcSMS::Bool=false, calcF::Bool=false, calcBoson::Bool=false, 
                            printBefore::Bool=true, bosonMass::Float64=0., levelSelection::LevelSelection=LevelSelection())
    Settings(calcNMS, calcSMS, calcF, calcBoson, printBefore, bosonMass, levelSelection)
end


# `Base.show(io::IO, settings::IsotopeShift.Settings)`  ... prepares a proper printout of the variable settings::IsotopeShift.Settings.
function Base.show(io::IO, settings::IsotopeShift.Settings) 
    println(io, "calcNMS:                  $(settings.calcNMS)  ")
    println(io, "calcSMS:                  $(settings.calcSMS)  ")
    println(io, "calcF:                    $(settings.calcF)  ")
    println(io, "calcBoson:                $(settings.calcBoson)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "bosonMass:                $(settings.bosonMass)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end



"""
`struct  IsotopeShift.Outcome`  
    ... defines a type to keep the outcome of a isotope-shift computation, such as the K and F parameters as well other 
        results.

    + level                     ::Level              ... Atomic level to which the outcome refers to.
    + Knms                      ::Float64            ... K_nms parameter
    + Ksms                      ::Float64            ... K_sms parameter
    + Fme                       ::Float64            ... F parameter from matrix element.
    + Fdensity                  ::Float64            ... F parameter from density.
    + Xboson                    ::Float64            ... X boson-field shift constant.
    + amplitudeKnms             ::Complex{Float64}   ... K_nms amplitude
    + amplitudeKsmsA            ::Complex{Float64}   ... K_sms,A amplitude
    + amplitudeKsmsB            ::Complex{Float64}   ... K_sms,B amplitude
    + amplitudeKsmsC            ::Complex{Float64}   ... K_sms,C amplitude
"""
struct Outcome 
    level                       ::Level 
    Knms                        ::Float64
    Ksms                        ::Float64
    Fme                         ::Float64
    Fdensity                    ::Float64
    Xboson                      ::Float64
    amplitudeKnms               ::Complex{Float64}
    amplitudeKsmsA              ::Complex{Float64}
    amplitudeKsmsB              ::Complex{Float64}
    amplitudeKsmsC              ::Complex{Float64}
end 


"""
`IsotopeShift.Outcome()`  ... constructor for an `empty` instance of Hfs.Outcome for the computation of isotope-shift properties.
"""
function Outcome()
    Outcome(Level(), 0., 0., 0., 0., 0.,   0., 0., 0., 0.)
end


# `Base.show(io::IO, outcome::IsotopeShift.Outcome)`  ... prepares a proper printout of the variable outcome::IsotopeShift.Outcome.
function Base.show(io::IO, outcome::IsotopeShift.Outcome) 
    println(io, "level:                   $(outcome.level)  ")
    println(io, "Knms:                    $(outcome.Knms)  ")
    println(io, "Ksms:                    $(outcome.Ksms)  ")
    println(io, "Fme:                     $(outcome.Fme)  ")
    println(io, "Fdensity:                $(outcome.Fdensity)  ")
    println(io, "Xboson:                  $(outcome.Xboson)  ")
    println(io, "amplitudeKnms:           $(outcome.amplitudeKnms)  ")
    println(io, "amplitudeKsmsA:          $(outcome.amplitudeKsmsA)  ")
    println(io, "amplitudeKsmsB:          $(outcome.amplitudeKsmsB)  ")
    println(io, "amplitudeKsmsC:          $(outcome.amplitudeKsmsC)  ")
end


"""
`IsotopeShift.amplitude(kind::String, rLevel::Level, sLevel::Level, nm::Nuclear.Model, grid::Radial.Grid)` 
    ... to compute either the  H^(NMS),  H^(SMS,A),  H^(SMS,B) or  H^(SMS,C) normal and specific mass-shift amplitudes
        <alpha_r J_r || H^(A)) || alpha_s J_s>  for a given pair of levels. A value::ComplexF64 is returned.
"""
function  amplitude(kind::String, rLevel::Level, sLevel::Level, nm::Nuclear.Model, grid::Radial.Grid)
    nr = length(rLevel.basis.csfs);    ns = length(sLevel.basis.csfs);    matrix = zeros(ComplexF64, nr, ns)
    printstyled("Compute mass-shift $(kind[1:9]) matrix of dimension $nr x $ns ... ", color=:light_green)
    #
    if  rLevel.parity != sLevel.parity   return( ComplexF64(0.) )   end
    #
    for  r = 1:nr
        for  s = 1:ns
            me = 0.
            if  rLevel.basis.csfs[r].parity  != rLevel.parity    ||  sLevel.basis.csfs[s].parity  != sLevel.parity  ||
                rLevel.parity != sLevel.parity    continue    
            end 
            #
            if      kind == "H^(NMS)   amplitude"
            #------------------------------------
                subshellList = rLevel.basis.subshells
                op = SpinAngular.OneParticleOperator(0, plus, true)
                wa = SpinAngular.computeCoefficients(op, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList)
                for  coeff in wa
                    ja = Basics.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                    jb = Basics.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                    tamp  = InteractionStrength.hamiltonian_nms(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], nm, grid)
                    me = me + coeff.T  * sqrt( ja + 1) * tamp  
                end 
            #
            elseif  kind == "H^(SMS,A) amplitude"
            #------------------------------------
                subshellList = rLevel.basis.subshells
                op = SpinAngular.TwoParticleOperator(0, plus, true)
                wa = SpinAngular.computeCoefficients(op, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList)
                for  coeff in wa
                    if  coeff.nu != 1  continue   end
                    ja = Basics.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                    jb = Basics.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                    tamp  = InteractionStrength.X_smsA(rLevel.basis.orbitals[coeff.a], rLevel.basis.orbitals[coeff.b],
                                                        sLevel.basis.orbitals[coeff.c], sLevel.basis.orbitals[coeff.d], nm, grid)
                    me = me + coeff.V * sqrt( ja + 1) * tamp  
                end
            #
            elseif  kind == "H^(SMS,B) amplitude"
            #------------------------------------
                subshellList = rLevel.basis.subshells
                op = SpinAngular.TwoParticleOperator(0, plus, true)
                wa = SpinAngular.computeCoefficients(op, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList)
                for  coeff in wa
                    if  coeff.nu != 1  continue   end
                    ja = Basics.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                    jb = Basics.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                    tamp  = InteractionStrength.X_smsB(rLevel.basis.orbitals[coeff.a], rLevel.basis.orbitals[coeff.b],
                                                        sLevel.basis.orbitals[coeff.c], sLevel.basis.orbitals[coeff.d], nm, grid)
                    me = me + coeff.V * sqrt( ja + 1) * tamp  
                end
            #
            elseif  kind == "H^(SMS,C) amplitude"
            #------------------------------------
                subshellList = rLevel.basis.subshells
                op = SpinAngular.TwoParticleOperator(0, plus, true)
                wa = SpinAngular.computeCoefficients(op, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList)
                for  coeff in wa
                    if  coeff.nu != 1  continue   end
                    ja = Basics.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                    jb = Basics.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                    tamp  = InteractionStrength.X_smsC(rLevel.basis.orbitals[coeff.a], rLevel.basis.orbitals[coeff.b],
                                                        sLevel.basis.orbitals[coeff.c], sLevel.basis.orbitals[coeff.d], nm, grid)
                    me = me + coeff.V * sqrt( ja + 1) * tamp  
                end
            #
            else    error("stop a")
            end
            #
            matrix[r,s] = me
        end
    end
    printstyled("done.\n", color=:light_green)
    amplitude = transpose(rLevel.mc) * matrix * sLevel.mc 
    #
    return( amplitude )
end


"""
`IsotopeShift.amplitude(kind::String, rLevel::Level, sLevel::Level, potential::Array{Float64,1}, grid::Radial.Grid)` 
    ... to compute the H^(field-shift) field-shift amplitude  <alpha_r J_r || H^(field-shift) || alpha_s J_s> or the
        H^(boson-field) shift amplitude  <alpha_r J_r || H^(boson-field) || alpha_s J_s>  for a given pair of levels.
        The potential has to provide delta V^(nuc) for the field-shift amplitudes and the effective potential for the
        boson-field shift. A value::ComplexF64 is returned.
"""
function  amplitude(kind::String, rLevel::Level, sLevel::Level, potential::Array{Float64,1}, grid::Radial.Grid)
    nr = length(rLevel.basis.csfs);    ns = length(sLevel.basis.csfs);    matrix = zeros(ComplexF64, nr, ns)
    printstyled("Compute field-shift $kind matrix of dimension $nr x $ns ... ", color=:light_green)
    #
    if  rLevel.parity != sLevel.parity   return( ComplexF64(0.) )   end
    #
    for  r = 1:nr
        for  s = 1:ns
            me = 0.
            if  rLevel.basis.csfs[r].parity  != rLevel.parity    ||  sLevel.basis.csfs[s].parity  != sLevel.parity  ||
                rLevel.parity != sLevel.parity    continue    
            end 
            #
            if      kind == "H^(field-shift) amplitude"
            #------------------------------------------
                subshellList = rLevel.basis.subshells
                op = SpinAngular.OneParticleOperator(0, plus, true)
                wa = SpinAngular.computeCoefficients(op, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList)
                for  coeff in wa
                    ja = Basics.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                    jb = Basics.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                    tamp  = InteractionStrength.fieldShift(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], potential, grid)
                    me = me + coeff.T  * sqrt( ja + 1) * tamp  
                end
            #
            elseif  kind == "H^(boson-field) amplitude"
            #------------------------------------------
                subshellList = rLevel.basis.subshells
                op = SpinAngular.OneParticleOperator(0, plus, true)
                wa = SpinAngular.computeCoefficients(op, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList)
                for  coeff in wa
                    ja = Basics.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                    jb = Basics.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                    tamp  = InteractionStrength.bosonShift(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], potential, grid)
                    me = me + coeff.T  * sqrt( ja + 1) * tamp  
                end
            #
            else    error("stop b")
            end
            #
            matrix[r,s] = me
        end
    end
    printstyled("done.\n", color=:light_green)
    amplitude = transpose(rLevel.mc) * matrix * sLevel.mc 
    #
    return( amplitude )
end



"""
`IsotopeShift.computeAmplitudesProperties(outcome::IsotopeShift.Outcome, nm::Nuclear.Model, grid::Radial.Grid, 
                                            settings::IsotopeShift.Settings) 
    ... to compute all amplitudes and properties for a given level; an outcome::IsotopeShift.Outcome is returned for 
        which the amplitudes and properties are now evaluated explicitly.
"""
function  computeAmplitudesProperties(outcome::IsotopeShift.Outcome, nm::Nuclear.Model, grid::Radial.Grid, settings::IsotopeShift.Settings)
    amplitudeKnms   = Ksms = F = Xboson = 0.
    amplitudeKsmsA  = amplitudeKsmsB = amplitudeKsmsC = amplitudeF = 0.
    if  settings.calcNMS
        amplitudeKnms  = IsotopeShift.amplitude("H^(NMS)   amplitude", outcome.level, outcome.level, nm, grid)
    end
    
    if  settings.calcSMS
        amplitudeKsmsA = IsotopeShift.amplitude("H^(SMS,A) amplitude", outcome.level, outcome.level, nm, grid)
        amplitudeKsmsB = IsotopeShift.amplitude("H^(SMS,B) amplitude", outcome.level, outcome.level, nm, grid)
        amplitudeKsmsC = IsotopeShift.amplitude("H^(SMS,C) amplitude", outcome.level, outcome.level, nm, grid)
        Ksms = (amplitudeKsmsA + amplitudeKsmsB + amplitudeKsmsC)
    end
    
    if  settings.calcF
        nmp       = Nuclear.Model(nm.Z, nm.mass+1.0)
        deltaPot  = Nuclear.nuclearPotential(nm, grid).Zr - Nuclear.nuclearPotential(nmp, grid).Zr
        deltaPot  = deltaPot / (nm.radius^2 - nmp.radius^2)
        lastPoint = Basics.lastPoint(deltaPot, 1.0e-9);   @show lastPoint
        #
        Fme       = IsotopeShift.amplitude("H^(field-shift) amplitude", outcome.level, outcome.level, deltaPot, grid)
        #
        density   = Basics.computeDensity(outcome.level, grid)
        wb        = zeros( lastPoint )
        wb[1:lastPoint] = - density[1:lastPoint] .* deltaPot[1:lastPoint] ./ grid.r[1:lastPoint]
        Fdensity  = RadialIntegrals.V0(wb, lastPoint, grid)
    end
    
    if  settings.calcBoson
        potential = zeros( grid.NoPoints );   malpha = settings.bosonMass / Defaults.getDefaults("alpha")
        for     i = 1:length(potential)    potential[i] = exp( -malpha*grid.r[i] ) / grid.r[i]    end
        Xboson    = IsotopeShift.amplitude("H^(boson-field) amplitude", outcome.level, outcome.level, potential, grid)
    end
    
    newOutcome = IsotopeShift.Outcome( outcome.level, amplitudeKnms, Ksms, Fme, Fdensity, Xboson,
                                        amplitudeKnms, amplitudeKsmsA, amplitudeKsmsB, amplitudeKsmsC )
    return( newOutcome )
end


"""
`IsotopeShift.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::IsotopeShift.Settings; output=true)` 
    ... to compute (as selected) the isotope-shift M and F parameters for the levels of the given multiplet and as 
        specified by the given settings. The results are printed in neat tables to screen but nothing is returned otherwise.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::IsotopeShift.Settings; output=true)
    println("")
    printstyled("IsotopeShift.computeOutcomes(): The computation of the isotope-shift parameters starts now ... \n", color=:light_green)
    printstyled("-------------------------------------------------------------------------------------------------- \n", color=:light_green)
    #
    outcomes = IsotopeShift.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    IsotopeShift.displayOutcomes(outcomes)    end
    # Calculate all amplitudes and requested properties
    newOutcomes = IsotopeShift.Outcome[]
    for  outcome in outcomes
        newOutcome = IsotopeShift.computeAmplitudesProperties(outcome, nm, grid, settings) 
        push!( newOutcomes, newOutcome)
    end
    # Print all results to screen
    IsotopeShift.displayResults(stdout, newOutcomes, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    IsotopeShift.displayResults(iostream, newOutcomes, settings)   end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`IsotopeShift.determineOutcomes(multiplet::Multiplet, settings::IsotopeShift.Settings)`  
    ... to determine a list of Outcomes's for the computation of the isotope-shift M and F parameters for the given multiplet. 
        It takes into account the particular selections and settings. An Array{IsotopeShift.Outcome,1} is returned. Apart from 
        the level specification, all physical properties are set to zero during the initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::IsotopeShift.Settings) 
    outcomes = IsotopeShift.Outcome[]
    @show settings.levelSelection
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            push!( outcomes, IsotopeShift.Outcome(level, 0., 0., 0., 0., 0.,   0., 0., 0., 0.) )
        end
    end
    return( outcomes )
end


"""
`IsotopeShift.displayOutcomes(outcomes::Array{IsotopeShift.Outcome,1})`  
    ... to display a list of levels that have been selected for the computations. A small neat table of all 
        selected levels and their energies is printed but nothing is returned otherwise.
"""
function  displayOutcomes(outcomes::Array{IsotopeShift.Outcome,1})
    nx = 43
    println(" ")
    println("  Selected IsotopeShift levels:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = Basics.LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`IsotopeShift.displayResults(stream::IO, outcomes::Array{Hfs.Outcome,1}, settings::IsotopeShift.Settings)`  
    ... to display the energies, M_ms and F-parameters, etc. for the selected levels. A neat table is printed but nothing 
        is returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{IsotopeShift.Outcome,1}, settings::IsotopeShift.Settings)
    nx = 133
    println(stream, " ")
    println(stream, "  IsotopeShift parameters and amplitudes:  Just take the difference K_i - K_f for transition lines.")
    println(stream, " ")
    println(stream, "  Boson mass = $(settings.bosonMass)")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  ";   sc = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=5)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5);  sc = sc * TableStrings.hBlank(45)
    sa = sa * TableStrings.center(11, "K_nms"; na=4)              
    sb = sb * TableStrings.center(11, "[GHz u]" ; na=4)
    sc = sc * TableStrings.center(11, "[a.u.]" ; na=4)
    sa = sa * TableStrings.center(11, "K_sms"; na=4)             
    sb = sb * TableStrings.center(11, "[GHz u]" ; na=4)
    sc = sc * TableStrings.center(11, "[a.u.]" ; na=4)
    sa = sa * TableStrings.center(11, "K_ms "; na=4)              
    sb = sb * TableStrings.center(11, "[GHz u]" ; na=4)
    sc = sc * TableStrings.center(11, "[a.u.]" ; na=4)
    sa = sa * TableStrings.center(11, " F [ME]";   na=4)              
    sb = sb * TableStrings.center(11, "[MHz/fm^2]" ; na=4)
    sa = sa * TableStrings.center(11, " F [dens]";   na=5)              
    sb = sb * TableStrings.center(11, "[MHz/fm^2]" ; na=5)
    sa = sa * TableStrings.center(11, "X^(boson)";   na=4)              
    sb = sb * TableStrings.center(11, "a.u."; na=4)
    println(stream, sa);    println(stream, sb);    ## println(stream, sc);    
    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy))            * "    "
        sa = sa * @sprintf("%.5e", outcome.Knms * 3609.4824)                                                      * "    "
        sa = sa * @sprintf("%.5e", outcome.Ksms * 3609.4824)                                                      * "    "
        sa = sa * @sprintf("%.5e", (outcome.Knms + outcome.Ksms) * 3609.4824)                                     * "    "
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic to Hz", outcome.Fme)      * 1.0e-6) * "    " 
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic to Hz", outcome.Fdensity) * 1.0e-6) * "    "
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic to Hz", outcome.Xboson) )      
        println(stream, sa )
        ## sb = TableStrings.hBlank(47)
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n\n")
    #
    #
    # Now printout the individual amplitudes in an extra table
    nx = 158
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(26, "H^(nms) amplitude"    ; na=3);             sb = sb * TableStrings.hBlank(35)
    sa = sa * TableStrings.center(26, "H^(sms,A) amplitude"  ; na=3);             sb = sb * TableStrings.hBlank(35)
    sa = sa * TableStrings.center(26, "H^(sms,B) amplitude"  ; na=3);             sb = sb * TableStrings.hBlank(35)
    sa = sa * TableStrings.center(26, "H^(sms,C) amplitude"  ; na=3);             sb = sb * TableStrings.hBlank(35)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy))          * "    "
        sa = sa * @sprintf("%.5e %s %.5e", outcome.amplitudeKnms.re,  " ", outcome.amplitudeKnms.im)  * "    "
        sa = sa * @sprintf("%.5e %s %.5e", outcome.amplitudeKsmsA.re, " ", outcome.amplitudeKsmsA.im) * "    "
        sa = sa * @sprintf("%.5e %s %.5e", outcome.amplitudeKsmsB.re, " ", outcome.amplitudeKsmsB.im) * "    "
        sa = sa * @sprintf("%.5e %s %.5e", outcome.amplitudeKsmsC.re, " ", outcome.amplitudeKsmsC.im) * "    "
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #
    # Now printout the transition isotope parameters
    nx = 144
    println(stream, " ")
    println(stream, "  IsotopeShift parameters for individual transitions:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  ";   sc = "  "
    sa = sa * TableStrings.center(14, "Transition";    na=2);                        sb = sb * TableStrings.hBlank(16)
    sa = sa * TableStrings.center(14, "I -- J^P -- F"; na=4);                        sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(14, "Energy"; na=5)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5);  sc = sc * TableStrings.hBlank(45)
    sa = sa * TableStrings.center(11, "K_nms"; na=5)              
    sb = sb * TableStrings.center(11, "[GHz u]" ; na=5)
    sc = sc * TableStrings.center(11, "[a.u.]" ; na=5)
    sa = sa * TableStrings.center(11, "K_sms"; na=5)             
    sb = sb * TableStrings.center(11, "[GHz u]" ; na=5)
    sc = sc * TableStrings.center(11, "[a.u.]" ; na=5)
    sa = sa * TableStrings.center(11, "K_ms "; na=5)              
    sb = sb * TableStrings.center(11, "[GHz u]" ; na=5)
    sc = sc * TableStrings.center(11, "[a.u.]" ; na=5)
    sa = sa * TableStrings.center(11, " F [ME]";   na=4)              
    sb = sb * TableStrings.center(11, "[MHz/fm^2]" ; na=4)
    sa = sa * TableStrings.center(11, " F [dens]";   na=4)              
    sb = sb * TableStrings.center(11, "[MHz/fm^2]" ; na=4)
    sa = sa * TableStrings.center(11, "X^(boson)";   na=4)              
    sb = sb * TableStrings.center(11, "a.u."; na=4)
    println(stream, sa);    println(stream, sb);    ## println(stream, sc);    
    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  ioutcome in outcomes
        for  foutcome in outcomes
        if ioutcome.level.index >= foutcome.level.index     continue    end
        sa  = "";    isym = LevelSymmetry( ioutcome.level.J, ioutcome.level.parity) 
                        fsym = LevelSymmetry( foutcome.level.J, foutcome.level.parity)
        sa = sa * TableStrings.center(14, TableStrings.levels_if(ioutcome.level.index, foutcome.level.index); na=2)
        sa = sa * TableStrings.center(18, string(isym) * "  --   " * string(fsym); na=2)
        en = abs(ioutcome.level.energy - foutcome.level.energy)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", en))                              * "     "
        sa = sa * @sprintf("%.5e", (ioutcome.Knms-foutcome.Knms) * 3609.4824)                                     * "    "
        sa = sa * @sprintf("%.5e", (ioutcome.Ksms-foutcome.Ksms) * 3609.4824)                                     * "     "
        sa = sa * @sprintf("%.5e", (ioutcome.Knms + ioutcome.Ksms - foutcome.Knms - foutcome.Ksms) * 3609.4824)   * "     "
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic to Hz", ioutcome.Fme-foutcome.Fme)  * 1.0e-6) * "    "      
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic to Hz", ioutcome.Fdensity-foutcome.Fdensity)  * 1.0e-6) * "    "      
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic to Hz", ioutcome.Xboson-foutcome.Xboson))      
        println(stream, sa )
        ## sb = TableStrings.hBlank(47)
        end
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n\n")

    return( nothing )
end

end # module

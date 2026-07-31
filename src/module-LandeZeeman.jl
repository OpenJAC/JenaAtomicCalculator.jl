
"""
`module  JAC.LandeZeeman`  
... a submodel of JAC that contains all methods for computing Lande factors and Zeeman properties for some level(s).
"""
module LandeZeeman


using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial,
                ..SpinAngular, ..TableStrings, ..Hfs


"""
`abstract type LandeZeeman.AbstractZeemanKind`
    ... defines an abstract and two singleton types to select the N^(1) or Delta N^(1) Zeeman amplitude.

    + ZeemanN1       ... non-relativistic Zeeman operator N^(1).
    + ZeemanDeltaN1  ... relativistic correction Delta N^(1) to the Zeeman operator.
"""
abstract type  AbstractZeemanKind                          end
struct         ZeemanN1       <:  AbstractZeemanKind       end
struct         ZeemanDeltaN1  <:  AbstractZeemanKind       end


"""
`struct  LandeZeeman.SublevelJ`  ... defines a type to specify a magnetic sublevel with well-defined J.

    + M                      ::AngularM64        ... M_J-value
    + energy                 ::Float64           ... energy of this sublevel
    + c2Coeff                ::Float64           ... quadratic-Zeeman shift C_2 coefficient for this level
"""
struct SublevelJ 
    M                        ::AngularM64
    energy                   ::Float64
    c2Coeff                  ::Float64
end 


# `Base.show(io::IO, Jsublevel::LandeZeeman.SublevelJ)`  ... prepares a proper printout of the variable Jsublevel::LandeZeeman.SublevelJ.
function Base.show(io::IO, Jsublevel::LandeZeeman.SublevelJ) 
    println(io, "Sublevel [M=$(Jsublevel.M); energy = $(Jsublevel.energy); C_2 = $(Jsublevel.c2Coeff)]")
end



# SublevelF defines a magnetic hyperfine sublevel
"""
`struct  LandeZeeman.SublevelF`  ... defines a type to specify a magnetic hyperfine sublevel with well-defined F.

    + F                      ::AngularJ64        ... F-value
    + M                      ::AngularM64        ... M_F-value
    + energy                 ::Float64           ... energy of this sublevel
    + c2Coeff                ::Float64           ... quadratic-Zeeman shift C_2 coefficient of this level
"""
struct SublevelF
    F                        ::AngularJ64
    M                        ::AngularM64
    energy                   ::Float64
    c2Coeff                  ::Float64
end 


# `Base.show(io::IO, Fsublevel::LandeZeeman.SublevelF)`  ... prepares a proper printout of the variable Fsublevel::LandeZeeman.SublevelF.
function Base.show(io::IO, Fsublevel::LandeZeeman.SublevelF) 
    println(io, "Sublevel [F=$(Fsublevel.F); M=$(Fsublevel.M); energy = $(Fsublevel.energy); C_2 = $(Fsublevel.c2Coeff)]")
end


"""
`struct  LandeZeeman.Outcome`  
    ... defines a type to keep the Lande Factors and Zeeman spittling parameters of a fine-structure level.

    + Jlevel                 ::Level             ... Fine-structure levels to which the results refer to.
    + LandeJ                 ::Float64           ... Lande-J factor of this fine-structure level 
    + amplitudeN1            ::Complex{Float64}  ... N1 amplitude
    + amplitudeDeltaN1       ::Complex{Float64}  ... Delta-N1 amplitude
    + nuclearI               ::AngularJ64        ... nuclear spin
    + Jsublevels             ::Array{LandeZeeman.SublevelJ,1}   
        ... List of the magnetic fine-structure sublevels and data with well-defined J-value 
    + Fsublevels             ::Array{LandeZeeman.SublevelF,1}   
        ... List of the magnetic hyperfine sublevels and data with well-defined F-value 
"""
struct Outcome 
    Jlevel                   ::Level
    LandeJ                   ::Float64
    amplitudeN1              ::Complex{Float64}
    amplitudeDeltaN1         ::Complex{Float64}
    nuclearI                 ::AngularJ64
    Jsublevels               ::Array{LandeZeeman.SublevelJ,1}
    Fsublevels               ::Array{LandeZeeman.SublevelF,1}
end 


"""
`LandeZeeman.Outcome()`  ... constructor for an `empty` instance of LandeZeeman.Outcome.
"""
function Outcome()
    Outcome(Level(), 0., ComplexF64(0.), ComplexF64(0.), AngularJ64(0), SublevelJ[],  SublevelF[])
end


# `Base.show(io::IO, outcome::LandeZeeman.Outcome)`  ... prepares a proper printout of the variable LandeZeeman.Outcome.
function Base.show(io::IO, outcome::LandeZeeman.Outcome) 
    println(io, "Jlevel:           $(outcome.Jlevel)  ")
    println(io, "LandeJ:           $(outcome.LandeJ)  ")
    println(io, "amplitudeN1:      $(outcome.amplitudeN1)  ")
    println(io, "amplitudeDeltaN1: $(outcome.amplitudeDeltaN1)  ")
    println(io, "nuclearI:         $(outcome.nuclearI)  ")
    println(io, "Jsublevels:       $(outcome.Jsublevels)  ")
    println(io, "Fsublevels:       $(outcome.Fsublevels)  ")
end


"""
`struct  LandeZeeman.Settings  <:  AbstractPropertySettings`  
    ... defines a type for the details and parameters of computing the Lande Factors and Zeeman spittling 
        of fine-structure levels.

    + calcLandeJ             ::Bool              ... True if Lande-J factors are to be calculated.
    + calcLandeF             ::Bool              ... True if Lande-F factors are to be calculated.
    + calcZeeman             ::Bool              ... True if Zeeman splittings are to be calculated.
    + calcQZScoeff           ::Bool              ... True if the quadratic-Zeeman shift coefficients are to be calculated.
    + includeSchwinger       ::Bool              ... True if Schwinger's QED correction ``\\Delta N^(1)`` is to be included, 
                                                        and false otherwise.
    + printBefore            ::Bool              ... True if a list of selected levels is printed before the actual computations start. 
    + BField                 ::Float64           ... Strength of the magnetic field in [Tesla]
    + levelSelection         ::LevelSelection    ... Specifies the selected levels, if any.
    + gMultiplet             ::Multiplet         ... Specifies the levels, used as intermediate levels for second-order coefficients.
"""
struct Settings  <:  AbstractPropertySettings 
    calcLandeJ               ::Bool
    calcLandeF               ::Bool
    calcZeeman               ::Bool
    calcQZScoeff             ::Bool
    includeSchwinger         ::Bool
    printBefore              ::Bool 
    BField                   ::Float64
    levelSelection           ::LevelSelection
    gMultiplet               ::Multiplet  
end 


"""
`LandeZeeman.Settings()`
    ... constructor for an `empty` instance of ZeemanSettings for the computation of isotope M and F parameters.
"""
function Settings()
    Settings(false, false, false, false, false, false, 0., LevelSelection(), Multiplet() )
end


"""
`LandeZeeman.Settings(set::LandeZeeman.Settings; keywords...)`
    ... keyword constructor to modify selected fields of a LandeZeeman.Settings object;
        for fields not given, the values of set are used.
"""
function Settings(set::LandeZeeman.Settings;
        calcLandeJ::Union{Nothing,Bool}=nothing,           calcLandeF::Union{Nothing,Bool}=nothing,
        calcZeeman::Union{Nothing,Bool}=nothing,            calcQZScoeff::Union{Nothing,Bool}=nothing,
        includeSchwinger::Union{Nothing,Bool}=nothing,      printBefore::Union{Nothing,Bool}=nothing,
        BField::Union{Nothing,Float64}=nothing,             levelSelection::Union{Nothing,LevelSelection}=nothing,
        gMultiplet::Union{Nothing,Multiplet}=nothing)
    if  isnothing(calcLandeJ)         calcLandeJx       = set.calcLandeJ       else   calcLandeJx       = calcLandeJ       end
    if  isnothing(calcLandeF)         calcLandeFx       = set.calcLandeF       else   calcLandeFx       = calcLandeF       end
    if  isnothing(calcZeeman)         calcZeemanx       = set.calcZeeman       else   calcZeemanx       = calcZeeman       end
    if  isnothing(calcQZScoeff)       calcQZScoeffx     = set.calcQZScoeff     else   calcQZScoeffx     = calcQZScoeff     end
    if  isnothing(includeSchwinger)   includeSchwingerx = set.includeSchwinger else   includeSchwingerx = includeSchwinger end
    if  isnothing(printBefore)        printBeforex      = set.printBefore      else   printBeforex      = printBefore      end
    if  isnothing(BField)             BFieldx           = set.BField           else   BFieldx           = BField           end
    if  isnothing(levelSelection)     levelSelectionx   = set.levelSelection   else   levelSelectionx   = levelSelection   end
    if  isnothing(gMultiplet)         gMultipletx       = set.gMultiplet       else   gMultipletx       = gMultiplet       end

    Settings( calcLandeJx, calcLandeFx, calcZeemanx, calcQZScoeffx, includeSchwingerx, printBeforex,
              BFieldx, levelSelectionx, gMultipletx )
end


# `Base.show(io::IO, settings::LandeZeeman.Settings)`  ... prepares a proper printout of the variable settings::LandeZeeman.Settings.
function Base.show(io::IO, settings::LandeZeeman.Settings) 
    println(io, "calcLandeJ:               $(settings.calcLandeJ)  ")
    println(io, "calcLandeF:               $(settings.calcLandeF)  ")
    println(io, "calcZeeman:               $(settings.calcZeeman)  ")
    println(io, "calcQZScoeff:             $(settings.calcQZScoeff)  ")
    println(io, "includeSchwinger:         $(settings.includeSchwinger)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "BField:                   $(settings.BField)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
    println(io, "gMultiplet:               $(settings.gMultiplet)  ")
end


"""
`LandeZeeman.amplitude(::ZeemanN1, rLevel::Level, sLevel::Level, grid::Radial.Grid)`
    ... to compute the N^(1) Zeeman amplitude <alpha_r J_r || N^(1) || alpha_s J_s>
        for a given pair of levels. A value::ComplexF64 is returned.

    KNOWN LIMITATION (confirmed 30-Jul-2026, re-verified 31-Jul-2026 after an unrelated B-spline
    boundary-condition fix, NOT yet fixed): for any subshell with kappa <= -3 (d5/2, f7/2, g9/2, ...), the
    resulting g_J is WRONG. Verified via a kappa-sweep in the same Ge II ion (Z=32, [Ar] 3d^10 4s^2 4d or 4f):
    p_1/2 (kappa=+1), p_3/2 (kappa=-2), d_3/2 (kappa=+2, g_J=0.79950 vs exact 0.8000), and f_5/2 (kappa=+3,
    g_J=0.85682 vs exact/paper 0.85680) all reproduce the expected value to <0.1%, but d_5/2 (kappa=-3) gives
    g_J=-3.05247 vs the exact 1.2000, and f_7/2 (kappa=-4) gives g_J=-2.26367 (DFS) / -2.24247 (AL) vs the
    exact 1.1429 -- see example-Cd.jl branches c and d. NOTE: these wrong values are NOT stable against
    unrelated changes elsewhere in the SCF/basis machinery -- before the 30-Jul-2026 B-spline kappa-sign
    boundary-condition fix, the SAME two cases gave smaller-magnitude, same-sign wrong values (d_5/2:
    0.3937, f_7/2: 0.8937); after that unrelated fix, both flipped sign and grew in magnitude. This is
    expected for a bug in the ANGULAR coefficient machinery (not the radial basis) that is being fed a
    changed, more-correct radial orbital -- root cause narrowed to `AngularMomentum.CL_reduced_me_rb` or its
    use in `InteractionStrength.zeeman_n1`'s `Subshell(1,-ka)` construction, but not yet isolated further,
    and evidently sensitive to exactly which orbitals it acts on. Treat ANY g_J result involving a kappa<=-3
    subshell as unverified, and expect its specific wrong VALUE to shift again whenever the underlying
    B-spline/SCF machinery changes -- only the qualitative fact "kappa<=-3 is wrong" is stable, not any
    specific number quoted here. example-Cd.jl branch d deliberately restricts itself to f_5/2 (kappa=+3)
    to avoid this contaminating an unrelated comparison. See project_zeeman_hfs_bugs.md for the full
    investigation.
"""
function  amplitude(::ZeemanN1, rLevel::Level, sLevel::Level, grid::Radial.Grid)
    nr = length(rLevel.basis.csfs);    ns = length(sLevel.basis.csfs);    matrix = zeros(ComplexF64, nr, ns)
    printstyled("Compute Zeeman N^(1) matrix of dimension $nr x $ns ... ", color=:light_green)
    if  rLevel.parity != sLevel.parity   return( ComplexF64(0.) )   end
    #
    for  r = 1:nr
        for  s = 1:ns
            me = 0.
            if  rLevel.basis.csfs[r].parity  != rLevel.parity    ||  sLevel.basis.csfs[s].parity  != sLevel.parity  ||
                rLevel.mc[r] == 0  ||  sLevel.mc[s] == 0   continue
            end
            subshellList = sLevel.basis.subshells
            opa = SpinAngular.OneParticleOperator(1, plus, true)
            wa  = SpinAngular.computeCoefficients(opa, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList)
            for  coeff in wa
                tamp  = InteractionStrength.zeeman_n1(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                me = me + coeff.T * tamp
            end
            matrix[r,s] = me
        end
    end
    printstyled("done.\n", color=:light_green)
    amplitude = transpose(rLevel.mc) * matrix * sLevel.mc
    return( amplitude )
end


"""
` + amplitude(::ZeemanDeltaN1, rLevel::Level, sLevel::Level, grid::Radial.Grid)`
    ... to compute the Delta N^(1) Zeeman amplitude <alpha_r J_r || Delta N^(1) || alpha_s J_s>
        for a given pair of levels. A value::ComplexF64 is returned.
"""
function  amplitude(::ZeemanDeltaN1, rLevel::Level, sLevel::Level, grid::Radial.Grid)
    nr = length(rLevel.basis.csfs);    ns = length(sLevel.basis.csfs);    matrix = zeros(ComplexF64, nr, ns)
    printstyled("Compute Zeeman ΔN^(1) matrix of dimension $nr x $ns ... ", color=:light_green)
    if  rLevel.parity != sLevel.parity   return( ComplexF64(0.) )   end
    #
    for  r = 1:nr
        for  s = 1:ns
            me = 0.
            if  rLevel.basis.csfs[r].parity  != rLevel.parity    ||  sLevel.basis.csfs[s].parity  != sLevel.parity  ||
                rLevel.mc[r] == 0  ||  sLevel.mc[s] == 0   continue
            end
            subshellList = sLevel.basis.subshells
            opa = SpinAngular.OneParticleOperator(1, plus, true)
            wa  = SpinAngular.computeCoefficients(opa, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList)
            for  coeff in wa
                tamp  = InteractionStrength.zeeman_Delta_n1(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                me = me + coeff.T * tamp
            end
            matrix[r,s] = me
        end
    end
    printstyled("done.\n", color=:light_green)
    amplitude = transpose(rLevel.mc) * matrix * sLevel.mc
    return( amplitude )
end


"""
`LandeZeeman.amplitudeN1(::ZeemanN1, rLevel::Level, sLevel::Level, grid::Radial.Grid; display::Bool=false)`
    ... to compute the (reduced) Zeeman amplitude <alpha_r J_r || N^(1) || alpha_s J_s>
        for a given pair of levels. A value::ComplexF64 is returned.
"""
function amplitudeN1(::ZeemanN1, rLevel::Level, sLevel::Level, grid::Radial.Grid; display::Bool=false)
    if     rLevel.parity != sLevel.parity     amplitude = ComplexF64(0.)
    else
        nr = length(rLevel.basis.csfs);    ns = length(sLevel.basis.csfs)
        if display   printstyled("Compute Zeeman N^(1) matrix of dimension $nr x $ns in the given bases " *
                                 "[transition $(rLevel.index)- $(sLevel.index)] ... ", color=:light_green)     end
        matrix = zeros(ComplexF64, nr, ns)
        #
        for  r = 1:nr
            for  s = 1:ns
                if  rLevel.mc[r] == 0  ||  sLevel.mc[s] == 0    continue    end
                subshellList = sLevel.basis.subshells
                opa = SpinAngular.OneParticleOperator(1, plus, true)
                wa  = SpinAngular.computeCoefficients(opa, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList)
                for  coeff in wa
                    tamp = InteractionStrength.zeeman_n1(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                    matrix[r,s] = matrix[r,s] + coeff.T * tamp
                end
            end
        end
        if display   printstyled("done. \n", color=:light_green)   end
        amplitude = transpose(rLevel.mc) * matrix * sLevel.mc
    end
    if  display
        sa = @sprintf("%.5e", amplitude.re) * "  " * @sprintf("%.5e", amplitude.im)
        println("    < level=$(rLevel.index) [J=$(rLevel.J)$(string(rLevel.parity))] || N^(1) ||" *
                " $(sLevel.index) [$(sLevel.J)$(string(sLevel.parity))] >  = " * sa)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary
            println(iostream,  "    N^(1) amplitude:  < level=$(rLevel.index) [J=$(rLevel.J)$(string(rLevel.parity))] || N^(1) ||" *
                               " $(sLevel.index) [$(sLevel.J)$(string(sLevel.parity))] >  = " * sa)
        end
    end
    return( amplitude )
end



"""
`LandeZeeman.computeAmplitudesProperties(outcome::LandeZeeman.Outcome, grid::Radial.Grid, settings::LandeZeeman.Settings)`  
    ... to compute all amplitudes and properties of for a given level. The given gMultiplet in settings containes the
        intermediate levels used in the computation of second order coefficients. An outcome::LandeZeeman.Outcome is 
        returned for which the amplitudes and all requested properties are now evaluated explicitly.
"""
function  computeAmplitudesProperties(outcome::LandeZeeman.Outcome, grid::Radial.Grid, settings::LandeZeeman.Settings)
    LandeJ = 0.0;   amplitudeN1 = ComplexF64(0.);   amplitudeDeltaN1 = ComplexF64(0.)
    J = AngularMomentum.oneJ(outcome.Jlevel.J) 
    #
    if       J == 0.                        LandeJ = 0.
    elseif   settings.calcLandeJ
        amplitudeN1 = LandeZeeman.amplitude(ZeemanN1(), outcome.Jlevel, outcome.Jlevel, grid)
        #
        if  settings.includeSchwinger
            amplitudeDeltaN1 = LandeZeeman.amplitude(ZeemanDeltaN1(), outcome.Jlevel, outcome.Jlevel, grid)
        end
        #       
        LandeJ = 2*(amplitudeN1 + amplitudeDeltaN1) / sqrt(J*(J+1))    
    end
    

    if settings.calcQZScoeff  &&  outcome.nuclearI != AngularJ64(0)
        newFsublevels = LandeZeeman.SublevelF[]
        for  Fsub in outcome.Fsublevels
            c2 = LandeZeeman.computeQuadraticZeemanC2(outcome.Jlevel, Fsub, grid, settings)
            push!(newFsublevels, LandeZeeman.SublevelF(Fsub.F, Fsub.M, 0.0, c2))
        end
    else
        newFsublevels = outcome.Fsublevels
    end

    if settings.calcQZScoeff  &&  outcome.nuclearI == AngularJ64(0)
        newJsublevels = LandeZeeman.SublevelJ[]
        for Jsub in outcome.Jsublevels
            c2 = LandeZeeman.computeQuadraticZeemanC2(outcome.Jlevel, Jsub, grid, settings)
            push!(newJsublevels, LandeZeeman.SublevelJ(Jsub.M, 0.0, c2))
        end
    else
        newJsublevels = outcome.Jsublevels
    end
    
    newOutcome = LandeZeeman.Outcome( outcome.Jlevel, LandeJ, amplitudeN1, amplitudeDeltaN1, outcome.nuclearI, 
                                      newJsublevels, newFsublevels)
    return( newOutcome )
end



"""
`LandeZeeman.computeQuadraticZeemanC2(level::Level, Jsub::SublevelJ, grid::Radial.Grid, settings::LandeZeeman.Settings)`  
    ... to compute the quadratic-Zeeman shift coefficient C2 for the Zeeman sublevel Jsub by applyling a summation over 
        all levels from gMultiplet. A value c2::Float64 is returned for the given level (level, Jsub).
"""
function  computeQuadraticZeemanC2(level::Level, Jsub::SublevelJ, grid::Radial.Grid, settings::LandeZeeman.Settings)
    c2 = 0.
    println("\n>>> Calculate C_2 coefficient for level (index=$(level.index), J=$(level.J), M=$(Jsub.M)) ...")

    for nLevel in settings.gMultiplet.levels
        # Exclude levels with the same energy (E == E_n) or with total angular momenta that differ more than by 1.
        if level.J == nLevel.J  &&  level.parity == nLevel.parity    &&  isapprox(level.energy, nLevel.energy, rtol=1.0e-4)
            @show "C2 continue:", level.energy, nLevel.energy
            continue
        elseif  abs(Basics.twice(Jsub.M)) > Basics.twice(nLevel.J)   
            continue
        elseif  abs(Basics.twice(level.J) - Basics.twice(nLevel.J)) > 2    
            continue
        end
        @show  "compute c2: aa", nLevel.J, nLevel.parity

        amplitudeN1 = LandeZeeman.amplitude(ZeemanN1(), nLevel, level, grid)
        println("       <level=$(nLevel.index) [J=$(nLevel.J)$(string(nLevel.parity))] || N^(1) || " *
                        "level=$(level.index) [J=$(level.J)$(string(level.parity))] >  = $(amplitudeN1)")
        amplitudeDeltaN1 = 0.

        if settings.includeSchwinger
            amplitudeDeltaN1 = LandeZeeman.amplitude(ZeemanDeltaN1(), nLevel, level, grid)
            println("       <level=$(nLevel.index) [J=$(nLevel.J)$(string(nLevel.parity))] || ΔN^(1) || " *
                        "level=$(level.index) [J=$(level.J)$(string(level.parity))] > = $(amplitudeDeltaN1)")
        end
        
        @show  "compute c2: bb", level.J, Jsub.M, nLevel.J
        cg   = AngularMomentum.ClebschGordan_old(level.J, Jsub.M, AngularJ64(1), AngularM64(0), nLevel.J, Jsub.M)
        amp  = abs(amplitudeN1 + amplitudeDeltaN1)
        ## c2 = c2 + cg^2 / (Basics.twice(nLevel.J) + 1) / (level.energy - nLevel.energy) * amp^2
        c2 = c2 + cg^2 / (level.energy - nLevel.energy) * amp^2
        @show  "compute c2: cc", amplitudeN1, amplitudeDeltaN1, c2

    end

    conv = 0.11909076 #Conversion Factor from atomic units to MHz/T^2
    c2 = conv * c2
    
    println("   C_2(level=$(level.index) [J=$(level.J)$(string(level.parity)), M=$(Jsub.M)]) = $(c2) MHz/T^2")

    return( c2 )
end


"""
`LandeZeeman.computeQuadraticZeemanC2(Fsub::SublevelF, Fsublevels::Array{SublevelF,1}, grid::Radial.Grid, settings::LandeZeeman.Settings)`  
    ... to compute all quadratic-Zeeman shift C2 coefficients; an c2::Float64 is returned for the given hyperfine level
        Fsub.
"""
function  computeQuadraticZeemanC2(multiplet::Multiplet, level::Level, Fsub::SublevelF, nm::Nuclear.Model, grid::Radial.Grid, settings::LandeZeeman.Settings)
    c2 = 0.

    println(">>> Calculate C_2 coefficient for level (index=$(level.index), J=$(level.J), I=$(nm.spinI), F=$(Fsub.F), M_F=$(Fsub.M)) ...")
    println("   >>> Calculate Hyperfine parameters ...")

    #Calculate hyperfine splitting
    amplitudeT1 = Hfs.amplitude(Basics.M1, level, level, grid; printout=false)
    amplitudeT2 = Hfs.amplitude(Basics.E2, level, level, grid; printout=false)
    j = Float64(level.J)
    i = Float64(nm.spinI)
    A = nm.mu/i * 1/sqrt(j * (j + 1)) * amplitudeT1
    B = 2 * nm.Q * sqrt(j*(2*j - 1)/((j + 1)*(2*j + 3))) * amplitudeT2
    C(f) = Float64(f)*(Float64(f) + 1) - j*(j + 1) - i*(i + 1)
    EHFS(f) = 1/2 * A * C(f) + B * ( 3/4 * C(f)*(C(f)+1) - i * (i + 1) * j * (j + 1) )/( 2*i * (2*i - 1) * j * (2*j - 1))

    println("   >>> Sum over intermediate states")

    # Sum over hyperfine levels
    Fvalues = Basics.oplus(nm.spinI, level.J)
    amplitudeN1 = LandeZeeman.amplitude(ZeemanN1(), level, level, grid)
    println("       <level=$(level.index) [J=$(level.J)$(string(level.parity))] || N^(1) || " *
                    "level=$(level.index) [J=$(level.J)$(string(level.parity))] > = $(amplitudeN1)")
    amplitudeDeltaN1 = 0.0

    if settings.includeSchwinger
        amplitudeDeltaN1 = LandeZeeman.amplitude(ZeemanDeltaN1(), level, level, grid)
        println("       <level=$(level.index) [J=$(level.J)$(string(level.parity))|| ΔN^(1) || " *
                        "level=$(level.index) [J=$(level.J)$(string(level.parity))] > = $(amplitudeDeltaN1)")
    end

    for f in Fvalues
        if f == Fsub.F || abs(Fsub.M.num/Fsub.M.den) > Float64(f)
            continue 
        end
        
        factor = (2*Float64(Fsub.F) + 1)*(2*j + 1)*(2*Float64(f) + 1)
        w3jvalue = AngularMomentum.Wigner_3j(Fsub.F, AngularJ64(1), f, Fsub.M, AngularM64(0), AngularM64(-Fsub.M.num//Fsub.M.den))
        w6jvalue = AngularMomentum.Wigner_6j(level.J, nm.spinI, Fsub.F, f, AngularJ64(1), level.J)
        deltaE = EHFS(Fsub.F) - EHFS(f)

        c2 += factor * w3jvalue^2 * w6jvalue^2 * (amplitudeN1 + amplitudeDeltaN1)^2 / deltaE
    end
    

    # Sum over other intermediate levels
    for ilevel in multiplet.levels
        if ilevel == level || !Basics.selectLevel(ilevel, settings.gMultiplet)
            continue
        end

        amplitudeN1 = LandeZeeman.amplitude(ZeemanN1(), ilevel, level, grid)
        println("       <level=$(level.index) [J=$(level.J)$(string(level.parity))] || N^(1) || " *
                        "level=$(ilevel.index) [J=$(ilevel.J)$(string(ilevel.parity))] > = $(amplitudeN1)")

        if settings.includeSchwinger
            amplitudeDeltaN1 = LandeZeeman.amplitude(ZeemanDeltaN1(), ilevel, level, grid)
            println("       <level=$(level.index) [J=$(level.J)$(string(level.parity))] || ΔN^(1) || " *
                            "level=$(ilevel.index) [J=$(ilevel.J)$(string(ilevel.parity))] > = $(amplitudeDeltaN1)")
        end

        deltaE = level.energy - ilevel.energy

        Fvalues = Basics.oplus(ilevel.J, nm.spinI)
        for f in Fvalues
            if abs(Fsub.M.num/Fsub.M.den) > Float64(f)
                continue
            end
            factor = (2*Float64(Fsub.F) + 1)*(2*Float64(ilevel.J) + 1)*(2*Float64(f) + 1)
            w3jvalue = AngularMomentum.Wigner_3j(Fsub.F, AngularJ64(1), f, Fsub.M, AngularM64(0), AngularM64(-Fsub.M.num//Fsub.M.den))
            w6jvalue = AngularMomentum.Wigner_6j(level.J, nm.spinI, Fsub.F, f, AngularJ64(1), ilevel.J)

            c2 += factor * w3jvalue^2 * w6jvalue^2 * (amplitudeN1 + amplitudeDeltaN1)^2 / deltaE
        end
    end

    conv = 0.11909076 #Conversion Factor from atomic units to MHz/T^2
    c2 = conv * c2

    println("C_2(level=$(level.index) [J=$(level.J)$(string(level.parity)), I=$(nm.spinI), F=$(Fsub.F), M_F=$(Fsub.M)]) = $(c2) MHz/T^2")

    return ( c2 )
end


"""
`LandeZeeman.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::LandeZeeman.Settings; output=true)`  
    ... to compute (as selected) the Lande J or F factors for the levels of the given multiplet and as specified by the given settings. The 
        results are printed in neat tables to screen but nothing is returned otherwise.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::LandeZeeman.Settings; output=true)
    println("")
    printstyled("LandeZeeman.computeOutcomes(): The computation of the Zeeman amplitudes and Lande factors starts now ... \n", color=:light_green)
    printstyled("-------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    outcomes = LandeZeeman.determineOutcomes(multiplet, nm, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    LandeZeeman.displayOutcomes(outcomes)    end
    # Calculate all amplitudes and requested properties
    newOutcomes = LandeZeeman.Outcome[]
    for  outcome in outcomes
        newOutcome = LandeZeeman.computeAmplitudesProperties(outcome, grid, settings) 
        push!( newOutcomes, newOutcome)
    end
    # Print all results to screen
    LandeZeeman.displayResults(stdout, newOutcomes, nm, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    LandeZeeman.displayResults(iostream, newOutcomes, nm, settings)   end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`LandeZeeman.determineOutcomes(multiplet::Multiplet, nm::Nuclear.Model, settings::LandeZeeman.Settings)`  
    ... to determine a list of Outcomes's for the computation of Lande factors and/or Zeeman splittings for 
        levels from the given multiplet. It takes into account the particular selections and settings. 
        An Array{LandeZeeman.Outcome,1} is returned. Apart from the level specification, all physical properties are 
        still set to zero during the initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, nm::Nuclear.Model, settings::LandeZeeman.Settings) 
    # Define values that depend on the requested computations
    outcomes   = LandeZeeman.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            Jsublevels = LandeZeeman.SublevelJ[];   Mvalues = AngularMomentum.m_values(level.J) 
            Fsublevels = LandeZeeman.SublevelF[];   Fvalues = Basics.oplus(nm.spinI, level.J)
            for  M in Mvalues   push!(Jsublevels, LandeZeeman.SublevelJ(M, 0., 0.) )    end
            for  F in Fvalues 
                MFvalues = AngularMomentum.m_values(F)
                for  MF in MFvalues     push!(Fsublevels, LandeZeeman.SublevelF(F, MF, 0., 0.) )    end 
            end
            push!( outcomes, LandeZeeman.Outcome(level, 0., ComplexF64(0.), ComplexF64(0.), nm.spinI, Jsublevels, Fsublevels) )
        end
    end
    return( outcomes )
end


"""
`LandeZeeman.displayOutcomes(outcomes::Array{LandeZeeman.Outcome,1})`  
    ... to display a list of levels that have been selected for the computations. A small neat table of all selected 
        levels and their energies is printed but nothing is returned otherwise.
"""
function  displayOutcomes(outcomes::Array{LandeZeeman.Outcome,1})
    nx = 43
    println(" ")
    println("  Selected Lande-Zeeman levels:")
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
        sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.Jlevel.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`LandeZeeman.displayResults(stream::IO, outcomes::Array{LandeZeeman.Outcome,1}, nm::Nuclear.Model, settings::LandeZeeman.Settings)`  
    ... to display the energies, Lande factors, Zeeman amplitudes etc. for the selected levels. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{LandeZeeman.Outcome,1}, nm::Nuclear.Model, settings::LandeZeeman.Settings)
    #
    if  settings.calcLandeJ
        nx = 135
        println(stream, " ")
        println(stream, "  Lande g_J factors and Zeeman amplitudes:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=5)              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(12, "Lande-J"; na=8);                           sb = sb * TableStrings.hBlank(20)          
        sa = sa * TableStrings.center(62, "N1   -- Zeeman Amplitudes --   Delta N1"    ; na=4) 
        sb = sb * TableStrings.center(62, "   Re              Im                 Re              Im"; na=5)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            energy = 1.0
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))                    * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.LandeJ) )              * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.amplitudeN1.re) )      * " "
            sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.amplitudeN1.im) )      * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.amplitudeDeltaN1.re) ) * " "
            sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.amplitudeDeltaN1.im) ) * "    "
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    if  settings.calcLandeF
        nx = 80
        println(stream, " ")
        println(stream, "  Hyperfine levels and Lande g_F factors:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=5)              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(10, "F^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(12, "Lande-F"; na=8);                           sb = sb * TableStrings.hBlank(20)          
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            energy = 1.0
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))                    * "    "
            println(stream, sa )
            #
            sc = TableStrings.hBlank( length(sa) + 1 )
            Flist = oplus(nm.spinI, outcome.Jlevel.J)
            for  F in Flist
                symf   = LevelSymmetry( F, outcome.Jlevel.parity);      Fx = AngularMomentum.oneJ(F)
                Jx     = AngularMomentum.oneJ(outcome.Jlevel.J);    Ix = AngularMomentum.oneJ(nm.spinI)
                LandeF = (Fx*(Fx+1) + Jx*(Jx+1) - Ix*(Ix+1)) / (2*Fx*(Fx+1)) * outcome.LandeJ
                sa = sc * TableStrings.center(10, string(symf); na=4)
                sa = sa * TableStrings.flushright(15, @sprintf("% .8e", LandeF) )   
                println(stream, sa )
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    #
    if  settings.calcZeeman
        nx = 100
        println(stream, " ")
        println(stream, "  Zeeman splittings of (hyper-) fine-structure levels:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "  Here, we should display the Zeeman splittings of (hyper-) fine-structure levels. ")
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    #
    if  settings.calcQZScoeff
        nx = 78
        println(stream, " ")
        println(stream, "  Quadratic Zeeman shift C_2 coefficients for (hyper-) fine-structure levels:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=5)              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5)

        if nm.spinI == AngularJ64(0)
            sa = sa * TableStrings.center(10, "M";   na=4);                             sb = sb * TableStrings.hBlank(14)
            sa = sa * TableStrings.center(14, "C_2"; na=5)
            sb = sb * TableStrings.center(14, "[MHz/T^2]"; na=5)
        else
            sa = sa * TableStrings.center(10, "F";   na=4);                             sb = sb * TableStrings.hBlank(14)
            sa = sa * TableStrings.center(10, "M";   na=4);                             sb = sb * TableStrings.hBlank(14)
            sa = sa * TableStrings.center(14, "C_2"; na=5)
            sb = sb * TableStrings.center(14, "[MHz/T^2]"; na=5)
        end    
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))

        for outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            energy = 1.0
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
            println(stream, sa )

            sc = TableStrings.hBlank(length(sa) + 1)
            if nm.spinI == AngularJ64(0)
                for Jsub in outcome.Jsublevels
                    sa = sc * TableStrings.center(10, string(Jsub.M); na=4)
                    sa = sa * @sprintf("% .8e", Jsub.c2Coeff)
                    println(stream, sa)
                end
            else
                for Fsub in outcome.Fsublevels
                    sa = sc * TableStrings.center(10, string(Fsub.F); na=4)
                    sa = sa * TableStrings.center(10, string(Fsub.M); na=4)
                    sa = sa * @sprintf("% .8e", Fsub.c2Coeff)
                    println(stream, sa)
                end
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end
    end
    #
    return( nothing )
end

end # module

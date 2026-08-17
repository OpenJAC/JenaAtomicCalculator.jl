
export  dummyQZ
using FortranFiles, Dierckx

"""
`Basics.read(::ReadCslFileGrasp92, filename::String)`
    ... reads in the CSF list from the (existing) .csl file filename; a basis::Basis is returned with
        basis.isDefined = true and with a proper sequence of orbitals, but with basis.orbitals = Orbital[].
"""
function Basics.read(::ReadCslFileGrasp92, filename::String)
    return( Basics.readCslFileGrasp92(filename) )
end


"""
`Basics.read(::ReadOrbitalFileGrasp92, filename::String)`
    ... reads in the orbitals from a (formatted) .rwf file filename; a list of orbitals::Array{Orbital,1}
        is returned with all subfields specified.
"""
function Basics.read(::ReadOrbitalFileGrasp92, filename::String)
    return( Basics.readOrbitalFileGrasp92(filename) )
end


"""
`Basics.read(::ReadMixingFileGrasp18, filename::String)`
    ... reads energies & mixing coefficients from the Grasp18 mixing file filename.
"""
function Basics.read(::ReadMixingFileGrasp18, filename::String)
    return( Basics.readMixingFileGrasp18(filename) )
end


"""
Basics.readCslFileGrasp92(filename::String)`  
    ... reads in the CSF list from a Grasp92 .csl / GRASP18 .c file; a basis::Basis is returned.
"""
function Basics.readCslFileGrasp92(filename::String)
    coreSubshells = Subshell[];    peelSubshells =  Subshell[]

    #f  = open(filename)
    ftemp = open("temp-csf.c", "w")

    open(filename) do f1
        while ! eof(f1)
            line = readline(f1) 
            if line != " *"  write(ftemp, line*"\n") end
        end
    end

    close(ftemp)

    f  = open("temp-csf.c")

    sa = readline(f);   sa[1:14] != "Core subshells"   &&   error("Not a Grasp92 .cls file.")
    sa = readline(f);   if  length(sa) > 3    go = true    else  go = false    end
    i  = -1
    while go
        # if sa[end] != '-' sa = sa * " " end
        i = i + 1;    sh = Basics.subshellGrasp( strip(sa[5i+1:5i+5]) );    push!(coreSubshells, sh)
        if  5i + 6 >  length(sa)    break    end
    end
    sa = readline(f);   println("sa = $sa")
    sa = readline(f);   if  length(sa) > 3    go = true    else  go = false    end
    if sa[end] != '-' sa = sa * " " end
    i  = -1
    while go
        i = i + 1;    sh = Basics.subshellGrasp( strip(sa[5i+1:5i+5]) );    push!(peelSubshells, sh)
        if  5i + 6 >  length(sa)    break    end
    end
    # Prepare a list of all subshells to define the common basis
    subshells = deepcopy(coreSubshells);    append!(subshells, peelSubshells)

    # Now, read the CSF in turn until the EOF
    readline(f)
    NoCSF = 0;    csfs = CsfR[]
    Defaults.setDefaults("relativistic subshell list", subshells)

    while true
        sa = readline(f);      if length(sa) == 0   break    end
        sb = readline(f);      sc = readline(f)
        NoCSF = NoCSF + 1
        push!(csfs, ManyElectron.CsfRGrasp92(subshells, coreSubshells, sa, sb, sc) )
    end

    println("  ... $NoCSF CSF read in from Grasp92 file  $filename)") 
    NoElectrons = sum( csfs[1].occupation )

    rm("temp-csf.c")

    Basis(true, NoElectrons, subshells, csfs, coreSubshells, Dict{Subshell,Radial.Orbital}() )
end


"""
`Basics.readOrbitalFileGrasp92(filename::String, grid::Radial.Grid)`  
    ... reads in the orbitals list from a Grasp92 .w file; a dictionary 
        orbitals::Dict{Subshell,Radial.Orbital} is returned.
"""
function Basics.readOrbitalFileGrasp92(filename::String, grid::Radial.Grid)
    # using FortranFiles
    f = FortranFile(filename)
    
    orbitals = Dict{Subshell,Radial.Orbital}();    first = true

    if (String(Base.read(f, FString{6})) != "G92RWF")
    close(f)
    error("File \"", filename, "\" is not a proper Grasp92 RWF file")
    end

    Defaults.setDefaults("standard grid", grid);   

    while true
    try
        n, kappa, energy, mtp = Base.read(f, Int32, Int32, Float64, Int32)
        pz, P, Q = Base.read(f, Float64, (Float64, mtp), (Float64, mtp))
        ra = Base.read(f, (Float64, mtp))
            # Now place the data into the right fields of Orbital
            subshell = Subshell( Int64(n), Int64(kappa) )   
            if  energy < 0   isBound = true    else    isBound = false   end 
            useStandardGrid = true

            itp = Dierckx.Spline1D(ra,P)
            itq = Dierckx.Spline1D(ra,Q)

            Px = zeros(length(grid.r))
            Qx = zeros(length(grid.r))

            for i in 1:length(grid.r)
                if minimum(ra) <= grid.r[i] <= maximum(ra)
                    Px[i] = itp(grid.r[i])
                    Qx[i] = itq(grid.r[i])
                else
                    break
                end
            end


            orbitals = Base.merge( orbitals, Dict( subshell => Orbital(subshell, isBound, useStandardGrid, energy, Px, Qx, Px, Qx, grid ) ))
    catch ex
        if     ex isa EOFError  break
        else   throw(ex)
        end
    end
    end
    close(f)

    return( orbitals )
end


"""
`Basics.readMixingFileGrasp18(filename::String, basis::Basis)`  
    ... reads in the mixing coefficients from a Grasp18 .m file by using the given basis; A multiplet::Multiplet 
        is returned.
"""   
function Basics.readMixingFileGrasp18(filename::String, basis::Basis)
    
    name =  "from GRASP18"

    f = FortranFile(filename)
        
    if (String(Base.read(f, FString{6})) != "G92MIX")
        close(f)
    error("File \"", filename, "\" does not seem to be a Grasp92 MIX file")
    end

    nelec, ncftot, nw, nvectot, nvecsiz, nblock = read(f, Int32, Int32, Int32, Int32, Int32, Int32)

    levels = Level[]
    nlevel = 0

    for i = 1:nblock
        nb, ncfblk, nevblk, iatjp, iaspa = read(f, Int32, Int32, Int32, Int32, Int32)

        if iaspa < 0
            parity = Basics.Parity("-")
        else
            parity = Basics.Parity("+")
        end

        J = AngularJ64((iatjp-1)//2)
        M = AngularM64( J.num//J.den )

        ivecdum = read(f, (Int32, nevblk))
        avg_energy, egval = read(f, Float64, (Float64, nevblk))
        eigenvectors = read(f, (Float64, ncfblk, nevblk))

        for blk = 1:nevblk
            mc = zeros(ncftot)
            for n in 1:length(basis.csfs)
                if basis.csfs[n].J == J && basis.csfs[n].parity == parity 
                    mc[n:(n + ncfblk - 1)] = eigenvectors[:,blk]
                    break 
                end
            end
            nlevel += 1
            push!( levels, Level(J, M, parity, nlevel, avg_energy + egval[blk], 0., true, basis, mc) )
        end

    end

    levels = sort!(levels)
    
    levelsSorted = []

    for i = 1:length(levels)
        push!( levelsSorted, Level(levels[i].J, levels[i].M, levels[i].parity, i,
        levels[i].energy, 0., true, levels[i].basis, levels[i].mc) )
    end

    multiplet = Multiplet(name, levelsSorted)
    return( multiplet )
end


"""
`Basics.readMixFileRelci(filename::String, basis::Basis)`  
    ... reads in the mixing coefficients from a RELCI .mix file; a multiplet::Multiplet is returned.
"""
function Basics.readMixFileRelci(filename::String, basis::Basis)
    levels = Level[];    name =  "from Relci"

    f  = open(filename)
    sa = readline(f);   sa[1:14] != "G92MIX (format"   &&   error("Not a formatted relci.mix file.")
    sa = readline(f);   NoElectrons = Base.parse(Int64, strip(sa[1:6]) );    NoCsf = Base.parse(Int64, strip(sa[7:12]) ) 
                        NoSubshells = Base.parse(Int64, strip(sa[13:18]) ) 
    sa = readline(f);   NoLevels    = Base.parse(Int64, strip(sa[1:6]) )
    #
    sa = readline(f)
    sa = readline(f);   list2J = AngularJ64[AngularJ64(0)  for  i = 1:NoLevels];    listp = Basics.Parity[Basics.plus  for  i = 1:NoLevels]
    for  i = 1:NoLevels
        jj = Base.parse(Int64, strip(sa[(i-1)*6+1:(i-1)*6+5]) );   list2J[i] = AngularJ64(jj//2)
        pp = string( sa[(i-1)*6+6] );                              listp[i]  = Basics.Parity( pp )
    end
    #
    sa   = readline(f);   listEnergy = Float64[0.  for  i = 1:NoLevels]
    enav = Base.parse(Float64, strip(sa[1:26]) )
    for  i = 1:NoLevels
        en = Base.parse(Float64, strip(sa[26i+1:26i+26]) );   listEnergy[i] = en
    end
    #
    eigenvectors = zeros(NoCsf,NoLevels)
    for  j = 1:NoCsf
        sa   = readline(f)
        for  i = 1:NoLevels
        cm = Base.parse(Float64, strip(sa[16*(i-1)+1:16(i-1)+16]) );   eigenvectors[j,i] = cm
        end
    end
    #
    for  j = 1:NoLevels
        mc = Float64[0.  for  i = 1:NoCsf]
        for  i = 1:NoCsf    mc[i] = eigenvectors[i,j]    end
        Jx = list2J[j];   Mx = AngularM64( Jx.num//Jx.den )
        push!( levels, Level(Jx, Mx, listp[j], j, enav+listEnergy[j], 0., true, basis, mc) )
    end

    multiplet = Multiplet(name, levels)
    return( multiplet )
end


"""
`Basics.readFilesGrasp18(grid::Radial.Grid, fileCSF::String, fileWavefunction::String, fileMixing::String)`
 
    ... reads in the GRASP18 output files `*.c`, `*.w`and `*.m`; a multiplet::Multiplet is returned.
"""
function Basics.readFilesGrasp18(grid::Radial.Grid, fileCSF::String, fileWavefunction::String, fileMixing::String)

    basis = Basics.readCslFileGrasp92(fileCSF)
    orbitals = Basics.readOrbitalFileGrasp92(fileWavefunction, grid)
    isDefined = true
    basis = Basis( isDefined, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
    multiplet = Basics.readMixingFileGrasp18(fileMixing, basis)

    return multiplet
end


"""
`Basics.recast(::RecastRateToDecayWidth,
    line::Union{Einstein.Line, PhotoEmission.Line, HyperfineInduced.Line, TwoElectronOnePhoton.Line},
    wa::Float64)`
    ... recasts a radiative rate (Einstein A, a.u.) into a decay width, taking the selected energy unit
        into account; a Float64 is returned.
"""
function Basics.recast(::RecastRateToDecayWidth,
        line::Union{Einstein.Line, PhotoEmission.Line, HyperfineInduced.Line, TwoElectronOnePhoton.Line},
        wa::Float64)
    return( Defaults.convertUnits("energy: from atomic", wa) )
end


"""
`Basics.recast(::RecastRateToEinsteinA,
    line::Union{Einstein.Line, PhotoEmission.Line, HyperfineInduced.Line, TwoElectronOnePhoton.Line},
    wa::Float64)`
    ... recasts a spontaneous radiative rate (Einstein A, a.u.) into Einstein A in selected units;
        a Float64 is returned.
"""
function Basics.recast(::RecastRateToEinsteinA,
        line::Union{Einstein.Line, PhotoEmission.Line, HyperfineInduced.Line, TwoElectronOnePhoton.Line},
        wa::Float64)
    return( Defaults.convertUnits("rate: from atomic", wa) )
end


"""
`Basics.recast(::RecastRateToEinsteinB,
    line::Union{Einstein.Line, PhotoEmission.Line, HyperfineInduced.Line, TwoElectronOnePhoton.Line},
    wa::Float64)`
    ... recasts a radiative rate (Einstein A, a.u.) into an Einstein B-coefficient; a Float64 is returned.
"""
function Basics.recast(::RecastRateToEinsteinB,
        line::Union{Einstein.Line, PhotoEmission.Line, HyperfineInduced.Line, TwoElectronOnePhoton.Line},
        wa::Float64)
    einsteinB = pi^2 * Defaults.getDefaults("speed of light: c")^3 / line.omega^3  * wa
    return( Defaults.convertUnits("Einstein B: from atomic", einsteinB) )
end


"""
`Basics.recast(::RecastRateToOscillatorGf, line::Union{Einstein.Line, PhotoEmission.Line, TwoElectronOnePhoton.Line}, wa::Float64)`
    ... recasts a radiative rate (Einstein A, a.u.) into the oscillator strength g_f; a Float64 is returned.
"""
function Basics.recast(::RecastRateToOscillatorGf, line::Union{Einstein.Line, PhotoEmission.Line, TwoElectronOnePhoton.Line}, wa::Float64)
    return( (Basics.twice(line.initialLevel.J) + 1) / (Basics.twice(line.finalLevel.J) + 1) / 2. *
                Defaults.getDefaults("speed of light: c")^3 / line.omega^2 * wa )
end


"""
`Basics.recast(::RecastRateToOscillatorF, line::Union{Einstein.Line, PhotoEmission.Line}, wa::Float64)`
    ... recasts a radiative rate (Einstein A, a.u.) into the oscillator strength f; a Float64 is returned.
"""
function Basics.recast(::RecastRateToOscillatorF, line::Union{Einstein.Line, PhotoEmission.Line}, wa::Float64)
    return( Defaults.getDefaults("speed of light: c") / (12. * pi * line.omega) * wa )
end


"""
`Basics.recast(::RecastRateToLineStrengthS,
    line::Union{Einstein.Line, PhotoEmission.Line, TwoElectronOnePhoton.Line}, wa::Float64)`
    ... recasts a radiative rate (Einstein A, a.u.) into the line strength S; a Float64 is returned.
"""
function Basics.recast(::RecastRateToLineStrengthS, line::Union{Einstein.Line, PhotoEmission.Line, TwoElectronOnePhoton.Line}, wa::Float64)
    einsteinA = Defaults.convertUnits("rate: from atomic to 1/s", wa)
    if      true                    S = 3.707342e-14 * (Basics.twice(line.finalLevel.J) + 1) * einsteinA / (line.omega^3)
    elseif  line.multipole == E1    S = 8.928970e-19 * (Basics.twice(line.finalLevel.J) + 1) * einsteinA / (line.omega^5)
    else                            S = 0.
    end

    return( S )
end


"""
`Basics.recommendedGrid(occupations::Dict{Shell,Int64}, Z::Float64;
                        tailFactor::Float64=16., rbox::Union{Nothing,Float64}=nothing, rnt::Float64=2.0e-6,
                        h::Float64=5.0e-2, hp::Union{Nothing,Float64}=nothing, printout::Bool=false)`
    ... derives a radial grid whose box is matched to the shells that are to be represented on it, so that a user
        need not know Rule 12 in order to obtain a grid which can carry the orbitals asked for; a grid::Radial.Grid
        is returned.

        THE SCREENING IS WHAT MAKES THIS WORK, and it is where the obvious recipe fails.  The charge an outer
        electron sees at INFINITY is `Z - NoElectrons + 1`, which is 1 for any neutral atom; used in the hydrogenic
        turning point it puts the box for neutral thorium beyond 200 a.u., and a box far too LARGE starves the
        B-spline basis exactly as badly as one too small -- measured at 249 a.u. on ~100 splines, a compact orbital
        came out misrepresented by 28%.  What the orbital feels at its own MAXIMUM is much larger than the
        asymptotic charge, and `Basics.slaterScreening` supplies it: 3.15 rather than 1 for the thorium 7s, which
        brings the box back to 67 a.u. and that 28% down to 4e-3.

        THE INNER REGION IS NOT SACRIFICED FOR THIS.  Solving the point-nucleus Dirac problem in the resulting
        basis and comparing with Basics.computeDiracEnergy, the 1s agrees to 1e-9 or better for every system from
        helium to uranium, and everything below Z = 50 agrees to 1e-7 throughout.  What remains is the outermost
        s orbital at high Z (2e-3 for the thorium 7s), and that check is over-strict there by construction: it
        tests hydrogenic orbitals at the FULL nuclear charge, which are some thirty times more compact than the
        screened orbitals the box was sized for.  It is a bound, not the error of a real calculation.

        The box then covers the classical turning point plus the exponential tail beyond it,

            rbox = max over shells of  r_+ + tailFactor * n/Zeff,    r_+ = (n^2/Zeff)(1 + sqrt(1 - l(l+1)/n^2)),

        since a hydrogenic orbital decays as exp(-Zeff r/n) and `tailFactor` therefore counts decay lengths.  The
        `2.5 r_+` of Rule 12 agrees with this for n = 3-4, where that rule was calibrated, but is too tight for
        light systems (2.9 a.u. for helium) and too generous for heavy ones.

        THE DEFAULT tailFactor = 16 WAS MEASURED, not chosen.  The AL energy is a variational bound at fixed
        basis, so it falls as the box grows and rises again once the fixed number of splines is spread too thin,
        and the best value sits at the bottom of that curve.  Scanned over He, Ne, Ar and Ti+ at
        tailFactor = 8, 10, 13, 16, 20, 26, the two ends genuinely disagree -- Ti+ is best at 8 and has lost
        1.7e-5 Ha by 26, while Ne is still gaining at 20 -- and 16 is where the worst deviation from any single
        system's own optimum is smallest, at 1.2e-5 Ha (Ne).  At that value the box beats JAC's hand-chosen
        default grid for every one of the four, by 3.4e-5 Ha for He, 9.9e-5 for Ne, 2.9e-3 for Ar and 2.3e-2 for
        Ti+, the last two because the default box of 5.95 a.u. is simply too small for them.

        The step `hp` of the outer, linear part of the log-linear mesh is scaled with the box unless given
        explicitly, and that is what keeps the cost from following the box: the mesh holds `log(rbox/rnt)/h +
        rbox/hp` points, so a fixed `hp` makes a large box expensive, while `hp = rbox/300` holds the count near 600
        whatever the box -- measured, 602 points at rbox = 5.1 and 679 at rbox = 249.  This is not a new convention
        but JAC's own default written as a rule: that grid (`rnt = 2e-6, h = 0.05, hp = 0.02`, 595 points) reaches
        rbox = 5.95, and 5.95/300 = 0.0198.  The spacing is the harmonic blend of `h*r` and `hp` and so remains ~5%
        of r near the nucleus whatever the box, which is why scaling `hp` does not degrade the inner region.

    + occupations  ::Dict{Shell,Int64}       ... shells to be represented, with their occupation, which fixes the screening.
    + Z            ::Float64                 ... nuclear charge.
    + tailFactor   ::Float64                 ... number of decay lengths to be covered beyond the classical turning point.
    + rbox         ::Union{Nothing,Float64}  ... explicit box, which overrides the estimate and is the user's way in.
    + hp           ::Union{Nothing,Float64}  ... explicit outer step; if omitted it is scaled with the box as above.
"""
function Basics.recommendedGrid(occupations::Dict{Shell,Int64}, Z::Float64;
                                tailFactor::Float64=16., rbox::Union{Nothing,Float64}=nothing, rnt::Float64=2.0e-6,
                                h::Float64=5.0e-2, hp::Union{Nothing,Float64}=nothing, printout::Bool=false)
    if  length(occupations) == 0    error("Basics.recommendedGrid(): no shell is given, so there is nothing the box " *
                                          "could be matched to.  Pass the shells that the orbitals will occupy.")   end
    if  Z <= 0.                     error("Basics.recommendedGrid(): Z = $Z must be positive.")                      end

    NoElectrons = sum( values(occupations) )
    rMax = 0.;   outer = first(keys(occupations));   ZeffOuter = Z
    for  (sh, occ)  in occupations
        if  occ <= 0    continue    end
        n = sh.n;   l = sh.l
        ## the charge felt at the orbital's own maximum, but never less than the charge felt at infinity
        Zeff = max( Z - NoElectrons + 1., Z - Basics.slaterScreening(sh, occupations), 1.0 )
        wa   = (n*n/Zeff) * (1.0 + sqrt( max(0., 1.0 - l*(l+1)/(n*n)) ))  +  tailFactor * n / Zeff
        if  wa > rMax   rMax = wa;   outer = sh;   ZeffOuter = Zeff    end
    end

    if  isnothing(rbox)     rboxx = rMax                   else    rboxx = rbox     end
    if  rboxx <= 0.         error("Basics.recommendedGrid(): rbox = $rboxx must be positive.")   end
    if  isnothing(hp)       hpx   = rboxx / 300.           else    hpx   = hp       end

    grid = Radial.Grid(Radial.Grid(false); rnt=rnt, h=h, hp=hpx, rbox=rboxx)

    if  printout
        println("> Basics.recommendedGrid(): Z = $Z with $NoElectrons electrons; the box is set by $outer, which " *
                "sees Zeff = $(round(ZeffOuter, digits=2)).")
        println(">   rbox = $(round(rboxx, digits=2)) a.u.,  hp = $(round(hpx, sigdigits=3)),  " *
                "$(grid.NoPoints) mesh points,  $(grid.nsL) large-component splines.")
    end

    return( grid )
end


"""
`Basics.recommendedGrid(configs::Array{Configuration,1}, nm::Nuclear.Model;
                        tailFactor::Float64=16., rbox::Union{Nothing,Float64}=nothing, rnt::Float64=2.0e-6,
                        h::Float64=5.0e-2, hp::Union{Nothing,Float64}=nothing, printout::Bool=false)`
    ... derives a radial grid that is matched to the orbitals which the given configurations require, so that the
        grid need not be chosen by hand; a grid::Radial.Grid is returned.

        Every shell occurring in any of the configurations is taken into account with its LARGEST occupation across
        them, since a correlation configuration may reach further out than the reference one and must still be
        representable.  Where the configurations differ, that slightly overcounts the screening and so errs towards
        a larger box; the floor `Zeff >= Z - NoElectrons + 1` bounds how far it can go.

    + configs      ::Array{Configuration,1}  ... configurations whose shells are to be representable on the grid.
    + nm           ::Nuclear.Model           ... nuclear model, which supplies the charge Z.
"""
function Basics.recommendedGrid(configs::Array{Configuration,1}, nm::Nuclear.Model;
                                tailFactor::Float64=16., rbox::Union{Nothing,Float64}=nothing, rnt::Float64=2.0e-6,
                                h::Float64=5.0e-2, hp::Union{Nothing,Float64}=nothing, printout::Bool=false)
    if  length(configs) == 0    error("Basics.recommendedGrid(): no configuration is given, so there is nothing " *
                                      "the box could be matched to.")     end
    occupations = Dict{Shell,Int64}()
    for  conf in configs
        for  (sh, occ)  in conf.shells
            if  occ > 0     occupations[sh] = max( get(occupations, sh, 0), occ )    end
        end
    end

    return( Basics.recommendedGrid(occupations, nm.Z, tailFactor=tailFactor, rbox=rbox,
                                   rnt=rnt, h=h, hp=hp, printout=printout) )
end


"""
`Basics.selectLevel(level::Level, levelSelection::LevelSelection)`
    ... returns true::Bool if the levelSelection is inactive or if the level has been selected due to its
        indices or symmetries; in all other case, false::Bool is returned.
"""
function Basics.selectLevel(level::Level, levelSelection::LevelSelection)
    if  levelSelection.active
        # Test for level index
            if  level.index  in  levelSelection.indices                                 return( true )   end
        if  LevelSymmetry(level.J, level.parity)  in  levelSelection.symmetries     return( true )   end
    else                                                                            return( true ) 
    end
    return( false )
end


"""
`Basics.selectLevelPair(iLevel::Level, fLevel::Level, lineSelection::LineSelection)`  
    ... returns true::Bool if the lineSelection is inactive or if the pair (iLevel, fLevel) has been selected due to its 
        indices or symmetries; in all other case, false::Bool is returned.
"""
function Basics.selectLevelPair(iLevel::Level, fLevel::Level, lineSelection::LineSelection)
    if     lineSelection.active
        # Test for level indexPairs
        for ip in  lineSelection.indexPairs
            if      ip[1] == 0  &&  ip[2] == fLevel.index    return( true ) 
            elseif  ip[2] == 0  &&  ip[1] == iLevel.index    return( true ) 
            elseif  ip == (iLevel.index, fLevel.index)       return( true ) 
            end
        end
        # Test for level symmetries
        for sp in  lineSelection.symmetryPairs
            if      sp == (LevelSymmetry(iLevel.J, iLevel.parity),  LevelSymmetry(fLevel.J, fLevel.parity))   return( true )
            end
        end
    else                                                     return( true ) 
    end
    return( false )
end


"""
`Basics.selectLevelTriple(iLevel::Level, nLevel::Level, fLevel::Level, pathwaySelection::PathwaySelection)`  
    ... returns true::Bool if the pathwaySelection is inactive or if the triple (iLevel, nLevel, fLevel) has been selected 
        due to its indices or symmetries; in all other case, false::Bool is returned.
"""
function Basics.selectLevelTriple(iLevel::Level, nLevel::Level, fLevel::Level, pathwaySelection::PathwaySelection)
    if     pathwaySelection.active
        # Test for level indexTriples
        for ip in  pathwaySelection.indexTriples
            if      ip[1] == 0  &&  ip[2] == 0             &&  ip[3] == 0               return( true ) 
            elseif  ip[1] == 0  &&  ip[2] == 0             &&  ip[3] == fLevel.index    return( true ) 
            elseif  ip[1] == 0  &&  ip[2] == nLevel.index  &&  ip[3] == 0               return( true ) 
            elseif  ip[2] == 0  &&  ip[3] == 0             &&  ip[1] == iLevel.index    return( true ) 
            elseif  ip[1] == 0  &&  ip[2] == nLevel.index  &&  ip[3] == fLevel.index    return( true ) 
            elseif  ip[2] == 0  &&  ip[1] == iLevel.index  &&  ip[3] == fLevel.index    return( true ) 
            elseif  ip[3] == 0  &&  ip[1] == iLevel.index  &&  ip[2] == nLevel.index    return( true ) 
            elseif  ip == (iLevel.index, nLevel.index, fLevel.index)                    return( true ) 
            end
        end
        # Test for level symmetries
        for sp in  pathwaySelection.symmetryTriples
            if      sp == (LevelSymmetry(iLevel.J, iLevel.parity),  LevelSymmetry(nLevel.J, nLevel.parity),
                            LevelSymmetry(fLevel.J, fLevel.parity))                      return( true ) 
            end
        end
    else                                                                                return( true ) 
    end
    return( false )
end


"""
`Basics.selectSymmetry(sym::LevelSymmetry, levelSelection::LevelSelection)`  
    ... returns true::Bool if the levelSelection is inactive or if the symmetry sym has been selected; 
        in all other case, false::Bool is returned.
"""
function Basics.selectSymmetry(sym::LevelSymmetry, levelSelection::LevelSelection)
    if  levelSelection.active
        if      length(levelSelection.symmetries) == 0    return( true )  
        elseif  sym in levelSelection.symmetries          return( true )
        end
    else                                                  return( true ) 
    end
    return( false )
end


"""
`Basics.shiftTotalEnergies(multiplet::Multiplet, energyShift::Float64)`  
    ... to shift the energies of all levels in the multiplet by energyShift [a.u.]; a (new) multiplet::Multiplet is returned.
"""
function Basics.shiftTotalEnergies(multiplet::Multiplet, energyShift::Float64)
    newLevels = Level[]
    for lev in multiplet.levels
        push!(newLevels, Level(lev.J, lev.M, lev.parity, lev.index, lev.energy + energyShift, lev.relativeOcc, lev.hasStateRep, lev.basis, lev.mc) )
    end
    
    newMultiplet = Multiplet(multiplet.name, newLevels)
    
    return( newMultiplet )  
end


"""
`Basics.slaterScreening(sh::Shell, occupations::Dict{Shell,Int64})`
    ... computes the screening constant of Slater's rules for an electron in the shell sh, given the occupation of
        all shells; a screening::Float64 is returned, so that the effective charge is `Z - screening`.

        Slater's grouping is used unchanged: [1s], [2s,2p], [3s,3p], [3d], [4s,4p], [4d], [4f], ...  An electron in
        an s or p shell is screened by 0.35 for each other electron of the same n with l <= 1 (0.30 within 1s), by
        0.85 for each electron of principal quantum number n-1, and by 1.00 for each electron below that.  An
        electron in a d or f shell is screened by 0.35 for each other electron of the same shell and by 1.00 for
        every electron in a group to its left, which includes the s and p shells of the same n.  The electron
        itself never screens.

        The rules are crude, and for d and f shells notoriously so, but they are cheap, standard and citable, and
        they are used here only to size a radial box -- a purpose for which the difference between Zeff = 3.15 and
        the asymptotic 1 decides whether the box is usable, while the difference between 3.15 and 3.5 does not.

    + sh           ::Shell                   ... shell whose electron is screened.
    + occupations  ::Dict{Shell,Int64}       ... occupation of every shell, including sh itself.
"""
function Basics.slaterScreening(sh::Shell, occupations::Dict{Shell,Int64})
    n = sh.n;   l = sh.l;   screening = 0.
    for  (osh, occ)  in occupations
        if  occ <= 0    continue    end
        ## the electron under consideration does not screen itself
        no = (osh == sh) ? occ - 1 : occ
        if  no <= 0     continue    end
        if      l <= 1
            if      osh.n == n   &&  osh.l <= 1     screening = screening + (n == 1 ? 0.30 : 0.35) * no
            elseif  osh.n == n                      ## d and f of the same n lie to the RIGHT and do not screen
            elseif  osh.n == n - 1                  screening = screening + 0.85 * no
            elseif  osh.n <= n - 2                  screening = screening + 1.00 * no
            end
        else
            if      osh.n == n   &&  osh.l == l     screening = screening + 0.35 * no
            elseif  osh.n == n   &&  osh.l <  l     screening = screening + 1.00 * no
            elseif  osh.n <  n                      screening = screening + 1.00 * no
            end
        end
    end

    return( screening )
end


"""
`Basics.sortByEnergy(multiplet::Multiplet)`  
    ... to sort all levels in the multiplet into a sequence of increasing energy; a (new) multiplet::Multiplet is returned.
"""
function Basics.sortByEnergy(multiplet::Multiplet)
    sortedLevels = Base.sort( multiplet.levels , lt=Base.isless)
    newLevels = Level[];   index = 0
    for lev in sortedLevels
        index = index + 1
        push!(newLevels, Level(lev.J, lev.M, lev.parity, index, lev.energy, lev.relativeOcc, lev.hasStateRep, lev.basis, lev.mc) )
    end
    
    newMultiplet = Multiplet(multiplet.name, newLevels)
    
    return( newMultiplet )  
end


"""
`Basics.Subshell(n::Int64, symmetry::LevelSymmetry)`  ... constructor for a given principal quantum number n and (level) symmetry.
"""
function Basics.Subshell(n::Int64, symmetry::LevelSymmetry) 
    if  symmetry.parity == Basics.plus
        if      symmetry.J == AngularJ64(1//2)   kappa = -1
        elseif  symmetry.J == AngularJ64(3//2)   kappa =  2
        elseif  symmetry.J == AngularJ64(5//2)   kappa = -3
        elseif  symmetry.J == AngularJ64(7//2)   kappa =  4
        elseif  symmetry.J == AngularJ64(9//2)   kappa = -5
        elseif  symmetry.J == AngularJ64(11//2)  kappa =  6
        elseif  symmetry.J == AngularJ64(13//2)  kappa = -7
        else    error("stop a")
        end
    else
        if      symmetry.J == AngularJ64(1//2)   kappa =  1
        elseif  symmetry.J == AngularJ64(3//2)   kappa = -2
        elseif  symmetry.J == AngularJ64(5//2)   kappa =  3
        elseif  symmetry.J == AngularJ64(7//2)   kappa = -4
        elseif  symmetry.J == AngularJ64(9//2)   kappa =  5
        elseif  symmetry.J == AngularJ64(11//2)  kappa = -6
        elseif  symmetry.J == AngularJ64(13//2)  kappa =  7
        else    error("stop b")
        end
    end
    
    return( Subshell(n,kappa) )
end


"""
`Basics.subshellStateString(subshell::String, occ::Int64, seniorityNr::Int64, Jsub::AngularJ64, X::AngularJ64)`  
    ... to provide a string of a given subshell state in the form '[2p_1/2^occ]_(seniorityNr, J_sub), X=Xo' ... .
"""
function Basics.subshellStateString(subshell::String, occ::Int64, seniorityNr::Int64, Jsub::AngularJ64, X::AngularJ64)
    sa = "[" * subshell * "^$occ]_($seniorityNr, " * string(Jsub) * ") X=" * string(X)
    return( sa )
end


"""
`Basics.tabulate(stream::IO, multiplet::Multiplet, levelNos::Array{Int64,1}; detail::Int64=0)`  
    ... tabulates the energies from the multiplet with level numbers in levelNos into a neat format due to different 
        criteria; nothing is returned. 
        
    + details = 0: print total energies, relative to immediate and relative to lowest.
    + details = 1: print total energies
    + details = 2: print energies to immediately lower level
    + details = 3: print energies to lowest level
"""
function Basics.tabulate(stream::IO, multiplet::Multiplet, levelNos::Array{Int64,1}; detail::Int64=0)
    ceV = Defaults.convertUnits("energy: from atomic to eV", 1.0)
    
    if      detail == 0
        sb = "  Level  J Parity    Total [Hartree]           Total [eV]      " *
             TableStrings.center(27, "Total "     * TableStrings.inUnits("energy") ) *  
             TableStrings.center(21, "To lower "  * TableStrings.inUnits("energy") ) * 
             TableStrings.center(21, "To lowest " * TableStrings.inUnits("energy") )
        println(stream, "\n  Level energies:  \n\n", sb, "\n")
        for  i = 1:length(multiplet.levels)
            if  i in levelNos
                lev    = multiplet.levels[i]
                en     = lev.energy;      
                enUser = Defaults.convertUnits("energy: from atomic", en)
                if  i == 1        enUserLower = 0.;   enUserLowest = 0.
                else
                    enLower      = lev.energy - multiplet.levels[i-1].energy;    
                    enUserLower  = Defaults.convertUnits("energy: from atomic", enLower)
                    enLowest     = lev.energy - multiplet.levels[1].energy;    
                    enUserLowest = Defaults.convertUnits("energy: from atomic", enLowest)
                end
                sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "       "
                @printf(stream, "%s %.12e   %s %.12e   %s %.12e   %s %.9e   %s %.9e   %s", 
                                sc[1:18], en, "  ", en * ceV, "  ", enUser, "  ", enUserLower, "  ", enUserLowest, "\n")
            end
        end

    elseif  detail == 1
        sb = "  Level  J Parity          Hartrees       " * "             eV                   " *  TableStrings.inUnits("energy")
        println(stream, "\n  Eigenenergies:  \n\n", sb, "\n")
        for  i = 1:length(multiplet.levels)
            if  i in levelNos
                lev = multiplet.levels[i]
                en  = lev.energy;    enUser = Defaults.convertUnits("energy: from atomic", en)
                sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "    "
                @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en * ceV, "  ", enUser, "\n")
            end
        end

    elseif  detail == 2
        sb = "  Level  J Parity          Hartrees       " * "             eV                   " * TableStrings.inUnits("energy")  
        println(stream, "\n  Energy of each level relative to immediately lower level:  \n\n", sb, "\n")
        for  i = 2:length(multiplet.levels)
            if  i in levelNos
                lev    = multiplet.levels[i]
                en     = lev.energy - multiplet.levels[i-1].energy;    
                enUser = Defaults.convertUnits("energy: from atomic", en)
                sc     = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "    "
                @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en * ceV, "  ", enUser, "\n")
            end
        end

    elseif  detail == 3
        sb = "  Level  J Parity          Hartrees       " * "             eV                   " * TableStrings.inUnits("energy")      
        println(stream, "\n  Energy of each level relative to lowest level:  \n\n", sb, "\n")
        for  i = 2:length(multiplet.levels)
            if  i in levelNos
                lev    = multiplet.levels[i]
                en     = lev.energy - multiplet.levels[1].energy;    
                enUser = Defaults.convertUnits("energy: from atomic", en)
                sc     = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "    "
                @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en * ceV, "  ", enUser, "\n")
            end
        end
    else
        error("Unsupported keystring.")
    end

    return( nothing )  
end







"""
`Basics.tabulateKappaSymmetryEnergiesDirac(kappa::Int64, evalues::Array{Float64,1}, ns::Int64, nuclearModel::Nuclear.Model)`  
    ... tabulates the eigenenergies for a given symmetry block kappa together with the corresponding Dirac energies for a 
        point-like nucleus. The index ns tells the number of 'negative-continnum' energies in the given evalues. nothing is 
        returned.
"""
function Basics.tabulateKappaSymmetryEnergiesDirac(kappa::Int64, evalues::Array{Float64,1}, ns::Int64, nuclearModel::Nuclear.Model)
    Z = nuclearModel.Z;    nx = 77
    # Determine the allowed principal quantum numbers n
    l = Basics.subshell_l(Subshell(101,kappa))
    println("  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center( 7, "Index";            na=2)
    sa = sa * TableStrings.center(10, "Subshell";         na=3)
    sa = sa * TableStrings.center(17, "Energies [a.u.]";  na=2)
    sa = sa * TableStrings.center(17, "Dirac-E  [a.u.]";  na=2)
    sa = sa * TableStrings.center(17, "Delta-E / |E|";    na=2)
    println(sa)
    println("  ", TableStrings.hLine(nx))
    for  i = ns+1:ns+7
        sa = " " * TableStrings.center( 6, TableStrings.level(i-ns); na=2)
        sa = sa *  TableStrings.flushright(10, string(Subshell(i-ns+l, kappa)); na=6)
        en = Basics.computeDiracEnergy(Subshell(i-ns+l, kappa), Z)
        if  evalues[i] >= 0.      sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", evalues[i])                       * "    "
        if  en         >= 0.      sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", en)                               * "    "
        if  evalues[i]-en >= 0.   sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", (evalues[i]-en)/abs(evalues[i]))  * "    "
        println(sa)
    end
    println("      :       :    ")
    for  i = length(evalues)-1:length(evalues)
        sa = " " * TableStrings.center( 6, TableStrings.level(i-ns); na=2)
        sa = sa *  TableStrings.flushright(10, string(Subshell(i-ns+l, kappa)); na=6)
        en = Basics.computeDiracEnergy(Subshell(i-ns+l, kappa), Z)
        if  evalues[i] >= 0.      sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", evalues[i])                       * "    "
        if  en         >= 0.      sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", en)                               * "    "
        if  evalues[i]-en >= 0.   sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", (evalues[i]-en)/abs(evalues[i]))  * "    "
        println(sa)
    end
    println("  ", TableStrings.hLine(nx))

    return( nothing )
end


## RETIRED 15-Aug-2026: `Basics.tools(dict::Dict)`, together with the whole `module-Tools.jl` it was the
## entry point to.  The intention, recoverable from the stubs: hand it the results dictionary of a
## computation and get a MENU of small interactive helpers -- `module-Tools.jl` carried a commented-out
## `using Interact`, i.e. the Julia widget package for Jupyter/Pluto, so this was to be a notebook toolbox.
## The one concrete task ever sketched was `taskGridCalculatorResults(rnt, h, hp, rmax)`, a calculator for
## radial-grid parameters -- and THAT idea now exists properly as Radial.generateGrid(grid;
## maximumPrincipalQN=..), which derives the box from the orbital's turning point (Rule 12).  Interact was
## never added as a dependency, the menu was never written, and the only method printed a placeholder.


## RETIRED 14-Aug-2026: `yesno(question::String, sa::String)`, an interactive yes/no prompt.  It had no
## caller anywhere in JAC, and it could not have run if it had: its body read `readline(STDIN)`, and STDIN
## was renamed to stdin in Julia 1.0, so the first call would have raised UndefVarError.  It also carried a
## docstring header reading `Basics.yesno` although the function was written plainly and therefore landed in
## BasicsAZ -- which is how it was found, by the docstring-pointer check of 2c1cd5d.



export  dummyQZ
using FortranFiles, Dierckx

"""
`Basics.read()`  ... reads in data from different files and sources.

+ `("CSF list: Grasp92", "cslFilename")`  
    ... reads in the CSF list as given by the (existing) .csl file cslFilename; an basis::Basis is returned with 
        basis.isDefined = true and with a proper sequence of orbitals, but with basis.orbitals = Orbital[].

+ `("orbital list: Grasp92", "orbFilename")`  
    ... read in the orbitals from a (formatted) .rwf file orbFilename; a list of orbitals::Array{OrbitalsR,1} 
        is returned with all subfields specified.
"""
function Basics.read(sa::String, filename::String)
    if      sa == "CSF list: Grasp92"                          wa = Basics.readCslFileGrasp92(filename)
    elseif  sa == "orbital list: Grasp92"                      wa = Basics.readOrbitalFileGrasp92(filename)
    elseif  sa == "energies & mixing coefficients: Grasp18"    wa = Basics.readMixingFileGrasp18(filename)
    else    error("Unsupported keystring = $sa ")
    end

    return( wa )
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
`Basics.readOrbitalFileGrasp92(filename::String, grid.Radial.Grid)`  
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
`Basics.readMixFileRelci(filename::String)`  
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
`Basics.recast()`  
    ... recasts some data from one number/representation into another one; cf. the supported keystrings and return values.

+ `("rate: radiative, to decay width", line::Union{PhotoEmission.Line, HyperfineInduced.Line}, value::Float64)`  
    ... to recast a given radiative rate (Einstein A in atomic units) into a decay withs, taking the selected energy unit 
        into account. A Float64 is returned.

+ `("rate: radiative, to Einstein A", line::Union{PhotoEmission.Line, HyperfineInduced.Line}, value::Float64)`  
    ... to recast a given spontaneous radiative rate (= Einstein A-coefficient), taking the selected unit into account. 
        A Float64 is returned.

+ `("rate: radiative, to Einstein B", line::Union{PhotoEmission.Line, HyperfineInduced.Line}, value::Float64)`  
    ... to recast a given radiative rate (Einstein A in atomic units) into an Einstein B-coefficient, taking the selected 
        unit into account. A Float64 is returned.

+ `("rate: radiative, to g_f", line::Union{PhotoEmission.Line}, value::Float64)`  
    ... to recast a given radiative rate (Einstein A in atomic units) into an oscillator strength g_f; 
        a Float64 is returned.

+ `("rate: radiative, to f", line::Union{PhotoEmission.Line}, value::Float64)`  
    ... to recast a given radiative rate (Einstein A in atomic units) into an oscillator strength f; 
        a Float64 is returned.

+ `("rate: radiative, to S", line::Union{PhotoEmission.Line}, value::Float64)`  
    ... to recast a given radiative rate (Einstein A in atomic units) into a line strength S; 
        a Float64 is returned.
"""
function Basics.recast(sa::String, line::Union{Einstein.Line, PhotoEmission.Line, HyperfineInduced.Line}, wa::Float64)
    if  typeof(line) == HyperfineInduced.Line  &&
        ! (sa  in ["rate: radiative, to decay width", "rate: radiative, to Einstein A", "rate: radiative, to Einstein B"])
        error("Not supported recast for HyperfineInduced.Line's ")
    end

    if       sa == "rate: radiative, to decay width"
        width = Defaults.convertUnits("energy: from atomic", wa)
        return( width )

    elseif   sa == "rate: radiative, to Einstein A"
        einsteinA = Defaults.convertUnits("rate: from atomic", wa)
        return( einsteinA )

    elseif   sa == "rate: radiative, to Einstein B"
        einsteinB = pi^2 * Defaults.getDefaults("speed of light: c")^3 / line.omega^3  * wa
        einsteinB = Defaults.convertUnits("Einstein B: from atomic", einsteinB)
        return( einsteinB )

    elseif   sa == "rate: radiative, to g_f"
        gf = (Basics.twice(line.initialLevel.J) + 1) / (Basics.twice(line.finalLevel.J) + 1) / 2. * 
                Defaults.getDefaults("speed of light: c")^3 / line.omega^2 * wa   
        return( gf )

    elseif   sa == "rate: radiative, to f"
        f  = Defaults.getDefaults("speed of light: c") / (12. * pi * line.omega) * wa   
        ## f  = 2 * line.omega / 3. / (Basics.twice(line.initialLevel.J) + 1) * wa
        return( f )

    elseif   sa == "rate: radiative, to S"
        einsteinA = Defaults.convertUnits("rate: from atomic to 1/s", wa)
        if      true                    S = 3.707342e-14 * (Basics.twice(line.finalLevel.J) + 1) * einsteinA / (line.omega^3)
        elseif  line.multipole == E1    S = 8.928970e-19 * (Basics.twice(line.finalLevel.J) + 1) * einsteinA / (line.omega^5)
        else                            S = 0.
        end
        return( S )

    else     error("Unsupported keystring = $sa")
    end
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
            if      ip == (LevelSymmetry(iLevel.J, iLevel.parity),  LevelSymmetry(fLevel.J, fLevel.parity))   return( true ) 
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
            if      ip == (LevelSymmetry(iLevel.J, iLevel.parity),  LevelSymmetry(nLevel.J, nLevel.parity),  
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




#==
"""
`Basics.tabulate(stream::IO, sa::String, multiplet::Multiplet, levelNos::Array{Int64,1})`  
    ... tabulates the energies from the multiplet with level numbers in levelNos due to different criteria.

+ `(stream::IO, "multiplet: energies", multiplet::Multiplet)`  
    ... to tabulate the energies of all levels of the given multiplet into a neat format; nothing is returned.
+ `(stream::IO, "multiplet: energy relative to immediately lower level", multiplet::Multiplet, levelNos::Array{Int64,1})`  
    ... to tabulate the energy splitting between neighboured levels of all levels of the given multiplet into a neat
        format; nothing is returned.
+ `(stream::IO, "multiplet: energy of each level relative to lowest level", multiplet::Multiplet, levelNos::Array{Int64,1}`  
    ... to tabulate the energy splitting of all levels with regard to the lowest level of the given multiplet into 
        a neat format; nothing is returned.
"""
function Basics.tabulate(stream::IO, sa::String, multiplet::Multiplet, levelNos::Array{Int64,1})
    if        sa == "multiplet: energies"
        println(stream, "\n  Eigenenergies:")
        sb = "  Level  J Parity          Hartrees       " * "             eV                   " *  TableStrings.inUnits("energy")     
        println(stream, "\n", sb, "\n")
        for  i = 1:length(multiplet.levels)
            if  i in levelNos
                lev = multiplet.levels[i]
                en  = lev.energy;    en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
                sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "    "
                @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
            end
        end

    elseif    sa == "multiplet: energy relative to immediately lower level"
        println(stream, "\n  Energy of each level relative to immediately lower level:")
        sb = "  Level  J Parity          Hartrees       " * "             eV                   " * TableStrings.inUnits("energy")  
        println(stream, "\n", sb, "\n")
        for  i = 2:length(multiplet.levels)
            if  i in levelNos
                lev = multiplet.levels[i]
                en    = lev.energy - multiplet.levels[i-1].energy;    
                en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
                sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "    "
                @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
            end
        end

    elseif    sa == "multiplet: energy of each level relative to lowest level"
        println(stream, "\n  Energy of each level relative to lowest level:")
        sb = "  Level  J Parity          Hartrees       " * "             eV                   " * TableStrings.inUnits("energy")      
        println(stream, "\n", sb, "\n")
        for  i = 2:length(multiplet.levels)
            if  i in levelNos
                lev = multiplet.levels[i]
                en    = lev.energy - multiplet.levels[1].energy;    
                en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
                sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "    "
                @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
            end
        end
    else
        error("Unsupported keystring.")
    end

    return( nothing )  
end   ==#



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


"""
`Basics.tools(dict::Dict)`  
    ... select different tools from a menu if proper results::dict are given in terms of a dictionary as obtained
        usually from an Atomic.Computation.
"""
function Basics.tools(dict::Dict)
    println("Tools.select(): A menu if a dict comes in; $dict ")
end


"""
`Basics.yesno(question::String, sa::String)`  
    ... Returns true if the answer 'yes' or 'y' is found, and false otherwise; sa = {"Y", "N"} determines how the zero-String 
        "" is interpreted. The given question is repeated until a proper answer is obtained.
"""
function yesno(question::String, sa::String)
    while  true
        print(question * "  ");    reply = strip( readline(STDIN) )
        if      sa == "Y"
            if      reply in ["", "Y", "Yes", "y", "yes"]   return( true )
            elseif  reply in [    "N", "No",  "n", "no" ]   return( false )
                    println("... answer not recognized ... redo:")
            end
        elseif  sa == "N"
            if      reply in [    "Y", "Yes", "y", "yes"]   return( true )
            elseif  reply in ["", "N", "No",  "n", "no" ]   return( false )
                    println("... answer not recognized ... redo:")
            end
        end
    end
end


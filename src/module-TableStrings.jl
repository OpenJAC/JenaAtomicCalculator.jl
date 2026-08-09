
"""
`module JAC.TableStrings`  
... a submodel of JAC that contains various methods to quickly printout the headers of various tables, 
    especially if these headers require several lines.
"""
module TableStrings


using Printf, ..Basics,  ..Defaults


"""
`TableStrings.center(n::Int64, sa::String)`  
    ... a string sa is centered within a string of length n, and optionally, na blanks are appended; a string is returned.
"""
function  center(n::Int64, sa::String; na::Int64=0) 
    sb = sa;   mb = n - length(sa)
    while mb > 2    sb = " "*sb*" ";            mb = mb - 2    end
    if    mb == 2   sb = " "*sb*" "     else    sb = " "*sb    end
    return( sb * repeat(" ", na) )
end


"""
`TableStrings.floatList(n::Int64, floats::Array{Float64,1})`  
    ... a list of Strings is returned, each containing up to n floats of length 6.
"""
function  floatList(n::Int64, floats::Array{Float64,1})
    sa = "";   wa = String[];   nc = 0
    for  flt  in  floats
        nc = nc + 1;   sa = sa * " " * @sprintf("%.3e", flt) * " "
        if   nc == n   push!(wa, sa);    sa = "";   nc = 0   end
    end
    if  sa != ""      push!(wa, sa)                         end
    
    return( wa )
end


"""
`TableStrings.flushleft(n::Int64, sa::String)`  
    ... a string sa is left-placed within a string of length n, and optionally, na blanks are appended; a string is returned.
"""
function  flushleft(n::Int64, sa::String; na::Int64=0) 
    sb = sa;   sb = sb * repeat(" ", n+na)
    return( sb[1:n+na] )
end


"""
`TableStrings.flushright(n::Int64, sa::String)`  
    ... a string sa is right-placed within a string of length n, and optionally, na blanks are appended; a string is returned.
"""
function  flushright(n::Int64, sa::String; na::Int64=0) 
    sb = sa;   sb = repeat(" ", n) * sb;   sb = sb[end-n+1:end]
    return( sb * repeat(" ", na) )
end


"""
`TableStrings.hBlank(n::Int64)`  ... a string of n 'blanks' is returned.
"""
function  hBlank(n::Int64) 
    return( repeat(" ", n) )
end


"""
`TableStrings.hLine(n::Int64)`  ... a string of n * "-" is returned.
"""
function  hLine(n::Int64) 
    return( repeat("-", n) )
end


"""
`TableStrings.inUnits(sa::String)`  ... a string is returned that displays the requested unit in the form: '[unit]'
"""
function  inUnits(sa::String) 
    if       sa == "energy"                         return( "[" * Defaults.GBL_ENERGY_UNIT * "]" )
    elseif   sa == "rate"                           return( "[" * Defaults.GBL_RATE_UNIT   * "]" )
    elseif   sa == "time"                           return( "[" * Defaults.GBL_TIME_UNIT   * "]" )
    elseif   sa == "cross section"                  return( "[" * Defaults.GBL_CROSS_SECTION_UNIT * "]" )
    elseif   sa == "differential cross section"     return( "[" * Defaults.GBL_EDIFF_CROSS_SECTION_UNIT * "]" )
    elseif   sa == "strength"                       return( "[" * Defaults.GBL_STRENGTH_UNIT * "]" )
    elseif   sa == "reduced strength"               return( "[" * Defaults.GBL_STRENGTH_UNIT * "^2" * "]" )
    else     error("stop a")
    end
end


"""
`TableStrings.level(i::Int64)`  ... a string "....i" is returned.
"""
function  level(i::Int64)
    si = string(i);   ni = length(si);    sa = repeat(" ", 5)
    return( sa[1:5-ni] * si )
end 


"""
`TableStrings.levels_if(i::Int64, f::Int64)`  ... a string "....i -- ....f" is returned
"""
function  levels_if(i::Int64, f::Int64)
    si = string(i);   ni = length(si);    sf = string(f);   nf = length(sf);    sa = repeat(" ", 5)
    return( sa[1:5-ni] * si * " -- " * sa[1:4-nf] * sf )
end 


"""
`TableStrings.levels_imf(i::Int64, m::Int64, f::Int64)`  ... a string "....i -- ....m -- ....f" is returned
"""
function  levels_imf(i::Int64, m::Int64, f::Int64)
    si = string(i);   ni = length(si);    sm = string(m);   nm = length(sm);    sf = string(f);   nf = length(sf);    
    sa = repeat(" ", 5)
    return( sa[1:5-ni] * si * " -- " * sa[1:5-nm] * sm * " -- " * sa[1:5-nf] * sf )
end 


"""
`TableStrings.gaugeMultipolesTupels(n::Int64, mpList::Array{Tuple{Basics.UseGauge,Array{Basics.EmMultipole,1}},1})`  
    ... a string of tupels Used gauge[E1,M1,E2,E1, ...] is returned.
"""
function  gaugeMultipolesTupels(n::Int64, mpList::Array{Tuple{Basics.UseGauge,Array{Basics.EmMultipole,1}},1}) 
    sa = "";   wa = String[]
    for tt in mpList
        sa = sa * string(tt[1])[1:3] * "[" 
        for  mp in tt[2]
            sa = sa * string(mp) * ","
        end
        sa = sa[1:end-1] * "], "
        if  length(sa) + 15 > n    push!(wa, sa[1:end-2]);    sa = ""    end
    end
    if  sa != ""    push!(wa, sa[1:end-2])    end
    return( wa )
end


"""
`TableStrings.kappaSymmetryTupelList(n::Int64, kappaList::Array{Tuple{Int64,LevelSymmetry},1})`  
    ... a list of Strings with maximal length n is returned; each string in this list comprises a number of 
        'shell -> symmetry -> shell' descriptors.
"""
function kappaSymmetryTupelList(n::Int64, kappaList::Array{Tuple{Int64,LevelSymmetry},1}) 
    sa = "";   wa = String[]
    for k in kappaList
        shin = Subshell(9, k[1])
        sa = sa * string(shin)[2:end] * "(" * string(k[2]) * ")" * ",  "
        if  length(sa) + 15 > n    push!(wa, sa[1:end-3]);    sa = ""    end
    end
    if  sa != ""    push!(wa, sa[1:end-3])    end
    return( wa )
end


"""
`TableStrings.kappaKappaSymmetryTupels(n::Int64, kappaList::Array{Tuple{Int64,Int64,LevelSymmetry},1})`  
    ... a list of Strings with maximal length n is returned; each string in this list comprises a number of 
        'shell -> symmetry -> shell' descriptors.
"""
function kappaKappaSymmetryTupels(n::Int64, kappaList::Array{Tuple{Int64,Int64,LevelSymmetry},1}) 
    sa = "";   wa = String[]
    for k in kappaList
        shin = Subshell(9, k[1]);    shout = Subshell(9, k[2])
        sa = sa * string(shin)[2:end] * "(" * string(k[3]) * ")" * string(shout)[2:end]  * ",  "
        if  length(sa) + 15 > n    push!(wa, sa[1:end-3]);    sa = ""    end
    end
    if  sa != ""    push!(wa, sa[1:end-3])    end
    return( wa )
end


"""
`TableStrings.kappaKappaKappaSymmetryTupels(n::Int64, kappaList::Array{Tuple{Int64,Int64,Int64,LevelSymmetry},1})`   
    ... a list of Strings with maximal length n is returned; each string in this list comprises a number of 
        '(shell -> shell) symmetry [shell]' descriptors.
"""
function kappaKappaKappaSymmetryTupels(n::Int64, kappaList::Array{Tuple{Int64,Int64,Int64,LevelSymmetry},1}) 
    sa = "";   wa = String[]
    for k in kappaList
        shin = Subshell(9, k[1]);    shout = Subshell(9, k[2]);    shauger = Subshell(9, k[3])
        sa = sa * "(" * string(shin)[2:end] * "->" * string(shout)[2:end] * ") " * string(k[4]) * "[" * string(shauger)[2:end]  * "],  "
        if  length(sa) + 18 > n    push!(wa, sa[1:end-3]);    sa = ""    end
    end
    if  sa != ""    push!(wa, sa[1:end-3])    end
    return( wa )
end


"""
`TableStrings.kappaMultipoleSymmetryTupels(n::Int64, kappaList::Array{Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry},1})` 
    ... a list of Strings with maximal length n is returned; each string in this list comprises a number of 
        'shell (multipole, gauge, symmetry)' descriptors.
"""
function  kappaMultipoleSymmetryTupels(n::Int64, kappaList::Array{Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry},1})
    sa = "";   wa = String[]
    for k in kappaList
        sh = Subshell(9, k[1])
        sa = sa * string(sh)[2:end] * " (" * string(k[2]) * ", " * string(k[3])[1:3] * "; " * string(k[4]) * "),  "
        if  length(sa) + 15 > n    push!(wa, sa[1:end-3]);    sa = ""    end
    end
    if  sa != ""    push!(wa, sa[1:end-3])    end
    return( wa )
end


"""
`TableStrings.kappaMultipoleSymmetryPhaseTupels(n::Int64, kappaList::Array{Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry,Float64},1})` 
    ... a list of Strings with maximal length n is returned; each string in this list comprises a number of 
        'shell (multipole, gauge, symmetry, phase)' descriptors.
"""
function  kappaMultipoleSymmetryPhaseTupels(n::Int64, kappaList::Array{Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry,Float64},1})
    sa = "";   wa = String[]
    for k in kappaList
        sh = Subshell(9, k[1])
        sa = sa * string(sh)[2:end] * " (" * string(k[2]) * ", " * string(k[3])[1:3] * "; " * string(k[4]) * "; " * @sprintf("%.4e", k[5]) * "),  "
        if  length(sa) + 15 > n    push!(wa, sa[1:end-3]);    sa = ""    end
    end
    if  sa != ""    push!(wa, sa[1:end-3])    end
    return( wa )
end


"""
`TableStrings.multipolesSymmetryTupels(n::Int64, mpList::Array{Tuple{EmMultipole,LevelSymmetry,EmMultipole,EmGauge},1})` 
    ... a list of Strings with maximal length n is returned; each string in this list comprises a number of 
        'gauge (multipole, symmetry, multipole)' descriptors.
"""
function  multipolesSymmetryTupels(n::Int64, mpList::Array{Tuple{EmMultipole,LevelSymmetry,EmMultipole,EmGauge},1})
    sa = "";   wa = String[]
    for k in mpList
        sa = sa * string(k[4])[1:3] * " (" * string(k[1]) * " > " * string(k[2]) * " > " * string(k[3]) * "),  "
        if  length(sa) + 15 > n    push!(wa, sa[1:end-3]);    sa = ""    end
    end
    if  sa != ""    push!(wa, sa[1:end-3])    end
    return( wa )
end


"""
`TableStrings.processSymmetryEnergyTupels(n::Int64, pList::Array{Tuple{Basics.AbstractProcess, Int64, LevelSymmetry, Float64},1}, sc::String)`  
        ... a list of Strings with maximal length n is returned; each string in this list comprises a number of 
        'sc(process: NoElectrons, symmetry, energy),  ' descriptors.
"""
function  processSymmetryEnergyTupels(n::Int64, pList::Array{Tuple{Basics.AbstractProcess, Int64, LevelSymmetry, Float64},1}, sc::String)
    sa = "";   wa = String[]
    for p in pList
        sa = sa * sc * "[" * string(p[1])[1:1] * ": " * string(p[2]) * ", " * string(p[3]) * ", " 
        sa = sa * @sprintf("%.4e", p[4]) * "],  "
        if  length(sa) + 20 > n    push!(wa, sa[1:end-3]);    sa = ""    end
    end
    if  sa != ""    push!(wa, sa[1:end-3])    end
    return( wa )
end


"""
`TableStrings.MMffTupels(n::Int64, pList::Array{Tuple{AngularM64, AngularM64, Float64, Float64},1}, sc::String)`  
        ... a list of Strings with maximal length n is returned; each string in this list comprises a number of 
        '[M, M, float, float],  ' descriptors.
"""
function  MMffTupels(n::Int64, pList::Array{Tuple{AngularM64, AngularM64, Float64, Float64},1}, sc::String)
    sa = "";   wa = String[]
    for p in pList
        sa = sa * sc * "[" * string(p[1]) * ", " * string(p[2]) * ", "  
        sa = sa * @sprintf("%.2e", p[3]) * ", " * @sprintf("%.2e", p[4]) * "],  "
        if  length(sa) + 20 > n    push!(wa, sa[1:end-3]);    sa = ""    end
    end
    if  sa != ""    push!(wa, sa[1:end-3])    end
    return( wa )
end


"""
`TableStrings.kappaSymmetryTupels(n::Int64, kappaList::Array{Tuple{Int64,LevelSymmetry},1})`  
    ... a string of shell (Symmetry) is returned.
"""
function  kappaSymmetryTupels(n::Int64, kappaList::Array{Tuple{Int64,LevelSymmetry},1})
    sa = ""
    for k in kappaList
        sh = Subshell(9, k[1])
        sa = sa * string(sh)[2:end] * " (" * string(k[2]) * "),  "
    end
    return( sa[1:end-3] )
end


"""
`TableStrings.multipoleGaugeTupels(n::Int64, mpList::Array{Tuple{Basics.EmMultipole,Basics.EmGauge},1})`  
    ... a string of tupels E1(Gauge), M1(Gauge) is returned.
"""
function  multipoleGaugeTupels(n::Int64, mpList::Array{Tuple{Basics.EmMultipole,Basics.EmGauge},1})
    sa = ""
    for mp in mpList
        sa = sa * string(mp[1]) * "(" * string(mp[2]) * "), "
    end
    return( sa[1:end-2] )
end


"""
`TableStrings.multipoleList(mpList::Array{EmMultipole,1})`  ... a string of multipoles is returned.
"""
function  multipoleList(mpList::Array{EmMultipole,1})
    sa = ""
    for mp in mpList
        sa = sa * string(mp) * ", "
    end
    return( sa[1:end-2] )
end


"""
`TableStrings.subshellList(n::Int64, subshells::Array{Subshell,1})`  
    ... a list of Strings is returned, each containing up to n subshells.
"""
function  subshellList(n::Int64, subshells::Array{Subshell,1})
    sa = "";   wa = String[];   nc = 0
    for  subsh  in  subshells
        nc = nc + 1;   sa = sa * " " * string(subsh) * "    "
        if   nc == n   push!(wa, sa);    sa = "";   nc = 0   end
    end
    if  sa != ""      push!(wa, sa)                         end
    
    return( wa )
end


"""
`TableStrings.subshellSymmetryTupels(n::Int64, subshellList::Array{Tuple{Subshell,LevelSymmetry},1})`  
    ... a string of subshell (Symmetry) is returned.
"""
function  subshellSymmetryTupels(n::Int64, subshellList::Array{Tuple{Subshell,LevelSymmetry},1})
    sa = ""
    for k in subshellList
        sa = sa * string(k[1]) * " (" * string(k[2]) * "),  "
    end
    return( sa[1:end-3] )
end


"""
`TableStrings.symmetries_if(isym::LevelSymmetry, fsym::LevelSymmetry)`  ... a string 'J^P --> J^P' is returned
"""
function  symmetries_if(isym::LevelSymmetry,  fsym::LevelSymmetry)
    sa = string(isym) * " --> " * string(fsym)
    return( sa )
end 


"""
`TableStrings.symmetries_imf(isym::LevelSymmetry, msym::LevelSymmetry, fsym::LevelSymmetry)`  ... a string 'J^P --> J^P --> J^P' is returned
"""
function  symmetries_imf(isym::LevelSymmetry, msym::LevelSymmetry, fsym::LevelSymmetry)
    sa = string(isym) * " --> " * string(msym) * " --> " * string(fsym)
    return( sa )
end 


"""
`TableStrings.twoPhotonGaugeTupels(n::Int64, mpList::Array{Tuple{AngularJ64, Basics.EmMultipole, LevelSymmetry, Basics.EmMultipole, Basics.EmGauge},1})
    ... a string of tupels (K-rank, multipole_1, Jsym, multipole_2, gauge) is returned.
"""
function  twoPhotonGaugeTupels(n::Int64, mpList::Array{Tuple{AngularJ64, Basics.EmMultipole, LevelSymmetry, Basics.EmMultipole, Basics.EmGauge},1}) 
    sa = "";   wa = String[];   mpReduced = String[]
    for tt in mpList
        sb = "(" * string(tt[1]) * ", " * string(tt[2]) * ", " * string(tt[3]) * ", " * string(tt[4]) * ", " * string(tt[5])[1:3] * ")"
        if    sb in mpReduced    else    push!(mpReduced, sb)    end
    end
    for sb in mpReduced
        sa = sa * sb * ", "  
        if  length(sa) + 20 > n    push!(wa, sa[1:end-2]);    sa = ""    end
    end
    if  sa != ""    push!(wa, sa[1:end-2])    end
    return( wa )
end


"""
`   + (n::Int64, mpList::Array{Tuple{String, Basics.EmMultipole, LevelSymmetry, Basics.EmMultipole, Basics.EmGauge},1})
    ... a string of tupels (Sik, multipole_1, Jsym, multipole_2, gauge) is returned.
"""
function  twoPhotonGaugeTupels(n::Int64, mpList::Array{Tuple{String, Basics.EmMultipole, LevelSymmetry, Basics.EmMultipole, Basics.EmGauge},1}) 
    sa = "";   wa = String[];   mpReduced = String[]
    for tt in mpList
        sb = "(" * string(tt[1]) * ", " * string(tt[2]) * ", " * string(tt[3]) * ", " * string(tt[4]) * ", " * string(tt[5])[1:3] * ")"
        if    sb in mpReduced    else    push!(mpReduced, sb)    end
    end
    for sb in mpReduced
        sa = sa * sb * ", "  
        if  length(sa) + 20 > n    push!(wa, sa[1:end-2]);    sa = ""    end
    end
    if  sa != ""    push!(wa, sa[1:end-2])    end
    return( wa )
end

end # module

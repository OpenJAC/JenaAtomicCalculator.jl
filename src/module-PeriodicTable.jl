
"""
`module  JAC.PeriodicTable`  
... a submodel of JAC that contains methods and data from the periodic table of elements. Although we here collect 
    data for various elements, no attempt will be made to set-up a full compilation. Instead, useful information
    and some semi-empirical data are compiled to simplify the use of the JAC package in the future.
"""
module PeriodicTable


"""
`abstract type PeriodicTable.AbstractEnergyData` 
    ... defines an abstract and a number of singleton types to deal with different data set for binding energies,
        ionization potentials, etc.

    + Larkins1977       ... F. Larkins' (1977) inner-shell (1s - 4p) subshell binding energies of the neutral atom;
                            Z = 1, ..., 45, with many (Z, subshell) combinations left undefined.
    + Nist2025          ... NIST's (2025) tabulated successive ionization potentials of the neutral atom and all its
                            ions, Z = 1, ..., 92; used by Empirical.ionizationPotential and by
                            Empirical.totalEnergy(data=Nist2025()); not supported by Empirical.bindingEnergy.
    + Williams2000      ... the subshell binding energies of Williams (2000), part of the LBNL X-Ray Data Booklet's
                            Table 1-1; Z = 1, ..., 36 only, shells 1s - 4p.
    + XrayDataBooklet   ... the LBNL X-Ray Data Booklet's Table 1-1 [J. A. Bearden & A. F. Burr, Rev. Mod. Phys. 39,
                            125 (1967); M. Cardona & L. Ley, Photoemission in Solids I (1978); J. C. Fuggle &
                            N. Mårtensson, J. Electron Spectrosc. Relat. Phenom. 21, 275 (1980)] -- the same lineage
                            as Williams2000, and numerically identical wherever both are tabulated, but with the
                            fuller Z = 1, ..., 92 and shells 1s - 6p coverage; the default data set for
                            Empirical.bindingEnergy/scaledBindingEnergy and the functions built on them.
"""
abstract type  AbstractEnergyData                                    end
struct     Larkins1977                   <:  AbstractEnergyData      end
struct     Nist2025                      <:  AbstractEnergyData      end
struct     Williams2000                  <:  AbstractEnergyData      end
struct     XrayDataBooklet               <:  AbstractEnergyData      end


export  AbstractEnergyData, Larkins1977, Nist2025, Williams2000, XrayDataBooklet


"""
`abstract type PeriodicTable.AbstractRecombinationData`
    ... defines an abstract and a number of singleton types to deal with different data sets of fitted
        recombination-rate-coefficient parameters.

    + Arnaud1985   ... M. Arnaud & R. Rothenflug, Astron. Astrophys. Suppl. Ser. 60, 425 (1985); their Table IIB
                        gives the dielectronic recombination fit parameters (ADI, BDI, T0, T1) for the H-like
                        sequence, but only for 4 specific ions: O VIII (Z=8), Mg XII (Z=12), Ca XX (Z=20),
                        Fe XXVI (Z=26) -- not a general Z-scaling formula. Their He-like sequence rate is instead
                        a genuine Z-scaling closed form (no table), which is implemented directly in
                        Empirical.dielectronicRecombinationPlasmaAlpha rather than stored here.
"""
abstract type  AbstractRecombinationData                             end
struct     Arnaud1985                    <:  AbstractRecombinationData end


export  AbstractRecombinationData, Arnaud1985


"""
`abstract type PeriodicTable.AbstractYieldData`
    ... defines an abstract and a number of singleton types to deal with different data sets for fluorescence,
        Auger and Coster-Kronig yields. These yields are dimensionless probabilities and are kept apart from the
        AbstractEnergyData above.

    + KrauseAdopted2016   ... applies the K-shell fluorescence yields as adopted by the xraylib project (2016)
"""
abstract type  AbstractYieldData                                     end
struct     KrauseAdopted2016             <:  AbstractYieldData       end


export  AbstractYieldData, KrauseAdopted2016

#################################################################################################################################
#################################################################################################################################


using  Printf, ..Basics, ..ManyElectron


"""
`PeriodicTable.getAtomicNumber(symbol::Symbol)`  
    ... returns the atomic number (nuclear charge Z) for the given element symbol.
"""
function getAtomicNumber(symbol::Symbol)
    for  (k,v) in PeriodicTable.data
        if  v.symbol == symbol    return( v.nuclearCharge )    end
    end
    
    error("Undefined element symbol in the periodic table.")
end


"""
`PeriodicTable.getData()`  
    ... provides various data from the periodic table of elements. 

+ `("mass", Z::Int64)`  or  `("mass", symbol::Symbol)`   
    ... to get mean mass value of the element with nuclear charge Z or symbol.

+ `("density", Z::Int64)`  or  `("density", symbol::Symbol)`   
    ... to get the density [g/cm^3] of the element with nuclear charge Z or symbol.

+ `("1st IP", Z::Int64)`  or  `("1st IP", symbol::Symbol)`   
    ... to get the first ionization potential [eV] of the element with nuclear charge Z or symbol.

+ `("ground configuration", Z::Int64)`  or  `("ground configuration", symbol::Symbol)`   
    ... to get the ground configuration of the neutral atom with nuclear charge Z or symbol.

+ `("atomic radius", Z::Int64)`  or  `("atomic radius", symbol::Symbol)`   
    ... to get the atomic radius of the neutral atom with nuclear charge Z or symbol.

+ `("polarizibility", Z::Int64)`  or  `("polarizibility", symbol::Symbol)`   
    ... to get the polarizibility of the neutral atom with nuclear charge Z or symbol.

+ `("ground configuration", Z::Int64)`  or  `("ground configuration", symbol::Symbol)`   
    ... to get the ground configuration of the neutral atom with nuclear charge Z or symbol.

+ `("chemical hardness", Z::Int64)`  or  `("chemical hardness", symbol::Symbol)`   
    ... to get the chemical hardness of the neutral atom with nuclear charge Z or symbol.
"""
function getData(sa::String, zsymbol::Union{Int64,Symbol})
    if  typeof(zsymbol) == Symbol   Z = PeriodicTable.getAtomicNumber(zsymbol)    else    Z = zsymbol     end

    if        sa == "mass"              return ( PeriodicTable.data[Z].mass )   
    elseif    sa == "density"           println("\nDensities are given [g/cm^3]: ...")
                                        return ( PeriodicTable.data[Z].density )   
    elseif    sa == "1st IP"            println("\nIonization potentials are given [eV]: ...")
                                        return ( PeriodicTable.data[Z].firstIP )   
    elseif    sa == "ground configuration"        
                                        return ( Configuration(PeriodicTable.data[Z].configuration) )   
    elseif    sa == "atomic radius"     println("\nAtomic radius is given [Angstroem]: ...")
                                        return ( PeriodicTable.data[Z].atomicRadius )   
    elseif    sa == "polarizibility"    println("\nAtomic polarizibility is given [10^-24 c.c]: ...")
                                        return ( PeriodicTable.data[Z].polarizibility )   
    elseif    sa == "chemical hardness" println("\nChemical hardness is given [eV]: ...")
                                        return ( PeriodicTable.data[Z].chemHardness )   
    else      error("Unsupported keystring:: $sa")
    end
end


"""
`PeriodicTable.getIsotopeData()`  
    ... provides various data from the periodic table of elements. Not yet implemented.

+ `("mass", Z::Int64)`  or  `("mass", symbol::Symbol)`   
    ... to get mean mass value of the element with nuclear charge Z or symbol

"""
function getIsotopeData(sa::String, zsymbol::Union{Int64,Symbol})

    if        sa in ["alpha", "fine-structure constant alpha"]          return (FINE_STRUCTURE_CONSTANT)   
    else      error("Unsupported keystring:: $sa")
    end
end



########################################################################################################################
########################################################################################################################
########################################################################################################################



"""
`struct  PeriodicTable.Isotope`  ... defines a type for collecting data for an individual isotope.

    + massNumber         ::Float64           ... Mass number of the isotope.
    + isotopeMass        ::Float64           ... Mass of the isotope.
    + spin               ::AngularJ64        ... Nuclear spin of the isotope
    + mu                 ::Float64           ... Magnetic dipole moment of the isotope in nuclear magnetons.
    + Q                  ::Float64           ... Electric quadrupole moment of the isotope.
"""
struct  Isotope 
    massNumber           ::Float64     
    isotopeMass          ::Float64  
    spin                 ::AngularJ64
    mu                   ::Float64   
    Q                    ::Float64 
end


"""
`PeriodicTable.Isotope()`  ... simple constructur for an empty instance of an isotope.
"""
function Isotope(i::Int64)
    Isotope( i, NaN, AngularJ64(0), NaN, NaN)
end


"""
`struct  PeriodicTable.Element`  ... defines a type for collecting data for each element.

    + nuclearCharge     ::Int64                 ... Nuclear charge
    + symbol            ::Symbol                ... Symbol of the element.
    + name              ::String                ... Name of the element.
    + mass              ::Float64               ... Mean mass number of the element.
    + density           ::Float64               ... Density of the element [g/cm^3].
    + firstIP           ::Float64               ... First ionization potential [eV].
    + configuration     ::String                ... Ground configuration of the element.
    + atomicRadius      ::Float64               ... Mean atomic radius of the neutral atom (Ghosh & Biswas, 2002).
    + polarizibility    ::Float64               ... Mean polarizibility of the neutral atom (Ghosh & Biswas, 2002).
    + chemHardness      ::Float64               ... Mean chemical hardness of the element (Ghosh & Biswas, 2002).
"""
struct Element
    nuclearCharge       ::Int64 
    symbol              ::Symbol
    name                ::String
    mass                ::Float64
    density             ::Float64 
    firstIP             ::Float64  
    configuration       ::String
    atomicRadius        ::Float64  
    polarizibility      ::Float64 
    chemHardness        ::Float64 
end 


"""
`PeriodicTable.Element(charge::Int64, symbol::Symbol)`  ... simple constructur for defining an element
"""
function Element(charge::Int64, symbol::Symbol)
    Element(charge, symbol, "not defined", 999.99, Isotope[])
end


# `Base.show(io::IO, element::PeriodicTable.Element)`  ... prepares a proper printout of the variable element::PeriodicTable.Element.
function Base.show(io::IO, element::PeriodicTable.Element) 
    println(io, "Element  $(element.name) ($(element.symbol), Z=$(element.nuclearCharge), mass=$(element.mass))  ")
    println(io, "configuration:            $(element.configuration)  ")
end


"""
`PeriodicTable.data::Dict{Int64,PeriodicTable.Element}`
    ... provides all available data for an element with given nuclear charge Z::Int64 by PeriodicTable.data[Z];
    cf. PeriodicTable.Element
"""
const data = Dict{Int64,PeriodicTable.Element}(
            1 => Element(  1, :H,  "hydrogen",         1.008,   0.09,  13.60,  "1s",			            0.5292,     0.6669,   13.588    ),
            2 => Element(  2, :He, "helium",           4.003,   0.18,  24.59,  "1s^2",			            0.3113,     0.1358,   22.383    ),
            3 => Element(  3, :Li, "lithum",           6.941,   0.53,   5.39,  "[He] 2s",			        1.6282,    19.4239,    4.4164   ),
            4 => Element(  4, :Be, "beryllium",        9.012,   1.85,   9.32,  "[He] 2s^2", 		        1.0855,     5.7558,    6.6244   ),
            5 => Element(  5, :B,  "boron",           10.811,   2.34,   8.30,  "[He] 2s^2 2p",		        0.8141,     2.428 ,    8.8328   ),
            6 => Element(  6, :C,  "carbon",          12.011,   2.26,  11.26,  "[He] 2s^2 2p^2",		    0.6513,     1.2432,   11.0407   ),
            7 => Element(  7, :N,  "nitrogen",        14.007,   1.25,  14.53,  "[He] 2s^2 2p^3",		    0.5427,     0.7193,   13.25     ),
            8 => Element(  8, :O,  "oxygen",          15.999,   1.43,  13.62,  "[He] 2s^2 2p^4",		    0.4652,     0.453 ,   15.4574   ),
            9 => Element(  9, :F,  "fluorine",        18.998,   1.70,  17.42,  "[He] 2s^2 2p^5",		    0.4071,     0.3034,   17.6634   ),
            10 => Element( 10, :Ne, "neon",            20.180,   0.90,  21.56,  "[Ne]",			            0.3618,     0.2131,   19.875    ),
            11 => Element( 11, :Na, "sodium",          22.990,   0.97,   5.14,  "[Ne] 3s",			        2.1649,    45.659 ,    3.3215   ),
            12 => Element( 12, :Mg, "magnesium",       24.305,   1.74,   7.65,  "[Ne] 3s^2", 		        1.6711,    21.0   ,    4.303    ),
            13 => Element( 13, :Al, "aluminium",       26.982,   2.70,   5.99,  "[Ne] 3s^2 3p",		        1.3607,    11.337 ,    5.2846   ),
            14 => Element( 14, :Si, "silicon",         28.086,   2.33,   8.15,  "[Ne] 3s^2 3p^2",		    1.1476,     6.8011,    6.2659   ),
            15 => Element( 15, :P,  "phosphorus",      30.974,   1.82,  10.49,  "[Ne] 3s^2 3p^3",		    0.9922,     4.3955,    7.2473   ),
            16 => Element( 16, :S,  "sulfur",          32.065,   2.07,  10.36,  "[Ne] 3s^2 3p^4",		    0.8738,     3.0023,    8.2293   ),
            17 => Element( 17, :Cl, "chlorine",        35.453,   3.21,  12.97,  "[Ne] 3s^2 3p^5",		    0.7807,     2.1412,    9.2107   ),
            18 => Element( 18, :Ar, "argon",           39.948,   1.78,  15.76,  "[Ar]",			            0.7056,     1.5808,    10.191   ),
            19 => Element( 19, :K,  "potassium",       39.098,   0.86,   4.34,  "[Ar] 4s",			        3.5598,   202.9969,     2.02    ),
            20 => Element( 20, :Ca, "calcium",         40.078,   1.55,   6.11,  "[Ar] 4s^2", 		        2.7479,    93.3717,     2.6162  ),
            21 => Element( 21, :Sc, "scandium",        44.956,   2.99,   6.56,  "[Ar] 3d 4s^2",		        2.6106,    80.0633,     2.7545  ),
            22 => Element( 22, :Ti, "titanium",        47.867,   4.54,   6.83,  "[Ar] 3d^2 4s^2",		    2.4861,    69.1712,     2.8924  ),
            23 => Element( 23, :V,  "vanadium",        50.942,   6.11,   6.75,  "[Ar] 3d^3 4s^2",		    2.3732,    60.1472,     3.03    ),
            24 => Element( 24, :Cr, "chromium",        51.996,   7.19,   6.77,  "[Ar] 3d^5 4s",		        2.2701,    52.6438,     3.1676  ),
            25 => Element( 25, :Mn, "manganese",       54.938,   7.43,   7.43,  "[Ar] 3d^5 4s^2",		    2.1754,    46.3265,     3.3055  ),
            26 => Element( 26, :Fe, "iron",            55.845,   7.87,   7.90,  "[Ar] 3d^6 4s^2",		    2.0885,    40.9939,     3.443   ),
            27 => Element( 27, :Co, "cobalt",          58.933,   8.90,   7.88,  "[Ar] 3d^7 4s^2",		    2.008 ,    36.4555,     3.5811  ),
            28 => Element( 28, :Ni, "nickel",          58.693,   8.90,   7.64,  "[Ar] 3d^8 4s^2",		    1.9337,    32.5372,     3.7187  ),
            29 => Element( 29, :Cu, "copper",          63.546,   8.96,   7.73,  "[Ar] 3d^10 4s",		        1.8648,    29.1769,     3.8561  ),
            30 => Element( 30, :Zn, "zinc",            65.390,   7.13,   9.39,  "[Ar] 3d^10 4s^2",		    1.8004,    26.2615,     3.994   ),
            31 => Element( 31, :Ga, "gallium",         69.723,   5.91,   6.00,  "[Ar] 3d^10 4s^2 4p",	    1.5663,    17.2917,     4.5909  ),
            32 => Element( 32, :Ge, "germanium",       72.640,   5.32,   7.90,  "[Ar] 3d^10 4s^2 4p^2",	    1.3862,    11.9864,     5.1874  ),
            33 => Element( 33, :As, "arsenic",         74.922,   5.72,   9.79,  "[Ar] 3d^10 4s^2 4p^3",	    1.2431,     8.6443,     5.7846  ),
            34 => Element( 34, :Se, "selenium",        78.960,   4.79,   9.75,  "[Ar] 3d^10 4s^2 4p^4",	    1.1269,     6.4397,     6.381   ),
            35 => Element( 35, :Br, "bromine",         79.904,   3.12,  11.81,  "[Ar] 3d^10 4s^2 4p^5",	    1.0305,     4.9244,     6.978   ),
            36 => Element( 36, :Kr, "krypton",         83.800,   3.75,  14.00,  "[Kr]",			            0.9493,     3.8497,     7.5748  ),
            37 => Element( 37, :Rb, "rubidium",        85.468,   1.63,   4.18,  "[Kr] 5s",			        4.8106,   500.9683,     1.4948  ),
            38 => Element( 38, :Sr, "strontium",       87.620,   2.54,   5.69,  "[Kr] 5s^2", 		        3.7135,   230.4462,     1.9364  ),
            39 => Element( 39, :Y,  "yttrium",         88.906,   4.47,   6.22,  "[Kr] 4d 5s^2",		        3.5278,   197.5715,     2.0383  ),
            40 => Element( 40, :Zr, "zirconium",       91.224,   6.51,   6.63,  "[Kr] 4d^2 5s^2",		    3.3598,   170.6683,     2.1402  ),
            41 => Element( 41, :Nb, "niobium",         92.906,   8.57,   6.76,  "[Kr] 4d^4 5s",		        3.2071,   148.4397,     2.2422  ),
            42 => Element( 42, :Mo, "molybdenum",      95.940,  10.22,   7.09,  "[Kr] 4d^5 5s",		        3.0677,   129.9126,     2.344   ),
            43 => Element( 43, :Tc, "technetium",      98.000,  11.50,   7.28,  "[Kr] 4d^5 5s^2",		    2.9398,   114.3315,     2.446   ),
            44 => Element( 44, :Ru, "ruthenium",      101.070,  12.37,   7.36,  "[Kr] 4d^7 5s",		        2.8222,   101.1523,     2.5479  ),
            45 => Element( 45, :Rh, "rhodium",        102.906,  12.41,   7.46,  "[Kr] 4d^8 5s",		        2.7137,    89.9286,     2.6498  ),
            46 => Element( 46, :Pd, "palladium",      106.420,  12.02,   8.34,  "[Kr] 4d^10",		        2.6132,    80.3028,     2.7517  ),
            47 => Element( 47, :Ag, "silver",         107.868,  10.50,   7.58,  "[Kr] 4d^10 5s",		        2.5199,    72.005 ,     2.8536  ),
            48 => Element( 48, :Cd, "cadmium",        112.411,   8.65,   8.99,  "[Kr] 4d^10 5s^2",		    2.433 ,    64.8095,     2.9555  ),
            49 => Element( 49, :In, "indium",         114.818,   7.31,   5.79,  "[Kr] 4d^10 5s^2 5p",	    2.1167,    42.6767,     3.3972  ),
            50 => Element( 50, :Sn, "tin",            118.710,   7.31,   7.34,  "[Kr] 4d^10 5s^2 5p^2",	    1.8732,    29.5777,     3.8388  ),
            51 => Element( 51, :Sb, "antimony",       121.760,   6.68,   8.61,  "[Kr] 4d^10 5s^2 5p^3",	    1.6799,    21.3335,     4.2805  ),
            52 => Element( 52, :Te, "tellurium",      127.600,   6.24,   9.01,  "[Kr] 4d^10 5s^2 5p^4",	    1.5228,    15.8906,     4.7221  ),
            53 => Element( 53, :I,  "iodine",         126.905,   4.93,  10.45,  "[Kr] 4d^10 5s^2 5p^5",	    1.3926,    12.1532,     5.0636  ),
            54 => Element( 54, :Xe, "xenon",          131.293,   5.90,  12.13,  "[Xe]",			            1.2828,     9.4993,     5.6056  ),
            55 => Element( 55, :Cs, "caesium",        132.906,   1.87,   3.89,  "[Xe] 6s",			        6.0615,  1002.1964,     1.1863  ),
            56 => Element( 56, :Ba, "barium",         137.327,   3.59,   5.21,  "[Xe] 6s^2", 		        4.6788,   460.9098,     1.5369  ),
            57 => Element( 57, :La, "lanthanum",      138.906,   6.15,   5.58,  "[Xe] 5d 6s^2",		        3.8102,   248.9177,     1.8873  ),
            58 => Element( 58, :Ce, "cerium",         140.116,   6.77,   5.54,  "[Xe] 4f 5d 6s^2",		    3.2133,   149.2883,     2.2378  ),
            59 => Element( 59, :Pr, "praseodymium",   140.908,   6.77,   5.47,  "[Xe] 4f^3 6s^2",		    2.778 ,    96.4738,     2.5885  ),
            60 => Element( 60, :Nd, "neodymium",      144.240,   7.01,   5.53,  "[Xe] 4f^4 6s^2",		    2.4468,    65.9186,     2.9389  ),
            61 => Element( 61, :Pm, "promethium",     145.000,   7.30,   5.58,  "[Xe] 4f^5 6s^2",		    2.1861,    47.0135,     3.2893  ),
            62 => Element( 62, :Sm, "samarium",       150.360,   7.52,   5.64,  "[Xe] 4f^6 6s^2",		    1.9756,    34.6984,     3.6398  ),
            63 => Element( 63, :Eu, "europium",       151.964,   5.24,   5.67,  "[Xe] 4f^7 6s^2",		    1.802 ,    26.3316,     3.9905  ),
            64 => Element( 64, :Gd, "gadolinium",     157.250,   7.90,   6.15,  "[Xe] 4f^7 5d 6s^2", 	    1.6565,    20.4544,     4.341   ),
            65 => Element( 65, :Tb, "terbium",        158.925,   8.23,   5.86,  "[Xe] 4f^9 6s^2",		    1.5328,    16.2057,     4.6913  ),
            66 => Element( 66, :Dy, "dysprosium",     162.500,   8.55,   5.94,  "[Xe] 4f^10 6s^2",		    1.4262,    13.0543,     5.0419  ),
            67 => Element( 67, :Ho, "holmium",        164.930,   8.80,   6.02,  "[Xe] 4f^11 6s^2",		    1.3335,    10.6707,     5.3924  ),
            68 => Element( 68, :Er, "erbium",         167.259,   9.07,   6.11,  "[Xe] 4f^12 6s^2",		    1.2521,     8.8334,     5.743   ),
            69 => Element( 69, :Tm, "thulium",        168.934,   9.32,   6.18,  "[Xe] 4f^13 6s^2",		    1.1801,     7.3955,     6.0934  ),
            70 => Element( 70, :Yb, "ytterbium",      173.040,   6.90,   6.25,  "[Xe] 4f^14 6s^2",		    1.1159,     6.253 ,     6.4439  ),
            71 => Element( 71, :Lu, "lutetium",       174.967,   9.84,   5.43,  "[Xe] 4f^14 5d 6s^2",	    1.0583,     5.3338,     6.7947  ),
            72 => Element( 72, :Hf, "hafnium",        178.490,  13.31,   6.83,  "[Xe] 4f^14 5d^2 6s^2",	    1.0079,     4.6075,     7.1344  ),
            73 => Element( 73, :Ta, "tantalum",       180.948,  16.65,   7.55,  "[Xe] 4f^14 5d^3 6s^2",	    0.9594,     3.9739,     7.4951  ),
            74 => Element( 74, :W,  "tungsten",       183.840,  19.35,   7.86,  "[Xe] 4f^14 5d^4 6s^2",	    0.9165,     3.4643,     7.8459  ),
            75 => Element( 75, :Re, "rhenium",        186.207,  21.04,   7.83,  "[Xe] 4f^14 5d^5 6s^2",	    0.8773,     3.0385,     8.1965  ),
            76 => Element( 76, :Os, "osmium",         190.230,  22.60,   8.44,  "[Xe] 4f^14 5d^6 6s^2",	    0.8413,     2.6796,     8.5472  ),
            77 => Element( 77, :Ir, "iridium",        192.217,  22.40,   8.97,  "[Xe] 4f^14 5d^7 6s^2",	    0.8182,     2.3756,     8.8973  ),
            78 => Element( 78, :Pt, "platinum",       195.078,  21.45,   8.96,  "[Xe] 4f^14 5d^9 6s",	    0.7776,     2.1175,     9.2474  ),
            79 => Element( 79, :Au, "gold",           196.967,  19.32,   9.23,  "[Xe] 4f^14 5d^10 6s",	    0.7492,     1.8924,     9.598   ),
            80 => Element( 80, :Hg, "mercury",        200.590,  13.55,  10.44,  "[Xe] 4f^14 5d^10 6s^2",	    3.0636,   129.646 ,     2.3456  ),
            81 => Element( 81, :Tl, "thallium",       204.383,  11.85,   6.11,  "[Xe] 4f^14 5d^10 6s^2 6p",  2.667 ,    85.3653,     2.6962  ),
            82 => Element( 82, :Pb, "lead",           207.200,  11.35,   7.42,  "[Xe] 4f^14 5d^10 6s^2 6p^2",2.3603,    59.1717,     3.0466  ),
            83 => Element( 83, :Bi, "bismuth",        208.980,   9.75,   7.29,  "[Xe] 4f^14 5d^10 6s^2 6p^3",2.1167,    42.6767,     3.3972  ),
            84 => Element( 84, :Po, "polonium",       209.000,   9.30,   8.42,  "[Xe] 4f^14 5d^10 6s^2 6p^4",1.9187,    31.7858,     3.7477  ),
            85 => Element( 85, :At, "astatine",       210.000,    NaN,   9.30,  "[Xe] 4f^14 5d^10 6s^2 6p^5",1.7546,    24.3079,     4.0983  ),
            86 => Element( 86, :Rn, "radon",          222.000,   9.73,  10.75,  "[Rn]",			            1.6164,    19.0046,     4.4487  ),
            87 => Element( 87, :Fr, "francium",       223.000,    NaN,   4.07,  "[Rn] 7s",			        7.2404,  1705.0   ,     0.9913  ),
            88 => Element( 88, :Ra, "radium",         226.000,   5.50,   5.28,  "[Rn] 7s^2", 		        5.5887,   785.4977,     1.2867  ),
            89 => Element( 89, :Ac, "actinium",       227.000,  10.07,   5.17,  "[Rn] 6d 7s^2",		        5.3091,   673.4033,     1.3544  ),
            90 => Element( 90, :Th, "thorium",        232.038,  11.72,   6.31,  "[Rn] 6d^2 7s^2",		    5.0569,   581.6815,     1.4222  ),
            91 => Element( 91, :Pa, "protactinium",   231.036,  15.40,   5.89,  "[Rn] 5f^2 6d 7s^2", 	    3.7042,   228.7156,     1.9413  ),
            92 => Element( 92, :U,  "uranium",        238.029,  18.95,   6.19,  "[Rn] 5f^3 6d 7s^2", 	    3.2177,   149.9164,     2.2348  ),
            93 => Element( 93, :Np, "neptunium",      237.000,  20.20,   6.27,  "[Rn] 5f^4 6d 7s^2", 	    2.8443,   103.5473,     2.5281  ),
            94 => Element( 94, :Pu, "plutonium",      244.000,  19.84,   6.03,  "[Rn] 5f^6 7s^2",		    2.3596,    59.1266,     3.0473  ),
            95 => Element( 95, :Am, "americium",      243.000,  13.67,   5.97,  "[Rn] 5f^7 7s^2",		    2.1525,    44.8789,     3.3407  ),
            96 => Element( 96, :Cm, "curium",         247.000,  13.50,   5.99,  "[Rn] 5f^7 6d 7s^2", 	    2.1097,    42.2547,     3.4084  ),
            97 => Element( 97, :Bk, "berkelium",      247.000,  14.78,   6.20,  "[Rn] 5f^9 7s^2",		    1.8308,    27.5918,     3.9277  ),
            98 => Element( 98, :Cf, "californium",    251.000,  15.10,   6.28,  "[Rn] 5f^10 7s^2",		    1.7035,    22.2453,     4.2212  ),
            99 => Element( 99, :Es, "einsteinium",    252.000,    NaN,   6.42,  "[Rn] 5f^11 7s^2",		    1.5928,    18.1843,     4.5146  ),
            100 => Element(100, :Fm, "fermium",        257.000,    NaN,   6.50,  "[Rn] 5f^12 7s^2",		    1.4956,    15.0542,     4.808   ),
            101 => Element(101, :Md, "mendelevium",    258.000,    NaN,   6.58,  "[Rn] 5f^13 7s^2",		    1.4096,    12.6038,     5.1013  ),
            102 => Element(102, :No, "nobelium",       259.000,    NaN,   6.65,  "[Rn] 5f^14 7s^2",		    1.3329,    10.6563,     5.3949  ),
            103 => Element(103, :Lr, "lawrencium",     262.000,    NaN,   4.90,  "[Rn] 5f^14 7s^2 7p",	    1.3164,    10.2631,     5.4629  ),
            104 => Element(104, :Rf, "rutherfordium",  261.000,    NaN,    NaN,  "[Rn] 5f^14 6d^2 7s^2",	       NaN,        NaN,        NaN  ),
            105 => Element(105, :Db, "dubnium",        262.000,    NaN,    NaN,  "[Rn] 5f^14 6d^3 7s^2",	       NaN,        NaN,        NaN  ),
            106 => Element(106, :Sg, "seaborgium",     266.000,    NaN,    NaN,  "[Rn] 5f^14 6d^4 7s^2",	       NaN,        NaN,        NaN  ),
            107 => Element(107, :Bh, "bohrium",        264.000,    NaN,    NaN,  "[Rn] 5f^14 6d^5 7s^2",	       NaN,        NaN,        NaN  ),
            108 => Element(108, :Hs, "hassium",        277.000,    NaN,    NaN,  "[Rn] 5f^14 6d^6 7s^2",	       NaN,        NaN,        NaN  ),
            109 => Element(109, :Mt, "meitnerium",     268.000,    NaN,    NaN,  "[Rn] 5f^14 6d^7 7s^2",	       NaN,        NaN,        NaN  ),
            110 => Element(110, :Ds, "darmstadtium",       NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^8 7s^2",	       NaN,        NaN,        NaN  ),
            111 => Element(111, :Rg, "roentgenium",        NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^9 7s^2",	       NaN,        NaN,        NaN  ),
            112 => Element(112, :Cn, "copernicium",        NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^10 7s^2",	       NaN,        NaN,        NaN  ),
            113 => Element(113, :Nh, "nihonium",           NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^8 7s^2 7p",	   NaN,        NaN,        NaN  ),
            114 => Element(114, :Fl, "flerovium",          NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^2 7s^2 7p^2",    NaN,        NaN,        NaN  ),
            115 => Element(115, :Mc, "moscovium",          NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^3 7s^2 7p^3",    NaN,        NaN,        NaN  ),
            116 => Element(116, :Lv, "livermorium",        NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^4 7s^2 7p^4",    NaN,        NaN,        NaN  ),
            117 => Element(117, :Ts, "tennessine",         NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^5 7s^2 7p^5",    NaN,        NaN,        NaN  ),
            118 => Element(118, :Og, "oganesson",          NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^6 7s^2 7p^6",    NaN,        NaN,        NaN  )	 )
                                                        
                                                            
"""
`PeriodicTable.isotopes::Dict{Int64,Array{PeriodicTable.Isotope,1}}`
    ... provides isotope data for an element with given nuclear charge Z::Int64 from PeriodicTable.isotopes[Z];
        cf. PeriodicTable.Isotope
"""
const isotopes = Dict{Int64,Array{PeriodicTable.Isotope,1}}(
            1 => Isotope[],
            2 => Isotope[ Isotope( 2, 2.015894,      AngularJ64(0), NaN, NaN),    Isotope( 2, 3.01602932265, AngularJ64(1//2), NaN, NaN), 
                            Isotope( 4, 4.00260325413, AngularJ64(0), NaN, NaN)],
            3 => Isotope[],
            4 => Isotope[],
            5 => Isotope[],
            6 => Isotope[],
            7 => Isotope[],
            8 => Isotope[],
            9 => Isotope[],
            10 => Isotope[],
            11 => Isotope[],
            12 => Isotope[],
            13 => Isotope[],
            14 => Isotope[],
            15 => Isotope[],
            16 => Isotope[],
            17 => Isotope[],
            18 => Isotope[],
            19 => Isotope[],
            20 => Isotope[],
            21 => Isotope[],
            22 => Isotope[],
            23 => Isotope[],
            24 => Isotope[],
            25 => Isotope[],
            26 => Isotope[],
            27 => Isotope[],
            28 => Isotope[],
            29 => Isotope[],
            30 => Isotope[],
            31 => Isotope[],
            32 => Isotope[],
            33 => Isotope[],
            34 => Isotope[],
            35 => Isotope[],
            36 => Isotope[],
            37 => Isotope[],
            38 => Isotope[],
            39 => Isotope[],
            40 => Isotope[],
            41 => Isotope[],
            42 => Isotope[],
            43 => Isotope[],
            44 => Isotope[],
            45 => Isotope[],
            46 => Isotope[],
            47 => Isotope[],
            48 => Isotope[],
            49 => Isotope[],
            50 => Isotope[], 
            51 => Isotope[],
            52 => Isotope[],
            53 => Isotope[],
            54 => Isotope[],
            55 => Isotope[],
            56 => Isotope[],
            57 => Isotope[],
            58 => Isotope[],
            59 => Isotope[],
            60 => Isotope[],
            61 => Isotope[],
            62 => Isotope[],
            63 => Isotope[],
            64 => Isotope[],
            65 => Isotope[],
            66 => Isotope[],
            67 => Isotope[],
            68 => Isotope[],
            69 => Isotope[],
            70 => Isotope[],
            71 => Isotope[],
            72 => Isotope[],
            73 => Isotope[],
            74 => Isotope[],
            75 => Isotope[],
            76 => Isotope[],
            77 => Isotope[],
            78 => Isotope[],
            79 => Isotope[],
            80 => Isotope[],
            81 => Isotope[],
            82 => Isotope[],
            83 => Isotope[],
            84 => Isotope[],
            85 => Isotope[],
            86 => Isotope[],
            87 => Isotope[],
            88 => Isotope[],
            89 => Isotope[],
            90 => Isotope[],
            91 => Isotope[],
            92 => Isotope[],
            93 => Isotope[],
            94 => Isotope[],
            95 => Isotope[],
            96 => Isotope[],
            97 => Isotope[],
            98 => Isotope[],
            99 => Isotope[],
            100 => Isotope[],
            101 => Isotope[],
            102 => Isotope[],
            103 => Isotope[],
            104 => Isotope[],
            105 => Isotope[],
            106 => Isotope[],
            107 => Isotope[],
            108 => Isotope[],
            109 => Isotope[],
            110 => Isotope[],
            111 => Isotope[],
            112 => Isotope[],
            113 => Isotope[],
            114 => Isotope[],
            115 => Isotope[],
            116 => Isotope[],
            117 => Isotope[],
            118 => Isotope[] )


"""
`PeriodicTable.bindingEnergies_Larkins1977(Z::Int64)`  
    ... to return the `stored values' of Larkins (1977) which are based on Sevier and others for an element with nuclear charge Z.
"""
function bindingEnergies_Larkins1977(Z::Int64)
    # This procedure 'stores' the binding energies of the inner-shell electrons 1s, ..., 3d for all elements from Z = 1, ..., 45;
    # all binding energies are given in eV
    # 
    #                        1s_1/2,   2s_1/2,   2p_1/2,   2p_3/2,   3s_1/2,   3p_1/2,   3p_3/2,   3d_3/2,   3d_5/2,   4s_1/2,   4p_1/2,   4p_3/2 
    if     Z ==  1   wa = [    13.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  2   wa = [    24.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  3   wa = [    54.8,      5.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  4   wa = [   112.1,      8.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] 
    elseif Z ==  5   wa = [   188.0,     12.6,      4.7,      4.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  6   wa = [   283.8,     18.0,      6.4,      6.4,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  7   wa = [   401.6,     24.4,      9.2,      9.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  8   wa = [   532.0,     28.5,      7.1,      7.1,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  9   wa = [   685.4,     34.0,      8.6,      8.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 10   wa = [   870.1,     48.5,     21.7,     21.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    #
    elseif Z == 11   wa = [  1072.1,     63.3,     31.1,     31.1,      1.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] 
    elseif Z == 12   wa = [  1305.0,     89.4,     51.4,     51.4,      2.1,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] 
    elseif Z == 13   wa = [  1559.6,    117.7,     73.2,     72.7,      7.0,      5.5,      5.5,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 14   wa = [  1838.9,    148.7,     99.5,     98.9,     10.6,      6.0,      6.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 15   wa = [  2145.5,    189.3,    136.2,    135.3,     16.2,      9.9,      9.9,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 16   wa = [  2472.0,    229.2,    165.4,    164.2,     15.8,      8.0,      8.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 17   wa = [  2822.4,    270.2,    201.6,    200.0,     17.5,      6.8,      6.8,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 18   wa = [  3206.0,    326.3,    250.7,    248.6,     29.2,     15.9,     15.8,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 19   wa = [  3607.4,    377.1,    296.3,    293.6,     33.9,     17.8,     17.8,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 20   wa = [  4038.1,    437.8,    350.0,    346.4,     43.7,     25.4,     25.4,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    #
    elseif Z == 21   wa = [  4492.8,    500.4,    406.7,    402.2,     53.8,     32.3,     32.3,      6.6,      6.6,     -1.0,     -1.0,     -1.0]        
    elseif Z == 22   wa = [  4966.4,    563.7,    461.5,    455.5,     60.3,     34.6,     34.6,      3.7,      3.7,     -1.0,     -1.0,     -1.0]  
    elseif Z == 23   wa = [  5465.1,    628.2,    520.5,    512.9,     66.5,     37.8,     37.8,      2.2,      2.2,     -1.0,     -1.0,     -1.0]
    elseif Z == 24   wa = [  5989.2,    694.6,    583.7,    574.5,     74.1,     42.5,     42.5,      2.3,      2.3,     -1.0,     -1.0,     -1.0]
    elseif Z == 25   wa = [  6539.0,    769.0,    651.4,    640.3,     83.9,     48.6,     48.6,      3.3,      3.3,     -1.0,     -1.0,     -1.0]
    elseif Z == 26   wa = [  7112.0,    846.1,    721.1,    708.1,     92.9,     54.0,     54.0,      3.6,      3.6,     -1.0,     -1.0,     -1.0]
    elseif Z == 27   wa = [  7708.9,    925.6,    793.6,    778.6,    100.7,     59.5,     59.5,      2.9,      2.9,     -1.0,     -1.0,     -1.0]
    elseif Z == 28   wa = [  8332.8,   1008.1,    871.9,    854.7,    111.8,     68.1,     68.1,      3.6,      3.6,     -1.0,     -1.0,     -1.0]
    elseif Z == 29   wa = [  8978.9,   1096.1,    951.0,    931.1,    119.8,     73.6,     73.6,      1.6,      1.6,     -1.0,     -1.0,     -1.0]
    elseif Z == 30   wa = [  9658.6,   1193.6,   1042.8,   1019.7,    135.9,     86.6,     86.6,      8.1,      8.1,     -1.0,     -1.0,     -1.0]
    #
    elseif Z == 31   wa = [ 10367.1,   1297.7,   1142.3,   1115.4,    158.1,    106.8,    102.9,     17.4,     17.4,      1.5,      0.8,      0.8]
    elseif Z == 32   wa = [ 11103.1,   1414.3,   1247.8,   1216.7,    180.0,    127.9,    120.8,     28.7,     28.7,      5.0,      2.3,      2.3]
    elseif Z == 33   wa = [ 11866.7,   1526.5,   1358.6,   1323.1,    203.5,    146.4,    140.5,     41.2,     41.2,      8.5,      2.5,      2.5]
    elseif Z == 34   wa = [ 12657.8,   1653.9,   1476.2,   1435.8,    231.5,    168.2,    161.9,     56.7,     56.7,     12.0,      5.6,      5.6]
    elseif Z == 35   wa = [ 13473.7,   1782.0,   1596.0,   1549.9,    256.5,    189.3,    181.5,     70.1,     69.0,     27.3,      5.2,      4.6]
    elseif Z == 36   wa = [ 14325.6,   1921.0,   1727.2,   1674.9,    292.1,    221.8,    214.5,     95.0,     93.8,     27.5,     14.7,     14.0]
    elseif Z == 37   wa = [ 15199.7,   2065.1,   1863.9,   1804.4,    322.1,    247.4,    238.5,    111.8,    110.3,     29.3,     14.8,     14.0]    
    elseif Z == 38   wa = [ 16104.6,   2216.3,   2006.8,   1939.6,    357.5,    279.8,    269.1,    135.0,    133.1,     37.7,     19.9,     19.9]    
    #
    ##  [keV]
    ##  39 Y   17.0384  2.3725  2.1555  2.0800  0.3936  0.3124  0.3003  0.1596  0.1574  0.0454  0.0256  0.0256  0.0024  0.0024
    ##  40 Zr  17.9976  2.5316  2.3067  2.2223  0.4303  0.3442  0.3305  0.1824  0.1800  0.0513  0.0287  0.0287  0.0030  0.0030
    ##  41 Nb  18.9856  2.6977  2.4647  2.3705  0.4684  0.3784  0.3630  0.2074  0.2046  0.0581  0.0339  0.0339  0.0032  0.0032
    ##  42 Mo  19.9995  2.8655  2.6251  2.5202  0.5046  0.4097  0.3923  0.2303  0.2270  0.0618  0.0348  0.0348  0.0018  0.0018
    ##  43 Tc  21.0440  3.0425  2.7932  2.6769  0.5440  0.4449  0.4250  0.2564  0.2529  0.0680  0.0389  0.0389  0.0020  0.0020
    ##  44 Ru  22.1172  3.2240  2.9669  2.8379  0.5850  0.4828  0.4606  0.2836  0.2794  0.0749  0.0431  0.0431  0.0020  0.0020
    ##  45 Rh  23.2199  3.4119  3.1461  3.0038  0.6271  0.5210  0.4962  0.3117  0.3070  0.0810  0.0479  0.0479  0.0025  0.0025        
    else  error("Data not available for Z = $Z")
    end
    return( wa )
end



"""
`PeriodicTable.bindingEnergies_Williams2000(Z::Int64)`  ... former  JAC.store_Williams2000(Z::Int64)
    ... to return the `stored values' of Williams et al. (2000), ... for the element with nuclear charge Z.
"""
function bindingEnergies_Williams2000(Z::Int64)
    # This procedure 'stores' the binding energies of the inner-shell electrons 1s, ..., 4p for all elements from Z = 1, ..., 36
    # (not 54, despite what an earlier version of this comment claimed -- verified against the actual branches below);
    # all binding energies are given in eV. For the fuller Z = 1, ..., 92 coverage, use XrayDataBooklet instead, which
    # follows the same underlying LBNL X-Ray Data Booklet Table 1-1 and is numerically identical wherever both are tabulated.
    # 
    #                        1s_1/2,   2s_1/2,   2p_1/2,   2p_3/2,   3s_1/2,   3p_1/2,   3p_3/2,   3d_3/2,   3d_5/2,   4s_1/2,   4p_1/2,   4p_3/2 
    if     Z ==  1   wa = [    13.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  2   wa = [    24.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  3   wa = [    54.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  4   wa = [   111.5,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] 
    elseif Z ==  5   wa = [   188.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  6 # wa = [   284.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
                     wa = [   284.2,     19.4,      7.0,      7.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] #ZS
    elseif Z ==  7   wa = [   409.9,     37.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  8   wa = [   543.1,     41.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  9 # wa = [   696.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
                     wa = [   696.7,     31.0,      9.0,      9.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] #ZS
    elseif Z == 10   wa = [   870.2,     48.5,     21.7,     21.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    #
    elseif Z == 11   wa = [  1070.8,     63.5,    30.65,    30.81,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] 
    elseif Z == 12   wa = [  1303.0,     88.7,    49.78,    49.50,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] 
    elseif Z == 13   wa = [  1559.6,    117.8,    72.95,    72.55,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 14   wa = [  1839.0,    149.7,    99.82,    99.42,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 15   wa = [  2145.5,    189.0,   136.0,    135.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 16   wa = [  2472.0,    230.9,   163.6,    162.5,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 17   wa = [  2822.4,    270.0,   202.0,    200.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 18   wa = [  3205.9,    326.3,   250.6,    248.4,      29.3,     15.9,     15.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 19   wa = [  3608.4,    378.6,   297.3,    294.6,      34.8,     18.3,     18.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 20   wa = [  4038.5,    438.4,   349.7,    346.2,      44.3,     25.4,     25.4,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    #
    elseif Z == 21   wa = [  4492.0,    498.0,   403.6,    398.7,      51.1,     28.3,     28.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 22   wa = [  4966.0,    560.9,   460.2,    453.8,      58.7,     32.6,     32.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 23   wa = [  5465.0,    626.7,   519.8,    512.1,      66.3,     37.2,     37.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 24   wa = [  5989.0,    696.0,   583.8,    574.1,      74.1,     42.2,     42.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 25   wa = [  6539.0,    769.1,   649.9,    638.7,      82.3,     47.2,     47.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 26   wa = [  7112.0,    844.6,   719.9,    706.8,      91.3,     52.7,     52.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 27   wa = [  7709.0,    925.1,   793.2,    778.1,     101.0,     58.9,     59.9,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 28   wa = [  8333.0,   1008.6,   870.0,    852.7,     110.8,     68.0,     66.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 29   wa = [  8979.0,   1096.7,   952.3,    932.7,     122.5,     77.3,     75.1,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 30   wa = [  9659.0,   1196.2,  1044.9,   1021.8,     139.8,     91.4,     88.6,     10.2,     10.1,     -1.0,     -1.0,     -1.0]
    #
    elseif Z == 31   wa = [ 10367.0,   1299.0,  1143.2,   1116.4,     159.5,    103.5,    100.0,     18.7,     18.7,     -1.0,     -1.0,     -1.0]
    elseif Z == 32   wa = [ 11103.0,   1414.6,  1248.1,   1217.0,     180.1,    124.9,    120.8,     29.8,     29.2,     -1.0,     -1.0,     -1.0]
    elseif Z == 33   wa = [ 11867.0,   1527.0,  1359.1,   1323.6,     204.7,    146.2,    141.2,     41.7,     41.7,     -1.0,     -1.0,     -1.0]
    elseif Z == 34   wa = [ 12658.0,   1652.0,  1474.3,   1433.9,     229.6,    166.5,    160.7,     55.5,     54.6,     -1.0,     -1.0,     -1.0]
    elseif Z == 35   wa = [ 13474.0,   1782.0,  1596.0,   1550.0,     257.0,    189.0,    182.0,     70.0,     69.0,     -1.0,     -1.0,     -1.0]
    elseif Z == 36   wa = [ 14326.0,   1921.0,  1730.9,   1678.4,     292.8,    222.2,    214.4,     95.0,     93.8,     27.5,     14.1,     14.1]
    #
    else  error("Data not available for Z = $Z")
    end
    return( wa )
end


"""
`PeriodicTable.bindingEnergies_XrayDataBooklet(Z::Int64)` 
    ... to return the `stored values' of Bearden et al. (1967), Cardona et al. (1978) and Fuggle et al. (1980) 
        for the element with nuclear charge Z; all binding energies are given in eV.
"""
function bindingEnergies_XrayDataBooklet(Z::Int64)
    # 
    #                        1s_1/2,   2s_1/2,   2p_1/2,   2p_3/2,   3s_1/2,   3p_1/2,   3p_3/2,   3d_3/2,   3d_5/2,   4s_1/2,   4p_1/2,   4p_3/2, 
    #                              4d_3/2,   4d_5/2,   4f_5/2,  4f_7/2,   5s_1/2,   5p_1/2,   5p_3/2,   5d_3/2,   5d_5/2,   6s_1/2,   6p_1/2,    6p_3/2
    if     Z ==  1   wa = [    13.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  2   wa = [    24.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  3   wa = [    54.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  4   wa = [   111.5,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0, 
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  5   wa = [   188.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  6   wa = [   284.2,     19.4,      7.0,      7.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ] #ZS: propagated from bindingEnergies_Williams2000's correction
    elseif Z ==  7   wa = [   409.9,     37.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  8   wa = [   543.1,     41.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  9   wa = [   696.7,     31.0,      9.0,      9.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0, 
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 10   wa = [   870.2,     48.5,     21.7,     21.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    #
    elseif Z == 11   wa = [  1070.8,     63.5,    30.65,    30.81,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0, 
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 12   wa = [  1303.0,     88.7,    49.78,    49.50,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0, 
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 13   wa = [  1559.6,    117.8,    72.95,    72.55,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 14   wa = [  1839.0,    149.7,    99.82,    99.42,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 15   wa = [  2145.5,    189.0,   136.0,    135.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 16   wa = [  2472.0,    230.9,   163.6,    162.5,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 17   wa = [  2822.4,    270.0,   202.0,    200.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 18   wa = [  3205.9,    326.3,   250.6,    248.4,      29.3,     15.9,     15.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 19   wa = [  3608.4,    378.6,   297.3,    294.6,      34.8,     18.3,     18.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 20   wa = [  4038.5,    438.4,   349.7,    346.2,      44.3,     25.4,     25.4,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    #
    elseif Z == 21   wa = [  4492.0,    498.0,   403.6,    398.7,      51.1,     28.3,     28.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 22   wa = [  4966.0,    560.9,   460.2,    453.8,      58.7,     32.6,     32.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 23   wa = [  5465.0,    626.7,   519.8,    512.1,      66.3,     37.2,     37.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 24   wa = [  5989.0,    696.0,   583.8,    574.1,      74.1,     42.2,     42.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 25   wa = [  6539.0,    769.1,   649.9,    638.7,      82.3,     47.2,     47.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 26   wa = [  7112.0,    844.6,   719.9,    706.8,      91.3,     52.7,     52.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 27   wa = [  7709.0,    925.1,   793.2,    778.1,     101.0,     58.9,     59.9,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 28   wa = [  8333.0,   1008.6,   870.0,    852.7,     110.8,     68.0,     66.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 29   wa = [  8979.0,   1096.7,   952.3,    932.7,     122.5,     77.3,     75.1,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 30   wa = [  9659.0,   1196.2,  1044.9,   1021.8,     139.8,     91.4,     88.6,     10.2,     10.1,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    #
    elseif Z == 31   wa = [ 10367.0,   1299.0,  1143.2,   1116.4,     159.5,    103.5,    100.0,     18.7,     18.7,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 32   wa = [ 11103.0,   1414.6,  1248.1,   1217.0,     180.1,    124.9,    120.8,     29.8,     29.2,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 33   wa = [ 11867.0,   1527.0,  1359.1,   1323.6,     204.7,    146.2,    141.2,     41.7,     41.7,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 34   wa = [ 12658.0,   1652.0,  1474.3,   1433.9,     229.6,    166.5,    160.7,     55.5,     54.6,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 35   wa = [ 13474.0,   1782.0,  1596.0,   1550.0,     257.0,    189.0,    182.0,     70.0,     69.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 36   wa = [ 14326.0,   1921.0,  1730.9,   1678.4,     292.8,    222.2,    214.4,     95.0,     93.8,     27.5,     14.1,     14.1,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    #
    elseif Z == 37   wa = [ 15200.0,   2065.0,  1864.0,   1804.0,     326.7,    248.7,    239.1,    113.0,    112.0,     30.5,     16.3,     15.3,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 38   wa = [ 16105.0,   2216.0,  2007.0,   1940.0,     358.7,    280.3,    270.0,    136.0,    134.2,     38.9,     21.3,     20.1,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 39   wa = [ 17038.0,   2373.0,  2156.0,   2080.0,     392.0,    310.6,    298.8,    157.7,    155.8,     43.8,     24.4,     23.1,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 40   wa = [ 17998.0,   2532.0,  2307.0,   2223.0,     430.3,    343.5,    329.8,    181.1,    178.8,     50.6,     28.5,     27.1,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 41   wa = [ 18986.0,   2698.0,  2465.0,   2371.0,     466.6,    376.1,    360.6,    205.0,    202.3,     56.4,     32.6,     30.8,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 42   wa = [ 20000.0,   2866.0,  2625.0,   2520.0,     506.3,    411.6,    394.0,    231.1,    227.9,     63.2,     37.6,     35.5,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 43   wa = [ 21044.0,   3043.0,  2793.0,   2677.0,     544.0,    447.6,    417.7,    257.6,    253.9,     69.5,     42.3,     39.9,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 44   wa = [ 22117.0,   3224.0,  2967.0,   2838.0,     586.1,    483.5,    461.4,    284.2,    280.0,     75.0,     46.3,     43.2,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 45   wa = [ 23220.0,   3412.0,  3146.0,   3004.0,     628.1,    521.3,    496.5,    311.9,    307.2,     81.4,     50.5,     47.3,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 46   wa = [ 24350.0,   3604.0,  3330.0,   3173.0,     671.6,    599.9,    532.3,    340.5,    335.2,     87.1,     55.7,     50.9,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 47   wa = [ 25514.0,   3806.0,  3524.0,   3351.0,     719.0,    603.8,    573.0,    374.0,    368.3,     97.0,     63.7,     58.3,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]  
    #
    elseif Z == 48   wa = [26711.0,    4018.0,  3727.0,   3538.0,     772.0,    652.6,    618.4,    411.9,    405.2,    109.8,     63.9,     63.9,
                                    11.7,     10.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 49   wa = [27940.0,    4238.0,  3938.0,   3730.0,     827.2,    703.2,    665.3,    451.4,    443.9,    122.9,     73.5,     73.5,
                                    17.7,     16.9,    -1.0,     -1.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 50   wa = [29200.0,    4465.0,  4156.0,   3929.0,     884.7,    756.5,    714.6,    493.2,    484.9,    137.1,     83.6,     83.6,
                                    24.9,     23.9,    -1.0,     -1.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 51   wa = [30491.0,    4698.0,  4380.0,   4132.0,     946.0,    812.7,    766.4,    537.5,    528.2,    153.2,     95.6,     95.6,
                                    33.3,     32.1,    -1.0,     -1.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 52   wa = [31814.0,    4939.0,  4612.0,   4341.0,    1006.0,    870.8,    820.0,    583.4,    573.0,    169.4,    103.3,    103.3,
                                    41.9,     40.4,    -1.0,     -1.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 53   wa = [33169.0,    5188.0,  4852.0,   4557.0,    1072.0,    931.0,    875.0,    630.8,    619.3,    186.0,    123.0,    123.0,
                                    50.6,     48.9,    -1.0,     -1.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 54   wa = [34561.0,    5453.0,  5107.0,   4786.0,    1148.7,   1002.1,    940.6,    689.0,    676.4,    213.2,    146.7,    145.5,
                                    69.5,     67.5,    -1.0,     -1.0,      23.3,     13.4,     12.1,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 55   wa = [35985.0,    5714.0,  5359.0,   5012.0,    1211.0,   1071.0,   1003.0,    740.5,    726.6,    232.3,    172.4,    161.3,
                                    79.8,     77.5,    -1.0,     -1.0,      22.7,     14.2,     12.1,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 56   wa = [37441.0,    5989.0,  5624.0,   5247.0,    1293.0,   1137.0,   1063.0,    795.7,    780.5,    253.5,    192.0,    178.6,
                                    92.6,     89.9,    -1.0,     -1.0,      30.3,     17.0,     14.8,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0] 
    elseif Z == 57   wa = [38925.0,    6266.0,  5891.0,   5483.0,    1362.0,   1209.0,   1128.0,    853.0,    836.0,    274.7,    205.8,    196.0,
                                    105.3,    102.5,    -1.0,     -1.0,      34.3,     19.3,     16.8,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 58   wa = [40443.0,    6549.0,  6164.0,   5723.0,    1436.0,   1274.0,   1187.0,    902.4,    883.8,    291.0,    223.2,    206.5,
                                    109.0,     -1.0,     0.1,      0.1,      37.8,     19.8,     17.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 59   wa = [41991.0,    6835.0,  6440.0,   5964.0,    1511.0,   1337.0,   1242.0,    948.3,    928.8,    304.5,    236.3,    217.6,
                                    115.1,    115.1,     2.0,      2.0,      37.4,     22.3,     22.3,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 60   wa = [43569.0,    7126.0,  6722.0,   6208.0,    1575.0,   1403.0,   1297.0,   1003.3,    980.4,    319.2,    243.3,    224.6,
                                    120.5,    120.5,     1.5,      1.5,      37.5,     21.1,     21.1,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 61   wa = [45184.0,    7428.0,  7013.0,   6459.0,      -1.0,   1471.0,   1357.0,   1052.0,   1027.0,     -1.0,    242.0,    242.0,
                                    120.0,    120.0,    -1.0,     -1.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 62   wa = [46834.0,    7737.0,  7312.0,   6716.0,    1723.0,   1541.0,   1420.0,   1110.9,   1083.4,    347.2,    265.6,    247.4,
                                    129.0,    129.0,     5.2,      5.2,      37.4,     21.3,     21.3,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 63   wa = [48519.0,    8052.0,  7617.0,   6977.0,    1800.0,   1614.0,   1481.0,   1158.6,   1127.5,    360.0,    284.0,    257.0,
                                    133.0,    127.7,     0.0,      0.0,      32.0,     22.0,     22.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 64   wa = [50239.0,    8376.0,  7930.0,   7243.0,    1881.0,   1688.0,   1544.0,   1221.9,   1189.6,    378.6,    286.0,    271.0,
                                    -1.0,    142.6,     8.6,      8.6,      36.0,     28.0,     21.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 65   wa = [51996.0,    8708.0,  8252.0,   7514.0,    1968.0,   1768.0,   1611.0,   1276.9,   1241.1,    396.0,    322.4,    284.1,
                                    150.5,    150.5,     7.7,      2.4,      45.6,     28.7,     22.6,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 66   wa = [53789.0,    9046.0,  8581.0,   7790.0,    2047.0,   1842.0,   1676.0,   1333.0,   1292.6,    414.2,    333.5,    293.2,
                                    153.6,    153.6,     8.0,      4.3,      49.9,     26.3,     26.3,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 67   wa = [55618.0,    9394.0,  8918.0,   8071.0,    2128.0,   1923.0,   1741.0,   1392.0,   1351.0,    432.4,    343.5,    308.2,
                                    160.0,    160.0,     8.6,      5.2,      49.3,     30.8,     24.1,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 68   wa = [57486.0,    9751.0,  9264.0,   8358.0,    2207.0,   2006.0,   1812.0,   1453.0,   1409.0,    449.8,    366.2,    320.2,
                                    167.6,    167.6,    -1.0,      4.7,      50.6,     31.4,     24.7,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 69   wa = [59390.0,   10116.0,  9617.0,   8648.0,    2307.0,   2090.0,   1885.0,   1515.0,   1468.0,    470.9,    385.9,    332.6,
                                    175.5,    175.5,    -1.0,      4.6,      54.7,     31.8,     25.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 70   wa = [61332.0,   10486.0,  9978.0,   8944.0,    2398.0,   2173.0,   1950.0,   1576.0,   1528.0,    480.5,    388.7,    339.7,
                                    191.2,    182.4,     2.5,      1.3,      52.0,     30.3,     24.1,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    #
    elseif Z == 71   wa = [63314.0,   10870.0, 10349.0,   9244.0,    2491.0,   2264.0,   2024.0,   1639.0,   1589.0,    506.8,    412.4,    359.2,
                                    206.1,    196.3,    89.0,      7.5,      57.3,     33.6,     26.7,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 72   wa = [65351.0,   11271.0, 10739.0,   9561.0,    2601.0,   2365.0,   2108.0,   1716.0,   1662.0,    538.0,    438.2,    380.7,
                                    220.0,    211.5,    15.9,     14.2,      64.2,     38.0,     29.9,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 73   wa = [67416.0,   11682.0, 11136.0,   9881.0,    2708.0,   2469.0,   2194.0,   1793.0,   1735.0,    563.4,    463.4,    400.9,
                                    237.9,    226.4,    23.5,     21.6,      69.7,     42.2,     32.7,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]		
    elseif Z == 74   wa = [69525.0,   12100.0, 11544.0,  10207.0,    2820.0,   2575.0,   2281.0,   1872.0,   1809.0,    594.1,    490.4,    423.6,
                                    255.9,    243.5,    33.6,     31.4,      75.6,     45.3,     36.8,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 75   wa = [71676.0,   12527.0, 11959.0,  10535.0,    2932.0,   2682.0,   2367.0,   1949.0,   1883.0,    625.4,    518.7,    446.8,
                                    273.9,    260.5,    42.9,     40.5,      83.0,     45.6,     34.6,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]	
    elseif Z == 76   wa = [73871.0,   12968.0, 12385.0,  10871.0,    3049.0,   2792.0,   2457.0,   2031.0,   1960.0,    658.2,    549,1,    470.7,
                                    293.1,    278.5,    53.4,     50.7,      84.0,     58.0,     44.5,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]	
    elseif Z == 77   wa = [76111.0,   13419.0, 12824.0,  11215.0,    3174.0,   2909.0,   2551.0,   2116.0,   2040.0,    691.1,    577.8,    495.8,
                                    311.9,    296.3,    63.8,     60.8,      95.2,     63.0,     48.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]	
    elseif Z == 78   wa = [78395.0,   13880.0, 13273.0,  11564.0,    3296.0,   3027.0,   2645.0,   2202.0,   2122.0,    725.4,    609.1,    519.4,
                                    331.6,    314.6,    74.5,     71.2,     101.7,     65.3,     51.7,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]	
    elseif Z == 79   wa = [80725.0,   14353.0, 13734.0,  11919.0,    3425.0,   3148.0,   2743.0,   2291.0,   2206.0,    762.1,    642.7,    546.3,
                                    353.2,    335.1,    87.6,     84.0,     107.2,     74.2,     57.2,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]	
    elseif Z == 80   wa = [83102.0,   14839.0, 14209.0,  12284.0,    3562.0,   3279.0,   2847.0,   2385.0,   2295.0,    802.2,    680.2,    576.6,
                                    378.2,    358.8,   104.0,     99.9,     127.0,     83.1,     64.5,     96.0,      7.8,     -1.0,     -1.0,      -1.0]	
    elseif Z == 81   wa = [85530.0,   15347.0, 14698.0,  12658.0,    3704.0,   3416.0,   2957.0,   2485.0,   2389.0,    846.2,    720.5,    609.5,
                                    405.7,    385.0,   122.2,    117.8,     136.0,     94.6,     73.5,     14.7,     12.5,     -1.0,     -1.0,      -1.0]
    elseif Z == 82   wa = [88005.0,   15861.0, 15200.0,  13035.0,    3851.0,   3554.0,   3066.0,   2586.0,   2484.0,    891.8,    761.9,    643.5,
                                    434.3,    412.2,   141.7,    136.9,     147.0,    106.4,     83.3,     20.7,     18.1,     -1.0,     -1.0,      -1.0]
    elseif Z == 83   wa = [90524.0,   16388.0, 15711.0,  13419.0,    3999.0,   3696.0,   3177.0,   2688.0,   2580.0,    939.0,    805.2,    678.8,
                                    464.0,    440.1,   162.3,    157.0,     159.3,    119.0,     92.6,     26.9,     23.8,     -1.0,     -1.0,      -1.0]
    elseif Z == 84   wa = [93105.0,   16939.0, 16244.0,  13814.0,    4149.0,   3854.0,   3302.0,   2798.0,   2683.0,    995.0,    851.0,    705.0,
                                    500.0,    473.0,   184.0,    184.0,     177.0,    132.0,    104.0,     31.0,     31.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 85   wa = [95730.0,   17493.0, 16785.0,  14214.0,    4317.0,   4008.0,   3426.0,   2909.0,   2787.0,    1042.0,   886.0,    740.0,
                                    533.0,    507.0,   210.0,    210.0,     195.0,    148.0,    115.0,     40.0,     40.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 86   wa = [98404.0,   18049.0, 17337.0,  14619.0,    4482.0,   4159.0,   3538.0,   3022.0,   2892.0,    1097.0,   929.0,    768.0,
                                    567.0,    541.0,   238.0,    238.0,     214.0,    164.0,    127.0,     48.0,     48.0,     26.0,     -1.0,      -1.0]	
    elseif Z == 87   wa = [101137.0,  18639.0, 17907.0,  15031.0,    4652.0,   4327.0,   3663.0,   3136.0,   3000.0,    1153.0,   980.0,    810.0,
                                    603.0,    577.0,   268.0,    268.0,     234.0,    182.0,    140.0,     58.0,     58.0,     34.0,     15.0,      15.0]
    elseif Z == 88   wa = [103922.0,  19237.0, 18484.0,  15444.0,    4822.0,   4490.0,   3792.0,   3248.0,   3105.0,    1208.0,  1058.0,    879.0,
                                    636.0,    603.0,   299.0,    299.0,     254.0,    200.0,    153.0,     68.0,     68.0,     44.0,     19.0,      19.0]
    elseif Z == 89   wa = [106755.0,  19840.0, 19083.0,  15871.0,    5002.0,   4656.0,   3909.0,   3370.0,   3219.0,    1269.0,  1080.0,    890.0,
                                    675.0,    639.0,   319.0,    319.0,     272.0,    215.0,    167.0,     80.0,     80.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 90   wa = [109651.0,  20472.0, 19693.0,  16300.0,    5182.0,   4830.0,   4046.0,   3491.0,   3332.0,    1330.0,  1168.0,    966.4,
                                    712.1,    675.2,   342.4,    333.1,     290.0,    229.0,    182.0,     92.5,     85.4,     41.4,     24.5,      16.6]
    elseif Z == 91   wa = [112601.0,  21105.0, 20314.0,  16733.0,    5367.0,   5001.0,   4174.0,   3611.0,   3442.0,    1387.0,  1224.0,   1007.0,
                                    743.0,    708.0,   371.0,    360.0,     310.0,    232.0,    232.0,     94.0,     94.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 92   wa = [115606.0,  21757.0, 20948.0,  17166.0,    5548.0,   5182.0,   4303.0,   3728.0,   3552.0,    1439.0,  1271.0,   1043.0, 
                                    778.3,    736.2,   388.2,    377.4,     321.0,    257.0,    192.0,    102.8,     94.2,     43.9,     26.8,      16.8] 
    #       
    else  error("Data not available for Z = $Z")
    end
    return( wa )
end

"""
`PeriodicTable.ionizationPotentials_Nist2025(Z::Int64)`  
    ... to return the `stored values' of the ionization potentials from NIST (2025) to readily estimate total and relative energies of atoms and ions.
"""
function ionizationPotentials_Nist2025(Z::Int64)
    # This procedure 'stores' the 1st, 2nd, 3rd, ... ionization potentials from NIST database 
    # (2025; http://physics.nist.gov/PhysRefData/ASD/ionEnergy.html) for all elements from Z = 1, ..., 92;
    # the q-th ionization potential leads from A^(q-1)+ --> A^q+; all ionization potentials are given in eV
    # 
    #       
	#      No.    IP               1st IP             2nd IP         3rd IP          4th IP         5th IP          6th             7th            8th            9th           10th
    #                            11/21/31 ..        12/22/32 ..    13/23/33 ..     14/24/34 ..    15/25/35 ..     16/26/36 ..     17/27/37 ..    18/28/38 ..    19/29/39 ..   20/30/40 ..
    if      Z ==  1   wa = [  13.598434599702]
	elseif  Z ==  2   wa = [     24.587389011,   54.4177655282]
	elseif  Z ==  3   wa = [      5.391714996,      75.6400970,   122.45435913]
	elseif  Z ==  4   wa = [         9.322699,        18.21115,     153.896205,   217.71858459]
	elseif  Z ==  5   wa = [         8.298019,        25.15483,       37.93059,     259.374379,   340.2260225]
	elseif  Z ==  6   wa = [       11.2602880,       24.383143,       47.88778,       64.49352,     392.09056,   489.99320779]
	elseif  Z ==  7   wa = [         14.53413,        29.60125,        47.4453,        77.4735,       97.8901,      552.06741,   667.0461377]
	elseif  Z ==  8   wa = [        13.618055,        35.12112,       54.93554,       77.41350,      113.8990,       138.1189,     739.32697,   871.4099138]
    elseif  Z ==  9   wa = [         17.42282,        34.97081,       62.70798,         87.175,       114.249,      157.16311,      185.1868,      953.8983,   1103.1175302]
    elseif  Z == 10   wa = [        21.564541,        40.96297,        63.4233,        97.1900,       126.247,        157.934,       207.271,      239.0970,      1195.8082,   1362.199256]
    elseif  Z == 11   wa = [       5.13907696,        47.28636,        71.6200,         98.936,       138.404,         172.23,       208.504,       264.192,        299.856,     1465.0992,  
                                  1648.702285]
    elseif  Z == 12   wa = [         7.646236,        15.03527,        80.1436,       109.2654,        141.33,         186.76,        225.02,       265.924,         327.99,       367.489,   
                                    1761.8049,     1962.663889]
    elseif  Z == 13   wa = [         5.985769,        18.82855,      28.447642,       119.9924,      153.8252,         190.49,        241.76,        284.64,         330.21,        398.65, 
                                      442.005,      2085.97693,    2304.140359]
    elseif  Z == 14   wa = [          8.15168,        16.34585,       33.49300,       45.14179,       166.767,        205.279,        246.57,        303.59,         351.28,        401.38, 
                                      476.273,         523.415,     2437.65805,    2673.177958]
    elseif  Z == 15   wa = [        10.486686,        19.76949,       30.20264,       51.44387,      65.02511,        220.430,        263.57,        309.60,         372.31,        424.40, 
                                       479.44,          560.62,        611.741,     2816.90868,   3069.842145]
    elseif  Z == 16   wa = [       10.3600167,        23.33788,          34.86,         47.222,       72.5945,        88.0529,       280.954,       328.794,         379.84,         447.7, 
                                       504.55,          564.41,         651.96,        706.994,    3223.78057,    3494.188518]
    elseif  Z == 17   wa = [        12.967633,        23.81364,          39.80,          53.24,         67.68,          96.94,      114.2013,       348.306,        400.851,         456.7, 
                                        530.0,          591.58,         656.30,         750.23,       809.198,     3658.34366,    3946.29179]
    elseif  Z == 18   wa = [       15.7596119,        27.62967,         40.735,          59.58,         74.84,         91.290,        124.41,      143.4567,         422.60,        479.96, 
                                        540.4,           619.0,          685.5,         755.13,         855.5,        918.375,    4120.66559,    4426.22407]
    elseif  Z == 19   wa = [       4.34066373,        31.62500,        45.8031,         60.917,         82.66,          99.44,        117.56,        154.87,       175.8174,        503.67, 
                                        565.6,           631.1,          714.7,          786.3,        860.92,         967.7,       1034.542,    4610.80714,     4934.04979]
    elseif  Z == 20   wa = [     6.1131549210,       11.871719,       50.91316,        67.2732,         84.34,        108.78,         127.71,        147.24,         188.54,       211.275, 
                                       591.60,           658.2,          728.6,          817.2,         894.0,         973.7,         1086.8,      1157.726,      5128.8576,    5469.86358]
    elseif  Z == 21   wa = [          6.56149,        12.79977,      24.756839,        73.4894,         91.95,        110.68,         137.99,        158.08,         180.03,        225.18, 
                                      249.798,          687.36,          757.7,          833.2,         926.5,        1008.6,         1093.5,        1213.1,       1287.957,     5674.9036, 
                                   6033.75643]
    elseif  Z == 22   wa = [         6.828120,         13.5755,       27.49171,       43.26717,        99.299,       119.533,         140.68,         170.5,          192.1,        215.92, 
                                       265.07,         291.500,         787.67,          864.0,         944.5,        1042.5,         1130.2,        1220.3,         1346.3,      1425.257, 
                                    6249.0226,      6625.81023]
    elseif  Z == 23   wa = [         6.746187,          14.634,        29.3111,         46.709,      65.28165,       128.125,         150.72,        173.55,          206.0,         230.5, 
                                        254.8,           308.5,        336.274,          896.0,         977.2,        1062.9,         1165.2,        1258.9,         1354.2,        1486.7, 
                                     1569.656,       6851.3112,     7246.12624]
    elseif  Z == 24   wa = [          6.76651,       16.486305,         30.959,          49.16,         69.46,       90.6349,         160.29,        184.76,          209.5,         244.5, 
                                        270.8,           296.7,          354.7,        384.163,        1011.6,        1097.2,         1188.0,        1294.8,         1394.5,        1495.1, 
                                       1634.1,        1721.183,      7481.8624,     7894.80289]
    elseif  Z == 25   wa = [        7.4340380,        15.63999,         33.668,          51.21,         72.41,        95.604,        119.203,         195.5,         221.89,        248.6, 
                                        286.1,           314.4,          343.6,         402.95,       435.172,        1133.7,         1224.1,        1320.3,         1430.9,       1537.2, 
                                       1643.2,          1788.7,       1879.873,      8140.4864,    8571.95438]
    elseif  Z == 26   wa = [        7.9024681,        16.19921,         30.651,          54.91,         75.00,        98.985,       124.9671,       151.060,          233.6,       262.10, 
                                        290.9,           330.8,          361.0,          392.2,         456.2,       489.312,         1262.7,        1357.8,           1460,       1575.6, 
                                       1687.0,          1798.4,         1950.4,       2045.759,     8828.1864,     9277.6886]
    elseif  Z == 27   wa = [          7.88101,         17.0844,          33.50,          51.27,         79.50,        102.00,         128.90,        157.80,         186.14,        275.4, 
                                       305.32,           336.1,          378.5,          410.0,         441.1,        511.96,        546.588,        1397.2,         1504.5,         1606, 
                                         1724,          1844.0,         1960.8,         2119.4,      2218.876,     9544.1817,     10012.1297]
    elseif  Z == 28   wa = [         7.639878,       18.168838,         35.187,          54.92,         76.06,         108.0,          132.0,         162.0,          193.2,        224.7, 
                                        319.5,           351.6,          384.5,          429.3,         462.8,         495.4,         571.07,       607.020,         1540.1,         1646, 
                                         1758,            1880,         2008.1,         2130.5,        2295.6,      2399.259,     10288.8848,    10775.3948]
    elseif  Z == 29   wa = [         7.726380,        20.29239,         36.841,          57.38,          79.8,         103.0,          139.0,         166.0,          198.0,        232.2, 
                                       265.33,           367.0,          401.0,          436.0,         483.1,         518.7,          552.8,         632.5,        670.608,       1690.5, 
                                         1800,            1918,           2044,         2179.4,        2307.3,        2479.1,       2586.954,    11062.4309,     11567.6237]
    elseif  Z == 30   wa = [         9.394197,        17.96439,       39.72330,         59.573,          82.6,         108.0,          133.9,         173.9,        203.0,          238.0, 
                                        274.4,           310.8,          417.6,          453.4,         490.6,         540.0,          577.8,         613.3,        697.5,        737.366, 
                                       1846.8,            1961,           2085,           2214,        2358.0,        2491.5,         2669.9,      2781.996,    11864.9401,    12388.9427]
    elseif  Z == 31   wa = [        5.9993020,        20.51514,       30.72576,         63.241,         86.01,         112.7,          140.8,         169.9,        211.0,          244.0, 
                                        280.0,           319.0,            356,          471.2,         508.8,         548.3,          599.8,         640.0,          677,          765.7, 
                                      807.308,          2010.0,           2129,           2258,          2391,        2543.9,         2683.0,          2868,     2984.426,     12696.5581, 
                                   13239.5029]
    elseif  Z == 32   wa = [         7.899435,       15.934610,        34.0576,        45.7155,        90.500,        115.90,          144.9,         176.4,        212.5,          252.1, 
                                          286,           326.0,            367,            407,         527.9,         567.3,          609.1,         662.8,        706.7,            744, 
                                        837.1,          880.44,         2178.2,           2304,          2439,          2575,         2737.1,        2881.9,         3074,       3194.293, 
                                   13557.4218,      14119.4457]
    elseif  Z == 33   wa = [          9.78855,         18.5892,         28.349,          50.15,         62.77,        121.19,          147.0,         180.0,        213.0,          247.0, 
                                        296.0,           333.0,            375,            418,           460,         587.6,          628.8,         672.9,        728.9,          774.0, 
                                          814,           911.7,         956.79,         2356.9,          2486,          2626,           2766,          2938,       3088.1,           3287, 
                                     3411.643,      14447.6799,    158028.9521]
    elseif  Z == 34   wa = [         9.752368,          21.196,         31.697,         42.947,         68.30,         81.83,        155.327,         184.0,        219.0,          255.0, 
                                        291.0,           342.9,            383,            426,           473,           517,          650.5,         693.4,        739.8,            798, 
                                        845.8,             887,          989.6,        1036.36,        2540.7,          2674,           2820,          2964,         3146,         3301.8, 
                                         3507,        3636.526,      15367.493,     15968.1075]
    elseif  Z == 35   wa = [         11.81381,          21.591,         34.871,         47.782,        59.595,        87.390,         103.03,        192.61,        224.0,          261.0, 
                                        301.0,           338.0,            393,            436,           481,           530,            577,         716.3,        761.0,          809.8, 
                                          870,           920.8,            963,         1070.6,       1119.17,        2731.4,           2869,          3021,         3169,           3361, 
                                       3523.1,            3735,       3868.986,      16317.014,    16937.1497]
    elseif  Z == 36   wa = [       13.9996055,        24.35984,         35.838,          50.85,         64.69,         78.49,         109.13,       125.802,        233.0,            268, 
                                          308,             350,            391,            446,           492,           540,            591,           640,      785.316,          831.6, 
                                        882.8,             945,          999.0,           1042,        1155.0,       1205.23,         2928.9,          3072,         3228,           3380, 
                                         3584,          3752.0,           3971,           4109.083, 17296.424,    17936.2405]
    elseif  Z == 37   wa = [        4.1771281,        27.28954,         39.247,         52.20,          68.44,          82.9,          98.67,        132.79,      150.628,         277.12, 
                                        313.1,           356.0,            400,            443,           502,           550,            601,           654,        706.0,            857, 
                                        905.3,           958.9,           1024,           1080,          1125,        1242.5,        1294.57,        3133.3,         3281,           3443, 
                                         3600,            3815,           3988,           4214,      4356.865,     18305.886,     18965.5484]
    elseif  Z == 38   wa = [       5.69486745,      11.0302765,       42.88353,         56.280,          70.7,          88.0,          104.0,        121.21,       158.33,         177.30, 
                                       324.07,           362.0,            408,            454,           499,           562,            612,           665,          722,            774, 
                                          932,           982.1,         1038.0,           1105,          1165,          1211,         1333.4,       1387.19,       3344.7,           3497, 
                                         3664,            3830,           4053,           4232,          4465,      4612.397,      19345.590,    20025.2673]
    elseif  Z == 39   wa = [          6.21726,         12.2236,       20.52441,        60.6072,         75.35,        91.390,         110.02,        128.12,       145.64,          185.7, 
                                      205.814,          374.04,            414,            463,           512,           559,            624,           677,          733,            790, 
                                          847,            1010,         1061.9,         1120.2,          1190,          1253,           1300,        1427.6,      1483.12,         3562.9,        
                                         3720,            3892,           4060,           4299,          4484,          4724,       4875.731,     20415.719,    21115.588]
    elseif  Z == 40   wa = [         6.634126,           13.13,         23.170,       34.41836,        80.348,         96.38,          112.0,         133.7,        153.0,         172.02, 
                                        214.9,         236.252,            426,            470,           520,           573,            622,           690,          745,            803, 
                                          863,             922,           1092,         1144.7,        1205.4,          1277,           1344,          1392,       1525.1,        1582.37, 
                                       3788.0,            3950,           4127,           4300,          4553,          4744,           4991,      5146.935,    21516.471,      22236.712]
    elseif  Z == 41   wa = [          6.75885,           14.32,          25.04,         37.611,       50.5728,       102.069,          119.1,         136.0,        159.2,          180.0, 
                                       200.28,           246.1,         268.59,          482.5,           530,           581,            636,           688,          758,            816, 
                                          877,             940,           1000,           1176,        1230.6,        1293.7,           1368,          1439,         1488,         1625.9, 
                                      1684.97,          4020.1,           4187,           4369,          4540,          4815,           5011,          5265,     5426.066,      22648.046, 
                                    23388.850]
    elseif  Z == 42   wa = [          7.09243,           16.16,          27.13,          40.33,        54.417,      68.82704,        125.638,         143.6,       164.12,          186.3, 
                                        209.3,          230.28,          279.1,         302.60,         544.0,           591,            646,           702,          758,            829, 
                                          890,             953,           1019,           1082,          1263,        1319.6,         1385.1,          1462,         1537,           1587, 
                                       1730.1,         1790.93,         4259.0,           4430,          4618,          4800,           5084,          5287,         5548,       5713.194, 
                                    23810.653,       24572.231]
    elseif  Z == 43   wa = [          7.11938,           15.26,          29.55,           41.0,          57.0,          72.0,           88.0,         150.0,        169.0,          189.9, 
                                        214.0,           239.0,         262.08,            311,        338.55,           604,            655,           713,          773,            829, 
                                          904,             968,           1032,           1102,          1166,          1354,         1411.6,        1479.5,         1559,           1638, 
                                         1689,            1838,        1900.28,         4505.0,          4681,          4874,           5060,          5361,         5570,           5838, 
                                     6008.391,       25004.531,      25787.047]
    elseif  Z == 44   wa = [          7.36050,           16.76,          28.47,           45.0,          59.0,          76.0,           93.0,         110.0,       178.41,          198.0, 
                                        219.9,           245.0,            271,         295.90,           348,        376.25,            670,           723,          784,            845, 
                                          905,             981,           1048,           1115,          1187,          1253,           1447,        1506.7,       1577.0,           1659, 
                                         1743,            1794,           1949,        2013.14,          4758,          4939,           5136,          5330,         5647,           5861, 
                                         6137,        6311.721,      26229.888,      27033.564]
    elseif  Z == 45   wa = [          7.45890,           18.08,          31.06,           42.0,          63.0,          80.0,           97.0,         115.1,        135.0,         207.51, 
                                        228.0,           252.1,            277,            306,        331.58,         389.3,         415.97,           739,          794,            857, 
                                          921,             984,           1061,           1131,          1202,          1274,           1344,          1544,       1604.9,         1677.6, 
                                         1763,            1851,           1903,           2063,       2129.22,          5018,           5203,          5406,         5600,           5940, 
                                         6161,            6444,       6623.262,      27486.979,     28312.031]
    elseif  Z == 46   wa = [         8.336839,           19.43,          32.93,           46.0,          61.0,          84.1,          101.0,         120.0,        141.0,          159.9, 
                                       238.57,           260.0,            286,            311,           342,        369.10,            427,        457.50,          810,            869, 
                                          933,            1000,           1065,           1145,          1218,          1290,           1366,          1438,         1644,         1706.2, 
                                       1781.3,            1869,           1962,           2016,          2181,       2248.87,           5284,          5475,         5683,           5880, 
                                         6242,            6469,           6759,       6943.097,     28776.032,     29622.678]
    elseif Z == 47    wa = [         7.576234,         21.4844,           34.8,           49.0,          65.0,          82.0,          106.0,         125.0,        145.1,          167.0, 
                                        188.0,          271.46,            294,            321,           347,           381,         408.43,           469,       500.87,            885, 
                                          946,            1013,           1082,           1149,          1231,          1308,           1382,          1460,         1535,           1747, 
                                       1810.5,          1888.0,           1979,           2077,          2131,          2302,        2371.99,          5558,         5753,           5966, 
                                         6170,            6551,           6785,           7082,      7271.298,     30097.314,      30965.780]
    elseif  Z == 48   wa = [         8.993820,       16.908313,         37.468,           51.0,          67.9,          87.0,          105.0,         130.1,        150.0,          173.0, 
                                        195.0,           218.0,            305,            329,           358,           385,            421,         452.6,          513,         546.19, 
                                          963,            1026,           1095,           1167,          1237,          1320,           1401,          1477,         1558,           1635, 
                                         1852,          1917.9,           1998,           2091,          2195,          2250,           2427,       2498.62,         5839,           6039, 
                                         6257,            6460,           6869,           7109,          7414,       7607.95,      31451.070,     32341.587]
    elseif  Z == 49   wa = [        5.7863558,        18.87041,       28.04415,          55.45,         69.3,           90.0,          109.0,         130.1,        156.0,          178.0, 
                                        201.0,           226.0,            249,            341,           368,           396,            425,           462,        497.1,            560, 
                                       593.38,            1043,           1109,           1181,          1255,          1328,           1413,          1496,         1575,           1659, 
                                         1738,            1961,         2028.5,           2111,          2207,          2317,           2373,          2555,      2628.77,           6126, 
                                         6331,            6554,           6770,           7196,          7442,          7754,        7953.14,     32837.593,    33750.404]
    elseif  Z == 50   wa = [         7.343918,        14.63307,         30.506,          40.74,         77.03,          94.0,          112.9,         135.0,        156.0,          184.0, 
                                        208.0,           232.0,            258,            282,           379,           407,            437,           466,          506,            537, 
                                          608,          642.35,           1127,           1195,          1269,          1347,           1421,          1508,         1596,           1676, 
                                         1763,            1844,           2074,         2142.1,          2227,          2326,           2443,          2499,         2687,        2762.49, 
                                         6421,            6631,           6859,           7080,          7531,          7790,           8103,       8306.95,    34257.148,      35192.501]
    elseif  Z == 51   wa = [         8.608389,          16.626,        25.3235,         43.804,         55.00,         99.51,          117.0,         139.0,        162.0,          185.0, 
                                        214.0,           238.0,            265,            292,           317,           420,            447,           479,          510,            552, 
                                          584,             657,         693.26,           1214,          1285,          1360,           1441,          1518,         1606,           1698, 
                                         1781,            1869,           1954,           2190,          2266,          2349,           2428,          2567,         2654,           2815, 
                                         2900,            6714,           6929,           7167,          7390,          7887,           8140,          8455,      8699.48,      35710.030, 
                                    36668.183]
    elseif  Z == 52   wa = [         9.009808,            18.6,          27.84,        37.4155,          59.3,          69.1,         124.20,         143.0,        167.0,          191.1, 
                                        215.0,           245.0,            272,            299,           328,           354,            461,           491,          522,            555, 
                                          599,             633,            709,         746.12,          1304,          1377,           1455,          1538,         1618,           1707, 
                                         1803,            1889,           1979,           2066,          2309,          2386,           2472,          2552,         2700,           2788, 
                                         2954,          3041.0,           7022,           7243,          7485,          7714,           8240,          8499,         8821,        9040.83, 
                                     37196.52,       38177.740]
    elseif  Z == 53   wa = [        10.451236,        19.13126,         29.570,         40.357,         51.52,          74.4,          87.61,        150.81,        171.0,          197.0, 
                                        220.9,           247.0,            279,            307,           335,           365,            393,           505,          535,            569, 
                                          601,             649,            683,            762,         800.8,          1397,           1472,          1553,         1639,           1720, 
                                         1812,            1911,           1999,           2093,          2181,          2431,           2510,          2598,         2680,           2836, 
                                         2926,            3096,         3185.5,           7337,          7563,          7811,           8044,          8601,         8867,           9196, 
                                      9421.10,        38717.00,      39721.549]
    elseif  Z == 54   wa = [       12.1298437,          20.975,          31.05,          42.20,          54.1,         66.703,          91.6,      105.9778,       179.84,          202.0, 
                                       229.02,           255.0,            281,            314,           343,            374,           404,           434,          549,            582, 
                                          616,             650,            700,            736,           818,            857,          1493,          1571,         1653,           1742, 
                                         1826,            1919,           2023,           2113,          2209,           2300,          2556,          2637,         2726,           2811, 
                                         2975,            3068,           3243,         3333.8,          7660,           7889,          8144,          8382,         8971,           9243, 
                                         9581,         9810.37,       40271.73,      41299.892]
    elseif  Z == 55   wa = [      3.893905727,        23.15745,         33.195,           43.0,          56.0,           69.1,          82.9,         110.1,       125.61,          213.3, 
                                        233.0,           261.0,            289,            316,           352,            382,           413,           445,          476,            597, 
                                          629,             666,            700,            753,           791,            875,         916.1,          1592,         1672,           1757, 
                                         1848,            1936,           2029,           2137,          2230,           2329,          2422,          2683,         2767,           2859, 
                                         2945,            3118,           3214,           3392,          3485,           7989,          8224,          8484,         8726,           9350, 
                                         9629,            9974,       10208.78,       41861.08,     42913.144]
    elseif  Z == 56   wa = [        5.2116646,       10.003826,        35.8438,           47.0,          58.0,           71.0,          86.0,         101.0,        130.5,         146.52, 
                                        241.0,           267.1,            296,            325,           354,            390,           422,           455,          488,            520, 
                                          646,             679,            717,            752,           809,            846,           935,        976.62,         1695,           1776, 
                                         1864,            1958,           2047,           2142,          2256,           2349,          2452,          2547,         2814,           2901, 
                                         2994,            3081,           3266,           3363,          3546,           3640,          8326,          8565,         8831,           9077, 
                                         9739,           10023,          10376,       10616.42,      43485.37,      44561.633]
	elseif  Z == 57   wa = [           5.5769,        11.18496,        19.1773,          49.95,          61.6,           74.0,          88.0,         105.0,        119.0,          151.4, 
	                                   168.77,           275.0,            303,            332,           364,            393,           431,           464,          498,            533, 
	                                 566, 696,             731,            770,            806,           865,            906,           995,       1039.09,         1800,           1884, 
	                                     1974,            2069,           2162,           2259,          2377,           2473,          2577,          2674,         2950,           3036,  
	                                     3133,            3222,           3416,           3515,          3704,           3800,          8669,          8914,         9184,           9437, 
	                                    10136,           10426,          10789,       11033.40,      45145.00,       46245.77]
	elseif  Z == 58   wa = [           5.5386,          10.956,        20.1974,         36.906,         65.55,           77.6,          91.0,         106.0,        125.0,          140.0, 
	                                    172.0,          192.24,            312,            340,           371,            403,           435,           472,          509,            543, 
	                                      579,             613,            749,            785,           824,            862,           924,           965,         1060,         1103.5, 
	                                     1908,            1994,           2087,           2185,          2280,           2378,          2500,          2600,         2706,           2806, 
	                                     3087,            3176,           3274,           3366,          3570,           3672,          3865,          3963,         9020,           9269, 
	                                     9545,            9803,          10542,          10840,         11210,       11459.85,      46840.31,      47965.89]
	elseif  Z == 59   wa = [           5.4702,          10.631,        21.6237,         38.981,         57.53,             82,            97,           112,          131,            148, 
	                                      162,             196,         217.02,            350,           378,            412,           445,           478,          516,            554, 
	                                      590,             627,            663,            803,           840,            880,           920,           985,         1028,           1124, 
	                                   1169.9,            2019,           2108,           2202,          2304,           2400,          2501,          2628,         2729,           2838, 
	                                     2941,            3227,           3319,           3419,          3512,           3729,          3832,          4030,         4130,           9378, 
                                         9632,            9913,          10175,          10959,         11262,          11641,      11895.89,      48571.71,     49722.44]
	elseif  Z == 60   wa = [          5.52475,          10.783,          22.09,          40.60,          60.0,             84,            99,           114,          136,            152, 
	                                      168,             195,            221,          243.0,           389,            420,           453,           489,          522,            562, 
	                                      602,             638,            678,            714,           859,            896,           939,           978,         1049,           1092, 
	                                     1191,         1238.42,           2134,           2224,          2321,           2425,          2525,          2627,         2758,           2861, 
	                                     2974,            3078,           3371,           3465,          3567,           3662,          3891,          3997,         4198,           4302, 
	                                     9742,           10002,          10288,          10555,         11384,          11694,         12082,      12341.66,     50339.59,       51515.78]
	elseif  Z == 61   wa = [          5.58187,          10.938,          22.44,          41.17,          61.7,             85,           101,           116,          138,            155, 
                                          174,             202,            229,            248,           269,            430,           462,           497,          534,            569, 
                                          609,             651,            689,            730,           767,            916,           956,           998,         1040,           1113, 
                                         1158,            1261,         1308.7,           2251,          2344,           2443,          2549,          2652,         2755,           2892, 
                                         2997,            3112,           3219,           3519,          3613,           3718,          3816,          4056,         4166,           4371, 
                                         4476,           10115,          10378,          10671,         10942,          11819,         12136,         12532,     12797.26,       52144.29, 
                                     53346.31]
	elseif  Z == 62   wa = [         5.643722,          11.078,          23.55,          41.64,          62.7,             87,           103,           118,          141,            158, 
                                          179,             208,            237,            257,           276,          306.5,           474,           506,          543,            581, 
                                          617,             658,            702,            742,           782,            822,           976,          1016,         1060,           1103, 
                                         1180,            1226,           1332,        1381.56,          2371,           2466,          2569,          2676,         2782,           2887, 
                                         3028,            3137,           3253,           3363,          3669,           3766,          3873,          3971,         4227,           4337,
                                         4548,            4655,          10494,          10762,         11060,          11337,         12264,         12588,        12992,       13262.85,
                                     53986.12,        55214.30]
	elseif  Z == 63   wa = [         5.670385,          11.240,          24.84,          42.94,          63.2,             89,           105,           120,          144,            161,  
                                          183,             213,            243,            263,           281,            311,         344.4,           518,          553,            590, 
                                          630,             667,            709,            755,           795,            838,           879,          1037,         1078,           1124, 
                                         1167,            1249,           1296,           1406,       1456.06,           2495,          2591,          2697,         2807,           2914, 
                                         3022,            3168,           3279,           3398,          3510,           3823,          3921,          4031,         4131,           4400,
                                         4513,            4729,           4838,          10880,         11153,          11457,         11739,         12718,        13050,          13462, 
                                     13738.58,        55865.93,       57120.64]
	elseif  Z == 64   wa = [          6.14980,          12.076,          20.54,          44.44,           64.8,            89,           106,           123,          144,            165, 
	                                      183,             213,            246,            268,             288,          319,           352,         384.4,          565,            601, 
	                                      639,             680,            719,            761,             810,          851,           895,           937,         1100,           1142, 
	                                     1189,            1233,           1321,           1368,            1481,       1532.3,          2621,          2720,         2827,           2941, 
	                                     3050,            3160,           3312,           3424,            3546,         3660,          3980,          4080,         4191,           4294,
                                         4578,            4693,           4914,           5025,           11273,        11552,         11861,         12147,        13183,          13521, 
	                                    13943,        14224.57,       57783.91,       59065.54]
	elseif  Z == 65   wa = [           5.8638,          11.513,          21.82,          39.33,            66.5,           90,           108,           125,          143,            168,
                                          186,             216,            250,            273,             294,          325,           358,           393,        426.6,            613,            
                                          651,             690,            732,            772,             816,          866,           909,           954,          997,           1165, 
                       	                 1208,            1256,           1301,           1393,            1443,         1559,        1610.4,          2750,         2852,           2961, 
                                         3078,            3189,           3300,           3458,            3573,         3698,          3814,          4139,         4242,           4355, 
                                         4460,            4760,           4877,           5103,            5217,        11673,         11957,         12272,        12563,          13658, 
                                        14003,           14434,       14721.02,       59741.12,         61050.1]
	elseif  Z == 66   wa = [         5.939061,          11.647,          22.89,          41.23,            62.1,           93,           110,           127,          152,            170,
                                          192,             224,            259,            279,             300,          332,           366,           399,          431,          464.9,
                                          664,             702,            743,            786,             827,          872,           924,           969,         1014,           1059,
                                         1232,            1275,           1325,           1371,            1468,         1520,          1638,        1691.7,         2882,           2987,
                                         3098,            3217,           3331,           3445,            3607,         3725,          3852,          3970,         4303,           4407,
                                         4523,            4629,           4945,           5066,            5296,         5412,         12081,         12370,        12690,          12986, 
                                        14144,           14495,          14936,       15228.06,        61736.62,     63073.23]
	elseif  Z == 67   wa = [           6.0215,          11.781,          22.79,          42.52,            63.9,           95,           112,           129,          155,            173, 
                                          197,             229,            263,            284,             305,          340,           373,           408,          441,            475, 
                                          510,             715,            755,            797,             842,          885,           929,           985,         1029,           1077, 
                                         1122,            1300,           1346,           1395,            1443,         1545,          1598,          1719,       1773.6,           3018,
                                         3125,            3238,           3359,           3476,            3592,         3760,          3880,          4009,         4131,           4469,
                                         4576,            4693,           4802,           5135,            5258,         5494,          5611,         12495,        12790,          13116,
                                        13417,           14639,          14998,          15448,        15745.77,     63772.42,      65137.13]
	elseif  Z == 68   wa = [           6.1077,          11.916,          22.70,          42.42,            65.1,           96,           114,           131,          158,            177,
                                          201,             235,            268,            290,             311,          345,           381,           415,          450,            486,
                                          520,             555,            770,            810,             853,          899,           943,           989,         1046,           1092,
                                         1142,            1188,           1370,           1416,            1468,         1516,          1625,          1678,         1803,         1858.5,
                                         3157,            3265,           3381,           3505,            3624,         3742,          3916,          4038,         4170,           4294, 
                                         4639,            4748,           4866,           4978,            5329,         5455,          5695,          5815,        12918,          13217,
                                        13548,           13855,          15146,          15511,           15971,     16274.56,      65848.24,      67241.48]
	elseif  Z == 69   wa = [         6.184402,          12.065,          23.66,          42.41,            65.4,           98,           116,           133,          160,            180,
                                          205,             239,            274,            295,             317,          352,           387,           424,          460,            496, 
                                          530,             570,            603,            825,             866,          911,           958,          1004,         1050,           1110,
                                         1157,            1207,           1255,           1442,            1490,         1542,          1591,          1706,         1761,           1889,
                                       1945.2,            3298,           3409,           3528,            3653,         3775,          3895,          4075,         4199,           4335,
                                         4461,            4812,           4922,           5044,            5157,         5527,          5656,          5901,         6023,          13347, 
                                        13651,           13988,          14300,          15663,           16036,        16510,      16814.34,      67965.25,     69387.45]
	elseif  Z == 70   wa = [         6.254160,       12.179185,         25.053,          43.61,            65.6,           99,           117,           135,          163,            182, 
                                          209,             244,            279,            301,             324,          360,           396,           431,          469,            505,
                                          540,             580,            610,            651,             882,          924,           971,          1019,         1065,           1114,  
                                         1175,            1224,           1275,           1324,            1516,         1564,          1618,          1668,         1789,           1845,
                                         1978,          2036.4,           3443,           3555,            3677,         3805,          3929,          4051,         4238,           4364,
                                         4502,            4630,           4988,           5101,            5224,         5339,          5731,          5860,         6111,           6236,
                                        13784,           14093,          14435,          14752,           16191,        16570,         17050,      17365.44,     70123.04,       71574.63]
    elseif  Z == 71   wa = [         5.425871,           14.13,        20.9594,         45.249,            66.8,           98,           117,           136,          159,            185, 
                                          205,             238,            276,            305,             328,          361,           399,           438,          476,            520,           
                                          560,             600,            630,            670,             713,          941,           985,          1032,         1081,           1130,
                                         1178,            1242,           1292,           1345,            1395,         1591,          1641,          1696,         1747,           1875,
                                         1933,            2067,         2125.5,           3590,            3706,         3828,          3960,          4086,         4211,           4403,
                                         4532,            4673,           4803,           5168,            5282,         5408,          5525,          5937,         6070,           6326,
                                         6452,           14228,          14542,          14890,           15211,        16730,         17120,         17610,     17928.05,       72322.87,
                                     73804.35]
	elseif  Z == 72   wa = [         6.825070,           14.61,          22.55,         33.370,           68.37,           98,           118,           137,          157,            187, 
                                          209,             230,            270,            310,             334,          359,           399,           440,          481,            520,
                                          570,             610,            650,            690,             730,          772,          1002,          1047,         1094,           1146,
                                         1195,            1245,           1311,           1362,            1417,         1467,          1669,          1719,         1776,           1827,
                                         1963,            2022,           2159,         2218.9,            3741,         3858,          3984,          4118,         4246,           4372, 
                                         4573,            4703,           4846,           4980,            5350,         5468,          5595,          5713,         6149,           6284, 
                                         6545,            6674,          14678,          14999,           15351,        15680,         17280,         17680,        18180,       18502.32,
                                     74565.91,        76077.70]
    elseif  Z == 73   wa = [         7.549571,            16.2,           23.1,           35.0,          48.272,        94.01,           119,           139,          159,            180, 
                                          213,             235,            262,            304,             338,          363,           396,           439,          482,            530,
                                          570,             610,            660,            700,             750,          790,           832,          1064,         1110,           1160, 
                                         1211,            1262,           1313,           1382,            1434,         1490,          1542,          1748,         1799,           1857,
                                         1910,            2053,           2113,           2254,          2314.7,       3898.7,          4014,          4143,         4278,           4410, 
                                         4537,            4745,           4877,           5024,            5159,         5537,          5655,          5785,         5907,           6364, 
                                         6502,            6769,           6900,          15137,           15461,        15820,         16150,         17840,        18250,          18760,
                                     19088.51,        76852.00,       78394.63]
	elseif  Z == 74   wa = [          7.86403,           16.37,           26.0,           38.2,            51.6,        64.77,        122.01,         141.2,        160.2,          179.0, 
                                        208.9,           231.6,          258.3,          290.7,           325.3,        361.9,         387.9,         420.7,        462.1,          502.6,
                                        543.4,           594.5,          640.6,          685.6,           734.1,        784.4,         833.4,         881.4,       1132.2,         1180.0, 
                                       1230.4,          1283.4,         1335.1,         1386.8,          1459.9,       1512.4,        1569.1,        1621.7,       1829.8,         1882.9, 
                                       1940.6,          1994.8,         2149.1,         2210.0,          2354.5,       2414.1,          4057,          4180,         4309,           4446,
                                         4578,            4709,           4927,           5063,            5209,         5348,          5719,          5840,         5970,           6093, 
                                         6596,            6735,           7000,           7130,           15566,        15896,         16252,         16588,        18476,          18872,
                                        19362,        19686.74,       79181.94,       80755.91]
	elseif  Z == 75   wa = [          7.83352,            16.6,           27.0,           39.1,            51.9,         67.0,         82.71,         144.4,          165,            187,
                                          208,             236,            268,            291,             330,          377,           403,           429,          476,            520, 
                                          570,             620,            670,            720,             760,          810,           860,           910,          953,           1194, 
                                         1242,            1294,           1349,           1402,            1454,         1530,          1583,          1641,         1696,           1912, 
                                         1966,            2025,           2080,           2240,            2302,         2450,        2514.5,          4214,         4335,           4468, 
                                         4609,            4745,           4877,           5099,            5236,         5388,          5528,          5919,         6042,           6176, 
                                         6300,            6810,           6952,           7230,            7366,        16080,         16410,         16780,        17120,          19000, 
                                        19420,           19950,       20297.40,       81556.58,        83162.41] 
    elseif  Z == 76   wa = [          8.43823,            17.0,           25.0,           41.0,            55.0,         70.1,          85.1,        102.02,        168.7,            190, 
                                          213,             235,            269,            298,             322,          367,           410,           436,          470,            520,
                                          570,             620,            670,            720,             770,          820,           870,           920,          970,           1015,
                                         1262,            1311,           1364,           1420,            1474,         1528,          1606,          1660,         1720,           1776,
                                         1996,            2052,           2112,           2168,            2336,         2400,          2552,        2615.5,         4374,           4510, 
                                         4635,            4779,           4917,           5052,            5280,         5421,          5575,          5717,         6115,           6240, 
                                         6376,            6503,           7039,           7185,            7468,         7610,         16560,         16900,        17270,          17620, 
                                        19600,           20030,          20570,       20920.60,        83976.18,     85614.42]
	elseif  Z == 77   wa = [          8.96702,            17.0,           28.0,           40.0,            57.0,         72.0,          89.0,         105.0,        122.7,          194.8, 
                                          217,             240,            264,            303,             329,          356,           407,           445,          472,            510, 
                                          560,             610,            670,            720,             770,          820,           870,           920,          980,           1030, 
                                         1080,            1331,           1381,           1436,            1493,         1548,          1603,          1684,         1739,           1801, 
                                         1857,            2083,           2139,           2201,            2258,         2435,          2500,          2656,       2720.4,           4540, 
                                         4668,            4806,           4952,           5092,            5229,         5466,          5609,          5765,         5910,           6315, 
                                         6441,            6580,           6708,           7274,            7421,         7710,          7850,         17040,        17390,          17770, 
                                        18120,           20210,          20650,          21200,        21556.60,     86442.44,       88113.6]
	elseif  Z == 78   wa = [          8.95883,           18.56,           29.0,           43.0,            56.0,         75.0,          91.0,         109.0,        126.0,          144.9, 
                                        220.4,             245,            269,            293,             332,          358,           392,           445,          479,            507, 
                                          550,             610,            660,            710,             760,          820,           870,           930,          980,           1040, 
                                         1090,            1140,           1402,           1454,            1509,         1567,          1624,          1680,         1763,           1821, 
                                         1883,            1941,           2171,           2228,            2291,         2350,          2536,          2603,         2762,         2827.8, 
                                         4715,            4839,           4980,           5128,            5270,         5410,          5654,          5800,         5959,           6106, 
                                         6517,            6646,           6787,           6918,            7512,         7660,          7960,          8100,        17540,          17890, 
                                        18280,           18630,          20840,          21280,           21840,      22205.7,       88955.1,      90659.84]
    elseif  Z == 79   wa = [         9.225554,          20.203,           30.0,           45.0,            60.0,         74.0,          94.0,         112.0,        130.1,          149.0, 
                                        168.2,           248.0,            275,            299,             324,          365,           392,           433,          487,            520, 
                                          550,             600,            650,            710,             760,          820,           870,           930,          990,           1040, 
                                         1100,            1150,           1210,           1475,            1527,         1584,          1644,          1702,         1758,           1845, 
                                         1904,            1967,           2026,           2261,            2320,         2383,          2443,          2640,         2708,           2870, 
                                       2941.0,            4888,           5013,           5156,            5307,         5452,          5594,          5846,         5994,           6156, 
                                         6305,            6724,           6854,           6997,            7130,         7760,          7910,          8210,         8360,          18040, 
                                        18400,           18790,          19150,          21470,           21920,        22500,       22868.1,       91515.8,     93254.62]
    elseif  Z == 80   wa = [        10.437504,        18.75688,          34.49,          48.55,           61.20,         76.6,          93.0,         113.9,        134.0,          153.0, 
                                        173.0,           192.7,          276.9,            307,             332,          357,           402,           429,          477,            530, 
                                          560,             590,            650,            710,             760,          820,           880,           930,          990,           1050, 
                                         1110,            1160,           1220,           1280,            1549,         1603,          1661,          1723,         1780,           1839, 
                                         1928,            1989,           2052,           2113,            2354,         2412,          2478,          2539,         2745,           2815, 
                                         2981,          3049.9,           5055,           5191,            5335,         5490,          5636,          5780,         6041,           6192, 
                                         6356,            6508,           6933,           7066,            7211,         7350,          8010,          8160,         8470,           8620, 
                                        18550,           18910,          19310,          19680,           22120,        22580,         23170,       23544.1,      94124.7,       95898.19]
    elseif  Z == 81   wa = [        6.1082873,         20.4283,        29.8520,          51.14,            62.6,         80.0,          97.9,         116.0,        135.0,          158.0, 
                                        177.0,           198.0,          218.3,          306.9,             340,          366,           392,           439,          467,            520, 
                                          570,             600,            640,            700,             760,          820,           880,           930,          990,           1060, 
                                         1110,            1170,           1230,           1290,            1350,         1625,          1681,          1740,         1802,           1862, 
                                         1920,            2014,           2075,           2140,            2202,         2447,          2508,          2574,         2635,           2854,
                                         2925,            3094,         3164.7,           5234,            5371,         5518,          5674,          5824,         5969,           6241, 
                                         6392,            6560,           6714,           7146,            7281,         7430,          7570,          8260,         8420,           8730, 
                                         8880,           19070,          19440,          19840,           20210,        22780,         23250,         23850,      24234.1,        96783.2, 
                                     98592.12]
    elseif  Z == 82   wa = [        7.4166799,       15.032499,        31.9373,       42.33256,            68.8,         82.9,         100.1,         120.0,        138.0,          158.0, 
                                        182.0,           203.0,            224,          245.1,           338.1,          374,           401,           427,          478,            507, 
                                          570,             610,            650,            690,             750,          810,           870,           930,          990,           1050, 
                                         1120,            1180,           1240,           1300,            1360,         1430,          1704,          1760,         1819,           1884, 
                                         1945,            2004,           2101,           2163,            2230,         2292,          2543,          2605,         2671,           2735, 
                                         2965,            3036,           3211,         3282.1,            5414,         5555,          5703,          5862,         6015,           6162, 
                                         6442,            6597,           6767,           6924,            7362,         7500,          7650,          7790,         8520,           8680, 
                                         9000,            9150,          19590,          19970,           20380,        20750,         23460,         23940,        24550,        24938.2, 
                                      99491.8,        101336.7]
	elseif  Z == 83   wa = [         7.285516,          16.703,       25.57075,          45.37,          54.856,         88.4,         103.0,         122.0,        143.0,          161.1, 
                                        183.0,           208.0,            229,            252,           272.6,        370.2,           409,           436,          464,            520, 
                                          550,             620,            660,            690,             750,          810,           870,           930,          990,           1060, 
                                         1120,            1180,           1250,           1310,            1380,         1440,          1500,          1784,         1840,           1902, 
                                         1967,            2029,           2090,           2190,            2253,         2321,          2385,          2641,         2703,           2771,
                                         2835,            3078,           3151,           3329,          3401.8,         5599,          5740,          5892,         6054,           6208, 
                                         6358,            6648,           6804,           6977,            7137,         7580,          7720,          7870,         8010,           8780, 
                                         8950,            9270,           9430,          20130,           20500,        20920,         21300,         24150,        24640,          25260, 
                                      25656.9,        102251.8,       104133.4]
	elseif  Z == 84   wa = [         8.418070,            19.3,           27.3,           36.0,            57.0,         69.1,         108.1,         125.0,        146.1,          166.0, 
                                        186.0,           209.0,            235,            257,             281,          304,           416,           444,          473,            502, 
                                          560,             590,            670,            700,             740,          800,           870,           930,          990,           1060, 
                                         1120,            1180,           1250,           1320,            1380,         1440,          1510,          1570,         1865,           1923, 
                                         1986,            2052,           2115,           2177,            2281,         2345,          2414,          2480,         2740,           2803, 
                                         2873,            2938,           3194,           3268,            3450,       3524.2,          5785,          5930,         6084,           6248, 
                                         6405,            6557,           6856,           7015,            7191,         7350,          7810,          7950,         8100,           8240, 
                                         9050,            9220,           9550,           9710,           20670,        21050,         21470,         21860,        24860,          25360, 
                                        25990,         26390.4,       105064.3,       106983.4]
	elseif  Z == 85   wa = [          9.31751,          17.880,          26.58,          39.65,           50.39,         72.0,          85.1,         130.1,        149.0,          169.0, 
                                        192.1,           212.0,            236,            263,             287,          311,           335,           452,          481,            510, 
                                          540,             600,            630,            720,             750,          790,           860,           920,          990,           1050, 
                                         1120,            1180,           1250,           1320,            1380,         1450,          1510,          1590,         1650,           1948, 
                                         2007,            2071,           2139,           2203,            2266,         2373,          2439,          2510,         2576,           2841, 
                                         2905,            2977,           3042,           3312,            3388,         3573,          3649,          5976,         6122,           6279, 
                                         6445,            6604,           6759,           7068,            7230,         7410,          7570,          8030,         8180,           8330, 
                                         8480,            9330,           9500,           9830,            9990,        21210,         21600,         22030,        22420,          25580,
                                        26090,           26730,        27139.0,       107930.0,        109887.2]
	elseif  Z == 86   wa = [         10.74850,           18.99,           29.4,           36.9,            52.9,         64.0,          88.0,         102.0,        154.0,          173.9, 
                                        195.0,           218.0,            240,            264,             293,          317,           342,           367,          488,            520, 
                                          550,             580,            640,            680,             760,          800,           850,           920,          980,           1050, 
                                         1110,            1180,           1250,           1310,            1390,         1460,          1520,          1590,         1660,           1720, 
                                         2033,            2094,           2158,           2227,            2293,         2357,          2467,          2535,         2606,           2674, 
                                         2944,            3010,           3082,           3149,            3433,         3510,          3699,          3777,         6169,           6318, 
                                         6476,            6646,           6807,           6964,            7283,         7450,          7630,          7800,         8260,           8410, 
                                         8570,            8710,           9610,           9780,           10120,        10290,         21770,         22160,        22600,          22990, 
                                        26310,           26830,          27490,        27903.1,        110846.3,     112842.2]
	elseif  Z == 87   wa = [        4.0727411,            22.4,           33.5,           39.1,            50.0,         67.0,          80.0,         106.0,        120.0,          179.0, 
                                        200.0,           222.1,            245,            269,             293,          324,           349,           375,          400,            530, 
                                          560,             590,            620,            690,             720,          810,           850,           910,          980,           1040, 
                                         1110,            1180,           1250,           1320,            1380,         1460,          1530,          1600,         1670,           1740, 
                                         1810,            2119,           2182,           2247,            2317,         2384,          2450,          2564,         2631,           2706, 
                                         2774,            3049,           3115,           3190,            3257,         3556,          3635,          3828,         3907,           6365, 
                                         6516,            6678,           6849,           7013,            7172,         7500,          7670,          7850,         8020,           8500, 
                                         8640,            8800,           8950,           9890,           10070,        10420,         10590,         22330,        22730,          23170, 
                                        23570,           27060,          27590,          28260,         28683.4,     113821.9,      115857.5]
	elseif  Z == 88   wa = [        5.2784239,        10.14718,           31.0,           41.0,            52.9,         64.0,          82.0,          97.0,        124.0,          140.0, 
                                        204.9,           227.0,            250,            274,             299,          324,           356,           382,          409,            435, 
                                          570,             600,            630,            660,             740,          770,           860,           900,          970,           1040, 
                                         1110,            1180,           1250,           1320,            1390,         1460,          1530,          1610,         1680,           1750, 
                                         1820,            1880,           2208,           2271,            2338,         2409,          2477,          2544,         2662,           2731, 
                                         2806,            2876,           3155,           3224,            3298,         3368,          3682,          3762,         3959,           4040, 
                                         6565,            6718,           6881,           7056,            7222,         7380,          7720,          7890,         8080,           8250, 
                                         8730,            8880,           9040,           9200,           10190,        10360,         10720,         10890,        22900,          23300, 
                                        23750,           24160,          27830,          28370,           29050,      29479.8,      116853.5,      118929.5]
	elseif  Z == 89   wa = [         5.380235,           11.75,         17.436,           44.8,            55.0,         67.0,          79.0,          98.9,        113.9,          143.9, 
	                                    161.1,           233.0,            255,            279,             305,          330,           355,           390,          416,            444, 
	                                      470,             610,            640,            670,             710,          780,           820,           920,          950,           1030, 
	                                     1100,            1170,           1240,           1310,            1380,         1460,          1530,          1610,         1680,           1750, 
	                                     1820,            1900,           1970,           2298,            2362,         2430,          2503,          2572,         2639,           2762, 
	                                     2833,            2908,           2980,           3264,            3334,         3409,          3479,          3811,         3893,           4093, 
	                                     4175,            6767,           6923,           7088,            7265,         7430,          7600,          7950,         8120,           8310, 
	                                     8480,            8970,           9120,           9290,            9440,        10480,         10660,         11030,        11200,          23480, 
	                                    23890,           24340,          24760,          28610,           29160,        29850,       30293.1,      119945.7,     122063.1]
	elseif  Z == 90  wa = [           6.30670,           12.10,          18.32,         28.648,            58.0,         69.1,          82.0,          95.0,        118.0,          133.0, 
	                                    165.0,           181.0,            262,            285,             310,          336,           362,           389,          424,            451, 
	                                      480,             508,            650,            680,             720,          750,           830,           870,          970,           1010, 
	                                     1090,            1160,           1240,           1310,            1380,         1460,          1530,          1600,         1680,           1760, 
	                                     1830,            1910,           1980,           2060,            2390,         2455,          2524,          2598,         2669,           2737, 
	                                     2864,            2935,           3013,           3086,            3375,         3445,          3522,          3593,         3943,           4025, 
	                                     4230,            4313,           6972,           7130,            7299,         7480,          7650,          7810,         8180,           8350, 
	                                     8550,            8720,           9220,           9370,            9540,         9690,         10790,         10970,        11340,          11510, 
                                        24060,           24480,          24940,          25360,           29410,        29970,         30680,       31122.8,     123091.0,       125250.3]
    else  error("Data not available for Z = $Z")
    end
    return( wa )
end

"""
`PeriodicTable.fluorescenceYields_KrauseAdopted2016(Z::Int64)`
    ... to provide the K-shell fluorescence yield omega_K of the neutral element with nuclear charge Z; this yield is
        the probability that a K-shell vacancy is filled radiatively rather than by an Auger transition, so that the
        associated Auger yield is simply a_K = 1 - omega_K. The yields are dimensionless and are tabulated for
        Z = 3, ..., 100; a value of -1. is returned for Z = 1, 2, where a K-shell vacancy cannot decay radiatively.
        An omega::Float64 is returned.

        Provenance (please read before relying on these numbers):
        The values are those adopted by the xraylib project [data/fluor_yield.dat, as of the K-shell revision of
        11-Apr-2016], which JAC follows here. Their history is layered and not fully documented upstream:
        + The base data of that file are due to M. Krause, C. Nestor Jr., C. Sparks Jr. and E. Ricci, "X-ray
          fluorescence cross sections for K and L X-rays of the elements", ORNL-5399, Oak Ridge National
          Laboratory (1978); this is the compilation behind M. O. Krause, J. Phys. Chem. Ref. Data 8, 307 (1979).
        + The K-shell values of that base data set were *replaced* in 2016 by the present ones. The xraylib
          documentation names newer sources (Campbell 2009; Ayri et al. 2021; Kaur et al. 2021) only for the
          L subshells, so that the origin of the revised K-shell values is not documented upstream.
        + The two sets differ appreciably: by factors 2 - 20 for Z <= 6 (Be: 3.3e-5 -> 6.9e-4), and in their trend
          for heavy elements, where the ORNL-5399 values peak near Z = 74 (omega_K = 0.983) and fall off to 0.969
          at Z = 92, whereas the adopted values rise monotonically to 0.972. The adopted values were preferred
          because they are those in current use and because the ORNL-5399 set places Be below Li (omega_K decreasing
          from Z = 3 to Z = 4), which is unphysical.
        These yields are semi-empirical and should be understood as good to a few percent at best; cf. also the
        review of Hubbell et al., J. Phys. Chem. Ref. Data 23, 339 (1994).
"""
function fluorescenceYields_KrauseAdopted2016(Z::Int64)
    # K-shell fluorescence yields omega_K [dimensionless] for Z = 3, ..., 100; -1. if no value is defined.
    wa = [
         -1.0,     -1.0,0.0002928,0.0006929, 0.001409, 0.002575, 0.004349, 0.006909,  0.01045,  0.01519,   ##  Z =  1 -  10
      0.02133,  0.02911,  0.03872,  0.05037,  0.06422,  0.08038,  0.09892,   0.1199,   0.1432,   0.1687,   ##  Z = 11 -  20
       0.1962,   0.2256,   0.2564,   0.2885,   0.3213,   0.3546,    0.388,   0.4212,   0.4538,   0.4857,   ##  Z = 21 -  30
       0.5166,   0.5464,   0.5748,   0.6019,   0.6275,   0.6517,   0.6744,   0.6956,   0.7155,    0.734,   ##  Z = 31 -  40
       0.7512,   0.7672,   0.7821,   0.7958,   0.8086,   0.8204,   0.8313,   0.8415,   0.8508,   0.8595,   ##  Z = 41 -  50
       0.8676,    0.875,   0.8819,   0.8883,   0.8942,   0.8997,   0.9049,   0.9096,    0.914,   0.9181,   ##  Z = 51 -  60
        0.922,   0.9255,   0.9289,    0.932,   0.9349,   0.9376,   0.9401,   0.9425,   0.9447,   0.9467,   ##  Z = 61 -  70
       0.9487,   0.9505,   0.9522,   0.9538,   0.9553,   0.9567,    0.958,   0.9592,   0.9604,   0.9615,   ##  Z = 71 -  80
       0.9625,   0.9634,   0.9643,   0.9652,   0.9659,   0.9667,   0.9674,    0.968,   0.9686,   0.9691,   ##  Z = 81 -  90
       0.9696,   0.9701,   0.9706,    0.971,   0.9713,   0.9717,    0.972,   0.9722,   0.9725,   0.9727]   ##  Z = 91 - 100

    if  Z < 1  ||  Z > 100   error("No K-shell fluorescence yield available for Z = $Z; the data cover Z = 3, ..., 100.")   end
    wb = wa[Z]
    if  wb == -1.   error("No K-shell fluorescence yield defined for Z = $Z; a K-shell vacancy cannot decay radiatively.")  end

    return( wb )
end


"""
`PeriodicTable.dielectronicRecombinationParameters_Arnaud1985(Z::Int64)`
    ... to return the H-like-sequence dielectronic recombination fit parameters (ADI, T0) of M. Arnaud &
        R. Rothenflug, Astron. Astrophys. Suppl. Ser. 60, 425 (1985), their Table IIB, for the bare, recombining
        nucleus of nuclear charge Z; the rate coefficient follows their Eq. (3) with BDI = 0 (as stated for all
        4 tabulated ions, so the {1 + BDI exp(-T1/T)} bracket reduces to 1):
            alpha_d(T) [cm^3/s] = ADI * T^(-3/2) * exp(-T0/T),   T, T0 in K, ADI in cm^3/s K^(3/2).
        Available for only 4 ions -- O VIII (Z=8), Mg XII (Z=12), Ca XX (Z=20), Fe XXVI (Z=26) -- since the
        paper tabulates specific ions rather than a general Z-scaling formula for this sequence. A
        (ADI, T0)::NTuple{2,Float64} is returned.

        Note on transcription: Table IIB's column headers carry small unit-exponent annotations that are easy to
        misread from a scanned 1985 paper. The values below were cross-checked by requiring T0 to scale
        consistently with each ion's H-like ionization potential (Table I of the same paper): T0/IP works out to
        55-60% across all 4 ions, a physically sensible and internally consistent ratio, which supports this
        reading -- but the values should still be treated as best-effort and spot-checked against an independent
        source before being relied on for precision work.
"""
function dielectronicRecombinationParameters_Arnaud1985(Z::Int64)
    #                Z    ADI [cm^3/s K^1.5]   T0 [K]
    if     Z ==  8   wa = (4.72e-2,            6.00e6)     ## O VIII
    elseif Z == 12   wa = (9.28e-2,            1.45e7)     ## Mg XII
    elseif Z == 20   wa = (2.68e-1,            3.74e7)     ## Ca XX
    elseif Z == 26   wa = (3.69e-1,            6.00e7)     ## Fe XXVI
    else   error("No H-like Arnaud (1985) dielectronic recombination parameters available for Z = $Z; " *
                 "only O VIII (Z=8), Mg XII (Z=12), Ca XX (Z=20), Fe XXVI (Z=26) are tabulated.")
    end

    return( wa )
end


end # module

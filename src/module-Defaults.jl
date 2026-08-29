
"""
`module  JAC.Defaults`  
... a submodel of JAC that contains all constants and global parameters of the program.
"""
module Defaults


using Dates,  JenaAtomicCalculator, ..Basics,  ..Radial,  ..Math
# 2014 CODATA recommended values, obtained from http://physics.nist.gov/cuu/Constants/

export  convertUnits, getDefaults,  setDefaults

# Dimensionless constants
# FINE_STRUCTURE_CONSTANT and INVERSE_FINE_STRUCTURE_CONSTANT are deliberately NOT const: both are
# overwritten together by setDefaults("fine structure constant: alpha", value) to support alpha-variation
# (q-factor) computations, which need to re-run the SCF/CI machinery at a slightly shifted alpha. Every
# other module reads these two globals by name (Defaults.FINE_STRUCTURE_CONSTANT / .INVERSE_...), so an
# update here propagates automatically, without touching those modules.
FINE_STRUCTURE_CONSTANT         = 7.297_352_566_4e-3
INVERSE_FINE_STRUCTURE_CONSTANT = 137.035_999_139

# Constants in SI units
const BOHR_RADIUS_SI                  = 0.529_177_210_67e-10
const BOLTZMANN_CONSTANT_SI           = 1.380_648_52e-23
const ELECTRON_MA_SI                  = 9.109_383_56e-31
## The muon is 206.77 times heavier than the electron and is otherwise the same particle: same charge, same
## spin, same Dirac equation.  In atomic units the electron mass is 1, so this ratio IS the muon mass.
const MUON_ELECTRON_MASS_RATIO        = 206.768_2830
const ELEMENTARY_CHARGE_SI            = 1.602_176_620_8e-19
const HARTREE_ENERGY_SI               = 4.359_744_650e-18
const PLANCK_CONSTANT_SI              = 6.626_070_040e-34
const PLANCK_CONSTANT_OVER_2_PI_SI    = 1.054_571_800e-34
const SPEED_OF_LIGHT_IN_VACUUM_SI     = 299_792_458

# Constants with eV instead of J
const BOLTZMANN_CONSTANT_EV           = 8.617_330_3e-5
const HARTREE_ENERGY_EV               = 27.211_386_02
const PLANCK_CONSTANT_EV              = 4.135_667_662e-15
const PLANCK_CONSTANT_OVER_2_PI_EV    = 6.582_119_514e-16

# Constants with other atomic units
const RYDBERG_IN_KAYSERS              = 109.73731534e3
const ELECTRON_MASS_IN_G              = 9.1093897e-28
const HBAR_IN_ERGS                    = 1.05457266e-27
const ELECTRON_CHARGE_IN_ESU          = 4.80320680e-10

# Constants in unified atomic mass units (u)
const ELECTRON_MASS_U                 = 5.485_799_090_70e-4
const NEUTRON_MASS_U                  = 1.008_664_915_88
const PROTON_MASS_U                   = 1.007_276_466_879

# Relationships between energy equivalents
const ELECTRON_VOLT_ATOMIC_MASS_UNIT_RELATIONSHIP = 1.073_544_110_5e-9

# Predefined conversion factors
const CONVERT_ENERGY_AU_TO_EV           = HARTREE_ENERGY_EV
const CONVERT_ENERGY_AU_TO_KAYSERS      = 2.0 * RYDBERG_IN_KAYSERS
const CONVERT_ENERGY_AU_TO_PER_SEC      = 6.57968974479e15
const CONVERT_INTENSITY_AU_TO_W_CM2     = 3.50944758e16
const CONVERT_TIME_AU_TO_SEC            = 2.418_884_254e-17
const CONVERT_CROSS_SECTION_AU_TO_BARN  = BOHR_RADIUS_SI^2 * 1.0e28
const CONVERT_RATE_AU_TO_PER_SEC        = (ELECTRON_MASS_IN_G/HBAR_IN_ERGS) * ((ELECTRON_CHARGE_IN_ESU^2/HBAR_IN_ERGS)^2)
const CONVERT_STRENGTH_AU_TO_BARN_EV    = CONVERT_CROSS_SECTION_AU_TO_BARN * HARTREE_ENERGY_EV
const CONVERT_STRENGTH_AU_TO_CM2_EV     = BOHR_RADIUS_SI^2 * 1.0e4 * HARTREE_ENERGY_EV
const CONVERT_LENGTH_AU_TO_FEMTOMETER   = BOHR_RADIUS_SI * 1.0e15

# Predefined coefficients for numerical integration
const FINITE_DIFFERENCE_NPOINTS         = 6
const fivePointCoefficients             = Array{Float64}([2 * 7, 32, 12, 32, 7]) * 2 / 45  # * grid.h     
const newtonCotesCoefficients           = fivePointCoefficients

n = FINITE_DIFFERENCE_NPOINTS     
    const weights = Array{Float64}(undef, 2*n + 1, 2*n + 1)
for i = 1:2*n + 1       
    weights[i,:] = Math.finiteDifferenceWeights(-n + i - 1, 2*n + 1, order=1)[2,:]     
end

#
# Global settings that can be (re-) defined by the user.
GBL_FRAMEWORK                = "relativistic"
GBL_CONT_POTENTIAL           = Basics.DFSField(0.42)
GBL_CONT_SOLUTION            = BsplineGalerkin()         ###  ContBessel(), ContSine(), AsymptoticCoulomb(), NonrelativisticCoulomb(), BsplineGalerkin()
GBL_CONT_NORMALIZATION       = AlokNorm()                ###  PureSineNorm(), CoulombSineNorm(), OngRussekNorm(), AlokNorm()
GBL_QED_HYDROGENIC_LAMBDAC   = [1.0,  1.0,  1.0,  1.0,  1.0]
GBL_QED_NUCLEAR_CHARGE       = 0.1
GBL_WARNINGS                 = String[]
## The warnings are collected so that a long job can be inspected AFTERWARDS -- anything printed to stdout
## during a ten-hour cascade has scrolled away long before it ends.  Two properties are needed for that to
## work (12-Aug-2026):  a message repeated by every line must appear ONCE, with a count, or the report is
## unreadable; and the collector must survive Threads.@threads, since push! onto a shared Vector can lose
## entries and can crash on a reallocation -- and the line loops of DielectronicRecombination and
## ImpactExcitation are already threaded.
GBL_WARNINGS_COUNT           = Dict{String,Int64}()
const GBL_WARNINGS_LOCK      = ReentrantLock()

GBL_ENERGY_UNIT              = "eV"
GBL_CROSS_SECTION_UNIT       = "barn"
GBL_EDIFF_CROSS_SECTION_UNIT = "barn/eV"
GBL_RATE_UNIT                = "1/s"
GBL_STRENGTH_UNIT            = "cm^2 eV"
GBL_TIME_UNIT                = "sec"

GBL_NUCLEAR_CHARGE           = 1.0
GBL_NUCLEAR_MODEL            = missing

GBL_PRINT_DATA               = false
GBL_PRINT_SUMMARY            = false
GBL_PRINT_TEST               = false
GBL_PRINT_DEBUG              = false

GBL_DATA_IOSTREAM            = nothing
GBL_DATAX_IOSTREAM           = nothing
GBL_DATAY_IOSTREAM           = nothing
GBL_SUMMARY_IOSTREAM         = nothing

GBL_STANDARD_GRID            = Radial.Grid(true, printout=false)        



"""
`Defaults.convertUnits()`  
    ... converts some data from one format/unit into another one; cf. the supported keystrings and return values.

+ `("cross section: from atomic to predefined unit", value::Float64)`  or  `("cross section: from atomic", value::Float64)`  
    ... to convert an cross section value from atomic to the predefined cross section unit; a Float64 is returned.

+ `("cross section: from barn to atomic unit", value::Float64)` 
    ... to convert an cross section value from barn atomic section unit; a Float64 is returned.

+ `("cross section: from atomic to barn", value::Float64)`  or  `("cross section: from atomic to Mbarn", value::Float64)`  or
    `("cross section: from atomic to cm^2", value::Float64)` 
    ... to convert an energy value from atomic to the speficied cross section unit; a Float64 is returned.

+ `("cross section: from predefined to atomic unit", value::Float64)`  or  `("cross section: to atomic", value::Float64)`
    ... to convert a cross section value from the predefined to the atomic cross section unit; a Float64 is returned.

+ `("einstein B: from atomic", value::Float64)`  
    ... to convert a Einstein B coefficient from atomic to the speficied energy units; a Float64 is returned.

+ `("density: from [g/cm^3] to atomic", value::Float64)`  
    ... to convert a mass density from [g/cm^3] to atomic units [u/a_o^3]; a Float64 is returned.

+ `("energy-diff. cross section: from atomic to predefined unit", value::Float64)`  or  
    `("energy-diff. cross section: from atomic", value::Float64)`  
    ... to convert an energy-diff. cross section value from atomic to the predefined energy-diff. cross section unit; 
        a Float64 is returned.

+ `("energy: from atomic to eV", value::Float64)`  or  `("energy: from atomic to Kayser", value::Float64)`    or
    `("energy: from atomic to Hz", value::Float64)`  or  `("energy: from atomic to Angstrom", value::Float64)`  or
    `("energy: from atomic to Ws", value::Float64)`
    ... to convert an energy value from atomic to the speficied energy unit; a Float64 is returned.

+ `("energy: from predefined to atomic unit", value::Float64)`  or  `("energy: to atomic", value::Float64)`... to convert an energy value 
                                                from the predefined to the atomic energy unit; a Float64 is returned.
+ `("energy: from eV to atomic", value::Float64)` ... to convert an energy value from eV to the atomic energy unit; a Float64 is returned.

+ `("energy: from wavelength [nm] to atomic", value::Float64)` ... to convert a wavelength [nm] to the atomic energy unit; a Float64 is returned.

+ `("intensity: from W/cm^2 to atomic", value::Float64)` ... to convert the intensity [in W/cm^2] to the atomic intensity unit; a Float64 is returned.

+ `("kinetic energy to wave number: atomic units", value::Float64)`  ... to convert a kinetic energy value (in a.u.) into a wave number
                                                k (a.u.); a Float64 is returned.

+ `("kinetic energy to wavelength: atomic units", value::Float64)`  ... to convert a kinetic energy value (in a.u.) into a wavelength (a.u.); 
                                                a Float64 is returned.

+ `("length: from fm to atomic", value::Float64)`  ... to convert a length value (in fm) into a.u.;  a Float64 is returned.
+ `("length: from atomic to fm", value::Float64)`  or  `("length: from atomic to cm", value::Float64)`  
    ... to convert a length value (in Bohr's a.u.) to the speficied length unit;  a Float64 is returned.

+ `("moment: from nuclear magneton to atomic", value::Float64)`  ... to convert a value in mu_nuc into a.u.;  a Float64 is returned.
+ `("moment: from nuclear magneton x fm^2 to atomic", value::Float64)`  ... to convert a value in [mu_nuc x fm^2] into a.u.;  a Float64 is returned.

+ `("rate: from atomic to predefined unit", value::Float64)`  or  ("rate: from atomic", value::Float64)  ... to convert a rate value 
                                            from atomic to the predefined rate unit; a Float64 is returned.

+ `("rate: from atomic to 1/s", value::Float64)`  ... to convert an rate value from atomic to the speficied rate unit; a Float64 is returned.

+ `("rate: from predefined to atomic unit", value::Float64)`  or  `("rate: to atomic", value::Float64)'... to convert a
                                                rate value from the predefined to the atomic rate unit; a Float64 is returned.


+ `("strength: from atomic to predefined unit", value::Float64)`  or  ("strength: from atomic", value::Float64)  ... to convert a (resonance) 
                                                strength value from atomic to the predefined rate unit; a Float64 is returned.
                                                
+ `("time: from atomic to predefined unit", value::Float64)`  or  ("time: from atomic", value::Float64)  ... to convert an time value 
                                            from atomic to the predefined time unit; a Float64 is returned.

+ `("time: from atomic to sec", value::Float64)`  or  `("time: from atomic to fs", value::Float64)'  or
    `("time: from atomic to as", value::Float64)`  ... to convert a time value from atomic to the speficied time unit; a Float64 is returned.

+ `("time: from predefined to atomic unit", value::Float64)`  or  `("time: to atomic", value::Float64)'  ... to convert a
                                            time value from the predefined to the atomic time unit; a Float64 is returned.

+ `("temperature: from Kelvin to (Hartree) units", value::Float64)`  ... to convert a temperature in Kelvin into atomic (Hartree) units; 
                                                a Float64 is returned.

+ `("temperature: from atomic to Kelvin", value::Float64)`  ... to convert an atomic (energy) unit into Kelvin; 1 Hartree = 315774.64 K; 
                                                a Float64 is returned.

+ `("wave number to total electron energy: atomic units", value::Float64)`  ... to convert a wavenumber (a.u.) into the total electron
                                            energy, including the rest energy; a Float64 is returned.

+ `("wave number to kinetic energy: atomic units", value::Float64)`  ... to convert a wavenumber (a.u.) into the kinetic energy; 
                                            a Float64 is returned.
"""
function convertUnits(sa::String, wa::Float64)
    global  CONVERT_ENERGY_AU_TO_EV, CONVERT_ENERGY_AU_TO_KAYSERS, CONVERT_ENERGY_AU_TO_PER_SEC, CONVERT_TIME_AU_TO_SEC,
            CONVERT_CROSS_SECTION_AU_TO_BARN,  CONVERT_RATE_AU_TO_PER_SEC,  CONVERT_LENGTH_AU_TO_FEMTOMETER 

    if       sa in ["cross section: from atomic to predefined unit", "cross section: from atomic"]
        if      Defaults.getDefaults("unit: cross section") == "a.u."   return( wa )
        elseif  Defaults.getDefaults("unit: cross section") == "barn"   return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN )
        elseif  Defaults.getDefaults("unit: cross section") == "Mbarn"  return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN * 1.0e-6 )
        elseif  Defaults.getDefaults("unit: cross section") == "cm^2"   return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN * 1.0e-24 )
        else    error("Defaults.convertUnits(): the current setting `unit: cross section` is currently $(Defaults.getDefaults("unit: cross section")), which this conversion does not handle. Valid: a.u., barn, Mbarn, cm^2.")
        end
    
    elseif   sa in ["cross section: from atomic to barn"]               return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN )
    elseif   sa in ["cross section: from atomic to Mbarn"]              return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN * 1.0e-6 )
    elseif   sa in ["cross section: from atomic to cm^2"]               return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN * 1.0e-24 )

    elseif   sa in ["cross section: from predefined to atomic unit", "cross section: to atomic"]
        if      Defaults.getDefaults("unit: cross section") == "a.u."   return( wa )
        elseif  Defaults.getDefaults("unit: cross section") == "barn"   return( wa / CONVERT_CROSS_SECTION_AU_TO_BARN )
        elseif  Defaults.getDefaults("unit: cross section") == "Mbarn"  return( wa / CONVERT_CROSS_SECTION_AU_TO_BARN * 1.0e6 )
        elseif  Defaults.getDefaults("unit: cross section") == "cm^2"   return( wa / CONVERT_CROSS_SECTION_AU_TO_BARN * 1.0e24 )
        else    error("Defaults.convertUnits(): the current setting `unit: cross section` is currently $(Defaults.getDefaults("unit: cross section")), which this conversion does not handle. Valid: a.u., barn, Mbarn, cm^2.")
        end

    elseif   sa in ["cross section: from barn to atomic unit"]          return( wa / CONVERT_CROSS_SECTION_AU_TO_BARN )

    elseif   sa in ["density: from [g/cm^3] to atomic"]                 wg = ELECTRON_MASS_U / ELECTRON_MASS_IN_G
        wcm = 1 / (100 * BOHR_RADIUS_SI);                               return( wa * wg / wcm^3 )  # 

    elseif   sa in ["Einstein B: from atomic"]                          return( wa )  # Einstein B not yet properly converted.
    
    elseif   sa in ["energy-diff. cross section: from atomic to predefined unit", "energy-diff. cross section: from atomic"]
        if      Defaults.getDefaults("unit: energy-diff. cross section") == "a.u."      return( wa )
        elseif  Defaults.getDefaults("unit: energy-diff. cross section") == "barn/eV"   
                                                                        return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN/CONVERT_ENERGY_AU_TO_EV )
        else    error("Defaults.convertUnits(): the current setting `unit: energy-diff. cross section` is currently $(Defaults.getDefaults("unit: energy-diff. cross section")), which this conversion does not handle. Valid: a.u., barn/eV.")
        end

    elseif   sa in ["energy: from atomic to predefined unit", "energy: from atomic"]
        if      Defaults.getDefaults("unit: energy") == "eV"            return( wa * CONVERT_ENERGY_AU_TO_EV )
        elseif  Defaults.getDefaults("unit: energy") == "Kayser"        return( wa * CONVERT_ENERGY_AU_TO_KAYSERS )
        elseif  Defaults.getDefaults("unit: energy") == "Hartree"       return( wa )
        elseif  Defaults.getDefaults("unit: energy") == "Hz"            return( wa * CONVERT_ENERGY_AU_TO_PER_SEC )
        elseif  Defaults.getDefaults("unit: energy") == "A"             return( 2pi * Defaults.INVERSE_FINE_STRUCTURE_CONSTANT / wa * 
                                                                                        Defaults.BOHR_RADIUS_SI * 1.0e10 )
        else    error("Defaults.convertUnits(): the current setting `unit: energy` is currently $(Defaults.getDefaults("unit: energy")), which this conversion does not handle. Valid: eV, Kayser, Hartree, Hz, A.")
        end
    
    elseif   sa in ["energy: from atomic to eV"]                        return( wa * CONVERT_ENERGY_AU_TO_EV )
    elseif   sa in ["energy: from atomic to Kayser"]                    return( wa * CONVERT_ENERGY_AU_TO_KAYSERS )
    elseif   sa in ["energy: from atomic to Hz"]                        return( wa * CONVERT_ENERGY_AU_TO_PER_SEC )
    elseif   sa in ["energy: from atomic to Angstrom"]                  return( 2pi * Defaults.INVERSE_FINE_STRUCTURE_CONSTANT / wa * 
                                                                                        Defaults.BOHR_RADIUS_SI * 1.0e10 )
    elseif   sa in ["energy: from atomic to Ws"]                        return( wa * 4.35974e-18 )

    elseif   sa in ["energy: from predefined to atomic unit", "energy: to atomic"]
        if      Defaults.getDefaults("unit: energy") == "eV"            return( wa / CONVERT_ENERGY_AU_TO_EV )
        elseif  Defaults.getDefaults("unit: energy") == "Kayser"        return( wa / CONVERT_ENERGY_AU_TO_KAYSERS )
        elseif  Defaults.getDefaults("unit: energy") == "Hartree"       return( wa )
        elseif  Defaults.getDefaults("unit: energy") == "Hz"            return( wa / CONVERT_ENERGY_AU_TO_PER_SEC )
        elseif  Defaults.getDefaults("unit: energy") == "A"             return( 2pi * Defaults.INVERSE_FINE_STRUCTURE_CONSTANT / wa * 
                                                                                        Defaults.BOHR_RADIUS_SI * 1.0e10 )
        else    error("Defaults.convertUnits(): the current setting `unit: energy` is currently $(Defaults.getDefaults("unit: energy")), which this conversion does not handle. Valid: eV, Kayser, Hartree, Hz, A.")
        end

    elseif   sa in ["energy: from eV to atomic"]                        return( wa / CONVERT_ENERGY_AU_TO_EV )
    elseif   sa in ["energy: from wavelength [nm] to atomic"]           return( 1.0e7 / (CONVERT_ENERGY_AU_TO_KAYSERS * wa) )

    elseif   sa in ["intensity: from W/cm^2 to atomic"]                 return( wa / CONVERT_INTENSITY_AU_TO_W_CM2 )
    elseif   sa in ["intensity: from atomic to W/cm^2"]                 return( wa * CONVERT_INTENSITY_AU_TO_W_CM2 )

    elseif  sa in ["kinetic energy to wave number: atomic units"]
        c = INVERSE_FINE_STRUCTURE_CONSTANT
        wb = sqrt( wa*wa/(c*c) + 2wa );                                 return( wb )

    elseif  sa in ["kinetic energy to wavelength: atomic units"]
        c = INVERSE_FINE_STRUCTURE_CONSTANT
        wb = sqrt( wa*wa/(c*c) + 2wa );                                 return( 2pi / wb )

    elseif  sa in ["length: from fm to atomic"]                         return (wa / CONVERT_LENGTH_AU_TO_FEMTOMETER )
    elseif  sa in ["length: from atomic to fm"]                         return (wa * CONVERT_LENGTH_AU_TO_FEMTOMETER )
    elseif  sa in ["length: from atomic to cm"]                         return (wa * CONVERT_LENGTH_AU_TO_FEMTOMETER * 1.0e-13 )

    elseif  sa in ["moment: from nuclear magneton to atomic"]           return (wa * 5.446170e-4 / 2. )
    elseif  sa in ["moment: from nuclear magneton x fm^2 to atomic"]    return (wa * 1.944690e-13 / 2. )

    elseif    sa in ["rate: from atomic to predefined unit", "rate: from atomic"]
        if       Defaults.getDefaults("unit: rate") == "1/s"            return( wa * CONVERT_RATE_AU_TO_PER_SEC )
        elseif   Defaults.getDefaults("unit: rate") == "a.u."           return( wa )
        else     error("Defaults.convertUnits(): the current setting `unit: rate` is currently $(Defaults.getDefaults("unit: rate")), which this conversion does not handle. Valid: 1/s, a.u..")
        end
    
    elseif   sa in ["rate: from atomic to 1/s"]                         return( wa * CONVERT_RATE_AU_TO_PER_SEC )

    elseif   sa in ["rate: from predefined to atomic unit", "rate: to atomic"]
        # "unit: rate", NOT "rate: time".  Until 29-Aug-2026 both lines read "rate: time", which is not a
        # settings key at all, so getDefaults raised "Unsupported keystring:: rate: time" BEFORE any conversion
        # was attempted -- convertUnits("rate: to atomic", ...) could never work, while the from-atomic
        # direction six lines above, which reads "unit: rate", always did.
        if      Defaults.getDefaults("unit: rate") == "1/s"             return( wa / CONVERT_RATE_AU_TO_PER_SEC )
        elseif  Defaults.getDefaults("unit: rate") == "a.u."            return( wa  )
        else    error("Defaults.convertUnits(): the current setting `unit: rate` is currently $(Defaults.getDefaults("unit: rate")), which this conversion does not handle. Valid: 1/s, a.u..")
        end

    elseif    sa in ["strength: from atomic to predefined unit", "strength: from atomic"]
        if       Defaults.getDefaults("unit: strength") == "barn eV"    return( wa * CONVERT_STRENGTH_AU_TO_BARN_EV )
        elseif   Defaults.getDefaults("unit: strength") == "cm^2 eV"    return( wa * CONVERT_STRENGTH_AU_TO_CM2_EV )
        elseif   Defaults.getDefaults("unit: strength") == "a.u."       return( wa )
        else     error("Defaults.convertUnits(): the current setting `unit: strength` is currently $(Defaults.getDefaults("unit: strength")), which this conversion does not handle. Valid: barn eV, cm^2 eV, a.u..")
        end
        
    elseif  sa in ["time: from atomic to predefined unit", "time: from atomic"]
        if      Defaults.getDefaults("unit: time") == "sec"             return( wa * CONVERT_TIME_AU_TO_SEC )
        elseif  Defaults.getDefaults("unit: time") == "fs"              return( wa * CONVERT_TIME_AU_TO_SEC * 1.0e15 )
        elseif  Defaults.getDefaults("unit: time") == "as"              return( wa * CONVERT_TIME_AU_TO_SEC * 1.0e18 )
        elseif  Defaults.getDefaults("unit: time") == "a.u."            return( wa )
        else    error("Defaults.convertUnits(): the current setting `unit: time` is currently $(Defaults.getDefaults("unit: time")), which this conversion does not handle. Valid: sec, fs, as, a.u..")
        end
    
    elseif   sa in ["time: from atomic to sec"]                         return( wa * CONVERT_TIME_AU_TO_SEC )
    elseif   sa in ["time: from atomic to fs"]                          return( wa * CONVERT_TIME_AU_TO_SEC * 1.0e15 )
    elseif   sa in ["time: from atomic to as"]                          return( wa * CONVERT_TIME_AU_TO_SEC * 1.0e18 )

    elseif  sa in ["time: from predefined to atomic unit", "time: to atomic"]
        if      Defaults.getDefaults("unit: time") == "sec"             return( wa / CONVERT_TIME_AU_TO_SEC )
        elseif  Defaults.getDefaults("unit: time") == "fs"              return( wa / (CONVERT_TIME_AU_TO_SEC * 1.0e15) )
        elseif  Defaults.getDefaults("unit: time") == "as"              return( wa / (CONVERT_TIME_AU_TO_SEC * 1.0e18) )
        elseif  Defaults.getDefaults("unit: time") == "a.u."            return( wa )
        else    error("Defaults.convertUnits(): the current setting `unit: time` is currently $(Defaults.getDefaults("unit: time")), which this conversion does not handle. Valid: sec, fs, as, a.u..")
        end
    
    elseif   sa in ["temperature: from Kelvin to (Hartree) units"]      return( wa / 315774.64 )
    elseif   sa in ["temperature: from atomic to Kelvin"]               return( wa * 315774.64 )

    elseif  sa in ["wave number to total electron energy: atomic units"]
        c = INVERSE_FINE_STRUCTURE_CONSTANT
        wb = sqrt( wa*wa/(c*c) + 2wa );                                 return( sqrt(wb*wb*c*c + c^4) )

    elseif  sa in ["wave number to kinetic energy: atomic units"]
        c = INVERSE_FINE_STRUCTURE_CONSTANT
        wb = sqrt( wa*wa/(c*c) + 2wa );                                 return( c + sqrt(wb*wb*c*c + c^2) )

    else
        error("Unsupported keystring = $sa")
    end
end


"""


+ `(string, values::Array{Float64,1})`  
    ... to convert for the same strings as above but for an list of values; a corresponding Array{Float64,1} is returned.
"""
function convertUnits(sa::String, values::Array{Float64,1})
    newValues = Float64[]
    for  wa  in  values   wb = convertUnits(sa, wa);     push!(newValues, wb)    end
    
    return( newValues )
end



"""
`Defaults.setDefaults()`  
    ... (re-) defines some 'standard' settings which are common to all the computations with the JAC module, and which can 
        be 'overwritten' by the user. --- An improper setting of some variable may lead to an error message, if recognized
        immediately. The following defaults apply if not specified otherwise by the user: the framework is 'relativistic', 
        energies are given in eV and cross sections in barn. Note that, internally, atomic units are used throughout for 
        all the computations within the program. nothing is returned if not indicated otherwise.

## + `("framework: relativistic")`  or  `("framework: non-relativistic")`   
##     ... to define a relativistic or non-relativistic framework for all subsequent computations. 

+ `("method: continuum, spherical Bessel")`  or  `("method: continuum, pure sine")`  or  
    `("method: continuum, asymptotic Coulomb")`  or  `("method: continuum, nonrelativistic Coulomb")`  or  
    `("method: continuum, Galerkin")`  
    ... to define a a method for the generation of the continuum orbitals as (pure) spherical Bessel, pure sine,
        asymptotic Coulomb, nonrelativistic Coulomb orbital or by means of the B-spline-Galerkin method.

+ `("method: normalization, pure sine")`  or  `("method: normalization, pure Coulomb")`  or  `("method: normalization, Ong-Russek")`   
    ... to define a method for the normalization of the continuum orbitals as asymptotically (pure) sine or Coulomb 
        functions, or following the procedure by Ong & Russek (1978).

+ `("method: Wigner symbols, exact")`  or  `("method: Wigner symbols, floating-point")`  or
    `("method: Wigner symbols, GSL")`
    ... to (pre-) define HOW the Wigner 3-j, 6-j and 9-j symbols are evaluated. All three agree to about one part in
        1e15, so this is a choice of cost and not of physics: `floating-point` evaluates in Float64, is 4.2x faster on a
        mid-size cascade with BIT-IDENTICAL results, and is the DEFAULT; `exact` evaluates in Rational{BigInt} and is
        what JAC used until 25-Aug-2026; and `GSL` is faster again but differs in the ninth digit, which would put
        approved reference data in question. See `AngularMomentum.AbstractWignerMethod`.

+ `("nuclear: charge", Z::Float64)`  or  `("nuclear: model", nm::Any")`
    ... to define the nuclear charge or nuclear model for the pedestrian approach to atomic computations.

+ `("fine structure constant: alpha", value::Float64)`
    ... to (re-) define the fine-structure constant alpha used throughout JAC, together with its reciprocal
        (the speed of light in atomic units); needed for alpha-variation (q-factor) computations, which
        re-run the SCF/CI machinery at a slightly shifted alpha. Defaults.FINE_STRUCTURE_CONSTANT and
        Defaults.INVERSE_FINE_STRUCTURE_CONSTANT are updated together so they never go out of sync.

+ `("QED model: Petersburg")`  or  `("QED model: Sydney")`   
    ... to define a model for the computation of the QED corrections following the work by Shabaev et al. (2011; Petersburg) 
        or Flambaum and Ginges (2004; Syney).

+ `("unit: energy", "eV")`  or  `("unit: energy", "Kayser")`  or  `("unit: energy", "Hartree")`  or  
    `("unit: energy", "Hz")`  or  `("unit: energy", "Hz")`  
    ... to (pre-) define the energy units for all further printouts and communications with the JAC module.

+ `("unit: cross section", "a.u.")`  or  `("unit: cross section", "barn")`  or  `("unit: cross section", "Mbarn")`  
    ... to (pre-) define the unit for the printout of cross sections.

+ `("unit: rate", "a.u.")`  or  `("unit: rate", "1/s")`  ... to (pre-) define the unit for the printout of rates.

+ `("unit: resonance strength", "a.u.")`  or  `("unit: resonance strength", "barn eV")`  or  
    `("unit: resonance strength", "cm^2 eV")`  ... to (pre-) define the unit for the printout of resonance strengths.

+ `("unit: time", "a.u.")`  or  `("unit: time", "sec")`  or  `("unit: time", "fs")`  or  `("unit: time", "as")`  
    ... to (pre-) define the unit for the printout and communications of times with the JAC module.
"""
function setDefaults(sa::String)
    global GBL_FRAMEWORK, GBL_CONT_SOLUTION, GBL_CONT_NORMALIZATION  

    if        sa == "framework: relativistic"                            GBL_FRAMEWORK           = "relativistic"    
    elseif    sa == "framework: non-relativistic"                        GBL_FRAMEWORK           = "non-relativistic"  
    elseif    sa == "method: continuum, spherical Bessel"                GBL_CONT_SOLUTION       = ContBessel()  
    elseif    sa == "method: continuum, pure sine"                       GBL_CONT_SOLUTION       = ContSine()  
    elseif    sa == "method: continuum, asymptotic Coulomb"              GBL_CONT_SOLUTION       = AsymptoticCoulomb()  
    elseif    sa == "method: continuum, nonrelativistic Coulomb"         GBL_CONT_SOLUTION       = NonrelativisticCoulomb()  
    elseif    sa == "method: continuum, Galerkin"                        GBL_CONT_SOLUTION       = BsplineGalerkin()   
    elseif    sa == "method: normalization, pure sine"                   GBL_CONT_NORMALIZATION  = PureSineNorm()   
    elseif    sa == "method: normalization, pure Coulomb"                GBL_CONT_NORMALIZATION  = CoulombSineNorm()  
    elseif    sa == "method: normalization, Ong-Russek"                  GBL_CONT_NORMALIZATION  = OngRussekNorm()
    elseif    sa == "method: normalization, Alok"                        GBL_CONT_NORMALIZATION  = AlokNorm()
    elseif    sa == "method: Wigner symbols, exact"                      selectWignerMethod(:exact)
    elseif    sa == "method: Wigner symbols, floating-point"             selectWignerMethod(:floating)
    elseif    sa == "method: Wigner symbols, GSL"                        selectWignerMethod(:gsl)
    else      error("Unsupported keystring:: $sa")
    end
    nothing
end


"""
`Defaults.selectWignerMethod(name::Symbol)`
    ... to switch how the Wigner symbols are evaluated, from the short name used by the `setDefaults` keystrings:
        `:exact`, `:floating` or `:gsl`. The choice itself is DISPATCHED on a singleton type inside
        `AngularMomentum`; this routine only turns a keystring into that type. Nothing is returned.
"""
function selectWignerMethod(name::Symbol)
    AM = JenaAtomicCalculator.AngularMomentum
    if        name == :exact      AM.setWignerMethod(AM.ExactWigner())
    elseif    name == :floating   AM.setWignerMethod(AM.FloatingWigner())
    elseif    name == :gsl        AM.setWignerMethod(AM.GslWigner())
    else      error("Unsupported Wigner method:: $name")
    end
    nothing
end




function setDefaults(sa::String, sb::String)
    if        sa == "unit: energy"
        units = ["eV", "Kayser", "Hartree", "Hz", "A"]
        !(sb in units)    &&    error("Currently supported energy units: $(units)")
        global GBL_ENERGY_UNIT = sb
    elseif    sa == "unit: cross section"
        units = ["a.u.", "barn", "Mbarn", "cm^2"]
        !(sb in units)    &&    error("Currently supported cross section units: $(units)")
        global GBL_CROSS_SECTION_UNIT = sb
    elseif    sa == "unit: rate"
        units = ["a.u.", "1/s"]
        !(sb in units)    &&    error("Currently supported rate units: $(units)")
        global GBL_RATE_UNIT = sb
    elseif    sa == "unit: strength"
        units = ["a.u.", "barn eV", "cm^2 eV"]
        !(sb in units)    &&    error("Currently supported resonance strength units: $(units)")
        global GBL_STRENGTH_UNIT = sb
    elseif    sa == "unit: time"
        units = ["a.u.", "sec", "fs", "as"]
        !(sb in units)    &&    error("Currently supported time units: $(units)")
        global GBL_TIME_UNIT = sb
    elseif    sa == "print data: open"
        global GBL_PRINT_DATA       = true
        global GBL_DATA_IOSTREAM    = open(sb, "w") 
        println(GBL_DATA_IOSTREAM, "Data file opened at $( string(now())[1:16] ): \n" *
                                    "=====================================   \n")
    elseif    sa == "print data: close"
        global GBL_PRINT_DATA       = false
        close(GBL_DATA_IOSTREAM)
    elseif    sa == "print data-X: open"
        global GBL_PRINT_DATAX       = true
        global GBL_DATAX_IOSTREAM    = open(sb, "w") 
        println(GBL_DATAX_IOSTREAM, "Data file opened at $( string(now())[1:16] ): \n" *
                                    "=====================================   \n")
    elseif    sa == "print data-X: close"
        global GBL_PRINT_DATAX      = false
        close(GBL_DATAX_IOSTREAM)
    elseif    sa == "print data-Y: open"
        global GBL_PRINT_DATAY      = true
        global GBL_DATAY_IOSTREAM   = open(sb, "w") 
        println(GBL_DATAY_IOSTREAM, "Data file opened at $( string(now())[1:16] ): \n" *
                                    "=====================================   \n")
    elseif    sa == "print data-Y: close"
        global GBL_PRINT_DATAY      = false
        close(GBL_DATAY_IOSTREAM)
    elseif    sa == "print summary: open"
        global GBL_PRINT_SUMMARY    = true
        global GBL_SUMMARY_IOSTREAM = open(sb, "w") 
        println(GBL_SUMMARY_IOSTREAM, "Summary file opened at $( string(now())[1:16] ): \n" *
                                    "========================================   \n")
    elseif    sa == "print summary: append"
        global GBL_PRINT_SUMMARY    = true
        global GBL_SUMMARY_IOSTREAM = open(sb, "a") 
        println(GBL_SUMMARY_IOSTREAM, "Summary file re-opened (to append) at $( string(now())[1:16]): \n" *
                                    "======================================================= \n")
    elseif    sa == "print summary: close"
        global GBL_PRINT_SUMMARY    = false
        close(GBL_SUMMARY_IOSTREAM)
    elseif    sa == "print test: open"
        global GBL_PRINT_TEST    = true
        JenaAtomicCalculator.JAC_TEST_IOSTREAM = open(sb, "w") 
        println(JenaAtomicCalculator.JAC_TEST_IOSTREAM, "Test report file opened at $( string(now())[1:16] ): \n" *
                                                        "============================================  \n")
    elseif    sa == "print test: close"
        global GBL_PRINT_TEST    = false
        close(GBL_TEST_IOSTREAM)
    else      error("Unsupported keystring:: $sa")
    end

    nothing
end


function setDefaults(sa::String, Z::Float64)

    if        sa == "nuclear: charge"      global GBL_NUCLEAR_CHARGE = Z
    elseif    sa == "fine structure constant: alpha"
        global FINE_STRUCTURE_CONSTANT         = Z
        global INVERSE_FINE_STRUCTURE_CONSTANT = 1.0 / Z
    else      error("Unsupported keystring:: $sa")
    end

    nothing
end


function setDefaults(sa::String, nm::Any)

    if        sa == "nuclear: model"      global GBL_NUCLEAR_MODEL = nm
    else      error("Unsupported keystring:: $sa")
    end

    nothing
end


"""
+ `("relativistic subshell list", subshells::Array{Subshell,1}; printout::Bool=true)`  
    ... to (pre-) define internally the standard relativistic subshell list on which the standard order of orbitals is based.
"""
function setDefaults(sa::String, subshells::Array{Subshell,1}; printout::Bool=true)
    if        sa == "relativistic subshell list"
        ## DEPRECATED as a user-facing setting (12-Aug-2026); it delegates to the internal setter below.
        ## See there for why this does not belong in a public configuration API.
        Defaults.setStandardSubshellList(subshells; printout=printout)
    else
        error("Unsupported keystring:: $sa")
    end

    nothing
end


"""
`Defaults.setStandardSubshellList(subshells::Array{Subshell,1}; printout::Bool=true)`
    ... sets the standard relativistic subshell list, which fixes the ORDER in which the subshells of a CSF are
        named when one is printed. Nothing is returned.

        THIS IS FOR DISPLAY ONLY, and it is internal. It is set in many places but read in exactly two, both of
        which produce output: `Base.string(csf::CsfR)` and `Basics.display(stream, orbitals, grid)`. No
        computation consults it -- every process module works from its own local subshell list, obtained from
        `Basics.generate(OrderedSubshellList(), basisA, basisB)` and passed explicitly to whatever needs it.

        IT IS NOT A USER SETTING, which is why it was demoted out of `Defaults.setDefaults` on 12-Aug-2026
        (the string key still delegates here, so existing scripts keep working). The correspondence to a CSF is
        POSITIONAL -- `csf.occupation[i]` belongs to `subshells[i]` -- and nothing validates the list against
        the CSFs it will label. A list that does not match therefore mislabels every printed CSF silently. The
        value is derived data with no physical content, so there is nothing a user answers by overriding it.
"""
function setStandardSubshellList(subshells::Array{Subshell,1}; printout::Bool=true)
    if printout    println("(Re-) Define a new standard subshell list.")    end
    global GBL_STANDARD_SUBSHELL_LIST = deepcopy(subshells)

    nothing
end


"""
+ `("standard grid", grid::Radial.Grid; printout::Bool=true)`  
    ... to (pre-) define internally the standard radial grid which is used to represent most orbitals.
"""
function setDefaults(sa::String, grid::Radial.Grid; printout::Bool=true)
    global GBL_STANDARD_GRID

    if        sa == "standard grid"
        if  printout    println("(Re-) Define the standard grid with $(grid.NoPoints) grid points.")    end
        GBL_STANDARD_GRID = grid
    else
    error("Unsupported keystring:: $sa")
    end

    nothing
end


"""
+ `("continuum: potential", scField::Basics.AbstractScField)`  
    ... to (re-) define the potential that is applied for the generation of the continuum orbitals.
"""
function setDefaults(sa::String, scField::Basics.AbstractScField)
    global GBL_CONT_POTENTIAL

    if        sa == "continuum: SCF potential"
        GBL_CONT_POTENTIAL = scField
    else
    error("Unsupported keystring:: $sa")
    end

    nothing
end


"""
+ `("QED: damped-hydrogenic", Znuc::Float64, wa::Array{Float64,1})`  
    ... to (re-) define the lambda-C damped overlap integrals of the lowest kappa-orbitals 
        [ wa_1s_1/2, wa_2p_1/2, wa_2p_3/2, wa_3d_3/2, wa_3d_5/2 ] for the (new) nuclear charge Znuc; 
        nothing is returned.
"""
function setDefaults(sa::String, Znuc::Float64, wa::Array{Float64,1})
    global GBL_QED_HYDROGENIC_LAMBDAC,  GBL_QED_NUCLEAR_CHARGE

    if        sa == "QED: damped-hydrogenic"
        if  length(wa) != 5  error("Defaults.setDefaults(\"QED: damped-hydrogenic\"): the array must hold five " *
                                   "damping parameters; got $(length(wa)).")   end
        println("Re-define the damped overlap integrals < a | e^{-r/lambda_C} | a > of the lowest kappa-orbitals for nuclear charge Z = " *
                "$(Znuc).")
        GBL_QED_NUCLEAR_CHARGE     = Znuc
        GBL_QED_HYDROGENIC_LAMBDAC = wa
    else
    error("Unsupported keystring:: $sa")
    end

    nothing
end

    
    
"""
`Defaults.getDefaults()`  
    ... gives/supplies different information about the (present) framework of the computation or about some 
        given data; cf. Defaults.setDefaults(). 

+ `("alpha")`  or  `("fine-structure constant alpha")`   
    ... to get the (current) value::Float64 of the fine-structure constant alpha.

+ `("electron mass: kg")`  or  `("electron mass: amu")`  
    ... to get the (current) value::Float64 of the electron mass in the specified unit.

+ `("framework")`  ... to give the (current) setting::String  of the overall framework.

+ `("electron rest energy")`  or  `("mc^2")`  ... to get the electron rest energy.

+ `("electron g-factor")`  ... to give the electron g-factor g_s = 2.00232.

+ `("nuclear: charge")`  or  `("nuclear: model")` 
    ... to get the nuclear charge or model, if it has been set before.

+ `("unit: energy")`  or  `("unit: cross section")`  or  `("unit: rate")`  or  `("unit: strength")`  or  `("unit: time")`  
    ... to get the corresponding (user-defined) unit::String for the current computations.

+ `("standard grid")`  
    ... to get the (current standard) grid::Array{Float64,1} to which all radial orbital functions usually refer.

+ `("speed of light: c")`  ... to get the speed of light in atomic units.
+ `("mass: muon")`  ... to get the muon mass in atomic units, i.e. in units of the electron mass.

+ `("summary flag/stream")`  
    ... to get the logical flag and stream for printing a summary file; a tupel (flag, iostream) is returned.
"""
function getDefaults(sa::String)
    global GBL_FRAMEWORK, GBL_CONTINUUM, GBL_ENERGY_UNIT, GBL_CROSS_SECTION_UNIT, GBL_RATE_UNIT, GBL_STRENGTH_UNIT, GBL_TIME_UNIT
    global ELECTRON_MASS_SI, ELECTRON_MASS_U, FINE_STRUCTURE_CONSTANT, GBL_STANDARD_GRID, 
            GBL_DATA_IOSTREAM, GBL_DATAX_IOSTREAM, GBL_DATAY_IOSTREAM, GBL_SUMMARY_IOSTREAM, GBL_PRINT_SUMMARY, GBL_TEST_IOSTREAM, GBL_TEST_SUMMARY

    if        sa in ["alpha", "fine-structure constant alpha"]          return (FINE_STRUCTURE_CONSTANT)   
    elseif    sa == "electron mass: kg"                                 return (ELECTRON_MASS_SI)   
    elseif    sa == "electron mass: amu"                                return (ELECTRON_MASS_U)   
    elseif    sa == "electron g-factor"                                 return (2.00232)   
    elseif    sa == "framework"                                         return (GBL_FRAMEWORK)
    elseif    sa == "method: continuum"                                 return (GBL_CONTINUUM)
    elseif    sa in ["electron rest energy", "mc^2"]                    return (1/(FINE_STRUCTURE_CONSTANT^2))   
    elseif    sa == "nuclear: charge"                                   return (GBL_NUCLEAR_CHARGE)
    elseif    sa == "nuclear: model"                                    return (GBL_NUCLEAR_MODEL)
    elseif    sa == "unit: energy"                                      return (GBL_ENERGY_UNIT)
    elseif    sa == "unit: cross section"                               return (GBL_CROSS_SECTION_UNIT)
    elseif    sa == "unit: energy-diff. cross section"                  return (GBL_EDIFF_CROSS_SECTION_UNIT)
    elseif    sa == "unit: rate"                                        return (GBL_RATE_UNIT)
    elseif    sa == "unit: strength"                                    return (GBL_STRENGTH_UNIT)
    elseif    sa == "unit: time"                                        return (GBL_TIME_UNIT)
    elseif    sa == "standard grid"                                     return (GBL_STANDARD_GRID)
    elseif    sa == "speed of light: c"                                 return (1.0/FINE_STRUCTURE_CONSTANT)
    elseif    sa == "mass: muon"                                        return (MUON_ELECTRON_MASS_RATIO)
    elseif    sa == "data flag/stream"                                  return ( (GBL_PRINT_DATA, GBL_DATA_IOSTREAM) )
    elseif    sa == "data-X flag/stream"                                return ( (GBL_PRINT_DATAX, GBL_DATAX_IOSTREAM) )
    elseif    sa == "data-Y flag/stream"                                return ( (GBL_PRINT_DATAY, GBL_DATAY_IOSTREAM) )
    elseif    sa == "summary flag/stream"                               return ( (GBL_PRINT_SUMMARY, GBL_SUMMARY_IOSTREAM) )
    elseif    sa == "test flag/stream"                                  return ( (GBL_PRINT_TEST,    JenaAtomicCalculator.JAC_TEST_IOSTREAM) )
    else      error("Unsupported keystring:: $sa")
    end
end


"""
+ `("ordered shell list: non-relativistic", n_max::Int64)`  
    ... to give an ordered list of non-relativistic shells::Array{Shell,1} up to the (maximum) principal number n_max.
        
+ `("ordered subshell list: relativistic", n_max::Int64)`   
    ... to give an ordered list of relativistic subshells::Array{Subshell,1} up to the (maximum) principal number n_max.
"""
function getDefaults(sa::String, n_max::Int64)

    !(1 <= n_max < 12)    &&    error("Unsupported value of n_max = $n_max")

    if        sa == "ordered shell list: non-relativistic"
        wa = Shell[]
        for  n = 1:n_max
        for  l = 0:n-1    push!(wa, Shell(n,l) )    end
        end
    elseif    sa == "ordered subshell list: relativistic"  
        wa = Subshell[]
        for  n = 1:n_max
            for l = 0:n_max-1
                if     l == 0   push!(wa, Subshell(n, -1) )
                else   j2 = 2l - 1;    kappa = Int64(  (j2+1)/2 );   push!(wa, Subshell(n, kappa) )
                    j2 = 2l + 1;    kappa = Int64( -(j2+1)/2 );   push!(wa, Subshell(n, kappa) )
                end
            end
        end
    else    error("Unsupported keystring:: $sa")
    end
    return( wa )
end


    
    
"""
`Defaults.warn()`  
    ... deals with warnings that occur during a run and session; it handles the global array GBL_WARNINGS.

+ `(AddWarning, warning::String)`  ... to add warning to the global array GBL_WARNINGS.

+ `(PrintWarnings)`  ... to print all warnings that are currently kept in the global array GBL_WARNINGS.

+ `(ResetWarnings)`  ... to reset the global array GBL_WARNINGS.
"""
function warn(wa::AbstractWarning)
    global GBL_WARNINGS, GBL_WARNINGS_COUNT

    if        wa == PrintWarnings()
        lock(GBL_WARNINGS_LOCK)
        try
            iostream = open("jac-warn.report", "w")
            println(iostream, " ")
            println(iostream, "\n\n",
                            "Warnings from the present sessions or run, printed at $( string(now())[1:16] ): \n",
                            "======================================================================= \n")
            if  length(GBL_WARNINGS) == 0    println(iostream, "++ none.")    end
            for sa in  GBL_WARNINGS
                nx = get(GBL_WARNINGS_COUNT, sa, 1)
                println(iostream, "++ " * sa * (nx > 1  ?  "   [$nx times]"  :  ""))
            end
            close(iostream)
        finally
            unlock(GBL_WARNINGS_LOCK)
        end
        #
    elseif    wa == ResetWarnings()
        lock(GBL_WARNINGS_LOCK)
        try     GBL_WARNINGS = String[];    GBL_WARNINGS_COUNT = Dict{String,Int64}()
        finally unlock(GBL_WARNINGS_LOCK)
        end
    else      error("Unsupported Warnings:: $wa")
    end

    nothing
end


function warn(wa::AbstractWarning, sa::String)
    global GBL_WARNINGS, GBL_WARNINGS_COUNT

    if        wa == AddWarning()
        ## Recorded once, but counted every time: a condition met by 500 lines is one entry saying so, not
        ## 500 entries.  The lock makes this safe to call from a threaded line loop; it is uncontended in
        ## the ordinary single-threaded case and costs nothing against the work that produced the warning.
        lock(GBL_WARNINGS_LOCK)
        try
            nx = get(GBL_WARNINGS_COUNT, sa, 0)
            GBL_WARNINGS_COUNT[sa] = nx + 1
            if  nx == 0     push!(GBL_WARNINGS, sa)     end
        finally
            unlock(GBL_WARNINGS_LOCK)
        end
    else      error("Unsupported Warnings:: $wa")
    end

    nothing
end

    
    
"""
`Defaults.saRatip()`  
    ... returns true or false if the spin-angular coefficients should be calculated by Fortran codes from RATIP/Grasp.
"""
saRatip() = false
saGG()    = true


end # module



"""
`module  JAC.StarkShift`
... a submodel of JAC that contains all methods for computing the static electric-dipole/-quadrupole
    moments (Theta^(1), Theta^(2)) and the static Stark shift of fine-structure levels in a uniform
    electric field.

    The Stark shift of a level |J,M> in a static field E is the direct physical use of the SAME
    static E1 scalar/tensor polarizability alpha_0(J)/alpha_2(J) that module-MultipolePolarizibility.jl
    computes: DeltaE(J,M) = -(1/2) alpha_0(J) E^2 - (1/2) alpha_2(J) E^2 [3M^2-J(J+1)] / [J(2J-1)]
    (standard second-order-perturbation-theory form; the Angel & Sandars alpha_2 normalization is
    defined precisely so this simple form applies -- see van Leeuwen & Hogervorst, Z. Phys. A 316, 149
    (1984), Eq. (5), the same paper module-MultipolePolarizibility.jl's docstring already references).
    Rather than re-deriving/re-validating this sum-over-states here, StarkShift.computeStarkShifts
    calls MultipolePolarizibility.computeScalarTensorPolarizability directly, on a caller-supplied,
    independently-generated StarkShift.Settings.gMultiplet (mirrors
    MultipolePolarizibility.Settings.gMultiplet/MultiPhotonDeExcitation.Settings.gMultiplet exactly --
    see examples/example-Cl.jl), and only adds the field-strength/M-dependent energy-shift conversion
    on top. Both Coulomb and Babushkin gauge are reported throughout (EmProperty), not calibrated onto
    a common scale, following the same convention already established for MultipolePolarizibility.

    THE AC SHIFT IS THE SAME EXPRESSION WITH alpha(omega) IN PLACE OF alpha(0), AND ONE FACTOR OF TWO
    THAT IS NOT A DETAIL.  For a field oscillating as E(t) = E_0 cos(omega t) the shift is averaged
    over the cycle, and <cos^2> = 1/2, so

        STATIC :  DeltaE(J,M) = -(1/2) alpha_0(0) E^2   - (1/2) alpha_2(0) E^2   [3M^2-J(J+1)]/[J(2J-1)]
        AC     :  DeltaE(J,M) = -(1/4) alpha_0(om) E_0^2 - (1/4) alpha_2(om) E_0^2 [3M^2-J(J+1)]/[J(2J-1)]

    where E in the first line is the CONSTANT field and E_0 in the second is the PEAK amplitude.  The
    omega -> 0 limit of the ac line is therefore HALF the static line, and that is correct physics
    rather than a discontinuity: a field held at E and a field swinging between +E and -E are not the
    same field.  `Settings.omega` selects between the two -- 0. means static, and nothing about the
    static path changes -- and `Settings.EField` is read as the constant field or as the peak
    amplitude accordingly.  THE HALVING IS ASSERTED IN THE TEST SUITE, since a factor of two silently
    adopted in one convention and quoted in the other is this quantity's classic error.

    WHAT THE AC SHIFT IS GOOD FOR.  -DeltaE at the trap frequency IS the optical-trap depth of that
    sublevel, so the M-dependence of the table is the differential light shift that a clock or a
    quantum-logic experiment has to control;  and a MAGIC WAVELENGTH is a zero of
    alpha_0(omega, upper) - alpha_0(omega, lower), which is read off by running two levels over a
    range of `omega`.  Both follow from the dynamic polarizability and need nothing further.

    THE BLACK-BODY SHIFT IS THE SAME EXPRESSION AVERAGED OVER THE PLANCK SPECTRUM, and it needs no
    further physics:  a level in thermal radiation at temperature T is shifted by

        DeltaE_BBR  =  -(1/2) <E^2> alpha_0(0) (1 + eta)        <E^2> = (4 pi^3/15) alpha_fs^3 T^4

    with <E^2> the MEAN-SQUARE field of the radiation -- already a time average, so there is NO second
    factor of two here, unlike the ac line above -- and eta the dynamic correction that comes from
    alpha_0(omega) not being flat across the spectrum.  The leading term is the universally quoted
    -(2/15)(alpha_fs pi)^3 T^4 alpha_0, which is the same number written differently.

    ONLY THE RANK-0 POLARIZABILITY ENTERS, and that is why a black-body shift is ONE number per level
    rather than a table over M.  Thermal radiation is isotropic and unpolarized, so its field
    correlation tensor is <E_i E_j> = (1/3) delta_ij <E^2>:  the rank-1 component averages to zero for
    want of a net handedness and the rank-2 component for want of a preferred axis.  `Settings.temperature`
    is read in KELVIN, to match `Settings.EField` being read in V/cm, and 0. -- the default -- means no
    black-body shift is computed and nothing about the other two paths changes.

    THE REGIME THAT MUST BE WATCHED IS THE RYDBERG ONE.  For a low-lying level at room temperature the
    Planck weight peaks around 2.8 kT, far BELOW every excitation energy, so alpha_0(omega) is flat over
    the whole spectrum and eta is a small correction -- measured on H 1s, +1.2e-04 at 300 K.  For a
    Rydberg level the excitation energies fall as 1/n^3 and the Planck weight sits ON them:  the integral
    is then a PRINCIPAL VALUE whose imaginary part is the BBR-INDUCED TRANSITION RATE, and neither is
    what a plain quadrature returns.  StarkShift.checkBlackbodyResonance measures the fraction of the
    Planck weight lying above the lowest excitation energy and warns when it is not negligible.

    THE POLE STRUCTURE IS INHERITED, NOT RE-CHECKED HERE.  Above the lowest excitation energy of the
    level the true polarizability is COMPLEX, and near a pole it rides on one energy denominator;
    `MultipolePolarizibility.checkFrequency` warns about both and runs on every call made from here,
    so a shift quoted close to a resonance carries its own warning.
"""
module StarkShift


using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..MultipoleMoment,
              ..MultipolePolarizibility, ..Nuclear, ..Radial, ..TableStrings


"""
`struct  StarkShift.SublevelJ`  ... defines a type to specify a Stark-shifted sublevel with well-defined J.

    + M                      ::AngularM64        ... M_J-value
    + energy                 ::EmProperty        ... Stark energy shift of this sublevel, in both the
                                                       Coulomb and Babushkin gauge (see module docstring).
"""
struct SublevelJ
    M                        ::AngularM64
    energy                   ::EmProperty
end


# `Base.show(io::IO, Jsublevel::StarkShift.SublevelJ)`  ... prepares a proper printout of the variable Jsublevel::StarkShift.SublevelJ.
function Base.show(io::IO, Jsublevel::StarkShift.SublevelJ) 
    println(io, "Sublevel [M=$(Jsublevel.M); energy = $(Jsublevel.energy)]")
end


"""
`struct  StarkShift.MagicWavelength`
    ... one magic wavelength of a PAIR of levels, i.e. one zero of alpha_0(omega, A) - alpha_0(omega, B). At such a
        frequency the two levels are shifted EQUALLY, so the transition between them is untouched by the trapping
        light to first order in the intensity -- which is what makes an optical trap compatible with a clock or a
        qubit. There are generically MANY of them, one or more between every pair of consecutive resonances, so
        StarkShift.magicWavelengths returns a list and never "the" magic wavelength.

    + omega          ::Float64   ... the frequency, in atomic units.
    + wavelength     ::Float64   ... the same frequency as a vacuum wavelength, in nm, which is how this quantity
                                       is universally quoted.
    + alphaCommon    ::Float64   ... the polarizability the two levels SHARE there, in a.u. IT IS THE TRAP DEPTH per
                                       unit intensity, and ITS SIGN DECIDES WHETHER THE TRAP TRAPS: a magic
                                       wavelength with alphaCommon < 0 shifts both levels UP and repels the atom
                                       from the focus, so it is magic and useless at the same time.
    + slope          ::Float64   ... d(alpha_A - alpha_B)/d omega at the zero, in a.u. per a.u. It says how
                                       SHARP the magic wavelength is: a small slope means a wide tolerance on the
                                       laser frequency, a large one means the differential shift returns at once.
    + aboveResonance ::Bool      ... true if omega lies at or above the lowest excitation energy of either level,
                                       where the true polarizability is COMPLEX and this real sum is a
                                       principal-value-like part with no width. Such a zero is still a real
                                       crossing of the two real parts, but it is not a usable magic wavelength
                                       without the widths.
"""
struct MagicWavelength
    omega            ::Float64
    wavelength       ::Float64
    alphaCommon      ::Float64
    slope            ::Float64
    aboveResonance   ::Bool
end


# `Base.show(io::IO, mw::StarkShift.MagicWavelength)`
#        ... prepares a proper printout of the variable mw::StarkShift.MagicWavelength.
function Base.show(io::IO, mw::StarkShift.MagicWavelength)
    sa = mw.aboveResonance ? "  [ABOVE a resonance: a principal value, not a polarizability]" : ""
    sb = mw.alphaCommon < 0.0 ? "  [alphaCommon < 0: this trap REPELS]" : ""
    println(io, "magic wavelength:  omega = $(mw.omega) a.u. = $(mw.wavelength) nm;  shared alpha_0 = " *
                "$(mw.alphaCommon) a.u.;  slope = $(mw.slope) a.u./a.u." * sa * sb)
end


"""
`struct  StarkShift.Outcome`
    ... defines a type to keep the Stark shift and electric-quadrupole-moment parameters of a fine-structure level.

    + Jlevel                 ::Level             ... Fine-structure levels to which the results refer to.
    + Theta1                 ::Float64           ... Electric-dipole moment of the atom in this level.
    + Theta2                 ::Float64           ... Electric-quadrupole moment of the atom in this level.
    + alpha0                 ::EmProperty        ... Static scalar (rank 0) E1 polarizability alpha_0(J), in
                                                       a.u., that the Stark shift below was computed from
                                                       (see module docstring); EmProperty(0.) if calcStarkshifts
                                                       was not requested.
    + alpha2                 ::EmProperty        ... Static tensor (rank 2) E1 polarizability alpha_2(J), in
                                                       a.u.; identically 0. for J<1 (see module docstring).
    + Jsublevels              ::Array{StarkShift.SublevelJ,1}
        ... List of the Stark-shifted fine-structure sublevels and data with well-defined J-value
    + bbrShift               ::EmProperty        ... Black-body radiation shift of the level at
                                                       Settings.temperature, in a.u.; 0. if no temperature was
                                                       given. It carries NO M-dependence, by the isotropy
                                                       argument in the module docstring.
    + magicWavelengths       ::Array{StarkShift.MagicWavelength,1}
                                                 ... The magic wavelengths of THIS level paired with
                                                       settings.magicPartnerIndex, empty unless
                                                       settings.calcMagicWavelengths. A magic wavelength belongs
                                                       to the PAIR, so it is stored on the outcome of the first
                                                       of the two and the partner's own outcome carries the same
                                                       list;  which is which is in the printed table.
    + bbrEta                 ::Float64           ... The dynamic correction eta of that shift, i.e. how far
                                                       alpha_0(omega) departs from alpha_0(0) across the Planck
                                                       spectrum. A small eta says the static polarizability was
                                                       sufficient; a large one says the spectrum reaches the
                                                       level's own excitation energies.
"""
struct Outcome
    Jlevel                   ::Level
    Theta1                   ::Float64
    Theta2                   ::Float64
    alpha0                   ::EmProperty
    alpha2                   ::EmProperty
    Jsublevels                ::Array{StarkShift.SublevelJ,1}
    bbrShift                 ::EmProperty
    bbrEta                   ::Float64
    magicWavelengths         ::Array{StarkShift.MagicWavelength,1}
end


"""
`StarkShift.Outcome()`  ... constructor for an `empty` instance of StarkShift.Outcome.
"""
function Outcome()
    Outcome(Level(), 0., 0., EmProperty(0.), EmProperty(0.), StarkShift.SublevelJ[], EmProperty(0.), 0.,
            StarkShift.MagicWavelength[])
end


# `Base.show(io::IO, outcome::StarkShift.Outcome)`  ... prepares a proper printout of the variable StarkShift.Outcome.
function Base.show(io::IO, outcome::StarkShift.Outcome)
    println(io, "Jlevel:           $(outcome.Jlevel)  ")
    println(io, "Theta1:           $(outcome.Theta1)  ")
    println(io, "Theta2:           $(outcome.Theta2)  ")
    println(io, "alpha0:           $(outcome.alpha0)  ")
    println(io, "alpha2:           $(outcome.alpha2)  ")
    println(io, "Jsublevels:       $(outcome.Jsublevels)  ")
    println(io, "bbrShift:         $(outcome.bbrShift)  ")
    println(io, "bbrEta:           $(outcome.bbrEta)  ")
end


"""
`struct  StarkShift.Settings  <:  AbstractPropertySettings`
    ... defines a type for the details and parameters of computing the Stark shift and electric-quadrupole-moments
        of fine-structure levels.

    + calcEDM                ::Bool              ... True if the EDM of selected levels is to be calculated.
    + calcEQM                ::Bool              ... True if the EQM of selected levels is to be calculated.
    + calcStarkshifts        ::Bool              ... True if the Stark energy shifts are to be computed.
    + calcMagicWavelengths   ::Bool              ... True if the MAGIC WAVELENGTHS of a PAIR of levels are to be
                                                       computed, i.e. the frequencies at which the two levels are
                                                       shifted equally and the transition between them is
                                                       insensitive to the trapping field. Needs gMultiplet and
                                                       magicPartnerIndex;  both are checked.
    + gMultiplet              ::Multiplet         ... A local, approximate Green multiplet of intermediate
                                                       (opposite-parity) ASF that the underlying alpha_0/alpha_2
                                                       sum-over-states runs over -- generated INDEPENDENTLY,
                                                       ahead of time, by the caller (see the module docstring
                                                       and examples/example-Cl.jl); only used if calcStarkshifts.
    + printBefore             ::Bool              ... True if a list of selected levels is printed before the actual computations start.
    + EField                 ::Float64           ... Strength of the electric field in [V/cm]: the CONSTANT
                                                       field if omega = 0., the PEAK amplitude E_0 otherwise
                                                       (see the module docstring on the factor of two).
    + omega                  ::Float64           ... Angular frequency of the driving field, in atomic units.
                                                       0., the default, means a STATIC field and leaves that
                                                       path bit-identical to what it was.
    + temperature            ::Float64           ... Temperature of a black-body radiation field, in KELVIN (the
                                                       same user-facing convention as EField in V/cm). 0., the
                                                       default, means no black-body shift is computed. It is
                                                       INDEPENDENT of EField and omega: a black-body shift is an
                                                       average over the whole Planck spectrum and does not read
                                                       either of them.
    + magicPartnerIndex      ::Int64             ... Index, IN THE COMPUTED MULTIPLET, of the level each selected
                                                       level is paired with for a magic wavelength. A magic
                                                       wavelength belongs to a PAIR and not to a level, and the
                                                       partner is always another level of the same atom, so an
                                                       index into the same multiplet names it -- the same
                                                       vocabulary levelSelection already uses. 0 means none.
    + magicOmegas            ::Array{Float64,1}  ... The frequency window [omegaMin, omegaMax] in atomic units
                                                       within which magic wavelengths are searched. EMPTY means
                                                       the documented defaults of StarkShift.magicWavelengths
                                                       (0.01 to 0.25 a.u.); any other length than 0 or 2 raises.
    + levelSelection         ::LevelSelection    ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings
    calcEDM                  ::Bool
    calcEQM                  ::Bool
    calcStarkshifts          ::Bool
    calcMagicWavelengths     ::Bool
    gMultiplet                ::Multiplet
    printBefore               ::Bool
    EField                   ::Float64
    omega                    ::Float64
    temperature              ::Float64
    magicPartnerIndex        ::Int64
    magicOmegas              ::Array{Float64,1}
    levelSelection           ::LevelSelection
end


"""
`StarkShift.Settings()`
    ... constructor for an `empty` instance of StarkShift.Settings.
"""
function Settings()
        Settings(false, false, false, false, Multiplet(), false, 0., 0., 0., 0, Float64[], LevelSelection() )
end


"""
`StarkShift.Settings(set::StarkShift.Settings;`

        calcEDM=.., calcEQM=.., calcStarkshifts=.., calcMagicWavelengths=.., gMultiplet=.., printBefore=..,
        EField=.., omega=.., temperature=.., magicPartnerIndex=.., magicOmegas=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::StarkShift.Settings.
"""
function Settings(set::StarkShift.Settings;
        calcEDM::Union{Nothing,Bool}=nothing,                   calcEQM::Union{Nothing,Bool}=nothing,
        calcStarkshifts::Union{Nothing,Bool}=nothing,           calcMagicWavelengths::Union{Nothing,Bool}=nothing,
        gMultiplet::Union{Nothing,Multiplet}=nothing,           printBefore::Union{Nothing,Bool}=nothing,
        EField::Union{Nothing,Float64}=nothing,                 omega::Union{Nothing,Float64}=nothing,
        temperature::Union{Nothing,Float64}=nothing,            magicPartnerIndex::Union{Nothing,Int64}=nothing,
        magicOmegas::Union{Nothing,Array{Float64,1}}=nothing,   levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  isnothing(calcEDM)           calcEDMx          = set.calcEDM          else   calcEDMx          = calcEDM          end
    if  isnothing(calcEQM)           calcEQMx          = set.calcEQM          else   calcEQMx          = calcEQM          end
    if  isnothing(calcStarkshifts)   calcStarkshiftsx  = set.calcStarkshifts  else   calcStarkshiftsx  = calcStarkshifts  end
    if  isnothing(calcMagicWavelengths)  calcMagicWavelengthsx = set.calcMagicWavelengths  else  calcMagicWavelengthsx = calcMagicWavelengths  end
    if  isnothing(gMultiplet)        gMultipletx       = set.gMultiplet       else   gMultipletx       = gMultiplet       end
    if  isnothing(printBefore)       printBeforex      = set.printBefore      else   printBeforex      = printBefore      end
    if  isnothing(EField)            EFieldx           = set.EField           else   EFieldx           = EField           end
    if  isnothing(omega)             omegax            = set.omega            else   omegax            = omega            end
    if  isnothing(temperature)       temperaturex      = set.temperature      else   temperaturex      = temperature      end
    if  isnothing(magicPartnerIndex)  magicPartnerIndexx = set.magicPartnerIndex  else  magicPartnerIndexx = magicPartnerIndex  end
    if  isnothing(magicOmegas)       magicOmegasx      = set.magicOmegas      else   magicOmegasx      = magicOmegas      end
    if  isnothing(levelSelection)    levelSelectionx   = set.levelSelection   else   levelSelectionx   = levelSelection   end

    Settings( calcEDMx, calcEQMx, calcStarkshiftsx, calcMagicWavelengthsx, gMultipletx, printBeforex, EFieldx,
              omegax, temperaturex, magicPartnerIndexx, magicOmegasx, levelSelectionx )
end


# `Base.show(io::IO, settings::StarkShift.Settings)`  ... prepares a proper printout of the variable settings::StarkShift.Settings.
function Base.show(io::IO, settings::StarkShift.Settings)
    println(io, "calcEDM:               $(settings.calcEDM)  ")
    println(io, "calcEQM:               $(settings.calcEQM)  ")
    println(io, "calcStarkshifts:       $(settings.calcStarkshifts)  ")
    println(io, "calcMagicWavelengths:  $(settings.calcMagicWavelengths)  ")
    println(io, "gMultiplet:            $(settings.gMultiplet)  ")
    println(io, "printBefore:           $(settings.printBefore)  ")
    println(io, "EField:                $(settings.EField)  ")
    println(io, "omega:                 $(settings.omega)  ")
    println(io, "temperature:           $(settings.temperature)  ")
    println(io, "magicPartnerIndex:     $(settings.magicPartnerIndex)  ")
    println(io, "magicOmegas:           $(settings.magicOmegas)  ")
    println(io, "levelSelection:        $(settings.levelSelection)  ")
end


"""
`StarkShift.computeAmplitudesProperties(outcome::StarkShift.Outcome, grid::Radial.Grid, settings::StarkShift.Settings)`
    ... to compute all amplitudes and properties of for a given level. An outcome::StarkShift.Outcome is returned for
        which the amplitudes and all requested properties are now evaluated explicitly.
"""
function  computeAmplitudesProperties(outcome::StarkShift.Outcome, grid::Radial.Grid, settings::StarkShift.Settings)
    edm = eqm = 0.
    if  settings.calcEDM    edm = StarkShift.computeTheta1(outcome.Jlevel, grid)    end
    if  settings.calcEQM    eqm = StarkShift.computeTheta2(outcome.Jlevel, grid)    end

    alpha0 = EmProperty(0.);   alpha2 = EmProperty(0.);   Jsublevels = outcome.Jsublevels
    if  settings.calcStarkshifts
        alpha0, alpha2, Jsublevels = StarkShift.computeStarkShifts(outcome.Jlevel, settings.gMultiplet, grid, settings)
    end

    ## THE BLACK-BODY SHIFT IS INDEPENDENT OF EField AND omega and is therefore gated on its own temperature
    ## alone, not on calcStarkshifts: a caller may want the thermal shift and no applied field at all.
    bbrShift = EmProperty(0.);   bbrEta = 0.
    if  settings.temperature > 0.0
        bbrShift, bbrEta = StarkShift.computeBlackbodyShift(outcome.Jlevel, settings.gMultiplet, grid,
                                                             settings.temperature)
    end

    newOutcome = StarkShift.Outcome( outcome.Jlevel, edm, eqm, alpha0, alpha2, Jsublevels, bbrShift, bbrEta,
                                     outcome.magicWavelengths)
    return( newOutcome )
end


"""
`StarkShift.AU_EFIELD_IN_VCM`
    ... the CODATA atomic unit of electric field, in V/cm (1 a.u. = 5.142206751e9 V/cm; verified
        against NIST/CODATA 31-Jul-2026, not from memory alone), used to convert
        StarkShift.Settings.EField (documented in V/cm) to atomic units.
"""
const AU_EFIELD_IN_VCM = 5.142206751e9


"""
`StarkShift.blackbodyFieldSquare(temperature::Float64)`  
    ... returns the mean-square electric field <E^2> of an isotropic black-body radiation field at `temperature`,
        given in KELVIN as StarkShift.Settings.temperature is, in ATOMIC UNITS. It is four pi times the
        Stefan-Boltzmann energy density u = (pi^2/15) alpha_fs^3 T^4, since u = <E^2>/(4 pi) for a field with
        <E^2> = <B^2>, so

            <E^2>  =  (4 pi^3/15) alpha_fs^3 T^4

        with T the temperature in Hartree. A value::Float64 is returned, and 0. for a non-positive temperature.

        THIS IS A COMPLETE EXTERNAL ANCHOR AND COSTS NOTHING TO CHECK, because no atomic quantity enters it: only
        the fine-structure constant does. At 300 K it gives E_rms = 831.945 V/m against the textbook 831.9 V/m,
        and the leading shift it implies for a level of alpha_0 = 1 a.u. is -0.008611 Hz. The same number by
        QUADRATURE over JAC's own `Distribution.photonNumberDensity(Distribution.PhotonPlanck(T), omega)` agrees
        to 1.3e-10, which is how the Planck normalisation used in StarkShift.computeBlackbodyShift was checked;
        that module is deliberately not a dependency here, since the closed form is exact and the spectral weight
        is one expression.
"""
function blackbodyFieldSquare(temperature::Float64)
    if  temperature <= 0.0    return( 0.0 )    end
    T     = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temperature)
    alpha = 1.0 / Defaults.getDefaults("speed of light: c")

    return( 4.0 * pi^3 / 15.0 * alpha^3 * T^4 )
end


"""
`StarkShift.checkBlackbodyResonance(temperature::Float64, energies::Array{Float64,1}, nPoints::Int64,
                                        xMax::Float64)`  
    ... inspects whether the Planck spectrum at `temperature` reaches the level's own excitation energies and
        WARNS LOUDLY, to stdout, when it does -- which is the one way a black-body shift comes out wrong while
        looking perfectly reasonable. It measures the FRACTION of the Planck weight x^3/(exp(x)-1) lying ABOVE
        the lowest excitation energy, since that is the part of the spectrum for which alpha_0(omega) is not a
        polarizability at all, and warns when that fraction exceeds 1e-3. Nothing is returned.

        WHY A FRACTION RATHER THAN A RULE OF THUMB. The weight peaks at x = 2.82 and is spent by x = 15, so
        "kT well below the first excitation energy" is the usual statement; but what decides whether a number is
        usable is how much weight actually sits above it, and that is one integral. For a low-lying level at
        room temperature the fraction is utterly negligible and nothing is printed. For a RYDBERG level the
        excitation energies fall as 1/n^3 and the fraction approaches one.

        WHAT IS WRONG IN THAT REGIME IS NOT PRECISION BUT THE QUANTITY. Above the lowest excitation energy the
        true alpha_0(omega) is COMPLEX, so the integral is a PRINCIPAL VALUE and its imaginary part is the
        BBR-INDUCED TRANSITION RATE -- a real, separate observable. A plain quadrature gives neither: it returns
        whatever the point nearest a pole happened to contribute, so the distance from the nearest quadrature
        point to the lowest excitation energy is reported too.
"""
function checkBlackbodyResonance(temperature::Float64, energies::Array{Float64,1}, nPoints::Int64, xMax::Float64)
    T   = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temperature)
    dEs = [dE  for dE in energies  if dE > 0.0]
    if  isempty(dEs)    return( nothing )    end
    dEmin = minimum(dEs);      xMin = dEmin / T
    if  xMin >= xMax    return( nothing )    end

    # the fraction of the Planck weight above the lowest excitation energy, by Simpson on [xMin, xMax]
    nf = 2000;    hf = (xMax - xMin)/nf;    wa = 0.0
    for  i = 0:nf
        x  = xMin + i*hf
        fx = x == 0.0 ? 0.0 : x^3 / (exp(x) - 1.0)
        wa = wa + ((i == 0 || i == nf) ? 1.0 : (isodd(i) ? 4.0 : 2.0)) * fx
    end
    frac = (hf/3.0 * wa) / (pi^4/15.0)
    if  frac <= 1.0e-3    return( nothing )    end

    h     = xMax/nPoints
    near  = minimum( [abs(T*i*h - dEmin)  for i = 0:nPoints] ) / dEmin
    printstyled("\n>>> THE PLANCK SPECTRUM AT T = " * @sprintf("%.1f", temperature) * " K REACHES THIS LEVEL'S OWN " *
                "EXCITATION ENERGIES: " * @sprintf("%.1f", 100*frac) * " % of the\n" *
                "    weight lies above the lowest one at " * @sprintf("%.6f", dEmin) * " a.u. There alpha_0(omega) " *
                "is COMPLEX, so the integral is a\n    PRINCIPAL VALUE and its imaginary part is the BBR-induced " *
                "transition rate -- and a plain quadrature\n    returns neither. The nearest quadrature point sits " *
                @sprintf("%.2f", 100*near) * " % away from that energy, so the number below is\n" *
                "    carried by whichever point came closest to the pole. DO NOT QUOTE IT.\n", color=:light_red)

    return( nothing )
end


"""
`StarkShift.polarizabilityTerms(level::Level, gMultiplet::Multiplet, grid::Radial.Grid)`  
    ... returns the per-perturber pieces from which the rank-0 E1 polarizability of `level` can be rebuilt at ANY
        frequency by arithmetic alone,

            alpha_0(omega) = sum_a strengths[a] * MultipolePolarizibility.frequencyFactor(energies[a], omega, false)

        with strengths[a] = c(L=1, K=0, J, J_a) |<J_a||O^(E1)||J>|^2 and energies[a] = W_a - W_J. The sum over states
        is therefore evaluated ONCE, and that is what makes a frequency quadrature or a root search affordable: a few
        hundred frequencies would otherwise be a few hundred sum-over-states evaluations, each rebuilding subshell
        unions and reduced matrix elements. A tuple (strengths::Array{Float64,1}, energies::Array{Float64,1}) is
        returned, both empty if no perturber of gMultiplet can be reached at all.

        AND THE ENERGIES ARE ALSO THE POLE POSITIONS, which is the second reason this function exists: a root search
        on a difference of polarizabilities has to know where the poles are, because a sign change ACROSS one is not
        a zero -- see StarkShift.magicWavelengths.

        THE REBUILD IS EXACT RATHER THAN AN APPROXIMATION, and that is measured instead of asserted: against a full
        sum-over-states evaluation it is BIT-IDENTICAL at omega = 0, 1e-3, 1e-2, 0.1 and 0.3 a.u. on the approved
        H 1s case (tools/probe-blackbodyShift.jl), the last being a frequency at which alpha_0 has already grown
        from 3.49 to 9.18 a.u.
"""
function polarizabilityTerms(level::Level, gMultiplet::Multiplet, grid::Radial.Grid)
    strengths = Float64[];      energies = Float64[]
    if  isempty(gMultiplet.levels)    return( (strengths, energies) )    end
    _, contributions = MultipolePolarizibility.computeMultipolePolarizabilities(E1, level, gMultiplet, grid)
    for  c in contributions
        push!( strengths, MultipolePolarizibility.rankCoefficient(1, 0, level.J, c.nuLevel.J) * abs2(c.reducedMe) )
        push!( energies,  -c.deltaE )             ## = W_J' - W_J
    end

    return( (strengths, energies) )
end


"""
`StarkShift.rebuiltPolarizability(strengths::Array{Float64,1}, energies::Array{Float64,1}, omega::Float64)`  
    ... returns the rank-0 E1 polarizability alpha_0(omega), in a.u., from the pieces that
        StarkShift.polarizabilityTerms extracted once, on the REAL frequency axis. It is the one place where those
        pieces are recombined, so that a quadrature and a root search cannot drift apart in how they do it. A
        value::Float64 is returned, and 0. for an empty set of terms.
"""
function rebuiltPolarizability(strengths::Array{Float64,1}, energies::Array{Float64,1}, omega::Float64)
    if  isempty(strengths)    return( 0.0 )    end

    return( sum( strengths[k] * MultipolePolarizibility.frequencyFactor(energies[k], omega, false)
                 for k = 1:length(strengths) ) )
end


"""
`StarkShift.computeBlackbodyShift(level::Level, gMultiplet::Multiplet, grid::Radial.Grid, temperature::Float64;
                                     nPoints::Int64=400, xMax::Float64=40.0)`  
    ... computes the black-body radiation shift of `level` at `temperature` (in Kelvin) by averaging the rank-0
        E1 polarizability over the Planck spectrum,

            DeltaE_BBR = -(1/2) <E^2> alpha_eff      alpha_eff = int x^3 alpha_0(Tx)/(e^x-1) dx / int x^3/(e^x-1) dx

        where x = omega/T and BOTH integrals are taken by the SAME quadrature on the SAME nodes -- see the comment
        at the loop, which is the one thing here that must not be tidied -- so that alpha_eff -> alpha_0(0) for a
        flat polarizability to machine precision rather than to the quadrature error, and the leading term is the
        universally quoted -(2/15)(alpha_fs pi)^3 T^4 alpha_0. There is NO second factor of two here, unlike the
        ac Stark shift: <E^2> is already a mean square. Only the RANK-0 polarizability enters, because thermal
        radiation is isotropic and unpolarized -- see the module docstring. A tuple
        (bbrShift::EmProperty, eta::Float64) is returned, with eta = alpha_eff/alpha_0(0) - 1.

        THE SUM OVER STATES IS EVALUATED ONCE, AT omega = 0, and the quadrature over omega is then arithmetic on
        its own stored terms -- the reduced matrix elements and energies that MultipolePolarizibility.Contribution
        carries -- recombined through MultipolePolarizibility.rankCoefficient and .frequencyFactor. That is the
        idiom module-DispersionCoefficient.jl uses for its Casimir-Polder integral, and it matters twice here: a
        few hundred quadrature points would otherwise mean a few hundred sum-over-states evaluations, and
        MultipolePolarizibility.checkFrequency would fire once per point, which is a guard that gets scrolled past.

        eta IS SMALL FOR A LOW-LYING LEVEL AND GROWS AS T^2, being first order in <omega^2>/dE^2 with <omega^2>
        proportional to T^2. Measured on the approved H 1s case: +1.15e-04 at 300 K, +1.28e-03 at 1000 K and
        +1.18e-02 at 3000 K. The law is far sharper on H 2p, whose smaller gap makes eta bigger, and that is the
        case the test suite asserts: the successive ratios at 1, 3, 10, 30, 100 and 300 K come out 9.000, 11.111,
        9.001, 11.120 and 9.062 against the exact (T2/T1)^2 of 9, 11.111, 9, 11.111 and 9 -- holding from
        eta = 3.8e-08 upwards, over two and a half decades. Where eta is NOT small the spectrum has reached the
        level's own excitation energies and StarkShift.checkBlackbodyResonance says so.
"""
function computeBlackbodyShift(level::Level, gMultiplet::Multiplet, grid::Radial.Grid, temperature::Float64;
                               nPoints::Int64=400, xMax::Float64=40.0)
    if  temperature <= 0.0  ||  isempty(gMultiplet.levels)    return( (EmProperty(0.), 0.) )    end
    T = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temperature)

    strengths, energies = StarkShift.polarizabilityTerms(level, gMultiplet, grid)
    if  isempty(strengths)    return( (EmProperty(0.), 0.) )    end
    alpha0 = StarkShift.rebuiltPolarizability(strengths, energies, 0.0)
    if  alpha0 == 0.0    return( (EmProperty(0.), 0.) )    end

    StarkShift.checkBlackbodyResonance(temperature, energies, nPoints, xMax)

    ## THE PLANCK WEIGHT IS NORMALISED BY ITS OWN QUADRATURE, on the SAME nodes, and NOT by the analytic
    ## pi^4/15.  The two differ by the quadrature error -- +2.57e-07 relative at these 400 nodes -- and dividing
    ## a Simpson numerator by an exact denominator puts that error straight into the answer, giving eta a FLOOR
    ## that is not physics: measured on H 2p, eta came out 5.97e-07 at 3 K where the T^2 law demands 3.4e-07, and
    ## the 3 K -> 10 K ratio was 6.8 against the required 11.1.  Sharing the quadrature makes alphaEff an exact
    ## weighted average of alpha_0 over the nodes, so it is EXACTLY alpha_0 for a flat polarizability and eta
    ## measures the dynamic departure alone.  Do not "simplify" this back to pi^4/15.
    h = xMax/nPoints;      wa = 0.0;      wb = 0.0
    for  i = 0:nPoints
        x  = i*h
        gx = x == 0.0 ? 0.0 : x^3 / (exp(x) - 1.0)
        w  = (i == 0 || i == nPoints) ? 1.0 : (isodd(i) ? 4.0 : 2.0)
        wa = wa + w * gx * StarkShift.rebuiltPolarizability(strengths, energies, T*x)
        wb = wb + w * gx
    end
    alphaEff = wa / wb                    ## the h/3 of the Simpson rule cancels in the ratio
    shift    = -0.5 * StarkShift.blackbodyFieldSquare(temperature) * alphaEff

    return( (EmProperty(shift, shift), alphaEff/alpha0 - 1.0) )
end


"""
`StarkShift.magicWavelengths(levelA::Level, levelB::Level, gMultiplet::Multiplet, grid::Radial.Grid;
                                omegaMin::Float64=0.01, omegaMax::Float64=0.25, nSamples::Int64=40,
                                tolerance::Float64=1.0e-12, residualTolerance::Float64=1.0e-6)`  
    ... finds every magic wavelength of the pair (`levelA`, `levelB`) between `omegaMin` and `omegaMax`, i.e. every
        zero of alpha_0(omega, A) - alpha_0(omega, B), and returns them as an
        Array{StarkShift.MagicWavelength,1} in ascending frequency -- empty if the range holds none. At such a
        frequency the trapping light shifts both levels equally, so the TRANSITION between them is untouched to
        first order in the intensity, which is what lets an optical trap coexist with a clock or a qubit. One
        gMultiplet serves both levels: each one's own selection rules decide which perturbers it reaches, so
        unreachable ones drop out on their own.

        A SIGN CHANGE IS NOT A ZERO, AND THAT IS THE WHOLE DIFFICULTY. alpha_0(omega) runs to +infinity as omega
        approaches an excitation energy from below and comes back from -infinity above it, so the difference jumps
        sign across EVERY resonance of EITHER level. A bisection let loose on such a sign change converges
        contentedly onto the POLE and reports a resonance as a magic wavelength -- which is the one error that
        looks entirely plausible, since the answer is a sharp, reproducible frequency.

        THE FIRST HALF OF THE REMEDY IS THAT THE POLES ARE KNOWN EXACTLY: they are the excitation energies that
        StarkShift.polarizabilityTerms returns. The range is CUT at every one of them, and each bracket is taken
        strictly inside a pole-free subinterval.

        AND THAT IS NECESSARY BUT NOT SUFFICIENT, which was MEASURED rather than foreseen -- the first version of
        this function did only the partition and returned two poles as magic wavelengths. THE POLES COME IN
        NEAR-DEGENERATE PAIRS: in hydrogen 2p -> 3s and 2p -> 3d sit 1e-7 apart, and the fine-structure partners of
        any subshell closer still. The sliver between such a pair is below numerical resolution, so the sign change
        inside it is the PAIR OF POLES and not a zero, and bisection converges onto one of them. Measured on
        H 1s against H 2p: the two spurious roots came back with |alpha_A - alpha_B|/alpha of 1.2e+01 and 1.1e+01
        and slopes of 1e12 and 1e13, against 3.5e-09 and 7.8e+04 for the one genuine root.

        SO EVERY CANDIDATE IS VERIFIED by its residual: |alpha_A - alpha_B| must be below `residualTolerance` times
        the scale of the two polarizabilities there, or it is discarded and counted. The two populations are seven
        orders apart, so the default 1e-6 is nowhere near either of them.

        AND EACH SUBINTERVAL IS SAMPLED RATHER THAN MERELY TESTED AT ITS ENDS. Between a pole of A and a pole of B
        the difference runs from -infinity to -infinity, so it holds an EVEN number of zeros -- commonly two, and
        an end-point test finds neither. `nSamples` points per subinterval is the resolution at which two zeros
        can still be separated; raise it if a known magic wavelength is missed.

        WHAT COMES BACK IS NOT AUTOMATICALLY USABLE, and MagicWavelength carries the two reasons why: a zero with
        `alphaCommon` negative is a trap that REPELS, and one with `aboveResonance` true sits where the true
        polarizability is complex and this real sum has no width.
"""
function magicWavelengths(levelA::Level, levelB::Level, gMultiplet::Multiplet, grid::Radial.Grid;
                          omegaMin::Float64=0.01, omegaMax::Float64=0.25, nSamples::Int64=40,
                          tolerance::Float64=1.0e-12, residualTolerance::Float64=1.0e-6)
    results = StarkShift.MagicWavelength[]
    sA, eA  = StarkShift.polarizabilityTerms(levelA, gMultiplet, grid)
    sB, eB  = StarkShift.polarizabilityTerms(levelB, gMultiplet, grid)
    if  isempty(sA)  ||  isempty(sB)
        printstyled("\n>>> ONE OF THE TWO LEVELS REACHES NO PERTURBER of this gMultiplet at all -- a magic " *
                    "wavelength needs an E1 response\n    from BOTH, and a set of one parity serves only one of " *
                    "them. No zero can be found and none is returned.\n", color=:light_red)
        return( results )
    end
    if  omegaMax <= omegaMin    return( results )    end
    dAlpha = omega -> StarkShift.rebuiltPolarizability(sA, eA, omega) -
                      StarkShift.rebuiltPolarizability(sB, eB, omega)

    # the poles inside the range, from BOTH levels, are where the range must be cut
    poles    = sort( unique( [dE  for dE in vcat(eA, eB)  if  omegaMin < dE < omegaMax] ) )
    edges    = vcat( omegaMin, poles, omegaMax )
    positive = [dE  for dE in vcat(eA, eB)  if dE > 0.0]
    lowest   = isempty(positive) ? Inf : minimum(positive)

    # a candidate survives only if its RESIDUAL is small on the scale of the polarizabilities there;  that is what
    # separates a zero from a pole caught between two near-degenerate resonances
    nRejected = 0
    keep      = function(omega::Float64)
        aA = StarkShift.rebuiltPolarizability(sA, eA, omega)
        aB = StarkShift.rebuiltPolarizability(sB, eB, omega)
        if  abs(aA - aB) > residualTolerance * max(abs(aA), abs(aB), 1.0)
            nRejected = nRejected + 1
            return( nothing )
        end
        push!( results, StarkShift.makeMagicWavelength(omega, sA, eA, dAlpha, lowest) )
        return( nothing )
    end

    for  i = 1:length(edges)-1
        delta = 1.0e-6 * (edges[i+1] - edges[i])
        lo    = edges[i] + delta;        hi = edges[i+1] - delta
        hi <= lo   &&   continue
        step  = (hi - lo) / nSamples
        xPrev = lo;      fPrev = dAlpha(lo)
        for  k = 1:nSamples
            xNext = lo + k*step;     fNext = dAlpha(xNext)
            if      fPrev == 0.0            keep(xPrev)
            elseif  fPrev * fNext < 0.0
                # a sign change strictly inside a pole-free subinterval;  plain bisection, then verified
                a = xPrev;   b = xNext;   fa = fPrev
                while  b - a > tolerance
                    m = 0.5*(a + b);     fm = dAlpha(m)
                    fm == 0.0   &&   (a = b = m;  break)
                    if  fa * fm < 0.0    b = m    else    a = m;  fa = fm    end
                end
                keep( 0.5*(a + b) )
            end
            xPrev = xNext;      fPrev = fNext
        end
    end
    if  nRejected > 0
        println(">>> magicWavelengths: $nRejected candidate(s) discarded as POLES rather than zeros -- a sign " *
                "change between two\n    near-degenerate resonances, whose residual |alpha_A - alpha_B| was not " *
                "small. $(length(results)) zero(s) kept.")
    end

    return( results )
end


"""
`StarkShift.makeMagicWavelength(omega::Float64, strengths::Array{Float64,1}, energies::Array{Float64,1},
                                    dAlpha::Function, lowest::Float64)`  
    ... assembles the StarkShift.MagicWavelength record for a zero already located at `omega`: it converts the
        frequency to a vacuum wavelength in nm, reads off the polarizability the two levels share there from
        `strengths`/`energies` (either level will do, since that is what a zero of the difference means), takes
        d(alpha_A - alpha_B)/d omega by a central difference on `dAlpha`, and compares omega with `lowest`, the
        lowest excitation energy of either level. A mw::StarkShift.MagicWavelength is returned.
"""
function makeMagicWavelength(omega::Float64, strengths::Array{Float64,1}, energies::Array{Float64,1},
                             dAlpha::Function, lowest::Float64)
    h     = 1.0e-6 * omega
    slope = (dAlpha(omega + h) - dAlpha(omega - h)) / (2*h)

    return( StarkShift.MagicWavelength(omega,
                Defaults.convertUnits("energy: from atomic to Angstrom", omega) / 10.0,
                StarkShift.rebuiltPolarizability(strengths, energies, omega), slope, omega >= lowest) )
end


"""
`StarkShift.computeStarkShifts(level::Level, gMultiplet::Multiplet, grid::Radial.Grid, settings::StarkShift.Settings)`
    ... computes the Stark shift of every M-sublevel of `level` in a field of strength settings.EField [V/cm] and
        angular frequency settings.omega [a.u.], from the level's own E1 scalar/tensor polarizability (taken from
        MultipolePolarizibility rather than re-derived here -- see module docstring):

            DeltaE(J,M) = - c alpha_0(omega) E^2  -  c alpha_2(omega) E^2 [3M^2-J(J+1)] / [J(2J-1)]

        with c = 1/2 for a STATIC field (omega = 0.) and c = 1/4 for a field oscillating at omega, where the
        second case is averaged over the cycle and E is then the PEAK amplitude.  The tensor term is undefined
        and set to 0. for J < 1, matching MultipolePolarizibility's own alpha_2 guard.  A tuple
        (alpha0::EmProperty, alpha2::EmProperty, Jsublevels::Array{StarkShift.SublevelJ,1}) is returned.

        AT omega = 0. THE STATIC ENTRY POINT IS CALLED, not the general one with omega = 0. passed through.  The
        two agree bit for bit -- that is asserted in `testModule_MultipolePolarizibility` -- but keeping the
        static call makes it visible here that adding the dynamic case moved no static number.
"""
function computeStarkShifts(level::Level, gMultiplet::Multiplet, grid::Radial.Grid, settings::StarkShift.Settings)
    if  settings.omega == 0.0
        alpha0, alpha2, _ = MultipolePolarizibility.computeScalarTensorPolarizability(level, gMultiplet, grid)
        cycleFactor = 0.5                ## a constant field
    else
        alphas, _ = MultipolePolarizibility.computeMultipolePolarizabilities(E1, level, gMultiplet, grid;
                                                                              omega = settings.omega)
        alpha0 = EmProperty(0.);     alpha2 = EmProperty(0.)
        for  (K, alphaK) in alphas
            if      K == 0   alpha0 = alphaK
            elseif  K == 2   alpha2 = alphaK
            end
        end
        cycleFactor = 0.25               ## a field oscillating at omega, averaged over the cycle
    end
    Efield  = settings.EField / StarkShift.AU_EFIELD_IN_VCM
    Efield2 = Efield^2
    Jd      = Basics.twice(level.J) / 2.0

    Jsublevels = StarkShift.SublevelJ[]
    for  M in AngularMomentum.m_values(level.J)
        Md = Basics.twice(M) / 2.0
        scalarShift = (-cycleFactor) * alpha0 * Efield2
        if  Jd*(2*Jd-1) == 0.0
            tensorShift = EmProperty(0.)
        else
            tensorFactor = (3*Md^2 - Jd*(Jd+1)) / (Jd*(2*Jd-1))
            tensorShift  = (-cycleFactor*tensorFactor) * alpha2 * Efield2
        end
        push!(Jsublevels, StarkShift.SublevelJ(M, scalarShift + tensorShift))
    end
    return  (alpha0, alpha2, Jsublevels)
end


"""
`StarkShift.computeTheta1(level::Level, grid::Radial.Grid)`  
    ... to compute the electric-dipole-moment [Theta^(1)] for the given level; an edm::Float64 is returned.
"""
function  computeTheta1(level::Level, grid::Radial.Grid)
    print(">>> Calculate the Theta^(1) (shift) coefficient for level (index=$(level.index), J=$(level.J)) ...")
    M = AngularM64(level.J.num, level.J.den)
    Theta1 = MultipoleMoment.emmStaticAmplitude(1, level, level, grid; display=false)
             ## AngularMomentum.ClebschGordan(level.J, M, AngularJ64(1), AngularM64(0), level.J, M) / sqrt(Basics.twice(level.J) + 1)
    J      = level.J.num/level.J.den
    w3j    = AngularMomentum.Wigner_3j(J, 2, J, J, 0, -J)
    Theta1 = Theta1 * (-1)^(2*J - 1 + 1) * w3j * sqrt(2*J+1)
         
    Theta1JG = MultipoleMoment.amplitude(1, level, grid, display=false)
    println("   Theta^(1) (level=$(level.index) [J=$(level.J)$(string(level.parity))]) = $(Theta1)  $(Theta1JG)  e a_o")
    return( Theta1 )
end



"""
`StarkShift.computeTheta2(level::Level, grid::Radial.Grid)`  
    ... to compute the electric-quadrupole-moment [Theta^(2)] for the given level; an eqm::Float64 is returned.
        Until July 2024, it remains unclear why the Clebsch-Gordan coefficient does not give a proper result
"""
function  computeTheta2(level::Level, grid::Radial.Grid)
    print(">>> Calculate the Theta^(2) (shift) coefficient for level (index=$(level.index), J=$(level.J)) ...")
    M = AngularM64(level.J.num, level.J.den)
    Theta2 = MultipoleMoment.emmStaticAmplitude(2, level, level, grid; display=false) 
             ## AngularMomentum.ClebschGordan(level.J, M, AngularJ64(2), AngularM64(0), level.J, M) / sqrt(Basics.twice(level.J) + 1)
    # Use the conversion factor by Jan Gilles (2024)
    J      = level.J.num/level.J.den
    w3j    = AngularMomentum.Wigner_3j(J, 2, J, J, 0, -J)
    Theta2 = Theta2 * (-1)^(2*J - 2 + 1) * w3j * sqrt(2*J+1)

    Theta2JG = MultipoleMoment.amplitude(2, level, grid, display=false)
    println("   Theta^(2) (level=$(level.index) [J=$(level.J)$(string(level.parity))]) = $(Theta2)  $(Theta2JG)  e a_o^2")
    return( Theta2 )
end


"""
`StarkShift.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::StarkShift.Settings; output=true)`  
    ... to compute (as selected) the EDM and EQM factors for the levels of the given multiplet and as specified by the given 
        settings. The results are returned in table (if required) or nothing is returned otherwise.
        The nuclear model is not used at present.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::StarkShift.Settings; output=true)
    println("")
    printstyled("StarkShift.computeOutcomes(): The computation of the electric multipole moments starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    outcomes = StarkShift.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    StarkShift.displayOutcomes(outcomes)    end
    # Calculate all amplitudes and requested properties
    newOutcomes = StarkShift.Outcome[]
    for  outcome in outcomes
        newOutcome = StarkShift.computeAmplitudesProperties(outcome, grid, settings) 
        push!( newOutcomes, newOutcome)
    end
    # A magic wavelength belongs to a PAIR of levels and not to one, so it cannot be formed in
    # computeAmplitudesProperties, which sees a single outcome:  the partner is named by its index in THIS
    # multiplet and is looked up here, where the multiplet is in hand.
    if  settings.calcMagicWavelengths
        newOutcomes = StarkShift.addMagicWavelengths(newOutcomes, multiplet, grid, settings)
    end
    # Print all results to screen
    StarkShift.displayResults(stdout, newOutcomes, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    StarkShift.displayResults(iostream, newOutcomes, settings)   end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`StarkShift.addMagicWavelengths(outcomes::Array{StarkShift.Outcome,1}, multiplet::Multiplet, grid::Radial.Grid,
                                settings::StarkShift.Settings)`
    ... returns the given outcomes with their `magicWavelengths` filled in, each selected level paired with the
        level `settings.magicPartnerIndex` of the SAME multiplet. This exists separately from
        StarkShift.computeAmplitudesProperties for one reason: a magic wavelength is a property of a PAIR of
        levels, so it cannot be formed from a single outcome, and the partner has to be looked up where the
        multiplet is in hand. An Array{StarkShift.Outcome,1} is returned.

        THREE THINGS ARE REFUSED RATHER THAN RETURNED EMPTY, because an empty list reads as "this pair has no
        magic wavelength" -- a physical statement -- where the truth is that the computation was not set up:
        a `gMultiplet` with no levels (the sum-over-states has nothing to run over), a `magicPartnerIndex` that
        names no level of the multiplet, and a `magicOmegas` of any length other than 0 or 2.

    + outcomes    ::Array{StarkShift.Outcome,1}  ... the outcomes whose amplitudes are already computed.
    + multiplet   ::Multiplet                    ... the multiplet the partner level is taken from.
    + grid        ::Radial.Grid                  ... the radial grid.
    + settings    ::StarkShift.Settings          ... carries magicPartnerIndex and magicOmegas.
"""
function addMagicWavelengths(outcomes::Array{StarkShift.Outcome,1}, multiplet::Multiplet, grid::Radial.Grid,
                             settings::StarkShift.Settings)
    if  isempty(settings.gMultiplet.levels)
        error("StarkShift.addMagicWavelengths(): calcMagicWavelengths = true needs a gMultiplet of "      *
              "intermediate, opposite-parity levels for the alpha_0(omega) sum-over-states to run over;  "  *
              "none was given.  See the module docstring and examples/example-Cl.jl for how it is built.")
    end
    partners = [ lev  for lev in multiplet.levels  if  lev.index == settings.magicPartnerIndex ]
    if  length(partners) != 1
        error("StarkShift.addMagicWavelengths(): magicPartnerIndex = $(settings.magicPartnerIndex) names "   *
              "no single level of this multiplet, whose indices are "                                        *
              "$([lev.index  for lev in multiplet.levels]).  A magic wavelength belongs to a PAIR, so the "  *
              "partner must be named.")
    end
    partner = partners[1]
    if      length(settings.magicOmegas) == 0    omegaMin, omegaMax = 0.01, 0.25
    elseif  length(settings.magicOmegas) == 2    omegaMin, omegaMax = settings.magicOmegas[1], settings.magicOmegas[2]
    else
        error("StarkShift.addMagicWavelengths(): magicOmegas must be empty, for the documented defaults of "  *
              "StarkShift.magicWavelengths, or hold exactly [omegaMin, omegaMax];  it holds "                 *
              "$(length(settings.magicOmegas)) entries.")
    end

    newOutcomes = StarkShift.Outcome[]
    for  outcome in outcomes
        if  outcome.Jlevel.index == partner.index
            # the partner paired with itself has no magic wavelength;  it keeps an empty list
            push!( newOutcomes, outcome );      continue
        end
        mws = StarkShift.magicWavelengths(outcome.Jlevel, partner, settings.gMultiplet, grid;
                                          omegaMin=omegaMin, omegaMax=omegaMax)
        push!( newOutcomes, StarkShift.Outcome(outcome.Jlevel, outcome.Theta1, outcome.Theta2, outcome.alpha0,
                                               outcome.alpha2, outcome.Jsublevels, outcome.bbrShift,
                                               outcome.bbrEta, mws) )
    end

    return( newOutcomes )
end


"""
`StarkShift.determineOutcomes(multiplet::Multiplet, settings::StarkShift.Settings)`  
    ... to determine a list of Outcomes's for the computation of electric-multipole-moments and energy shifts for 
        levels from the given multiplet. It takes into account the particular selections and settings. 
        An Array{StarkShift.Outcome,1} is returned. Apart from the level specification, all physical properties are 
        still set to zero during the initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::StarkShift.Settings) 
    # Define values that depend on the requested computations
    outcomes   = StarkShift.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            Jsublevels = StarkShift.SublevelJ[];   Mvalues = AngularMomentum.m_values(level.J)
            for  M in Mvalues   push!(Jsublevels, StarkShift.SublevelJ(M, EmProperty(0.)) )    end
            push!( outcomes, StarkShift.Outcome(level, 0., 0., EmProperty(0.), EmProperty(0.), Jsublevels,
                                                 EmProperty(0.), 0., StarkShift.MagicWavelength[]) )
        end
    end
    return( outcomes )
end


"""
`StarkShift.displayOutcomes(outcomes::Array{StarkShift.Outcome,1})`  
    ... to display a list of levels that have been selected for the computations. A small neat table of all selected 
        levels and their energies is printed but nothing is returned otherwise.
"""
function  displayOutcomes(outcomes::Array{StarkShift.Outcome,1})
    nx = 43
    println(" ")
    println("  Selected levels for EMM and Stark-shift computations:")
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
`StarkShift.displayResults(stream::IO, outcomes::Array{StarkShift.Outcome,1}, settings::StarkShift.Settings)`  
    ... to display the energies, Lande factors, Zeeman amplitudes etc. for the selected levels. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{StarkShift.Outcome,1}, settings::StarkShift.Settings)
    #
    if  settings.calcEDM  ||  settings.calcEQM
        nx = 82
        println(stream, " ")
        println(stream, "  Atomic electric-multipole-moment (EMM) shift coefficients:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=5)              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(14, "Theta^(1)/EDM"; na=5)                          
        sb = sb * TableStrings.center(14, "[e a_o]"; na=5)        
        sa = sa * TableStrings.center(14, "Theta^(2)/EQM"; na=5)                          
        sb = sb * TableStrings.center(14, "[e a_o^2]"; na=5)        
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            energy = outcome.Jlevel.energy
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))      * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.Theta1) )              * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.Theta2) )              * "    "
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    if  settings.calcStarkshifts
        nx1 = 106
        println(stream, " ")
        if  settings.omega == 0.0
            println(stream, "  Static scalar/tensor polarizabilities entering the Stark shift (Coulomb and " *
                            "Babushkin gauge):")
        else
            println(stream, "  Scalar/tensor polarizabilities at omega = $(settings.omega) a.u. entering the ac " *
                            "Stark shift, i.e. the cycle-averaged shift of a field of PEAK amplitude EField " *
                            "(Coulomb and Babushkin gauge):")
        end
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx1))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(15, "alpha_0 [Coul]"; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
        sa = sa * TableStrings.center(15, "alpha_0 [Bab]" ; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
        sa = sa * TableStrings.center(15, "alpha_2 [Coul]"; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
        sa = sa * TableStrings.center(15, "alpha_2 [Bab]" ; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx1))
        #
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha0.Coulomb) )      * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha0.Babushkin) )    * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha2.Coulomb) )      * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha2.Babushkin) )    * "    "
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(nx1))
        #
        nx2 = 92
        println(stream, " ")
        println(stream, "  Stark-shifts (energy-shifts) of the atomic sublevels at E = $(settings.EField) V/cm:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx2))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(8,  "M";     na=4);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(17, "Delta-E [Coul]"; na=3);                    sb = sb * TableStrings.center(17, "[Hz]"; na=3)
        sa = sa * TableStrings.center(17, "Delta-E [Bab]" ; na=3);                    sb = sb * TableStrings.center(17, "[Hz]"; na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx2))
        #
        for  outcome in outcomes
            sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            for  Jsublevel in outcome.Jsublevels
                sa  = "  "
                sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
                sa = sa * TableStrings.center(10, string(sym); na=4)
                sa = sa * TableStrings.center(8,  string(Jsublevel.M); na=4)
                sa = sa * TableStrings.flushright(17, @sprintf("%.6e", Defaults.convertUnits("energy: from atomic to Hz", Jsublevel.energy.Coulomb)) )   * "    "
                sa = sa * TableStrings.flushright(17, @sprintf("%.6e", Defaults.convertUnits("energy: from atomic to Hz", Jsublevel.energy.Babushkin)) ) * "    "
                println(stream, sa )
            end
        end
        println(stream, "  ", TableStrings.hLine(nx2))
    end
    #
    if  settings.temperature > 0.0
        nx3 = 104
        println(stream, " ")
        println(stream, "  Black-body radiation shift at T = $(settings.temperature) K, i.e. the rank-0 " *
                        "polarizability averaged over the Planck")
        println(stream, "  spectrum. It carries NO M-dependence: thermal radiation is isotropic and unpolarized, " *
                        "so only rank 0 survives.")
        println(stream, "  eta is the dynamic correction -- how far alpha_0(omega) departs from alpha_0(0) across " *
                        "the spectrum:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx3))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                   sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                   sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(17, "Delta-E [Coul]"; na=3);          sb = sb * TableStrings.center(17, "[Hz]"; na=3)
        sa = sa * TableStrings.center(17, "Delta-E [Bab]";  na=3);          sb = sb * TableStrings.center(17, "[Hz]"; na=3)
        sa = sa * TableStrings.center(15, "eta";            na=3);          sb = sb * TableStrings.center(15, " "   ; na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx3))
        #
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * TableStrings.flushright(17, @sprintf("%.6e",
                          Defaults.convertUnits("energy: from atomic to Hz", outcome.bbrShift.Coulomb)) )   * "    "
            sa = sa * TableStrings.flushright(17, @sprintf("%.6e",
                          Defaults.convertUnits("energy: from atomic to Hz", outcome.bbrShift.Babushkin)) ) * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%+.3e", outcome.bbrEta) )                       * "    "
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(nx3))
    end
    #
    if  settings.calcMagicWavelengths  &&  any(oc -> !isempty(oc.magicWavelengths), outcomes)
        nx4 = 112
        println(stream, " ")
        println(stream, "  Magic wavelengths, each level paired with level $(settings.magicPartnerIndex):")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx4))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                   sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=3);                   sb = sb * TableStrings.hBlank(13)
        sa = sa * TableStrings.center(16, "omega [a.u.]"; na=2);            sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(16, "lambda [nm]";  na=2);            sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(16, "alpha_common"; na=2);            sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(14, "d(dalpha)/domega"; na=2)
        println(stream, sa);   println(stream, sb);   println(stream, "  ", TableStrings.hLine(nx4))
        #
        for  outcome in outcomes
            isempty(outcome.magicWavelengths)  &&  continue
            sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            for  (k, mw)  in  enumerate(outcome.magicWavelengths)
                sa = "  "
                sa = sa * TableStrings.center(10, k == 1 ? TableStrings.level(outcome.Jlevel.index) : ""; na=2)
                sa = sa * TableStrings.center(10, k == 1 ? string(sym) : ""; na=3)
                sa = sa * TableStrings.flushright(16, @sprintf("%.8e", mw.omega) )           * "  "
                sa = sa * TableStrings.flushright(16, @sprintf("%.6f",  mw.wavelength) )     * "  "
                sa = sa * TableStrings.flushright(16, @sprintf("%+.6e", mw.alphaCommon) )    * "  "
                sa = sa * TableStrings.flushright(14, @sprintf("%+.2e", mw.slope) )
                # the two reasons a zero may not be usable, carried by MagicWavelength itself
                if  mw.alphaCommon < 0.    sa = sa * "   [alpha < 0: this trap REPELS]"                        end
                if  mw.aboveResonance      sa = sa * "   [ABOVE a resonance: a principal value, no width]"     end
                println(stream, sa )
            end
        end
        println(stream, "  ", TableStrings.hLine(nx4))
    end
    #
    return( nothing )
end

end # module

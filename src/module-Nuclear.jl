
"""
`module JAC.Nuclear`  
... a submodel of JAC that contains procedures for defining a Nuclear.Model and for calculating various nuclear
    potentials.
"""
module Nuclear

using  QuadGK, ..Basics,  ..Defaults, ..Radial, ..Math
export Model, AbstractNuclearModel, PointNucleus, UniformNucleus, FermiNucleus, ThreeParameterFermiNucleus


"""
`abstract type  Nuclear.AbstractNuclearModel`
    ... defines an abstract type to distinguish the shapes of the nuclear charge distribution:

    + struct PointNucleus     ... a point-like nucleus, rho(r) = Z delta(r).
    + struct UniformNucleus   ... a homogeneously charged sphere of radius R.
    + struct FermiNucleus     ... a two-parameter Fermi distribution,
                                  rho(r) = rho_0 / (1 + exp((r-c)/a)), with the skin thickness a fixed by
                                  Nuclear.fermiA and c determined from the requested rms radius.
    + struct ThreeParameterFermiNucleus
                              ... a three-parameter Fermi distribution,
                                  rho(r) = rho_0 (1 + w r^2/c^2) / (1 + exp((r-c)/a)); see below.

        Until 12-Aug-2026 the model was carried as a String and compared with ==, which made every new
        distribution another branch and left the accepted spellings inconsistent ("Fermi", but "point" and
        "uniform").  Dispatching on the type instead follows Basics.AbstractScField and its singletons, and
        makes a further distribution one new subtype plus one new method of Nuclear.nuclearPotential.
"""
abstract type  AbstractNuclearModel                            end
struct     PointNucleus          <:  AbstractNuclearModel      end
struct     UniformNucleus        <:  AbstractNuclearModel      end
struct     FermiNucleus          <:  AbstractNuclearModel      end


"""
`struct  Nuclear.ThreeParameterFermiNucleus  <:  Nuclear.AbstractNuclearModel`
    ... a three-parameter Fermi distribution of the nuclear charge,

            rho(r) = rho_0 (1 + w r^2/c^2) / (1 + exp((r-c)/a)),

        which relaxes the flat interior that the two-parameter form imposes by construction.  Here w > 0
        describes a central bump and w < 0 a CENTRAL DEPRESSION.  Two quite different mechanisms produce
        the latter: in medium-mass nuclei it is a shell effect (the 34Si bubble), whereas in the superheavy
        region it is driven by the Coulomb repulsion pushing protons outward, which makes it robust rather
        than model-dependent.  This is the form tabulated by the elastic electron-scattering compilations
        (de Vries, de Jager and de Vries, At. Data Nucl. Data Tables 36 (1987) 495) and used in the analysis
        of muonic atoms, where the muon orbits largely INSIDE the nucleus and hence samples the shape of the
        charge distribution rather than only its second moment.

    + w        ::Float64    ... dimensionless shape parameter; w = 0 recovers Nuclear.FermiNucleus exactly.
    + a        ::Float64    ... skin-thickness parameter [fm].  Unlike the two-parameter model, which always
                                uses the fixed Nuclear.fermiA, this is a field: the published parameter sets
                                come as triples (w, c, a) with a different a for every nucleus, so a fixed a
                                would make the model impossible to compare against any tabulated nucleus.
                                The one-argument constructor defaults it to Nuclear.fermiA.

        NOTE on the parametrisation itself: for w < 0 the profile 1 + w r^2/c^2 turns negative beyond
        r = c/sqrt(|w|).  For the tabulated magnitudes (|w| <~ 0.1) that lies beyond 3c, where the Fermi
        factor has already fallen to ~exp(-2c/a), so the contribution is numerically irrelevant.  It is a
        property of the three-parameter form, however, and not an oversight here.
"""
struct  ThreeParameterFermiNucleus  <:  AbstractNuclearModel
    w        ::Float64
    a        ::Float64
end


"""
`Nuclear.ThreeParameterFermiNucleus(w::Float64)`
    ... constructor for a three-parameter Fermi nucleus with the shape parameter w and the skin thickness
        left at JAC's standard Nuclear.fermiA.
"""
ThreeParameterFermiNucleus(w::Float64) = ThreeParameterFermiNucleus(w, Nuclear.fermiA)


"""
`Nuclear.name(model::Nuclear.AbstractNuclearModel)`
    ... returns a short name::String of the given nuclear model, for printout.
"""
name(model::Nuclear.PointNucleus)                 = "Point"
name(model::Nuclear.UniformNucleus)               = "Uniform"
name(model::Nuclear.FermiNucleus)                 = "Fermi"
name(model::Nuclear.ThreeParameterFermiNucleus)   = "Fermi-3p"


"""
`struct  Nuclear.Model`  ... defines a type for the nuclear model, i.e. for its form and parameters.

    + Z        ::Float64         ... nuclear charge
    + model    ::AbstractNuclearModel
                                 ... shape of the nuclear charge distribution: PointNucleus(),
                                     UniformNucleus() or FermiNucleus().  This was a String until
                                     12-Aug-2026, with the spellings "Fermi", "point" and "uniform" -- a
                                     capitalisation the docstring itself got wrong, so that
                                     Nuclear.Model(Z, "Point") raised for anyone following the manual.
                                     A wrong model is now a MethodError at the call site instead.
    + mass     ::Float64         ... atomic mass
    + radius   ::Float64         ... (root-mean square) radius of a uniform or Fermi-distributed nucleus
    + spinI    ::AngularJ64      ... nuclear spin I, must be >= 0
    + mu       ::Float64         ... magnetic dipole moment [nuclear magnetons]
    + Q        ::Float64         ... electric quadrupole moment [barn]
    + Omega    ::Float64         ... magnetic octupole moment [nuclear magnetons x barn]
"""
struct  Model
    Z          ::Float64
    model      ::AbstractNuclearModel
    mass       ::Float64
    radius     ::Float64
    spinI      ::AngularJ64      
    mu         ::Float64          
    Q          ::Float64
    Omega      ::Float64
end


"""
`Nuclear.Model(Z::Real)`  
    ... to specify a Fermi-type nucleus with charge Z, and where the nuclear spin and nuclear moments are all set to zero.
"""
function Model(Z::Real)
    Z < 0.1  &&  error("Z must be >= 0.1")
    model    = FermiNucleus()
    mass     = 2*Z + 0.005*Z^2
    radius   = Nuclear.rrmsRadius(mass)
    spinI    = AngularJ64(0)
    mu       = 0.
    Q        = 0.
    Omega    = 0.

    Model(Z, model, mass, radius, spinI, mu, Q, Omega) 
end


"""
`Nuclear.Model(Z::Real, M::Float64)`  
    ... to specify a Fermi-type nucleus with charge Z and mass M, and where the nuclear spin and nuclear moments are all set to zero.
"""
function Model(Z::Real, M::Float64)
    Z < 0.1  &&  error("Z must be >= 0.1")
    model    = FermiNucleus()
    mass     = M
    radius   = Nuclear.rrmsRadius(mass)
    spinI    = AngularJ64(0)
    mu       = 0.
    Q        = 0.
    Omega    = 0.

    Model(Z, model, mass, radius, spinI, mu, Q, Omega) 
end


"""
`Nuclear.Model(Z::Real, model::Nuclear.AbstractNuclearModel)`  
    ... to specify a nucleus with charge Z and the given model -- PointNucleus(), UniformNucleus() or
        FermiNucleus() -- and where the nuclear spin and nuclear moments are all set to zero.
"""
function Model(Z::Real, model::AbstractNuclearModel)
    Z < 0.1  &&  error("Z must be >= 0.1")
    if  model isa PointNucleus
        mass     = 1000*Z
        radius   = 0.
    else
        mass     = 2*Z + 0.005*Z^2
        radius   = Nuclear.rrmsRadius(mass)
    end
    spinI    = AngularJ64(0)
    mu       = 0.
    Q        = 0.
    Omega    = 0.

    Model(Z, model, mass, radius, spinI, mu, Q, Omega) 
end


"""
`Nuclear.Model(nm::Nuclear.Model;`
    
            Z=..,         model=..,         mass=..,        radius=..,     
            spinI=..,     mu=..,            Q=..,           Omega=..)
    ... constructor for re-defining a nuclear model nm::Nuclear.Model.
"""
function Model(nm::Nuclear.Model;            Z::Union{Nothing,Float64}=nothing,          model::Union{Nothing,AbstractNuclearModel}=nothing,         
    mass::Union{Nothing,Float64}=nothing,    radius::Union{Nothing,Float64}=nothing,     spinI::Union{Nothing,AngularJ64}=nothing,  
    mu::Union{Nothing,Float64}=nothing,      Q::Union{Nothing,Float64}=nothing,          Omega::Union{Nothing,Float64}=nothing)

    if  isnothing(Z)           Zx          = nm.Z           else   Zx          = Z          end 
    if  isnothing(model)       modelx      = nm.model       else   modelx      = model      end 
    if  isnothing(mass)        massx       = nm.mass        else   massx       = mass       end 
    if  isnothing(radius)      radiusx     = nm.radius      else   radiusx     = radius     end 
    if  isnothing(spinI)       spinIx      = nm.spinI       else   spinIx      = spinI      end 
    if  isnothing(mu)          mux         = nm.mu          else   mux         = mu         end 
    if  isnothing(Q)           Qx          = nm.Q           else   Qx          = Q          end 
    if  isnothing(Omega)       Omegax      = nm.Omega       else   Omegax      = Omega      end 
    
    Model(Zx, modelx, massx, radiusx, spinIx, mux, Qx, Omegax)
end


# `Base.show(io::IO, m::Model)`  ... prepares a proper printout of the variable  m::Model.
function Base.show(io::IO, m::Model) 
    print(io, "$(Nuclear.name(m.model)) nuclear model for Z = $(m.Z) with mass = $(m.mass), " *
              "radius R = $(m.radius) fm and ")
    print(io, "nuclear spin I = $(m.spinI), dipole moment mu = $(m.mu), quadrupole moment Q = $(m.Q) and octupole moment Omega = $(m.Omega).")
end

        
"""
`struct  Nuclear.Isomer`  
    ... defines a type for modeling isomeric levels that are involved in hyperfine-induced transitions and structure.
        It assumes that the nuclear charge, model and mass, radius is defined by an associated nm::Nuclear.Model.

    + spinI         ::AngularJ64   ... nuclear spin I >= 0 of the isomeric nuclear level, could be the ground level.
    + parity        ::Parity       ... parity of the isomeric nuclear level
    + energy        ::Float64      ... nuclear excitation energy of the isomeric level; 0. if nuclear ground level [in user-specified units]
    + mu            ::Float64      ... magnetic dipole moment [Bohr magnetons].
    + Q             ::Float64      ... electric quadrupole moment.
    + Omega         ::Float64      ... magnetic octupole moment [Bohr magnetons x fm^2].
    + multipoleM    ::Array{EmMultipole,1}
        ... multipoles of the <Ia || M^(multipole) || Ib > nuclear matrix elements.
    + elementM      ::Array{Float64,1}
        ... (real) values of the <Ia || M^(multipole) || Ib > nuclear matrix elements in [a.u.], one for each entry
            of multipoleM and IN THE SAME ORDER; the two arrays must therefore have equal length.

        These last two fields carry the NUCLEAR transition data: they are what makes a hyperfine-induced nuclear
        transition possible at all, and what mixes two nuclear levels of the same F. Use
        Nuclear.reducedTransitionAmplitude to obtain elementM from a reduced transition probability B(multipole),
        rather than converting by hand -- that conversion carries a Weisskopf unit and a convention factor, and
        has already been got wrong once in an application script.
"""
struct  Isomer
    spinI           ::AngularJ64
    parity          ::Parity  
    energy          ::Float64 
    mu              ::Float64
    Q               ::Float64
    Omega           ::Float64
    multipoleM      ::Array{EmMultipole,1}
    elementM        ::Array{Float64,1}  
end


"""
`Nuclear.Isomer()`  ... constructor for an `empty` instance of Nuclear.Isomer.

        Note (06-Aug-2026): multipoleM and elementM are both left EMPTY. They previously read `[E1]` and
        `Float64[]`, i.e. of unequal length, which breaks the pairing the two fields are supposed to express and
        would give an out-of-bounds the moment anything iterated them together. An empty isomer has no known
        nuclear transition, which is the honest default.
"""
function Isomer()
    Isomer( AngularJ64(0), Basics.plus, 0., 0., 0., 0., EmMultipole[], Float64[])
end


"""
`Nuclear.Isomer(isomer::Nuclear.Isomer;`

            spinI=..,        parity=..,       energy=..,       mu=..,
            Q=..,            Omega=..,        multipoleM=..,   elementM=..)

    ... constructor for re-defining an isomer::Nuclear.Isomer, i.e. the standard JAC keyword copy-constructor.
        Its absence is why every application script builds an Isomer from positional arguments, and why every one
        of them is now wrong against the eight-field struct.
"""
function Isomer(isomer::Nuclear.Isomer;
    spinI::Union{Nothing,AngularJ64}=nothing,               parity::Union{Nothing,Parity}=nothing,
    energy::Union{Nothing,Float64}=nothing,                 mu::Union{Nothing,Float64}=nothing,
    Q::Union{Nothing,Float64}=nothing,                      Omega::Union{Nothing,Float64}=nothing,
    multipoleM::Union{Nothing,Array{EmMultipole,1}}=nothing, elementM::Union{Nothing,Array{Float64,1}}=nothing)

    if  isnothing(spinI)        spinIx      = isomer.spinI       else   spinIx      = spinI       end
    if  isnothing(parity)       parityx     = isomer.parity      else   parityx     = parity      end
    if  isnothing(energy)       energyx     = isomer.energy      else   energyx     = energy      end
    if  isnothing(mu)           mux         = isomer.mu          else   mux         = mu          end
    if  isnothing(Q)            Qx          = isomer.Q           else   Qx          = Q           end
    if  isnothing(Omega)        Omegax      = isomer.Omega       else   Omegax      = Omega       end
    if  isnothing(multipoleM)   multipoleMx = isomer.multipoleM  else   multipoleMx = multipoleM  end
    if  isnothing(elementM)     elementMx   = isomer.elementM    else   elementMx   = elementM    end
    #
    if  length(multipoleMx) != length(elementMx)
        error("\n\nNuclear.Isomer():  STOP -- multipoleM and elementM must have the same length, but are \n"  *
              "$multipoleMx  and  $elementMx.\nEach nuclear matrix element belongs to exactly one multipole.\n")
    end
    Isomer(spinIx, parityx, energyx, mux, Qx, Omegax, multipoleMx, elementMx)
end


"""
`Nuclear.weisskopfUnit(mp::EmMultipole, A::Int64)`
    ... to return one Weisskopf unit for the multipole mp and mass number A, in ATOMIC units, i.e. the single-particle
        estimate against which reduced nuclear transition probabilities B(mp) are conventionally quoted. A
        value::Float64 is returned.

        With R = 1.2 A^(1/3) fm,

            B_W(EL) = 1/(4 pi) * [3/(L+3)]^2 * R^(2L)      [e^2 fm^(2L)]
            B_W(ML) = 10/pi   * [3/(L+3)]^2 * R^(2L-2)     [mu_N^2 fm^(2L-2)]

        so that B_W(M1) = 1.79 mu_N^2, independent of A, as it must be. The result is converted into atomic units
        here, so that a caller never has to mix fm with a_0 or nuclear magnetons with atomic ones -- which is
        precisely where the factors have gone missing before.
"""
function  weisskopfUnit(mp::EmMultipole, A::Int64)
    L      = mp.L
    Rfm    = 1.2 * A^(1/3)
    fmToAu = Defaults.convertUnits("length: from fm to atomic", 1.0)
    if      mp.electric
        wa = 1.0/(4pi) * (3.0/(L+3))^2 * Rfm^(2L)
        return( wa * fmToAu^(2L) )                       ## e^2 a_0^(2L) = a.u. for EL
    else
        wa = 10.0/pi  * (3.0/(L+3))^2 * Rfm^(2L-2)
        muNToAu = Defaults.convertUnits("moment: from nuclear magneton to atomic", 1.0)
        return( wa * muNToAu^2 * fmToAu^(2L-2) )         ## mu_N^2 a_0^(2L-2)
    end
end


"""
`Nuclear.reducedTransitionAmplitude(mp::EmMultipole, B::Float64, A::Int64, spinI::AngularJ64;
                            inWeisskopfUnits::Bool=true)`
    ... to convert a reduced nuclear transition probability B(mp) into the reduced nuclear matrix element
        <I_f || M^(mp) || I_i> that Nuclear.Isomer.elementM expects, in atomic units; spinI is the spin of the
        INITIAL (decaying) nuclear state. A value::Float64 is returned. If inWeisskopfUnits, B is taken in
        Weisskopf units and converted with Nuclear.weisskopfUnit; otherwise B is taken as already in atomic units.

        THE CONVENTION, corrected 06-Aug-2026, and it is worth reading before changing anything here. The
        reduced transition probability is DEFINED as

            B(sigma L; I_i -> I_f)  =  |<I_f || M^(sigma L) || I_i>|^2 / (2 I_i + 1)

        so the reduced matrix element is simply sqrt( (2 I_i + 1) * B ), which is what is implemented.

        THIS IS NOT WHAT THE MANUSCRIPT SAYS. examples/papers/b26.pra-hyperfine-induced.tex gives, for the 235U
        E3 case, <I_g || W^(E3) || I_e> = sqrt( 8 pi / 7 * B(E3) ), i.e. sqrt( 8 pi/(2L+1) * B ), and that form
        was implemented here first. It is wrong for this purpose, and not merely by a constant: since B already
        carries the 1/(2 I_i + 1), the manuscript's W has the initial-state degeneracy folded into it and is
        therefore NOT a reduced matrix element. Two things then break.

          * elementM is contracted with a 6-j symbol in Hfs.computeInteractionAmplitudeM, where nothing but a
            true reduced matrix element is admissible.
          * the radiation field factor in HyperfineInduced.amplitude acquired a spurious (2 I_i + 1) dependence
            to compensate -- and a field factor cannot depend on the spin of the emitting state, since it
            describes the photon, not the nucleus. That is what exposed the error.

        The manuscript needs a corresponding correction; it is marked in cyan in the -claude.tex copy.

        UNITS. The returned matrix element is in JAC's own moment convention, i.e. magnetic multipoles carry
        mu_N = 5.446170e-4/2 a.u. ("moment: from nuclear magneton to atomic"), the SAME convention as
        Nuclear.Isomer.mu and hence as the diagonal hyperfine matrix elements. It is deliberately NOT converted
        to the Gaussian convention mu_B = alpha/2 that the multipole RADIATION formula needs: that factor of
        alpha belongs to the 1/c of the magnetic multipole radiation operator and is applied in
        HyperfineInduced.amplitude, not to the nuclear moment itself.

        Doing this by hand has already gone wrong: apps/apps-wuwang-hfs-induced/job-a-uranium.jl computes
        `sqrt(8/7. * 5.38e-27)` -- without the pi, and from a B value that differs from the one in the text.
"""
function  reducedTransitionAmplitude(mp::EmMultipole, B::Float64, A::Int64, spinI::AngularJ64;
                                     inWeisskopfUnits::Bool=true)
    Bau = inWeisskopfUnits  ?  B * Nuclear.weisskopfUnit(mp, A)  :  B
    if  Bau < 0.   error("\n\nNuclear.reducedTransitionAmplitude():  STOP -- a reduced transition probability "  *
                         "cannot be negative, but B = $B was given.\n")   end
    return( sqrt( (Basics.twice(spinI) + 1) * Bau ) )
end


# `Base.show(io::IO, isomer::Nuclear.Isomer)`  ... prepares a proper printout of the variable isomer:Nuclear.Isomer.
function Base.show(io::IO, isomer::Nuclear.Isomer) 
    println(io, "spinI:          $(isomer.spinI)  ")
    println(io, "parity:         $(isomer.parity)  ")
    println(io, "energy:         $(isomer.energy)  ")
    println(io, "mu:             $(isomer.mu)  ")
    println(io, "Q:              $(isomer.Q)  ")
    println(io, "Omega:          $(isomer.Omega)  ")
    println(io, "multipoleM:     $(isomer.multipoleM)  ")
    println(io, "elementM:       $(isomer.elementM)  ")
end

        
#################################################################################################################################
#################################################################################################################################



"""
`Nuclear.fermiA`  ... provides a value::Float64 for the fermi_a parameter.
"""
const fermiA = 2.3/(4 * log(3))     ## made const 12-Aug-2026: it is read inside the r_rho/rr_rho
                                    ## integrands, where a non-const global is boxed at every call.


"""
`Nuclear.rrmsRadius(A::Float64)`  
    ... provides a value::Float64 for the root-mean-squared radius R [fm] of a uniformly-distributed charge
        density of a nucleus with mass number A.  Renamed from Nuclear.Rrms on 12-Aug-2026 to follow JAC's
        camelCase rule for function names.
"""
function rrmsRadius(A::Float64)
    return( 0.836 * A^(1/3) + 0.57 )
end


"""
`Nuclear.fermiRrms(b::Float64)`  
    ... provides a value::Float64 for the root-mean-squared radius of a fermi-distributed charge 
        density for given fermi_a and fermi_b parameters.
"""
function fermiRrms(b::Float64)

    return( sqrt(12 * fermiA^2 * Math.polylogExp(b/fermiA, 5) / Math.polylogExp(b/fermiA, 3)) )
end


"""
`Nuclear.threeParameterFermiRrms(c::Float64, a::Float64, w::Float64)`
    ... provides a value::Float64 for the root-mean-square radius [fm] of the three-parameter Fermi charge
        density rho(r) = rho_0 (1 + w r^2/c^2) / (1 + exp((r-c)/a)), in closed form.

        The moments of the Fermi factor are analytic,

            M_n = Int_0^Inf r^n / (1 + exp((r-c)/a)) dr = - n! a^(n+1) Li_(n+1)(-exp(c/a)),

        and Math.polylogExp(x, s) is exactly Li_s(-exp(x)) (it wraps GSL's complete Fermi-Dirac integral).
        With the weight 1 + w r^2/c^2 folded in,

            <r^2> = [M_4 + (w/c^2) M_6] / [M_2 + (w/c^2) M_4],

        which is the expression below.  For w = 0 it collapses identically to 12 a^2 L_5 / L_3, i.e. to
        Nuclear.fermiRrms, and that is used as a regression test of this whole model.
"""
function threeParameterFermiRrms(c::Float64, a::Float64, w::Float64)
    L3 = Math.polylogExp(c/a, 3);   L5 = Math.polylogExp(c/a, 5);   L7 = Math.polylogExp(c/a, 7)
    ## u collects the dimensionless combination in which w always enters; M_(n+2)/M_n carries a^2 (n+1)(n+2)
    u   = w * a^2 / c^2
    num = 24 * L5 + 720 * u * L7
    den =  2 * L3 +  24 * u * L5
    ## The L_k are negative, so num and den are both negative for w = 0 and the ratio is positive.  Once
    ## |u| = |w| a^2/c^2 grows to O(1) -- i.e. once c falls towards a -- the w-term can outweigh the leading
    ## one and den crosses zero; the "distribution" is then no longer a sensible charge density and <r^2> is
    ## meaningless rather than merely inaccurate.  Say so instead of returning a NaN.
    if  den == 0.  ||  num/den <= 0.
        error("Nuclear.threeParameterFermiRrms(): the three-parameter Fermi form is not meaningful for " *
              "c = $c fm, a = $a fm, w = $w -- the weight (1 + w r^2/c^2) outweighs the Fermi factor " *
              "(|w| a^2/c^2 = $(abs(u)) is not small) and the second moment ceases to be positive.")
    end
    return( sqrt( a^2 * num / den ) )
end


"""
`Nuclear.computeThreeParameterFermiC(R::Float64, a::Float64, w::Float64)`
    ... computes a value::Float64 for the half-density radius c [fm] of a three-parameter Fermi nucleus with
        root-mean-square radius R, skin thickness a and shape parameter w.

        This BISECTS rather than iterating the fixed-point step used by Nuclear.computeFermiBParameter.  That
        step, b <- b - (fermiRrms(b) - R) with unit damping, converges at a rate set by d(fermiRrms)/db, and
        that derivative vanishes as R approaches the smallest radius the form can represent -- which is why
        it stalls on helium (see the comment in computeFermiBParameter).  Bisection cannot stall: it either
        converges or reports an unreachable R, which is the distinction that matters to the caller.

        The bracket is grown outwards from c0 = sqrt(5/3) R -- the sharp-sphere relation R = sqrt(3/5) c --
        rather than being fixed in advance, because threeParameterFermiRrms is monotone in c ONLY while
        c stays comfortably above a.  That was measured, not assumed: for a = 0.523 fm the closed form is
        smooth and increasing for c >~ 2.5 fm, and for smaller c it is non-monotone and eventually singular,
        since |w| a^2/c^2 then reaches O(1).  cFloor below marks that boundary.  Every physically meaningful
        nucleus has c ~ 1.2 R with R >~ 2 fm and so sits far above it.
"""
function computeThreeParameterFermiC(R::Float64, a::Float64, w::Float64)
    R <= 0.  &&  error("Nuclear.computeThreeParameterFermiC(): R_rms = $R fm must be positive.")
    cFloor = 4.0 * a          ## below this the w-weight competes with the Fermi factor; see the docstring
    c0     = sqrt(5/3) * R    ## sharp-sphere estimate, always slightly above the true c
    ##
    ## Grow the bracket outwards from c0 until it straddles R.
    cUpp = max(c0, cFloor);   nUpp = 0
    while  threeParameterFermiRrms(cUpp, a, w) < R   &&   nUpp < 60
        cUpp = 1.5 * cUpp;    nUpp = nUpp + 1
    end
    cLow = cUpp
    while  threeParameterFermiRrms(cLow, a, w) > R
        cLow = cLow / 1.5
        if  cLow < cFloor
            error("Nuclear.computeThreeParameterFermiC(): R_rms = $R fm is not representable by a three-" *
                  "parameter Fermi distribution with a = $a fm and w = $w.  The half-density radius would " *
                  "have to fall below $cFloor fm, where the weight (1 + w r^2/c^2) competes with the Fermi " *
                  "factor itself; use UniformNucleus() for such a small radius.")
        end
    end
    ##
    tolAccept = 1.0e-12
    for  i = 1:200
        cMid = 0.5 * (cLow + cUpp)
        if      cUpp - cLow < tolAccept                                 return( cMid )
        elseif  threeParameterFermiRrms(cMid, a, w) > R                 cUpp = cMid
        else                                                            cLow = cMid
        end
    end
    return( 0.5 * (cLow + cUpp) )
end


"""
`Nuclear.computeFermiBParameter(R::Float64)`  
    ... computes a value::Float64 for the fermi_b parameter for a fermi-distributed nuclear with 
        root-mean square radius R.
"""
function computeFermiBParameter(R::Float64)

    function eq(b :: Float64)
        return fermiRrms(b) - R
    end

    b = R
    diff = eq(b)

    eps = 1E-14
    count = 0
    max = 100

    while  abs(diff) > eps   &&   count < max
        b = b -  diff
        
        diff = eq(b)
        count = count + 1
    end
    ## Signal non-convergence rather than returning whatever the last iterate happened to be (12-Aug-2026).
    ## The loop above exits on EITHER the target eps or the iteration cap, and the caller could not tell the
    ## two apart, so a b-parameter the fixed-point step never resolved was handed back silently and used to
    ## build a potential.  The iteration itself is unchanged.
    ##
    ## THE TEST IS AGAINST A PHYSICALLY MEANINGFUL RESIDUAL, NOT AGAINST eps, and the distinction is not
    ## academic: helium (R_rms = 1.8993 fm) stalls at a residual of 5.9e-8 fm and never reaches eps = 1e-14
    ## at all.  That was found the first time this check ran, and it is worth knowing -- the loop's stated
    ## target has never been attainable there.  The step is b <- b - (fermiRrms(b) - R) with unit damping,
    ## and its convergence rate follows d(fermiRrms)/db, which vanishes as R_rms approaches the smallest
    ## radius the 2-parameter form can represent (~1.86 fm); helium sits just above that edge.  A residual
    ## of 6e-8 fm is 3e-8 relative, far below any experimental nuclear radius, so such a result is perfectly
    ## usable and must not be rejected.  tolAccept is what separates "stalled but fine" from "wrong".
    tolAccept = 1.0e-6
    if  abs(diff) > tolAccept
        error("Nuclear.computeFermiBParameter(): no convergence for R_rms = $R fm -- the fixed-point " *
              "iteration stopped after $count steps with a residual of $diff fm, above the acceptable " *
              "$tolAccept fm.  The 2-parameter Fermi form cannot represent an arbitrarily small radius " *
              "(its floor is about 1.86 fm); for a smaller R_rms use the uniform nuclear model instead.")
    end
    
    return( b )
end


"""
`Nuclear.fermiDistributedNucleus(Rrms::Float64, Z::Float64, grid::Radial.Grid)`  
    ... computes the effective, radial-dependent charge Z(r) for a Fermi-distributed nucleus with rms 
        radius R and nuclear charge Z. The full nuclear potential is then given by V_nuc = - Z(r)/r; 
        a potential::Radial.Potential is returned.
"""
function fermiDistributedNucleus(Rrms::Float64, Z::Float64, grid::Radial.Grid)

    zznew = zeros(Float64, grid.NoPoints)

    if  Z < 1.2
        # With JAC's fixed Fermi skin-thickness parameter fermiA, the 2-parameter Fermi charge distribution
        # cannot represent an rms radius below about 1.86 fm (computeFermiBParameter diverges to negative b
        # for smaller targets) -- far above the physical charge radius of hydrogen (~0.88 fm). An earlier
        # version of this code silently substituted b = computeFermiBParameter(1.89) here, ignoring whatever
        # Rrms was actually requested; that silent substitution has been removed (25-Jul-2026) in favor of
        # this explicit error, since it was giving quietly wrong isotope-shift/HFS results for Z=1.
        error("The Fermi nuclear model cannot represent Z = $Z (only reachable in practice for hydrogen, " *
              "Z=1): with JAC's fixed Fermi skin-thickness parameter (fermiA = $(round(fermiA, digits=4)) fm), " *
              "the 2-parameter Fermi charge distribution has a minimum representable rms radius of about " *
              "1.86 fm, far above hydrogen's physical charge radius (~0.88 fm). Use " *
              "Nuclear.Model(Z, UniformNucleus(), mass, Rrms, ...) instead for Z=1.")
    else
        b = computeFermiBParameter(Rrms);    b < 0  &&  error("Inappropriate R_rms radius.")
    end
    
    ## The two integrands are defined AFTER fermiA_au/fermiC_au (moved here 12-Aug-2026).  They used to be
    ## written above, closing over variables that were only assigned thirty lines further down: legal Julia,
    ## since a closure captures the variable rather than its value, but it boxes both captures and reads as
    ## an error.  Same values, same order of evaluation.
    fermiA_au = Defaults.convertUnits("length: from fm to atomic", fermiA)
    fermiC_au = Defaults.convertUnits("length: from fm to atomic", b)
    function  r_rho(r::Float64)   r   / (1.0 + exp( (r-fermiC_au)/fermiA_au ) )   end
    function  rr_rho(r::Float64)  r^2 / (1.0 + exp( (r-fermiC_au)/fermiA_au ) )   end
    N = 1. / QuadGK.quadgk(rr_rho, 0., 1.3, rtol=1.0e-10)[1]
    #
    #
    for  i = 2:length(zznew)   
        zznew[i] = 1 / grid.r[i] * quadgk(rr_rho, 0.0, 1.0e-3, grid.r[i], rtol=1.0e-8)[1] 
        zznew[i] = zznew[i] + quadgk(r_rho, grid.r[i], 1.0, 1.0e+1, 1.0e+3, rtol=1.0e-8)[1]
        zznew[i] = Z * N * zznew[i] * grid.r[i]
    end 

    potential = Radial.Potential("nuclear-potential: Fermi-distributed", zznew, deepcopy(grid))
    return( potential )
end


"""
`Nuclear.threeParameterFermiNucleus(Rrms::Float64, Z::Float64, model::Nuclear.ThreeParameterFermiNucleus, grid::Radial.Grid)`
    ... computes the effective, radial-dependent charge Z(r) for a nucleus whose charge density follows the
        three-parameter Fermi form rho(r) = rho_0 (1 + w r^2/c^2) / (1 + exp((r-c)/a)) with rms radius Rrms
        and nuclear charge Z.  The full nuclear potential is then given by V_nuc = - Z(r)/r;
        a potential::Radial.Potential is returned.

        The structure follows Nuclear.fermiDistributedNucleus exactly -- the same split of the Coulomb
        integral into an inner and an outer part, the same quadrature tolerances -- with both integrands and
        the normalisation carrying the extra weight (1 + w r^2/c^2).  For w = 0 the two therefore agree to
        the quadrature tolerance, which is checked as a regression test.
"""
function threeParameterFermiNucleus(Rrms::Float64, Z::Float64, model::Nuclear.ThreeParameterFermiNucleus, grid::Radial.Grid)

    zznew = zeros(Float64, grid.NoPoints)
    c     = computeThreeParameterFermiC(Rrms, model.a, model.w)

    fermiA_au = Defaults.convertUnits("length: from fm to atomic", model.a)
    fermiC_au = Defaults.convertUnits("length: from fm to atomic", c)
    w         = model.w
    ## The shape weight; for w < 0 it turns negative beyond r = c/sqrt(|w|), which for the tabulated
    ## magnitudes lies beyond 3c where the Fermi factor has already fallen to ~exp(-2c/a).  See the docstring
    ## of Nuclear.ThreeParameterFermiNucleus: this belongs to the parametrisation, and is not clamped here.
    function  shape(r::Float64)   1.0 + w * (r/fermiC_au)^2                                    end
    function  r_rho(r::Float64)   r   * shape(r) / (1.0 + exp( (r-fermiC_au)/fermiA_au ) )     end
    function  rr_rho(r::Float64)  r^2 * shape(r) / (1.0 + exp( (r-fermiC_au)/fermiA_au ) )     end
    N = 1. / QuadGK.quadgk(rr_rho, 0., 1.3, rtol=1.0e-10)[1]
    #
    for  i = 2:length(zznew)
        zznew[i] = 1 / grid.r[i] * quadgk(rr_rho, 0.0, 1.0e-3, grid.r[i], rtol=1.0e-8)[1]
        zznew[i] = zznew[i] + quadgk(r_rho, grid.r[i], 1.0, 1.0e+1, 1.0e+3, rtol=1.0e-8)[1]
        zznew[i] = Z * N * zznew[i] * grid.r[i]
    end

    potential = Radial.Potential("nuclear-potential: three-parameter Fermi", zznew, deepcopy(grid))
    return( potential )
end


"""
`Nuclear.pointNucleus(Z::Float64, grid::Radial.Grid)`
    ... computes the effective, radial-dependent charge Z(r) for a point-like nucleus with nuclear charge Z. 
        The full nuclear potential is then given by V_nuc = - Z(r)/r; a potential::Radial.Potential is returned.
"""
function pointNucleus(Z::Float64, grid::Radial.Grid)

    zz = Z * ones(Float64, grid.NoPoints)

    potential = Radial.Potential("nuclear-potential: point-like", zz, deepcopy(grid))
    return( potential )
end


"""
`Nuclear.uniformNucleus(R::Float64, Z::Float64, grid::Radial.Grid)`  
    ... computes the effective, radial-dependent charge Z(r) for a uniformly-distributed nucleus with radius 
        R [in fm] and nuclear charge Z. The full nuclear potential is then given by V_nuc = - Z(r)/r outside 
        the nucleus; a potential::Radial.Potential is returned.
"""
function uniformNucleus(R::Float64, Z::Float64, grid::Radial.Grid)

    zz = zeros(Float64, grid.NoPoints);   R_au = Defaults.convertUnits("length: from fm to atomic", R)

    for i = 1:grid.NoPoints
        if     grid.r[i] <= R_au    zz[i] = Z / (2 * R_au) * (3. - grid.r[i]^2/R_au^2) * grid.r[i]
        ## if     grid.r[i] <= R_au    zz[i] = Z / (2 * R_au) * (grid.r[i]^2/R_au^2 - 3.) * grid.r[i] * grid.r[i]
        else   zz[i] = Z
        end
    end

    potential = Radial.Potential("nuclear-potential: uniformly-distributed", zz, deepcopy(grid))
    return( potential )
end


"""
`Nuclear.nuclearPotential(nm::Nuclear.Model, grid::Radial.Grid)`  
    ... computes the effective, radial-dependent charge Z(r) for a nucleus with radius R and nuclear charge Z. 
        The full nuclear potential is then given by V_nuc = - Z(r)/r; a potential::Radial.Potential 
        is returned.
"""
function nuclearPotential(nm::Nuclear.Model, grid::Radial.Grid)
    return( Nuclear.nuclearPotential(nm.model, nm, grid) )
end

nuclearPotential(::Nuclear.PointNucleus,   nm::Nuclear.Model, grid::Radial.Grid) =
    Nuclear.pointNucleus(nm.Z, grid)
nuclearPotential(::Nuclear.UniformNucleus, nm::Nuclear.Model, grid::Radial.Grid) =
    Nuclear.uniformNucleus(nm.radius, nm.Z, grid)
nuclearPotential(::Nuclear.FermiNucleus,   nm::Nuclear.Model, grid::Radial.Grid) =
    Nuclear.fermiDistributedNucleus(nm.radius, nm.Z, grid)
nuclearPotential(model::Nuclear.ThreeParameterFermiNucleus, nm::Nuclear.Model, grid::Radial.Grid) =
    Nuclear.threeParameterFermiNucleus(nm.radius, nm.Z, model, grid)


"""
`Nuclear.nuclearPotentialDH(nm::Nuclear.Model, grid::Radial.Grid, lambda::Float64)`  
    ... computes the effective, radial-dependent charge Z(r) for a nucleus with radius R, nuclear charge Z 
        and for a Debye-Hueckel screening exp(-lambda r). The full Debye-Hueckel nuclear potential is then 
        given by V_nuc = - Z(r)/r; a potential::Radial.Potential is returned.
""" 
function nuclearPotentialDH(nm::Nuclear.Model, grid::Radial.Grid, lambda::Float64)
    pot = nuclearPotential(nm, grid)
    Zr  = pot.Zr
    for  i = 1:length(Zr)   Zr[i] = Zr[i] *exp(-lambda * pot.grid.r[i])    end
    
    potential = Radial.Potential(pot.name* "+ Debey-Hueckel screening", Zr, deepcopy(pot.grid))
    return( potential )
end

end # module


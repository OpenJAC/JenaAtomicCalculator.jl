#
# Validation of the BLACK-BODY RADIATION SHIFT added to module-StarkShift.jl on 20-Oct-2026.
#
# The quantity:  a level in isotropic thermal radiation at temperature T is shifted by
#
#     DeltaE_BBR = -(1/2) <E^2> alpha_0(0) (1 + eta),      <E^2> = (4 pi^3/15) alpha_fs^3 T^4
#
# whose leading term is the universally quoted -(2/15)(alpha_fs pi)^3 T^4 alpha_0.  THERE IS NO SECOND FACTOR
# OF TWO here, unlike the ac Stark shift in the same module: <E^2> is already a mean square, whereas the ac
# line is written with the PEAK amplitude and picks up <cos^2> = 1/2.  Getting that wrong is this quantity's
# classic error, so it is worth saying twice.
#
# ONLY THE RANK-0 POLARIZABILITY ENTERS.  Thermal radiation is isotropic and unpolarized, so its field
# correlation tensor is <E_i E_j> = (1/3) delta_ij <E^2>: the rank-1 component averages away for want of a net
# handedness, the rank-2 for want of a preferred axis.  That is why a black-body shift is ONE number per level
# and not a table over M -- and why clock papers quote one number.
#
# CASE A  <E^2> BY TWO ROUTES.  The closed form against a quadrature over JAC's own
#         Distribution.photonNumberDensity(Distribution.PhotonPlanck(T), omega), which is what fixes the Planck
#         normalisation used in the module.  THE 300 K VALUE IS A COMPLETE EXTERNAL ANCHOR and costs nothing to
#         check, because no atomic quantity enters it -- only the fine-structure constant does.
# CASE B  THE DYNAMIC CORRECTION eta on the approved H 1s case, at three temperatures, both by an independent
#         quadrature written out here and through StarkShift.computeBlackbodyShift.  eta must grow as T^2.
# CASE C  THE RESONANCE GUARD.  It must stay silent for a low-lying level at room temperature and speak when
#         the Planck spectrum reaches the level's own excitation energies.
#
# RESULTS, 20-Oct-2026 (this file is the validation record):
#
#   CASE A   <E^2>(300 K) = 2.61752213e-18 a.u., i.e. E_rms = 831.945 V/m against the TEXTBOOK 831.9 V/m.  The
#            quadrature over Distribution.photonNumberDensity gives the same to 1.3e-10.  The implied leading
#            shift for a level of alpha_0 = 1 a.u. is -0.008611 Hz at 300 K, and -(1/2)<E^2> and
#            -(2/15)(alpha_fs pi)^3 T^4 agree to all six printed digits, as they must, being the same expression.
#   CASE B   H 1s, perturbers 2p/3p/4p in the nuclear field (the approved case, alpha_0 = 3.491829):
#                 300 K   DeltaE = -0.030072 Hz     eta = +1.153e-04
#                1000 K   DeltaE = -3.716979 Hz     eta = +1.282e-03
#                3000 K   DeltaE = -304.241942 Hz   eta = +1.181e-02
#            The module reproduces the independent quadrature DIGIT FOR DIGIT at all three.  And eta follows T^2
#            to three digits -- 1.153e-04 * (10/3)^2 = 1.28e-03 against the measured 1.282e-03 -- which is what
#            it must do, since eta is first order in <omega^2>/dE^2 and <omega^2> goes as T^2.  That is the
#            check that the WEIGHT in the quadrature is right and not merely its normalisation.
#   CASE C   silent at 300 K.  At 50000 K (kT = 0.158 a.u. against a lowest excitation energy of 0.375 a.u.) it
#            reports 74.4 % of the Planck weight above that energy and the nearest quadrature point 1.34 % away
#            from it, and says the number must not be quoted.  eta there is -8.97e-01, i.e. the "correction" has
#            eaten the quantity -- exactly the shape of failure the guard exists to label.
#
# WHAT IS NOT COMPUTED, AND IT IS THE HALF THE RYDBERG REGIME NEEDS.  Above the lowest excitation energy
# alpha_0(omega) is COMPLEX.  The real part of the integral is then a PRINCIPAL VALUE, which a plain quadrature
# does not give, and the imaginary part is the BBR-INDUCED TRANSITION RATE, which is a separate observable and
# is not computed anywhere in JAC.  For a Rydberg level the excitation energies fall as 1/n^3 and the Planck
# weight sits ON them, so BOTH are needed there and NEITHER is available -- the guard is what keeps that from
# being discovered by reading a wrong number.
#
# Run as:  julia --project=. tools/probe-blackbodyShift.jl
#

using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf

const MP = MultipolePolarizibility
const SS = StarkShift

const AU_EFIELD_IN_VM = 5.14220675e11      ## atomic unit of electric field, in V/m


# <E^2> by QUADRATURE over JAC's own Planck photon number density:  <E^2> = int 4 pi omega n(omega) domega,
# written out here so that the module's closed form is checked against something independent of it
function fieldSquareByQuadrature(temperature::Float64; nx::Int64=4000, xMax::Float64=60.0)
    T = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temperature)
    h = xMax/nx;    wa = 0.0
    for  i = 0:nx
        x  = i*h;    om = T*x
        fx = x == 0.0 ? 0.0 :
             4.0*pi * om * Distribution.photonNumberDensity(Distribution.PhotonPlanck(T), om)
        wa = wa + ((i == 0 || i == nx) ? 1.0 : (isodd(i) ? 4.0 : 2.0)) * fx
    end

    return( T * h/3.0 * wa )
end


# the black-body shift by an INDEPENDENT quadrature, taking alpha_0(omega) from the module one point at a time
function shiftByQuadrature(level::Level, gMultiplet::Multiplet, grid::Radial.Grid, temperature::Float64;
                           nx::Int64=400, xMax::Float64=40.0)
    T  = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temperature)
    a0 = om -> redirect_stdout(devnull) do
             alphas, _ = MP.computeMultipolePolarizabilities(E1, level, gMultiplet, grid; omega=om)
             alphas[1][2].Babushkin   end
    h  = xMax/nx;    wa = 0.0
    for  i = 0:nx
        x  = i*h
        fx = x == 0.0 ? 0.0 : x^3 * a0(T*x) / (exp(x) - 1.0)
        wa = wa + ((i == 0 || i == nx) ? 1.0 : (isodd(i) ? 4.0 : 2.0)) * fx
    end
    alphaEff = (h/3.0 * wa) / (pi^4/15.0)

    return( (-0.5 * SS.blackbodyFieldSquare(temperature) * alphaEff, alphaEff/a0(0.0) - 1.0) )
end


function caseA()
    println("\nCASE A   <E^2> BY TWO ROUTES, and the 300 K value is a complete external anchor.")
    alphaFs = 1.0 / Defaults.getDefaults("speed of light: c")
    @printf("\n   %-10s %-20s %-14s %-20s %-12s\n", "T [K]", "<E^2> closed [a.u.]", "E_rms [V/m]", "<E^2> quadrature",
            "rel. diff")
    for  TK in [300.0, 1000.0, 3000.0]
        ec = SS.blackbodyFieldSquare(TK);     eq = fieldSquareByQuadrature(TK)
        @printf("   %-10.0f %-20.8e %-14.3f %-20.8e %-12.2e\n", TK, ec, sqrt(ec)*AU_EFIELD_IN_VM, eq,
                abs(eq-ec)/ec)
    end
    println("\n   THE TEXTBOOK VALUE AT 300 K IS 831.9 V/m.")
    T300 = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 300.0)
    @printf("\n   for alpha_0 = 1 a.u. at 300 K:  -(1/2)<E^2> = %.6f Hz    -(2/15)(alpha_fs pi)^3 T^4 = %.6f Hz\n",
            Defaults.convertUnits("energy: from atomic to Hz", -0.5*SS.blackbodyFieldSquare(300.0)),
            Defaults.convertUnits("energy: from atomic to Hz", -2.0/15.0*(alphaFs*pi)^3*T300^4))

    return( nothing )
end


function approvedH1s()
    ni   = Nuclear.Model(1.0, PointNucleus())
    grid = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=4.0e-2, rbox=80.0)
    scf  = Basics.NuclearField();    asf = AsfSettings(AsfSettings(); scField=scf)
    gMp  = redirect_stdout(devnull) do
        generate(Representation("np perturbers", ni, grid,
                   [Configuration("2p"), Configuration("3p"), Configuration("4p")],
                   MeanFieldMultiplet(MeanFieldSettings(scf))), output=true)[ResultKeys.MeanFieldMultiplet]   end
    lev  = redirect_stdout(devnull) do
        perform(Atomic.Computation(Atomic.Computation(), name="H 1s", grid=grid, nuclearModel=ni,
                   configs=[Configuration("1s")], asfSettings=asf); output=true)[ResultKeys.Multiplet].levels[1]  end

    return( (lev, gMp, grid) )
end


function caseB(lev, gMp, grid)
    println("\nCASE B   THE DYNAMIC CORRECTION eta on the approved H 1s case, by two independent quadratures.")
    @printf("\n   %-10s %-20s %-14s %-20s %-14s\n", "T [K]", "DeltaE module [Hz]", "eta module",
            "DeltaE probe [Hz]", "eta probe")
    etas = Float64[]
    for  TK in [300.0, 1000.0, 3000.0]
        sh, et = redirect_stdout(devnull) do;  SS.computeBlackbodyShift(lev, gMp, grid, TK)  end
        sq, eq = shiftByQuadrature(lev, gMp, grid, TK)
        push!( etas, et )
        @printf("   %-10.0f %-20.6f %-14.3e %-20.6f %-14.3e\n", TK,
                Defaults.convertUnits("energy: from atomic to Hz", sh.Babushkin), et,
                Defaults.convertUnits("energy: from atomic to Hz", sq), eq)
    end
    @printf("\n   the T^2 law:  eta(300) * (1000/300)^2 = %.3e against the measured eta(1000) = %.3e\n",
            etas[1]*(1000.0/300.0)^2, etas[2])
    @printf("                 eta(1000) * (3000/1000)^2 = %.3e against the measured eta(3000) = %.3e\n",
            etas[2]*9.0, etas[3])

    return( nothing )
end


function caseC(lev, gMp, grid)
    println("\nCASE C   THE RESONANCE GUARD.  H 1s has its lowest excitation energy at 0.375 a.u.")
    println("\n   --- at 300 K (kT = 9.5e-04 a.u.):  the guard must stay silent")
    sh, et = redirect_stdout(devnull) do;   SS.computeBlackbodyShift(lev, gMp, grid, 300.0)   end
    @printf("   silent;  DeltaE = %+.6f Hz,  eta = %+.3e\n",
            Defaults.convertUnits("energy: from atomic to Hz", sh.Babushkin), et)
    println("\n   --- at 50000 K (kT = 0.158 a.u.):  the guard must speak")
    sh2, et2 = SS.computeBlackbodyShift(lev, gMp, grid, 50000.0)
    @printf("\n   eta = %+.3e -- the 'correction' has eaten the quantity, which is the point of the warning\n", et2)

    return( nothing )
end


function main()
    println("\n", "="^118)
    println("THE BLACK-BODY RADIATION SHIFT -- validation record, 20-Oct-2026")
    println("="^118)
    caseA()
    lev, gMp, grid = approvedH1s()
    caseB(lev, gMp, grid)
    caseC(lev, gMp, grid)
    println("\n", "="^118)

    return( nothing )
end

main()

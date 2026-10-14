using Printf

"""
`module  JAC.ExportAdas`
    ... a submodel of JAC that writes an `ExportDataSet.DataSet` as a BASIC ADAS adf04 file (ADAS User
        Manual, Appendix A/ADF04, "resolved specific ion data collections", fetched 13-Oct-2026).

        SCOPE, DELIBERATELY NARROW. Only the "basic file structure" the manual itself shows (its Table B4c,
        a hydrogen example): a single parent, LS-style level labels, ITYP=3 ("Upsilons as a function of
        T_e", the usual form -- and exactly the quantity `ImpactExcitation.computeEffStrengths` already
        computes), and plain electron-impact excitation data (CCODE blank). This writer does NOT attempt
        dielectronic recombination, charge-exchange or electron-impact-ionisation transition lines (CCODE
        r/R, h/H, i/I/s/S), unresolved dielectronic satellites (CCODE l/L), multi-parent ionisation data,
        Eissner configuration notation, or the optional error block -- the adf04 specification covers all
        of these, and a caller needing them should write them by hand following the same manual appendix.
        This module exists to show the shape of the integration (JAC's own computed levels, A-values and
        effective collision strengths map onto a real external format with no new physics required), not to
        be a complete ADAS writer.

        COLUMN LAYOUT, HONESTLY STATED. Field widths below are matched to the manual's own worked example
        (Table B4c) rather than independently re-derived from its abstract Fortran format-code listing,
        which did not extract cleanly from the source PDF. The output LOOKS like a real adf04 file and
        carries the right numbers in the right order; byte-exact compliance with the `xxdata_04` reading
        routine's column positions has not been verified against that routine itself.
"""
module ExportAdas

using Printf, ..ExportDataSet


"""
`ExportAdas.formatExponential(x::Float64)`
    ... to format a number the way adf04's data lines do: a 2-decimal mantissa, a sign, and a 2-digit
        exponent, with NO 'e' or 'd' letter (e.g. 6.27e8 -> "6.27+08", 0.26 -> "2.60-01"), matching the
        `fm.n` / `e8.2`-without-letter convention the ADAS manual describes for STRG4/STRG5.
        A value::String of exactly 7 characters is returned.
"""
function  formatExponential(x::Float64)
    x == 0.0   &&   return( "0.00+00" )
    expo = floor(Int64, log10(abs(x)))
    mant = x / 10.0^expo
    if   round(abs(mant), digits=2) >= 10.0    mant = mant / 10.0;   expo = expo + 1   end
    sgn  = expo >= 0 ? "+" : "-"

    return( @sprintf("%4.2f%s%02d", mant, sgn, abs(expo)) )
end


"""
`ExportAdas.writeAdf04(stream::IO, data::ExportDataSet.DataSet)`
    ... to write `data` as a BASIC adf04 file onto `stream` -- see this module's own docstring for exactly
        what "basic" excludes. Nothing is returned.
"""
function  writeAdf04(stream::IO, data::ExportDataSet.DataSet)
    sym = data.element * ("+" ^ data.ionCharge)
    @printf(stream, "%-3s%2d%10d%10d%12.1f\n", sym, data.ionCharge, data.nuclearCharge, data.ionCharge+1,
                     data.ionizationPotential)

    for  lev in data.levels
        @printf(stream, "%5d  %-16s(%1d)%1d( %3.1f)%16.1f\n",
                lev.index, lev.label, lev.multiplicity, lev.L, lev.J, lev.energy)
    end
    println(stream, "   -1")

    @printf(stream, "%5.2f%5d      ", 1.0, 3)
    for  temp in data.temperatures    print(stream, " ", formatExponential(temp))    end
    println(stream)

    for  tr in data.transitions
        @printf(stream, "%4d %4d %s", tr.upperIndex, tr.lowerIndex, formatExponential(tr.aValue))
        for  ups in tr.upsilons    print(stream, " ", formatExponential(ups))    end
        println(stream)
    end
    println(stream, "   -1")

    return( nothing )
end

end # module

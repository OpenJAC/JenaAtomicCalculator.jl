using Printf

"""
`module  JAC.ExportChianti`
    ... a submodel of JAC that writes an `ExportDataSet.DataSet` as CHIANTI-style `.elvlc` (energy levels)
        and `.wgfa` (radiative data) files, in the classic (CHIANTI 1 through 7.1) fixed-column layout.

        SCOPE, DELIBERATELY NARROW, the same way `ExportAdas` is. Only `.elvlc` and `.wgfa` are written;
        `.scups`/`.splups` (electron-impact collisional data) is NOT, because CHIANTI stores that data
        Burgess & Tully (1992) SCALED -- a transition-type-dependent reduced-temperature transformation on
        top of the raw Upsilon(T) this module's own `TransitionRecord`s already carry, not a reformatting
        of it. This module exists to show the shape of the integration, not to be a complete CHIANTI writer;
        see [[project_challenges_list]] S24 for what remains.

        COLUMN LAYOUT, HONESTLY STATED, same caveat as `ExportAdas`: field widths follow the documented
        Ch7 column-position table (fetched 14-Oct-2026), not independent verification against a real
        ChiantiPy/fiasco parse.

        THE gf (WEIGHTED OSCILLATOR STRENGTH) COLUMN is DERIVED, not supplied: gf = g_upper * A * lambda^2
        / 6.670e15 (lambda in Angstrom, A in s^-1), from the standard A-f relation (IUPAC Gold Book) with
        gf = g_lower*f_absorption and A the emission rate. Cross-checked once (see `example-Ra.jl`) against
        the well-known He-like "w-line" oscillator strength (f ~ 0.7 across the isoelectronic sequence).
"""
module ExportChianti

using Printf, ..ExportDataSet

const  RYDBERG_PER_CM = 1.0 / 109737.31568160


"""
`ExportChianti.writeElvlc(stream::IO, data::ExportDataSet.DataSet)`
    ... to write `data`'s level list as a CHIANTI-style `.elvlc` file onto `stream`. JAC's own levels are
        THEORETICAL, not measured, so the observed-energy columns are written as -1.0, CHIANTI's own
        convention for "not available" (see this module's own docstring). Nothing is returned.
"""
function  writeElvlc(stream::IO, data::ExportDataSet.DataSet)
    for  lev in data.levels
        eRy = lev.energy * RYDBERG_PER_CM
        @printf(stream, "%3d %-22s%2d%3d%3.1f%6d%15.3f%15.6f%15.3f%15.6f\n",
                lev.index, lev.label, lev.multiplicity, lev.L, lev.J, round(Int64, 2*lev.J+1),
                -1.0, -1.0, lev.energy, eRy)
    end
    println(stream, " -1")

    return( nothing )
end


"""
`ExportChianti.writeWgfa(stream::IO, data::ExportDataSet.DataSet)`
    ... to write `data`'s transition list as a CHIANTI-style `.wgfa` file onto `stream`, deriving the gf
        column as described in this module's own docstring. Nothing is returned.
"""
function  writeWgfa(stream::IO, data::ExportDataSet.DataSet)
    levelOf = Dict(lev.index => lev for lev in data.levels)

    for  tr in data.transitions
        upper = levelOf[tr.upperIndex];   lower = levelOf[tr.lowerIndex]
        deltaE = upper.energy - lower.energy
        deltaE <= 0.0   &&   continue
        lambdaAngstrom = 1.0e8 / deltaE
        gUpper         = 2*upper.J + 1
        gf             = gUpper * tr.aValue * lambdaAngstrom^2 / 6.670e15

        @printf(stream, "%5d%5d%15.3f%13.4e%15.4e\n", tr.lowerIndex, tr.upperIndex, lambdaAngstrom, gf, tr.aValue)
    end
    println(stream, " -1")

    return( nothing )
end

end # module

using Printf

"""
`module  JAC.NistCalibration`
    ... a submodel of JAC that replaces selected level energies of a `ManyElectron.Multiplet` by externally
        supplied (typically NIST ASD) measured values, producing a `CalibratedMultiplet` -- the SAME CSFs,
        orbitals and CI mixing coefficients as the input, only the `.energy` field of matched levels changed.
        No process module (`Einstein`, `PhotoEmission`, `DielectronicRecombination`, `AutoIonization`, ...)
        ever uses or is used by this module: the dependency runs one way only. A caller who wants calibrated
        rates or cross sections extracts `cm.multiplet` (an ORDINARY `Multiplet`, indistinguishable in type
        from any other) and passes it to whichever property-computation function they would have used anyway
        -- `Einstein.computeLines(cm.multiplet, grid, settings)`, `IsotopeShift.computeOutcomes(cm.multiplet,
        ...)`, and so on, with no change needed in any of those functions.

        WHY A WRAPPER TYPE, NOT A BARE Multiplet. Returning a bare `Multiplet` would make it impossible to
        tell, from the type alone, whether a given multiplet was calibrated or not, and would lose the
        provenance (which NIST snapshot, which levels matched) the moment the caller passed it on.
        `CalibratedMultiplet` keeps that information attached and inspectable for as long as the caller keeps
        the wrapper around; `cm.multiplet` carries none of it once extracted, which is why `NistCalibration.
        summary(cm)` exists -- to make citing the snapshot in a results report one line rather than a thing to
        remember from scratch.

        MATCHING IS BY LEVEL INDEX ONLY, NEVER BY ENERGY. This is deliberate, for two reasons. First, an
        automatic derivation of an LS-term label from an arbitrary (possibly jj-coupled, possibly heavy-
        element) JAC level is its own small research question -- the same reasoning `ExportDataSet` is built
        on. Second, and more importantly, matching by index rather than by energy proximity is what makes
        re-applying calibration to an already-calibrated `Multiplet` idempotent rather than a silent
        corruption: a second pass matches the same index to the same NIST entry and sets the same energy
        again, regardless of how far the first pass moved anything.

        DOUBLE CALIBRATION, TWO LAYERS OF DEFENCE. `calibrate` only accepts a plain `Multiplet` as its first
        argument; calling it on a `CalibratedMultiplet` directly hits an explicit method below that raises a
        clear error rather than silently re-shifting anything. If a caller works around that by extracting
        `cm.multiplet` and calibrating it again, the result is UNCHANGED from the first pass -- not corrupted
        -- because matching is index-based (see above) and the reference level's energy (see
        `groundLevelIndex` below) is not itself altered by a normal calibration unless explicitly given its
        own, necessarily-zero, NIST entry.

        THE REFERENCE LEVEL IS NAMED, NOT GUESSED. Early design considered auto-detecting "the ground level"
        as whichever input level has the lowest energy -- rejected, because if calibration ever changes which
        level that is, a second pass would anchor to a DIFFERENT level than the first, breaking the very
        idempotence this module depends on. `groundLevelIndex` is therefore a required argument: the caller
        names the level whose ORIGINAL (ab initio) energy anchors every NIST-relative energy supplied.
"""
module NistCalibration

using Printf, ..ManyElectron, ..Defaults


"""
`struct  NistCalibration.NistEntry`
    ... one measured level, supplied by the caller (never fetched live -- see this module's own docstring on
        reproducibility).

    + levelIndex   ::Int64    ... the JAC level index (Level.index) this entry refers to.
    + energy       ::Float64  ... measured excitation energy relative to the multiplet's OWN ground level
                                    (see `groundLevelIndex` in `calibrate`), in cm^-1 -- NIST ASD's own
                                    convention.
    + label        ::String   ... free text, e.g. a configuration/term label or an NIST ASD line reference.
"""
struct  NistEntry
    levelIndex     ::Int64
    energy         ::Float64
    label          ::String
end


"""
`struct  NistCalibration.CalibratedLevel`
    ... the calibration record for one level, matched or not.

    + index        ::Int64    ... the JAC level index.
    + matched      ::Bool     ... True if a NistEntry with this index was found (and not excluded).
    + rawEnergy    ::Float64  ... the ORIGINAL (ab initio) level energy, in JAC's own atomic units.
    + nistEnergy   ::Float64  ... the calibrated level energy, in atomic units; NaN if not matched.
    + label        ::String   ... the matching NistEntry's own label, or "" if not matched.
"""
struct  CalibratedLevel
    index          ::Int64
    matched        ::Bool
    rawEnergy      ::Float64
    nistEnergy     ::Float64
    label          ::String
end


"""
`struct  NistCalibration.CalibratedMultiplet`
    ... the result of `NistCalibration.calibrate`: a `Multiplet` with calibrated energies, together with the
        provenance that made it. See this module's own docstring for why this is a distinct type rather than
        a bare `Multiplet`.

    + multiplet      ::ManyElectron.Multiplet         ... same CSFs/orbitals/mixing as the input; `.energy`
                                                            replaced for matched levels.
    + levelRecords   ::Array{CalibratedLevel,1}        ... one record per level, matched or not.
    + nistSource     ::String                          ... a dated citation for the data used, e.g.
                                                            "NIST ASD, snapshot 2026-09-01".
"""
struct  CalibratedMultiplet
    multiplet        ::ManyElectron.Multiplet
    levelRecords     ::Array{CalibratedLevel,1}
    nistSource       ::String
end


"""
`NistCalibration.calibrate(multiplet::ManyElectron.Multiplet, nistTable::Array{NistCalibration.NistEntry,1},
                            nistSource::String; groundLevelIndex::Int64, excludeLevels::Array{Int64,1}=Int64[])`
    ... to replace the energy of every level in `multiplet` that has a matching entry in `nistTable` (by
        level index, not energy -- see this module's own docstring) by the measured value, anchored to
        `multiplet`'s own level `groundLevelIndex`'s ORIGINAL energy. Levels in `excludeLevels` are left
        ab initio even if a matching entry exists. Every level, matched or not, gets a `CalibratedLevel`
        record. A cm::CalibratedMultiplet is returned.
"""
function  calibrate(multiplet::ManyElectron.Multiplet, nistTable::Array{NistEntry,1}, nistSource::String;
                     groundLevelIndex::Int64, excludeLevels::Array{Int64,1}=Int64[])
    groundMatches = [lev for lev in multiplet.levels if lev.index == groundLevelIndex]
    isempty(groundMatches)   &&   error("NistCalibration.calibrate(): no level with index $groundLevelIndex " *
                                         "in the given multiplet.")
    groundEnergy = groundMatches[1].energy
    auPerKayser  = 1.0 / Defaults.convertUnits("energy: from atomic to Kayser", 1.0)

    entryFor = Dict(e.levelIndex => e for e in nistTable)

    newLevels = ManyElectron.Level[];   records = CalibratedLevel[]
    for  lev in multiplet.levels
        if  haskey(entryFor, lev.index)  &&  !(lev.index in excludeLevels)
            entry     = entryFor[lev.index]
            newEnergy = groundEnergy + entry.energy * auPerKayser
            push!( newLevels, ManyElectron.Level(lev.J, lev.M, lev.parity, lev.index, newEnergy, lev.relativeOcc,
                                                  lev.hasStateRep, lev.basis, lev.mc) )
            push!( records, CalibratedLevel(lev.index, true, lev.energy, newEnergy, entry.label) )
        else
            push!( newLevels, lev )
            push!( records, CalibratedLevel(lev.index, false, lev.energy, NaN, "") )
        end
    end

    return( CalibratedMultiplet(ManyElectron.Multiplet(multiplet.name, newLevels), records, nistSource) )
end


"""
`NistCalibration.calibrate(cm::NistCalibration.CalibratedMultiplet, args...; kwargs...)`
    ... REFUSES: re-calibrating an already-calibrated multiplet directly is not supported, so that the
        common accidental double-calibration is a clear, immediate error rather than a silent re-shift (see
        this module's own docstring). Nothing is returned; an error is always raised.
"""
function  calibrate(cm::CalibratedMultiplet, args...; kwargs...)
    error("NistCalibration.calibrate(): the argument is already a CalibratedMultiplet (source: " *
          "\"$(cm.nistSource)\"). Re-calibrating it directly is not supported -- start from the original, " *
          "uncalibrated Multiplet if you need to calibrate again, or use cm.multiplet explicitly if you " *
          "intend to re-apply calibration on purpose (this is idempotent, see the module docstring, but is " *
          "not done implicitly).")
end


"""
`NistCalibration.summary(cm::NistCalibration.CalibratedMultiplet)`
    ... to build a one-line, ready-to-print citation of `cm`'s provenance, for pasting into a results report
        header -- see this module's own docstring for why this cannot be carried automatically through a
        later rate or cross-section computation. A value::String is returned.
"""
function  summary(cm::CalibratedMultiplet)
    nMatched = count( r -> r.matched, cm.levelRecords );   nTotal = length(cm.levelRecords)

    return( "Levels calibrated against $(cm.nistSource): $nMatched of $nTotal matched." )
end


"""
`NistCalibration.displayResults(stream::IO, cm::NistCalibration.CalibratedMultiplet)`
    ... to print a table of `cm`'s level records (raw and calibrated energy side by side, never one in place
        of the other) onto `stream`. Nothing is returned.
"""
function  displayResults(stream::IO, cm::CalibratedMultiplet)
    nx = 90
    println(stream, " ");   println(stream, "  ", summary(cm));   println(stream, " ")
    println(stream, "  ", "-"^nx)
    println(stream, "  ", rpad("Level",8), rpad("Matched?",10), rpad("Raw energy [a.u.]",20),
                     rpad("NIST energy [a.u.]",20), "Label")
    println(stream, "  ", "-"^nx)
    for  r in cm.levelRecords
        nistStr = r.matched ? @sprintf("%.10f", r.nistEnergy) : "--"
        println(stream, "  ", rpad(string(r.index),8), rpad(r.matched ? "yes" : "no", 10),
                         rpad(@sprintf("%.10f", r.rawEnergy), 20), rpad(nistStr, 20), r.label)
    end
    println(stream, "  ", "-"^nx)

    return( nothing )
end

end # module

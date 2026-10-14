
"""
`module  JAC.ExportDataSet`
    ... a submodel of JAC that defines a small, format-agnostic container for one ion's spectroscopic data
        (levels, radiative A-values, effective electron-impact collision strengths), meant to sit between
        JAC's own computed outcomes and a specific external format writer such as `ExportAdas`. It carries
        NO physics and NO dependence on JAC's own `Level`/`Multiplet`/`Outcome` types: the caller builds the
        records from whatever computation produced them, picking the term label, (2S+1), L by hand for each
        level. This is deliberate -- a general, automatic derivation of an LS term label from an arbitrary
        (possibly jj-coupled, possibly heavy-element) JAC level is its own small research question, and
        bolting a best-effort version onto an export writer would silently mislabel exactly the levels where
        it matters, which is a worse failure than asking the caller for it explicitly.
"""
module ExportDataSet


"""
`struct  ExportDataSet.LevelRecord`
    ... one atomic level, in the plain form the external formats this module feeds want it: an LS-style
        term label (not derived here -- see the module's own docstring) and an excitation energy.

    + index          ::Int64     ... level index, 1-based, as the external format will reference it.
    + label          ::String    ... a short configuration/term label; free text, at the caller's choice.
    + multiplicity   ::Int64     ... 2S+1.
    + L              ::Int64     ... total orbital angular momentum quantum number.
    + J              ::Float64   ... total angular momentum quantum number.
    + energy         ::Float64   ... excitation energy relative to the ground level, in cm^-1.
"""
struct  LevelRecord
    index           ::Int64
    label           ::String
    multiplicity    ::Int64
    L               ::Int64
    J               ::Float64
    energy          ::Float64
end


"""
`struct  ExportDataSet.TransitionRecord`
    ... the radiative and electron-impact-excitation data for one pair of levels.

    + upperIndex     ::Int64              ... index of the upper level, referring to a LevelRecord.index.
    + lowerIndex     ::Int64              ... index of the lower level, referring to a LevelRecord.index.
    + aValue         ::Float64            ... Einstein A-value, in s^-1; 0.0 if not computed or forbidden.
    + upsilons       ::Array{Float64,1}   ... effective collision strengths, aligned one-to-one with the
                                               owning DataSet's own `temperatures` array.
"""
struct  TransitionRecord
    upperIndex      ::Int64
    lowerIndex      ::Int64
    aValue          ::Float64
    upsilons        ::Array{Float64,1}
end


"""
`struct  ExportDataSet.DataSet`
    ... the complete, format-agnostic record for one ion that an external-format writer consumes.

    + element               ::String                        ... element symbol, e.g. "H", "Nd".
    + ionCharge              ::Int64                         ... spectroscopic ion charge (0 for neutral).
    + nuclearCharge           ::Int64                        ... nuclear charge Z.
    + ionizationPotential      ::Float64                     ... ionisation potential of this ion stage, in cm^-1.
    + levels                   ::Array{LevelRecord,1}        ... the level list.
    + temperatures              ::Array{Float64,1}           ... electron temperatures, in Kelvin, that every
                                                                   TransitionRecord's `upsilons` are tabulated at.
    + transitions                ::Array{TransitionRecord,1} ... the transition list.
"""
struct  DataSet
    element             ::String
    ionCharge           ::Int64
    nuclearCharge       ::Int64
    ionizationPotential ::Float64
    levels              ::Array{LevelRecord,1}
    temperatures        ::Array{Float64,1}
    transitions         ::Array{TransitionRecord,1}
end

end # module

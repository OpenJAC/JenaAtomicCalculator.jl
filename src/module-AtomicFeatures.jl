
"""
`module  JAC.AtomicFeatures`  
... a submodel of JAC that contains all methods for computing the features for a neural-network training. 
    It assumes that the electronic and level structure of atoms and ions can be encoded by (large, real)
    x-vectors (feature vectors) that are based on a proper list of shells, subshells, orbitals as well as
    the configuration-interaction and coupling information of leading CSF.
"""
module AtomicFeatures


using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals, ..TableStrings

export  AtomicModel


"""
`struct  AtomicFeatures.AtomicModel`  
    ... defines a type for specifying the atomic model that is used in the training, test and use of an NN. It contains 
        all necessary data but not the underlying "recipi/scheme", how these data are generated and utilized in dealing 
        with the NN. At present, it is assumed that each NN is based on a set of shells (subshells) and associated
        orbitals, from which are features (descriptors) can be generated/extracted.
        

    + nMax              ::Int64                      ... maximum n-shell in the feature generation/extraction.
    + lMax              ::Int64                      ... maximum l-shell in the feature generation/extraction.
    + grid              ::Radial.Grid                ... The radial grid to be used for the computation.
    + asfSettings       ::ManyElectron.AsfSettings   ... AsfSettings used for simplified atomic-structure computations.
    + nuclearModel      ::Nuclear.Model              ... nuclear model used to define the orbitals.
    + shells            ::Array{Shell,1}             ... List of shells of model considered.
    + subshells         ::Array{Subshell,1}          ... List of subshells; this is redundant but simplifies the feature extraction. 
    + orbitals01        ::Dict{Subshell, Orbital}    ... Set of associated orbitals for charge state q=1.
    + orbitals02        ::Dict{Subshell, Orbital}    ... Set of associated orbitals for charge state q=2.
    + orbitals03        ::Dict{Subshell, Orbital}    ... Set of associated orbitals for charge state q=3.
"""
struct  AtomicModel
    nMax                ::Int64
    lMax                ::Int64  
    grid                ::Radial.Grid 
    asfSettings         ::ManyElectron.AsfSettings 
    nuclearModel        ::Nuclear.Model 
    shells              ::Array{Shell,1}  
    subshells           ::Array{Subshell,1} 
    orbitals01          ::Dict{Subshell, Orbital}
    orbitals02          ::Dict{Subshell, Orbital}
    orbitals03          ::Dict{Subshell, Orbital}
end


"""
`AtomicFeatures.AtomicModel()`  ... constructor for an 'empty' instance of AtomicFeatures.AtomicModel.
"""
function AtomicModel()
    AtomicModel( 0, 0, Radial.Grid(), ManyElectron.AsfSettings(), Nuclear.Model(1.), Shell[], Subshell[],
                 Dict{Subshell,Orbital}(), Dict{Subshell,Orbital}(), Dict{Subshell,Orbital}() )
end


# `Base.show(io::IO, am::AtomicFeatures.AtomicModel)`  ... prepares a proper printout of the variable am::AtomicFeatures.AtomicModel
function Base.show(io::IO, am::AtomicFeatures.AtomicModel) 
    println(io, "nMax:              $(am.nMax)  ")
    println(io, "lMax:              $(am.lMax)  ")
    println(io, "grid:              $(am.grid)  ")
    println(io, "asfSettings:       $(am.asfSettings)  ")
    println(io, "nuclearModel:      $(am.nuclearModel)  ")
    println(io, "shells:            $(am.shells)  ")
    println(io, "subshells:         $(am.subshells)  ")
    println(io, "orbitals01:        $(am.orbitals01)  ")
    println(io, "orbitals02:        $(am.orbitals02)  ")
    println(io, "orbitals03:        $(am.orbitals03)  ")
end


"""
`struct  AtomicFeatures.XyVector`
    ... defines a type for specifying a particular LSJ-coupled atomic level in terms of its LSJ-quantum
        numbers, a computed energy as well as a generate feature x-vector.

    + L                 ::AngularJ64          ... total L of the LSJ-level
    + S                 ::AngularJ64          ... total S of the LSJ-level
    + J                 ::AngularJ64          ... total J of the LSJ-level
    + parity            ::Basics.Parity       ... total parity J of the LSJ-level
    + energy            ::Float64             ... total energy, compute with the given atomic model.
    + nistEnergy        ::Float64             ... excitation energy, if available in NIST, and zero otherwise.
    + xVector           ::Array{Float64,1}    ... x-vector, generated with the given atomic model. 
"""
struct  XyVector
    L                   ::AngularJ64 
    S                   ::AngularJ64
    J                   ::AngularJ64 
    parity              ::Basics.Parity 
    energy              ::Float64  
    nistEnergy          ::Float64 
    xVector             ::Array{Float64,1} 
end 


# `Base.show(io::IO, v::AtomicFeatures.XyVector)`  ... prepares a proper printout of the variable v::AtomicFeatures.XyVector
function Base.show(io::IO, v::AtomicFeatures.XyVector) 
    println(io, "L:           $(v.L)  ")
    println(io, "S:           $(v.S)  ")
    println(io, "J:           $(v.J)  ")
    println(io, "parity:      $(v.parity)  ")
    println(io, "energy:      $(v.energy)  ")
    println(io, "nistEnergy:  $(v.nistEnergy)  ")
    println(io, "xVector:     $(v.xVector)  ")
end


#################################################################################################################################
#################################################################################################################################
## Feature-extraction primitives: each returns a fixed-length Vector{Float64} for a given, fixed shell/subshell list, so that
## every configuration/level built from the same AtomicModel produces an x-vector of the same total length.


"""
`AtomicFeatures.extractShellOccupations(shells::Array{Shell,1}, conf::Configuration)`
    ... returns, for each shell in the fixed `shells` list, the occupation number of that shell in `conf`
        (0.0 if the shell does not occur in conf). A Vector{Float64} of length(shells) is returned.
"""
function extractShellOccupations(shells::Array{Shell,1}, conf::Configuration)
    x = zeros(Float64, length(shells))
    for  (i, shell)  in  enumerate(shells)
        if  haskey(conf.shells, shell)   x[i] = Float64(conf.shells[shell])   end
    end
    return  x
end


"""
`AtomicFeatures.extractMeanOccupationNumbers(subshells::Array{Subshell,1}, level::Level)`
    ... returns, for each subshell in the fixed `subshells` list, its mean occupation number for the given
        level, obtained by averaging the occupation of each CSF in level.basis over the CI mixture (weight
        |level.mc[r]|^2). A Vector{Float64} of length(subshells) is returned.
"""
function extractMeanOccupationNumbers(subshells::Array{Subshell,1}, level::Level)
    x = zeros(Float64, length(subshells))
    basisIndex = Dict{Subshell,Int64}( s => k  for (k,s) in enumerate(level.basis.subshells) )
    for  r = 1:length(level.basis.csfs)
        w = abs(level.mc[r])^2;   csf = level.basis.csfs[r]
        for  (i, subshell)  in  enumerate(subshells)
            if  haskey(basisIndex, subshell)   x[i] = x[i] + w * csf.occupation[ basisIndex[subshell] ]   end
        end
    end
    return  x
end


"""
`AtomicFeatures.extractRkExpectation(subshells::Array{Subshell,1}, k::Int64, orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid)`
    ... returns <r^k> for the orbital of each subshell in the fixed `subshells` list, using
        RadialIntegrals.rkDiagonal(k, orb, orb, grid). Subshells not present in `orbitals` get 0.0.
        `grid` must be passed explicitly (the atomic model's grid) rather than relied upon via an
        orbital's own `.grid` field, which is not reliably populated for mean-field-generated
        orbitals -- the same convention every other JAC process module already follows.
        A Vector{Float64} of length(subshells) is returned.
"""
function extractRkExpectation(subshells::Array{Subshell,1}, k::Int64, orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid)
    x = zeros(Float64, length(subshells))
    for  (i, subshell)  in  enumerate(subshells)
        if  haskey(orbitals, subshell)
            orb   = orbitals[subshell]
            x[i]  = RadialIntegrals.rkDiagonal(k, orb, orb, grid)
        end
    end
    return  x
end


"""
`AtomicFeatures.extractFkIntegrals(subshells::Array{Subshell,1}, k::Int64, orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid)`
    ... returns the direct Slater integral F^k(a,b) = R^k(abab) for every unique pair (a,b) with a<=b (by
        position in the fixed `subshells` list), 0.0 if either orbital is missing from `orbitals`. The pair
        order is fixed by the `subshells` list, so every configuration built from the same list produces a
        (sub-)vector of the same length(subshells)*(length(subshells)+1)/2. `grid` must be passed explicitly
        (see AtomicFeatures.extractRkExpectation).
"""
function extractFkIntegrals(subshells::Array{Subshell,1}, k::Int64, orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid)
    n = length(subshells);   x = Float64[]
    for  i = 1:n
        for  j = i:n
            if  haskey(orbitals, subshells[i])  &&  haskey(orbitals, subshells[j])
                a = orbitals[subshells[i]];   b = orbitals[subshells[j]]
                push!(x, RadialIntegrals.SlaterRk(k, a, b, a, b, grid))
            else
                push!(x, 0.0)
            end
        end
    end
    return  x
end


"""
`AtomicFeatures.extractGkIntegrals(subshells::Array{Subshell,1}, k::Int64, orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid)`
    ... returns the exchange Slater integral G^k(a,b) = R^k(abba) for every unique pair (a,b) with a<=b (by
        position in the fixed `subshells` list), 0.0 if either orbital is missing from `orbitals`. Mirrors
        AtomicFeatures.extractFkIntegrals; same fixed pair order and (sub-)vector length; `grid` must be
        passed explicitly (see AtomicFeatures.extractRkExpectation).
"""
function extractGkIntegrals(subshells::Array{Subshell,1}, k::Int64, orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid)
    n = length(subshells);   x = Float64[]
    for  i = 1:n
        for  j = i:n
            if  haskey(orbitals, subshells[i])  &&  haskey(orbitals, subshells[j])
                a = orbitals[subshells[i]];   b = orbitals[subshells[j]]
                push!(x, RadialIntegrals.SlaterRk(k, a, b, b, a, grid))
            else
                push!(x, 0.0)
            end
        end
    end
    return  x
end


"""
`AtomicFeatures.extractIntermediateQN(nLeading::Int64, level::Level)`
    ... returns a fixed-length feature block describing the nLeading CSFs of `level` with the largest
        |mixing coefficient|^2 (idiom: sortperm(abs.(level.mc).^2, rev=true)), padded with zeros if level
        has fewer than nLeading CSFs. For each of the nLeading CSFs, four numbers are appended: [mixing
        coefficient, sum of subshell seniority numbers, the cumulative (total) subshellJ, the cumulative
        (total) subshellX]. This is a simplified, first-cut coupling summary (not a full per-subshell
        breakdown), matching the exploratory, script-like character requested for this feature. A
        Vector{Float64} of length 4*nLeading is returned.
"""
function extractIntermediateQN(nLeading::Int64, level::Level)
    x   = zeros(Float64, 4*nLeading)
    idx = sortperm(abs.(level.mc).^2, rev=true)
    for  p = 1:min(nLeading, length(idx))
        r         = idx[p];   csf = level.basis.csfs[r]
        seniority = Float64(sum(csf.seniorityNr))
        J         = Basics.twice(csf.subshellJ[end]) / 2.0
        X         = Basics.twice(csf.subshellX[end]) / 2.0
        x[4*(p-1)+1] = level.mc[r]
        x[4*(p-1)+2] = seniority
        x[4*(p-1)+3] = J
        x[4*(p-1)+4] = X
    end
    return  x
end

end # module


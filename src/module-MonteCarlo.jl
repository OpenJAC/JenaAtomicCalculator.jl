
# September 2025
#
# This code/module has been written without any tests so far; apart from syntax errors, a number of
# logical errors need to be expected. Nonetheless, it provides a rather complete set of functions in order 
# collect useful histograms of different kind and resolution from the stepwise cascade decay of atoms 
# and ions. Not all the functions have been realized in full details yet.

# At present, this module is not (yet) included into the JAC code.
# It might be useful to think about useful termination conditions for cases, in which the full tree
# need not to be analyzed for each individual shot. These conditions need to be formulated with proper care.


"""
`module  JAC.MonteCarlo`  
    ... a submodel of JAC that contains all methods for dealing with the (Monte-Carlo) simulations of coincidence
    spectra and maps. This module is based on (the definition of) a large set of MonteCarlos.Level's and the 
    idea of Monte-Carlo shots, which follow the probabilities to come along different charge states and atomic
    levels. The coincidence maps are then obtained by following a great numbers of shots in order to construct the 
    requested n-dimensional histograms. Each shot gives rise to one or several (particle) events, i.e. the emission
    of electrons, photons, etc. with certain energies (properties). In fact, this Monte-Carlo approach appears as
    the method of choice to extract useful one- and multi-dimensional histograms (abundance maps of coincidence
    events).
    
    This module provides a rather general machinery to extract different cuts (maps) through coincidence data.
    It enables one to extract such one- or multi-dimensional maps from cascade simultations with little adaption
    of input and specialized functions. Here, functions are provided to make single Monte-Carlo shots through a 
    tree (list) of MonteCarlo.Level's, by starting from a pre-set of such levels. Other functions enable the user 
    to initalize a histogram of given kind or to fill such histograms by analyzing the evens from single shots.
    
    A Monte-Carlo tree follows the idea of well-defined parent-daughter relations: Each parent level has none, 
    one or several daughters, i.e. subsequent levels to which it can decay. A Monte-Carlo shop ends, if no 
    daughter occurs, or ti proceeds to one of the daughters due to some given probabilities. While each level
    can act as parent, its parents themselve remain largely unknown and do not affect the subsequent steps.
    
    All levels are based  on the charge state, energy, parity, and total J of levels, a total decay rate as well 
    as a number of daughter levels, which are encoded by their level index, the emitted particle type and the 
    accumulated probability, i.e. the sum of all probabilities, including those for the current daughter in the list. 
    This encoding with accumulated probabilities simplifies the use of individual shots as well as the extraction 
    of useful events.
"""
module MonteCarlo


using Base.Threads, Distributed, Printf, ProgressMeter, SpecialFunctions,
        ..AngularMomentum, ..Defaults, ..TableStrings


"""
`abstract type MonteCarlo.AbstractEvent` 
    ... defines an abstract type to distinguish different types of events that can be extracted in or without
        coincidence. These individual events are collected from individual MonteCarlo shots through the 
        levels and later distributed to other (concidence) events, if needed. The following events concrete 
        data types are supported by this module:
    
    + struct MonteCarlo.ParticleEvent  
        ... to refer to the emission of an particle (electron, photon, ...) with certain energy.
    + struct MonteCarlo.Particle2Event  
        ... to refer to the coincidence emission of two particles (electrons, photons, ...) with certain energies.
"""
abstract type  AbstractEvents           end


"""
`struct  MonteCarlo.ParticleEvent <: MonteCarlo.AbstractEvent`  
    ... to refer to the emission of an particle (electron, photon, ...) with certain energy.

    + kind            ::Symbol            ... to indicate the emission of an electron (:e), photon (:p), ...
    + energy          ::Float64           ... List of shells from which excitations are to be considered.
"""
struct  ParticleEvent            <: MonteCarlo.AbstractEvent
    kind              ::Symbol 
    energy            ::Float64   
end 


"""
`MonteCarlo.ParticleEvent()`  ... constructor for an 'empty' instance of a MonteCarlo.ParticleEvent
"""
function  ParticleEvent()
    ParticleEvent(:x, 0.)
end


# `Base.show(io::IO, event::MonteCarlo.ParticleEvent)`  ... prepares a proper printout of event::MonteCarlo.ParticleEvent.
function Base.show(io::IO, event::MonteCarlo.ParticleEvent) 
    println(io, "kind:             $(event.kind)  ")
    println(io, "energy:           $(event.energy)  ")
end


"""
`struct  MonteCarlo.Particle2Event <: MonteCarlo.AbstractEvent`  
    ... to refer to the emission of two particles (electrons, photons, ...) with certain energies.
        The kind and energies are insignificant for the definition of the tuples, although the 
        sequence is the same in both tuples

    + kinds           ::Tuple{Symbol, Synbol}  .
        .. tuple of symbols of the two particles to indicate the emission of an electron (:e), photon (:p), ...
    + energies        ::Tuple{Float64, Float64} ... tuple of energies of the two particles.
"""
struct  ParticleEvent            <: MonteCarlo.AbstractEvent
    kinds             ::Tuple{Symbol, Synbol}
    energyies         ::Tuple{Float64, Float64}
end 


"""
`MonteCarlo.Particle2Event()`  ... constructor for an 'empty' instance of a MonteCarlo.Particle2Event
"""
function  Particle2Event()
    Particle2Event( (:x :x), (0.0, 0.0) )
end


# `Base.show(io::IO, event::MonteCarlo.Particle2Event)`  ... prepares a proper printout of event::MonteCarlo.Particle2Event.
function Base.show(io::IO, event::MonteCarlo.Particle2Event) 
    println(io, "kinds:            $(event.kinds)  ")
    println(io, "energies:         $(event.energies)  ")
end

"""
`abstract type MonteCarlo.AbstractHistogram` 
    ... defines an abstract type to distinguish and to collect data for different types of histograms with or without
        coincidence. The concrete histograms refer to 1-, 2- or 3-dimensional count maps with given xValuesm
        yValues and zValues (for the three different dimension). Each dimension x, y, z has a kind an a list
        of N+1 values to distinguish N intervals of different length. A histogram is filled by comparing the 
        values of different events with the kind of dimension and the values with the corresponding intervals
        of the dimension. Although we always add 1, the histograms are  stored in Array{Float64,n} arrays.
        The dimensions x, y, z refer (in this sequence) to different kinds of values, for which a Cartesian histrogram
        can be constructed. Concrete histograms are:
    
    + struct MonteCarlo.Histogram1D  
        ... to refer to a 1-dim histogram of particle events (electrons, photons, ...) with values in N intervals.
    + struct MonteCarlo.Histogram2D  
        ... to refer to a 2-dim histogram of particle events (electrons, photons, ...) with values in N intervals.
"""
abstract type  AbstractHistogram          end


"""
`struct  MonteCarlo.Histogram1D     <: MonteCarlo.AbstractHistogram`  
    ... to comprise the definition and (count) data of a 1-dimensional histogram for values
        xLower <= x < xUpper with stepsize xDelta = (xUpper - xLower) /xN.

    + xKind           ::Symbol            ... to indicate the emission of an electron (:e), photon (:p), ...
    + xN              ::Int64             ... Number of bins Nx.
    + xLower          ::Float64           ... Lower x-value of the histogram.
    + xUpper          ::Float64           ... Upper x-value of the histogram.
    + xValues         ::Array{Float64,1}  ... List of Nx + 1 values to define N intervals, into which the 
                                              events are counted.
    + counts          ::Array{Float64,1}  
        ... Number of counts in the Nx bins of the histogram; the i-th bin takes counts for 
            xValues[i] <= x < xValues[i+1], and the xValues[Nx+1] has no associated bin.
"""
struct  Histogram1D                 <: MonteCarlo.AbstractHistogram
    xKind             ::Symbol 
    xN                ::Int64   
    xLower            ::Float64 
    xUpper            ::Float64 
    xValues           ::Array{Float64,1}  
    counts            ::Array{Float64,1}  
end 


"""
`MonteCarlo.Histogram1D()`  ... constructor for an 'empty' instance of a MonteCarlo.Histogram1D
"""
function  Histogram1D()
    Histogram1D(:x, 0, 0., 0., Float64[], Float64[] )
end


# `Base.show(io::IO, histogram::MonteCarlo.Histogram1D)`  ... prepares a proper printout of histogram::MonteCarlo.Histogram1D.
function Base.show(io::IO, histogram::MonteCarlo.Histogram1Dt) 
    println(io, "xKind:         $(histogram.xKind)  ")
    println(io, "xN:            $(histogram.xN)  ")
    println(io, "xLower:        $(histogram.xLower)  ")
    println(io, "xUpper:        $(histogram.xUpper)  ")
    println(io, "xValues:       $(histogram.xValues)  ")
    println(io, "counts:        $(histogram.counts)  ")
end


"""
`struct  MonteCarlo.Histogram2D     <: MonteCarlo.AbstractHistogram`  
    ... to comprise the definition and (count) data of a 2-dimensional histogram for values
        xLower <= x < xUpper with stepsize xDelta = (xUpper - xLower) /xN   and
        yLower <= y < yUpper with stepsize yDelta = (yUpper - yLower) /yN.

    + xKind           ::Symbol            ... to indicate the emission of an electron (:e), photon (:p), ...
    + xN              ::Int64             ... Number of bins Nx.
    + xLower          ::Float64           ... Lower x-value of the histogram.
    + xUpper          ::Float64           ... Upper x-value of the histogram.
    + xValues         ::Array{Float64,1}  ... List of Nx + 1 values to define Nx intervals, into which the 
                                              events are counted.
    + yKind           ::Symbol            ... to indicate the emission of an electron (:e), photon (:p), ...
    + yN              ::Int64             ... Number of bins Ny.
    + yLower          ::Float64           ... Lower y-value of the histogram.
    + yUpper          ::Float64           ... Upper y-value of the histogram.
    + yValues         ::Array{Float64,1}  ... List of Ny + 1 values to define Ny intervals, into which the 
                                              events are counted.
    + counts          ::Array{Float64,2}  
        ... Number of counts in the Nx x Ny bins of the histogram; the (i,j)-th bin takes counts for 
            xValues[i] <= x < xValues[i+1] and yValues[j] <= y < yValues[j+1]; no bins are associate with
            xValues[Nx+1] and yValues[Ny+1], respectively.
"""
struct  Histogram2D                 <: MonteCarlo.AbstractHistogram
    xKind             ::Symbol 
    xN                ::Int64   
    xLower            ::Float64 
    xUpper            ::Float64 
    xValues           ::Array{Float64,1}  
    yKind             ::Symbol 
    yN                ::Int64   
    yLower            ::Float64 
    yUpper            ::Float64 
    yValues           ::Array{Float64,1}  
    counts            ::Array{Float64,2}  
end 


"""
`MonteCarlo.Histogram2D()`  ... constructor for an 'empty' instance of a MonteCarlo.Histogram2D
"""
function  Histogram1D()
    Histogram1D(:x, 0, 0., 0., Float64[], :y, 0, 0., 0., Float64[], zero(2,2) )
end


# `Base.show(io::IO, histogram::MonteCarlo.Histogram1D)`  ... prepares a proper printout of histogram::MonteCarlo.Histogram1D.
function Base.show(io::IO, histogram::MonteCarlo.Histogram1D) 
    println(io, "xKind:         $(histogram.xKind)  ")
    println(io, "xN:            $(histogram.xN)  ")
    println(io, "xLower:        $(histogram.xLower)  ")
    println(io, "xUpper:        $(histogram.xUpper)  ")
    println(io, "xValues:       $(histogram.xValues)  ")
    println(io, "yKind:         $(histogram.yKind)  ")
    println(io, "yN:            $(histogram.yN)  ")
    println(io, "yLower:        $(histogram.yLower)  ")
    println(io, "yUpper:        $(histogram.yUpper)  ")
    println(io, "yValues:       $(histogram.yValues)  ")
    println(io, "counts:        $(histogram.counts)  ")
end

#######################################################################################################################
#######################################################################################################################


"""
`struct  MonteCarlo.DaughterId`  
    ... defines a type to refer to and identify the daughter levels of a given level in the Monte-Carlo tree.

    + index          ::Int64       ... Index of the daughter level in the given Monte-Carlo tree.
    + kind           ::Symbol      
        ... symbol to indicate the emission process [electron (:e), photon (:p), ...] that lead to this daughter.
    + prob           ::Float64     
        ... Total probability of all daughters, including the present one; this total probability sums up
            the probabilities of all previous daughters and, hence, readily supports which daughter is reached
            in the next step. For the last daughter, this probability must always be 1.
"""
struct  DaughterId
    index            ::Int64 
    kind             ::Symbol      
    prob             ::Float64     
end 


"""
`MonteCarlo.DaughterId()`  ... constructor for an 'empty' instance of a MonteCarlo.DaughterId
"""
function  DaughterId()
    DaughterId( 0, :x, 0. )
end


# `Base.show(io::IO, daughter::MonteCarlo.DaughterId)`  ... prepares a proper printout of daughter::MonteCarlo.DaughterId.
function Base.show(io::IO, daughter::MonteCarlo.DaughterId) 
    println(io, "indexd:            $(daughter.index)  ")
    println(io, "kind:              $(daughter.kind)  ")
    println(io, "prob:              $(daughter.prob)  ")
end


"""
`struct  MonteCarlo.Level`  
    ... defines a (simple) type to order all cascade level into a structured Monte-Carlo tree; each level is uniquely 
        identified by the charge state, its total angular momentum, parity, energy and rate, and is associated with 
        a number of daughter-Id's to indicate, into which further levels the atom can decay to. The total angular 
        momentum, parity and rate of the level is often obsolete to follow individual shots through the cascade tree 
        but have been introduced to identify the physical origin of the level.
        
        An accumulated probability can be assigned to each level, if appropriate and if a list of such levels 
        are treated together. Such an accumulated probability, if a list of initial levels need to be specified, from
        which the individual shots may start from
        
        A cascade branch ends if not further daughter can be reached from a given level. As typical for a cascade, 
        there is no information about the parent level, from which the level was reached. In the future, it might be 
        useful to assign an origin::Symbol to each level, if the physical origin of levels need to be distinguished.
        
        The overall Monte-Carlo tree is tree::Vector{MonteCarlo.Level}. Care has to be taken that the tree is properly 
        linked internally before it can be used to follow individual shots through the cascade. The probabilities of the 
        daughters are always handled as partially-summed probabilities of the previous daughters, including the present
        daughter. This convention simplifies the use of such a cascade tree by keepig the size for its storage moderate.

    + totalJ         ::AngularJ64           ... Total angular momentum J of the level
    + parity         ::Basics.Parity        ... Total parity of the level.
    + charge         ::Int64                ... Charge state of the atom/ion.
    + energy         ::Float64              ... Total energy of the level; a shift of all levels has no particular meaning.
    + rate           ::Float64              ... Total decay rate of the level, if needed.
    + prob           ::Float64              ... accumulated probability, if needed, including the current level
    + daughters      ::Array{DaughterId,1}  ... Id's of all daughter levels, including their kind and partially-summed 
                                                probabilities.
"""
struct  Level
    totalJ           ::AngularJ64  
    parity           ::Basics.Parity  
    charge           ::Int64 
    energy           ::Float64 
    rate             ::Float64
    daughters        ::Array{DaughterId,1}  
end 


"""
`MonteCarlo.Level()`  ... constructor for an 'empty' instance of a MonteCarlo.Level
"""
function  Level()
    Level( AngularJ64(0), Basics.plus, 0, 0., 0., DaughterId[] )
end


# `Base.show(io::IO, level::MonteCarlo.Level)`  ... prepares a proper printout of level::MonteCarlo.Level.
function Base.show(io::IO, level::MonteCarlo.Level) 
    println(io, "totalJ:         $(level.totalJ)  ")
    println(io, "parity:         $(level.parity)  ")
    println(io, "charge:         $(level.charge)  ")
    println(io, "energy:         $(level.energy)  ")
    println(io, "rate:           $(level.rate)  ")
    println(io, "daughters:      $(level.daughters)  ")
end

#######################################################################################################################
#######################################################################################################################


"""
`MonteCarlo.addEvents!(counts::Array{Float64,1}, xValues::Array{Float64,1}, hEvents::Array{ParticleEvent,1})`  
    ... to add to counts the additional count(s) from the given events due to their energy, which needs to be 
        compared with xValues. The array counts
        is modified but nothing is returned otherwise.
"""
function addEvents!(counts::Array{Float64,1}, xValues::Array{Float64,1}, hEvents::Array{ParticleEvent,1})
    for  event  in  hEvents
        for  n  = 1:length(xValues)-1
            if  xValue[n] <= event.energy < xValue[n+1]   counts[n] = counts[n] + 1;   break   end
        end 
    end 
    
    return( nothing )
end


"""
`MonteCarlo.addEvents!(counts::Array{Float64,2}, xValues::Array{Float64,1}, yValues::Array{Float64,1}, 
                       hEvents::Array{Particle2Event,1})`  
    ... to add to counts the additional count(s) from the given events due to their energies, 
        which needs to be compared with xValues and yValues. The array counts is modified but nothing is 
        returned otherwise.
"""
function addEvents!(counts::Array{Float64,2}, xValues::Array{Float64,1}, yValues::Array{Float64,1}, 
                    hEvents::Array{Particle2Event,1})
    for  event  in  hEvents
        for  nx  = 1:length(xValues)-1,  ny  = 1:length(yValues)-1
            if  xValue[nx] <= event.energies[1] < xValue[nx+1]   && 
                yValue[ny] <= event.energies[2] < yValue[ny+1]   
                counts[nx, ny] = counts[nx, ny] + 1;   break   
            end
        end 
    end 
    
    return( nothing )
end


"""
`MonteCarlo.defineLevelTree(... line-lists, ...)`  
    ... to convert the data from a Cascade computation (or any other suitable data) into the format of a Monte-Carlo
        level tree, from which the desired histograms and coincidence maps can be readily extracted. A 
        tree::Array{MonteCarlo.Level,1} is returned for which all internal data and links to daughters are specified 
        properly. This tree is the basis for all data evaluations.
"""
function defineLevelTree()
    tree = MonteCarlo.Level[]
                
    return( tree )
end


"""
`MonteCarlo.determineRandomIndex(accumProbs::Array{Float64}, rn::Float64)`  
    ... returns the index in the vector of accumulated probabilities which is selected by the 
        random number rn. It is assumed that all random numbers are 0<= rn <= 1 and that the 
        last element of accumProbs[end] = 1. An ndx::Int64 is returned.
"""
function determineRandomIndex(accumProbs::Array{Float64}, rn::Float64)
    for  n = 1:length(accumProbs)   if   rn <= accumProbs[n]    return(n)   end    end
    error("stop a")
end


"""
`MonteCarlo.extractEvents(histrogram::Histogram1D, pEvents::Array{ParticleEvent,1})`  
    ... to extract from pEvent, a list of particle events from a shot through the level tree, those events, which can 
        contribute in their kind to the given histogram. A (reduced) list of hEvents::Array{ParticleEvent,1} is 
        returned. For 1D histograms this reduction is always simple.
"""
function extractEvents(histrogram::Histogram1D, pEvents::Array{ParticleEvent,1})
    hEvents = MonteCarlo.ParticleEvent[]
    for pEvent in pEvents   if pEvent.kind == histogram.xKind   push!(hEvents, pEvent)   end
                
    return( hEvents )
end


"""
`MonteCarlo.extractEvents(histrogram::Histogram2D, pEvents::Array{ParticleEvent,1})`  
    ... to extract from pEvent, a list of particle events from a shot through the level tree, those events, which can 
        contribute in their kinds to the given histogram. A (reduced) list of hEvents::Array{Particle2Event,1} is 
        returned. All pEvents are considered as coincident but need to be brought in the right sequence for the 
        histogram.
"""
function extractEvents(histrogram::Histogram2D, pEvents::Array{ParticleEvent,1})
    hEvents = MonteCarlo.Particle2Event[]
    for (p, pEvent)  in  enumerate(pEvents), (q, qEvent)  in  enumerate(pEvents)
        if p == q    continue   end
        if      (pEvent.kind, qEvent.kind) == histogram.kinds   
            push!(hEvents, Particle2Event( (pEvent.kind, qEvent.kind), (pEvent.energy, qEvent.energy) ) )
        elseif  (qEvent.kind, pEvent.kind) == histogram.kinds 
            push!(hEvents, Particle2Event( (qEvent.kind, pEvent.kind), (qEvent.energy, pEvent.energy) ) )
        end
    end
                
    return( hEvents )
end



"""
`MonteCarlo.initializeHistogram(histogram::MonteCarlo.Histogram1D)`  
    ... to initialize a 1D histogram from the given x-bounds and number of bins in the histogram.
        A newHistogram::MonteCarlo.Histogram1D is returned in which the xValues and counts arrays have 
        proper dimension to be filled due to the analysis of a sufficient large number of Monte-Carlo
        shots. A report about this initialization is provided as well.
"""
function initializeHistogram(histogram::MonteCarlo.Histogram1D)
    xDelta = (histogram.xUpper - histogram.xLower) / histogram.xN
    if  xDelta <= zero   error("Improper values for histogram initialization: xUpper=$(histogram.xUpper), " *
                               "xLower=$(histogram.xLower), xN=$(histogram.xN) ")
    end
    #
    xValues = Float64;   wa = 0.;   for  i = 1:histogram.xN   wa = wa + xDelta;    push!(xValues, wa)   end
    counts  = zeros(histogram.xN)
    
    newHistogram = MonteCarlo.Histogram1D(histogram.xKind, histogram.xN, histogram.xLower, histogram.xUpper,
                                          xValues, counts )
    
    if      histogram.xKind == :e   xKind = "electrons"
    elseif  histogram.xKind == :p   xKind = "photons"
    else    error("Unknown kind of histogram: xKind=$(histogram.xKind)")
    ends
    println(">> A 1D histogram has been initialized for $xKind with Nx=$(histogram.xN) bins and the bounds " *
            "xLower =$(histogram.xLower) : xUpper=$(histogram.xUpper)" *
            "\n   xValues = $(xValues[1:2]), ..., $(xValues[end-1:end]) ")
                
    return( newHistogram )
end


"""
`MonteCarlo.initializeHistogram(histogram::MonteCarlo.Histogram2D)`  
    ... to initialize a 2D histogram from the given x,y-bounds and numbers of x,y-bins in the histogram.
        A newHistogram::MonteCarlo.Histogram2D is returned in which the xValues/yValues and counts arrays have 
        proper dimension to be filled due to the analysis of a sufficient large number of Monte-Carlo
        shots. A report about this initialization is provided as well.
"""
function initializeHistogram(histogram::MonteCarlo.Histogram2D)
    xDelta = (histogram.xUpper - histogram.xLower) / histogram.xN
    yDelta = (histogram.yUpper - histogram.yLower) / histogram.yN
    if      xDelta <= zero   error("Improper values for histogram initialization: xUpper=$(histogram.xUpper), " *
                                   "xLower=$(histogram.xLower), xN=$(histogram.xN) ")
    elseif  xDelta <= zero   error("Improper values for histogram initialization: yUpper=$(histogram.yUpper), " *
                                   "yLower=$(histogram.yLower), yN=$(histogram.yN) ")
    end
    #
    xValues = Float64;   wa = 0.;   for  i = 1:histogram.xN   wa = wa + xDelta;    push!(xValues, wa)   end
    yValues = Float64;   wa = 0.;   for  i = 1:histogram.yN   wa = wa + yDelta;    push!(yValues, wa)   end
    counts  = zeros(histogram.xN, histogram.yN)
    
    newHistogram = MonteCarlo.Histogram1D(histogram.xKind, histogram.xN, histogram.xLower, histogram.xUpper, xValues, 
                                          histogram.yKind, histogram.yN, histogram.yLower, histogram.yUpper, yValues, 
                                          counts )
    
    if      histogram.xKind == :e   xKind = "electrons"
    elseif  histogram.xKind == :p   xKind = "photons"
    else    error("Unknown kind of histogram: xKind=$(histogram.xKind)")
    ends
    if      histogram.yKind == :e   yKind = "electrons"
    elseif  histogram.yKind == :p   yKind = "photons"
    else    error("Unknown kind of histogram: yKind=$(histogram.yKind)")
    ends
    
    println(">> A 2D histogram has been initialized for ($xKind, $yKind) with " *
            "\n   Nx=$(histogram.xN) bins and bounds " *
            "xLower =$(histogram.xLower) : xUpper=$(histogram.xUpper), " *
            "xValues = $(xValues[1:2]), ..., $(xValues[end-1:end]) "     *
            "\n   Ny=$(histogram.yN) bins and bounds " *
            "yLower =$(histogram.yLower) : yUpper=$(histogram.yUpper), " *
            "yValues = $(yValues[1:2]), ..., $(yValues[end-1:end]) " )
                
    return( newHistogram )
end


"""
`MonteCarlo.initializeStartLevels(tree::Array{MonteCarlo.Level,1})`  
    ... to initialize the set of starting levels from which the individual Monte-Carlo shots may start from.
        Like for daughter levels, an accumulated probability (including the current level) is used to select
        a starting level. Moreover, all indices refer to the given Monte-Carlo tree. A list of
        startLevels::Array{MonteCarlo.Level,1} is returned, which are typical a subset of the given Monte-Carlo 
        tree.
"""
function initializeStartLevels(tree::Array{MonteCarlo.Level,1})
    startLevels = MonteCarlo.Level[]
    # First collect proper start levels/indices and determine their (individual) probabilities.
                
    # Complete the initialization by computing the accumulated probabilities for all startLevels.
                
    return( startLevels )
end


"""
`MonteCarlo.runMultipleShots(tree::Array{MonteCarlo.Level,1}, startLevels::Array{MonteCarlo.Level,1}, nShots::Int64,
                             histogram::Histogram1D)`  
    ... to run nShots (individual) single shot through the given tree by starting from a randomly selected
        startLevel. The particle-events from each shot are analyzed how they might contributed to the requested
        1D histogram  and the counts are added accordingly.  A newHistogram::Histogram1D is returned in which 
        all useful counts from the nShots are added properly.
"""
function runMultipleShots(tree::Array{MonteCarlo.Level,1}, startLevels::Array{MonteCarlo.Level,1}, nShots::Int64,
                          histogram::Histogram1D)
    # Extract the starting probabilities to later find the right startLevel
    startProbs = Float64[];    for  level in startLevels   push!(startProbs, level.prob)   end 
    counts     = zeros(histogram.xN)
    
    for ns = 1:nShots
        # Run nShots by first randomly selecting a startLevel
        ndx     = MonteCarlo.determineRandomIndex(startProbs, rand() )
        pEvents = MonteCarlo.runSingleShot(tree, startLevels[ndx])
        
        # Analyze the events for their use in the histogram; only events of proper particle type are useful
        hEvents = MonteCarlo.extractEvents(histrogram, pEvents)
        MonteCarlo.extractEvents!(counts, hEvents)
    end
    
    newHistogram = MonteCarlo.Histogram1D(histogram.xKind, histogram.xN, histogram.xLower, histogram.xUpper,
                                          histogram.xValues, counts )
    
    if      histogram.xKind == :e   xKind = "electrons"
    elseif  histogram.xKind == :p   xKind = "photons"
    else    error("Unknown kind of histogram: xKind=$(histogram.xKind)")
    ends
    println(">> A 1D histogram has been updated after $nShots shots for $xKind with Nx=$(histogram.xN) bins and the bounds " *
            "xLower =$(histogram.xLower) : xUpper=$(histogram.xUpper)" *
            "\n   xValues = $(histogram.xValues[1:2]), ..., $(histogram.xValues[end-1:end]) ")
                
    return( newHistogram )
end


"""
`MonteCarlo.runMultipleShots(tree::Array{MonteCarlo.Level,1}, startLevels::Array{MonteCarlo.Level,1}, nShots::Int64,
                             histogram::Histogram2D)`  
    ... to run nShots (individual) single shot through the given tree by starting from a randomly selected
        startLevel. The particle-events from each shot are analyzed how they might contributed to the requested
        2D histogram  and the counts are added accordingly.  A newHistogram::Histogram2D is returned in which 
        all useful counts from the nShots are added properly.
"""
function runMultipleShots(tree::Array{MonteCarlo.Level,1}, startLevels::Array{MonteCarlo.Level,1}, nShots::Int64,
                          histogram::Histogram2D)
    # Extract the starting probabilities to later find the right startLevel
    startProbs = Float64[];    for  level in startLevels   push!(startProbs, level.prob)   end 
    counts     = zeros(histogram.xN, histogram.yN)
    
    for ns = 1:nShots
        # Run nShots by first randomly selecting a startLevel
        ndx     = MonteCarlo.determineRandomIndex(startProbs, rand() )
        pEvents = MonteCarlo.runSingleShot(tree, startLevels[ndx])
        
        # Analyze the events for their use in the histogram; only events of proper particle type are useful
        hEvents = MonteCarlo.extractEvents(histrogram, pEvents)
        MonteCarlo.addEvents!(counts, hEvents)
    end
    
    newHistogram = MonteCarlo.Histogram2D(histogram.xKind, histogram.xN, histogram.xLower, histogram.xUpper,
                                          histogram.xValues, 
                                          histogram.yKind, histogram.yN, histogram.yLower, histogram.yUpper,
                                          histogram.yValues, counts )
    
    if      histogram.xKind == :e   xKind = "electrons"
    elseif  histogram.xKind == :p   xKind = "photons"
    else    error("Unknown kind of histogram: xKind=$(histogram.xKind)")
    ends
    if      histogram.yKind == :e   yKind = "electrons"
    elseif  histogram.yKind == :p   yKind = "photons"
    else    error("Unknown kind of histogram: yKind=$(histogram.yKind)")
    ends
    println(">> A 2D histogram has been updated after $nShots shots for ($xKind, $yKind) with " *
            "\n   Nx=$(histogram.xN) bins and bounds " *
            "xLower =$(histogram.xLower) : xUpper=$(histogram.xUpper), " *
            "xValues = $(xValues[1:2]), ..., $(xValues[end-1:end]) "     *
            "\n   Ny=$(histogram.yN) bins and bounds " *
            "yLower =$(histogram.yLower) : yUpper=$(histogram.yUpper), " *
            "yValues = $(yValues[1:2]), ..., $(yValues[end-1:end]) " )
                
    return( newHistogram )
end


"""
`MonteCarlo.runSingleShot(tree::Array{MonteCarlo.Level,1}, startLevel::MonteCarlo.Level)`  
    ... to run a single shot through the given tree. A set of events::Array{MonteCarlo.ParticleEvent,1} is returned, 
        which summarize the emitted electrons/photons from the given shot. These particle events need to be 
        further analyzed to extract the desired (coincidence) events from the cascade.
"""
function runSingleShot(tree::Array{MonteCarlo.Level,1}, , startLevel::MonteCarlo.Level)
    events = MonteCarlo.ParticleEvent[]
    # Collect the list of events by following the parent-daughter relations in the tree until the tree ends
    # or some termination condition applies. Each individual step usually contributes one particle event
    # to the events.
    # First, randomly select a startLevel and (randomly) follow through the tree until it terminates.
                
    return( events )
end


#######################################################################################################################
#######################################################################################################################
# This function must be called within the cascade module to perform the corresponding simulations by
# using the associated Cascade.Scheme.


"""
`MonteCarlo.runSimulation( atomic line data, histogram::MonteCarlo.AbstractHistogram, 
                           scheme::Cascade.MonteCarloScheme)`  
    ... perform the Monte-Carlo simulation of a one- or multi-dimensional histogram of different kind and 
        purpose. A histogram::MonteCarlo.Histogram1D or histogram::MonteCarlo.Histogram2D is returned, in dependence 
        of the given input from scheme::Cascade.MonteCarloScheme.
        
        The following tasks need to be solved by this function
        (1) tree        = MonteCarlo.defineLevelTree(... line-lists, ...)
        (2) startLevels = MonteCarlo.initializeStartLevels(tree::Array{MonteCarlo.Level,1})
        (3) histogram   = MonteCarlo.runMultipleShots(tree, startLevels, nShots::Int64, histogram)
        (4) return final histogram and more.
"""
function runSimulation()
    histogram = MonteCarlo.Histogram1D()
                
    return( histogram )
end

end # module




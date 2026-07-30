
"""
`module  JAC.DeepLearning`  
... a submodel of JAC that contains all methods for dealing with the training, test and use of neutral
    networks for atomic applications. This module distinguishes between different NN and associated 
    applications, whose details are compiled in the corresponding inc files.
    
    Each NN on atomic data/information is typically based on a particular AtomicModel which comprises
    a set of shells (subshells) and orbitals as well as a "recipi" how the atomic features are constructed/
    extracted from this information. This may include simplified atomic-structure computations.
    To extract the information from a (trained) NN, the particular requests are formulated & processes
    by an DeepLearning.Application(), which can be "preform"(ed) like other Computation's and Estimate's
    within the JAC toolbox.
"""
module DeepLearning


using Printf, ..AngularMomentum, ..AtomicFeatures, ..AtomicState, ..Basics, ..Defaults, ..Hamiltonian, ..LSjj, ..ManyElectron,
      ..Nuclear, ..Radial, ..TableStrings

export  AbstractNeuralNetwork, AbstractNeuralNetworkRequest, LevelEstimationRequest, 
        Application


"""
`abstract type DeepLearning.AbstractNeuralNetwork` 
    ... defines an abstract type to distinguish different neutral networks that were trained on atomic data; see als
    
    + struct LevelEstimationNNn4  
        ... to characterize the neutral network that is based on a shell list with nMax = 4.
"""
abstract type  AbstractNeuralNetwork                                   end
struct   LevelEstimationNNn4  <:  DeepLearning.AbstractNeuralNetwork   end


"""
`abstract type DeepLearning.AbstractNeuralNetworkRequest` 
    ... defines an abstract type to distinguish different requests on neutral networks that were trained on atomic 
        data. These requests distinguish between different applications of such network; see als
    
    + struct LevelEstimationRequest  
        ... to characterize the configuration and details, for which level estimates need to be made (and
            compare with available NIST data).
"""
abstract type  AbstractNeuralNetworkRequest       end



"""
`struct  DeepLearning.Application`  
    ... defines a type for a cascade computation, i.e. for the computation of a whole photon excitation, photon ionization and/or 
        decay cascade. The -- input and control -- data from this computation can be modified, adapted and refined to the practical needs, 
        and before the actual computations are carried out explictly. Initially, this struct just contains the physical meta-data about the 
        cascade, that is to be calculated, but a new instance of the same Cascade.Computation gets later enlarged in course of the 
        computation in order to keep also wave functions, level multiplets, etc.

    + name               ::String                          ... A name for the cascade
    + atomicModel        ::AtomicFeatures.AtomicModel      ... Atomic model used for feature generation & extraction.
    + request            ::AbstractNeuralNetworkRequest    ... User request for extracting data from the NN.
"""
struct  Application
    name                 ::String
    atomicModel          ::AtomicFeatures.AtomicModel
    request              ::AbstractNeuralNetworkRequest 
end 


"""
`DeepLearning.Application()`  ... constructor for an 'default' instance of a DeepLearning.Application.
"""
function Application()
    Application( "Default deep-learning application", AtomicFeatures.AtomicModel(), LevelEstimationRequest() )
end


"""
`DeepLearning.Application(app::DeepLearning.Application;`
    
            name=..,               atomicModel=..,          request=..)
            
    ... constructor for re-defining the application::DeepLearning.Application.
"""
function Application(app::DeepLearning.Application;                              
    name::Union{Nothing,String}=nothing,                                                         
    atomicModel::Union{Nothing,AtomicFeatures.AtomicModel}=nothing,   request::Union{Nothing,AbstractNeuralNetworkRequest}=nothing)
    
    if  isnothing(name)                   namex                 = app.name                     else  namex = name                              end
    if  isnothing(atomicModel)            atomicModelx          = app.atomicModel               else  atomicModelx = atomicModel                end
    if  isnothing(request)                requestx              = app.request                  else  requestx = request                        end
    
    Application(namex, atomicModelx, requestx)
end


# `Base.string(applic::DeepLearning.Application)`  ... provides a String notation for the variable applic::DeepLearning.Application.
function Base.string(applic::DeepLearning.Application)
    if     typeof(applic.request) == DeepLearning.LevelEstimationRequest     sb = "request for level estimation"
    else   error("unknown typeof(applic.request)")
    end

    sa = "Deep learning application  $(applic.name)  with a $sb and the atomic model $(applic.atomicModel)  as well as "
    sa = sa * "for Z = $(applic.atomicModel.nuclearModel.Z) and for the configurations: \n "

    return( sa )
end


# `Base.show(io::IO, applic::DeepLearning.Application)`  ... prepares a proper printout applic::DeepLearning.Application.
function Base.show(io::IO, applic::DeepLearning.Application)
    sa = Base.string(applic)
    print(io, sa, "\n")
end


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

"""
`Basics.run(applic::DeepLearning.Application)`  
    ... to set-up and run a deep-learning application that is based on a particular atomic model and some given request.
        The results of all individual steps are printed to screen but nothing is returned otherwise.

`Basics.run(applic::DeepLearning.Application; output::Bool=true)`   
    ... to run the same but to return the complete output in a dictionary;  the particular output depends on the type 
        and specifications of the deep-learning application but can easily accessed by the keys of this dictionary.
"""
function Basics.run(applic::DeepLearning.Application; output::Bool=true)
    DeepLearning.run(applic.request, applic::DeepLearning.Application, output=output)
end


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

include("module-DeepLearning-inc-LevelEstimates.jl")

end # module


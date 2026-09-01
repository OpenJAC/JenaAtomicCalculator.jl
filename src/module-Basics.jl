
"""
`module  JAC.Basics`  
	... a submodel of JAC that contains many basic types/struct that are specific to the JAC module; this module also defines
	    a number of basic functions/methods that are later extended. We here provide proper docstrings for all abstract and
	    conrete types (struct) in order to allow the user to easy understand and access the individual fields of each type/struct 
	    definition.
"""
module Basics


using Printf, QuadGK

include("module-Basics-inc-first.jl")
include("module-Basics-inc-abstract.jl")
include("module-Basics-inc-second.jl")

export  add, analyze, compute, diagonalize, generate, integrate, perform, recast, run,
        tabulate       

export  checkConfigurations,  displayConfiguration,  displayConfigurations,  extractConfiguration,  extractConfigurations,  
        extractFromConfiguration,  extractFromConfigurations,  generateConfigurations
        
# Functions/methods that are later added to the module Basics
function add                                                    end
function addZerosToCsfR                                         end
function analyze                                                end
function analyzeConvergence                                     end
function analyzeGrid                                            end
function checkConfigurations                                    end
function compute                                                end
function computeDensity                                         end
function computeDiracEnergy                                     end
function computeMeanSubshellOccupation                          end
function computeMultipletForGreenApproach                       end
function computePotential                                       end
function determineEnergySharings                                end
function determineHoleShells                                    end
function determineMeanEnergy                                    end
function determineNearestPoints                                 end
function determineNonorthogonalShellOverlap                     end
function determinePolarizationLambda                            end
function determinePolarizationVector                            end
function determineSelectedLines                                 end
function determineSelectedPathways                              end
function diagonalize                                            end
function diracDelta                                             end
function display                                                end
function displayConfigurationThemes                             end
function displayConfiguration                                   end
function displayConfigurations                                  end
function displayLevels                                          end
function displayMeanEnergies                                    end
function displayOrbitalOverlap                                  end
function displayOrbitalProperties                               end
function expandOrbital                                          end
function extractConfiguration                                   end
function extractConfigurations                                  end
function extractFromConfiguration                               end
function extractFromConfigurations                              end
function extractMeanEnergy                                      end
function extractMeanOccupation                                  end
function extractNonrelativisticShellList                        end
function extractOpenShellQNfromCsfR                             end
function extractRelativisticConfigurationFromCsfR               end
function extractRelativisticSubshellList                        end
function extractRydbergSubshellList                             end
function extractShellList                                       end
function extractShellOccupationFromCsfR                         end
function extractSubshellList                                    end
function extractValenceShell                                    end
function FermiDirac                                             end
function fixEigenvectorPhase!                                   end
function generate                                               end
function generateBasis                                          end
function generateConfigurations                                 end
function generateConfigurationsForExcitationScheme              end
function generateCsfRs                                          end
function generateFieldCoordinates                               end
function generateLevelWithExtraElectron                         end
function generateLevelWithExtraSubshell                         end
function generateLevelWithExtraSubshells                        end
function generateLevelWithExtraTwoElectrons                     end
function generateLevelWithSymmetryReducedBasis                  end
function generateMeshCoordinates                                end
function generateOrbitalsForBasis                               end
function printRasStepDiagnostic                                 end
function generateOrbitalSuperposition                           end
function generateShellList                                      end
function generateSubshellList                                   end
function generateSpectrumLorentzian                             end
function hasSubshell                                            end
"""
`Basics.integrate()`
    ... generic function for the numerical integration of a (radial) function on a grid; see the methods of
        Basics.integrate for the integration rules that are available.
"""
function integrate                                              end
function interpolateOnGridGrasp92                               end
function isSimilar                                              end
function isSymmetric                                            end
function isStandardSubshellList                                 end
function isViolated                                             end
function isZero                                                 end
function lastPoint                                              end
function merge                                                  end
function modifyLevelMixing                                      end
function perform                                                end
function performCI                                              end
function plot                                                   end
function read                                                   end
function readCslFileGrasp92                                     end
function readFilesGrasp18                                  		end
function readOrbitalFileGrasp92                                 end
function readMixFileRelci                                       end
function readMixingFileGrasp18                                  end
function recast                                                 end
function recommendedGrid                                        end
function run                                                    end
function selectLevel                                            end
function selectLevelPair                                        end
function selectLevelTriple                                      end
function selectSymmetry                                         end
function shiftTotalEnergies                                     end
function slaterScreening                                        end
function sortByEnergy                                           end
function subshellStateString                                    end
function tabulate                                               end
function tabulateKappaSymmetryEnergiesDirac                     end

end  ## module

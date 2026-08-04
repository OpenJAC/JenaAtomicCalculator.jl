#
println("Am) Apply & test procedures to establish a configuration-based language and to deal with electron configurations.")

# This file tests many different features by first comparing the new procedures in module-BasicsAZ-inc-configurations.jl 
# to previously developed functions; it also test new features by using the special 
# configuration-themes  <:  AbstractConfigurationTheme. This file does not strictly follow the list of themes as displayed
# by ? Basics.displayConfigurationthemes

include("test-code.jl")

heConf   = Configuration("[He]");    arConf = Configuration("[Ar]");    xeConf =  Configuration("[Xe]")
heShells = Basics.generateShellList(2, 4, [0,1,2,3]);          heSubshells = Basics.generateSubshellList(heShells)
arShells = Basics.generateShellList(4, 6, [0,1,2,3,4,5]);      arSubshells = Basics.generateSubshellList(arShells)
xeShells = Basics.generateShellList(6, 8, [0,1,2,3,4,5,6]);    xeSubshells = Basics.generateSubshellList(xeShells)
allConfs = Basics.merge([heConf], [arConf], [xeConf])

if      false
    # Basic operations
    # ----------------
    # Merge three configuration lists
    wa = Basics.mergeConfigurations([heConf], [arConf], [xeConf])
         Basics.displayConfigurations(stdout, wa, longForm=true, details="He + Ar + Xe")
    #
elseif  false
    # Add, excite or remove electrons from some reference configurations
    # ------------------------------------------------------------------
    fromShells = Basics.generateShellList(5, 5, 2)
    intoShells = Basics.generateShellList(6, 7, 5)
    wa = Basics.generateConfigurations(AddElectrons(1, intoShells), [arConf])
         Basics.displayConfigurations(stdout, wa, longForm=true, details="1 added electron")
    wb = Basics.generateConfigurations(ExciteElectrons(1, fromShells, intoShells), [xeConf])
         Basics.displayConfigurations(stdout, wa, longForm=true, details="1 excited electron")
    wc = Basics.generateConfigurations(RemoveElectrons(1, fromShells), [xeConf])
         Basics.displayConfigurations(stdout, wa, longForm=true, details="1 removed electron")
    #
elseif  false
    # Add, excite or remove electrons due to given atomic processes
    # -------------------------------------------------------------
    fromShells = [Shell("3s"), Shell("3p")]
    toShells   = [Shell("3d")]
    intoShells = Basics.generateShellList(4, 5, 3)
    conf       = Configuration("1s 2s^2 2p^6 3s^2 3p^4 3d^5")
    wa = Basics.generateConfigurations(ForAutoIonization(), [conf])
         Basics.displayConfigurations(stdout, wa, longForm=true, details="ForAutoIonization")
    wb = Basics.generateConfigurations(ForDielectronicCapture(fromShells, toShells, intoShells), [arConf])
         Basics.displayConfigurations(stdout, wb, longForm=true, details="ForDielectronicCapture")
    wc = Basics.generateConfigurations(ForElectronCapture(intoShells), [arConf])
         Basics.displayConfigurations(stdout, wc, longForm=true, details="ForElectronCapture")
    wd = Basics.generateConfigurations(ForPhotoEmission(), [conf])
         Basics.displayConfigurations(stdout, wd, longForm=true, details="ForPhotoEmission")
    we = Basics.generateConfigurations(ForPhotoIonization(), [arConf])
         Basics.displayConfigurations(stdout, we, longForm=true, details="ForPhotoIonization")
    wf = Basics.generateConfigurations(ForPhotoRecombination(intoShells), [arConf])
         Basics.displayConfigurations(stdout, wf, longForm=true, details="ForPhotoRecombination")
    wg = Basics.generateConfigurations(ForRasExcitations(true, true, false, false, fromShells, intoShells), [arConf])
         Basics.displayConfigurations(stdout, wg, longForm=true, details="ForRasExcitations")
    #
elseif  true
    # Add, excite or remove electrons due to given atomic processes
    # -------------------------------------------------------------
    conf        = Configuration("1s^2 2s^2 2p^2")
    fromShells  = [Shell("2s"), Shell("2p")]
    toShells    = [Shell("2p"), Shell("3s")]
    intoShells  = Basics.generateShellList(4, 4, 3)
    decayShells = Basics.generateShellList(2, 3, 3)
    wa = Basics.generateConfigurations(ForDielectronicRecombination(fromShells, toShells, intoShells, decayShells), [conf])
         Basics.displayConfigurations(stdout, wa[1], longForm=true, details="ForDielectronicRecombination: intermediate confs")
         Basics.displayConfigurations(stdout, wa[2], longForm=true, details="ForDielectronicRecombination: final confs")
    #
elseif  false
    # Add, excite or remove electrons due to given atomic processes
    # -------------------------------------------------------------
    fromShells  = [Shell("3s"), Shell("3p")]
    toShells    = [Shell("3d")]
    intoShells  = Basics.generateShellList(5, 5, 3)
    decayShells = Basics.generateShellList(3, 4, 3)
    conf       = Configuration("1s 2s^2 2p^6")
    wa = Basics.generateConfigurations(ForStepwiseDecay(2), [conf])
         Basics.displayConfigurations(stdout, wa, longForm=true, details="ForStepwiseDecay")
    wb = Basics.generateConfigurations(ForHollowIons(2, intoShells, decayShells), [conf])
         Basics.displayConfigurations(stdout, wb, longForm=true, details="ForHollowIons")
         Basics.displayConfigurations(stdout, ByNumber([7, 8, 9, 10, 11]), wb, longForm=false, details="ForHollowIons")
    #
elseif  false
    #
    # Generate and extract configurations by applying additional restrictions
    # -----------------------------------------------------------------------
    fromShells = Basics.generateShellList(5, 5, 2)
    intoShells = Basics.generateShellList(6, 7, 5)
    restrictions = AbstractConfigurationRestriction[ RestrictParity(Basics.plus), RestrictToShellDoubles(7,2) ]
    ## restrictions = AbstractConfigurationRestriction[ RestrictParity(Basics.plus) ]
    thema        = ExciteElectrons(1, fromShells, intoShells)
    themeA       = Basics.RestrictExcitations(restrictions)
    themeB       = Basics.RestrictExcitations(thema, restrictions)
    wa           = Basics.generateConfigurations(ExciteElectrons(1, fromShells, intoShells), [xeConf])
                   Basics.displayConfigurations(stdout, wa, details="ExciteElectrons")
    wb           = Basics.extractConfigurations(themeA, wa)
                   Basics.displayConfigurations(stdout, wb, details="... extract by RestrictExcitations")
    wc           = Basics.generateConfigurations(themeB, [xeConf])
                   Basics.displayConfigurations(stdout, wc, details="... generate by RestrictExcitations")
elseif  false
    # Generate relativistic configurations
    # ------------------------------------
    conf = Configuration("[Ar] 4s 4p")
    wa = Basics.generateConfigurations(RelativisticConfigurations(), arConf);                           @show wa
    wb = Basics.generateConfigurations(RelativisticConfigurations(), conf);                             @show wb
    wc = Basics.generateConfigurations(RelativisticConfigurations(), xeConf);                           @show wc
    #
elseif  false
    # Extract some particular single configuration
    # -------------------------------------------- 
    name        = "Oxygen 1s^2 2s^2 2p^4 ground configuration"
    refConfigs  = [Configuration("[He] 2s^2 2p^4"), Configuration("[He] 2s 2p^5"), ]
    mfSettings  = MeanFieldSettings()
    #
    va          = Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, MeanFieldBasis(mfSettings) )
    vb          = generate(va, output=true)
    basis       = vb["mean-field basis"]

    wa = Basics.extractConfiguration(GroundConfiguration(29.,24));                                      @show wa
    wb = Basics.extractConfiguration(GroundConfiguration(29.,29));                                      @show wb
    wc = Basics.extractConfiguration(GroundConfiguration(54.,32));                                      @show wc
    wd = Basics.extractConfiguration(FromBasis(), basis, basis.csfs[1]);                                @show wd
    we = Basics.extractConfiguration(FromBasis(), basis, basis.csfs[8]);                                @show we
    #
elseif  false
    # Extract information from a single (relativistic) configuration
    # --------------------------------------------------------------
    wa = Atomic.Computation(Atomic.Computation(), name="Oxygen", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(18.), 
                            configs=[Configuration("[He] 2s^2 2p^4"), Configuration("[He] 2s 2p^5")])
    wb = perform(wa, output=true)
    mp = wb["multiplet:"]
    #
    wc = Basics.extractConfiguration(LeadingConfiguration(), mp.levels[1]);                             @show wc
    wd = Basics.extractConfiguration(LeadingConfiguration(), mp.levels[9]);                             @show wd
    we = Basics.extractConfiguration(LeadingConfigurationR(), mp.levels[1]);                            @show we
    wf = Basics.extractConfiguration(LeadingConfigurationR(), mp.levels[9]);                            @show wf
    #
elseif  false 
    # Extract selected configurations from large configuration lists
    # --------------------------------------------------------------
         Basics.displayConfigurations(stdout, allConfs, longForm=true, details="allConfs")
    wa = Basics.extractConfigurations(ByNumber([18]), allConfs)
         Basics.displayConfigurations(stdout, wa, longForm=true, details="ByNumber([18])")
    wb = Basics.extractConfigurations(ByNumber([2, 18]), allConfs)
         Basics.displayConfigurations(stdout, wb, longForm=true, details="ByNumber([2, 18])")
    #
elseif  false
    # Extract selected configurations from given lists or information about them
    # -------------------------------------------------------------------------- 
    wa = Basics.extractConfigurations(FromBasis(), basis);                                              @show wa
    wb = Basics.extractConfigurations(ByNumber([18]), allConfs);                                        @show wb
    wc = Basics.extractConfigurations(ByParity(Basics.plus), allConfs);                                 @show wc
    wd = Basics.extractConfigurations(ByParity(Basics.minus), allConfs);                                @show wd
    wc = Basics.extractConfigurations(TotalAM(false, [AngularJ64(0)]), allConfs);                       @show wc
    wd = Basics.extractConfigurations(TotalAM(false, [AngularJ64(3//2)]), allConfs);                    @show wd
    #
elseif  false
    # Extract several (relativistic) configurations
    # --------------------------------------------- 
    name        = "Oxygen 1s^2 2s^2 2p^4 ground configuration"
    refConfigs  = [Configuration("[He] 2s^2 2p^4"), Configuration("[He] 2s 2p^5"), ]
    mfSettings  = MeanFieldSettings()
    #
    va          = Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, MeanFieldBasis(mfSettings) )
    vb          = generate(va, output=true)
    basis       = vb["mean-field basis"]
    #
    wa = Basics.extractConfigurations(RelativisticConfigurations(), basis);                             @show wa
    wb = Basics.extractConfigurations(RelativisticConfigurations(), basis, AngularJ64(0));              @show wb
    wc = Basics.extractConfigurations(RelativisticConfigurations(), basis, AngularJ64(1));              @show wc
    wd = Basics.extractConfigurations(TotalAM(false, [AngularJ64(0)]), refConfigs);                     @show wd
    we = Basics.extractConfigurations(TotalAM(false, [AngularJ64(2)]), refConfigs);                     @show we
    #
elseif  false
    # Extract information from a single (relativistic) configuration
    # --------------------------------------------------------------
    name        = "Oxygen 1s^2 2s^2 2p^4 ground configuration"
    refConfigs  = [Configuration("[He] 2s^2 2p^4"), Configuration("[He] 2s 2p^5"), ]
    mfSettings  = MeanFieldSettings()
    #
    va          = Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, MeanFieldBasis(mfSettings) )
    vb          = generate(va, output=true)
    basis       = vb["mean-field basis"]
    #
    wa = Basics.extractConfigurations(RelativisticConfigurations(), basis);                             @show wa
    wb = Basics.extractFromConfiguration(ClosedSubshells(), wa[1]);                                     @show wb
    wc = Basics.extractFromConfiguration(OpenSubshells(), wa[1]);                                       @show wc
    conf = Configuration("[Ar] 3d^2 4s 4p")
    we = Basics.extractFromConfiguration(OpenShells(), conf);                                           @show we
    wf = Basics.extractFromConfiguration(OpenSubshells(), conf);                                        @show wf
    wg = Basics.extractFromConfiguration(TotalAM(), conf);                                              @show wg
    wh = Basics.extractFromConfiguration(ClosedSubshells(), "[Ne]");                                    @show wh
    wj = Basics.extractFromConfiguration(ClosedSubshells(), "[Xe]");                                    @show wj
    #
elseif  false
    # Extract parity, is-occupied and multiplicity from a single standard or relativistic configuration
    # -------------------------------------------------------------------------------------------------
    conf = Configuration("[Ar] 4p");   coreConf = Configuration("[Ne]");   innerConf = Configuration("[He] 2s^2 2p^4 3s 3p^6")
    wa = Basics.generateConfigurations(RelativisticConfigurations(), conf);                             @show wa
    wb = Basics.extractFromConfiguration(ClosedCore(), conf);                                           @show wb
    wc = Basics.extractFromConfiguration(GetParity(), conf);                                            @show wc
    wd = Basics.extractFromConfiguration(GetParity(), wa[1]);                                           @show wd
    we = Basics.extractFromConfiguration(GetParity(), wa[2]);                                           @show we
    wf = Basics.extractFromConfiguration(IsOccupied(), conf, Shell("3s"));                              @show wf
    wg = Basics.extractFromConfiguration(IsOccupied(), conf, Shell("5p"));                              @show wg
    wh = Basics.extractFromConfiguration(TotalAM(true, AngularJ64[]), conf);                            @show wh
    wj = Basics.extractFromConfiguration(TotalAM(false, AngularJ64[]), conf);                           @show wj
    wk = Basics.extractFromConfiguration(Multiplicity(), conf);                                         @show wk
    wl = Basics.extractFromConfiguration(ValenceOccupation(), conf, coreConf);                          @show wl
    #
elseif  false
    # Generate relativistic configurations and extract information from them.
    # -----------------------------------------------------------------------
    conf = Configuration("[Ar] 4s 4p")
    wa = Basics.extractFromConfiguration(ClosedShells(), xeConf);                                       @show wa
    wb = Basics.extractFromConfiguration(ClosedSubshells(), xeConf);                                    @show wb
    wc = Basics.extractFromConfigurations(ClosedShells(), allConfs);                                    @show wc
    #
elseif  false
    # Extract information from a list of configurations    
    # -------------------------------------------------
    fromShells = [Shell("3s"), Shell("3p")]
    intoShells = Basics.generateShellList(4, 5, 3)
    wa = Basics.generateConfigurations(ForRasExcitations(true, true, false, false, fromShells, intoShells), [arConf])
         Basics.displayConfigurations(stdout, wa, longForm=true, details="ForRasExcitations")
    wb = Basics.extractFromConfigurations(MeanOccupation(), wa);                                        @show wb
    wc = Basics.extractFromConfigurations(MeanOccupation(), [xeConf]);                                  @show wc
    wd = Basics.extractFromConfigurations(NumberOfElectrons(), wa);                                     @show wd
    we = Basics.extractFromConfigurations(ClosedCore(), wa);                                            @show we
    wf = Basics.extractFromConfigurations(OccupationDifference(), wa[1], wa[2]);                        @show wf
    wg = Basics.extractFromConfigurations(OccupationDifference(), wa[1], wa[11]);                       @show wg
    conf = Configuration("[Ar] 4p")
    wh = Basics.generateConfigurations(RelativisticConfigurations(), conf);                             @show wh
    wj = Basics.extractFromConfigurations(OccupationDifference(), wh[1], wh[2]);                        @show wj
    #
elseif  false
    # Display configurations with additional information 
    # --------------------------------------------------
    fromShells = Basics.generateShellList(5, 5, 2)
    intoShells = Basics.generateShellList(6, 7, 5)
    wa = Basics.generateConfigurations(AddElectrons(2, intoShells), [xeConf])
    Basics.displayConfigurations(stdout, wa, longForm=true)
    Basics.displayConfigurations(stdout, wa, details = "xenon")
    #
    Basics.displayConfiguration(stdout, Basics.MeanConfiguration(), wa, details = "xenon")
    Basics.displayConfiguration(stdout, Basics.FineStructure(), wa[15], details = "xenon")
    Basics.displayConfiguration(stdout, Basics.FineStructure(), wa[45], details = "xenon")
    Basics.displayConfiguration(stdout, Basics.FineStructure(), wa[75], details = "xenon")
    #
end
    

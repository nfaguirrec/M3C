program M3C
	use GOptions_
	use String_
	use CommandLineParser_
	use BlocksIFileParser_
	
	use GOptionsM3C_
	use FragmentsDB_
	use FragmentsList_
	use Reactor_
	use MarkovChain_
	
! 	use RMoleculeListJL_
! 	use RMJLReactor_
! 	use RMJL_
	
	implicit none
	
	type(CommandLineParser) :: programOptions
	type(BlocksIFileParser) :: iParser
	
	type(String) :: inputFileName
	
	type(FragmentsList) :: reactives
	
! 	type(Reactor) :: reactorJL
! 	type(MarkovChain) :: RMJLMethod
	
	type(Reactor) :: reactorMethod
	type(MarkovChain) :: MarkovChainMethod
	
	logical :: calculateMaxVib
	
	integer :: elapsedTime(3)
	integer :: status
	character(128) :: hostname
	
	!------------------------------------------------------------------------------
	
! 	call RigidMolecule_test()
! 	call RMoleculeListJL_test()
! 	call FragmentsDB_test()
! 	call RMJLReactor_test()
! 	stop

	write(*,*) "--------------------------------------"
	write(*,*) "M3C version: ", __DATE__, " - ", __TIME__
	write(*,*) "--------------------------------------"
	
	call programOptions.init()
! 	status = hostnm(hostname)
	
	inputFileName = programOptions.getString( "-i" )
	write(*,*) "Input file = ", trim(inputFileName.fstr)
	write(*,*) ""
	
	call iParser.init( inputFileName.fstr )
	call iParser.showContent()
	
	call GOptions_timer.init()
	call GOptions_timer.start()
	write(*,*) "---------------------------------------------------------------------------"
! 	write(*,*) "START TIME: ", trim(GOptions_timer.startDate()), "  @"//trim(hostname)
	write(*,*) "START TIME: ", trim(GOptions_timer.startDate())
	write(*,*) "---------------------------------------------------------------------------"
	
	!-------------------------------------------------------------------
	! Loading the global options
	!-------------------------------------------------------------------
	GOptions_zero = iParser.getReal( "GOPTIONS:zero", def=1d-12 )
	
	GOptionsM3C_systemRadius = iParser.getReal( "GOPTIONS:systemRadius", def=10.0_8*angs )*angs
	GOptionsM3C_randomWalkStepRadius = iParser.getReal( "GOPTIONS:randomWalkStepRadius", def=2.0_8*angs )*angs
	GOptionsM3C_useWeightedWalkStep = iParser.getLogical( "GOPTIONS:useWeightedWalkStep", def=.false. )
	GOptionsM3C_overlappingRadius = iParser.getReal( "GOPTIONS:overlappingRadius", def=0.0_8 )*angs
	GOptionsM3C_useZPECorrection = iParser.getLogical( "GOPTIONS:useZPECorrection", def=.false. )
	GOptionsM3C_useRandomWalkers = iParser.getLogical( "GOPTIONS:useRandomWalkers", def=.true. )
	GOptionsM3C_useLCorrection = iParser.getLogical( "GOPTIONS:useLCorrection", def=.false. )
	GOptionsM3C_useLDOSContrib = iParser.getLogical( "GOPTIONS:useLDOSContrib", def=.false. )
	GOptionsM3C_useLWeightContrib = iParser.getLogical( "GOPTIONS:useLWeightContrib", def=.false. )
	GOptionsM3C_useLReference = iParser.getLogical( "GOPTIONS:useLReference", def=.false. )
	GOptionsM3C_gammaLCorrection = iParser.getReal( "GOPTIONS:gammaLCorrection", def=3.0_8 )
	GOptionsM3C_useSpinConservationRules = iParser.getLogical( "GOPTIONS:useSpinConservationRules", def=.false. )
	GOptions_printLevel = iParser.getInteger( "GOPTIONS:printLevel", def=1 )
	GOptions_debugLevel = iParser.getInteger( "GOPTIONS:debugLevel", def=1 )
	
	write(*,*)
	write(*,"(A40,E15.2)") "GOptions:zero = ", GOptions_zero
	write(*,"(A40,F15.5,A)") "GOptions:systemRadius = ", GOptionsM3C_systemRadius/angs, " A"
	write(*,"(A40,F15.5,A)") "GOptions:randomWalkStepRadius = ", GOptionsM3C_randomWalkStepRadius/angs, " A"
	write(*,"(A40,F15.5,A)") "GOptions:overlappingRadius = ", GOptionsM3C_overlappingRadius/angs, " A"
	write(*,"(A40,L5)") "GOptions:useWeightedWalkStep = ", GOptionsM3C_useWeightedWalkStep
	write(*,"(A40,L5)") "GOptions:useLCorrection = ", GOptionsM3C_useLCorrection
	write(*,"(A40,L5)") "GOptions:useLDOSContrib = ", GOptionsM3C_useLDOSContrib
	write(*,"(A40,L5)") "GOptions:useLWeightContrib = ", GOptionsM3C_useLWeightContrib
	write(*,"(A40,L5)") "GOptions:useLReference = ", GOptionsM3C_useLReference
	write(*,"(A40,F5.0)") "GOptions:gammaLCorrection = ", GOptionsM3C_gammaLCorrection
	write(*,"(A40,L)") "useSpinConservationRules = ", GOptionsM3C_useSpinConservationRules
	write(*,"(A40,I5)") "GOptions:printLevel = ", GOptions_printLevel
	write(*,"(A40,I5)") "GOptions:debugLevel = ", GOptions_debugLevel
	write(*,*)
	
	calculateMaxVib = iParser.getLogical( "FRAGMENTS_DATABASE"//":maxVib", def=.false. )
	
	!-------------------------------------------------------------------
	! Loading the mass table
	!-------------------------------------------------------------------
	call FragmentsDB_instance.fromInputFile( iParser )
	
	!-------------------------------------------------------------------
	! Checking for System radius optimization block
	!-------------------------------------------------------------------
	call reactives.executeRadiusOptimization( iParser )
	
	!-------------------------------------------------------------------
	! Checking for System radius optimization block
	!-------------------------------------------------------------------
! 	call reactives.executeRVCurve( iParser )
		
	!-------------------------------------------------------------------
	! Checking for Reactor test block
	!-------------------------------------------------------------------
! 	call reactorJL.execute( iParser )
		
	!-------------------------------------------------------------------
	! Checking for MarkovChain calculation
	!-------------------------------------------------------------------
! 	call RMJLMethod.execute( iParser )
	
	!-------------------------------------------------------------------
	! Checking for Reactor test block
	!-------------------------------------------------------------------
	if( iParser.isThereBlock( "REACTOR" ) ) then
		call reactorMethod.execute( iParser )
		
	else if( calculateMaxVib ) then
		call reactorMethod.executeMinFragmentationEnergy( iParser )
! 		call reactorMethod.executeGenerateAllChannels( iParser )
		
	!-------------------------------------------------------------------
	! Checking for MarkovChain calculation
	!-------------------------------------------------------------------
	else if( iParser.isThereBlock( "MARKOV_CHAIN" ) ) then
		call MarkovChainMethod.execute( iParser )
		
	end if
	
	elapsedTime = GOptions_timer.elapsed()
	write(*,*) "---------------------------------------------------------------------------"
	write(*,*) "ELAPSED TIME: ", elapsedTime(1), "h", elapsedTime(2), "m", elapsedTime(3), "s"
	write(*,*) "END TIME: ", trim(GOptions_timer.currentDate())
	write(*,*) "---------------------------------------------------------------------------"
end program M3C

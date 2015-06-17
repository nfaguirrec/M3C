program M3C
	use String_
	use CommandLineParser_
	use BlocksIFileParser_
	
	use GOptions_
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
	
	!------------------------------------------------------------------------------
	
! 	call RigidMolecule_test()
! 	call RMoleculeListJL_test()
! 	call FragmentsDB_test()
! 	call RMJLReactor_test()
! 	stop
	
	call programOptions.init()
	
	inputFileName = programOptions.getString( "-i" )
	write(*,*) "Input file = ", trim(inputFileName.fstr)
	write(*,*) ""
	
	call iParser.init( inputFileName.fstr )
	call iParser.showContent()
	
	call GOptions_timer.init()
	call GOptions_timer.start()
	write(*,*) "---------------------------------------------------------------------------"
	write(*,*) "START TIME: ", trim(GOptions_timer.startDate())
	write(*,*) "---------------------------------------------------------------------------"
	
	!-------------------------------------------------------------------
	! Loading the global options
	!-------------------------------------------------------------------
	GOptions_zero = iParser.getReal( "GOPTIONS:zero", def=1d-12 )
	
	GOptions_systemRadius = iParser.getReal( "GOPTIONS:systemRadius", def=10.0_8*angs )*angs
	GOptions_randomWalkStepRadius = iParser.getReal( "GOPTIONS:randomWalkStepRadius", def=2.0_8*angs )*angs
	GOptions_useWeightedWalkStep = iParser.getLogical( "GOPTIONS:useWeightedWalkStep", def=.false. )
	GOptions_overlappingRadius = iParser.getReal( "GOPTIONS:overlappingRadius", def=0.0_8 )*angs
	GOptions_useZPECorrection = iParser.getLogical( "GOPTIONS:useZPECorrection", def=.false. )
	GOptions_useRandomWalkers = iParser.getLogical( "GOPTIONS:useRandomWalkers", def=.true. )
	GOptions_useLCorrection = iParser.getLogical( "GOPTIONS:useLCorrection", def=.false. )
	GOptions_useLDOSContrib = iParser.getLogical( "GOPTIONS:useLDOSContrib", def=.false. )
	GOptions_useLWeightContrib = iParser.getLogical( "GOPTIONS:useLWeightContrib", def=.false. )
	GOptions_useLReference = iParser.getLogical( "GOPTIONS:useLReference", def=.false. )
	GOptions_gammaLCorrection = iParser.getReal( "GOPTIONS:gammaLCorrection", def=3.0_8 )
	GOptions_printLevel = iParser.getInteger( "GOPTIONS:printLevel", def=1 )
	GOptions_debugLevel = iParser.getInteger( "GOPTIONS:debugLevel", def=1 )
	
	write(*,*)
	write(*,"(A40,E15.2)") "GOptions:zero = ", GOptions_zero
	write(*,"(A40,F15.5,A)") "GOptions:systemRadius = ", GOptions_systemRadius/angs, " A"
	write(*,"(A40,F15.5,A)") "GOptions:randomWalkStepRadius = ", GOptions_randomWalkStepRadius/angs, " A"
	write(*,"(A40,F15.5,A)") "GOptions:overlappingRadius = ", GOptions_overlappingRadius/angs, " A"
	write(*,"(A40,L5)") "GOptions:useWeightedWalkStep = ", GOptions_useWeightedWalkStep
	write(*,"(A40,L5)") "GOptions:useLCorrection = ", GOptions_useLCorrection
	write(*,"(A40,L5)") "GOptions:useLDOSContrib = ", GOptions_useLDOSContrib
	write(*,"(A40,L5)") "GOptions:useLWeightContrib = ", GOptions_useLWeightContrib
	write(*,"(A40,L5)") "GOptions:useLReference = ", GOptions_useLReference
	write(*,"(A40,F5.0)") "GOptions:gammaLCorrection = ", GOptions_gammaLCorrection
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

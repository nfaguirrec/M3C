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
	
	implicit none
	
	type(CommandLineParser) :: programOptions
	type(BlocksIFileParser) :: iParser
	
	type(String) :: inputFileName
	
	type(FragmentsList) :: reactives
	
	type(Reactor) :: reactorMethod
	type(MarkovChain) :: MarkovChainMethod
	
	logical :: calculateMaxVib
	
	integer :: elapsedTime(3)
	integer :: status
	character(128) :: hostname
	
	type(String) :: sBuffer
	
	!------------------------------------------------------------------------------
	
	call programOptions.init()
! 	status = hostnm( hostname )
	
	write(*,*) "--------------------------------------"
	write(*,*) "M3C version: ", __DATE__, " - ", __TIME__
! 	write(*,*) "   hostname: "//trim(hostname)
	write(*,*) "--------------------------------------"
	
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
	
	GOptionsM3C_systemRadius = iParser.getReal( "GOPTIONS:systemRadius", def=10.0_8*angs )*angs
	GOptionsM3C_randomWalkStepRadius = iParser.getReal( "GOPTIONS:randomWalkStepRadius", def=2.0_8*angs )*angs
	GOptionsM3C_useRandomWalkers = iParser.getLogical( "GOPTIONS:useRandomWalkers", def=.false. )
	GOptionsM3C_useWeightedWalkStep = iParser.getLogical( "GOPTIONS:useWeightedWalkStep", def=.false. )
	GOptionsM3C_overlappingRadius = iParser.getReal( "GOPTIONS:overlappingRadius", def=0.0_8 )*angs
	
	sBuffer = iParser.getString( "GOPTIONS:radiusType", def="COVALENT" )
	select case( trim(sBuffer.fstr) )
		case( "COVALENT" )
			GOptionsM3C_radiusType = AtomicElementsDB_COVALENT_RADIUS
		case( "VAN_DER_WAALS" )
			GOptionsM3C_radiusType = AtomicElementsDB_VANDERWAALS_RADIUS
		case default
			call GOptions_error( "GOPTIONS:radiusType = "//trim(sBuffer.fstr)//" is not implemented. Available values: COVALENT, VAN_DER_WAALS", "M3C.main()" )
	end select
	
	GOptionsM3C_useZPECorrection = iParser.getLogical( "GOPTIONS:useZPECorrection", def=.false. )
	GOptionsM3C_useSpinConservationRules = iParser.getLogical( "GOPTIONS:useSpinConservationRules", def=.false. )
	GOptionsM3C_angularMomentumCouplingScheme = iParser.getString( "GOPTIONS:angularMomentumCouplingScheme", def="JJ" )
	GOptionsM3C_structureSamplingMethod = iParser.getString( "GOPTIONS:structureSamplingMethod", def="RANDOM" )
	GOptionsM3C_totalJ(3) = iParser.getReal( "GOPTIONS:totalJ", def=0.0_8 )
	GOptions_printLevel = iParser.getInteger( "GOPTIONS:printLevel", def=1 )
	GOptions_debugLevel = iParser.getInteger( "GOPTIONS:debugLevel", def=1 )
	
	write(*,*)
	write(*,"(A40,E15.2)") "zero = ", GOptions_zero
	write(*,"(A40,F15.5,A)") "systemRadius = ", GOptionsM3C_systemRadius/angs, " A"
	write(*,"(A40,L5)") "useRandomWalkers = ", GOptionsM3C_useRandomWalkers
	write(*,"(A40,F15.5,A)") "randomWalkStepRadius = ", GOptionsM3C_randomWalkStepRadius/angs, " A"
	write(*,"(A40,L5)") "useWeightedWalkStep = ", GOptionsM3C_useWeightedWalkStep
	write(*,"(A40,F15.5,A)") "overlappingRadius = ", GOptionsM3C_overlappingRadius/angs, " A"
	
	select case( GOptionsM3C_radiusType )
		case( AtomicElementsDB_COVALENT_RADIUS )
			write(*,"(A40,A15)") "radiusType = ", "COVALENT"
		case( AtomicElementsDB_VANDERWAALS_RADIUS )
			write(*,"(A40,A15)") "radiusType = ", "VAN_DER_WAALS"
	end select
	
	write(*,"(A40,L5)") "useZPECorrection = ", GOptionsM3C_useZPECorrection
	write(*,"(A40,L)") "useSpinConservationRules = ", GOptionsM3C_useSpinConservationRules
	write(*,"(A40,A5)") "angularMomentumCouplingScheme = ", trim(GOptionsM3C_angularMomentumCouplingScheme.fstr)
	write(*,"(A40,A10)") "structureSamplingMethod = ", trim(GOptionsM3C_structureSamplingMethod.fstr)
	write(*,"(A40,3F15.5,A)") "totalJ = ", GOptionsM3C_totalJ, " a.u."
	write(*,"(A40,I5)") "printLevel = ", GOptions_printLevel
	write(*,"(A40,I5)") "debugLevel = ", GOptions_debugLevel
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

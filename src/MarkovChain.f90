!>
!! @brief
!!
module MarkovChain_
	use String_
	use StringList_
	use RealList_
	use RandomUtils_
	use UnitsConverter_
	use StringIntegerPair_
	use StringRealPair_
	use StringIntegerMap_
	use StringRealMap_
	use RealHistogram_
	use StringHistogram_
	use BlocksIFileParser_
	
	use GOptions_
	use FragmentsList_
	use FragmentsDB_
	use Reactor_
	
	implicit none
	private
	
	type, public :: MarkovChain
		type(String) :: task
		real(8) :: burnInFraction
		integer :: numberOfEvents
		integer :: numberOfExperiments
		real(8) :: excitationEnergy
		integer :: historyFileFrequency = -1
		type(String) :: geometryHistoryFilePrefix
		integer :: freqBlockingCheck = 4
		type(String) :: tracking
		
		type(RealHistogram), allocatable :: iTemperatureHistogram(:) ! One item for each experiment
		type(RealHistogram), allocatable :: entropyHistogram(:)      ! One item for each experiment
		
		type(StringHistogram) :: reactorAcceptedHistogram    ! One for all experiments
		type(StringHistogram) :: reactorRejectedHistogram    ! One for all experiments
		type(StringHistogram) :: reactorStatusHistogram      ! One for all experiments
		
		type(StringHistogram) :: transitionHistogram         ! One for all experiments. Probability to change from one channel to another
		type(StringHistogram) :: transitionDetHistogram      ! One for all experiments. Probability to change from one channel to another
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Histogramas
		type(StringIntegerMap), allocatable :: speciesHistogram(:) ! One item for each experiment
		type(StringIntegerMap), allocatable :: speciesDetHistogram(:) ! One item for each experiment
		
		type(StringIntegerMap), allocatable :: channelHistogram(:) ! One item for each experiment
		type(StringIntegerMap), allocatable :: channelDetHistogram(:) ! One item for each experiment
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Strings de history
		type(StringList), private :: energyHistory
		type(StringList), private :: weightHistory
		type(StringList), private :: JHistory
		type(StringList), private :: LHistory
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		contains
			generic :: init => initDefault
			generic :: assignment(=) => copy
			
			procedure :: initDefault
			procedure :: copy
			final :: destroy
			procedure, private :: clearHistograms
			procedure, private :: initHistograms
			procedure, NOPASS :: showAverHistogram
			procedure :: run
			procedure :: saveEnergyHistory
			procedure :: saveWeightHistory
			procedure :: saveJHistory
			procedure :: saveLHistory
			procedure :: saveHistograms
			
			procedure :: execute
	end type MarkovChain
	
	contains
	
	!>
	!! @brief Constructor
	!!
	subroutine initDefault( this )
		class(MarkovChain) :: this
		
	end subroutine initDefault
	
	!>
	!! @brief Copy constructor
	!!
	subroutine copy( this, other )
		class(MarkovChain), intent(out) :: this
		class(MarkovChain), intent(in) :: other
		
		call GOptions_error( &
			"This method is not implemented", &
			"Reactor.copyReactor()" &
		)
	end subroutine copy
	
	!>
	!! @brief Destructor
	!!
	subroutine destroy( this )
		type(MarkovChain) :: this
		
		call this.clearHistograms()
	end subroutine destroy
	
	!>
	!! @brief
	!!
	subroutine clearHistograms( this )
		class(MarkovChain) :: this
		
		integer :: nExp
		
		if( allocated(this.iTemperatureHistogram) ) then
			deallocate( this.iTemperatureHistogram )
		end if
		
		if( allocated(this.entropyHistogram) ) then
			deallocate( this.entropyHistogram )
		end if
		
		call this.reactorAcceptedHistogram.clear()
		call this.reactorRejectedHistogram.clear()
		call this.reactorStatusHistogram.clear()
		
		call this.transitionHistogram.clear()
		call this.transitionDetHistogram.clear()
		
		if( allocated(this.speciesHistogram) ) then
			do nExp=1,this.numberOfExperiments
				call this.speciesHistogram(nExp).clear()
			end do
			
			deallocate( this.speciesHistogram )
		end if
		
		if( allocated(this.speciesHistogram) ) then
			do nExp=1,this.numberOfExperiments
				call this.speciesDetHistogram(nExp).clear()
			end do
			
			deallocate( this.speciesDetHistogram )
		end if
			
		if( allocated(this.channelHistogram) ) then
			do nExp=1,this.numberOfExperiments
				call this.channelHistogram(nExp).clear()
			end do
			
			deallocate( this.channelHistogram )
		end if
			
		if( allocated(this.channelDetHistogram) ) then
			do nExp=1,this.numberOfExperiments
				call this.channelDetHistogram(nExp).clear()
			end do
			
			deallocate( this.channelDetHistogram )
		end if
	end subroutine clearHistograms
	
	!>
	!! @brief
	!!
	subroutine initHistograms( this )
		class(MarkovChain) :: this
		
		integer :: nExp
		
		allocate( this.iTemperatureHistogram(this.numberOfExperiments) )
		allocate( this.entropyHistogram(this.numberOfExperiments) )
		allocate( this.speciesHistogram(this.numberOfExperiments) )
		allocate( this.speciesDetHistogram(this.numberOfExperiments) )
		allocate( this.channelHistogram(this.numberOfExperiments) )
		allocate( this.channelDetHistogram(this.numberOfExperiments) )
		
		call this.reactorAcceptedHistogram.initStringHistogram()
		call this.reactorRejectedHistogram.initStringHistogram()
		call this.reactorStatusHistogram.initStringHistogram()
		
		call this.transitionHistogram.initStringHistogram()
		call this.transitionDetHistogram.initStringHistogram()
		
		do nExp=1,this.numberOfExperiments
			call this.iTemperatureHistogram(nExp).initRealHistogram()
			call this.entropyHistogram(nExp).initRealHistogram()
			call this.speciesHistogram(nExp).init()
			call this.speciesDetHistogram(nExp).init()
			call this.channelHistogram(nExp).init()
			call this.channelDetHistogram(nExp).init()
		end do
	end subroutine initHistograms
	
	!>
	!! @brief
	!!
	subroutine run( this, reactives )
		class(MarkovChain) :: this
		type(FragmentsList), intent(in) :: reactives
		
		type(Reactor) :: react
		
		integer :: nTimesBlocked ! si se bloquea 2 veces consecutivas fuerza a random
		integer :: historyFrequency
		character(10), allocatable :: taskTokens(:)
		character(10) :: currentTask
		real(8) :: p, Pi
		integer :: n, i, j, nExp, iBuffer
		type(String) :: sBuffer
		character(100) :: geometryFileName
		character(2) :: origin
		
		call this.energyHistory.clear()
		call this.weightHistory.clear()
		
		call this.clearHistograms()
		call this.initHistograms()
		
		call this.task.split( taskTokens, "," )
		
		if( this.historyFileFrequency == -1 ) then
			historyFrequency = max(this.numberOfEvents,this.numberOfEvents/1000)
		else
			historyFrequency = this.historyFileFrequency 
		end if
		
		if( this.tracking /= "energy" .and. this.tracking /= "weight" .and. this.tracking /= "none" ) then
			call GOptions_error( &
				"This tracking variable is not available ("//trim(this.tracking.fstr)//")", &
				"MMMCMethod.run()", &
				"The tracing variables available are: energy, weight, none" &
				)
		end if
		
		if( this.tracking == "energy" ) then
		
			write(6,"(A)") "#------------------------------------"
			write(6,"(A)") "# ENERGY HISTORY"
			write(6,"(A)") "#------------------------------------"
			write(6,"(A1,3X,5A15,5X,A)") "#", "trans", "intermol", "vib", "rot", "tot", "formula"
			write(6,"(A1,3X,5A15,5X,A)") "#", "eV", "eV", "eV", "eV", ""
			write(6,"(A1,3X,5A15,5X,A)") "#", "-------", "--------", "-------", "-----", "-------"

		else if( this.tracking == "weight" ) then
		
			write(6,"(A)") "#------------------------------------"
			write(6,"(A)") "# WEIGHT HISTORY"
			write(6,"(A)") "#------------------------------------"
			write(6,"(A1,3X,6A15,5X,A)") "#", "LnWe", "LnWv", "LnWn", "LnWr", "LnWt", "LnW", "formula"
			write(6,"(A1,3X,6A15,5X,A)") "#", "arb.", "arb.", "arb.", "arb.", "arb.", "arb.", ""
			write(6,"(A1,3X,6A15,5X,A)") "#", "-------", "-------", "--------", "--------", "--------", "-------", "-------"
			
		end if
		
! 		call reactives.initialGuessFragmentsList()

		do nExp=1,this.numberOfExperiments
		
			nTimesBlocked = 0
			
			call react.init( reactives, this.excitationEnergy )
			
			if( .not. this.geometryHistoryFilePrefix.isEmpty() ) then
				geometryFileName = trim(this.geometryHistoryFilePrefix.fstr)//"-"//trim(FString_fromInteger(nExp))//".xyz"
				call react.reactives.save(geometryFileName)
			end if
			
			n = 1
			do while( n <= this.numberOfEvents )
				do i=1,size(taskTokens)
					if( n == this.numberOfEvents+1 ) exit
					
					currentTask = trim(adjustl(taskTokens(i)))
					call react.setType( currentTask )
					
! 					if( n==1 ) then
! 						if( this.tracking == "energy" ) then
! 							sBuffer = react.reactives.energyHistoryLine( origin )
! 						else if( this.tracking == "weight" ) then
! 							sBuffer = react.reactives.weightHistoryLine( origin )
! 						end if
! 					end if

					call react.run()
					
					! Si la energía cinetica es negativa
					if( .not. react.state ) then
						
! 						if( n == 1 ) then
! 							call react.reactives.initialGuessFragmentsList()
! 							write(*,*) "Trying another initial guess"
! 							cycle
! 						end if
						
						! Se fuerza el uso de RigidMoleculeList.randomCenters() solo para la siguiente llamada
						! a RigidMoleculeList.changeGeometry()
						if( nTimesBlocked == this.freqBlockingCheck ) then
							call GOptions_info( "Blocked. Forcing centers completely random", "MarkovChain" )
							
							react.reactives.forceRandomCenters = .true.
							nTimesBlocked = 0
						else
							nTimesBlocked = nTimesBlocked + 1
						end if
							
						call GOptions_info( "Step rejected ( Negative energy )", "MarkovChain" )
						
						origin = "e"//trim(currentTask)
						
						call this.reactorRejectedHistogram.append( FString_toString( trim(currentTask) ) )
						call this.reactorStatusHistogram.append( FString_toString( "e.REJECTED(E<0) " ) )
					else
					
						nTimesBlocked = 0  ! El bloqueo debe ser consecutivo, así que si no pasa por negative energy se cuenta nuevamente
						
						Pi = react.products.LnW()-react.reactives.LnW()
						
						if( Pi > 0.0_8 ) then
							if( trim(react.reactives.label( details=.false. )) /= trim(react.products.label( details=.false. )) ) then
								sBuffer = trim(react.reactives.label( details=.false. ))//"-->"//trim(react.products.label( details=.false. ))
								call this.transitionHistogram.add( sBuffer )
							end if
						
							if( trim(react.reactives.label( details=.true. )) /= trim(react.products.label( details=.true. )) ) then
								sBuffer = trim(react.reactives.label( details=.true. ))//"-->"//trim(react.products.label( details=.true. ))
								call this.transitionDetHistogram.add( sBuffer )
							end if
							
							react.reactives = react.products
							
							call GOptions_info( &
							  "Step accepted logPi="//trim(adjustl(FString_fromReal(Pi,"(F10.3)")))//"  ( logPi > 0 )", "MarkovChain" )
							
							origin = "a"//trim(currentTask)
							
							call this.reactorAcceptedHistogram.append( FString_toString( trim(currentTask) ) )
							call this.reactorStatusHistogram.append( FString_toString( "a.ACCEPTED      " ) )
						else
							p = log( RandomUtils_uniform( [0.0_8, 1.0_8] ) )
							
							if( p <= Pi ) then
								if( trim(react.reactives.label( details=.false. )) /= trim(react.products.label( details=.false. )) ) then
									sBuffer = trim(react.reactives.label( details=.false. ))//"-->"//trim(react.products.label( details=.false. ))
									call this.transitionHistogram.add( sBuffer )
								end if
							
								if( trim(react.reactives.label( details=.true. )) /= trim(react.products.label( details=.true. )) ) then
									sBuffer = trim(react.reactives.label( details=.true. ))//"-->"//trim(react.products.label( details=.true. ))
									call this.transitionDetHistogram.add( sBuffer )
								end if

								react.reactives = react.products
								
								call GOptions_info( &
								  "Step accepted logP="//trim(adjustl(FString_fromReal(p,"(F10.3)")))// &
								  ",logPi="//trim(adjustl(FString_fromReal(Pi,"(F10.3)")))//"  ( logP <= logPi )", "MarkovChain" )
								
								origin = "p"//trim(currentTask)
								
								call this.reactorAcceptedHistogram.append( FString_toString( trim(currentTask) ) )
								call this.reactorStatusHistogram.append( FString_toString( "p.ACCEPTED(p<PI)" ) )
							else
								call GOptions_info( &
								  "Step rejected logP="//trim(adjustl(FString_fromReal(p,"(F10.3)")))// &
								  ",logPi="//trim(adjustl(FString_fromReal(Pi,"(F10.3)")))//"  ( logP > logPi )", "MarkovChain" )
								
								origin = "r"//trim(currentTask)
								
								call this.reactorRejectedHistogram.append( FString_toString( trim(currentTask) ) )
								call this.reactorStatusHistogram.append( FString_toString( "r.REJECTED      " ) )
							end if
							
						end if
						
					end if
					
						
					if( n > this.burnInFraction*this.numberOfEvents ) then
					
						if( mod(n,historyFrequency) == 0 ) then
						
							if( this.tracking == "energy" ) then
								sBuffer = react.reactives.energyHistoryLine( origin )
							else if( this.tracking == "weight" ) then
								sBuffer = react.reactives.weightHistoryLine( origin )
							end if
								
							if( this.tracking /= "none" ) then
								write(6,"(A)") trim(sBuffer.fstr)
							end if
							
							call this.weightHistory.append( react.reactives.weightHistoryLine( origin ) )
							call this.energyHistory.append( react.reactives.energyHistoryLine( origin ) )
							
							if( .not. this.geometryHistoryFilePrefix.isEmpty() ) then
								geometryFileName = trim(this.geometryHistoryFilePrefix.fstr)//"-"//trim(FString_fromInteger(nExp))//".xyz"
								call react.reactives.save(geometryFileName, append=.true.)
							end if
							
						end if
					
						call this.iTemperatureHistogram(nExp).append( react.reactives.iTemperature() )
						call this.entropyHistogram(nExp).append( react.reactives.LnW() )
						
						do j=1,react.reactives.nMolecules()
							sBuffer = react.reactives.clusters(j).label( details=.false. )
							call this.speciesHistogram(nExp).set( sBuffer, this.speciesHistogram(nExp).at( sBuffer, defaultValue=0 )+1 )
							
							sBuffer = react.reactives.clusters(j).label( details=.true. )
							call this.speciesDetHistogram(nExp).set( sBuffer, this.speciesDetHistogram(nExp).at( sBuffer, defaultValue=0 )+1 )
						end do
						
						call this.JHistory.append( react.reactives.JHistoryLine() )
						call this.LHistory.append( react.reactives.LHistoryLine() )
						
						sBuffer = react.reactives.label( details=.false. )
						call this.channelHistogram(nExp).set( sBuffer, this.channelHistogram(nExp).at( sBuffer, defaultValue=0 )+1 )
						
						sBuffer = react.reactives.label( details=.true. )
						call this.channelDetHistogram(nExp).set( sBuffer, this.channelDetHistogram(nExp).at( sBuffer, defaultValue=0 )+1 )
						
					end if
					
					n = n+1
				end do
			end do
			
			sBuffer = ""
! 			if( this.tracking == "energy" ) then
				
				call this.weightHistory.append( sBuffer )
				call this.weightHistory.append( sBuffer )
				
! 			else if( this.tracking == "weight" ) then
			
				call this.energyHistory.append( sBuffer )
				call this.energyHistory.append( sBuffer )
				
! 			end if
			write(6,"(A)") ""
			write(6,"(A)") ""
			
			call this.JHistory.append( sBuffer )
			call this.JHistory.append( sBuffer )
			
			call this.LHistory.append( sBuffer )
			call this.LHistory.append( sBuffer )
		end do
		
		if( GOptions_printLevel >= 1 ) then
			write(*,"(A)") ""
		end if
		
		if( this.tracking == "energy" ) then
			call this.saveWeightHistory()
		else if( this.tracking == "weight" ) then
			call this.saveEnergyHistory()
		end if
			
		call this.saveHistograms()
		
	end subroutine run
	
	!>
	!! @brief
	!!
	subroutine showAverHistogram( histogramVec, unit, ebklFileName, excitationEnergy )
		type(StringIntegerMap), intent(in) :: histogramVec(:)
		integer, optional, intent(in) :: unit
		character(*), optional, intent(in) :: ebklFileName
		real(8), optional, intent(in) :: excitationEnergy
		
		integer :: effUnit
		
		type(StringIntegerMap) :: averHistogram
		class(StringIntegerMapIterator), pointer :: iter, ptr
		integer :: n, i, nExp, numberOfExperiments
		type(String) :: sBuffer
		type(StringIntegerPair) :: pair
		real(8) :: aver, desv, stdError
		type(OFStream) :: ebklFile
		real(8) :: normConst
		integer :: maxNameLength
		
		effUnit = 6
		if( present(unit) ) effUnit = unit
		
		! Abre el fichero para almacenar el bloque de energía si así fué seleccionado
		if( present(ebklFileName) .and. present(excitationEnergy) ) then
			call ebklFile.open( ebklFileName )
		end if
		
		call averHistogram.init()
		numberOfExperiments = size(histogramVec)
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! averHistogram corresponde a un histograma generado
		! por la totalidad de los datos
		normConst = 0.0_8
		maxNameLength = 0
		do nExp=1,numberOfExperiments
			iter => histogramVec(nExp).begin
			do while( associated(iter) )
				pair = histogramVec(nExp).pair( iter )
				call averHistogram.set( pair.first, 1 )
				
				if( nExp == 1 ) then
					normConst = normConst + real(pair.second, 8)
				end if
				
				if( pair.first.length() > maxNameLength ) then
					maxNameLength = pair.first.length()
				end if
				
				iter => iter.next
			end do
		end do
		
		normConst = 1.0_8/normConst
		
		write(effUnit,"(A1,10X,A<maxNameLength-1>,<numberOfExperiments>I10,5X,2A10)") "#", "item", ( i, i=1,numberOfExperiments ), "aver", "desv"
		write(effUnit,"(A1,10X,A<maxNameLength-1>,<numberOfExperiments>A10,5X,2A10)") "#", "----", ( "-----", i=1,numberOfExperiments ), "----", "----"
		
		iter => averHistogram.begin
		do while( associated(iter) )
			pair = averHistogram.pair( iter )
			
			aver = 0.0_8
			do nExp=1,numberOfExperiments
				sBuffer = pair.first
				if( histogramVec(nExp).find( sBuffer, ptr ) ) then
						aver = aver + real( histogramVec(nExp).at( ptr ), 8 )
				end if
			end do
			aver = aver/real( numberOfExperiments, 8 )
			
			desv = 0.0_8
			do nExp=1,numberOfExperiments
				sBuffer = pair.first
				if( histogramVec(nExp).find( sBuffer, ptr ) ) then
						desv = desv + ( real( histogramVec(nExp).at( ptr ), 8 ) - aver )**2
				end if
			end do
			desv = sqrt( desv/real( numberOfExperiments, 8 ) )
			
! 			stdError = 1.96_8*desv/sqrt( real( numberOfExperiments, 8 ) )
			
			do nExp=1,numberOfExperiments
				sBuffer = pair.first
				if( histogramVec(nExp).find( sBuffer, ptr ) ) then
					if( nExp == 1 ) then
						write(effUnit,"(10X,A<maxNameLength>,F10.3)",advance="no") pair.first.fstr, normConst*histogramVec(nExp).at( ptr )
					else
						write(effUnit,"(F10.3)",advance="no") normConst*histogramVec(nExp).at( ptr )
					end if
				else
					if( nExp == 1 ) then
						write(effUnit,"(10X,A<maxNameLength>,F10.3)",advance="no") pair.first.fstr, 0
					else
						write(effUnit,"(F10.3)",advance="no") 0
					end if
				end if
			end do
			
! 			write(effUnit,"(5X,2F10.3)") normConst*aver, normConst*stdError
			write(effUnit,"(5X,2F10.3)") normConst*aver, normConst*desv
			
			if( present(ebklFileName) .and. present(excitationEnergy) ) then
				write(ebklFile.unit,"(A)") &
					trim(adjustl(FString_fromReal(excitationEnergy,"(F10.5)")))//"#" &
					//trim(pair.first.fstr)//"#" &
					//trim(adjustl(FString_fromReal(normConst*aver,"(F10.5)")))//":" &
! 					//trim(adjustl(FString_fromReal(normConst*stdError,"(F10.5)")))
					//trim(adjustl(FString_fromReal(normConst*desv,"(F10.5)")))
			end if
			
			iter => iter.next
		end do
		
		call averHistogram.clear()
		
		if( present(ebklFileName) .and. present(excitationEnergy) ) then
			call ebklFile.close()
		end if
	end subroutine showAverHistogram
	
	!>
	!! @brief Save the energy for the accepted steps in a file
	!!
	subroutine saveEnergyHistory( this, oFileName )
		class(MarkovChain), intent(in) :: this
		character(*), optional, intent(in) :: oFileName
		
		integer :: unit
		
		type(OFStream) :: oFile
		class(StringListIterator), pointer :: iter
		type(String) :: strBuffer
		integer :: i
		
		unit = 6
		if( present(oFileName) ) then
			call oFile.open( oFileName )
			unit = oFile.unit
		end if
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# ENERGY HISTORY"
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A1,3X,5A15,5X,A)") "#", "trans", "intermol", "vib", "rot", "tot", "formula"
		write(unit,"(A1,3X,5A15,5X,A)") "#", "eV", "eV", "eV", "eV", ""
		write(unit,"(A1,3X,5A15,5X,A)") "#", "-------", "--------", "-------", "-----", "-------"
		
		iter => this.energyHistory.begin
		do while( associated(iter) )
			
			strBuffer = this.energyHistory.at(iter)
			
			if( strBuffer.isEmpty() ) then
				write(unit,*) ""
			else
				write(unit,"(A)",advance="no") trim(strBuffer.fstr)
				write(unit,*)
			end if
			
			iter => iter.next
		end do
		
		write(unit,"(A)") ""
		
		if( present(oFileName) ) then
			call oFile.close()
		end if
	end subroutine saveEnergyHistory
	
	!>
	!! @brief Save the energy for the accepted steps in a file
	!!
	subroutine saveWeightHistory( this, oFileName )
		class(MarkovChain), intent(in) :: this
		character(*), optional, intent(in) :: oFileName
		
		integer :: unit
		
		type(OFStream) :: oFile
		class(StringListIterator), pointer :: iter
		type(String) :: strBuffer
		integer :: i
		
		unit = 6
		if( present(oFileName) ) then
			call oFile.open( oFileName )
			unit = oFile.unit
		end if
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# WEIGHT HISTORY"
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A1,3X,6A15,5X,A)") "#", "LnWe", "LnWv", "LnWn", "LnWr", "LnWt", "LnW", "formula"
		write(unit,"(A1,3X,6A15,5X,A)") "#", "arb.", "arb.", "arb.", "arb.", "arb.", "arb.", ""
		write(unit,"(A1,3X,6A15,5X,A)") "#", "-------", "-------", "--------", "--------", "--------", "-------", "-------"
		
		iter => this.weightHistory.begin
		do while( associated(iter) )
			
			strBuffer = this.weightHistory.at(iter)
			
			if( strBuffer.isEmpty() ) then
				write(unit,*) ""
			else
				write(unit,"(A)",advance="no") trim(strBuffer.fstr)
				write(unit,*)
			end if
			
			iter => iter.next
		end do
		
		write(unit,"(A)") ""
		
		if( present(oFileName) ) then
			call oFile.close()
		end if
	end subroutine saveWeightHistory
	
	!>
	!! @brief Save the energy for the accepted steps in a file
	!!
	subroutine saveJHistory( this, oFileName )
		class(MarkovChain), intent(in) :: this
		character(*), optional, intent(in) :: oFileName
		
		integer :: unit
		
		type(OFStream) :: oFile
		class(StringListIterator), pointer :: iter
		type(String) :: strBuffer
		
		unit = 6
		if( present(oFileName) ) then
			call oFile.open( oFileName )
			unit = oFile.unit
		end if
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# J HISTORY"
		write(unit,"(A)") "#------------------------------------"
		
		iter => this.JHistory.begin
		do while( associated(iter) )
			
			strBuffer = this.JHistory.at(iter)
			
			if( strBuffer.isEmpty() ) then
				write(unit,*) ""
			else
				write(unit,"(A)",advance="no") trim(strBuffer.fstr)
				write(unit,*)
			end if
			
			iter => iter.next
		end do
		
		write(unit,"(A)") ""
		
		if( present(oFileName) ) then
			call oFile.close()
		end if
	end subroutine saveJHistory
	
	!>
	!! @brief
	!!
	subroutine saveLHistory( this, oFileName )
		class(MarkovChain), intent(in) :: this
		character(*), optional, intent(in) :: oFileName
		
		integer :: unit
		
		type(OFStream) :: oFile
		class(StringListIterator), pointer :: iter
		type(String) :: strBuffer
		
		unit = 6
		if( present(oFileName) ) then
			call oFile.open( oFileName )
			unit = oFile.unit
		end if
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# L HISTORY"
		write(unit,"(A)") "#------------------------------------"
		
		iter => this.LHistory.begin
		do while( associated(iter) )
			
			strBuffer = this.LHistory.at(iter)
			
			if( strBuffer.isEmpty() ) then
				write(unit,*) ""
			else
				write(unit,"(A)",advance="no") trim(strBuffer.fstr)
				write(unit,*)
			end if
			
			iter => iter.next
		end do
		
		write(unit,"(A)") ""
		
		if( present(oFileName) ) then
			call oFile.close()
		end if
	end subroutine saveLHistory
	
	!>
	!! @brief
	!!
	subroutine saveHistograms( this, oFileName, genEbkl )
		class(MarkovChain), intent(in) :: this
		character(*), optional, intent(in) :: oFileName
		logical, optional, intent(in) :: genEbkl
		
		logical :: EffgenEbkl
		
		character(:), allocatable :: ebklFileName
		type(OFStream) :: oFile
		integer :: unit
		type(RealHistogram) :: histBuffer
		
		class(StringRealMapIterator), pointer :: simIter
		type(StringRealPair) :: srPair
		
		integer :: i
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Estos valores son temporales para el
		! cálculo de la temperatura
		class(RealListIterator), pointer :: iter
		real(8) :: temperature, stdTemperature
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		unit = 6
		if( present(oFileName) ) then
			call oFile.open( oFileName )
			unit = oFile.unit
		end if
		
		EffgenEbkl = .false.
		if( present(genEbkl) ) EffgenEbkl = genEbkl
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# Channels histogram"
		write(unit,"(A)") "#------------------------------------"
		
		if( EffgenEbkl ) then
			ebklFileName = "E_"//trim(adjustl(FString_fromReal(this.excitationEnergy/eV,"(F10.5)")))//".eblkC"
			call showAverHistogram( this.channelHistogram, unit, ebklFileName, this.excitationEnergy/eV )
		else
			call showAverHistogram( this.channelHistogram, unit )
		end if
		
		write(unit,*) ""
		write(unit,*) ""
		
		if( EffgenEbkl ) then
			ebklFileName = "E_"//trim(adjustl(FString_fromReal(this.excitationEnergy/eV,"(F10.5)")))//".eblkCd"
			call showAverHistogram( this.channelDetHistogram, unit, ebklFileName, this.excitationEnergy/eV )
		else
			call showAverHistogram( this.channelDetHistogram, unit )
		end if
		
		write(unit,*) ""
		write(unit,*) ""
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# Species histogram"
		write(unit,"(A)") "#------------------------------------"
		
		if( EffgenEbkl ) then
			ebklFileName = "E_"//trim(adjustl(FString_fromReal(this.excitationEnergy/eV,"(F10.5)")))//".eblkS"
			call showAverHistogram( this.speciesHistogram, unit, ebklFileName, this.excitationEnergy/eV )
		else
			call showAverHistogram( this.speciesHistogram, unit )
		end if
		
		write(unit,*) ""
		write(unit,*) ""
		
		if( EffgenEbkl ) then
			ebklFileName = "E_"//trim(adjustl(FString_fromReal(this.excitationEnergy/eV,"(F10.5)")))//".eblkSd"
			call showAverHistogram( this.speciesDetHistogram, unit, ebklFileName, this.excitationEnergy/eV )
		else
			call showAverHistogram( this.speciesDetHistogram, unit )
		end if
		
		write(unit,*) ""
		write(unit,*) ""
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# Temperature (eV)"
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A1,9X,<this.numberOfExperiments>I10,5X,2A10)") "#", ( i, i=1,this.numberOfExperiments ), "aver", "desv"
		write(unit,"(A1,9X,<this.numberOfExperiments>A10,5X,2A10)") "#", ( "-----", i=1,this.numberOfExperiments ), "----", "----"
		
		call histBuffer.initRealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( 1.0_8/this.iTemperatureHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,F10.3)",advance="no") 1.0_8/this.iTemperatureHistogram(i).mean()/eV
			else
				write(unit,"(F10.3)",advance="no") 1.0_8/this.iTemperatureHistogram(i).mean()/eV
			end if
		end do
		
		write(unit,"(5X,2F10.3)") histBuffer.mean(), histBuffer.stdev()
		
		write(unit,*) ""
		write(unit,*) ""
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# Entropy"
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A1,9X,<this.numberOfExperiments>I10,5X,2A10)") "#", ( i, i=1,this.numberOfExperiments ), "aver", "desv"
		write(unit,"(A1,9X,<this.numberOfExperiments>A10,5X,2A10)") "#", ( "-----", i=1,this.numberOfExperiments ), "----", "----"
		
		call histBuffer.initRealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( this.entropyHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,F10.3)",advance="no") this.entropyHistogram(i).mean()
			else
				write(unit,"(F10.3)",advance="no") this.entropyHistogram(i).mean()
			end if
		end do
		
		write(unit,"(5X,2F10.3)") histBuffer.mean(), histBuffer.stdev()
		
		write(unit,*) ""
		write(unit,*) ""
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# Markov chain statistics"
		write(unit,"(A)") "#------------------------------------"
		
		write(unit,"(A1,A)") "#", " Reactor type (ACCEPTED)"
		call this.reactorAcceptedHistogram.build()
		call this.reactorAcceptedHistogram.densityBegin( simIter )
		do while( associated(simIter) )
			srPair = this.reactorAcceptedHistogram.pair( simIter )
			write(unit,"(A20,F15.5)") srPair.first.fstr, srPair.second
			
			simIter => simIter.next
		end do
		
		write(unit,*) ""
		write(unit,"(A1,A)") "#", " Reactor type (REJECTED)"
		call this.reactorRejectedHistogram.build()
		call this.reactorRejectedHistogram.densityBegin( simIter )
		do while( associated(simIter) )
			srPair = this.reactorRejectedHistogram.pair( simIter )
			write(unit,"(A20,F15.5)") srPair.first.fstr, srPair.second
			
			simIter => simIter.next
		end do
		
		write(unit,*) ""
		write(unit,"(A1,A)") "#", " Reactor status"
		call this.reactorStatusHistogram.build()
		call this.reactorStatusHistogram.densityBegin( simIter )
		do while( associated(simIter) )
			srPair = this.reactorStatusHistogram.pair( simIter )
			write(unit,"(A20,F15.5)") srPair.first.fstr, srPair.second
			
			simIter => simIter.next
		end do
		
		write(unit,"(5X)")
		
		write(unit,*) ""
		write(unit,*) ""
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# Transition statistics"
		write(unit,"(A)") "#------------------------------------"
		
		write(unit,"(A1,A49,A15)") "#", "reaction", "prob."
		write(unit,"(A1,A49,A15)") "#", "--------", "-----"
		call this.transitionHistogram.build()
		call this.transitionHistogram.densityBegin( simIter )
		do while( associated(simIter) )
			srPair = this.transitionHistogram.pair( simIter )
			write(unit,"(A50,F15.5)") srPair.first.fstr, srPair.second
			
			simIter => simIter.next
		end do
		
		write(unit,*) ""
		write(unit,*) ""
		
		write(unit,"(A1,A49,A15)") "#", "reaction", "prob."
		write(unit,"(A1,A49,A15)") "#", "--------", "-----"
		call this.transitionDetHistogram.build()
		call this.transitionDetHistogram.densityBegin( simIter )
		do while( associated(simIter) )
			srPair = this.transitionDetHistogram.pair( simIter )
			write(unit,"(A50,F15.5)") srPair.first.fstr, srPair.second
			
			simIter => simIter.next
		end do
		
		write(unit,*) ""
		write(unit,*) ""
		
		if( present(oFileName) ) then
			call oFile.close()
		end if
	end subroutine saveHistograms
	
	!>
	!! @brief
	!!
	subroutine execute( this, iParser )
		class(MarkovChain) :: this
		type(BlocksIFileParser), intent(in) :: iParser
		
		type(FragmentsList) :: reactives
		character(20), allocatable :: reactiveTokens(:)
		type(String) :: energyHistoryFile
		type(String) :: weightHistoryFile
		type(String) :: JHistoryFile
		type(String) :: LHistoryFile
		type(String) :: histogramFile
		logical :: genEbkl
		
		type(String) :: strReactives
		type(String) :: sBuffer
		integer :: i, iBuffer
		
! 		if( size(FragmentsDB_instance.clusters) >= 40 ) then
! 			write(*,"(A)") "### ERROR ###: MarkovChain.execute(). The fragments database limit of 40 have been exceeded."
! 			write(*,"(A)") "               Please contact us:"
! 			write(*,"(A)") "                      Dr. Nestor F. Aguirre ( nestor.aguirre@uam.es )"
! 			write(*,"(A)") "                      Dr. Sergio Díaz-Tendero ( sergio.diaztendero@uam.es )"
! 			write(*,"(A)") "                      Prof. M. Paul-Antoine Hervieux ( Paul-Antoine.Hervieux@ipcms.unistra.fr )"
! 			write(*,"(A)") "                      Prof. Fernando Martín ( fernando.martin@uam.es )"
! 			write(*,"(A)") "                      Prof. Manuel Alcamí ( manuel.alcami@uam.es )"
! 			stop
! 		end if
		
		sBuffer = iParser.getString( "MARKOV_CHAIN:reactives" )
		
		call sBuffer.split( reactiveTokens, ":" )
		
		if( reactiveTokens(1) == "file" ) then
			call reactives.loadXYZ( reactiveTokens(2) )
		else
			strReactives = FragmentsDB_instance.extendFragmentsListName( sBuffer.fstr )
			call strReactives.split( reactiveTokens, "+" )
			
			call reactives.init( size(reactiveTokens) )
			do i=1,size(reactiveTokens)
				iBuffer = FragmentsDB_instance.getIdFromName( reactiveTokens(i) )
				call reactives.set( i, FragmentsDB_instance.clusters(iBuffer) )
			end do
		end if
		
		call this.init()
		
		this.burnInFraction = iParser.getReal( "MARKOV_CHAIN:burnInFraction", def=0.1_8 )
		this.excitationEnergy = iParser.getReal( "MARKOV_CHAIN:excitationEnergy", def=10.0_8 )*eV
		this.numberOfEvents = iParser.getInteger( "MARKOV_CHAIN:numberOfEvents", def=1000 )
		this.numberOfExperiments = iParser.getInteger( "MARKOV_CHAIN:numberOfExperiments", def=5 )
		iBuffer = max(10,this.numberOfEvents/1000)
		this.historyFileFrequency = iParser.getInteger( "MARKOV_CHAIN:historyFileFrequency", def=iBuffer )
		this.task = iParser.getString( "MARKOV_CHAIN:task", def="V" )
		this.geometryHistoryFilePrefix = iParser.getString( "MARKOV_CHAIN:geometryHistoryFilePrefix", def="" )
		this.freqBlockingCheck = iParser.getInteger( "MARKOV_CHAIN:freqBlockingCheck", def=4 )
		this.tracking = iParser.getString( "MARKOV_CHAIN:tracking", def="none" )
		
		write(*,*)
		write(*,"(A40,A)") "reactives = ", strReactives.fstr
		write(*,"(A40,F15.5)") "burnInFraction = ", this.burnInFraction
		write(*,"(A40,F15.5,A)") "excitationEnergy = ", this.excitationEnergy/eV, " eV"
		write(*,"(A40,I15)") "numberOfEvents = ", this.numberOfEvents
		write(*,"(A40,I15)") "numberOfExperiments = ", this.numberOfExperiments
		write(*,"(A40,A)") "task = ", trim(this.task.fstr)
		write(*,"(A40,A)") "geometryHistoryFilePrefix = ", trim(this.geometryHistoryFilePrefix.fstr)
		write(*,"(A40,I15)") "freqBlockingCheck = ", this.freqBlockingCheck
		write(*,"(A40,A)") "tracking = ", trim(this.tracking.fstr)
		write(*,*)
		
		call this.run( reactives )
		
		energyHistoryFile = iParser.getString( "MARKOV_CHAIN:energyHistoryFile", def="#@NONE@#" )
		if( trim(energyHistoryFile.fstr) /= "#@NONE@#" ) then
			call this.saveEnergyHistory( energyHistoryFile.fstr )
		end if
		
		weightHistoryFile = iParser.getString( "MARKOV_CHAIN:weightHistoryFile", def="#@NONE@#" )
		if( trim(weightHistoryFile.fstr) /= "#@NONE@#" ) then
			call this.saveWeightHistory( weightHistoryFile.fstr )
		end if
		
		JHistoryFile = iParser.getString( "MARKOV_CHAIN:JHistoryFile", def="#@NONE@#" )
		if( trim(JHistoryFile.fstr) /= "#@NONE@#" ) then
			call this.saveJHistory( JHistoryFile.fstr )
		end if
		
		LHistoryFile = iParser.getString( "MARKOV_CHAIN:LHistoryFile", def="#@NONE@#" )
		if( trim(LHistoryFile.fstr) /= "#@NONE@#" ) then
			call this.saveLHistory( LHistoryFile.fstr )
		end if
		
		histogramFile = iParser.getString( "MARKOV_CHAIN:histogramFile", def="#@NONE@#" )
		genEbkl = iParser.getLogical( "MARKOV_CHAIN:genEbkl", def=.false. )
		if( trim(histogramFile.fstr) /= "#@NONE@#" ) then
			call this.saveHistograms( histogramFile.fstr, genEbkl )
		end if
		
	end subroutine execute
end module MarkovChain_

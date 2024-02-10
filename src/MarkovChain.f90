!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                   !!
!! This file is part of M3C project                                                  !!
!!                                                                                   !!
!! Copyright (c) 2019-2020 by authors                                                !!
!! Authors:                                                                          !!
!!                         * Néstor F. Aguirre (2019-2020)                           !!
!!                           nfaguirrec@gmail.com                                    !!
!!                                                                                   !!
!! Copyright (c) 2013-2016 Departamento de Química                                   !!
!!                         Universidad Autónoma de Madrid                            !!
!!                         All rights reserved.                                      !!
!!                                                                                   !!
!!                         * Néstor F. Aguirre (2013-2016)                           !!
!!                           nestor.aguirre@uam.es                                   !!
!!                         * Sergio Díaz-Tendero (2013-2016)                         !!
!!                           sergio.diaztendero@uam.es                               !!
!!                         * M. Paul-Antoine Hervieux (2013-2015)                    !!
!!                           Paul-Antoine.Hervieux@ipcms.unistra.fr                  !!
!!                         * Manuel Alcamí (2013-2016)                               !!
!!                           manuel.alcami@uam.es                                    !!
!!                         * Fernando Martín (2013-2016)                             !!
!!                           fernando.martin@uam.es                                  !!
!!                                                                                   !!
!!  Redistribution and use in source and binary forms, with or without               !!
!!  modification, are permitted provided that the following conditions are met:      !!
!!                                                                                   !!
!!  1. Redistributions of source code must retain the above copyright notice, this   !!
!!     list of conditions and the following disclaimer.                              !!
!!  2. Redistributions in binary form must reproduce the above copyright notice,     !!
!!     this list of conditions and the following disclaimer in the documentation     !!
!!     and/or other materials provided with the distribution.                        !!
!!  3. Neither the name of the copyright holders nor the names of its contributors   !!
!!     may be used to endorse or promote products derived from this software         !!
!!     without specific prior written permission.                                    !!
!!                                                                                   !!
!!  The copyright holders provide no reassurances that the source code provided      !!
!!  does not infringe any patent, copyright, or any other intellectual property      !!
!!  rights of third parties.  The copyright holders disclaim any liability to any    !!
!!  recipient for claims brought against recipient by any third party for            !!
!!  infringement of that parties intellectual property rights.                       !!
!!                                                                                   !!
!!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  !!
!!  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    !!
!!  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           !!
!!  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR  !!
!!  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES   !!
!!  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     !!
!!  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      !!
!!  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       !!
!!  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    !!
!!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     !!
!!                                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! @brief
!!
module MarkovChain_
	use String_
	use StringList_
	use RealList_
	use RandomUtils_
	use UnitsConverter_
	use IOStream_
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
	
	use GOptionsM3C_
	
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
		
		type(OFStream) :: energyHistoryFile
		type(OFStream) :: weightHistoryFile
		type(OFStream) :: JHistoryFile
		type(OFStream) :: LHistoryFile
		
		type(RealHistogram), allocatable :: iTemperatureHistogram(:) ! One item for each experiment
		
		type(RealHistogram), allocatable :: translationalEnergyHistogram(:)   ! One item for each experiment
		type(RealHistogram), allocatable :: intermolEnergyHistogram(:)        ! One item for each experiment
		type(RealHistogram), allocatable :: vibrationalEnergyHistogram(:)     ! One item for each experiment
		type(RealHistogram), allocatable :: rotationalEnergyHistogram(:)      ! One item for each experiment
		
		type(RealHistogram), allocatable :: electronicWeightHistogram(:)      ! One item for each experiment
		type(RealHistogram), allocatable :: vibrationalWeightHistogram(:)     ! One item for each experiment
		type(RealHistogram), allocatable :: combinatorialWeightHistogram(:)   ! One item for each experiment
		type(RealHistogram), allocatable :: rotationalWeightHistogram(:)      ! One item for each experiment
		type(RealHistogram), allocatable :: translationalWeightHistogram(:)   ! One item for each experiment
		type(RealHistogram), allocatable :: totalWeightHistogram(:)           ! One item for each experiment
		
		type(StringHistogram) :: reactorAcceptedHistogram    ! One for all experiments
		type(StringHistogram) :: reactorRejectedHistogram    ! One for all experiments
		type(StringHistogram) :: reactorStatusHistogram      ! One for all experiments
		
		type(StringHistogram) :: transitionHistogram         ! One for all experiments. Probability to change from one channel to another
		type(StringHistogram) :: transitionDetHistogram      ! One for all experiments. Probability to change from one channel to another
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Histogramas
		type(StringIntegerMap), allocatable :: nFragsHistogram(:) ! One item for each experiment
		
		type(StringIntegerMap), allocatable :: speciesHistogram(:) ! One item for each experiment
		type(StringIntegerMap), allocatable :: speciesDetHistogram(:) ! One item for each experiment
		
		type(StringIntegerMap), allocatable :: channelHistogram(:) ! One item for each experiment
		type(StringIntegerMap), allocatable :: channelDetHistogram(:) ! One item for each experiment
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Strings de history
! 		type(StringList), private :: energyHistory
! 		type(StringList), private :: weightHistory
! 		type(StringList), private :: JHistory
! 		type(StringList), private :: LHistory
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
! 			procedure :: saveEnergyHistory
! 			procedure :: saveWeightHistory
! 			procedure :: saveJHistory
! 			procedure :: saveLHistory
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
		
		if( allocated(this.iTemperatureHistogram) ) deallocate( this.iTemperatureHistogram )
		
		if( allocated(this.translationalEnergyHistogram) ) deallocate( this.translationalEnergyHistogram )
		if( allocated(this.intermolEnergyHistogram) ) deallocate( this.intermolEnergyHistogram )
		if( allocated(this.vibrationalEnergyHistogram) ) deallocate( this.vibrationalEnergyHistogram )
		if( allocated(this.rotationalEnergyHistogram) ) deallocate( this.rotationalEnergyHistogram )
		
		if( allocated(this.electronicWeightHistogram) ) deallocate( this.electronicWeightHistogram )
		if( allocated(this.vibrationalWeightHistogram) ) deallocate( this.vibrationalWeightHistogram )
		if( allocated(this.combinatorialWeightHistogram) ) deallocate( this.combinatorialWeightHistogram )
		if( allocated(this.rotationalWeightHistogram) ) deallocate( this.rotationalWeightHistogram )
		if( allocated(this.translationalWeightHistogram) ) deallocate( this.translationalWeightHistogram )
		if( allocated(this.totalWeightHistogram) ) deallocate( this.totalWeightHistogram )
		
		call this.reactorAcceptedHistogram.clear()
		call this.reactorRejectedHistogram.clear()
		call this.reactorStatusHistogram.clear()
		
		call this.transitionHistogram.clear()
		call this.transitionDetHistogram.clear()
		
		if( allocated(this.nFragsHistogram) ) then
			do nExp=1,this.numberOfExperiments
				call this.nFragsHistogram(nExp).clear()
			end do
			
			deallocate( this.nFragsHistogram )
		end if
		
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
		
		allocate( this.translationalEnergyHistogram(this.numberOfExperiments) )
		allocate( this.intermolEnergyHistogram(this.numberOfExperiments) )
		allocate( this.vibrationalEnergyHistogram(this.numberOfExperiments) )
		allocate( this.rotationalEnergyHistogram(this.numberOfExperiments) )
		
		allocate( this.electronicWeightHistogram(this.numberOfExperiments) )
		allocate( this.vibrationalWeightHistogram(this.numberOfExperiments) )
		allocate( this.combinatorialWeightHistogram(this.numberOfExperiments) )
		allocate( this.rotationalWeightHistogram(this.numberOfExperiments) )
		allocate( this.translationalWeightHistogram(this.numberOfExperiments) )
		allocate( this.totalWeightHistogram(this.numberOfExperiments) )
		
		allocate( this.nFragsHistogram(this.numberOfExperiments) )
		allocate( this.speciesHistogram(this.numberOfExperiments) )
		allocate( this.speciesDetHistogram(this.numberOfExperiments) )
		allocate( this.channelHistogram(this.numberOfExperiments) )
		allocate( this.channelDetHistogram(this.numberOfExperiments) )
		
		call this.reactorAcceptedHistogram.initStringHistogram( algorithm=Histogram_RUNNING )
		call this.reactorRejectedHistogram.initStringHistogram( algorithm=Histogram_RUNNING )
		call this.reactorStatusHistogram.initStringHistogram( algorithm=Histogram_RUNNING )
		
		call this.transitionHistogram.initStringHistogram( algorithm=Histogram_RUNNING )
		call this.transitionDetHistogram.initStringHistogram( algorithm=Histogram_RUNNING )
		
		do nExp=1,this.numberOfExperiments
			this.iTemperatureHistogram(nExp) = RealHistogram( algorithm=Histogram_RUNNING )
			
			this.translationalEnergyHistogram(nExp) = RealHistogram( algorithm=Histogram_RUNNING )
			this.intermolEnergyHistogram(nExp) = RealHistogram( algorithm=Histogram_RUNNING )
			this.vibrationalEnergyHistogram(nExp) = RealHistogram( algorithm=Histogram_RUNNING )
			this.rotationalEnergyHistogram(nExp) = RealHistogram( algorithm=Histogram_RUNNING )
			
			this.electronicWeightHistogram(nExp) = RealHistogram( algorithm=Histogram_RUNNING )
			this.vibrationalWeightHistogram(nExp) = RealHistogram( algorithm=Histogram_RUNNING )
			this.combinatorialWeightHistogram(nExp) = RealHistogram( algorithm=Histogram_RUNNING )
			this.rotationalWeightHistogram(nExp) = RealHistogram( algorithm=Histogram_RUNNING )
			this.translationalWeightHistogram(nExp) = RealHistogram( algorithm=Histogram_RUNNING )
			this.totalWeightHistogram(nExp) = RealHistogram( algorithm=Histogram_RUNNING )
			
			this.nFragsHistogram(nExp) = StringIntegerMap()
			this.speciesHistogram(nExp) = StringIntegerMap()
			this.speciesDetHistogram(nExp) = StringIntegerMap()
			this.channelHistogram(nExp) = StringIntegerMap()
			this.channelDetHistogram(nExp) = StringIntegerMap()
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
		character(100), allocatable :: taskTokens(:)
		character(100), allocatable :: bufferTokens(:)
		character(100) :: currentTask
		integer :: taskMult
		real(8) :: p, Pi
		integer :: n, i, j, k, nExp, iBuffer
		type(String) :: sBuffer
		character(100) :: geometryFileName
		character(2) :: origin
		
		!!--------------------------------------------------------------------
		!! Este bloque expande los operadores n*(XXX) a XXX,XXX,XXX,...,XXX
		!! @todo El codigo funciona bien, pero no me gusta.
		!!--------------------------------------------------------------------
		character(100), allocatable :: taskTokens2(:)
		
		call this.task.split( taskTokens, "()" )
		
		do i=1,size(taskTokens)
			sBuffer = trim(taskTokens(i))
			call sBuffer.split( taskTokens2, ",*" )
			
			if( FString_isInteger( trim(taskTokens2(size(taskTokens2))) ) ) then
				taskMult = FString_toInteger( taskTokens2(size(taskTokens2)) )
				currentTask = trim(adjustl( taskTokens(i+1) ))
				
				sBuffer = ""
				do j=1,taskMult
					sBuffer = trim(sBuffer.fstr)//trim(currentTask)//","
				end do
				
				call this.task.replace( trim(FString_fromInteger(taskMult))//"*("//trim(currentTask)//")", trim(sBuffer.fstr) )
			end if
		end do
		!!--------------------------------------------------------------------
		
! 		call this.energyHistory.clear()
! 		call this.weightHistory.clear()
		
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
			write(6,"(A1,3X,5A15,5X,A)") "#", "trans", "intermol", "vib", "rot", "tot", "label"
			write(6,"(A1,3X,5A15,5X,A)") "#", "eV", "eV", "eV", "eV", ""
			write(6,"(A1,3X,5A15,5X,A)") "#", "-------", "--------", "-------", "-----", "-------"

		else if( this.tracking == "weight" ) then
		
			write(6,"(A)") "#------------------------------------"
			write(6,"(A)") "# WEIGHT HISTORY"
			write(6,"(A)") "#------------------------------------"
			write(6,"(A1,3X,6A15,5X,A)") "#", "LnWe", "LnWv", "LnWn", "LnWr", "LnWt", "LnW", "label"
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
					
					call FString_split( taskTokens(i), bufferTokens, "*" )
					
					if( FString_isNumeric( bufferTokens(1) ) ) then
						taskMult = FString_toInteger( bufferTokens(1) )
						currentTask = trim(adjustl( bufferTokens(2) ))
					else
						taskMult = 1
						currentTask = trim(adjustl( bufferTokens(1) ))
					end if
					
					if( allocated(bufferTokens) ) deallocate( bufferTokens )
					
! 					currentTask = trim(adjustl(taskTokens(i)))
					call react.setType( currentTask )
					
					do k=1,taskMult
						if( n == this.numberOfEvents+1 ) exit
						
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
								
								if( GOptionsM3C_useRandomWalkers ) then
									GOptionsM3C_randomWalkStepRadius = GOptionsM3C_randomWalkStepRadius*0.9_8
									call GOptions_info( &
										"Blocked. Reducing randomWalkStepRadius in 1% ("//FString_fromReal(GOptionsM3C_randomWalkStepRadius/angs)//" A)", &
										"MarkovChain" )
								else
									call GOptions_info( "Blocked. Forcing centers completely random", "MarkovChain" )
									react.reactives.forceRandomCenters = .true.
								end if
								
								nTimesBlocked = 0
							else
								nTimesBlocked = nTimesBlocked + 1
							end if
								
							call GOptions_info( "Step rejected ( Negative energy )", "MarkovChain" )
							
							origin = "e"//trim(currentTask)
							
							call this.reactorRejectedHistogram.add( FString_toString( trim(currentTask) ) )
							call this.reactorStatusHistogram.add( FString_toString( "e.REJECTED(E<0) " ) )
						else
						
							nTimesBlocked = 0  ! El bloqueo debe ser consecutivo, así que si no pasa por negative energy se cuenta nuevamente
							
! 							if( react.replaceTS ) then
! 								react.replaceTS = .false.
! 								write(*,*) "Hola = ", react.productsTS.LnW(), react.reactives.LnW()
! 								Pi = react.productsTS.LnW()-react.reactives.LnW()
! 							else
! 								Pi = react.products.LnW()-react.reactives.LnW()
! 							end if

							if( react.replaceTS ) then
								react.replaceTS = .false.
								
								select case( trim(GOptionsM3C_TSModel.fstr) )
									case( "NONE" )
										Pi = react.products.LnW()-react.reactives.LnW()  ! No correction by TS
									case( "EARLY" )
										Pi = react.productsTS.LnW()-react.reactives.LnW()
										
! 										write(*,"(A,4F10.5)") "W,WTS,Wreact,Pi =", react.products.LnW(), react.productsTS.LnW(), react.reactives.LnW(), Pi
									case( "LATE" )
										Pi = react.productsTS2.LnW()+react.TS.LnWv()-react.reactives.LnW()
										
! 										write(*,"(A,4F10.5)") "W,WTS,Wreact,Pi =", react.products.LnW(), react.productsTS2.LnW()+react.TS.LnWv(), react.reactives.LnW(), Pi
									case default
										write(*,"(A)") "### ERROR ###: MarkovChain.run(). Unknown TS model ("//trim(GOptionsM3C_TSModel.fstr)//")"
										write(*,"(A)") "               Available options: NONE, EARLY, LATE. Default NONE"
										stop
								end select
							else
								Pi = react.products.LnW()-react.reactives.LnW()  ! Normal channel without TS. No correction by TS needed
							end if
							
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
								
								call this.reactorAcceptedHistogram.add( FString_toString( trim(currentTask) ) )
								call this.reactorStatusHistogram.add( FString_toString( "a.ACCEPTED      " ) )
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
									
									call this.reactorAcceptedHistogram.add( FString_toString( trim(currentTask) ) )
									call this.reactorStatusHistogram.add( FString_toString( "p.ACCEPTED(p<PI)" ) )
								else
									call GOptions_info( &
									"Step rejected logP="//trim(adjustl(FString_fromReal(p,"(F10.3)")))// &
									",logPi="//trim(adjustl(FString_fromReal(Pi,"(F10.3)")))//"  ( logP > logPi )", "MarkovChain" )
									
									origin = "r"//trim(currentTask)
									
									call this.reactorRejectedHistogram.add( FString_toString( trim(currentTask) ) )
									call this.reactorStatusHistogram.add( FString_toString( "r.REJECTED      " ) )
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
								
								if( this.weightHistoryFile.isOpen() ) then
									sBuffer = react.reactives.weightHistoryLine( origin )
									write( this.weightHistoryFile.unit, "(A)" ) trim(sBuffer.fstr)
								end if
								
								if( this.energyHistoryFile.isOpen() ) then
									sBuffer = react.reactives.energyHistoryLine( origin )
									write( this.energyHistoryFile.unit, "(A)" ) trim(sBuffer.fstr)
								end if
								
								if( .not. this.geometryHistoryFilePrefix.isEmpty() ) then
									geometryFileName = trim(this.geometryHistoryFilePrefix.fstr)//"-"//trim(FString_fromInteger(nExp))//".xyz"
									call react.reactives.save(geometryFileName, append=.true.)
								end if
								
							end if
						
							call this.iTemperatureHistogram(nExp).add( react.reactives.iTemperature() )
							
							call this.translationalEnergyHistogram(nExp).add( react.reactives.translationalEnergy() )
							call this.intermolEnergyHistogram(nExp).add( react.reactives.intermolEnergy() )
							call this.vibrationalEnergyHistogram(nExp).add( react.reactives.vibrationalEnergy() )
							call this.rotationalEnergyHistogram(nExp).add( react.reactives.rotationalEnergy() )
							
							call this.electronicWeightHistogram(nExp).add( react.reactives.LnWe() )
							call this.vibrationalWeightHistogram(nExp).add( react.reactives.LnWv() )
							call this.combinatorialWeightHistogram(nExp).add( react.reactives.LnWn() )
							call this.rotationalWeightHistogram(nExp).add( 0.5_8*react.reactives.LnDiagI()+react.reactives.logVtheta_ )
							call this.translationalWeightHistogram(nExp).add( react.reactives.LnLambda()-0.5_8*react.reactives.LnDiagI()-react.reactives.logVtheta_ )
							call this.totalWeightHistogram(nExp).add( react.reactives.LnW() )
							
							sBuffer = trim(FString_fromInteger( react.reactives.nMolecules() ))
							call this.nFragsHistogram(nExp).set( sBuffer, this.nFragsHistogram(nExp).at( sBuffer, defaultValue=0 )+1 )
							
							do j=1,react.reactives.nMolecules()
								sBuffer = react.reactives.clusters(j).label( details=.false. )
								call this.speciesHistogram(nExp).set( sBuffer, this.speciesHistogram(nExp).at( sBuffer, defaultValue=0 )+1 )
								
								sBuffer = react.reactives.clusters(j).label( details=.true. )
								call this.speciesDetHistogram(nExp).set( sBuffer, this.speciesDetHistogram(nExp).at( sBuffer, defaultValue=0 )+1 )
							end do
							
							if( this.JHistoryFile.isOpen() ) then
								sBuffer = react.reactives.JHistoryLine()
								write( this.JHistoryFile.unit, "(A)" ) trim(sBuffer.fstr)
							end if
							
							if( this.LHistoryFile.isOpen() ) then
								sBuffer = react.reactives.LHistoryLine()
								write( this.LHistoryFile.unit, "(A)" ) trim(sBuffer.fstr)
							end if
							
							sBuffer = react.reactives.label( details=.false. )
							call this.channelHistogram(nExp).set( sBuffer, this.channelHistogram(nExp).at( sBuffer, defaultValue=0 )+1 )
							
							sBuffer = react.reactives.label( details=.true. )
							call this.channelDetHistogram(nExp).set( sBuffer, this.channelDetHistogram(nExp).at( sBuffer, defaultValue=0 )+1 )
							
						end if
						
						n = n+1
					end do
				end do
			end do
			
			sBuffer = ""
! 			if( this.tracking == "energy" ) then
				
				if( this.weightHistoryFile.isOpen() ) then
					write( this.weightHistoryFile.unit, "(A)" ) trim(sBuffer.fstr)
					write( this.weightHistoryFile.unit, "(A)" ) trim(sBuffer.fstr)
				end if
				
! 			else if( this.tracking == "weight" ) then
			
				if( this.energyHistoryFile.isOpen() ) then
					write( this.energyHistoryFile.unit, "(A)" ) trim(sBuffer.fstr)
					write( this.energyHistoryFile.unit, "(A)" ) trim(sBuffer.fstr)
				end if
				
! 			end if

			if( this.tracking /= "none" ) then
				write(6,"(A)") ""
				write(6,"(A)") ""
			end if
			
			if( this.JHistoryFile.isOpen() ) then
				write( this.JHistoryFile.unit, "(A)" ) trim(sBuffer.fstr)
				write( this.JHistoryFile.unit, "(A)" ) trim(sBuffer.fstr)
			end if
			
			if( this.LHistoryFile.isOpen() ) then
				write( this.LHistoryFile.unit, "(A)" ) trim(sBuffer.fstr)
				write( this.LHistoryFile.unit, "(A)" ) trim(sBuffer.fstr)
			end if
		end do
		
		if( GOptions_printLevel >= 1 ) then
			write(*,"(A)") ""
		end if
		
! 		if( this.tracking == "energy" ) then
! 			call this.saveWeightHistory()
! 		else if( this.tracking == "weight" ) then
! 			call this.saveEnergyHistory()
! 		end if
			
		call this.saveHistograms()
		
		if( allocated(taskTokens) ) deallocate( taskTokens )
		if( allocated(bufferTokens) ) deallocate( bufferTokens )
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
		
		averHistogram = StringIntegerMap()
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
	
! 	!>
! 	!! @brief Save the energy for the accepted steps in a file
! 	!!
! 	subroutine saveEnergyHistory( this, oFileName )
! 		class(MarkovChain), intent(in) :: this
! 		character(*), optional, intent(in) :: oFileName
! 		
! 		integer :: unit
! 		
! 		type(OFStream) :: oFile
! 		class(StringListIterator), pointer :: iter
! 		type(String) :: strBuffer
! 		integer :: i
! 		
! 		unit = 6
! 		if( present(oFileName) ) then
! 			call oFile.open( oFileName )
! 			unit = oFile.unit
! 		end if
! 		
! 		write(unit,"(A)") "#------------------------------------"
! 		write(unit,"(A)") "# ENERGY HISTORY"
! 		write(unit,"(A)") "#------------------------------------"
! 		write(unit,"(A1,3X,5A15,5X,A)") "#", "trans", "intermol", "vib", "rot", "tot", "label"
! 		write(unit,"(A1,3X,5A15,5X,A)") "#", "eV", "eV", "eV", "eV", ""
! 		write(unit,"(A1,3X,5A15,5X,A)") "#", "-------", "--------", "-------", "-----", "-------"
! 		
! 		iter => this.energyHistory.begin
! 		do while( associated(iter) )
! 			
! 			strBuffer = this.energyHistory.at(iter)
! 			
! 			if( strBuffer.isEmpty() ) then
! 				write(unit,*) ""
! 			else
! 				write(unit,"(A)",advance="no") trim(strBuffer.fstr)
! 				write(unit,*)
! 			end if
! 			
! 			iter => iter.next
! 		end do
! 		
! 		write(unit,"(A)") ""
! 		
! 		if( present(oFileName) ) then
! 			call oFile.close()
! 		end if
! 	end subroutine saveEnergyHistory
	
! 	!>
! 	!! @brief Save the energy for the accepted steps in a file
! 	!!
! 	subroutine saveWeightHistory( this, oFileName )
! 		class(MarkovChain), intent(in) :: this
! 		character(*), optional, intent(in) :: oFileName
! 		
! 		integer :: unit
! 		
! 		type(OFStream) :: oFile
! 		class(StringListIterator), pointer :: iter
! 		type(String) :: strBuffer
! 		integer :: i
! 		
! 		unit = 6
! 		if( present(oFileName) ) then
! 			call oFile.open( oFileName )
! 			unit = oFile.unit
! 		end if
! 		
! 		write(unit,"(A)") "#------------------------------------"
! 		write(unit,"(A)") "# WEIGHT HISTORY"
! 		write(unit,"(A)") "#------------------------------------"
! 		write(unit,"(A1,3X,6A15,5X,A)") "#", "LnWe", "LnWv", "LnWn", "LnWr", "LnWt", "LnW", "label"
! 		write(unit,"(A1,3X,6A15,5X,A)") "#", "arb.", "arb.", "arb.", "arb.", "arb.", "arb.", ""
! 		write(unit,"(A1,3X,6A15,5X,A)") "#", "-------", "-------", "--------", "--------", "--------", "-------", "-------"
! 		
! 		iter => this.weightHistory.begin
! 		do while( associated(iter) )
! 			
! 			strBuffer = this.weightHistory.at(iter)
! 			
! 			if( strBuffer.isEmpty() ) then
! 				write(unit,*) ""
! 			else
! 				write(unit,"(A)",advance="no") trim(strBuffer.fstr)
! 				write(unit,*)
! 			end if
! 			
! 			iter => iter.next
! 		end do
! 		
! 		write(unit,"(A)") ""
! 		
! 		if( present(oFileName) ) then
! 			call oFile.close()
! 		end if
! 	end subroutine saveWeightHistory
	
! 	!>
! 	!! @brief Save the energy for the accepted steps in a file
! 	!!
! 	subroutine saveJHistory( this, oFileName )
! 		class(MarkovChain), intent(in) :: this
! 		character(*), optional, intent(in) :: oFileName
! 		
! 		integer :: unit
! 		
! 		type(OFStream) :: oFile
! 		class(StringListIterator), pointer :: iter
! 		type(String) :: strBuffer
! 		
! 		unit = 6
! 		if( present(oFileName) ) then
! 			call oFile.open( oFileName )
! 			unit = oFile.unit
! 		end if
! 		
! 		write(unit,"(A)") "#------------------------------------"
! 		write(unit,"(A)") "# J HISTORY"
! 		write(unit,"(A)") "#------------------------------------"
! 		
! 		iter => this.JHistory.begin
! 		do while( associated(iter) )
! 			
! 			strBuffer = this.JHistory.at(iter)
! 			
! 			if( strBuffer.isEmpty() ) then
! 				write(unit,*) ""
! 			else
! 				write(unit,"(A)",advance="no") trim(strBuffer.fstr)
! 				write(unit,*)
! 			end if
! 			
! 			iter => iter.next
! 		end do
! 		
! 		write(unit,"(A)") ""
! 		
! 		if( present(oFileName) ) then
! 			call oFile.close()
! 		end if
! 	end subroutine saveJHistory
	
! 	!>
! 	!! @brief
! 	!!
! 	subroutine saveLHistory( this, oFileName )
! 		class(MarkovChain), intent(in) :: this
! 		character(*), optional, intent(in) :: oFileName
! 		
! 		integer :: unit
! 		
! 		type(OFStream) :: oFile
! 		class(StringListIterator), pointer :: iter
! 		type(String) :: strBuffer
! 		
! 		unit = 6
! 		if( present(oFileName) ) then
! 			call oFile.open( oFileName )
! 			unit = oFile.unit
! 		end if
! 		
! 		write(unit,"(A)") "#------------------------------------"
! 		write(unit,"(A)") "# L HISTORY"
! 		write(unit,"(A)") "#------------------------------------"
! 		
! 		iter => this.LHistory.begin
! 		do while( associated(iter) )
! 			
! 			strBuffer = this.LHistory.at(iter)
! 			
! 			if( strBuffer.isEmpty() ) then
! 				write(unit,*) ""
! 			else
! 				write(unit,"(A)",advance="no") trim(strBuffer.fstr)
! 				write(unit,*)
! 			end if
! 			
! 			iter => iter.next
! 		end do
! 		
! 		write(unit,"(A)") ""
! 		
! 		if( present(oFileName) ) then
! 			call oFile.close()
! 		end if
! 	end subroutine saveLHistory
	
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
		write(unit,"(A)") "# NFragments histogram"
		write(unit,"(A)") "#------------------------------------"
		
		if( EffgenEbkl ) then
			ebklFileName = "E_"//trim(adjustl(FString_fromReal(this.excitationEnergy/eV,"(F10.5)")))//".eblkN"
			call showAverHistogram( this.nFragsHistogram, unit, ebklFileName, this.excitationEnergy/eV )
		else
			call showAverHistogram( this.nFragsHistogram, unit )
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
		write(unit,"(A)") "# Temperature (eV)"
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A1,9X,<this.numberOfExperiments>I15,5X,2A15)") "#", ( i, i=1,this.numberOfExperiments ), "aver", "desv"
		write(unit,"(A1,9X,<this.numberOfExperiments>A15,5X,2A15)") "#", ( "-----", i=1,this.numberOfExperiments ), "----", "----"
		
		histBuffer = RealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( 1.0_8/this.iTemperatureHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,F15.5)",advance="no") 1.0_8/this.iTemperatureHistogram(i).mean()/eV
			else
				write(unit,"(F15.5)",advance="no") 1.0_8/this.iTemperatureHistogram(i).mean()/eV
			end if
		end do
		
		write(unit,"(5X,2F15.5)") histBuffer.mean()/eV, histBuffer.stdev()/eV
		
		write(unit,*) ""
		write(unit,*) ""
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# Energy components (eV)"
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A1,9X,A5,<this.numberOfExperiments>I15,5X,2A15)") "#", "  ", ( i, i=1,this.numberOfExperiments ), "aver", "desv"
		write(unit,"(A1,9X,A5,<this.numberOfExperiments>A15,5X,2A15)") "#", "  ", ( "-----", i=1,this.numberOfExperiments ), "----", "----"
		
		histBuffer = RealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( this.translationalEnergyHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,A5,F15.5)",advance="no") "trans", this.translationalEnergyHistogram(i).mean()/eV
			else
				write(unit,"(F15.5)",advance="no") this.translationalEnergyHistogram(i).mean()/eV
			end if
		end do
		
		write(unit,"(5X,2F15.5)") histBuffer.mean()/eV, histBuffer.stdev()/eV
		
		histBuffer = RealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( this.intermolEnergyHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,A5,F15.5)",advance="no") "intermol", this.intermolEnergyHistogram(i).mean()/eV
			else
				write(unit,"(F15.5)",advance="no") this.intermolEnergyHistogram(i).mean()/eV
			end if
		end do
		
		write(unit,"(5X,2F15.5)") histBuffer.mean()/eV, histBuffer.stdev()/eV
		
		histBuffer = RealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( this.vibrationalEnergyHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,A5,F15.5)",advance="no") "vib", this.vibrationalEnergyHistogram(i).mean()/eV
			else
				write(unit,"(F15.5)",advance="no") this.vibrationalEnergyHistogram(i).mean()/eV
			end if
		end do
		
		write(unit,"(5X,2F15.5)") histBuffer.mean()/eV, histBuffer.stdev()/eV
		
		histBuffer = RealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( this.rotationalEnergyHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,A5,F15.5)",advance="no") "rot", this.rotationalEnergyHistogram(i).mean()/eV
			else
				write(unit,"(F15.5)",advance="no") this.rotationalEnergyHistogram(i).mean()/eV
			end if
		end do
		
		write(unit,"(5X,2F15.5)") histBuffer.mean()/eV, histBuffer.stdev()/eV
		
		write(unit,*) ""
		write(unit,*) ""
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# Weight components (arb.)"
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A1,9X,A5,<this.numberOfExperiments>I15,5X,2A15)") "#", "  ", ( i, i=1,this.numberOfExperiments ), "aver", "desv"
		write(unit,"(A1,9X,A5,<this.numberOfExperiments>A15,5X,2A15)") "#", "  ", ( "-----", i=1,this.numberOfExperiments ), "----", "----"
		
		histBuffer = RealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( this.electronicWeightHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,A5,F15.5)",advance="no") "LnWe", this.electronicWeightHistogram(i).mean()
			else
				write(unit,"(F15.5)",advance="no") this.electronicWeightHistogram(i).mean()
			end if
		end do
		
		write(unit,"(5X,2F15.5)") histBuffer.mean(), histBuffer.stdev()
		
		histBuffer = RealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( this.vibrationalWeightHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,A5,F15.5)",advance="no") "LnWv", this.vibrationalWeightHistogram(i).mean()
			else
				write(unit,"(F15.5)",advance="no") this.vibrationalWeightHistogram(i).mean()
			end if
		end do
		
		write(unit,"(5X,2F15.5)") histBuffer.mean(), histBuffer.stdev()
		
		histBuffer = RealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( this.combinatorialWeightHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,A5,F15.5)",advance="no") "LnWn", this.combinatorialWeightHistogram(i).mean()
			else
				write(unit,"(F15.5)",advance="no") this.combinatorialWeightHistogram(i).mean()
			end if
		end do
		
		write(unit,"(5X,2F15.5)") histBuffer.mean(), histBuffer.stdev()
		
		histBuffer = RealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( this.rotationalWeightHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,A5,F15.5)",advance="no") "LnWr", this.rotationalWeightHistogram(i).mean()
			else
				write(unit,"(F15.5)",advance="no") this.rotationalWeightHistogram(i).mean()
			end if
		end do
		
		write(unit,"(5X,2F15.5)") histBuffer.mean(), histBuffer.stdev()
		
		histBuffer = RealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( this.translationalWeightHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,A5,F15.5)",advance="no") "LnWt", this.translationalWeightHistogram(i).mean()
			else
				write(unit,"(F15.5)",advance="no") this.translationalWeightHistogram(i).mean()
			end if
		end do
		
		write(unit,"(5X,2F15.5)") histBuffer.mean(), histBuffer.stdev()
		
		histBuffer = RealHistogram()
		
		do i=1,this.numberOfExperiments
			call histBuffer.add( this.totalWeightHistogram(i).mean() )
			
			if( i == 1 ) then
				write(unit,"(10X,A5,F15.5)",advance="no") "LnW", this.totalWeightHistogram(i).mean()
			else
				write(unit,"(F15.5)",advance="no") this.totalWeightHistogram(i).mean()
			end if
		end do
		
		write(unit,"(5X,2F15.5)") histBuffer.mean(), histBuffer.stdev()
		
		write(unit,*) ""
		write(unit,*) ""
		
		write(unit,"(A)") "#------------------------------------"
		write(unit,"(A)") "# Markov chain statistics"
		write(unit,"(A)") "#------------------------------------"
		
		write(unit,"(A1,A)") "#", " Reactor type (ACCEPTED)"
! 		call this.reactorAcceptedHistogram.build()
		call this.reactorAcceptedHistogram.densityBegin( simIter )
		do while( associated(simIter) )
			srPair = this.reactorAcceptedHistogram.pair( simIter )
			write(unit,"(A20,F15.5)") srPair.first.fstr, srPair.second
			
			simIter => simIter.next
		end do
		
		write(unit,*) ""
		write(unit,"(A1,A)") "#", " Reactor type (REJECTED)"
! 		call this.reactorRejectedHistogram.build()
		call this.reactorRejectedHistogram.densityBegin( simIter )
		do while( associated(simIter) )
			srPair = this.reactorRejectedHistogram.pair( simIter )
			write(unit,"(A20,F15.5)") srPair.first.fstr, srPair.second
			
			simIter => simIter.next
		end do
		
		write(unit,*) ""
		write(unit,"(A1,A)") "#", " Reactor status"
! 		call this.reactorStatusHistogram.build()
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
		
		write(unit,"(A1,A15,5X,A)") "#", "prob.", "reaction"
		write(unit,"(A1,A15,5X,A)") "#", "-----", "--------"
! 		call this.transitionHistogram.build()
		call this.transitionHistogram.densityBegin( simIter )
		do while( associated(simIter) )
			srPair = this.transitionHistogram.pair( simIter )
			write(unit,"(1X,F15.5,5X,A)") srPair.second, srPair.first.fstr
			
			simIter => simIter.next
		end do
		
		write(unit,*) ""
		write(unit,*) ""
		
		write(unit,"(A1,A15,5X,A)") "#", "prob.", "reaction"
		write(unit,"(A1,A15,5X,A)") "#", "-----", "--------"
! 		call this.transitionDetHistogram.build()
		call this.transitionDetHistogram.densityBegin( simIter )
		do while( associated(simIter) )
			srPair = this.transitionDetHistogram.pair( simIter )
			write(unit,"(1X,F15.5,5X,A)") srPair.second, srPair.first.fstr
			
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
				iBuffer = FragmentsDB_instance.getIdClusterFromLabel( reactiveTokens(i) )
				call reactives.set( i, FragmentsDB_instance.clusters(iBuffer) )
			end do
		end if
		
		if( allocated(reactiveTokens) ) deallocate( reactiveTokens )
		
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
		
		sBuffer = iParser.getString( "MARKOV_CHAIN:energyHistoryFile", def="#@NONE@#" )
		if( trim(sBuffer.fstr) /= "#@NONE@#" ) then
			call this.energyHistoryFile.init( trim(sBuffer.fstr) )
			
			write(this.energyHistoryFile.unit,"(A)") "#------------------------------------"
			write(this.energyHistoryFile.unit,"(A)") "# ENERGY HISTORY"
			write(this.energyHistoryFile.unit,"(A)") "#------------------------------------"
			write(this.energyHistoryFile.unit,"(A1,3X,5A15,5X,A)") "#", "trans", "intermol", "vib", "rot", "tot", "label"
			write(this.energyHistoryFile.unit,"(A1,3X,5A15,5X,A)") "#", "eV", "eV", "eV", "eV", ""
			write(this.energyHistoryFile.unit,"(A1,3X,5A15,5X,A)") "#", "-------", "--------", "-------", "-----", "-------"
			
			write(*,"(A40,A)") "energyHistoryFile = ", trim(sBuffer.fstr)
		end if
		
		sBuffer = iParser.getString( "MARKOV_CHAIN:weightHistoryFile", def="#@NONE@#" )
		if( trim(sBuffer.fstr) /= "#@NONE@#" ) then
			call this.weightHistoryFile.init( trim(sBuffer.fstr) )
			
			write(this.weightHistoryFile.unit,"(A)") "#------------------------------------"
			write(this.weightHistoryFile.unit,"(A)") "# WEIGHT HISTORY"
			write(this.weightHistoryFile.unit,"(A)") "#------------------------------------"
			write(this.weightHistoryFile.unit,"(A1,3X,6A15,5X,A)") "#", "LnWe", "LnWv", "LnWn", "LnWr", "LnWt", "LnW", "label"
			write(this.weightHistoryFile.unit,"(A1,3X,6A15,5X,A)") "#", "arb.", "arb.", "arb.", "arb.", "arb.", "arb.", ""
			write(this.weightHistoryFile.unit,"(A1,3X,6A15,5X,A)") "#", "-------", "-------", "--------", "--------", "--------", "-------", "-------"
			
			write(*,"(A40,A)") "weightHistoryFile = ", trim(sBuffer.fstr)
		end if
		
		sBuffer = iParser.getString( "MARKOV_CHAIN:JHistoryFile", def="#@NONE@#" )
		if( trim(sBuffer.fstr) /= "#@NONE@#" ) then
			call this.JHistoryFile.init( trim(sBuffer.fstr) )
			
			write(this.JHistoryFile.unit,"(A)") "#------------------------------------"
			write(this.JHistoryFile.unit,"(A)") "# J HISTORY"
			write(this.JHistoryFile.unit,"(A)") "#------------------------------------"
			
			write(*,"(A40,A)") "JHistoryFile = ", trim(sBuffer.fstr)
		end if
		
		sBuffer = iParser.getString( "MARKOV_CHAIN:LHistoryFile", def="#@NONE@#" )
		if( trim(sBuffer.fstr) /= "#@NONE@#" ) then
			call this.LHistoryFile.init( trim(sBuffer.fstr) )
			
			write(this.LHistoryFile.unit,"(A)") "#------------------------------------"
			write(this.LHistoryFile.unit,"(A)") "# L HISTORY"
			write(this.LHistoryFile.unit,"(A)") "#------------------------------------"
			
			write(*,"(A40,A)") "LHistoryFile = ", trim(sBuffer.fstr)
		end if
		
		write(*,*)
		
		call this.run( reactives )
		
		histogramFile = iParser.getString( "MARKOV_CHAIN:histogramFile", def="#@NONE@#" )
		genEbkl = iParser.getLogical( "MARKOV_CHAIN:genEbkl", def=.false. )
		if( trim(histogramFile.fstr) /= "#@NONE@#" ) then
			call this.saveHistograms( histogramFile.fstr, genEbkl )
		end if
		
		if( this.energyHistoryFile.isOpen() ) then
			call this.energyHistoryFile.close()
		end if
		
		if( this.weightHistoryFile.isOpen() ) then
			call this.weightHistoryFile.close()
		end if
		
		if( this.JHistoryFile.isOpen() ) then
			call this.JHistoryFile.close()
		end if
		
		if( this.LHistoryFile.isOpen() ) then
			call this.LHistoryFile.close()
		end if
		
	end subroutine execute
end module MarkovChain_

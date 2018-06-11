!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                   !!
!! This file is part of M3C project                                                  !!
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
module FragmentsListBase_
	use GOptions_
	use String_
	use Math_
	use Matrix_
	use SpecialMatrix_
	use StringIntegerMap_
	use StringIntegerPair_
	use UnitsConverter_
	use IOStream_
	use RandomUtils_
	use RandomSampler_
	use AtomicElementsDB_
	use Atom_
	use Molecule_
	use BlocksIFileParser_
	use RealList_
	
	use GOptionsM3C_
	use Fragment_
	use FragmentsDB_
	
	implicit none
	private
	
	type, abstract, public :: FragmentsListBase
		type(Fragment), allocatable :: clusters(:)    !< Fragment list
		integer, allocatable :: idSorted(:)            !< Position of the molecule sorted by mass. It is calculated in updateFormula procedure
													   !<            this.clusters(i) ----> this.clusters( this.idSorted(i) )
		
		real(8) :: diagInertiaTensor(3)    !< Main inertia tensor components [ Ixx, Iyy, Izz ]. It is calculated into updateInertiaTensor procedure
		type(Matrix) :: inertiaAxes        !< Inertia axes considering all clusters
		
		integer :: nTrials_    !< Number of trials to generate the current geometry
		
		real(8) :: vibrationalEnergy_    !< It is calculated into LnWi procedure
		real(8) :: intermolEnergy_       !< It is calculated into changeGeometryFragmentsListBase procedure
		
		real(8) :: kineticEnergy_       !< Translational energy internally calculated by difference with the reactor energy
		real(8), private :: reactorEnergy_    !< Total energy into the reactor which is externally chosen
		
		real(8) :: logGFactor_    !< Undistinguished particles factor log( 1(Na!*Nb!...*Ni!) )
		integer, private :: ft_            !< Translational number of degrees of freedom
		integer, private :: fl_            !< Global number of rotational degrees of freedom (rigid-body)
		real(8) :: logVfree_      !< This is calculated only for the first time
		real(8) :: logVtheta_     !< This is calculated only for the first time
		real(8) :: logVJ_         !< This is calculated only for the first time
		
		real(8), private :: LnWe_    !< Electronic degeneration weight
		real(8), private :: LnWv_    !< Vibrational weight
		real(8), private :: LnWn_    !< Combinatorial weight
		
		character(200), private :: label_    !< Formula given by the user, for example C1
		character(200), private :: dlabel_   !< Detailed label given by the user, for example C1(sl)
		logical, private :: testLabel_       !< Checks that actualization of the label is made only one time
		
		logical :: forceRandomCenters    !< Fuerza a utilizar randomCenters en el siguiente llamado a changeGeometryFragmentsListBase
		logical :: forceInitializing     !< Fuerza a utilizar initialGuessFragmentsListBase
		
		integer :: totalComposition( AtomicElementsDB_nElems ) !< [ n1, n2, ..., nN ] n=numberOfAtomsWithZ, pos = atomicNumber 
		
		real(8) :: L_(3)          ! Orbital angular momentum
		
		integer, allocatable :: currentProducts(:,:) !< This allows to get the current potential energy surface @see updateIntermolecularPotential
		
		logical :: state !< True if there was not any problem generating its state
		
		contains
			generic :: init => initFragmentsListBase
			generic :: assignment(=) => copyFragmentsListBase
			
			procedure, NON_OVERRIDABLE :: initFragmentsListBase
			procedure, NON_OVERRIDABLE :: copyFragmentsListBase
			procedure, NON_OVERRIDABLE :: destroyFragmentsListBase
			procedure, NON_OVERRIDABLE :: save
			procedure, NON_OVERRIDABLE, private :: saveXYZ
			procedure, NON_OVERRIDABLE :: loadXYZ
			
			procedure, NON_OVERRIDABLE :: nMolecules
			procedure, NON_OVERRIDABLE :: nAtoms
			procedure, NON_OVERRIDABLE :: charge
			procedure, NON_OVERRIDABLE :: mass
			procedure, private :: updateFormula
			procedure, NON_OVERRIDABLE :: label
			generic :: set => setFromFragment
			procedure, NON_OVERRIDABLE :: setFromFragment
			
			procedure, NON_OVERRIDABLE :: setFrozen
			
			procedure, NON_OVERRIDABLE :: kineticEnergy
			procedure, NON_OVERRIDABLE :: setReactorEnergy
			procedure, NON_OVERRIDABLE :: reactorEnergy
			procedure, NON_OVERRIDABLE :: vibrationalEnergy
			procedure, NON_OVERRIDABLE :: intermolEnergy
			
			procedure, NON_OVERRIDABLE :: initialGuessFragmentsListBase
			procedure, private :: randomCenters
			procedure, private :: randomCentersByRandomWalkStep
			procedure, NON_OVERRIDABLE :: atomicOverlapping
			procedure, NON_OVERRIDABLE :: changeGeometryFragmentsListBase
			procedure, NON_OVERRIDABLE :: interpolateGeometryFragmentsListBase
			procedure, NON_OVERRIDABLE :: changeOrientationsFragmentsListBase
			procedure, NON_OVERRIDABLE :: changeVibrationalEnergyFragmentsListBase
			
			procedure, NON_OVERRIDABLE :: electronicEnergy
			procedure :: internalEnergy
			procedure :: totalEnergy
			procedure, NON_OVERRIDABLE :: ft
			procedure, NON_OVERRIDABLE :: fl
			procedure, NON_OVERRIDABLE :: fr
			
			procedure :: energyHistoryLine
			procedure :: weightHistoryLine
			procedure, NON_OVERRIDABLE :: JHistoryLine
			procedure, NON_OVERRIDABLE :: LHistoryLine
			
			procedure :: LnW !abstract y probablamente para quitar
			procedure, NON_OVERRIDABLE :: LnWe
			procedure, NON_OVERRIDABLE :: LnWv
			procedure, NON_OVERRIDABLE :: LnWn
			procedure, NON_OVERRIDABLE :: updateLnWe
			procedure, NON_OVERRIDABLE :: updateLnWv
			procedure, NON_OVERRIDABLE :: updateLnWn
			procedure :: showLnWComponents !abstract, revisar si no se puede reemplazar por una llamada a weightHistoryLine
			
			procedure, private :: updateLogGFactor ! Gibbs factor
			procedure, private :: updateLogVfree
			procedure, private :: updateLogVtheta
			procedure, private :: updateInertiaTensor
			procedure, private :: updateIntermolecularPotential
			procedure :: updateKineticEnergy
			
			procedure :: iTemperature
			procedure, NON_OVERRIDABLE :: centerOfMass
			procedure, NON_OVERRIDABLE :: buildInertiaTensor
			procedure, NON_OVERRIDABLE :: orient
			procedure, NON_OVERRIDABLE :: radius
			procedure, NON_OVERRIDABLE :: nFragments
			
			procedure, NON_OVERRIDABLE :: spinRange
			procedure, NON_OVERRIDABLE :: spinAvailable
			
			procedure, NON_OVERRIDABLE :: executeRadiusOptimization
	end type FragmentsListBase
	
	contains
	
	!>
	!! @brief Default constructor
	!!
	!! @param[in] nMolecules Number of clusters that will be contained
	!!
	subroutine initFragmentsListBase( this, nMolecules )
		class(FragmentsListBase) :: this 
		integer, intent(in) :: nMolecules
		
		if( allocated(this.clusters) ) deallocate(this.clusters)
		allocate( this.clusters(nMolecules) )
		
		if( allocated(this.idSorted) ) deallocate(this.idSorted)
		allocate( this.idSorted(nMolecules) )
		
		this.idSorted = -1.0_8
		
		this.diagInertiaTensor = 0.0_8
		this.inertiaAxes = SpecialMatrix_identity(3,3)
		
		this.nTrials_ = -1
		
		this.vibrationalEnergy_ = 0.0_8
		this.intermolEnergy_ = 0.0_8
		
		this.kineticEnergy_ = 0.0_8
		this.reactorEnergy_ = 0.0_8
		
		this.logGFactor_ = 0.0_8
		this.logVfree_ = 0.0_8
		this.logVtheta_ = 0.0_8
		this.ft_ = -1
		this.fl_ = -1
		
		this.LnWe_ = 0.0_8
		this.LnWv_ = 0.0_8
		this.LnWn_ = 0.0_8
		
		this.label_ = ""
		this.dlabel_ = ""
		this.testLabel_ = .false.
		
		this.forceRandomCenters = .true.
! 		this.forceInitializing = .true.
		this.forceInitializing = .false.
		
		this.totalComposition = -1
		
		this.L_ = 0.0_8
		
		if( allocated(this.currentProducts) ) deallocate(this.currentProducts)
		allocate(this.currentProducts(nMolecules,nMolecules) )
		
		this.state = .true.
	end subroutine initFragmentsListBase
	
	!>
	!! @brief Copy constructor
	!!
	subroutine copyFragmentsListBase( this, other )
		class(FragmentsListBase), intent(out) :: this
		class(FragmentsListBase), intent(in) :: other
		
		integer :: i
		
		if( allocated(this.clusters) ) deallocate(this.clusters)
		allocate( this.clusters( size(other.clusters) ) )
		
		do i=1,size(other.clusters)
			this.clusters(i) = other.clusters(i)
		end do
		
		if( allocated(this.idSorted) ) deallocate(this.idSorted)
		allocate( this.idSorted( size(other.idSorted) ) )
		
		this.idSorted = other.idSorted
		
		this.diagInertiaTensor = other.diagInertiaTensor
		this.inertiaAxes = other.inertiaAxes
		
		this.nTrials_ = other.nTrials_
		
		this.vibrationalEnergy_ = other.vibrationalEnergy_
		this.intermolEnergy_ = other.intermolEnergy_
		
		this.kineticEnergy_ = other.kineticEnergy_
		this.reactorEnergy_ = other.reactorEnergy_
		
		this.logGFactor_ = other.logGFactor_
		this.logVfree_ = other.logVfree_
		this.logVtheta_ = other.logVtheta_
		this.ft_ = other.ft_
		this.fl_ = other.fl_
		
		this.LnWe_ = other.LnWe_
		this.LnWv_ = other.LnWv_
		this.LnWn_ = other.LnWn_
		
		this.label_ = other.label_
		this.dlabel_ = other.dlabel_
		this.testLabel_ = other.testLabel_
		
		this.forceRandomCenters = other.forceRandomCenters
		this.forceInitializing = other.forceInitializing
		
		this.totalComposition = other.totalComposition
		
		this.L_ = other.L_
		
		if( allocated(this.currentProducts) ) deallocate(this.currentProducts)
		allocate( this.currentProducts( size(other.currentProducts,dim=1), size(other.currentProducts,dim=2) ) )
		this.currentProducts = other.currentProducts
		
		this.state = other.state
	end subroutine copyFragmentsListBase
	
	!>
	!! @brief Destructor
	!!
	subroutine destroyFragmentsListBase( this )
		class(FragmentsListBase) :: this
		
		if( allocated(this.clusters) ) deallocate(this.clusters)
		if( allocated(this.idSorted) ) deallocate(this.idSorted)
		
		if( allocated(this.currentProducts) ) deallocate(this.currentProducts)
	end subroutine destroyFragmentsListBase
	
	!>
	!! @brief Save the molecule to file
	!!
	subroutine save( this, fileName, format, append )
		class(FragmentsListBase), intent(in) :: this
		character(*), optional, intent(in) :: fileName
		integer, optional, intent(in) :: format
		logical, optional, intent(in) :: append
		
		integer :: effFormat
		
		effFormat = XYZ
		if( present(format) ) effFormat = format
		
		select case ( effFormat )
			case( XYZ )
				call this.saveXYZ( fileName, append )
		end select
	end subroutine save
	
	!>
	!! @brief Save the molecule to file
	!!
	subroutine saveXYZ( this, fileName, append )
		class(FragmentsListBase), intent(in) :: this
		character(*), intent(in) :: fileName
		logical, optional, intent(in) :: append
		
		logical :: effAppend
		
		type(OFStream) :: ofile
		integer :: i, j, molCounter
		
		effAppend = .false.
		if( present(append) ) effAppend = append
		
		call ofile.init( fileName, append=effAppend )
		
		write(ofile.unit,*) this.nAtoms()
		write(ofile.unit,*) "FragmentsListBase"
		
		molCounter = 1
		do i=1,this.nMolecules()
			do j=1,this.clusters(i).nAtoms()
! 				write(ofile.unit,'(A5,3F20.8,I10)') &
! 					this.clusters(i).atoms(j).symbol, &
! 					this.clusters(i).atoms(j).x/angs, &
! 					this.clusters(i).atoms(j).y/angs, &
! 					this.clusters(i).atoms(j).z/angs, &
! 					molCounter
					
				write(ofile.unit,'(A5,3F20.8)') &
					this.clusters(i).atoms(j).symbol, &
					this.clusters(i).atoms(j).x/angs, &
					this.clusters(i).atoms(j).y/angs, &
					this.clusters(i).atoms(j).z/angs
					
				if( j == this.clusters(i).nAtoms() ) then
					molCounter = molCounter + 1
				end if
			end do
		end do
		
		call ofile.close()
	end subroutine saveXYZ
	
	!>
	!! @brief Load the molecule from file
	!!
	subroutine loadXYZ( this, fileName )
		class(FragmentsListBase) :: this
		character(*), intent(in) :: fileName
		
		type(IFStream) :: ifile
		type(String) :: buffer
		character(1000), allocatable :: tokens(:)
		integer :: i
		integer :: nMols
		integer :: idClus
		real(8) :: x, y, z
		type(String) :: name
		
		if( allocated(this.clusters) ) deallocate(this.clusters)
		
		call ifile.init( trim(fileName) )
		
		if( ifile.numberOfLines < 3 ) then
			write(*,*) "### ERROR ### bad format in XYZ file"
			stop
		end if
		
		buffer = ifile.readLine()
		call buffer.split( tokens, " " )
		
		if( buffer.length() == 0 .or. size(tokens)<1 ) then
			call GOptions_error( &
				"Inconsistent XYZ file, line "//FString_fromInteger(ifile.currentLine), &
				"FragmentsListBase.loadXYZ()", &
				"Check the XYZ format file for FragmentsList" &
				)
		end if
		
		nMols = FString_toInteger(tokens(1))
		
		call this.init( nMols )
		
		buffer = ifile.readLine() ! It loads the name, but in this case it's not important
			
		do i=1,nMols
			if( ifile.eof() ) then
				call GOptions_error( &
					"Inconsistent number of clusters in XYZ file, line "//FString_fromInteger(ifile.currentLine), &
					"FragmentsListBase.loadXYZ()", &
					"Check the XYZ format file for FragmentsList" &
					)
			end if
			
			buffer = ifile.readLine()
			
			call buffer.split( tokens, " " )
			
			if( buffer.length() == 0 .or. size(tokens)<4 ) then
				call GOptions_error( &
					"Inconsistent XYZ file (4 tokens are neccesary), line "//FString_fromInteger(ifile.currentLine), &
					"FragmentsListBase.loadXYZ()", &
					"Check the XYZ format file for FragmentsList" &
					)
				stop
			end if
			
			name = tokens(1)
			x = FString_toReal( tokens(2) )*angs
			y = FString_toReal( tokens(3) )*angs
			z = FString_toReal( tokens(4) )*angs
			
			idClus = FragmentsDB_instance.getIdFromName( name.fstr )
			
			call this.setFromFragment( i, FragmentsDB_instance.clusters( idClus ) )
			call this.clusters(i).setCenter( [ x, y, z ] )
		end do
		
		deallocate( tokens )
		call ifile.close()
		
		call this.orient()
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Ya no es necesario generar centros aleatorios
		! pues se han leido del fichero de gemetría
		this.forceRandomCenters = .false.
		
	end subroutine loadXYZ
	
	!>
	!! @brief Save the molecule to file
	!!
	pure function nMolecules( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		integer :: output
		
		output = size( this.clusters )
	end function nMolecules
	
	!>
	!! @brief Returns the number of atoms
	!!
	pure function nAtoms( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		integer :: output
		
		integer :: i
		
		output = 0
		do i=1,this.nMolecules()
			output = output + this.clusters(i).nAtoms()
		end do
	end function nAtoms
	
	!>
	!! @brief Returns the total charge
	!!
	pure function charge( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		integer :: output
		
		integer :: i
		
		output = 0
		do i=1,this.nMolecules()
			output = output + this.clusters(i).charge
		end do
	end function charge
	
	!>
	!! @brief Returns the total mass
	!!
	pure function mass( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		real(8) :: output
		
		integer :: i
		
		output = 0.0_8
		do i=1,this.nMolecules()
			output = output + this.clusters(i).mass()
		end do
	end function mass
	
	!>
	!! @brief
	!!
	subroutine updateFormula( this )
		class(FragmentsListBase) :: this
		
		integer :: i
		real(8), allocatable :: massVec(:) ! @todo Hay que hacer el Math_sort para enteros
		
		allocate( massVec(size(this.clusters)) )
		
		this.totalComposition = 0
		do i=1,size(this.clusters)
		
			massVec(i) = this.clusters(i).mass()/amu &
					+ this.clusters(i).charge/10.0 &
						+ this.clusters(i).multiplicity/100.0_8 

! 			massVec(i) = 10000000*this.clusters(i).mass()/amu + 1000000*this.clusters(i).nAtoms() + 100*this.clusters(i).charge &
! 					+ this.clusters(i).multiplicity

			this.totalComposition = this.totalComposition + this.clusters( i ).composition
		end do
			
		call Math_sort( massVec, this.idSorted )
		
		this.label_ = ""
		this.dlabel_ = ""
		do i=1,size(this.clusters)
			if( i /= size(this.clusters) ) then
				this.label_ = trim(this.label_)//trim(this.clusters( this.idSorted(i) ).label( details=.false. ))//"+"
				this.dlabel_ = trim(this.dlabel_)//trim(this.clusters( this.idSorted(i) ).label( details=.true. ))//"+"
			else
				this.label_ = trim(this.label_)//trim(this.clusters( this.idSorted(i) ).label( details=.false. ))
				this.dlabel_ = trim(this.dlabel_)//trim(this.clusters( this.idSorted(i) ).label( details=.true. ))
			end if
		end do
		
		this.testLabel_ = .true.
		
		deallocate( massVec )
	end subroutine updateFormula
	
	!>
	!! @brief
	!!
	function label( this, details ) result( output )
		class(FragmentsListBase) :: this 
		character(200) :: output
		logical, optional :: details
		
		logical :: effDetails
		
		effDetails = .true.
		if( present(details) ) effDetails = details
		
		if( .not. this.testLabel_ ) then
			call this.updateFormula()
		end if
		
		if( effDetails ) then
			output = this.dlabel_
		else
			output = this.label_
		end if
	end function label
	
	!>
	!! @brief
	!!
	subroutine setFromFragment( this, pos, clus )
		class(FragmentsListBase) :: this
		integer, intent(in) :: pos
		class(Fragment), intent(in) :: clus
		
		if( GOptions_printLevel >= 3 ) then
			write(*,"(<GOptions_indentLength*2>X,A)") "Added cluster "//trim(clus.label())
		end if
		
		if( size(this.clusters) >= pos ) then
			this.clusters(pos) = clus
		else
			write(*,*) "### ERROR ### FragmentsListBase.setFromFragment: ( pos > size )"
			stop
		end if
		
		this.testLabel_ = .false.
	end subroutine setFromFragment
	
	!>
	!! @brief
	!!
	subroutine setFrozen( this, frozen )
		class(FragmentsListBase) :: this
		logical, intent(in) :: frozen
		
		integer :: i
		
		do i=1,this.nMolecules()
			this.clusters(i).frozen = frozen
		end do
	end subroutine setFrozen
	
	!>
	!! @brief
	!!
	subroutine setKineticEnergy( this, kineticEnergy )
		class(FragmentsListBase) :: this
		real(8), intent(in) :: kineticEnergy
		
		this.kineticEnergy_ = kineticEnergy
	end subroutine setKineticEnergy
	
	!>
	!! @brief
	!!
	pure function kineticEnergy( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		real(8) :: output
		
		output = this.kineticEnergy_
	end function kineticEnergy
	
	!>
	!! @brief
	!!
	subroutine setReactorEnergy( this, reactorEnergy )
		class(FragmentsListBase) :: this
		real(8), intent(in) :: reactorEnergy
		
		this.reactorEnergy_ = reactorEnergy
		
		call this.updateKineticEnergy()
	end subroutine setReactorEnergy
	
	!>
	!! @brief
	!!
	pure function reactorEnergy( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		real(8) :: output
		
		output = this.reactorEnergy_
	end function reactorEnergy
	
	!>
	!! @brief
	!!
	pure function vibrationalEnergy( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		real(8) :: output
		
		output = this.vibrationalEnergy_
	end function vibrationalEnergy
	
	!>
	!! @brief
	!!
	pure function intermolEnergy( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		real(8) :: output
		
		output = this.intermolEnergy_
	end function intermolEnergy
	
	!>
	!! @brief
	!!
	subroutine initialGuessFragmentsListBase( this )
		class(FragmentsListBase) :: this
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Esto debe ser los primero que haga esta subrutina
		! ya que los métodos change llaman esta función cuando
		! forceInitializing=.true. y si este no se cambia
		! se generaría un ciclo infinito
		this.forceInitializing = .false.
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		call this.orient()
		
		call this.updateFormula()
		
		call this.updateLogGFactor()
		
		call this.updateLnWn()
		call this.updateLnWe()
		
! 		call this.changeGeometryFragmentsListBase()
! 		call this.changeOrientationsFragmentsListBase()
		
		call this.updateLogVfree() ! El metodo 1 necesita this.nTrials_ de changeGeometryFragmentsListBase
		call this.updateLogVtheta()
		
! 		call this.changeVibrationalEnergyFragmentsListBase()
	end subroutine initialGuessFragmentsListBase
	
	!>
	!! @brief
	!! @param[in]  maxIter
	!!             Maximum allowed number of iterations (default = 1000000)
	!!
	subroutine randomCenters( this, maxIter )
		class(FragmentsListBase) :: this
		integer, optional, intent(in) ::  maxIter
		
		integer :: effMaxIter
		
		integer :: nTrials_
		type(RandomSampler) :: rs
		real(8), allocatable :: sample(:,:)
		real(8) :: rVec1(3), rVec2(3), cm(3)
		integer :: i, j, n
		logical :: overlap
		
		this.state = .true.
		
		effMaxIter = 1000
		if( present(maxIter) ) effMaxIter = maxIter
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_subsection( "Random centers "//trim(this.label()), indent=2 )
			call GOptions_valueReport( "maxIter", effMaxIter, indent=2 )
		end if
		
! Testing overlap
#define OVERLAPPING(i,j) this.clusters(i).radius( type=GOptionsM3C_radiusType )+this.clusters(j).radius( type=GOptionsM3C_radiusType )-GOptionsM3C_overlappingRadius > norm2( rVec2-rVec1 )
! ! Sampling on "x" axis, which is the axis with inertia moment equal to cero
! #define RVEC_X(i) [ sample(1,i), 0.0_8, 0.0_8 ]
! ! Sampling on "x-y" axes in polar coordinates
! #define RVEC_XY(i) [ sample(1,i)*cos(sample(2,i)), sample(1,i)*sin(sample(2,i)), 0.0_8 ]
! ! Sampling on "x-y-z" axes in spherical coordinates
! #define RVEC_XYZ(i) [ sample(1,i)*sin(sample(2,i))*cos(sample(3,i)), sample(1,i)*sin(sample(2,i))*sin(sample(3,i)), sample(1,i)*cos(sample(2,i)) ]

! Sampling on "x" axis, which is the axis with inertia moment equal to cero
#define RVEC_X(i) [ sample(1,i), 0.0_8, 0.0_8 ]
! Sampling on "x-y" axes in polar coordinates
#define RVEC_XY(i) [ sample(1,i), sample(2,i), 0.0_8 ]
! Sampling on "x-y-z" axes in spherical coordinates
#define RVEC_XYZ(i) [ sample(1,i), sample(2,i), sample(3,i) ]
		
		allocate( sample(3,this.nMolecules()) )
		
		call rs.init( nDim=3 )
! 		call rs.setRange( 1, [0.0_8,GOptionsM3C_systemRadius] )     ! r in (0,Rsys)
! 		call rs.setRange( 2, [0.0_8,MATH_PI] )        ! theta in (0,pi)
! 		call rs.setRange( 3, [0.0_8,2.0_8*MATH_PI] )  ! phi in (0,2pi)
		call rs.setRange( 1, [-GOptionsM3C_systemRadius,GOptionsM3C_systemRadius] )  ! x in (-Rsys,Rsys)
		call rs.setRange( 2, [-GOptionsM3C_systemRadius,GOptionsM3C_systemRadius] )  ! y in (-Rsys,Rsys)
		call rs.setRange( 3, [-GOptionsM3C_systemRadius,GOptionsM3C_systemRadius] )  ! z in (-Rsys,Rsys)
		
		overlap = .false.
		
		select case( this.nMolecules() )
			case( 1 )
				call this.clusters(1).setCenter( [0.0_8, 0.0_8, 0.0_8] )
				nTrials_ = 1
				
			case( 2 )
				
				nTrials_ = 0
				do n=1,effMaxIter
					call rs.uniform( sample )
					
					rVec1 = RVEC_X(1)
					rVec2 = RVEC_X(2)
					
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					! Se verifica que a los centros no sobrelapen ni que
					! se salgan del radio del sistema
					overlap = OVERLAPPING(1,2) .or. norm2(rVec1) > GOptionsM3C_systemRadius .or. norm2(rVec2) > GOptionsM3C_systemRadius
					
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					! Se verifica que a ajustar el centro de masas el sistema
					! no se salga del volumen de simulación
					cm = 0.0_8
					do i=1,this.nMolecules()
						cm = cm + this.clusters(i).mass()*RVEC_X(i)
					end do
					cm = cm/this.mass()
					
					do i=1,this.nMolecules()
						if( norm2(RVEC_X(i)-cm) > GOptionsM3C_systemRadius ) then
							overlap = .true.
						end if
					end do
					
					nTrials_ = nTrials_ + 1
					if( .not. overlap ) exit
				end do
				
				do i=1,this.nMolecules()
					call this.clusters(i).setCenter( RVEC_X(i) )
				end do
				
			case( 3 )
				
				nTrials_ = 0
				do n=1,effMaxIter
					call rs.uniform( sample )
					
					overlap = .false.
					
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					! Se verifica que a los centros no sobrelapen ni que
					! se salgan del radio del sistema
					do i=1,this.nMolecules()-1
						do j=i+1,this.nMolecules()
							
							rVec1 = RVEC_XY(i)
							rVec2 = RVEC_XY(j)
							
							overlap = OVERLAPPING(i,j) .or. norm2(rVec1) > GOptionsM3C_systemRadius .or. norm2(rVec2) > GOptionsM3C_systemRadius
							
							if( overlap ) exit
						end do
						
						if( overlap ) exit
					end do
					
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					! Se verifica que a ajustar el centro de masas el sistema
					! no se salga del volumen de simulación
					cm = 0.0_8
					do i=1,this.nMolecules()
						cm = cm + this.clusters(i).mass()*RVEC_XY(i)
					end do
					cm = cm/this.mass()
					
					do i=1,this.nMolecules()
						if( norm2(RVEC_XY(i)-cm) > GOptionsM3C_systemRadius ) then
							overlap = .true.
						end if
					end do
					
					nTrials_ = nTrials_ + 1
					if( .not. overlap ) exit
				end do
				
				do i=1,this.nMolecules()
					call this.clusters(i).setCenter( RVEC_XY(i) )
				end do
				
			case default
				
				nTrials_ = 0
				do n=1,effMaxIter
					call rs.uniform( sample )
					
					overlap = .false.
					
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					! Se verifica que a los centros no sobrelapen ni que
					! se salgan del radio del sistema
					do i=1,this.nMolecules()-1
						rVec1 = RVEC_XYZ(i)
						
						overlap = norm2(rVec1) > GOptionsM3C_systemRadius
						if( overlap ) exit
						
						do j=i+1,this.nMolecules()
							
							rVec2 = RVEC_XYZ(j)
							
							overlap = OVERLAPPING(i,j) .or. norm2(rVec2) > GOptionsM3C_systemRadius
							
							if( overlap ) exit
						end do
						
						if( overlap ) exit
					end do
					
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					! Se verifica que a ajustar el centro de masas el sistema
					! no se salga del volumen de simulación
					cm = 0.0_8
					do i=1,this.nMolecules()
						cm = cm + this.clusters(i).mass()*RVEC_XYZ(i)
					end do
					cm = cm/this.mass()
					
					do i=1,this.nMolecules()
						if( norm2(RVEC_XYZ(i)-cm) > GOptionsM3C_systemRadius ) then
							overlap = .true.
						end if
					end do
					
					nTrials_ = nTrials_ + 1
					if( .not. overlap ) exit
				end do
				
				do i=1,this.nMolecules()
					call this.clusters(i).setCenter( RVEC_XYZ(i) )
				end do
				
		end select
#undef OVERLAPPING
#undef RVEC_X
#undef RVEC_XY
#undef RVEC_XYZ
			
		this.nTrials_ = nTrials_
		deallocate( sample )
			
		if( overlap ) then
			this.state = .false.
! 			call GOptions_error( &
! 				"Maximum number of iterations reached"//" (n = "//trim(FString_fromInteger(n))//")", &
! 				"FragmentsListBase.randomCenters()", &
! 				"Consider to increase GOptions:systemRadius" &
! 				)
		end if
		
	end subroutine randomCenters

	!>
	!! @brief
	!! @param[in]  maxIter
	!!             Maximum allowed number of iterations (default = 100000)
	!!
	subroutine randomCentersByRandomWalkStep( this, maxIter )
		class(FragmentsListBase) :: this
		integer, optional, intent(in) ::  maxIter
		
		integer :: effMaxIter
		
		integer :: nTrials_
		type(RandomSampler) :: rs
		real(8), allocatable :: sample(:,:)
		real(8) :: rVec1(3), rVec2(3), dR(3)
		integer :: i, j, n, lastCase
		logical :: overlap
		
		this.state = .true.
		
		effMaxIter = 1000
		if( present(maxIter) ) effMaxIter = maxIter
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_paragraph( "RANDOM WALK STEP", indent=2 )
			call GOptions_valueReport( "radius", GOptionsM3C_randomWalkStepRadius/angs, "A", indent=2 )
			call GOptions_valueReport( "systemRadius", GOptionsM3C_systemRadius/angs, "A", indent=2 )
			call GOptions_valueReport( "maxIter", effMaxIter, indent=2 )
		end if
		
! Testing overlap
#define OVERLAPPING(i,j) this.clusters(i).radius( type=GOptionsM3C_radiusType )+this.clusters(j).radius( type=GOptionsM3C_radiusType )-GOptionsM3C_overlappingRadius > norm2( rVec2-rVec1 )
! Sampling on "x" axis, which is the axis with inertia moment equal to cero (r)
#define RVEC_X(i) [ sample(1,i), 0.0_8, 0.0_8 ]
! Sampling on "x-y" axes in polar coordinates (r,phi)
#define RVEC_XY(i) [ sample(1,i)*cos(sample(3,i)), sample(1,i)*sin(sample(3,i)), 0.0_8 ]
! Sampling on "x-y-z" axes in spherical coordinates (r,theta,phi)
#define RVEC_XYZ(i) [ sample(1,i)*sin(sample(2,i))*cos(sample(3,i)), sample(1,i)*sin(sample(2,i))*sin(sample(3,i)), sample(1,i)*cos(sample(2,i)) ]

#define SCALE_BY_MASS(i) this.clusters(i).mass()/this.mass()
! #define SCALE_BY_MASS(i) max( 0.0_8, log(1.01_8*this.mass()-this.clusters(i).mass())/log(this.mass()) )

		allocate( sample(3,this.nMolecules()) )
		
		call rs.init( nDim=3 )
		call rs.setRange( 1, [0.0_8,GOptionsM3C_randomWalkStepRadius] )     ! r in (0,dR)
		call rs.setRange( 2, [0.0_8,MATH_PI] )        ! theta in (0,pi)
		call rs.setRange( 3, [0.0_8,2.0_8*MATH_PI] )  ! phi in (0,2pi)
		
		overlap = .false.
		
		select case( this.nMolecules() )
			case( 1 )
				lastCase = 1
				call this.clusters(1).setCenter( [0.0_8, 0.0_8, 0.0_8] )
				
			case( 2 )
				lastCase = 2
				! Se genera un configuración cercana no solapante
				! y dentro del radio del sistema
				do n=1,effMaxIter
					call rs.uniform( sample )
					
					if( GOptionsM3C_useWeightedWalkStep ) then
						rVec1 = this.clusters(1).center()+RVEC_X(1)*( 1.0_8-SCALE_BY_MASS(1) )
						rVec2 = this.clusters(2).center()+RVEC_X(2)*( 1.0_8-SCALE_BY_MASS(2) )
					else
						rVec1 = this.clusters(1).center()+RVEC_X(1)
						rVec2 = this.clusters(2).center()+RVEC_X(2)
					end if
						
					overlap = OVERLAPPING(1,2) .or. norm2(rVec1) > GOptionsM3C_systemRadius .or. norm2(rVec2) > GOptionsM3C_systemRadius
					
					if( .not. overlap ) exit
				end do
				
				do i=1,this.nMolecules()
					if( GOptionsM3C_useWeightedWalkStep ) then
						call this.clusters(i).setCenter( this.clusters(i).center()+RVEC_X(i)*( 1.0_8-SCALE_BY_MASS(i) ) )
					else
						call this.clusters(i).setCenter( this.clusters(i).center()+RVEC_X(i) )
					end if
				end do
				
			case( 3 )
				lastCase = 3
				! Se genera un configuración cercana no solapante
				! y dentro del radio del sistema
				do n=1,effMaxIter
					call rs.uniform( sample )
					
					overlap = .false.
					
					do i=1,this.nMolecules()-1
						do j=i+1,this.nMolecules()
							
							if( GOptionsM3C_useWeightedWalkStep ) then
								rVec1 = this.clusters(i).center()+RVEC_XY(i)*( 1.0_8-SCALE_BY_MASS(i) )
								rVec2 = this.clusters(j).center()+RVEC_XY(j)*( 1.0_8-SCALE_BY_MASS(j) )
							else
								rVec1 = this.clusters(i).center()+RVEC_XY(i)
								rVec2 = this.clusters(j).center()+RVEC_XY(j)
							end if
							
							overlap = OVERLAPPING(i,j) .or. norm2(rVec1) > GOptionsM3C_systemRadius .or. norm2(rVec2) > GOptionsM3C_systemRadius
							
							if( overlap ) exit
						end do
						
						if( overlap ) exit
					end do
					
					if( .not. overlap ) exit
				end do
				
				do i=1,this.nMolecules()
					if( GOptionsM3C_useWeightedWalkStep ) then
						call this.clusters(i).setCenter( this.clusters(i).center()+RVEC_XY(i)*( 1.0_8-SCALE_BY_MASS(i) ) )
					else
						call this.clusters(i).setCenter( this.clusters(i).center()+RVEC_XY(i) )
					end if
				end do
				
			case default
				lastCase = 4
				! Se genera un configuración cercana no solapante
				! y dentro del radio del sistema
				do n=1,effMaxIter
					call rs.uniform( sample )
					
					overlap = .false.
					
					do i=1,this.nMolecules()-1
						do j=i+1,this.nMolecules()
							
							if( GOptionsM3C_useWeightedWalkStep ) then
								rVec1 = this.clusters(i).center()+RVEC_XYZ(i)*( 1.0_8-SCALE_BY_MASS(i) )
								rVec2 = this.clusters(j).center()+RVEC_XYZ(j)*( 1.0_8-SCALE_BY_MASS(j) )
							else
								rVec1 = this.clusters(i).center()+RVEC_XYZ(i)
								rVec2 = this.clusters(j).center()+RVEC_XYZ(j)
							end if
							
							overlap = OVERLAPPING(i,j) .or. norm2(rVec1) > GOptionsM3C_systemRadius .or. norm2(rVec2) > GOptionsM3C_systemRadius
							
							if( overlap ) exit
						end do
						
						if( overlap ) exit
					end do
					
					if( .not. overlap ) exit
				end do
				
				do i=1,this.nMolecules()
					if( GOptionsM3C_useWeightedWalkStep ) then
						call this.clusters(i).setCenter( this.clusters(i).center()+RVEC_XYZ(i)*( 1.0_8-SCALE_BY_MASS(i) ) )
					else
						call this.clusters(i).setCenter( this.clusters(i).center()+RVEC_XYZ(i) )
					end if
				end do
				
		end select
#undef OVERLAPPING
#undef RVEC_X
#undef RVEC_XY
#undef RVEC_XYZ
#undef SCALE_BY_MASS
		
		deallocate(sample)
		
		if( overlap ) then
			this.state = .false.
! 			call GOptions_error( &
! 				"Maximum number of iterations reached"//" (n = "//trim(FString_fromInteger(n))//")", &
! 				"FragmentsListBase.randomCentersByRandomWalkStep() case="//trim(FString_fromInteger(lastCase)), &
! 				"Consider to increase GOptions:randomWalkStepRadius ( "//trim(this.label(details=.true.))//" )" &
! 				)
		end if
	end subroutine randomCentersByRandomWalkStep
	
	!>
	!! @brief Checks for atomic overlapping
	!!
	pure function atomicOverlapping( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		logical :: output
		
		integer :: i, j
		
		if( .not. GOptionsM3C_checkAtomicOverlapping ) then
			output = .false.
			return
		end if
		
		output = .false.
		do i=1,this.nMolecules()-1
			do j=i+1,this.nMolecules()
				
				output = this.clusters(i).overlapping( this.clusters(j), &
						overlappingRadius=GOptionsM3C_atomicOverlappingRadius, radiusType=GOptionsM3C_radiusType )
				
				if( output ) return
			end do
		end do
		
	end function atomicOverlapping
	
	!>
	!! @brief Builds a random intermolecular energy by building a random
	!!        geometry for the clusters. The center of mass is fixed as
	!!        the origin and the main inertia directions as axes of the
	!!        coordinates system.
	!!
	subroutine changeGeometryFragmentsListBase( this )
		class(FragmentsListBase) :: this
		
		logical :: effInitialize
		
		real(8) :: rVec1(3), rVec2(3)
		integer :: i, j, n, m
		real(8) :: centerOfMass(3)
		logical :: check
		
		real(8), allocatable :: geomBackup(:,:)
		
		allocate( geomBackup(3,this.nMolecules()) )
		
		do i=1,this.nMolecules()
			geomBackup(i,:) = this.clusters(i).center()
		end do
		
		this.state = .true.
		
! 		if( this.forceInitializing ) then
! 			call this.initialGuessFragmentsListBase()
! 			return
! 		end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Por omisión se generan centros aleatorios
		! la primera vez que entra a esta función
		if( this.forceRandomCenters .or. ( .not. GOptionsM3C_useRandomWalkers ) ) then
			
			write(*,"(A)",advance="no") "Random ... "
			check = .true.
			do n=1,1000
				if( GOptions_printLevel >= 2 ) write(*,*) "Changing geometry avoiding atomic overlapping (step="//trim(FString_fromInteger(n))//")"
				
				call this.randomCenters()
				write(*,"(A)",advance="no") ":"
				
				if( .not. this.atomicOverlapping() ) then
					check = .false.
					write(*,"(A)",advance="no") " OK"
					exit
				end if
			end do
			
! 			this.forceRandomCenters = .false.  !@todo forceRandomCenters debe desaparecer. Si no hay una configuración simplemente this.state=.false.
			if( check ) then
				this.state = .false.
				
				call this.save("hola.xyz")
				read(*,*)
				
				do i=1,this.nMolecules()
					call this.clusters(i).setCenter( geomBackup(i,:) )
				end do
			end if
		else
			do m=1,10
				check = .true.
				do n=1,1000
					if( GOptions_printLevel >= 2 ) write(*,*) "Changing geometry avoiding atomic overlapping (trial="//trim(trim(FString_fromInteger(m)))//",step="//trim(FString_fromInteger(n))//")"
					
					call this.randomCentersByRandomWalkStep()
					
					if( this.radius() < GOptionsM3C_systemRadius .and. .not. this.atomicOverlapping() ) then
						check = .false.
						exit
					end if
				end do
				
! 				if( check ) then
! 					call this.randomCenters()
! 				else
! 					exit
! 				end if
			end do
				
			if( check ) then
				this.state = .false.
! 				call GOptions_error( &
! 				"Maximum number of iterations reached"//" (n = "//trim(FString_fromInteger(n-1))//", m="//trim(FString_fromInteger(m-1))//")", &
! 				"FragmentsListBase.changeGeometryFragmentsListBase()", &
! 				"Consider to increase GOptions:systemRadius" &
! 				)
			end if
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Se verifica que a ajustar el centro de masas el sistema
			! no se salga del volumen de simulación
			centerOfMass = this.centerOfMass()
			do i=1,this.nMolecules()
				if( norm2(this.clusters(i).center()-centerOfMass) > GOptionsM3C_systemRadius ) then
					call this.randomCenters()
					exit
				end if
			end do
		end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se establece como origen el centro de masas y se
		! actualiza el tensor de inercia
		centerOfMass = this.centerOfMass()
		do i=1,this.nMolecules()
			call this.clusters(i).setCenter( this.clusters(i).center() - centerOfMass )
		end do
		
		call this.updateInertiaTensor()
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_paragraph( "Geometry", indent=2 )
			
			call GOptions_valueReport( "CM", this.centerOfMass(), "A", indent=2 )
			write(*,*) ""
			
			write(*,"(<GOptions_indentLength*2>X,A20,3X,3A10,3X,A10)") "id", "X", "Y", "Z", "R"
			write(*,"(<GOptions_indentLength*2>X,A33,2A10,3X,A10)") "A", "A", "A", "A"
			do i=1,this.nMolecules()
				write(*,"(<GOptions_indentLength*2>X,A20,3X,3F10.5,3X,F10.5)") trim(this.clusters(i).label()), &
						this.clusters(i).center()/angs, this.clusters(i).radius( type=GOptionsM3C_radiusType )/angs
			end do
			write(*,*) ""
		end if
		
		call this.updateIntermolecularPotential()
		
		deallocate( geomBackup )
	end subroutine changeGeometryFragmentsListBase
	
	!>
	!! @brief Interpolates the geometry from other
	!!
	subroutine interpolateGeometryFragmentsListBase( this, other )
		class(FragmentsListBase) :: this
		class(FragmentsListBase) :: other
		
		logical :: check
		integer :: i, n
		real(8) :: centerOfMass(3)
		
		! @todo No tengo claro que esto deba ir aqui
! 		if( this.forceInitializing ) then
! 			call this.initialGuessFragmentsListBase()
! 			return
! 		end if

		if( this.nFragments() == other.nFragments() ) then
		
			do i=1,this.nMolecules()
				call this.clusters(i).setCenter( other.clusters(i).center() )
			end do
			
			if( GOptions_printLevel >= 2 ) write(*,*) "Geometry interpolated ... OK"
			
! 			call this.changeOrientationsFragmentsListBase()
			
			! @todo No se por que tengo este bloque. Pero check nunca es asignado
! 			if( check ) then
! 				call GOptions_error( &
! 				"Maximum number of iterations reached"//" (n = "//trim(FString_fromInteger(n-1))//")", &
! 				"FragmentsListBase.interpolateGeometryFragmentsListBase()", &
! 				"Consider to change the initial geometry" &
! 				)
! 			end if
		
		else if( this.nFragments() > other.nFragments() ) then
			
			call GOptions_error( &
				"Case this.nFragments() > other.nFragments() is not implemented yet", &
				"FragmentsListBase.geometryFrom()" &
				)
				
		else if( this.nFragments() < other.nFragments() ) then
		
			call GOptions_error( &
				"Case this.nFragments() < other.nFragments() is not implemented yet", &
				"FragmentsListBase.geometryFrom()" &
				)
				
		end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se establece como origen el centro de masas y se
		! actualiza el tensor de inercia
		centerOfMass = this.centerOfMass()
		do i=1,this.nMolecules()
			call this.clusters(i).setCenter( this.clusters(i).center() - centerOfMass )
		end do
		
		call this.updateInertiaTensor()
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_paragraph( "Geometry", indent=2 )
			
			call GOptions_valueReport( "CM", this.centerOfMass(), "A", indent=2 )
			write(*,*) ""
			
			write(*,"(<GOptions_indentLength*2>X,A20,3X,3A10,3X,A10)") "id", "X", "Y", "Z", "R"
			write(*,"(<GOptions_indentLength*2>X,A33,2A10,3X,A10)") "A", "A", "A", "A"
			do i=1,this.nMolecules()
				write(*,"(<GOptions_indentLength*2>X,A20,3X,3F10.5,3X,F10.5)") trim(this.clusters(i).label()), &
						this.clusters(i).center()/angs, this.clusters(i).radius( type=GOptionsM3C_radiusType )/angs
			end do
			write(*,*) ""
		end if
		
		call this.updateIntermolecularPotential()
	end subroutine interpolateGeometryFragmentsListBase
	
	!>
	!! @brief 
	!!
	subroutine changeOrientationsFragmentsListBase( this )
		class(FragmentsListBase) :: this
		
		integer :: i, n
		logical :: check
		
		this.state = .true.
		
! 		if( this.forceInitializing ) then
! 			call this.initialGuessFragmentsListBase()
! 			return
! 		end if
		
		if( this.nMolecules() > 1 ) then
			
			if( GOptions_printLevel >= 2 ) write(*,"(A)", advance="no") "Checking overlaping"
			
			check = .true. 
			do n=1,1000
				
				if( GOptions_printLevel >= 2 ) write(*,*) "Changing orientations avoiding atomic overlapping (step="//trim(FString_fromInteger(n))//")"
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Caso átomo - molécula lineal
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				if( this.nMolecules() == 2 .and. &
				( ( this.clusters(1).isLineal() .and. this.clusters(2).nAtoms() == 1 ) .or. &
				( this.clusters(2).isLineal() .and. this.clusters(1).nAtoms() == 1 ) ) ) then
					do i=1,this.nMolecules()
						call this.clusters(i).changeOrientation( force1D=.true. )
					end do
					
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Caso átomo - molécula
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				else if( this.nMolecules() == 2 .and. &
				( this.clusters(1).nAtoms() == 1 .or. this.clusters(2).nAtoms() == 1 ) ) then
					do i=1,this.nMolecules()
						call this.clusters(i).changeOrientation( force2D=.true. )
					end do
					
				else
					do i=1,this.nMolecules()
						call this.clusters(i).changeOrientation()
					end do
				end if
				
				if( GOptions_printLevel >= 2 ) then
					write(*,*) this.atomicOverlapping()
! 					call this.save("hola"//trim(FString_fromInteger(n))//".xyz")
				end if
				
				if( .not. this.atomicOverlapping() ) then
					check = .false.
					exit
				end if
			end do
			
			if( check ) then
				this.state = .false.
			end if
		end if
		
		call this.updateLogVtheta()
	end subroutine changeOrientationsFragmentsListBase
	
	!>
	!! @brief
	!!
	subroutine changeVibrationalEnergyFragmentsListBase( this )
		class(FragmentsListBase) :: this
		
		integer :: i
		
! 		if( this.forceInitializing ) then
! 			call this.initialGuessFragmentsListBase()
! 			return
! 		end if
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_subsection( "Random vibrational energies "//trim(this.label()), indent=2 )
		end if
		
		this.vibrationalEnergy_ = 0.0_8
		do i=1,this.nMolecules()
			call this.clusters(i).changeVibrationalEnergy()
			this.vibrationalEnergy_ = this.vibrationalEnergy_ + this.clusters(i).vibrationalEnergy_
		end do
			
		if( GOptions_printLevel >= 3 ) then
			call GOptions_paragraph( "Vibrational energy summary", indent=2 )
			
			if( GOptionsM3C_useZPECorrection ) then
				write(*,"(<GOptions_indentLength*2>X,A15,2A10)") "id", "Evib", "maxEvib"
				write(*,"(<GOptions_indentLength*2>X,A25,A10)") "eV", "eV"
				do i=1,this.nMolecules()
					write(*,"(<GOptions_indentLength*2>X,A15,2F10.5)") &
						trim(this.clusters(i).label()), this.clusters(i).vibrationalEnergy_/eV, this.clusters(i).maxEvib/eV
				end do
			else
				write(*,"(<GOptions_indentLength*2>X,A15,3A10)") "id", "Evib", "ZPE", "maxEvib"
				write(*,"(<GOptions_indentLength*2>X,A25,2A10)") "eV", "eV", "eV"
				do i=1,this.nMolecules()
					write(*,"(<GOptions_indentLength*2>X,A15,3F10.5)") &
						trim(this.clusters(i).label()), this.clusters(i).vibrationalEnergy_/eV, this.clusters(i).ZPE/eV, this.clusters(i).maxEvib/eV
				end do
			end if
			
			write(*,*) ""
			call GOptions_valueReport( "Evib", this.vibrationalEnergy_/eV, "eV", indent=2 )
		end if
		
		call this.updateLnWv()
	end subroutine changeVibrationalEnergyFragmentsListBase
	
	!>
	!! @brief Returns the electronic energy
	!!
	function electronicEnergy( this ) result( output )
		class(FragmentsListBase) :: this
		real(8) :: output
		
		integer :: i
		
		output = 0.0_8
		do i=1,this.nMolecules()
			output = output + this.clusters(i).electronicEnergy
		end do
	end function electronicEnergy
	
	!>
	!! @brief Return the internal energy
	!!
	function internalEnergy( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		real(8) :: output
		
		call GOptions_error( &
			"This method is abstract type so that you should always override ", &
			"FragmentsListBase.internalEnergy()" &
			)
	end function internalEnergy
	
	!>
	!! @brief Returns the total energy
	!!
	function totalEnergy( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		real(8) :: output
		
		call GOptions_error( &
			"This method is abstract type so that you should always override ", &
			"FragmentsListBase.totalEnergy()" &
			)
	end function totalEnergy
	
	!>
	!! @brief Returns number of vibrational degrees of freedom
	!!
	pure function ft( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		integer :: output
		
		output = this.ft_
	end function ft
	
	!>
	!! @brief
	!!
	pure function fl( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		integer :: output
		
		output = this.fl_
	end function fl
	
	!>
	!! @brief Number of rotational degrees of freedom
	!!
	pure function fr( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		integer :: output
		
		integer :: i
		
		output = 0
		do i=1,this.nMolecules()
			output = output + this.clusters( this.idSorted(i) ).fr()
		end do
	end function fr
	
	!>
	!! @brief
	!!
	function energyHistoryLine( this, prefix ) result( output )
		class(FragmentsListBase) :: this
		character(*), optional, intent(in) :: prefix
		type(String) :: output
		
		call GOptions_error( &
			"This method is abstract type so that you should always override ", &
			"FragmentsListBase.energyHistoryLine()" &
			)
	end function energyHistoryLine
	
	!>
	!! @brief
	!!
	function weightHistoryLine( this, prefix ) result( output )
		class(FragmentsListBase) :: this
		character(*), optional, intent(in) :: prefix
		type(String) :: output
		
		call GOptions_error( &
			"This method is abstract type so that you should always override ", &
			"FragmentsListBase.weightHistoryLine()" &
			)
	end function weightHistoryLine
	
	!>
	!! @brief
	!!
	function JHistoryLine( this, prefix ) result( output )
		class(FragmentsListBase) :: this
		type(String) :: output
		character(1), optional, intent(in) :: prefix
		
		character(1000) :: line
		character(1) :: prefixEff
		integer :: i
		
		prefixEff = ""
		if( present(prefix) ) prefixEff = prefix
		
! #define JVal(i) trim(adjustl(FString_fromReal(norm2(this.clusters(i).J_),"(F10.5)")))
#define JVal(i) trim(adjustl(FString_fromInteger(int(norm2(this.clusters(i).J_)),"(I10)")))
		write(line,"(1X,A1,2X,<this.nMolecules()>A)") trim(prefixEff), &
			( trim(this.clusters(this.idSorted(i)).label())//"#"//JVal(this.idSorted(i))//"  ", i=1,this.nMolecules() )
#undef JVal
		
		output = line
	end function JHistoryLine
	
	!>
	!! @brief
	!!
	function LHistoryLine( this, prefix ) result( output )
		class(FragmentsListBase) :: this
		type(String) :: output
		character(1), optional, intent(in) :: prefix
		
		character(1000) :: line
		character(1) :: prefixEff
		integer :: i
		
		prefixEff = ""
		if( present(prefix) ) prefixEff = prefix
		
! #define LVal trim(adjustl(FString_fromReal(norm2(this.L_),"(F10.5)")))
#define LVal trim(adjustl(FString_fromInteger(int(norm2(this.L_)),"(I10)")))
		write(line,"(1X,A1,2X,A)") trim(prefixEff), &
			trim(this.label())//"#"//LVal//"  "
#undef LVal
		
		output = line
	end function LHistoryLine
	
	!>
	!! @brief Total weight
	!!
	function LnW( this ) result( output )
		class(FragmentsListBase) :: this
		real(8) :: output
		
		call GOptions_error( &
			"This method is abstract type so that you should always override ", &
			"FragmentsListBase.LnW()" &
			)
	end function LnW
	
	!>
	!! @brief
	!!
	function LnWe( this ) result( output )
		class(FragmentsListBase) :: this
		real(8) :: output
		
		output = this.LnWe_
	end function LnWe
	
	!>
	!! @brief
	!!
	function LnWv( this ) result( output )
		class(FragmentsListBase) :: this
		real(8) :: output
		
		output = this.LnWv_
	end function LnWv
	
	!>
	!! @brief
	!!
	function LnWn( this ) result( output )
		class(FragmentsListBase) :: this
		real(8) :: output
		
		output = this.LnWn_
	end function LnWn

	!>
	!! @brief
	!!
	subroutine updateLnWe( this )
		class(FragmentsListBase) :: this
		
		integer :: i
		
		this.LnWe_ = 0.0_8
		do i=1,this.nMolecules()
			this.LnWe_ = this.LnWe_ + this.clusters(i).LnWe()
		end do
			
		if( GOptions_printLevel >= 3 ) then
			call GOptions_valueReport( "LnWe", this.LnWe_, indent=2 )
		end if
	end subroutine updateLnWe
	
	!>
	!! @brief Number of ordered partitions of the cluster
	!!
	subroutine updateLnWv( this )
		class(FragmentsListBase) :: this
		
		integer :: i
		
		this.LnWv_ = 0.0_8
		do i=1,this.nMolecules()
			this.LnWv_ = this.LnWv_ + this.clusters(i).LnWv()
		end do
			
		if( GOptions_printLevel >= 3 ) then
			call GOptions_valueReport( "LnWv", this.LnWv_, indent=2 )
		end if
	end subroutine updateLnWv
	
	!>
	!! @brief Number of ordered partitions of the cluster
	!!
	subroutine updateLnWn( this )
		class(FragmentsListBase) :: this
		
! 		this.LnWn_ = -log( real( Math_fact(this.nMolecules()), 8) ) &
! 				+log( real( Math_comb(this.nAtoms()-1,this.nMolecules()-1), 8) ) &
! 				+log( real( Math_comb(this.charge()+this.nMolecules()-1,this.nMolecules()-1), 8) )
		this.LnWn_ = this.logGFactor_
			
		if( GOptions_printLevel >= 3 ) then
			call GOptions_valueReport( "LnWn", this.LnWn_, indent=2 )
		end if
	end subroutine updateLnWn

	!>
	!! @brief
	!!
	subroutine showLnWComponents( this )
		class(FragmentsListBase) :: this
		
		call GOptions_error( &
			"This method is abstract type so that you should always override ", &
			"FragmentsListBase.showLnWComponents()" &
			)
	end subroutine showLnWComponents
	
	!>
	!! @brief
	!!
	subroutine updateLogGFactor( this )
		class(FragmentsListBase) :: this
		
		type(StringIntegerMap) :: counts
		integer :: i
		type(String) :: label
		
		class(StringIntegerMapIterator), pointer :: iter
		type(StringIntegerPair) :: pair
		
		call counts.init()
		
		do i=1,this.nMolecules()
			label = this.clusters(i).label( details=.true. )
			
			call counts.set( label, counts.at( label, defaultValue=0 )+1 )
		end do
		
		this.logGFactor_ = 0.0_8
		iter => counts.begin
		do while( associated(iter) )
			pair = counts.pair( iter )
			this.logGFactor_ = this.logGFactor_ - log( 1.0_8*Math_fact(pair.second) )
			
			iter => iter.next
		end do
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_valueReport( "-log(N!)", this.logGFactor_, indent=2 )
		end if
	end subroutine updateLogGFactor
	
	!>
	!! @brief
	!!
	subroutine updateLogVtheta( this )
		class(FragmentsListBase) :: this
		
		integer :: i
		
		this.logVtheta_ = 0.0_8
		do i=1,this.nMolecules()
			this.logVtheta_ = this.logVtheta_ + this.clusters(i).logVtheta_
		end do
			
		if( GOptions_printLevel >= 3 ) then
			write(*,*) ""
			call GOptions_valueReport( "logVtheta", this.logVtheta_, indent=2 )
		end if
	end subroutine updateLogVtheta
	
	!>
	!! @brief
	!!
	subroutine updateLogVfree( this, maxIter )
		class(FragmentsListBase) :: this
		integer, optional, intent(in) ::  maxIter
		
		integer :: method = 1
		integer :: effMaxIter
		
		real(8) :: Vsys ! Volume of the system
		real(8) :: probability ! Probability to generate a non-overlaping configuration
		integer :: nSuccesses
		type(RandomSampler) :: rs
		real(8), allocatable :: sample(:,:)
		real(8) :: rVec1(3), rVec2(3)
		integer :: n, i, j
		logical :: overlap
		
		effMaxIter = 100000
		if( present(maxIter) ) effMaxIter = maxIter
		
		this.state = .true.
		
! Testing overlap
#define OVERLAPPING(i,j) this.clusters(i).radius( type=GOptionsM3C_radiusType )+this.clusters(j).radius( type=GOptionsM3C_radiusType )-GOptionsM3C_overlappingRadius > norm2( rVec2-rVec1 )
! Convertion from spherical to cartesian coordinates
#define RVEC(i) [ sample(1,i)*sin(sample(2,i))*cos(sample(3,i)), sample(1,i)*sin(sample(2,i))*sin(sample(3,i)), sample(1,i)*cos(sample(2,i)) ]
		
		if( method == 1 ) then
			
			this.logVfree_ = 0.0_8
			do i=1,this.nMolecules()
				this.logVfree_ = this.logVfree_ &
					+ 3.0_8*log( (GOptionsM3C_systemRadius-this.clusters(i).radius( type=GOptionsM3C_radiusType ))/(2.0_8*Math_PI) )+log( 4.0_8*Math_PI/3.0_8 )
					
				if( GOptionsM3C_systemRadius < this.clusters(i).radius( type=GOptionsM3C_radiusType ) ) then
					write(*,*) "### ERROR ### FragmentsListBase.updateLogVfree: Rsys < Ri then Log( Rsys-Ri ) = NaN"
					stop
				end if
			end do
			
			if( this.nTrials_ /= -1 ) then
				this.logVfree_ = this.logVfree_ - log(real(this.nTrials_,8))
			end if
								
		else if( method == 2 ) then
		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Cálculo del volumen del sistema
			Vsys = 0.0_8
			
			! Dos clusters solo tienen un grado de libertad
			! por lo cual el volumen de integración es un línea
			if( this.nMolecules() == 2 ) then
				Vsys = 2.0_8*GOptionsM3C_systemRadius
			! Tres clusters se mueven sobre un plano por lo
			! cual el volumen de integración es un círculo
			else if( this.nMolecules() == 3 ) then
				Vsys = Math_PI*GOptionsM3C_systemRadius**2
			! Para el resto de casos el sistema se mueve dentro
			! de un volumen esférico
			else if( this.nMolecules() > 3 ) then
				Vsys = 4.0_8*Math_PI*GOptionsM3C_systemRadius**3/3.0_8
			end if
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			allocate( sample(3,this.nMolecules()) )
			
			call rs.init( nDim=3 )
			call rs.setRange( 1, [0.0_8,GOptionsM3C_systemRadius] )   ! r
			call rs.setRange( 2, [0.0_8,MATH_PI] )                 ! theta
			call rs.setRange( 3, [0.0_8,2.0_8*MATH_PI] )           ! phi
			
			nSuccesses = 0
			do n=1,effMaxIter
				
				if( GOptions_printLevel >= 3 ) write(*,*) "Calculating LogVfree avoiding fragments overlapping (step="//trim(FString_fromInteger(n))//")"
				
				call rs.uniform( sample )
				
				overlap = .false.
				
				do i=1,this.nMolecules()-1
					do j=i+1,this.nMolecules()
						
						rVec1 = RVEC(i)
						rVec2 = RVEC(j)
							
						overlap = OVERLAPPING(i,j) .or. norm2(rVec1) > GOptionsM3C_systemRadius .or. norm2(rVec2) > GOptionsM3C_systemRadius
						
						if( overlap ) exit
					end do
					
					if( overlap ) exit
				end do
				
				if( .not. overlap ) then
					nSuccesses = nSuccesses + 1
				end if
				
			end do
			
			deallocate( sample )
		
			if( nSuccesses == 0 ) then
				this.state = .false.
! 				call GOptions_error( &
! 					"Maximum number of iterations reached"//" (n = "//trim(FString_fromInteger(effMaxIter))//")", &
! 					"FragmentsListBase.updateLogVfree(method=2)", &
! 					"Consider to increase GOptions:systemRadius" &
! 					)
			end if
			
			! Probabilidad de generar un configuración no solapante
			probability = real(nSuccesses,8)/real(effMaxIter,8)
			
			this.logVfree_ = this.nMolecules()*log(probability*Vsys)
			
			! Esto no estoy seguro si es de aquí, pero fijo no es del método 1
			this.logVfree_ = this.logVfree_ - this.ft_*log( 2.0_8*Math_PI )
		end if
		
#undef OVERLAPPING
#undef RVEC
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_valueReport( "Rsys", GOptionsM3C_systemRadius/angs, "A", indent=2 )
			
			if( method == 2 ) then
				call GOptions_valueReport( "Vsys", Vsys/angs**3, "A^3", indent=2 )
				call GOptions_valueReport( "Ps", probability, indent=2 )
			end if
			
			call GOptions_valueReport( "logVfree", this.logVfree_, "nlog(a.u.^3)", indent=2 )
		end if
		
	end subroutine updateLogVfree
	
	!>
	!! @brief Builds the intertia tensor and updates the number of
	!!        tranlational degrees of freedom (ft)
	!!
	subroutine updateInertiaTensor( this )
		class(FragmentsListBase) :: this
		
		type(Matrix) :: Im
		real(8), allocatable :: eValues(:)
		integer :: i
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Construyendo los ejes de inercia por diagonalización
		! hay una fase de pi/2 que no se puede controlar. En
		! la clase Molecule de scift, esto se ha solucionado
		! rotando gradualmente los ejes de inercia, sin embargo
		! aquí no hay control sobre los ejes, así que creo que
		! no hay otra solución para este caso.
		call this.buildInertiaTensor( Im )
		call Im.eigen( eValues=eValues, eVecs=this.inertiaAxes )
		this.diagInertiaTensor = eValues(1:3)
		deallocate( eValues )
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		this.ft_ = 0
		if( this.nMolecules() == 1 ) then
			this.ft_ = 0
			this.fl_ = 0
		else if( abs( this.diagInertiaTensor(1) ) < 1e-3 ) then ! Este valor detecta a partir 179.988 deg.
			this.ft_ = 3*this.nMolecules()-5
			this.fl_ = 2
		else
			this.ft_ = 3*this.nMolecules()-6
			this.fl_ = 3
		end if
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_paragraph( "Inertia tensor", indent=2 )
			
			call GOptions_valueReport( "ft", this.ft_, indent=2 )
			call GOptions_valueReport( "fl", this.fl_, indent=2 )
			write(*,*) ""
			write(*,"(<GOptions_indentLength*2>X,3A20)") "Ixx", "Iyy", "Izz"
			write(*,"(<GOptions_indentLength*2>X,3A20)") "amu*A^2", "amu*A^2", "amu*A^2"
			write(*,"(<GOptions_indentLength*2>X,3F20.5)") this.diagInertiaTensor/amu/angs**2
			write(*,*) ""
		end if
	end subroutine updateInertiaTensor
	
	!>
	!! @brief
	!!
	subroutine updateIntermolecularPotential( this )
		class(FragmentsListBase) :: this
		
		integer :: i, j
		real(8) :: Rij, rvij, twoFragContrib, oneFragContrib
		integer :: idProd
		real(8) :: rBuffer
		
		this.intermolEnergy_ = 0.0_8
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_paragraph( "Intermolecular energy contributions", indent=2 )
			
			write(*,"(<GOptions_indentLength*2>X,A)") "One body contributions"
			
			write(*,"(<GOptions_indentLength*2>X,A30,3X,A15)") "id", "E"
			write(*,"(<GOptions_indentLength*2>X,33X,A15)") "eV"
		end if
		
		oneFragContrib = 0.0_8
		do i=1,this.nMolecules()
			oneFragContrib = oneFragContrib + this.clusters(i).electronicEnergy
			
			if( GOptions_printLevel >= 3 ) then
				write(*,"(<GOptions_indentLength*2>X,A30,3X,F15.5)") &
					trim(this.clusters(i).label()), &
					this.clusters(i).electronicEnergy/eV
			end if
		end do
		
! 		if( this.nMolecules() < 2 ) ... @todo Creo que el siguiente bloque va solo si esto se satisface

		if( GOptions_printLevel >= 3 ) then
			write(*,"(<GOptions_indentLength*2>X,A)") "Two body contributions"
			write(*,"(<GOptions_indentLength*2>X,A30,3X,A10,3X,A10,3X,A15)") "id1--id2", "r12", "E"
			write(*,"(<GOptions_indentLength*2>X,A43,3X,A10,3X,A15)") "A", "eV"
		end if
		
		this.currentProducts = 0
		twoFragContrib = 0.0_8
		do i=1,this.nMolecules()-1
			do j=i+1,this.nMolecules()
				
				! @todo Esto hay que pensarlo bien, pero tambien habría que cambiarlo en
				!       la generación de las geometrías
! 				rvij = norm2( this.clusters(i).centerOfMass()-this.clusters(j).centerOfMass() )
				rvij = norm2( this.clusters(i).center()-this.clusters(j).center() )
				
				! @todo Esto solo ocurren en algunos casos concretos y no se porque. Esto ya deberia estar filtrado desde la generación de geometría. Cuando ocurre hay desviaciones de solo ~0.05A
! 				if( this.clusters(i).radius( type=GOptionsM3C_radiusType )+this.clusters(j).radius( type=GOptionsM3C_radiusType )-GOptionsM3C_overlappingRadius > rvij ) then
! 					call GOptions_error( &
! 						 "Overlapping configuration found", &
! 						 "FragmentsListBase.updateIntermolecularPotential()", &
! 						 "( R1+R2-S = "//&
! 						 trim(FString_fromReal((this.clusters(i).radius( type=GOptionsM3C_radiusType )+this.clusters(j).radius( type=GOptionsM3C_radiusType )-GOptionsM3C_overlappingRadius)/angs,"(F5.3)"))//" A ) > " &
! 						 //"( rij = "//trim(FString_fromReal(rvij/angs,"(F5.3)"))//" A ) " &
! 						 //trim(this.clusters(i).label())//" --- "//trim(this.clusters(j).label()) &
! 					)
! 				end if
				
				rBuffer = FragmentsDB_instance.potential( this.clusters(i).id, this.clusters(j).id, rvij )
				
				if( GOptions_printLevel >= 3 ) then
					write(*,"(<GOptions_indentLength*2>X,A30,3X,F10.5,3X,F15.5)") &
						trim(this.clusters(i).label())//"--"//trim(this.clusters(j).label()), &
						rvij/angs, rBuffer/eV
				end if
				
				twoFragContrib = twoFragContrib + rBuffer
			end do
		end do
		
		this.intermolEnergy_ = oneFragContrib + twoFragContrib - FragmentsDB_instance.energyReference()
		
		if( GOptions_printLevel >= 3 ) then
			write(*,*) ""
			write(*,"(<GOptions_indentLength*2>X,A30,3X,F15.5)") "V^(1)", oneFragContrib/eV
			write(*,"(<GOptions_indentLength*2>X,A30,16X,F15.5)") "V^(2)", twoFragContrib/eV
			write(*,"(<GOptions_indentLength*2>X,A30,16X,F15.5)") "E^0", FragmentsDB_instance.energyReference()/eV
			write(*,*) ""
			call GOptions_valueReport( "V(r)", this.intermolEnergy_/eV, "eV", indent=2 )
		end if
	end subroutine updateIntermolecularPotential

! 	!>
! 	!! @brief
! 	!!
! 	subroutine updateIntermolecularPotential( this )
! 		class(FragmentsListBase) :: this
! 		
! 		integer :: i, j
! 		real(8) :: Rij, rvij, twoFragContrib
! 		integer :: idProd
! 		
! 		this.intermolEnergy_ = 0.0_8
! 		
! 		if( this.nMolecules() < 2 ) return
! 		
! 		if( GOptions_printLevel >= 3 ) then
! 			call GOptions_paragraph( "Intermolecular energy contributions", indent=2 )
! 			
! 			write(*,"(<GOptions_indentLength*2>X,A30,3X,A10,3X,A10,3X,A10)") "id1--id2", "r12", "E"
! 			write(*,"(<GOptions_indentLength*2>X,A43,3X,A10,3X,A10)") "A", "eV"
! 		end if
! 		
! 		do i=1,this.nMolecules()-1
! 			do j=i+1,this.nMolecules()
! 				
! 				rvij = norm2( this.clusters(i).center()-this.clusters(j).center() )
! 				
! 				if( this.clusters(i).radius( type=GOptionsM3C_radiusType )+this.clusters(j).radius( type=GOptionsM3C_radiusType )-GOptions_overlappingRadius > rvij ) then
! 					call GOptions_error( &
! 						 "Overlapping configuration found", &
! 						 "Fragment.updateIntermolecularPotential()", &
! 						 "( R1+R2 = "//&
! 						 trim(FString_fromReal((this.clusters(i).radius( type=GOptionsM3C_radiusType )+this.clusters(j).radius( type=GOptionsM3C_radiusType ))/angs,"(F5.3)"))//" A ) > " &
! 						 //"( rij = "//trim(FString_fromReal(rvij/angs,"(F5.3)"))//" A ) " &
! 						 //trim(this.clusters(i).label())//" --- "//trim(this.clusters(j).label()) &
! 					)
! 				end if
! 				
! 				twoFragContrib = FragmentsDB_instance.potential( this.clusters(i).id, this.clusters(j).id, rvij )
! 				
! 				if( GOptions_printLevel >= 3 ) then
! 					write(*,"(<GOptions_indentLength*2>X,A30,3X,F10.5,3X,F10.5)") &
! 						trim(this.clusters(i).label())//"--"//trim(this.clusters(j).label()), &
! 						rvij/angs, twoFragContrib/eV
! 				end if
! 				
! 				this.intermolEnergy_ = this.intermolEnergy_ + twoFragContrib
! 			end do
! 		end do
! 		
! 		if( GOptions_printLevel >= 3 ) then
! 			write(*,*) ""
! 			call GOptions_valueReport( "V(r)", this.intermolEnergy_/eV, "eV", indent=2 )
! 		end if
! 	end subroutine updateIntermolecularPotential
	
	!>
	!! @brief Returns the total energy
	!!
	subroutine updateKineticEnergy( this )
		class(FragmentsListBase) :: this
		
		this.kineticEnergy_ = this.reactorEnergy() - this.internalEnergy()
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_valueReport( "updated K", this.kineticEnergy()/eV, indent=2 )
		end if
	end subroutine updateKineticEnergy
	
	!>
	!! @brief Inverse of the temperature
	!!
	function iTemperature( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		real(8) :: output
		
		call GOptions_error( &
			"This method is abstract type so that you should always override ", &
			"FragmentsListBase.iTemperature()" &
			)
	end function iTemperature
	
	!>
	!! @brief
	!!
	function centerOfMass( this ) result( cm )
		class(FragmentsListBase), intent(in) :: this
		real(8) :: cm(3)
		
		integer :: i
		
		cm = 0.0_8
		
		do i=1,this.nMolecules()
			cm = cm + this.clusters(i).mass()*this.clusters(i).centerOfMass()
		end do
		
		cm = cm/this.mass()
	end function centerOfMass
	
	!>
	!! @brief Builds the inertia tensor with the current geometry
	!! around of the origin of the coordinates system by considering
	!! each cluster as a punctual mass
	!!
	subroutine buildInertiaTensor( this, Im, center )
		class(FragmentsListBase), intent(in) :: this
		type(Matrix), intent(out) :: Im
		real(8), optional, intent(in) :: center(3)
		
		real(8) :: effCenter(3)
		
		integer :: i
		real(8) :: iCenterOfMass(3)
		real(8), allocatable :: X(:), Y(:), Z(:), m(:)
		
		effCenter = 0.0_8
		if( present(center) ) effCenter = center
		
		if ( this.nMolecules() == 1 ) then
			Im = this.clusters(1).inertiaTensor( center )
			return
		end if
		
		allocate( X(this.nMolecules()) )
		allocate( Y(this.nMolecules()) )
		allocate( Z(this.nMolecules()) )
		allocate( m(this.nMolecules()) )
		
		do i=1,this.nMolecules()
			iCenterOfMass = this.clusters(i).centerOfMass()
			X(i) = iCenterOfMass(1) - effCenter(1)
			Y(i) = iCenterOfMass(2) - effCenter(2)
			Z(i) = iCenterOfMass(3) - effCenter(3)
			m(i) = this.clusters(i).mass()
		end do
		
		call Im.init(3,3)
		
		call Im.set( 1, 1,  sum( m*(Y**2+Z**2) ) )
		call Im.set( 1, 2, -sum( m*X*Y ) )
		call Im.set( 1, 3, -sum( m*X*Z ) )
		call Im.set( 2, 1, -sum( m*Y*X ) )
		call Im.set( 2, 2,  sum( m*(X**2+Z**2) ) )
		call Im.set( 2, 3, -sum( m*Y*Z ) )
		call Im.set( 3, 1, -sum( m*Z*X ) )
		call Im.set( 3, 2, -sum( m*Z*Y ) )
		call Im.set( 3, 3,  sum( m*(X**2+Y**2) ) )
		
		deallocate( X )
		deallocate( Y )
		deallocate( Z )
		deallocate( m )
	end subroutine buildInertiaTensor
	
! 	!>
! 	!! @brief Orienta el FragmentsListBase en dirección de los ejes seleccionados "axes"
! 	!!        asumiendo que el centro de masas ya ha sido previamente seleccioando
! 	!!
! 	subroutine orient( this, axes )
! 		class(FragmentsListBase) :: this
! 		type(Matrix), intent(in) :: axes
! 		
! 		type(Matrix) :: Rot, r, Im
! 		real(8), allocatable :: eValues(:)
! 		real(8) :: rThetaPhi(3)
! 		integer :: i
! 		
! 		rThetaPhi = Math_cart2Spher( axes.data(:,3) )
! 		Rot = SpecialMatrix_rotation( rThetaPhi(3), rThetaPhi(2), 0.0_8 )
! 		
! 		call r.columnVector( 3, values=axes.data(:,1) )
! 		r = Rot*r
! 		rThetaPhi = Math_cart2Spher( r.data(:,1) )
! 		Rot = SpecialMatrix_zRotation( rThetaPhi(3) )*Rot
! 		
! 		do i=1,this.nMolecules()
! 			call r.columnVector( 3, values=this.clusters(i).center() )
! 			r = Rot*r
! 			call this.clusters(i).setCenter( r.data(:,1) )
! 		end do
! 		
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 		! Construyendo los ejes de inercia por diagonalización
! 		! hay una fase de pi/2 que no se puede controlar. En
! 		! la clase Molecule de scift, esto se ha solucionado
! 		! rotando gradualmente los ejes de inercia, sin embargo
! 		! aquí no hay control sobre los ejes, así que creo que
! 		! no hay otra solución para este caso.
! 		!
! 		! Se recalculan los ejes principales y el tensor de inercia
! 		! de tal forma que si "axes" son los ejes de inercia, los
! 		! nuevos ejes deberían ser la matriz unitaria
! 		call this.buildInertiaTensor( Im )
! 		call Im.eigen( eValues=eValues, eVecs=this.inertiaAxes )
! 		this.diagInertiaTensor = eValues(1:3)
! 		deallocate( eValues )
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	end subroutine orient

	!>
	!! @brief Orient molecule
	!!
	!! Orient molecule such that origin is centre of mass, and axes are eigenvectors of inertia tensor
	!! primary axis   Z
	!! primary plane YZ
	!! 
	!!
	subroutine orient( this )
		class(FragmentsListBase) :: this
		
		real(8) :: centerOfMass(3)
		type(Matrix) :: Im, Vm, Rot, r
		real(8), allocatable :: Imd(:)
		real(8) :: rThetaPhi(3)
		integer :: i
		
		real(8), allocatable :: myDiagTensor(:)
		
		! Este pequeño bloque debería hacerse solo una vez
		! En Molecule de scift, esta implementado con un flag que se
		! activa una vez se ha hecho la primera vez, pero acá no lo
		! he hecho por pereza
		centerOfMass = 0.0_8
		do i=1,this.nMolecules()
			centerOfMass = centerOfMass + &
				this.clusters(i).mass()*this.clusters(i).center()
		end do
		centerOfMass = centerOfMass/this.mass()
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! 1) It choose the origin as the center of mass
! 		centerOfMass = this.centerOfMass()
		do i=1,this.nMolecules()
			call this.clusters(i).setCenter( this.clusters(i).center() - centerOfMass )
		end do
		
		! Esto es para comprobar
		centerOfMass = 0.0_8
		do i=1,this.nMolecules()
			centerOfMass = centerOfMass + &
				this.clusters(i).mass()*this.clusters(i).center()
		end do
		centerOfMass = centerOfMass/this.mass()
		
! 		if( this.nMolecules() == 1 ) then
! 			call this.clusters(1).orient()
! 			
! 			this.diagInertiaTensor = [ this.clusters(1).diagInertiaTensor.get(1,1), &
! 									   this.clusters(1).diagInertiaTensor.get(2,2), &
! 									   this.clusters(1).diagInertiaTensor.get(3,3) ]
! 			
! 			this.inertiaAxes = this.clusters(1).inertiaAxes()
! 			
! 			this.ft_ = 0
! 			this.fr_ = 0
! 			this.fv_ = this.clusters(1).fv()
! 			
! 			return
! 		end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! 2) It builds the inertia tensor which is diagonalizes
		!    in order to get the principal axes
		call this.buildInertiaTensor( Im )
		call Im.eigen( eVecs=Vm, eValues=myDiagTensor )
		
		this.diagInertiaTensor = myDiagTensor(1:3)
		deallocate(myDiagTensor)
		
		rThetaPhi = Math_cart2Spher( Vm.data(:,3) )
		Rot = SpecialMatrix_rotation( rThetaPhi(3), rThetaPhi(2), 0.0_8 )
		
		call r.columnVector( 3, values=Vm.data(:,1) )
		r = Rot*r
		rThetaPhi = Math_cart2Spher( r.data(:,1) )
		Rot = SpecialMatrix_zRotation( rThetaPhi(3) )*Rot
		
		do i=1,this.nMolecules()
			call r.columnVector( 3, values=this.clusters(i).center() )
			r = Rot*r
			call this.clusters(i).setCenter( r.data(:,1) )
		end do
		
		this.inertiaAxes.data(:,1) = [ 1.0_8, 0.0_8, 0.0_8 ]
		this.inertiaAxes.data(:,2) = [ 0.0_8, 1.0_8, 0.0_8 ]
		this.inertiaAxes.data(:,3) = [ 0.0_8, 0.0_8, 1.0_8 ]
		
! 		this.ft_ = 0
! 		this.fr_ = 0
! 		if( this.nMolecules() == 1 ) then
! 			this.ft_ = 0
! 			this.fr_ = 0
! 			
! 			this.diagInertiaTensor = [ this.clusters(1).diagInertiaTensor.get(1,1), &
! 									   this.clusters(1).diagInertiaTensor.get(2,2), &
! 									   this.clusters(1).diagInertiaTensor.get(3,3) ]
! 		else
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			! Construyendo los ejes de inercia por diagonalización
! 			! hay una fase de pi/2 que no se puede controlar. En
! 			! la clase Molecule de scift, esto se ha solucionado
! 			! rotando gradualmente los ejes de inercia, sin embargo
! 			! aquí no hay control sobre los ejes, así que creo que
! 			! no hay otra solución para este caso.
! ! 			call this.buildInertiaTensor( Im )
! ! 			call Im.eigen( eValues=eValues, eVecs=this.inertiaAxes )
! ! 			this.diagInertiaTensor = eValues(1:3)
! ! 			deallocate( eValues )
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			
! 			if( abs( this.diagInertiaTensor(1) ) < 1e-3 ) then ! Este valor detecta a partir 179.988 deg.
! 				this.ft_ = 3*this.nMolecules()-5
! 				this.fr_ = 2
! 			else
! 				this.ft_ = 3*this.nMolecules()-6
! 				this.fr_ = 3
! 			end if
! 		end if
! 		
! 		this.fv_ = 0
! 		do i=1,this.nMolecules()
! 			this.fv_ = this.fv_ + this.clusters(i).fv()
! 		end do
		
! 		if( GOptions_printLevel >= 3 ) then
! 			call GOptions_paragraph( "Inertia tensor", indent=2 )
! 			
! 			call GOptions_valueReport( "ft", this.ft_, indent=2 )
! 			call GOptions_valueReport( "fr", this.fr_, indent=2 )
! 			call GOptions_valueReport( "fv", this.fv_, indent=2 )
! 			write(*,*) ""
! 			write(*,"(<GOptions_indentLength*2>X,3A20)") "Ixx", "Iyy", "Izz"
! 			write(*,"(<GOptions_indentLength*2>X,3A20)") "amu*A^2", "amu*A^2", "amu*A^2"
! 			write(*,"(<GOptions_indentLength*2>X,3F20.5)") this.diagInertiaTensor/amu/angs**2
! 			write(*,*) ""
! 		end if
	end subroutine orient
	
	!>
	!! @brief Return the radius of the system in atomic units
	!!        defined as half the largest distance between
	!!        two fragments
	!!
	real(8) function radius( this ) result ( output )
		class(FragmentsListBase), intent(in) :: this
		
		real(8) :: rij
		real(8) :: cRadius1, cRadius2
		integer :: i, j
		integer :: si, sj
		
		output = 0.0_8
		
		if( this.nFragments() > 1 ) then
			
			si = -1
			sj = -1
			do i=1,this.nFragments()-1
				do j=i+1,this.nFragments()
					rij = norm2( this.clusters(j).center()-this.clusters(i).center() )
					
					if( rij > output ) then
						output = rij
						si = i
						sj = j
					end if
				end do
			end do
			
			cRadius1 = this.clusters(si).radius( type=GOptionsM3C_radiusType )
			cRadius2 = this.clusters(sj).radius( type=GOptionsM3C_radiusType )
			
			output = output/2.0_8 + max( cRadius1, cRadius2 )
		else
			output = this.clusters(1).radius( type=GOptionsM3C_radiusType )
		end if
	end function radius
	
	!>
	!! @brief
	!!
	pure function nFragments( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		integer :: output
		
		output = size(this.clusters)
	end function nFragments
	
	!>
	!! @brief
	!!
	function spinRange( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		real(8) :: output(2)
		
		real(8) :: S
		real(8) :: minS, maxS
		logical :: firstTime
		
		integer :: i
		
		firstTime = .true.
		minS = 0.0_8
		maxS = 0.0_8
		do i=1,this.nMolecules()
			S = (this.clusters(i).multiplicity-1.0_8)/2.0_8
			
			if( S < 0.0_8 ) cycle
			
			maxS = maxS + S
			minS = minS + merge( 1.0_8, -1.0_8, firstTime )*S
			firstTime = .false.
			
! 			write(6,*) "Spin = ", S, maxS, minS
		end do
		minS = merge( 0.0_8, abs(minS), minS <= 0.0_8 )
! 		write(6,*) "check Spin = ", S, maxS, minS
		
		output = [ minS, maxS ]
	end function spinRange
	
	!>
	!! @brief
	!!
	function spinAvailable( this ) result( output )
		class(FragmentsListBase), intent(in) :: this
		type(RealList) :: output
		
		real(8) :: S, Si, Sj
		integer :: i, j
		
		call output.init()
		
		if( this.nMolecules() == 1 ) then
			S = (this.clusters(1).multiplicity-1.0_8)/2.0_8
			call output.append( S )
		else
			do i=1,this.nMolecules()-1
				Si = (this.clusters(i).multiplicity-1.0_8)/2.0_8
				
				if( Si < 0.0_8 ) Si=0.0_8
				
				do j=i+1,this.nMolecules()
					Sj = (this.clusters(j).multiplicity-1.0_8)/2.0_8
					
					if( Sj < 0.0_8 ) Sj=0.0_8
					
					S = abs(Si-Sj)
					do while( int(2.0*S) <= int(2.0*(Si+Sj)) )
						call output.append( S )
						
						S = S + 1.0_8
					end do
				end do
			end do
		end if
		
		if( output.size() == 0 ) then
			write(*,*) "### ERROR ### FragmentsListBase.spinAvailable.size() == 0", this.nMolecules(), S, Si, Sj
		end if
	end function spinAvailable
	
	!>
	!! @brief
	!!
	subroutine executeRadiusOptimization( this, iParser )
		class(FragmentsListBase) :: this
		type(BlocksIFileParser), intent(in) :: iParser
		
		type(String) :: strReactives
		character(20), allocatable :: reactiveTokens(:)
		real(8) :: r, rMin, rMax, rStep
		real(8), allocatable :: logIData(:), nTrials_Data(:)
		real(8) :: averLogI, averNTrials, stdevLogI, stdevNTrials
		integer :: nExp
		type(String) :: rSysOptOFile
		type(OFStream) :: oFile
		
		type(String) :: sBuffer
		integer :: iBuffer
		integer :: i, k
		
		if( .not. iParser.isThereBlock( "SYSROPT" ) ) then
			return
		end if
		
		sBuffer = iParser.getString( "SYSROPT:reactives" )
		
		call sBuffer.split( reactiveTokens, ":" )
		
		if( reactiveTokens(1) == "file" ) then
			call this.loadXYZ( reactiveTokens(2) )
		else
			strReactives = FragmentsDB_instance.extendFragmentsListName( sBuffer.fstr )
			call strReactives.split( reactiveTokens, "+" )
			
			call this.init( size(reactiveTokens) )
			do i=1,size(reactiveTokens)
				iBuffer = FragmentsDB_instance.getIdFromName( reactiveTokens(i) )
				call this.set( i, FragmentsDB_instance.clusters(iBuffer) )
			end do
		end if
		
		if( allocated(reactiveTokens) ) deallocate( reactiveTokens )
		
		rMin = iParser.getReal( "SYSROPT:rMin", def=2.0_8 )*angs
		rMax = iParser.getReal( "SYSROPT:rMax", def=10.0_8 )*angs
		rStep = iParser.getReal( "SYSROPT:rStep", def=2.0_8 )*angs
		nExp = iParser.getInteger( "SYSROPT:nExp", def=10 )
		rSysOptOFile = iParser.getString( "SYSROPT:outputFile", def="#@NONE@#" )
		
		if( trim(rSysOptOFile.fstr) /= "#@NONE@#" ) then
			call oFile.init( rSysOptOFile.fstr )
		end if
		
		allocate( nTrials_Data(nExp) )
		allocate( logIData(nExp) )
		
		write(*,*) ""
		write(*,"(A)") "----------------------------------------------------------"
		write(*,"(A)") " SYSTEM RADIUS OPTIMIZATION"
		write(*,*) ""
		write(*,"(15X,3A10)") "R", "nTrials_", "log(I1*I2*I3)"
		write(*,"(A15,3A10)") "", "A", "", ""

		GOptionsM3C_systemRadius = rMax
		do while( GOptionsM3C_systemRadius >= rMin )
			
			do i=1,nExp
				this.forceRandomCenters = .true.
				call this.changeGeometryFragmentsListBase()
				
				logIData(i) = 0.0_8
				do k=1,3
					if( abs(this.diagInertiaTensor(k)) > 1e4*GOptions_zero ) then
						logIData(i) = logIData(i) + log( this.diagInertiaTensor(k) )
					end if
				end do
				
				nTrials_Data(i) = real(this.nTrials_,8)
			end do
			
			averNTrials = sum(nTrials_Data)/nExp
			stdevNTrials = sqrt( sum((nTrials_Data-averNTrials)**2)/real(nExp-1,8) )
			
			averLogI = sum(logIData)/nExp
			stdevNTrials = sqrt( sum((logIData-averLogI)**2)/real(nExp-1,8) )
			
			if( trim(rSysOptOFile.fstr) /= "#@NONE@#" ) then
				write(oFile.unit,"(15X,F10.5,5X,2F10.5,5X,2F15.5)") GOptionsM3C_systemRadius/angs, &
					averNTrials, stdevNTrials, averLogI, stdevNTrials
			end if
				
			write(*,"(15X,F10.5,5X,2F10.5,5X,2F15.5)") GOptionsM3C_systemRadius/angs, &
				averNTrials, stdevNTrials, averLogI, stdevNTrials
			
			GOptionsM3C_systemRadius = GOptionsM3C_systemRadius - rStep
		end do
		write(*,"(A)") "----------------------------------------------------------"
		
		if( trim(rSysOptOFile.fstr) /= "#@NONE@#" ) then
			call oFile.close()
		end if
		
		deallocate( nTrials_Data )
		deallocate( logIData )
		
		stop
	end subroutine executeRadiusOptimization
		
end module FragmentsListBase_
	

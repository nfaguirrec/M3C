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
module Fragment_
	use GOptions_
	use Math_
	use IOStream_
	use UnitsConverter_
	use String_
	use Matrix_
	use SpecialMatrix_
	use RandomUtils_
	use Molecule_
	use StringList_
	use RealVector_
	
	use GOptionsM3C_
	
	implicit none
	private
	
	public :: &
		Fragment_getLabelFragments, &
		Fragment_test
	
	type, public, extends( Molecule ):: Fragment
		
		integer :: id
		integer :: charge
		integer :: multiplicity
		integer :: degL
		integer :: sigmaSym
		real(8) :: electronicEnergy
		real(8), allocatable :: vibFrequencies(:)
		logical :: isTransitionState
		real(8) :: maxEvib
		real(8) :: ZPE
		real(8) :: maxJ
		character(:), allocatable, private :: fileName
		
		logical :: frozen
		
		real(8) :: vibrationalEnergy_
		real(8) :: J_(3)
		
		real(8) :: logVtheta_
		real(8) :: logVj_
		
		real(8), private :: LnWe_
		real(8), private :: LnWv_
		
		character(100), private :: label_
		character(100), private :: dlabel_
		logical, private :: testLabel_
		
		contains
			generic :: init => initDefault, fromMassTableRow
			generic :: assignment(=) => copyFragment
			
			procedure :: initDefault
			procedure :: fromMassTableRow
			procedure :: copyFragment
			final :: destroyFragment
			
			procedure, private :: updateLabel
			procedure :: label
			
! 			procedure, private :: loadRXYZ
			procedure :: loadRXYZ
			procedure, private :: loadMOLDEN
			
			procedure :: changeVibrationalEnergy
			procedure :: changeOrientation
! 			procedure :: changeAngularMomentum ! No lo he borrado, porque podría ser importante para el reactor JL
			
			procedure :: LnW
			procedure :: LnWe
			procedure :: LnWv
			procedure, private :: updateLnWe
			procedure, private :: updateLnWv
			
			procedure, private :: updateLogVtheta
			procedure, private :: updateLogVJ
			
			procedure :: showLnWComponents
			
			procedure :: spinAvailable
! 			procedure :: rotSymNumber
	end type Fragment
	
	contains
	
	!>
	!! @brief Constructor
	!!
	subroutine initDefault( this )
		class(Fragment) :: this 
		
		this.id = -1
		this.charge = 0
		this.multiplicity = 1
		this.degL = 1
		this.sigmaSym = 1
		this.electronicEnergy = 0.0_8
		
		if( allocated(this.vibFrequencies) ) deallocate( this.vibFrequencies )
		this.isTransitionState = .false.
		this.maxEvib = 0.0_8
		this.ZPE = 0.0_8
		this.maxJ = 0.0_8
		if( allocated(this.fileName) ) deallocate(this.fileName)
		
		this.frozen = .false.
		
		this.vibrationalEnergy_ = 0.0_8
		this.J_(3) = 0.0_8
		this.logVtheta_ = 0.0_8
		this.logVJ_ = 0.0_8
		
		this.LnWe_ = 0.0_8
		this.LnWv_ = 0.0_8
		
		this.label_ = ""
		this.dlabel_ = ""
		this.testLabel_ = .false.
	end subroutine initDefault
	
	!>
	!! @brief Constructor
	!!
	subroutine fromMassTableRow( this, strRow, id, store )
		class(Fragment) :: this 
		character(*), intent(in) :: strRow
		integer, optional :: id
		character(*), optional :: store
		
		real(8) :: rBuffer
		character(100), allocatable :: tokens(:)
		character(100) :: extension
		
		type(String) :: effStore
		
		effStore = "."
		if( present(store) ) effStore = store
		
		call this.initDefault()
		
		if( GOptions_printLevel >= 4 ) then
			call GOptions_paragraph( "Loading Fragment data from mass table row", indent=3 )
		end if
		
		if( allocated(this.fileName) ) deallocate( this.fileName )
		
		this.id = -1
		if( present(id) ) this.id = id
		
		call FString_split( strRow, tokens, " " )
		
		if( size(tokens) >= 7 ) then
			this.name = trim(adjustl(tokens(1)))
			this.charge = FString_toInteger( tokens(2) )
			this.multiplicity = FString_toInteger( tokens(3) )
			this.degL = FString_toInteger( tokens(4) )
			this.sigmaSym = FString_toInteger( tokens(5) )
			this.fileName = trim(effStore.fstr)//"/"//trim(tokens(6))
			this.electronicEnergy = FString_toReal( tokens(7) )*eV
		else
			call GOptions_error( &
				"Mass table row incomplete, this row should have at least 7 parameters", &
				"Fragment.fromMassTableRow()", &
				trim(strRow) &
			)
		end if
		
		deallocate( tokens )
		
		! Las frecuencias en cm-1 son leidas del fichero de geometría
		! enriquecido y son convertidas a a.u.
		call FString_split( this.fileName, tokens, "." )
		extension = ""
		if( size( tokens ) > 1 ) extension = trim(tokens(size(tokens)))
		deallocate( tokens )
		
		if( trim(extension) == "xyz" ) then
			call this.loadRXYZ( this.fileName, loadName=.false. )
		else if( trim(extension) == "rxyz" ) then
			call this.loadRXYZ( this.fileName, loadName=.false. )
		else if( trim(extension) == "molden" ) then
			call this.loadMOLDEN( this.fileName, loadName=.false. )
		else
			call GOptions_error( &
				"Unknown format in Geometry-Frequency file "//trim(this.fileName), &
				"SMolecule.fromMassTableRow()", &
				trim(strRow) &
			)
		end if
		
		if( any( this.vibFrequencies < 0.0_8 ) ) then
			this.isTransitionState = .true.
		end if
		
		! Zero point energy
		this.ZPE = 0.0_8
		if( this.nAtoms() > 1 ) then
			this.ZPE = sum( merge(this.vibFrequencies, 0.0_8,this.vibFrequencies>0.0_8) )/2.0_8
			
			if( GOptions_printLevel >= 4 ) then
				call GOptions_valueReport( "ZPE", this.ZPE, "eV", indent=3 )
			end if
		end if
			
		if( GOptionsM3C_useZPECorrection ) then
			this.electronicEnergy = this.electronicEnergy + this.ZPE
		end if
			
! 		if( GOptions_printLevel >= 4 ) then
! 			call GOptions_valueReport( "Eelec", this.electronicEnergy/eV, "eV", indent=3 )
! 		end if
				
		! The maxEvib and maxJ values, are updated in ClusterDB
		this.maxEvib = 0.0_8
		this.maxJ = 0
				
		this.LnWe_ = 0.0_8
		this.LnWv_ = 0.0_8
		
		! Orienta la molécula, de tal forma que el eje
		! Z corresponde con el eje principal de mayor 
		! inercia y el eje X con el de menor. Hay que
		! resaltar que en el caso de una molécula lineal
		! el eje X tiene un momento de inercia nulo.
		call this.orient()
				
		if( this.fv() /= size(this.vibFrequencies) ) then
			call GOptions_error( &
				"Mass table row incomplete, this row should have at least 7 parameters", &
				"Fragment.fromMassTableRow()", &
				"Inconsistent vibrational degrees of freedom, " &
					//trim(FString_fromInteger(size(this.vibFrequencies)))//" should be " &
					//trim(FString_fromInteger(this.fv()))//" ("//trim(this.fileName)//")" &
			)
		end if
		
! 		if( GOptions_printLevel >= 4 ) then
			write(IO_STDOUT,"(A)") ""
! 			write(IO_STDOUT,"(4X,A22,I15)") "id = ", this.id
			write(IO_STDOUT,"(4X,A22,A)") "file name = ", this.fileName
			write(IO_STDOUT,"(4X,A22,A)") "name = ", this.name
			write(IO_STDOUT,"(4X,A22,A)") "formula = ", trim(this.chemicalFormula())
			write(IO_STDOUT,"(4X,A23,3F20.5,A)") "Moments of inertia = [", &
					this.diagInertiaTensor.get(1,1)/amu/angs**2, &
					this.diagInertiaTensor.get(2,2)/amu/angs**2, &
					this.diagInertiaTensor.get(3,3)/amu/angs**2, &
					"  ]   amu*angs**2"
			write(IO_STDOUT,"(4X,A23,3F20.5,A)") "Moments of inertia = [", &
					this.diagInertiaTensor.get(1,1), &
					this.diagInertiaTensor.get(2,2), &
					this.diagInertiaTensor.get(3,3), &
					"  ]   a.u."
! 			write(IO_STDOUT,"(4X,A23,3F14.5,A)") "Moments of inertia = [", &
! 					this.diagInertiaTensor.get(1,1)/amu/12.01_8, &
! 					this.diagInertiaTensor.get(2,2)/amu/12.01_8, &
! 					this.diagInertiaTensor.get(3,3)/amu/12.01_8, &
! 					"  ]   amu*angs**2/mC"
			
			write(IO_STDOUT,"(4X,A22,F15.5,A)")  "            Radius = ", this.radius( type=GOptionsM3C_radiusType )/angs, "   A"
			write(IO_STDOUT,"(4X,A22,F15.7,A)")  "             Eelec = ", this.electronicEnergy/eV, "   eV"
			write(IO_STDOUT,"(4X,A22,F15.7,A)")  "             Eelec = ", this.electronicEnergy, "   a.u."
			write(IO_STDOUT,"(4X,A22,F15.7,A)")  "              Mass = ", this.mass()/amu, "   amu"
			if( this.nAtoms() > 1 ) then
! 				write(IO_STDOUT,"(4X,A22,F15.7,A)")  "   aver. vib. freq = ", product(this.vibFrequencies)**(1.0_8/size(this.vibFrequencies))/eV, "   eV"
				write(IO_STDOUT,"(4X,A22,F15.7,A)")  "               ZPE = ", this.ZPE/eV, "   eV"
! 				write(IO_STDOUT,"(4X,A22,F15.7,A)")  "               ZPE = ", this.ZPE, "   a.u."
			end if
			if( this.isTransitionState ) then
				write(IO_STDOUT,"(4X,A22,L)")        " isTransitionState = ", this.isTransitionState
			end if
			write(IO_STDOUT,"(4X,A23,2I5,A)")    "          (fr, fv) = (", this.fr(), this.fv(), "  )"
! 			write(IO_STDOUT,"(A)") ""
! 		end if
		
		! Esto no depende de nada, así que se puede evaluar desde el principio
		call this.updateLnWe()
		call this.updateLogVtheta() ! Esto solo depende de que orient() se haya llamado antes
	end subroutine fromMassTableRow
	
	!>
	!! @brief Copy constructor
	!!
	subroutine copyFragment( this, other )
		class(Fragment), intent(inout) :: this
		class(Fragment), intent(in) :: other
		
		if( allocated(this.fileName) ) deallocate(this.fileName)
		
		call this.copyMolecule( other )
		
		this.id = other.id
		this.charge = other.charge
		this.multiplicity = other.multiplicity
		this.degL = other.degL
		this.sigmaSym = other.sigmaSym
		this.electronicEnergy = other.electronicEnergy
		
		if( allocated(this.vibFrequencies) ) deallocate(this.vibFrequencies)
		
		allocate( this.vibFrequencies( size(other.vibFrequencies) ) )
		this.vibFrequencies = other.vibFrequencies
		this.isTransitionState = other.isTransitionState
		
		this.maxEvib = other.maxEvib
		this.ZPE = other.ZPE
		this.maxJ = other.maxJ
		this.fileName = other.fileName
		
		this.frozen = other.frozen
		
		this.vibrationalEnergy_ = other.vibrationalEnergy_
		this.J_ = other.J_
		this.logVtheta_ = other.logVtheta_
		this.logVJ_ = other.logVJ_
		
		this.LnWe_ = other.LnWe_
		this.LnWv_ = other.LnWv_
		
		this.label_ = other.label_
		this.dlabel_ = other.dlabel_
		this.testLabel_ = other.testLabel_
		
	end subroutine copyFragment
	
	!>
	!! @brief Destructor
	!!
	subroutine destroyFragment( this )
		type(Fragment) :: this
		
		if( allocated(this.vibFrequencies) ) deallocate(this.vibFrequencies)
	end subroutine destroyFragment
	
	!>
	!! @brief
	!!
	subroutine updateLabel( this )
		class(Fragment) :: this 
		
		character(100), allocatable :: tokens(:)
		
		this.dlabel_ = this.name
		
		call FString_split( this.name, tokens, "(" )
		
		this.label_ = trim(tokens(1))
		
		deallocate( tokens )
	end subroutine updateLabel
	
	!>
	!! @brief
	!!
	function label( this, details ) result( output )
		class(Fragment) :: this 
		character(100) :: output
		logical, optional :: details
		
		logical :: effDetails
		
		effDetails = .true.
		if( present(details) ) effDetails = details
		
		if( .not. this.testLabel_ ) then
			call this.updateLabel()
			this.testLabel_ = .true.
		end if
		
		if( effDetails ) then
			output = this.dlabel_
		else
			output = this.label_
		end if
	end function label
	
	!>
	!! @brief Load data from Rich XYZ format file
	!!
	subroutine loadRXYZ( this, fileName, loadName )
		class(Fragment) :: this
		character(*), intent(in) :: fileName
		logical, optional, intent(in) :: loadName
		
		type(IFStream) :: ifile
		type(String) :: buffer
		character(1000), allocatable :: tokens(:)
		integer :: i, nItems
		
		call ifile.init( trim(fileName) )
		
		call this.loadXYZ( ifile, loadName )
		
		do while( .not. ifile.eof() )
			buffer = ifile.readLine()
			if( buffer.length() /= 0 ) then
				call buffer.split( tokens, " " )
				
				if( trim(tokens(1)) == "FREQUENCIES" ) then
					nItems = FString_toInteger(tokens(2))
					
					if( allocated(this.vibFrequencies) ) deallocate(this.vibFrequencies)
					allocate( this.vibFrequencies(nItems) )
					
					do i=1,nItems
						buffer = ifile.readLine()
						
						this.vibFrequencies(i) = buffer.toReal()*cm1
					end do
				end if
			end if
		end do
		
		call ifile.close()
		if( allocated(tokens) ) deallocate(tokens)
	end subroutine loadRXYZ
	
	!>
	!! @brief Load data from MOLDEN format file
	!!
	subroutine loadMOLDEN( this, fileName, loadName )
		class(Fragment) :: this
		character(*), intent(in) :: fileName
		logical, optional, intent(in) :: loadName
		
		type(IFStream) :: ifile
		type(String) :: buffer
		character(1000), allocatable :: tokens(:)
		integer :: i, nItems
		
		logical :: advance
		type(StringList) :: frequencyBlock
		type(StringListIterator), pointer :: iter
		
		call ifile.init( trim(fileName) )
		
		call this.loadGeomMOLDEN( ifile, loadName )
		
		advance = .true.
		do while( .not. ifile.eof() )
			if( advance ) then
				buffer = ifile.readLine()
				call buffer.split( tokens, " " )
				advance = .true.
			end if
			
			!-----------------------------------
			! search for Frequencies data
			!-----------------------------------
			if( tokens(1) == "[FREQ]" ) then
				do while( .not. ifile.eof() )
					buffer = ifile.readLine()
					call buffer.split( tokens, " " )
					
					if( index( tokens(1), "[" ) == 1 ) then
						advance = .false.
						exit
					end if
					
					call frequencyBlock.append( buffer )
				end do
			end if
			
			if( frequencyBlock.size() > 0 ) exit
				
			if( .not. advance ) then
				buffer = ifile.readLine()
				call buffer.split( tokens, " " )
			end if
		end do
		
		if( this.nAtoms() > 1 .and. frequencyBlock.size() == 0 ) then
			call GOptions_error( &
				"Frequencies not found in MOLDEN file ("//trim(ifile.name)//")", &
				"Fragment.loadMOLDEN()", &
				"Fragment = "//trim(this.name) &
			)
		end if
		
		call ifile.close()
		if( allocated(tokens) ) deallocate( tokens )
		
		! Check for zero frequencies
		nItems=0
		iter => frequencyBlock.begin
		do while( associated(iter) )
			if( abs(iter.data.toReal()) > 1.0_8 ) then
				nItems = nItems + 1
			end if
			
			iter => iter.next
		end do
		
		if( allocated(this.vibFrequencies) ) deallocate(this.vibFrequencies)
		allocate( this.vibFrequencies( nItems ) )
		
		i=1
		iter => frequencyBlock.begin
		do while( associated(iter) )
			if( abs(iter.data.toReal()) > 1.0_8 ) then
				this.vibFrequencies(i) = iter.data.toReal()*cm1
				i = i + 1
			end if
			
			iter => iter.next
		end do
	end subroutine loadMOLDEN
	
	!>
	!! @brief Builds a random vibrational energy
	!!
	subroutine changeVibrationalEnergy( this )
		class(Fragment) :: this
		
		if( this.frozen ) then
			this.vibrationalEnergy_ = 0.0_8
		else
			if( GOptionsM3C_useZPECorrection ) then
				this.vibrationalEnergy_ = RandomUtils_uniform( [0.0_8,this.maxEvib] )
			else
				this.vibrationalEnergy_ = RandomUtils_uniform( [this.ZPE,this.maxEvib] )
			end if
		end if
		
		if( GOptions_printLevel >= 4 ) then
			call GOptions_paragraph( "Change vibrational energy "//trim(this.name), indent=3 )
			call GOptions_valueReport( "ZPE", this.ZPE/eV, "eV", indent=3 )
			call GOptions_valueReport( "maxEvib", this.maxEvib/eV, "eV", indent=3 )
			call GOptions_valueReport( "Evib", this.vibrationalEnergy_/eV, "eV", indent=3 )
		end if
		
		call this.updateLnWv()
	end subroutine changeVibrationalEnergy
	
	!>
	!! @brief Random choosing of the orientation
	!!
	subroutine changeOrientation( this, force1D, force2D, value )
		class(Fragment) :: this
		logical, optional, intent(in) :: force1D
		logical, optional, intent(in) :: force2D
		integer, optional, intent(in) :: value(3)
		
		type(Matrix) :: effSfAxes
		logical :: effForce1D
		logical :: effForce2D
		
		effForce1D = .false.
		effForce2D = .false.
		if( present(force1D) ) effForce1D = force1D
		if( present(force2D) ) effForce2D = force2D
		
		if( present(value) ) then
			call GOptions_error( &
				"This method has not implemented yet", &
				"Fragment.changeOrientation( value= )", &
				"Fragment = "//trim(this.name) &
			)
		else
			call this.rotate( random=.true., force1D=effForce1D, force2D=effForce2D )
		end if
		
		if( GOptions_printLevel >= 4 ) then
			call GOptions_subsection( "Change orientation "//trim(this.name), indent=3 )
		end if
		
		call this.updateLogVtheta( effForce1D, effForce2D )
		
	end subroutine changeOrientation
	
! 	!>
! 	!! @brief Random choosing of the angular momentum J
! 	!!
! 	subroutine changeAngularMomentum( this, force1D, force2D, value )
! 		class(Fragment) :: this
! 		logical, optional, intent(in) :: force1D
! 		logical, optional, intent(in) :: force2D
! 		integer, optional, intent(in) :: value(3)
! 		
! 		type(Matrix) :: effSfAxes
! 		logical :: effForce1D
! 		logical :: effForce2D
! 		
! 		real(8) :: J, theta, phi
! 		integer :: Jmax
! 		
! 		effForce1D = .false.
! 		effForce2D = .false.
! 		if( present(force1D) ) effForce1D = force1D
! 		if( present(force2D) ) effForce2D = force2D
! 		
! #define Ic(k) this.diagInertiaTensor.get(k,k)
! 		if( this.nAtoms() == 1 .or. this.frozen ) then
! 		
! 			this.J_ = 0.0_8   ! Jx = Jy = Jz = 0
! 			
! 		else if( present(value) ) then
! 		
! 			this.J_ = value
! 			
! 		else
! 			J = RandomUtils_uniform( [0,this.maxJ] )
! 			theta = RandomUtils_uniform( [0.0_8,Math_PI] )
! 			phi = RandomUtils_uniform( [0.0_8,2.0_8*Math_PI] )
! 			
! 			if( effForce1D ) then
! 				this.J_(1) = 0.0_8
! 				this.J_(2) = J
! 				this.J_(3) = 0.0_8
! 				
! 			else if( effForce2D ) then
! 				this.J_(1) = 0.0_8
! 				this.J_(2) = J*cos(phi)
! 				this.J_(3) = J*sin(phi)
! 				
! 			else
! 				select case( this.fr() )
! 					case(3)
! 						this.J_(1) = J*sin(theta)*cos(phi)
! 						this.J_(2) = J*sin(theta)*sin(phi)
! 						this.J_(3) = J*cos(theta)
! 						
! 					case(2)
! 						this.J_(1) = 0.0_8
! 						this.J_(2) = J*cos(phi)
! 						this.J_(3) = J*sin(phi)
! 						
! 					case default
! 						call GOptions_error( &
! 							"Inconsistent number of rotational degrees of freedom (fr="//FString_fromInteger(this.fr())//")", &
! 							"Fragment.changeAngularMomentum()", &
! 							"Fragment = "//trim(this.name) &
! 						)
! 				end select
! 			end if
! 		end if
! #undef Ic
! 		
! 		if( GOptions_printLevel >= 4 ) then
! 			call GOptions_paragraph( "Change angular momentum "//trim(this.name), indent=3 )
! 			call GOptions_valueReport( "J", this.J_, "a.u.", indent=3 )
! 			call GOptions_valueReport( "|J|", norm2(this.J_), "a.u.", indent=3 )
! 		end if
! 		
! 		call this.updateLogVJ( effForce1D, effForce2D )
! 	end subroutine changeAngularMomentum
	
	!>
	!! @brief Molecular weight
	!!
	pure function LnW( this ) result( output )
		class(Fragment), intent(in) :: this
		real(8) :: output
		
		output = this.LnWe_ + this.LnWv_
	end function LnW
	
	!>
	!! @brief Electronic weight
	!!
	pure function LnWe( this ) result( output )
		class(Fragment), intent(in) :: this
		real(8) :: output
		
		output = this.LnWe_
	end function LnWe

	!>
	!! @brief Vibrational weight
	!!
	pure function LnWv( this ) result( output )
		class(Fragment), intent(in) :: this
		real(8) :: output
		
		output = this.LnWv_
	end function LnWv
	
	!>
	!! @brief Electronic weight
	!!
	subroutine updateLnWe( this )
		class(Fragment) :: this
		
		this.LnWe_ = log( real(this.multiplicity*this.degL,8) )
		
		if( GOptions_printLevel >= 4 ) then
			call GOptions_valueReport( "LnWe", this.LnWe_, indent=3 )
		end if
	end subroutine updateLnWe
	
	!>
	!! @brief Rotational volume
	!!
	subroutine updateLogVtheta( this, force1D, force2D )
		class(Fragment) :: this
		logical, optional, intent(in) :: force1D
		logical, optional, intent(in) :: force2D
		
		logical :: effForce1D
		logical :: effForce2D
		
		effForce1D = .false.
		if( present(force1D) ) effForce1D = force1D
		
		! @todo Esto hay que implementarlo
		effForce2D = .false.
		if( present(force2D) ) effForce2D = force2D
		
		this.logVtheta_ = 0.0_8
		
		! @todo Hay que estar seguro que en el caso que no se considere la rotación
		!       el elemento de volumen se mantiene o desaparece
		if( this.nAtoms() /= 1 .and. .not. this.frozen ) then
			if( effForce1D ) then
			
				if( GOptions_printLevel >= 4 ) then
					write(IO_STDOUT,"(A)") "** Forcing 1D"
				end if
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Caso átomo-diátomo
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				this.logVtheta_ = -log(2.0_8*Math_PI) - log(real(this.sigmaSym,8)) + log(2.0_8)
				
			else if( effForce2D ) then
			
				if( GOptions_printLevel >= 4 ) then
					write(IO_STDOUT,"(A)") "** Forcing 2D"
				end if
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Caso átomo-molécula
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				this.logVtheta_ = -2.0_8*log(2.0_8*Math_PI) - log(real(this.sigmaSym,8)) + log(4.0_8*Math_PI)
			else
				this.logVtheta_ = -this.fr()*log(2.0_8*Math_PI) - log(real(this.sigmaSym,8))
				
				select case( this.fr() )
					case(3)
						this.logVtheta_ = this.logVtheta_ + log(8.0_8*Math_PI**2)
					case(2)
						this.logVtheta_ = this.logVtheta_ + log(4.0_8*Math_PI)
					case(0)
						this.logVtheta_ = this.logVtheta_ + log(2.0_8)
					case default
						call GOptions_error( &
							"Inconsistent number of rotational degrees of freedom (fr="//FString_fromInteger(this.fr())//")", &
							"Fragment.updateLogVtheta()", &
							"Fragment = "//trim(this.name) &
						)
				end select
			end if
			
		end if

		
		if( GOptions_printLevel >= 4 ) then
			call GOptions_valueReport( "fr", this.fr(), indent=3 )
			call GOptions_valueReport( "logVtheta", this.logVtheta_, indent=3 )
		end if
		
	end subroutine updateLogVtheta
	
	!>
	!! @brief Angular momentum volume
	!!
	subroutine updateLogVJ( this, force1D, force2D )
		class(Fragment) :: this
		logical, optional, intent(in) :: force1D
		logical, optional, intent(in) :: force2D
		
		logical :: effForce1D
		logical :: effForce2D
		
		effForce1D = .false.
		if( present(force1D) ) effForce1D = force1D
		
		! @todo Esto hay que implementarlo
		effForce2D = .false.
		if( present(force2D) ) effForce2D = force2D
		
		this.logVJ_ = 0.0_8
		
		! @todo Hay que estar seguro que en el caso que no se considere la rotación
		!       el elemento de volumen se mantiene o desaparece
! 		if( this.nAtoms() /= 1 .and. .not. this.frozen ) then
! 			if( .not. effForce1D ) then
				this.logVJ_ = - this.fr()*log(2.0_8*Math_PI) - log(real(this.sigmaSym,8))
				
				select case( this.fr() )
					case(3)
						this.logVJ_ = this.logVJ_ + log(8.0_8*Math_PI**2)
					case(2)
						this.logVJ_ = this.logVJ_ + log(4.0_8*Math_PI)
					case default
						this.logVJ_ = 0.0_8
				end select
! 			else
! 				if( GOptions_printLevel >= 4 ) then
! 					write(IO_STDOUT,"(A)") "** Forcing 1D"
! 				end if
! 
! 				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 				! Caso átomo-diátomo
! 				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 				this.logVJ_ = - log(2.0_8*Math_PI) - log(real(this.sigmaSym,8)) + log(2.0_8)
! 			end if
! 		end if
		
		if( GOptions_printLevel >= 4 ) then
			call GOptions_valueReport( "fr", this.fr(), indent=3 )
			call GOptions_valueReport( "logVJ", this.logVJ_, indent=3 )
		end if
		
	end subroutine updateLogVJ
	
	!>
	!! @brief Vibrational weight
	!!
	subroutine updateLnWv( this )
		class(Fragment) :: this
		
		integer :: i, eff_fv
		real(8) :: ssum
		
		if( this.nAtoms() == 1 .or. this.frozen ) then
			this.LnWv_ = 0.0_8
			
			if( GOptions_printLevel >= 4 ) then
				write(*,*) "frozen"
			end if
		else
			ssum = 0.0_8
			eff_fv = 0
			do i=1,size(this.vibFrequencies)
				if( this.vibFrequencies(i) > 0.0_8 ) then
					ssum = ssum + log(this.vibFrequencies(i))
					eff_fv = eff_fv + 1
				end if
			end do
			
			this.LnWv_ = real(eff_fv-1,8)*log(this.vibrationalEnergy_) &
						-log(Gamma(real(eff_fv,8))) - ssum
			
! 			this.LnWv_ = real(this.fv()-1,8)*log(this.vibrationalEnergy_) &
! 						-log(Gamma(real(this.fv(),8))) - sum(log(this.vibFrequencies))

! 			this.LnWv_ = real(this.fv()-1,8)*log(this.vibrationalEnergy_) &
! 						-log(Gamma(real(this.fv(),8))) - real(this.fv(),8)*log( averVibFreq ) ! Estilo Sergio
		end if
		
		if( GOptions_printLevel >= 4 ) then
			call GOptions_valueReport( "LnWv", this.LnWv_, indent=3 )
		end if
	end subroutine updateLnWv
	
	!>
	!! @brief
	!!
	subroutine showLnWComponents( this )
		class(Fragment) :: this
		
		if( GOptions_printLevel >= 4 ) then
			write(6,"(A)") "C         "//&
				trim(FString_fromReal(this.LnWe_,"(F10.5)"))// &
				trim(FString_fromReal(this.LnWv_,"(F10.5)"))// &
				trim(FString_fromReal(this.LnWe_+this.LnWv_,"(F10.5)"))// &
				"     "//trim(this.name)
		end if
	end subroutine showLnWComponents
	
	!>
	!! @brief
	!!
	function spinAvailable( this ) result( output )
		class(Fragment), intent(in) :: this
		type(RealVector) :: output
		
		real(8) :: S, Si, Sj
		integer :: i, j
		
		call output.init()
		
		S = (this.multiplicity-1.0_8)/2.0_8
		call output.append( S )
	end function spinAvailable
	
!       Real*8 Function RotSN(Linear,Group)
!       Implicit Real*8 (A-H,O-Z)
! C
! C ROTational Symmetry Number
! C     Returns the rotational symmetry number for a given point group,
! C     as follows:
! C                   POINT GROUP                    ROTATION SYMM.NUMBER
! C
! C C1,CI, CS, C(inf)v                                         1
! C C2,C2V,C2H           D(inf)H                               2
! C C3,C3V,C3H                          S6                     3
! C C4,C4V,C4H           D2,D2D,D2H                            4
! C C6,C6V,C6H           D3,D3D,D3H                            6
! C                      D4,D4D,D4H                            8
! C                      D6,D6D,D6H     T ,TD                 12
! C                                     OH                    24
! C
! C Input
! C     Group  : Point group, given as a string
! C
! C Output:
! C     Linear : Indicates if the molecule is linear
! C
! C     Input
!       Character Group*(*)
! C     Output
!       Logical Linear
! C
!       Linear = (Group(1:2).eq.'C*').or.(Group(1:2).eq.'D*')
!       RotSN = 1.0D0
!       If((Group(1:2).eq.'C1').or.(Group(1:2).eq.'CI')) Return
!       If((Group(1:2).eq.'CS').or.(Group(1:2).eq.'C*')) Return
!       RotSN = 2.0D0
!       If((Group(1:2).eq.'C2').or.(Group(1:2).eq.'D*')) Return
!       RotSN = 3.0D0
!       If((Group(1:2).eq.'C3').or.(Group(1:2).eq.'S6')) Return
!       RotSN = 4.0D0
!       If((Group(1:2).eq.'C4').or.(Group(1:2).eq.'D2')) Return
!       RotSN = 6.0D0
!       If((Group(1:2).eq.'C6').or.(Group(1:2).eq.'D3')) Return
!       RotSN = 8.0D0
!       If(Group(1:2).eq.'D4') Return
!       RotSN = 12.0D0
!       If(Group(1:2).eq.'D6'.or.Group(1:1).eq.'T')Return
!       RotSN = 24.0D0
!       If(Group(1:2).eq.'OH ') Return
!       RotSN = 0.0D0
!       Return
!       End

!       Subroutine NewLab(IOut,NReps,Group,ChrLab)
!       Implicit Real*8(A-H,O-Z)
! C
! C NEW symmetry LABels
! C     Returns the number of irreducible representations, as well as
! C     their label Group symmetry given in Group
! C
! C Input:
! C     Group  : Group symmetry (character string)
! C
! C Output:
! C     NReps  : Number of irreducible representations
! C     ChrLab : (NReps) Labels of the irreducible representations
! C
! C     Dimensions
!       Integer MAXREP
!       Parameter (MAXREP=15)
! C     Input
!       Character Group*(*)
! C     Output
!       Character ChrLab(MAXREP)*4
! C     Local
!       Integer i
! C
!       Do 100 i = 1, MAXREP
!   100   ChrLab(i) = '    '
! C
!       If(Index(Group,'OH  ').ne.0) then
!         NReps = 10
!         ChrLab( 1) = 'A1G'
!         ChrLab( 2) = 'A1U'
!         ChrLab( 3) = 'A2G'
!         ChrLab( 4) = 'A2U'
!         ChrLab( 5) = 'EG'
!         ChrLab( 6) = 'EU'
!         ChrLab( 7) = 'T1G'
!         ChrLab( 8) = 'T1U'
!         ChrLab( 9) = 'T2G'
!         ChrLab(10) = 'T2U'
!         ChrLab(11) = '??'
!       else if(Index(Group,'O   ').ne.0) then
!         NReps = 5
!         ChrLab(1) = 'A1'
!         ChrLab(2) = 'A2'
!         ChrLab(3) = 'E'
!         ChrLab(4) = 'T1'
!         ChrLab(5) = 'T2'
!         ChrLab(6) = '??'
!       else if(Index(Group,'TD  ').ne.0) then
!         NReps = 5
!         ChrLab(1) = 'A1'
!         ChrLab(2) = 'A2'
!         ChrLab(3) = 'E'
!         ChrLab(4) = 'T1'
!         ChrLab(5) = 'T2'
!         ChrLab(6) = '??'
!       else if(Index(Group,'T   ').ne.0) then
!         NReps = 3
!         ChrLab(1) = 'A'
!         ChrLab(2) = 'E'
!         ChrLab(3) = 'T'
!         ChrLab(4) = '??'
!       else if(Index(Group,'D*  ').ne.0 .or.
!      $        Index(Group,'D*H ').ne.0) then
!         NReps = 8
!         ChrLab(1) = 'SGG'
!         ChrLab(2) = 'SGU'
!         ChrLab(3) = 'PIG'
!         ChrLab(4) = 'PIU'
!         ChrLab(5) = 'DLTG'
!         ChrLab(6) = 'DLTU'
!         ChrLab(7) = 'PHIG'
!         ChrLab(8) = 'PHIU'
!         ChrLab(9) = '??'
!       else if(Index(Group,'D2H').ne.0) then
!         NReps = 8
!         ChrLab(1) = 'AG'
!         ChrLab(2) = 'AU'
!         ChrLab(3) = 'B1G'
!         ChrLab(4) = 'B1U'
!         ChrLab(5) = 'B2G'
!         ChrLab(6) = 'B2U'
!         ChrLab(7) = 'B3G'
!         ChrLab(8) = 'B3U'
!         ChrLab(9) = '??'
!       else if(Index(Group,'D3H').ne.0) then
!         NReps = 6
!         ChrLab(1) = 'A1'''
!         ChrLab(2) = 'A1"'
!         ChrLab(3) = 'A2'''
!         ChrLab(4) = 'A2"'
!         ChrLab(5) = 'E'''
!         ChrLab(6) = 'E"'
!         ChrLab(7) = '??'
!       else if(Index(Group,'D4H').ne.0) then
!         NReps = 10
!         ChrLab( 1) = 'A1G'
!         ChrLab( 2) = 'A1U'
!         ChrLab( 3) = 'A2G'
!         ChrLab( 4) = 'A2U'
!         ChrLab( 5) = 'B1G'
!         ChrLab( 6) = 'B1U'
!         ChrLab( 7) = 'B2G'
!         ChrLab( 8) = 'B2U'
!         ChrLab( 9) = 'EG'
!         ChrLab(10) = 'EU'
!         ChrLab(11) = '??'
!       else if(Index(Group,'D5H').ne.0) then
!         NReps = 8
!         ChrLab(1) = 'A1'''
!         ChrLab(2) = 'A1"'
!         ChrLab(3) = 'A2'''
!         ChrLab(4) = 'A2"'
!         ChrLab(5) = 'E1'''
!         ChrLab(6) = 'E1"'
!         ChrLab(7) = 'E2'''
!         ChrLab(8) = 'E2"'
!         ChrLab(9) = '??'
!       else if(Index(Group,'D6H').ne.0) then
!         NReps = 12
!         ChrLab( 1) = 'A1G'
!         ChrLab( 2) = 'A1U'
!         ChrLab( 3) = 'A2G'
!         ChrLab( 4) = 'A2U'
!         ChrLab( 5) = 'B1G'
!         ChrLab( 6) = 'B1U'
!         ChrLab( 7) = 'B2G'
!         ChrLab( 8) = 'B2U'
!         ChrLab( 9) = 'E1G'
!         ChrLab(10) = 'E1U'
!         ChrLab(11) = 'E2G'
!         ChrLab(12) = 'E2U'
!         ChrLab(13) = '??'
!       else if(Index(Group,'D2D').ne.0) then
!         NReps = 5
!         ChrLab(1) = 'A1'
!         ChrLab(2) = 'A2'
!         ChrLab(3) = 'B1'
!         ChrLab(4) = 'B2'
!         ChrLab(5) = 'E'
!         ChrLab(6) = '??'
!       else if(Index(Group,'D3D').ne.0) then
!         NReps = 6
!         ChrLab(1) = 'A1G '
!         ChrLab(2) = 'A1U'
!         ChrLab(3) = 'A2G'
!         ChrLab(4) = 'A2U'
!         ChrLab(5) = 'EG'
!         ChrLab(6) = 'EU'
!         ChrLab(7) = '??'
!       else if(Index(Group,'D4D').ne.0) then
!         NReps = 7
!         ChrLab(1) = 'A1'
!         ChrLab(2) = 'A2'
!         ChrLab(3) = 'B1'
!         ChrLab(4) = 'B2'
!         ChrLab(5) = 'E1'
!         ChrLab(6) = 'E2'
!         ChrLab(7) = 'E3'
!         ChrLab(8) = '??'
!       else if(Index(Group,'D2 ').ne.0) then
!         NReps = 4
!         ChrLab(1) = 'A'
!         ChrLab(2) = 'B1'
!         ChrLab(3) = 'B2'
!         ChrLab(4) = 'B3'
!         ChrLab(5) = '??'
!       else if(Index(Group,'D3 ').ne.0) then
!         NReps = 3
!         ChrLab(1) = 'A1'
!         ChrLab(2) = 'A2'
!         ChrLab(3) = 'E'
!         ChrLab(4) = '??'
!       else if(Index(Group,'D4 ').ne.0) then
!         NReps = 5
!         ChrLab(1) = 'A1'
!         ChrLab(2) = 'A2'
!         ChrLab(3) = 'B1'
!         ChrLab(4) = 'B2'
!         ChrLab(5) = 'E'
!         ChrLab(6) = '??'
!       else if(Index(Group,'D6 ').ne.0) then
!         NReps = 6
!         ChrLab(1) = 'A1'
!         ChrLab(2) = 'A2'
!         ChrLab(3) = 'B1'
!         ChrLab(4) = 'B2'
!         ChrLab(5) = 'E1'
!         ChrLab(6) = 'E2'
!         ChrLab(7) = '??'
!       else if(Index(Group,'C*V ').ne.0 .or.
!      $        Index(Group,'C*  ').ne.0) then
!         NReps = 4
!         ChrLab(1) = 'SG'
!         ChrLab(2) = 'PI'
!         ChrLab(3) = 'DLTA'
!         ChrLab(4) = 'PHI'
!         ChrLab(5) = '??'
!       else if(Index(Group,'C2H').ne.0) then
!         NReps = 4
!         ChrLab(1) = 'AG'
!         ChrLab(2) = 'AU'
!         ChrLab(3) = 'BG'
!         ChrLab(4) = 'BU'
!         ChrLab(5) = '??'
!       else if(Index(Group,'C3H').ne.0) then
!         NReps = 4
!         ChrLab(1) = 'A'''
!         ChrLab(2) = 'A"'
!         ChrLab(3) = 'E'''
!         ChrLab(4) = 'E"'
!         ChrLab(5) = '??'
!       else if(Index(Group,'C4H').ne.0) then
!         NReps = 6
!         ChrLab(1) = 'AG'
!         ChrLab(2) = 'AU'
!         ChrLab(3) = 'BG'
!         ChrLab(4) = 'BU'
!         ChrLab(5) = 'EG'
!         ChrLab(6) = 'EU'
!         ChrLab(7) = '??'
!       else if(Index(Group,'C6H').ne.0) then
!         NReps = 8
!         ChrLab(1) = 'AG'
!         ChrLab(2) = 'AU'
!         ChrLab(3) = 'BG'
!         ChrLab(4) = 'BU'
!         ChrLab(5) = 'E1G'
!         ChrLab(6) = 'E1U'
!         ChrLab(7) = 'E2G'
!         ChrLab(8) = 'E2U'
!         ChrLab(9) = '??'
!       else if(Index(Group,'C2V').ne.0) then
!         NReps = 4
!         ChrLab(1) = 'A1'
!         ChrLab(2) = 'A2'
!         ChrLab(3) = 'B1'
!         ChrLab(4) = 'B2'
!         ChrLab(5) = '??'
!       else if(Index(Group,'C3V').ne.0) then
!         NReps = 3
!         ChrLab(1) = 'A1'
!         ChrLab(2) = 'A2'
!         ChrLab(3) = 'E'
!         ChrLab(4) = '??'
!       else if(Index(Group,'C4V').ne.0) then
!         NReps = 5
!         ChrLab(1) = 'A1'
!         ChrLab(2) = 'A2'
!         ChrLab(3) = 'B1'
!         ChrLab(4) = 'B2'
!         ChrLab(5) = 'E'
!         ChrLab(6) = '??'
!       else if(Index(Group,'C5V').ne.0) then
!         NReps = 4
!         ChrLab(1) = 'A1'
!         ChrLab(2) = 'A2'
!         ChrLab(3) = 'E1'
!         ChrLab(4) = 'E2'
!         ChrLab(5) = '??'
!       else if(Index(Group,'C6V').ne.0) then
!         NReps = 6
!         ChrLab(1) = 'A1'
!         ChrLab(2) = 'A2'
!         ChrLab(3) = 'B1'
!         ChrLab(4) = 'B2'
!         ChrLab(5) = 'E1'
!         ChrLab(6) = 'E2'
!         ChrLab(7) = '??'
!       else if(Index(Group,'CS  ').ne.0) then
!         NReps = 2
!         ChrLab(1) = 'A'''
!         ChrLab(2) = 'A"'
!         ChrLab(3) = '??'
!       else if(Index(Group,'CI  ').ne.0) then
!         NReps = 2
!         ChrLab(1) = 'AG'
!         ChrLab(2) = 'AU'
!         ChrLab(3) = '??'
!       else if(Index(Group,'C2 ').ne.0) then
!         NReps = 2
!         ChrLab(1) = 'A'
!         ChrLab(2) = 'B'
!         ChrLab(3) = '??'
!       else if(Index(Group,'C3 ').ne.0) then
!         NReps = 2
!         ChrLab(1) = 'A'
!         ChrLab(2) = 'E'
!         ChrLab(3) = '??'
!       else if(Index(Group,'C4 ').ne.0) then
!         NReps = 3
!         ChrLab(1) = 'A'
!         ChrLab(2) = 'B'
!         ChrLab(3) = 'E'
!         ChrLab(4) = '??'
!       else if(Index(Group,'C6 ').ne.0) then
!         NReps = 4
!         ChrLab(1) = 'A'
!         ChrLab(2) = 'B'
!         ChrLab(3) = 'E1'
!         ChrLab(4) = 'E2'
!         ChrLab(5) = '??'
!       else
!         NReps = 1
!         ChrLab(1) = 'A1'
!         ChrLab(2) = '??'
!       endIf
!       Return
!       End

	!>
	!! @brief
	!!
	subroutine Fragment_getLabelFragments( clusters, label, dlabel, idSorted )
		type(Fragment), allocatable, intent(in) :: clusters(:)
		character(200), intent(out) :: label
		character(200), intent(out) :: dlabel
		integer, allocatable, intent(out), optional :: idSorted(:)
		
		integer :: i
		real(8), allocatable :: massVec(:) ! @todo Hay que hacer el Math_sort para enteros
		integer, allocatable :: effIdSorted(:)
		
		allocate( massVec(size(clusters)) )
		allocate( effIdSorted(size(clusters)) )
		
		do i=1,size(clusters)
		
			massVec(i) = clusters(i).mass()/amu &
					+ clusters(i).charge/10.0 &
						+ clusters(i).multiplicity/100.0_8 

! 			massVec(i) = 10000000*clusters(i).mass()/amu + 1000000*clusters(i).nAtoms() + 100*clusters(i).charge &
! 					+ clusters(i).multiplicity

		end do
			
		call Math_sort( massVec, effIdSorted )
		
		label = ""
		dlabel = ""
		do i=1,size(clusters)
			if( i /= size(clusters) ) then
				label = trim(label)//trim(clusters( effIdSorted(i) ).label( details=.false. ))//"+"
				dlabel = trim(dlabel)//trim(clusters( effIdSorted(i) ).label( details=.true. ))//"+"
			else
				label = trim(label)//trim(clusters( effIdSorted(i) ).label( details=.false. ))
				dlabel = trim(dlabel)//trim(clusters( effIdSorted(i) ).label( details=.true. ))
			end if
		end do
		
		if( present(idSorted) ) then
			if( allocated(idSorted) .and. size(idSorted) == size(effIdSorted) ) then
				idSorted = effIdSorted
			else
				if( allocated(idSorted) ) deallocate(idSorted)
				allocate( idSorted(size(clusters)) )
				idSorted = effIdSorted
			end if
		end if
		
		deallocate( massVec )
		deallocate( effIdSorted )
	end subroutine Fragment_getLabelFragments
	
	!>
	!! @brief Test method
	!!
	subroutine Fragment_test()
		character(:), allocatable :: fstr
		type(Fragment) :: cluster, cluster2
		type(Matrix) :: Im, Rot
		integer :: i
		
		fstr = "      slC4    0  1  0    2   prueba.xyz     -4129.754238  0.4159"
		
		write(IO_STDOUT,*) "================================"
		call cluster.init( fstr )
! 		call cluster.show()
		
! 		fstr = "      slC5    0  1  0    2   prueba.xyz     -5165.261895  0.4159"
		fstr = "      slC5    0  1  0    2   prueba.xyz     -5165.261895  0.4159"
		
		do i=1,1000000
		write(*,*) i
		
		call cluster.init( fstr )
! 		call cluster.show()
		
		cluster.maxEvib = 0.4561*eV
		cluster.maxJ = 100
		
		call cluster.changeVibrationalEnergy()
		call cluster.changeOrientation()
! 		write(*,*) "logVtheta = ", cluster.logVtheta_
! 		write(*,*) "LnWe = ", cluster.LnWe()
! 		write(*,*) "LnWv = ", cluster.LnWv()
! 		write(*,*) " LnW = ", cluster.LnW()
		
		call cluster.setCenter( [-1.5_8, 3.8_8, 2.8_8] )
		call cluster.rotate( random=.true. )
! 		write(IO_STDOUT,"(A,3F5.1)") "newCenter = ", cluster.center()
		
		cluster = cluster2
		
		end do
		
! 		call cluster.save("salida.xyz")
	end subroutine Fragment_test
	
end module Fragment_

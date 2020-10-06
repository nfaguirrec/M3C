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
module FragmentsList_
	use GOptions_
	use AtomicElementsDB_
	use String_
	use Math_
	use Matrix_
	use SpecialMatrix_
	use RandomSampler_
	
	use GOptionsM3C_
	use FragmentsListBase_
	
	implicit none
	private
	
	public :: &
		FragmentsList_test
	
	type, public, extends( FragmentsListBase ) :: FragmentsList
		real(8) :: rotationalEnergy_      !< It is calculated into updateRotationalEnergy procedure
		real(8) :: E_totJ                !< Experimental en versión 1.9
		
		real(8), private :: LnLambda_     !< Translational weight
		real(8), private :: LnDiagI_     !< Contiene el log del producto de la diagonal de los tensores de inercia efectivos
		
		contains
			generic :: init => initFragmentsList
			generic :: assignment(=) => copyFragmentsList
			
			procedure :: initFragmentsList
			procedure :: copyFragmentsList
			final :: destroyFragmentsList
			
			procedure :: initialGuessFragmentsList
			procedure :: changeGeometry
			procedure :: changeVibrationalEnergy
			procedure :: changeOrientations
			
			procedure :: totalEnergy
			procedure :: internalEnergy
			procedure :: translationalEnergy
			procedure :: rotationalEnergy
			
			procedure :: energyHistoryLine
			procedure :: weightHistoryLine
				
			procedure :: LnW
			procedure :: LnLambda
			procedure :: LnDiagI
			procedure, private :: updateLambda
! 			procedure :: showLnWComponents
			
			procedure, private :: updateDiagInertiaTensor
			procedure, private :: updateDiagInertiaTensorJJ
			procedure, private :: updateDiagInertiaTensorJJL
			
			procedure :: iTemperature
	end type FragmentsList
	
	contains
	
	!>
	!! @brief Default constructor
	!!
	!! @param[in] nMolecules Number of clusters that will be contained
	!!
	subroutine initFragmentsList( this, nMolecules )
		class(FragmentsList) :: this 
		integer, intent(in) :: nMolecules
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_subsection( "INITIALIZING CLUSTERLIST", indent=2 )
		end if
		
		call this.initFragmentsListBase( nMolecules )
		
		this.rotationalEnergy_ = 0.0_8
		this.E_totJ = 0.0_8
		this.LnLambda_ = 0.0_8
		this.LnDiagI_ = 0.0_8
	end subroutine initFragmentsList
	
	!>
	!! @brief Copy constructor
	!!
	subroutine copyFragmentsList( this, other )
		class(FragmentsList), intent(inout) :: this
		class(FragmentsList), intent(in) :: other
		
		call this.copyFragmentsListBase( other )
		
		this.rotationalEnergy_ = other.rotationalEnergy_
		this.E_totJ = other.E_totJ
		this.LnLambda_ = other.LnLambda_
		this.LnDiagI_ = other.LnDiagI_
	end subroutine copyFragmentsList
	
	!>
	!! @brief Destructor
	!!
	subroutine destroyFragmentsList( this )
		type(FragmentsList) :: this
		
		call this.destroyFragmentsListBase()
	end subroutine destroyFragmentsList	

	!>
	!! @brief
	!!
	subroutine initialGuessFragmentsList( this )
		class(FragmentsList) :: this
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_section( "BUILDING INITIAL CONFIGURATION "//trim(this.label()), indent=2 )
			write(IO_STDOUT,*) ""
		end if
		
		this.forceInitializing = .false.
		
		! Este frozen, no estoy seguro si debería ir aquí, pero para no
		! escribir más código lo dejo así por el momento
! 		call this.setFrozen( .true. )
		
		call this.initialGuessFragmentsListBase() ! Inicializa geometría (incluidas orientaciones) y parte vibracional
! 		call this.updateKineticEnergy()
! 		
! 		call this.changeAngularMomenta() ! Inicializa el momento angular
! 		call this.updateLogVJ()
! 		call this.updateKineticEnergy()
! 		
! 		call this.updateLambda()
		
! 		call this.setFrozen( .false. )

		call this.changeVibrationalEnergy()
		call this.changeGeometry()
		call this.changeOrientations()
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_section( "END BUILDING INITIAL CONFIGURATION "//trim(this.label()), indent=2 )
		end if
	end subroutine initialGuessFragmentsList
	
	!>
	!! @brief
	!!
	subroutine changeGeometry( this )
		class(FragmentsList) :: this
		
		integer :: j
		
		if( this.forceInitializing ) then
			call this.initialGuessFragmentsList()
			return
		end if
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_section( "CHANGE GEOMETRY "//trim(this.label()), indent=2 )
			write(IO_STDOUT,*) ""
		end if

		call this.changeGeometryFragmentsListBase()
		
		!-----------------------------------------------------------------------------------------
		! @todo Probablemente este bloque debería ir dentro de changeGeometryFragmentsListBase
		!       ya que es un potencial centrifugo. Al igula habría que cambiar V(r) --> V(r)+E_totJ
		if( norm2(GOptionsM3C_totalJ(:)) > 1d-5 ) then
			this.E_totJ = 0.0_8
			do j=3,4-this.fl(),-1
				this.E_totJ = this.E_totJ + 0.5_8*GOptionsM3C_totalJ(j)**2/this.diagInertiaTensor( j )
			end do
			
! 			write(*,*) this.E_totJ
		end if
		!-----------------------------------------------------------------------------------------
		
		call this.updateKineticEnergy()
		call this.updateLambda()
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_section( "END CHANGE GEOMETRY "//trim(this.label()), indent=2 )
		end if
	end subroutine changeGeometry
	
	!>
	!! @brief
	!!
	subroutine changeVibrationalEnergy( this )
		class(FragmentsList) :: this
		
		if( this.forceInitializing ) then
			call this.initialGuessFragmentsList()
			return
		end if
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_section( "CHANGE VIBRATIONAL ENERGY "//trim(this.label()), indent=2 )
		end if

! 		call this.changeVibrationalEnergyFragmentsListBase()
! 		call this.updateKineticEnergy()
! 		call this.updateLambda()

		call this.changeVibrationalEnergyFragmentsListBase()
		call this.updateKineticEnergy()
		call this.updateLambda()
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_section( "END CHANGE VIBRATIONAL ENERGY "//trim(this.label()), indent=2 )
		end if
	end subroutine changeVibrationalEnergy
	
	!>
	!! @brief 
	!!
	subroutine changeOrientations( this )
		class(FragmentsList) :: this
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_section( "CHANGE ORIENTATIONS "//trim(this.label()), indent=2 )
		end if
		
		call this.changeOrientationsFragmentsListBase()
		call this.updateDiagInertiaTensor()
		call this.updateLambda()
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_section( "END CHANGE ORIENTATIONS "//trim(this.label()), indent=2 )
		end if
	end subroutine changeOrientations
	
	!>
	!! @brief Total weight
	!!
	function LnW( this ) result( output )
		class(FragmentsList) :: this
		real(8) :: output
		
		output = this.LnWn() + this.LnWe() + this.LnWv() + this.LnLambda_
	end function LnW
	
	!>
	!! @brief
	!!
	function LnLambda( this ) result( output )
		class(FragmentsList) :: this
		real(8) :: output
		
		output = this.LnLambda_
	end function LnLambda
	
	!>
	!! @brief
	!!
	function LnDiagI( this ) result( output )
		class(FragmentsList) :: this
		real(8) :: output
		
		output = this.LnDiagI_
	end function LnDiagI

	! Creo que esto ya no lo uso
! 	subroutine showLnWComponents( this )
! 		class(FragmentsList) :: this
! 		
! 		integer :: i
! 		
! 		write(6,*) ""
! 		write(6,"(A)") &
! 			trim(FString_fromReal(this.LnWe(),"(F10.5)"))// &
! 			trim(FString_fromReal(this.LnWv(),"(F10.5)"))// &
! 			trim(FString_fromReal(this.LnWn(),"(F10.5)"))// &
! 			trim(FString_fromReal(this.LnDiagI_,"(F10.5)"))// &
! 			trim(FString_fromReal(this.LnLambda_,"(F10.5)"))// &
! 			trim(FString_fromReal(this.LnWe()+this.LnWv()+this.LnWn()+this.LnLambda_,"(F10.5)"))// &
! 			"     "//trim(this.label())
! 	end subroutine showLnWComponents
	
	!>
	!! @brief Returns the total energy
	!!
	function totalEnergy( this ) result( output )
		class(FragmentsList), intent(in) :: this
		real(8) :: output
		
		output = this.kineticEnergy() + this.internalEnergy()
	end function totalEnergy
	
	!>
	!! @brief Return the internal energy
	!!
	function internalEnergy( this ) result( output )
		class(FragmentsList), intent(in) :: this
		real(8) :: output
		
		output = this.vibrationalEnergy_  + this.intermolEnergy_ + this.E_totJ
	end function internalEnergy
	
	!>
	!! @brief Return the translacional energy
	!!
	pure function translationalEnergy( this ) result( output )
		class(FragmentsList), intent(in) :: this
		real(8) :: output
		
		output = this.kineticEnergy() - this.rotationalEnergy_
	end function translationalEnergy
	
	!>
	!! @brief Return the translacional energy
	!!
	pure function rotationalEnergy( this ) result( output )
		class(FragmentsList), intent(in) :: this
		real(8) :: output
		
		output = this.rotationalEnergy_
	end function rotationalEnergy
	
	!>
	!! @brief
	!!
	function energyHistoryLine( this, prefix ) result( output )
		class(FragmentsList) :: this
		type(String) :: output
		character(*), optional, intent(in) :: prefix
		
		character(1000) :: line
		character(100) :: prefixEff
		
		prefixEff = ""
		if( present(prefix) ) prefixEff = prefix
		
		write(line,"(1X,A2,1X,5F15.5,5X,A)") &
			trim(prefixEff), &
			this.translationalEnergy()/eV, this.intermolEnergy_/eV, &
			this.vibrationalEnergy_/eV, &
			this.rotationalEnergy_/eV, &
			this.totalEnergy()/eV, trim(this.label())
			
		output = line
	end function energyHistoryLine
	
	!>
	!! @brief
	!!
	function weightHistoryLine( this, prefix ) result( output )
		class(FragmentsList) :: this
		type(String) :: output
		character(*), optional, intent(in) :: prefix
		
		character(1000) :: line
		character(100) :: prefixEff
		
		prefixEff = ""
		if( present(prefix) ) prefixEff = prefix
		
		write(line,"(1X,A2,1X,6F15.5,5X,A)") &
			trim(prefixEff), &
			this.LnWe(), this.LnWv(), &
! 			this.LnWn(), 0.0_8, this.LnLambda_, &
			this.LnWn(), 0.5*this.LnDiagI_+this.logVtheta_, this.LnLambda_-0.5*this.LnDiagI_-this.logVtheta_, &
			this.LnW(), trim(this.label())
			
		output = line
	end function weightHistoryLine
	
	!>
	!! @brief Actualiza la energía rotacional
	!!
	subroutine updateDiagInertiaTensor( this )
		class(FragmentsList) :: this
		
		integer :: j
		
		select case( trim(GOptionsM3C_angularMomentumCouplingScheme.fstr) )
			case( "JJ" )
				call this.updateDiagInertiaTensorJJ()
			case( "JJL" )
				call this.updateDiagInertiaTensorJJL()
			case default
				call GOptions_error( &
					"Unknown angular momentum coupling scheme"//" ("//trim(GOptionsM3C_angularMomentumCouplingScheme.fstr)//")", &
					"FragmentsListBase.updateDiagInertiaTensor()", &
					"Implemented methods: JJ, JJL" &
					)
		end select
		
	end subroutine updateDiagInertiaTensor
	
	!>
	!! @brief Actualiza la energía rotacional
	!!
	subroutine updateDiagInertiaTensorJJ( this )
		class(FragmentsList) :: this
		
		integer :: i, j, n
		integer :: mu, nu, effMu, effNu
		type(Matrix) :: invIn    !< Tensor de inercia de N proyectado sobre los ejes de i y su inversa
		type(Matrix) :: invIi    !< Tensor de inercia de i proyectado sobre los ejes de i (diagonal) y su inversa
		type(Matrix) :: invIt    !< Tensor de inercia de i proyectado sobre los ejes de i (diagonal) y su inversa
		type(Matrix) :: RotMu, RotNu      !< Matrices de rotación
		type(Matrix) :: invBigI
		
		type(Matrix) :: Ui
		type(Matrix) :: invBi
		real(8) :: maxErot
		integer :: fr_sf
		integer :: molRef
		
		if( this.forceInitializing ) then
			call this.initialGuessFragmentsList()
! 			return
		end if
		
		maxErot = this.kineticEnergy()
		if( maxErot < 0.0_8 ) return
		
		this.LnDiagI_ = 0.0_8
		
		n = this.nMolecules()
		if( n == 1 ) return
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Calcula el numero efectivo de fragmentos no atomicos
		! e inicializa la matriz I grande
		effMu = 0  ! << Comienza en cero porque es solo para contar
		effNu = 0  ! << Comienza en cero porque es solo para contar
		do mu=1,n-1
			if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
				effNu = 0  ! << Comienza en cero porque es solo para contar
				do nu=1,n-1
					if( this.clusters( this.idSorted(nu) ).fr() /= 0 ) then
						effNu = effNu + 1
					end if
				end do
				effMu = effMu + 1
			end if
		end do
		
		if( GOptions_debugLevel >= 3 ) then
			call GOptions_subsection( "Angular momentum coupling JJ --> "//trim(this.label()), indent=2 )
		end if
		
		if( effNu > 1 .and. effMu /= effNu ) then
			write(*,"(A,2I8)") "FragmentsList.updateDiagInertiaTensorJJ(). effMu /= effNu, ", effMu, effNu
			stop 
		end if
		
		if( effMu < 1 ) return
		
		call invBigI.init( 3*effMu, 3*effNu, 0.0_8 )
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se obtiene la matriz inversa del tensor
		! de inercia de n. Es diagonal en body fix
		call invIn.init( 3, 3, 0.0_8 )
		if( this.clusters( this.idSorted(n) ).fr() /= 0 ) then
			do j=3,4-this.clusters( this.idSorted(n) ).fr(),-1
				invIn.data( j, j ) = 1.0_8/this.clusters( this.idSorted(n) ).diagInertiaTensor.data( j, j )
			end do
		end if
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Creación de los Ij
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		molRef = n
		
		effMu = 1
		do mu=1,n-1
		
			if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Matriz de rotación que permite transformar a los ejes
				! de N en los ejes de i, pasando por el sistem fix
				! RotMu = R_n^T*R_mu
				RotMu = SpecialMatrix_rotationTransform( &
					this.clusters( this.idSorted(mu) ).inertiaAxes(), &
					this.clusters( this.idSorted(molRef) ).inertiaAxes() )
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				effNu = 1
				do nu=1,n-1
					
					if( this.clusters( this.idSorted(nu) ).fr() /= 0 ) then
					
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! Matriz de rotación que permite transformar a los ejes
						! de N en los ejes de i, pasando por el sistem fix
						! RotNu = R_n^T*R_nu
						RotNu = SpecialMatrix_rotationTransform( &
							this.clusters( this.idSorted(nu) ).inertiaAxes(), &
							this.clusters( this.idSorted(molRef) ).inertiaAxes() )
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						
						if( mu == nu ) then
						
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							! Se obtiene la matriz inversa del tensor
							! de inercia de i
							call invIi.init( 3, 3, 0.0_8 )
							do j=3,4-this.clusters( this.idSorted(mu) ).fr(),-1
								invIi.data( j, j ) = 1.0_8/this.clusters( this.idSorted(mu) ).diagInertiaTensor.data( j, j )
							end do
							
							if( GOptions_debugLevel >= 3 ) then
								write(*,*) "Base(", mu, ",",  nu, ")"
								call invIi.show( formatted=.true., precision=10 )
! 								write(*,*) "Base(", mu, ",",  nu, ") inverse"
! 								Ii = invIi.inverse()
! 								call Ii.show( formatted=.true., precision=8 )
								write(*,*) ""
							end if
							
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							! Se obtiene la matriz inversa del tensor
							! de inercia efectivo.
							invIt =  invIi + RotMu*invIn*RotMu.transpose()
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							
							if( GOptions_debugLevel >= 3 ) then
								write(*,*) "JJ(", effMu, ",",  effNu, ") = "
								call invIt.show( formatted=.true., precision=10 )
								write(*,*) ""
							end if
							
							do i=1,3
								do j=1,3
									call invBigI.set( 3*(effMu-1)+i, 3*(effNu-1)+j, invIt.get(i,j) )
								end do
							end do
							
						else
						
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							! Se obtiene la matriz de acoplamieto
							invIt =  RotMu*invIn*RotNu.transpose()
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							
							if( GOptions_debugLevel >= 3 ) then
								write(*,*) "JJ(", effMu, ",",  effNu, ")"
								call invIt.show( formatted=.true., precision=10 )
								write(*,*) ""
							end if
							
							do i=1,3
								do j=1,3
									call invBigI.set( 3*(effMu-1)+i, 3*(effNu-1)+j, invIt.get(i,j) )
								end do
							end do
						
						end if
						
						effNu = effNu + 1
					end if
					
				end do
				
				effMu = effMu + 1
			end if
			
		end do
		
		if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
			write(*,*) ""
			write(*,*) "Total invI matrix"
			call invBigI.show( formatted=.true., precision=10 )
			write(*,*) ""
		end if
		
		call invBigI.eigen( eVals=invBi, eVecs=Ui )
		
		fr_sf = this.fr() - this.clusters( this.idSorted(n) ).fr()
		
		if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
			write(*,*) "effective fr = ", fr_sf
			write(*,*) ""
			write(*,*) "Total diag(invI) matrix"
			call invBi.show( formatted=.true., precision=10 )
			write(*,*) ""
			call GOptions_subsection( "LnDiagI_ contributions", indent=3 )
			write(*,"(10X,A5,A20,A20)") "i", "diag(invI)_i", "log(I_i)"
			write(*,"(10X,A5,A20,A20)") "---", "---------", "---------------"
		end if
			
		this.LnDiagI_ = 0.0_8
		do effMu=invBi.nRows,invBi.nRows-fr_sf+1,-1
			if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
				write(*,"(10X,I5,F20.8,F20.5)") effMu, invBi.get(effMu,effMu), -log(invBi.get(effMu,effMu))
			end if
			
			this.LnDiagI_ = this.LnDiagI_ - log(invBi.get(effMu,effMu))
		end do
		
		if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
			write(*,"(10X,5X,20X,A20)") "---------------"
			write(*,"(5X,A10,20X,F20.5)") "Total", this.LnDiagI_
			write(*,*) ""
			stop
		end if
		
	end subroutine updateDiagInertiaTensorJJ
	
	!>
	!! @brief Actualiza la energía rotacional
	!!
	subroutine updateDiagInertiaTensorJJL( this )
		class(FragmentsList) :: this
		
		integer :: i, j, n
		integer :: mu, nu, effMu, effNu, nMu, nNu
		type(Matrix) :: invIn    !< Tensor de inercia de N proyectado sobre los ejes de i y su inversa
		type(Matrix) :: invIi    !< Tensor de inercia de i proyectado sobre los ejes de i (diagonal) y su inversa
		type(Matrix) :: invIt    !< Tensor de inercia de i proyectado sobre los ejes de i (diagonal) y su inversa
		type(Matrix) :: RotMu, RotNu      !< Matrices de rotación
		type(Matrix) :: invBigI
		
		type(Matrix) :: Ui, Jprime
		type(Matrix) :: invBi
		
		real(8) :: maxErot
		
		real(8) :: randNumber
		integer :: randDirection
		type(Matrix) :: Jmu, vecL
		
		if( this.forceInitializing ) then
			call this.initialGuessFragmentsList()
! 			return
		end if
		
		maxErot = this.kineticEnergy()
		if( maxErot < 0.0_8 ) return
		
		this.LnDiagI_ = 0.0_8
		
		n = this.nMolecules()
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Calcula el numero efectivo de fragmentos no atomicos
		! e inicializa la matriz I grande
		effMu = 0  ! << Comienza en cero porque es solo para contar
		effNu = 0  ! << Comienza en cero porque es solo para contar
		do mu=1,n
			if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
				effNu = 0  ! << Comienza en cero porque es solo para contar
				do nu=1,n
					if( this.clusters( this.idSorted(nu) ).fr() /= 0 ) then
						effNu = effNu + 1
					end if
				end do
				effMu = effMu + 1
			end if
		end do
		
		nMu = effMu
		nNu = effNu
		
		! @TODO Agregar los grados de libertad de los productos del TS y sumarlos a nMu y nNu
		
		if( GOptions_debugLevel >= 3 ) then
			call GOptions_subsection( "Angular momentum coupling JJL --> "//trim(this.label()), indent=2 )
		end if
		
		if( nNu > 1 .and. nMu /= nNu ) then
			write(*,"(A,2I8)") "FragmentsList.updateDiagInertiaTensorJJL(). nMu /= nNu, ", nMu, nNu
			stop 
		end if
		
		if( nMu < 1 ) return
		
		call invBigI.init( 3*nMu, 3*nNu, 0.0_8 )
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se obtiene la matriz inversa del tensor
		! de inercia de n. Es diagonal en body fix
		call invIn.init( 3, 3, 0.0_8 )
		if ( n == 1 ) then
			do j=3,4-this.clusters(1).fr(),-1
				invIn.data( j, j ) = 1.0_8/this.clusters(1).diagInertiaTensor.get( j, j )
			end do
		else
			do j=3,4-this.fl(),-1
				invIn.data( j, j ) = 1.0_8/this.diagInertiaTensor( j )
			end do
		end if
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Creación de los Ij
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		effMu = 1
		do mu=1,n
		
			if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Matriz de rotación que permite transformar a los ejes
				! de N en los ejes de i, pasando por el sistem fix
				! RotMu = R_n^T*R_mu
				RotMu = SpecialMatrix_rotationTransform( &
					this.clusters( this.idSorted(mu) ).inertiaAxes(), &
					this.inertiaAxes )
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				effNu = 1
				do nu=1,n
					
					if( this.clusters( this.idSorted(nu) ).fr() /= 0 ) then
					
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! Matriz de rotación que permite transformar a los ejes
						! de N en los ejes de i, pasando por el sistem fix
						! RotNu = R_n^T*R_nu
						RotNu = SpecialMatrix_rotationTransform( &
							this.clusters( this.idSorted(nu) ).inertiaAxes(), &
							this.inertiaAxes )
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						
						if( mu == nu ) then
						
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							! Se obtiene la matriz inversa del tensor
							! de inercia de i
							call invIi.init( 3, 3, 0.0_8 )
							do j=3,4-this.clusters( this.idSorted(mu) ).fr(),-1
								invIi.data( j, j ) = 1.0_8/this.clusters( this.idSorted(mu) ).diagInertiaTensor.data( j, j )
							end do
							
							if( GOptions_debugLevel >= 3 ) then
								write(*,*) "Base(", mu, ",",  nu, ")"
								call invIi.show( formatted=.true., precision=10 )
								write(*,*) ""
							end if
							
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							! Se obtiene la matriz inversa del tensor
							! de inercia efectivo.
							invIt =  invIi + RotMu*invIn*RotMu.transpose()
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							
							if( GOptions_debugLevel >= 3 ) then
								write(*,*) "JJ(", effMu, ",",  effNu, ") = "
								call invIt.show( formatted=.true., precision=10 )
								write(*,*) ""
							end if
							
							do i=1,3
								do j=1,3
									call invBigI.set( 3*(effMu-1)+i, 3*(effNu-1)+j, invIt.get(i,j) )
								end do
							end do
							
						else
						
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							! Se obtiene la matriz de acoplamieto
							invIt =  RotMu*invIn*RotNu.transpose()
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							
							if( GOptions_debugLevel >= 3 ) then
								write(*,*) "JJ(", effMu, ",",  effNu, ")"
								call invIt.show( formatted=.true., precision=10 )
								write(*,*) ""
							end if
							
							do i=1,3
								do j=1,3
									call invBigI.set( 3*(effMu-1)+i, 3*(effNu-1)+j, invIt.get(i,j) )
								end do
							end do
						
						end if
						
						effNu = effNu + 1
					end if
					
				end do
				
				effMu = effMu + 1
			end if
			
		end do
		
		! @TODO Aca hay que agregar los valores que faltan a la matriz grande que vienen de los productos del TS
		
		if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
			write(*,*) ""
			write(*,*) "Total invI matrix"
			call invBigI.show( formatted=.true., precision=10 )
			write(*,*) ""
		end if
		
		call invBigI.eigen( eVals=invBi, eVecs=Ui )
		
		if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
			write(*,*) "fr = ", this.fr()
			write(*,*) ""
			write(*,*) "Total diag(invI) matrix"
			call invBi.show( formatted=.true., precision=10 )
			write(*,*) ""
			call GOptions_subsection( "LnDiagI_ contributions", indent=3 )
			write(*,"(10X,A5,A20,A20)") "i", "diag(invI)_i", "log(I_i)"
			write(*,"(10X,A5,A20,A20)") "---", "---------", "---------------"
		end if
			
		this.LnDiagI_ = 0.0_8
		do effMu=invBi.nRows,invBi.nRows-this.fr()+1,-1
			if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
				write(*,"(10X,I5,F20.8,F20.5)") effMu, invBi.get(effMu,effMu), -log(invBi.get(effMu,effMu))
			end if
			
			this.LnDiagI_ = this.LnDiagI_ - log(invBi.get(effMu,effMu))
		end do
		
		if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
			write(*,"(10X,5X,20X,A20)") "---------------"
			write(*,"(10X,5X,20X,F20.5)") this.LnDiagI_
			write(*,*) ""
		end if
		
		!--------------------------------------------------------------
		! Muestreo de la energía rotacional
		! Este codigo no necesita la dirección ni identidad de los Jmu
		!--------------------------------------------------------------
		do while( .true. )
			
			this.rotationalEnergy_ = 0.0_8
			
			do effMu=invBi.nRows,invBi.nRows-this.fr()+1,-1
				call random_number( randNumber ) ! [0-1]
				
				this.rotationalEnergy_ = this.rotationalEnergy_ + randNumber**2*maxErot
			end do
			
			if( this.rotationalEnergy_ <= maxErot ) exit
		end do
		
! 		call invBigI.eigenNotSorted( eVals=invBi, eVecs=Ui )  ! Este era un intento de diagonalizar la matriz sin ordenar los valores ni vectores propios
! 		call Jprime.columnVector( 3*nMu )
! 		
! 		do while( .true. )
! 			
! 			this.rotationalEnergy_ = 0.0_8
! 			
! 			do i=1,3*nMu
! 			
! 				if( .not. Math_isInf( 1.0_8/invBi.get(i,i) ) ) then
! 					call random_number( randNumber ) ! [0-1]
! 					randDirection = merge( 1, 0, randNumber>0.5_8 )  ! [0,1]
! 					
! 					call random_number( randNumber ) ! [0-1]
! 					
! 					call Jprime.set( i, 1, (-1.0_8)**randDirection*randNumber*sqrt( 2.0_8*maxErot/invBi.get( i, i ) ) )
! 					
! 					this.rotationalEnergy_ = this.rotationalEnergy_ + 0.5_8*Jprime.get( i, 1 )**2*invBi.get( i, i )
! 					
! 					write(*,"(I3,3F15.7)") i, invBi.get( i, i ), Jprime.get( i, 1 ), this.rotationalEnergy_
! 				end if
! 				write(*,*) ">>"
! 			end do
			
! 			El problema que hay es que despues de la diagonalización
! 			no hay forma de saber que vector de momento angular o porción
! 			de la matriz invBi, pertenece a que molécula, ya que esta
! 			información se pierde durante la diagonalización
! 			write(*,*) "Total diag(invI) matrix (non-sort)"
! 			call invBi.show( formatted=.true., precision=10 )
! 			write(*,*) "MMMM0", this.rotationalEnergy_, maxErot, trim(this.label()), nMu
! 			
! 			effMu = 1
! 			do mu=1,n
! 				if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
! 					
! 					call random_number( randNumber ) ! [0-1]
! 					randDirection = merge( 1, 0, randNumber>0.5_8 )  ! [0,1]
! 					
! 					call random_number( randNumber ) ! [0-1]
! 					
! 					do i=3,4-this.clusters( this.idSorted(mu) ).fr(),-1
! 						call Jprime.set( 3*(effMu-1)+i, 1, (-1.0_8)**randDirection*randNumber*sqrt( 2.0_8*maxErot/invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) )
! 						this.rotationalEnergy_ = this.rotationalEnergy_ + 0.5_8*Jprime.get( 3*(effMu-1)+i, 1 )**2*invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i )
! 						
! 						write(*,"(3I3,3F15.7)") mu, i, 3*(effMu-1)+i, invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i ), Jprime.get( 3*(effMu-1)+i, 1 ), this.rotationalEnergy_
! 					end do
! 							
! 					effMu = effMu + 1
! 				end if
! 			end do
! 
! 			write(*,*) "MMMM1", this.rotationalEnergy_, maxErot, trim(this.label()), nMu
! 			if( this.rotationalEnergy_ <= maxErot ) exit
! 			
! 		end do
! 		
! 		write(*,*) "Jprime"
! 		call Jprime.show( formatted=.true. )
! 		
! 		Jprime = Ui.transpose()*Jprime
! 		
! 		write(*,*) "Jprime2"
! 		call Jprime.show( formatted=.true. )
! 		
! 		call vecL.columnVector( 3 )
! 		effMu = 1
! 		do mu=1,n
! 			if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
! 				
! 				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 				! Matriz de rotación que permite transformar a los ejes
! 				! de N en los ejes de i, pasando por el sistem fix
! 				! RotMu = R_n^T*R_mu
! 				RotMu = SpecialMatrix_rotationTransform( &
! 					this.clusters( this.idSorted(mu) ).inertiaAxes(), &
! 					this.inertiaAxes )
! 				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 				
! 				call Jmu.columnVector( 3 )
! 				do i=3,4-this.clusters( this.idSorted(mu) ).fr(),-1
! 					call Jmu.set( i, 1, Jprime.get( 3*(effMu-1)+i, 1 ) )
! 				end do
! 				
! 				! Contribución a L
! 				this.clusters( this.idSorted(mu) ).J_ = Jmu.data(:,1)
! 				vecL = vecL - RotMu.transpose()*Jmu
! 				
! 				write(*,"(A,I2,3F10.5)") "Jmu ", mu, this.clusters( this.idSorted(mu) ).J_
! 						
! 				effMu = effMu + 1
! 			end if
! 		end do
! 		
! 		this.L_ = vecL.data(:,1)
! 		write(*,"(A,3F10.5)")  "L ", this.L_
	end subroutine updateDiagInertiaTensorJJL
	
	!>
	!! @brief
	!!
! 	subroutine updateRotationalEnergy( this )
! 		class(FragmentsList) :: this
! 		
! 		if( GOptions_printLevel >= 3 ) then
! 			write(*,"(10X,5X,20X,A20)") "----------"
! 			write(*,"(10X,5X,20X,F20.5)") weight
! 			write(*,*) ""
! 		
! 			call GOptions_subsection( "Angular momenta", indent=3 )
! 			write(*,"(10X,A30)") "           J(au)"
! 			write(*,"(10X,A30)")  "-------------------------"
! 		end if
! 
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 		! Genera los momentos angulares en sf
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 		Ein = 0.0_8
! 		effMu = 1  !< effMu da las coordenadas para localizar un valor en la matriz bigI
! 		do mu=1,n-1
! 			if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
! 			
! 				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 				! Matriz de rotación que permite transformar a los ejes
! 				! de N en los ejes de i, pasando por el sistem fix
! 				! RotMu = R_n^T*R_mu
! 				RotMu = SpecialMatrix_rotationTransform( &
! 					this.clusters( this.idSorted(effMu) ).inertiaAxes(), &
! 					this.clusters( this.idSorted(n) ).inertiaAxes() )
! 				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 				
! 				do while( .true. )
! 					tmpErot = 0.0_8
! 					do i=1,3
! ! 						if( ( .not. Math_isInf( 1.0_8/invBigI.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) ) .and. 1.0_8/invBigI.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) > maxIval ) then
! ! 						if( ( invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) > 0.0_8 .and. abs( invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) > tol ) then
! 						if( ( .not. Math_isInf( 1.0_8/invBigI.get(effMu,effMu) ) ) .and. 1.0_8/invBigI.get(effMu,effMu) > maxIval ) then
! 
! 							write(*,*) "Entro", mu
! 							
! 							call random_number( randNumber ) ! [0,1]
! 							call random_number( randDirection ) ! [0,1]
! 							randDirection = 0.0_8+2.0_8*randDirection  ! [0:2]
! 							
! 							this.clusters( this.idSorted(mu) ).J_(i) = &
! 								(-1.0_8)**int(randDirection)*sqrt( 2.0_8*maxErot/abs( invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) )*randNumber
! 							
! 							tmpErot = tmpErot + 0.5_8*this.clusters( this.idSorted(mu) ).J_(i)**2*invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i )
! 							
! 						end if
! 					end do
! 					
! 					if( tmpErot <= maxErot ) exit
! 				end do
! 				Ein = tmpErot
! 				
! 				if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
! 					write(*,"(10X,3F10.5,5X,A)") this.clusters( this.idSorted(mu) ).J_, trim(this.clusters( this.idSorted(mu) ).label())
! 				end if
! 				
! 				! @todo Hace falta una linea transformadolo J a bf (La energia si se calcula en sf, tal como esta arriba)
! 				
! 				! Momento angular para el fragmento mu-esimo
! 				call Ji.columnVector( 3, values=this.clusters( this.idSorted(mu) ).J_ )
! 				
! 				! Contribución a Jn
! 				Jn = Jn - RotMu*Ji
! 				
! 				effMu = effMu + 1
! 			end if
! 		end do
! 		
! 		this.clusters(this.idSorted(n)).J_ = Jn.data(:,1)
! 		
! 		if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
! 			write(*,"(10X,3F10.5,5X,A,A)") this.clusters( this.idSorted(n) ).J_, trim(this.clusters( this.idSorted(n) ).label()), "   <<<< ref"
! 		end if
! 		
! 		this.rotationalEnergy_ = Ein
! 		
! 		if( GOptions_debugLevel >= 3 ) then
! 			write(*,"(A,A,3F10.1)") trim(this.clusters( this.idSorted(n) ).label()), ", Jcal = ", this.clusters( this.idSorted(n) ).J_
! 			write(*,"(A,F10.5)") "Energy = ", Ein
! 			stop
! 		end if
! 		
! 		call this.updateLambda() !! <<<<<<<<<<< OOJJJOOOOO esto no me gusta aquí
! 	end subroutine updateRotationalEnergy
	
	!>
	!! @brief Translational weight
	!!
	subroutine updateLambda( this )
		class(FragmentsList) :: this
		
		integer :: i, n, s
		real(8) :: logMu, Et
		integer :: j
		real(8) :: randNumber
		integer :: sTS
		real(8) :: logTSt, Et_TS, iEt_TS
		character(100), allocatable :: tokens(:)
		character(100) :: fsBuffer
		real(8) :: reducedMass
		
		n = this.nMolecules()
		
		logTSt = 0.0_8
		Et_TS = 0.0_8
		sTS = 0
		do i=1,n
			if( this.clusters(i).isTransitionState .and. this.kineticEnergy() > 0.0_8 ) then
			
				call this.clusters(i).vibFrequenciesData(i).split( tokens,  ";" )
				fsBuffer = tokens(1)
				call FString_split( fsBuffer, tokens, "=" )
				if( tokens(1) == "rMass" ) then
					reducedMass = FString_toReal( tokens(2) )
				else
					write(*,*) "### ERROR ### FragmentsList.updateLambda: rMass undefined for cluster"
					stop
				end if
				deallocate( tokens )
				
				logTSt = &
					logTSt &
					+ 0.5_8*log(2.0_8*Math_PI) &
					- log( Gamma(0.5_8) ) &
					+ 0.5_8*log(reducedMass)
					
! 				write(*,"(A,3F15.8)") "logTStA = ", logTSt, log(reducedMass), reducedMass
					
				call random_number( randNumber ) ! [0-1]
				iEt_TS = randNumber*( this.kineticEnergy() - Et_TS )  ! Rand [0:E1]
				Et_TS = Et_TS + iEt_TS
				
				logTSt = logTSt - 0.5*log(iEt_TS)
				
! 				write(*,"(A,2F15.8)") "logTSt = ", logTSt, Et_TS
				
				sTS = sTS + 1
			end if
		end do
		
		Et = this.kineticEnergy() - Et_TS
! 		Et = this.kineticEnergy()
		
		if( Et < 0.0_8 ) then
			this.LnLambda_ = 0.0_8
! 		else if( this.nMolecules() == 1 ) then
! ! 			this.LnLambda_ = this.logVfree_ + this.logVtheta_
! 			this.LnLambda_ = 0.0_8
		else
			select case( trim(GOptionsM3C_angularMomentumCouplingScheme.fstr) )
				case( "JJ" )
					s = this.ft() + this.fr() - this.clusters( this.idSorted(n) ).fr()
				case( "JJL" )
					s = this.ft() + this.fr()
				case default
					call GOptions_error( &
						"Unknown angular momentum coupling scheme"//" ("//trim(GOptionsM3C_angularMomentumCouplingScheme.fstr)//")", &
						"FragmentsListBase.updateLambda()", &
						"Posible implemented values: JJ, JL" &
						)
			end select
			
			logMu = 0.0_8
			do i=1,n
				logMu = logMu + log(this.clusters(i).mass())
			end do
			logMu = logMu - log( this.mass() )
			
			if ( n == 1 ) then	
				this.LnLambda_ = &
					logTSt &
					+ 0.5_8*this.clusters(1).fr()*log(2.0_8*Math_PI) &
					- log( Gamma(0.5_8*this.clusters(1).fr()) ) &
					+ this.logVtheta_ &
					+ 0.5_8*this.LnDiagI_
			else
				this.LnLambda_ = &
					logTSt &
					+ 0.5_8*s*log(2.0_8*Math_PI) &
					- log( Gamma(0.5_8*s) ) &
					+ 1.5*logMu &
					+ this.logVfree_ + this.logVtheta_ &
					+ 0.5_8*this.LnDiagI_ &
					+ (0.5_8*s-1.0_8)*log(Et)
			end if
				
		end if
		
		if( ( Math_isNaN(this.LnLambda_) .or. Math_isInf(this.LnLambda_) ) .and. GOptions_printLevel >= 2 ) then
			write(*,*) ""
			call GOptions_valueReport( "reactorEnergy", this.reactorEnergy()/eV, "eV", indent=2 )
			call GOptions_valueReport( "internalEnergy", this.internalEnergy()/eV, "eV", indent=2 )
			call GOptions_valueReport( "kineticEnergy", this.kineticEnergy()/eV, "eV", indent=2 )
			call GOptions_valueReport( "vibrationalEnergy", this.vibrationalEnergy_/eV, "eV", indent=2 )
			call GOptions_valueReport( "rotationalEnergy", this.rotationalEnergy_/eV, "eV", indent=2 )
			call GOptions_valueReport( "used Et", Et/eV, "eV", indent=2 )
			call GOptions_valueReport( "logVfree", this.logVfree_, indent=2 )
			call GOptions_valueReport( "logVtheta", this.logVtheta_, indent=2 )
			call GOptions_valueReport( "logVJ", this.logVJ_, indent=2 )
			call GOptions_valueReport( "0.5*logMu", 0.5*logMu, indent=2 )
			call GOptions_valueReport( "LnLambda_", this.LnLambda_, indent=2 )
		end if
	end subroutine updateLambda
		
	!>
	!! @brief Inverse of the temperature
	!!
	function iTemperature( this ) result( output )
		class(FragmentsList), intent(in) :: this
		real(8) :: output
		
		integer :: s
		integer :: n
		
		n = this.nMolecules()
		
		select case( trim(GOptionsM3C_angularMomentumCouplingScheme.fstr) )
			case( "JJ" )
				s = this.ft() + this.fr() - this.clusters( this.idSorted(n) ).fr()
			case( "JJL" )
				s = this.ft() + this.fr()
			case default
				call GOptions_error( &
					"Unknown angular momentum coupling scheme"//" ("//trim(GOptionsM3C_angularMomentumCouplingScheme.fstr)//")", &
					"FragmentsListBase.updateLambda()", &
					"Posible implemented values: JJ, JL" &
					)
		end select
		
		output = (0.5_8*s-1.0_8)/this.kineticEnergy()
	end function iTemperature

	!>
	!! @brief Test method
	!!
	subroutine FragmentsList_test()
! 		use FragmentsDB_
! 		character(:), allocatable :: fstr
! 		
! 		type(FragmentsList) :: clist
! 		type(RigidMolecule) :: cluster
! 		type(Matrix) :: Im
! 		real(8), allocatable :: eVal(:)
! 		integer :: i, nTrials
! 		real(8) :: rBuffer
! 		
! 		call clist.init(4)
! 		
! 		write(*,*) "================================"
! 		fstr = "    tC1    0  F    3  1    0   C1.xyz             -1028.031662  0.000000"
! 		
! 		call cluster.init( fstr, 1 )
! 		call cluster.show()
! 		write(*,*) "radius = ", cluster.radius( type=GOptionsM3C_radiusType )
! 		call clist.set( 1, cluster )
! 		
! 		write(*,*) "================================"
! 		fstr = "    tlC2    0  T    3  1    2   C2T-linear.xyz     -2061.703744  0.209800       2*tC1    300"
! 		
! 		call cluster.init( fstr, 2 )
! 		call cluster.show()
! 		write(*,*) "radius = ", cluster.radius( type=GOptionsM3C_radiusType )
! 		call clist.set( 2, cluster )
! 		write(*,*) "================================"
! 		fstr = "    tcC4    0  F    3  0    4   C4T-cyclic.xyz     -4129.378872  0.089866    tC1,slC3    300"
! 		
! 		call cluster.init( fstr, 3 )
! 		call cluster.show()
! 		write(*,*) "radius = ", cluster.radius( type=GOptionsM3C_radiusType )
! 		call clist.set( 3, cluster )
! 		write(*,*) "================================"
! 		fstr = "    tcC5    0  F    3  0    1   C5T-cyclic.xyz     -5162.267819  0.093690   slC2,slC3    300"
! 		
! 		call cluster.init( fstr, 4 )
! 		call cluster.show()
! 		write(*,*) "radius = ", cluster.radius( type=GOptionsM3C_radiusType )
! 		call clist.set( 4, cluster )
! 		write(*,*) "================================"
! 		write(*,*) ""
! 		call clist.show( formatted=.true. )
! 		
! 		clist.kineticEnergy = 3.0_8*eV
! 		
! 		call clist.LnWi()
! 		call clist.LnWn()
! 		call clist.LnWp()
! 		call clist.logVtheta()
! 		
! 		write(*,*) "------------------------------------"
! 		write(*,*) " Checking reactivity"
! 		write(*,*) "------------------------------------"
! 		
! 		call clist.init(2)
! 		
! 		fstr = "    tC1    0  F    3  1    0   C1.xyz             -1028.031662  0.000000"
! 		call cluster.init( fstr, 1 )
! 		call clist.set( 1, cluster )
! 		
! 		fstr = "    tlC2    0  T    3  1    2   C2T-linear.xyz     -2061.703744  0.209800       2*tC1    300"
! 		call cluster.init( fstr, 2 )
! 		call clist.set( 2, cluster )
! 		
! 		call clist.show( formatted=.true. )
! 		clist.kineticEnergy = 3.0_8*eV
! 		rBuffer = clist.LnW()
! 		
! ! 		call clist.diagInertiaTensor( Im )
! ! 		
! ! 		call Im.show( formatted=.true. )
! ! 		call Im.eigen( eValues=eVal )
! ! 		write(*,"(A,3F20.5)") "eValues = ", eVal
! ! 		write(*,*)
! ! 		do i=1,clist.nMolecules()
! ! 			write(*,*) i, clist.clusters(i).mass()
! ! 		end do
! 		
! 		call clist.save( "superMol.xyz" )

		integer :: i, j
		type(String), allocatable :: massTable(:)
		type(String), allocatable :: potential(:)
		character(:), allocatable :: forbidden
		character(100), allocatable :: tokens(:), items(:)
		real(8) :: r, rBuffer
		type(FragmentsList) :: clist, reactives, products
		type(String) :: store
		
		allocate( massTable(3) )
		massTable(1) = "       C(s)    0  1  2    0   C1.xyz         -1026.581828"
		massTable(2) = "      C9(s)    0  1  0    2   prueba.xyz     -2062.282893    0.458"
		massTable(3) = "      C9(t)    0  1  0    2   prueba.xyz     -3097.388207    0.458"
		
! 		allocate( massTable(16) )

		!                     Label    Z  M  L  SYM         geomFile            Eelec   Chann.Vib  Jmax
! 		massTable(1)  = "       sC1    0  1  2    0   C1.xyz             -1026.581828"
! 		massTable(2)  = "       tC1    0  3  1    0   C1.xyz             -1028.031662"
! 		massTable(3)  = "      slC2    0  1  0    2   C2S-linear.xyz     -2062.282893       2*tC1   300"
! 		massTable(4)  = "      tlC2    0  3  1    2   C2T-linear.xyz     -2061.703744       2*tC1   300"
! 		massTable(5)  = "      slC3    0  1  0    2   C3S-linear.xyz     -3097.388207    tC1,slC2   400"
! 		massTable(6)  = "      tlC3    0  3  1    2   C3T-linear.xyz     -3095.315621    tC1,slC2   400"
! 		massTable(7)  = "      scC3    0  1  0    2   C3S-cyclic.xyz     -3096.460007    tC1,slC2   300"
! 		massTable(8)  = "      tcC3    0  3  0    6   C3T-cyclic.xyz     -3096.641729    tC1,slC2   300"
! 		massTable(9)  = "      slC4    0  1  0    2   C4S-linear.xyz     -4129.754238    tC1,slC3  1000"
! 		massTable(10) = "      tlC4    0  3  0    2   C4T-linear.xyz     -4129.134441    tC1,slC3  1000"
! 		massTable(11) = "      scC4    0  1  0    4   C4S-cyclic.xyz     -4130.294751    tC1,slC3   300"
! 		massTable(12) = "      tcC4    0  3  0    4   C4T-cyclic.xyz     -4129.378872    tC1,slC3   300"
! 		massTable(13) = "      slC5    0  1  0    2   C5S-linear.xyz     -5165.261895   slC2,slC3  1000"
! 		massTable(14) = "      tlC5    0  3  1    2   C5T-linear.xyz     -5162.952781   slC2,slC3  1000"
! 		massTable(15) = "      scC5    0  1  0    1   C5S-cyclic.xyz     -5160.722344   slC2,slC3   300"
! 		massTable(16) = "      tcC5    0  3  0    1   C5T-cyclic.xyz     -5162.267819   slC2,slC3   300"
		
		! En esta clase no se debería utilizar la base de datos
! 		store="."
! 		call FragmentsDB_instance.fromMassTable( massTable, store=store )
! 		
! 		GOptionsM3C_angularMomentumCouplingScheme = "JJL"
! 		
! 		call reactives.init( 3 )
! 		reactives.clusters(1) = FragmentsDB_instance.clusters(3)
! 		reactives.clusters(2) = FragmentsDB_instance.clusters(3)
! 		reactives.clusters(3) = FragmentsDB_instance.clusters(3)
! 		
! 		do i=1,1000000
! 			call reactives.changeGeometry()
! 			products = reactives
! 		end do
		
! 		write(*,*) "getEelecFromName('slC3,tcC4,sC1') = ", (-3097.388207+0.049062-4129.378872+0.089866-1026.581828), " eV"
! 		write(*,*) RigidMoleculeDatabase_instance.getEelecFromName( "slC3,tcC4,sC1" )/eV, " eV"
! 		
! 		write(*,*) "getEelecFromName('slC3,2*tcC4,3*sC1') = ", (-3097.388207+0.049062+2.0*(-4129.378872+0.089866)+3.0*(-1026.581828)), " eV"
! 		write(*,*) RigidMoleculeDatabase_instance.getEelecFromName( "slC3,2*tcC4,3*sC1" )/eV, " eV"
! 		
! 		forbidden = "  slC2,slC4,slC4  -->  tlC2,slC3,tlC2,tlC2  "
! 		call FString_split( forbidden, tokens, "->" )
! 		
! 		do i=1,size(tokens)
! 			call FString_split( tokens(i), items, "," )
! 			
! 			do j=1,size(items)
! 				write(*,*) i, trim(adjustl(items(j))), RigidMoleculeDatabase_instance.getIdClusterFromLabel( items(j) )
! 			end do
! 			write(*,*) ""
! 		end do
! 		
! 		allocate( potential(3) )
! 		potential(1) = "  tC1+tC1->tlC2    5.0   MORSE(1.5,1.0,1.0,-1.5)  "
! 		potential(2) = "  slC2+slC3->scC5  1.5   MORSE(10.5,1.0,5.0,-10.5)"
! 		potential(3) = "  slC4+tC1->slC5   1.5   MORSE(0.5,1.0,3.0,-0.5)  "
! 		
! 		call RigidMoleculeDatabase_instance.setPotentialTable( potential )
		
! 		write(*,*) "------------------------------------"
! 		write(*,*) " Checking reactivity"
! 		write(*,*) "------------------------------------"
! 		
! 		GOptionsM3C_printRigidMoleculeData = .false.
! 		
! 		call clist.init(3)
! 		call clist.setFromRigidMolecule( 1, RigidMoleculeDatabase_instance.clusters(1) )
! 		call clist.setFromRigidMolecule( 2, RigidMoleculeDatabase_instance.clusters(2) )
! 		call clist.setFromRigidMolecule( 3, RigidMoleculeDatabase_instance.clusters(3) )
! 		call clist.initialGuessFragmentsList()
! 		
! 		clist.reactorEnergy = 10.0_8*eV + clist.electronicEnergy()
! 		
! 		call clist.changeGeometry()
! 		
! 		call clist.changeVibrationalEnergy()
! 		call clist.changeAngularMomenta()
! 		
! 		rBuffer = clist.LnW()

		write(*,*) "------------------------------------"
		write(*,*) " Load geometry from file"
		write(*,*) "------------------------------------"
! 		call clist.loadXYZ( "prueba.xyz" )
! 		call clist.save()
		
	end subroutine FragmentsList_test
	
end module FragmentsList_

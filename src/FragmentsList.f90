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
		real(8) :: rotationalEnergy      !< It is calculated into updateRotationalEnergy procedure
		
		real(8), private :: LnLambda     !< Translational weight
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
			
			procedure :: energyHistoryLine
			procedure :: weightHistoryLine
				
			procedure :: LnW
			procedure, private :: updateLambda
! 			procedure :: showLnWComponents
			
			procedure, private :: N2LCorrection
			procedure, private :: updateDiagInertiaTensor
			procedure, private :: updateDiagInertiaTensorJJ
			procedure, private :: updateDiagInertiaTensorJL
			
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
		
		this.rotationalEnergy = 0.0_8
		this.LnLambda = 0.0_8
		this.LnDiagI_ = 0.0_8
	end subroutine initFragmentsList
	
	!>
	!! @brief Copy constructor
	!!
	subroutine copyFragmentsList( this, other )
		class(FragmentsList), intent(inout) :: this
		class(FragmentsList), intent(in) :: other
		
		call this.copyFragmentsListBase( other )
		
		this.rotationalEnergy = other.rotationalEnergy
		this.LnLambda = other.LnLambda
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
			write(STDOUT,*) ""
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
		
		if( this.forceInitializing ) then
			call this.initialGuessFragmentsList()
			return
		end if
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_section( "CHANGE GEOMETRY "//trim(this.label()), indent=2 )
			write(STDOUT,*) ""
		end if

		call this.changeGeometryFragmentsListBase()
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
		
		output = this.LnWn() + this.LnWe() + this.LnWv() + this.LnLambda
	end function LnW

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
! 			trim(FString_fromReal(this.LnLambda,"(F10.5)"))// &
! 			trim(FString_fromReal(this.LnWe()+this.LnWv()+this.LnWn()+this.LnLambda,"(F10.5)"))// &
! 			"     "//trim(this.label())
! 	end subroutine showLnWComponents
	
	!>
	!! @brief Returns the total energy
	!!
	function totalEnergy( this ) result( output )
		class(FragmentsList) :: this
		real(8) :: output
		
		output = this.kineticEnergy() + this.vibrationalEnergy_ + this.intermolEnergy_
	end function totalEnergy
	
	!>
	!! @brief Return the internal energy
	!!
	function internalEnergy( this ) result( output )
		class(FragmentsList) :: this
		real(8) :: output
		
		output = this.vibrationalEnergy_  + this.intermolEnergy_
	end function internalEnergy
	
	!>
	!! @brief Return the translacional energy
	!!
	function translationalEnergy( this ) result( output )
		class(FragmentsList) :: this
		real(8) :: output
		
		output = this.kineticEnergy() - this.rotationalEnergy
	end function translationalEnergy
	
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
			this.rotationalEnergy/eV, &
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
! 			this.LnWn(), 0.0_8, this.LnLambda, &
			this.LnWn(), this.LnDiagI_, this.LnLambda, &   ! << Solo por visualización del valor de LnDiagI_
			this.LnW(), trim(this.label())
			
		output = line
	end function weightHistoryLine
	
	!>
	!! @brief Actualiza la energía rotacional
	!!
	subroutine updateDiagInertiaTensor( this )
		class(FragmentsList) :: this
		
! 		select case( trim(GOptionsM3C_angularMomentumCouplingScheme.fstr) )
! 			case( "JJ" )
! 				call this.updateDiagInertiaTensorJJ()
! 			case( "JL" )
				call this.updateDiagInertiaTensorJL()
! 			case default
! 				call GOptions_error( &
! 					"Unknown angular momentum coupling scheme"//" ("//trim(GOptionsM3C_angularMomentumCouplingScheme.fstr)//")", &
! 					"FragmentsListBase.updateDiagInertiaTensor()", &
! 					"Posible implemented values: JJ, JL" &
! 					)
! 		end select
	end subroutine updateDiagInertiaTensor
	
	!>
	!! @brief Actualiza la energía rotacional
	!!
	subroutine N2LCorrection( this )
		class(FragmentsList) :: this
		
		integer :: n
		real(8) :: valI_L, Re, mu, r, m, m_1, m_2
		
		type(Matrix) :: invIi, invBi, Ui, I_L
		real(8) :: Il, molRef
		integer :: j
		type(Matrix) :: RotMu
		
		n = this.nMolecules()
	
		! Contribución orbital para el caso de 1 atomo y una molecula lineal
		if( n==2 .and. GOptionsM3C_useLCorrection ) then
! 			if( this.clusters( this.idSorted(1) ).fr() == 0 .and. this.clusters( this.idSorted(2) ).fr() /= 0 ) then
! 				if( this.clusters( this.idSorted(1) ).fr() == 0 .and. this.clusters( this.idSorted(2) ).nAtoms() == 2 ) then
! 				
! 	!                               Re = 0.5*this.clusters( this.idSorted(2) ).radius()  ! El valor exacto es radius-rcov(1)-rcov(n), el valor que está es para una diatómica
! 	!                               mu = this.clusters( this.idSorted(2) ).mass()
! 						Re = norm2(this.clusters( this.idSorted(2) ).atoms(1).r-this.clusters( this.idSorted(2) ).atoms(2).r)
! 						
! 						m_1 = AtomicElementsDB_instance.atomicMass( this.clusters( this.idSorted(2) ).atoms(1).symbol )
! 						m_2 = AtomicElementsDB_instance.atomicMass( this.clusters( this.idSorted(2) ).atoms(2).symbol )
! 						mu = m_1*m_2/(m_1+m_2)
! 						
! 						m = this.clusters( this.idSorted(1) ).mass()
! 						r = norm2(this.clusters(1).center()-this.clusters(2).center())
! 												
! 						valI_L = 1.0_8/( 1.0_8/( 2.0_8*mu*Re**2  ) + 1.0_8/( 2.0_8*m*r**2 ) )
! 												
! 						this.LnDiagI_ = this.LnDiagI_ + 1.0_8*log(sqrt(2.0_8*valI_L))

! 				Re = 0.5*this.clusters( this.idSorted(2) ).radius()  ! El valor exacto es radius-rcov(1)-rcov(n), el valor que está es para una diatómica
! 				mu = this.clusters( this.idSorted(2) ).mass()
! 
! 				r = norm2(this.clusters(1).center()-this.clusters(2).center())
! 				m = this.clusters( this.idSorted(1) ).mass()*this.clusters( this.idSorted(2) ).mass()/this.mass()
! 				
! 				valI_L = ( mu*Re**2 * m*r**2 )/( mu*Re**2  + m*r**2 )
! 				
! 				this.LnDiagI_ = this.LnDiagI_ + log(valI_L)

				!--------------------------------------
				
				! mu es el que está rotando y molRef el que está en el centro
! 				if( this.clusters( this.idSorted(1) ).fr() > this.clusters( this.idSorted(2) ).fr() ) then
! 				      mu = 2
! 				      molRef = 1
! 				else
				      mu = 1
				      molRef = 2
! 				end if

				r = norm2(this.clusters(mu).center()-this.clusters(molRef).center())
				Il = this.clusters( this.idSorted(mu) ).mass()*this.clusters( this.idSorted(molRef) ).mass()*r**2/this.mass()
				
! 				! El valor que sale para valIl es el mismo que con el tensor de inercia
! ! 				call this.buildInertiaTensor( I_L, this.clusters( this.idSorted(molRef) ).center() )
! ! ! 				call this.buildInertiaTensor( I_L )
! ! 				write(*,*) "Il = ", log(Il)
! ! 				write(*,*) "I_L = "
! ! 				call I_L.show( formatted=.true. )
! ! 				I_L.data = log(I_L.data)
! ! 				write(*,*) "log(I_L) = "
! ! 				call I_L.show( formatted=.true. )
! ! 				stop
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Calculo del tensor de inerciaq que está en el centro
				call invIi.init( 3, 3, 0.0_8 )
				do j=3,4-this.clusters( this.idSorted(molRef) ).fr(),-1
				      invIi.data( j, j ) = 1.0_8/this.clusters( this.idSorted(molRef) ).diagInertiaTensor.data( j, j )
				end do

				call invIi.eigen( eVals=invBi, eVecs=Ui )
! 				
! 				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 				! Matriz de rotación que permite transformar a los ejes
! 				! de N en los ejes de i, pasando por el sistem fix
! 				! RotMu = R_n^T*R_mu
! 				RotMu = SpecialMatrix_rotationTransform( &
! 					this.inertiaAxes(), &
! 					this.clusters( this.idSorted(molRef) ).inertiaAxes() )
! 				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 				
! 				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 				! Se obtiene la matriz inversa del tensor
! 				! de inercia efectivo.
! 				invIt =  invIi + RotMu*invIn*RotMu.transpose()
! 				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				        
! 				this.LnIm_ = 0.0_8
! 				do j=3,4-this.clusters( this.idSorted(mu) ).fr(),-1
! 				      this.LnDiagI_ = this.LnDiagI_ - log(2.0_8*invBi.get(j,j))
! 				end do
! 				this.LnDiagI_ = this.LnDiagI_ + log(Math_PI**(this.clusters( this.idSorted(mu) ).fr()/2.0)) - log(gamma(this.clusters( this.idSorted(mu) ).fr()/2.0))

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Calculo de la corrección
				valI_L = 0.0_8
				do j=3,4-this.clusters( this.idSorted(molRef) ).fr(),-1
				      valI_L = valI_L + (- log( 1.0_8 + 1.0_8/( invBi.get(j,j)*Il ) ) - log(invBi.get(j,j)))/3.0_8
				end do
! 				valI_L = - log( 1.0_8 + 1.0_8/( invBi.get(3,3)*Il ) ) - log(invBi.get(3,3))
				
				this.LnDiagI_ = this.LnDiagI_ + valI_L
				
				if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
					write(*,*) ""
					write(*,*) trim(this.label())
					write(*,*) mu, molRef, this.clusters( this.idSorted(molRef) ).fr(), invBi.get(1,1), invBi.get(2,2), invBi.get(3,3), Il
					write(*,*) log( 1.0_8 + 1.0_8/( invBi.get(1,1)*Il ) ), log( 1.0_8 + 1.0_8/( invBi.get(2,2)*Il ) ), log( 1.0_8 + 1.0_8/( invBi.get(3,3)*Il ) )
					write(*,*) log(invBi.get(1,1)), log(invBi.get(2,2)), log(invBi.get(3,3))
					write(*,"(5X,A10,20X,F20.5)") "L-corr", valI_L
					write(*,*) ""
				end if
! 			end if
		end if
		
	end subroutine N2LCorrection
	
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
		integer :: fr_sf
		integer :: molRef
		
! 		type(Matrix) :: I_L, invI_L
! 		integer :: fl_sf
		
		if( this.forceInitializing ) then
			call this.initialGuessFragmentsList()
! 			return
		end if
		
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
		
		call this.N2LCorrection()
		
		if( effNu > 1 .and. effMu /= effNu ) then
			write(*,"(A,2I8)") "FragmentsList.updateDiagInertiaTensorJJ(). effMu /= effNu, ", effMu, effNu
			stop 
		end if
		
		if( effMu < 1 ) return
		
		call invBigI.init( 3*effMu, 3*effNu, 0.0_8 )
! 		call invBigI.init( 3*effMu+3, 3*effNu+3, 0.0_8 )
		
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
		
! 		call this.buildInertiaTensor( I_L, this.clusters( this.idSorted(molRef) ).center() )
! 		invI_L = I_L.inverse()
! 		invI_L.data = invI_L.data + 1e-8
! 			
! 		if( GOptions_debugLevel >= 3 ) then
! 			write(*,*) "I_L"
! 			call I_L.show( formatted=.true., precision=10 )
! 			write(*,*) ""
! 			
! 			write(*,*) "inv I_L"
! 			call invI_L.show( formatted=.true., precision=10 )
! 			write(*,*) ""
! 		end if
! 		
! 		select case( n )
! 			case( 1 )
! 				fl_sf = 0
! 			case( 2 )
! 				fl_sf = 1
! ! 			case( 3 )
! ! 				fl_sf = 2
! 			case default
! ! 				fl_sf = 3
! 				fl_sf = 1
! 		end select
! 
! 		do i=3,4-fl_sf,-1
! 			do j=3,4-fl_sf,-1
! 				call invBigI.set( 3*(effMu-1)+i, 3*(effNu-1)+j, invI_L.get(i,j) )
! 			end do
! 		end do
		
		if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
			write(*,*) ""
			write(*,*) "Total invI matrix"
			call invBigI.show( formatted=.true., precision=10 )
			write(*,*) ""
		end if
		
		call invBigI.eigen( eVals=invBi, eVecs=Ui )
		
		fr_sf = 0
		do i=1,n-1
			fr_sf = fr_sf + this.clusters( this.idSorted(i) ).fr()
		end do
		
		if( GOptions_printLevel >= 3 .or. GOptions_debugLevel >= 3 ) then
			write(*,*) "effective fr = ", fr_sf
! 			write(*,*) "effective fl = ", fl_sf
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
! 		do effMu=invBi.nRows,invBi.nRows-fr_sf-fl_sf+1,-1
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
	subroutine updateDiagInertiaTensorJL( this )
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
		integer :: fr_sf
		
		if( this.forceInitializing ) then
			call this.initialGuessFragmentsList()
! 			return
		end if
		
		this.LnDiagI_ = 0.0_8
		
		n = this.nMolecules()
		if( n == 1 ) return
		
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
		
		if( GOptions_debugLevel >= 3 ) then
			call GOptions_subsection( "Angular momentum coupling JL --> "//trim(this.label()), indent=2 )
		end if
		
		if( effNu > 1 .and. effMu /= effNu ) then
			write(*,"(A,2I8)") "FragmentsList.updateDiagInertiaTensorJL(). effMu /= effNu, ", effMu, effNu
			stop 
		end if
		
		if( effMu < 1 ) return
		
		call invBigI.init( 3*effMu, 3*effNu, 0.0_8 )
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se obtiene la matriz inversa del tensor
		! de inercia de n. Es diagonal en body fix
		call invIn.init( 3, 3, 0.0_8 )
		do j=3,4-this.fl(),-1
				invIn.data( j, j ) = 1.0_8/this.diagInertiaTensor( j )
		end do
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
! 							invIt =  invIi + RotMu*invIn*RotMu.transpose()
							invIt =  invIi + RotMu.transpose()*invIn*RotMu
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
! 							invIt =  RotMu*invIn*RotNu.transpose()
							invIt =  RotMu.transpose()*invIn*RotNu
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
		
		fr_sf = 0
		do i=1,n
			fr_sf = fr_sf + this.clusters( this.idSorted(i) ).fr()
		end do
		
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
			write(*,"(10X,5X,20X,F20.5)") this.LnDiagI_
			write(*,*) ""
		end if
		
		! Lo he probado en varios sistemas y no funciona. No se si
		! al cambiar el radio del sistema tenga algún efecto
! 		do j=3,4-this.fl(),-1
! 			this.LnDiagI_ = this.LnDiagI_ - log(this.diagInertiaTensor( j ))
! 		end do
		
	end subroutine updateDiagInertiaTensorJL
	
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
! 		this.rotationalEnergy = Ein
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
		
		integer :: fr
		
		n = this.nMolecules()
		
		Et = this.kineticEnergy()
		
		if( Et < 0.0_8 ) then
			this.LnLambda = 0.0_8
! 			write(*,"(A,A,F10.5)") "E < 0", "   "//trim(this.label()), this.LnLambda
		else if( this.nMolecules() == 1 ) then
! 			this.LnLambda = this.logVfree_ + this.logVtheta_
			this.LnLambda = 0.0_8
! 			write(*,"(A,A,F10.5)") "n = 0", "   "//trim(this.label()), this.LnLambda
		else
			logMu = 0.0_8
			fr = 0
			do i=1,n
				logMu = logMu + log(this.clusters(i).mass())
				
! 				if( i /= n ) then
					fr = fr + this.clusters( this.idSorted(i) ).fr()
! 				end if
			end do
! 			logMu = logMu - log( this.mass() )
			
! 			if( n==2 .and. GOptionsM3C_useLCorrection ) then
! 				s = this.ft() + this.fl() + fr + 1
! 			else
				s = this.ft() + this.fl() + fr
! 			end if
			
			this.LnLambda = \
				0.5_8*s*log(2.0_8*Math_PI) &
! 				- (3.0_8+this.fl())*log(2.0_8*Math_PI) &
				- log( Gamma(0.5_8*s) ) &
				+ 1.5*logMu &
				+ this.logVfree_ + this.logVtheta_ &
				+ 0.5_8*fr*log(2.0_8) &
				+ 0.5_8*this.LnDiagI_ &
				+ (0.5_8*s-1.0_8)*log(Et)
				
! 			write(*,"(A,4I5,A,F10.5)") "s = ", s, this.ft(), this.fl(), fr, "   "//trim(this.label()), this.LnLambda
		end if
		
		if( ( Math_isNaN(this.LnLambda) .or. Math_isInf(this.LnLambda) ) .and. GOptions_printLevel >= 2 ) then
			write(*,*) ""
			call GOptions_valueReport( "reactorEnergy", this.reactorEnergy()/eV, "eV", indent=2 )
			call GOptions_valueReport( "internalEnergy", this.internalEnergy()/eV, "eV", indent=2 )
			call GOptions_valueReport( "kineticEnergy", this.kineticEnergy()/eV, "eV", indent=2 )
			call GOptions_valueReport( "vibrationalEnergy", this.vibrationalEnergy_/eV, "eV", indent=2 )
			call GOptions_valueReport( "rotationalEnergy", this.rotationalEnergy/eV, "eV", indent=2 )
			call GOptions_valueReport( "used Et", Et/eV, "eV", indent=2 )
			call GOptions_valueReport( "logVfree", this.logVfree_, indent=2 )
			call GOptions_valueReport( "logVtheta", this.logVtheta_, indent=2 )
			call GOptions_valueReport( "logVJ", this.logVJ_, indent=2 )
			call GOptions_valueReport( "0.5*logMu", 0.5*logMu, indent=2 )
			call GOptions_valueReport( "LnLambda", this.LnLambda, indent=2 )
		end if
	end subroutine updateLambda
		
	!>
	!! @brief Inverse of the temperature
	!!
	function iTemperature( this ) result( output )
		class(FragmentsList), intent(in) :: this
		real(8) :: output
		
		real(8) :: Et, Ed
		
		integer :: fr, i
		
		! @todo Hay que hacer una función que mantenga fr en memoria, lo importante es
		!       aclarar que solo tiene n-1 contribuciones
		fr = 0
		do i=1,this.nMolecules()-1
			fr = fr + this.clusters( this.idSorted(i) ).fr()
		end do
		
		Et = this.kineticEnergy()
		
		output = 0.5_8*( real( 3*this.nMolecules()-3-2 , 8 )/Et + fr/Ed  )
	end function iTemperature

	!>
	!! @brief Test method
	!!
	subroutine FragmentsList_test()
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
! 		write(*,*) "radius = ", cluster.radius()
! 		call clist.set( 1, cluster )
! 		
! 		write(*,*) "================================"
! 		fstr = "    tlC2    0  T    3  1    2   C2T-linear.xyz     -2061.703744  0.209800       2*tC1    300"
! 		
! 		call cluster.init( fstr, 2 )
! 		call cluster.show()
! 		write(*,*) "radius = ", cluster.radius()
! 		call clist.set( 2, cluster )
! 		write(*,*) "================================"
! 		fstr = "    tcC4    0  F    3  0    4   C4T-cyclic.xyz     -4129.378872  0.089866    tC1,slC3    300"
! 		
! 		call cluster.init( fstr, 3 )
! 		call cluster.show()
! 		write(*,*) "radius = ", cluster.radius()
! 		call clist.set( 3, cluster )
! 		write(*,*) "================================"
! 		fstr = "    tcC5    0  F    3  0    1   C5T-cyclic.xyz     -5162.267819  0.093690   slC2,slC3    300"
! 		
! 		call cluster.init( fstr, 4 )
! 		call cluster.show()
! 		write(*,*) "radius = ", cluster.radius()
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
		type(FragmentsList) :: clist
		
		allocate( massTable(3) )
		massTable(1) = "       sC1    0  1  2    0   C1.xyz             -1026.581828"
		massTable(2) = "      slC2    0  1  0    2   C2S-linear.xyz     -2062.282893       2*sC1   300"
		massTable(3) = "      slC3    0  1  0    2   C3S-linear.xyz     -3097.388207    sC1,slC2   400"
		
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
! 		call RigidMoleculeDatabase_instance.init( massTable )
		
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
! 				write(*,*) i, trim(adjustl(items(j))), RigidMoleculeDatabase_instance.getIdFromName( items(j) )
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
		call clist.loadXYZ( "initialGeom.xyz" )
		call clist.save()
		
	end subroutine FragmentsList_test
	
end module FragmentsList_

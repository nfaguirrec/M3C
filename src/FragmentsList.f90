!>
!! @brief
!!
module FragmentsList_
	use String_
	use Math_
	use Matrix_
	use SpecialMatrix_
	use RandomSampler_
	
	use GOptions_
	use FragmentsListBase_
	
	implicit none
	private
	
	public :: &
		FragmentsList_test
	
	type, public, extends( FragmentsListBase ) :: FragmentsList
		real(8) :: rotationalEnergy      !< It is calculated into updateRotationalEnergy procedure
		
! 		real(8), private :: LnWt         !< Translational weight
! 		real(8), private :: LnWr         !< Rotational diagonal weight
		real(8) :: LnWt         !< Translational weight
		real(8) :: LnWr         !< Rotational diagonal weight
		real(8), private :: LnIm_         !< Contiene el log del producto de la diagonal de los tensores de inercia efectivos
		                                  !< Se calculan al corregir la energía rotacional, así que aprovecho este lugar para
		                                  !< obtenerlos, ya que su calculo es probablemente lo maś costoso
		
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
			procedure, private :: updateLnWt
! 			procedure :: showLnWComponents
			
			procedure, private :: updateRotationalEnergy
			procedure, private :: updateRotationalEnergyJnFull
			procedure, private :: updateRotationalEnergyJnFull2
			procedure, private :: updateRotationalEnergyJn
			
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
		this.LnWt = 0.0_8
		this.LnWr = 0.0_8
		this.LnIm_ = 0.0_8
	end subroutine initFragmentsList
	
	!>
	!! @brief Copy constructor
	!!
	subroutine copyFragmentsList( this, other )
		class(FragmentsList), intent(inout) :: this
		class(FragmentsList), intent(in) :: other
		
		call this.copyFragmentsListBase( other )
		
		this.rotationalEnergy = other.rotationalEnergy
		this.LnWt = other.LnWt
		this.LnWr = other.LnWr
		this.LnIm_ = other.LnIm_
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
! 		call this.updateLnWt()
		
! 		call this.setFrozen( .false. )

		call this.changeVibrationalEnergy()
		call this.changeOrientations()
		call this.changeGeometry()
		call this.updateRotationalEnergy()
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_section( "END BUILDING INITIAL CONFIGURATION "//trim(this.label()), indent=2 )
		end if
		
		this.forceInitializing = .false.
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
		call this.updateRotationalEnergy()
		call this.updateKineticEnergy()
		call this.updateLnWt()
		
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
			call GOptions_section( "CHANGE VIBRATIONAL ENERGY"//trim(this.label()), indent=2 )
			write(STDOUT,*) ""
		end if

! 		call this.changeVibrationalEnergyFragmentsListBase()
! 		call this.updateKineticEnergy()
! 		call this.updateLnWt()

		call this.changeVibrationalEnergyFragmentsListBase()
		call this.updateRotationalEnergy()
		call this.updateKineticEnergy()
		call this.updateLnWt()
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_section( "END CHANGE VIBRATIONAL ENERGY"//trim(this.label()), indent=2 )
		end if
	end subroutine changeVibrationalEnergy
	
	!>
	!! @brief 
	!!
	subroutine changeOrientations( this )
		class(FragmentsList) :: this
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_subsection( "Random orientations "//trim(this.label()), indent=2 )
		end if
		
		call this.changeOrientationsFragmentsListBase()
		call this.updateRotationalEnergy()
	end subroutine changeOrientations
	
	!>
	!! @brief Total weight
	!!
	function LnW( this ) result( output )
		class(FragmentsList) :: this
		real(8) :: output
		
		output = this.LnWe() + this.LnWv() + this.LnWn() +  this.LnWr + this.LnWt
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
! 			trim(FString_fromReal(this.LnWr,"(F10.5)"))// &
! 			trim(FString_fromReal(this.LnWt,"(F10.5)"))// &
! 			trim(FString_fromReal(this.LnWe()+this.LnWv()+this.LnWn()+this.LnWt,"(F10.5)"))// &
! 			"     "//trim(this.label())
! 	end subroutine showLnWComponents
	
	!>
	!! @brief Returns the electronic energy
	!!
	function electronicEnergy( this ) result( output )
		class(FragmentsList) :: this
		real(8) :: output
		
		integer :: i
		
		output = 0.0_8
		do i=1,this.nMolecules()
			output = output + this.clusters(i).electronicEnergy
		end do
	end function electronicEnergy
	
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
		
! 		output = this.vibrationalEnergy_  + this.rotationalEnergy + this.intermolEnergy_
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
			this.LnWn(), this.LnWr, this.LnWt, &
			this.LnW(), trim(this.label())
			
		output = line
	end function weightHistoryLine
	
	!>
	!! @brief Actualiza la energía rotacional
	!!
	subroutine updateRotationalEnergy( this )
		class(FragmentsList) :: this
		
! 		if( GOptions_useLReference ) then
! 			call updateRotationalEnergyL( this )
! 		else
! 			call updateRotationalEnergyJn( this )
			call updateRotationalEnergyJnFull( this )
! 			call this.updateRotationalEnergyJnFull2()
! 		end if
		
	end subroutine updateRotationalEnergy
	
	!>
	!! @brief Actualiza la energía rotacional
	!!
	subroutine updateRotationalEnergyJnFull( this )
		class(FragmentsList) :: this
		
		real(8) :: maxErot
		integer :: i, j, n
		integer :: mu, nu, effMu, effNu
		type(Matrix) :: In, invIn    !< Tensor de inercia de N proyectado sobre los ejes de i y su inversa
		type(Matrix) :: Ii, invIi    !< Tensor de inercia de i proyectado sobre los ejes de i (diagonal) y su inversa
		type(Matrix) :: Il, invIl
		type(Matrix) :: It, invIt      !< Tensor de inercia de i proyectado sobre los ejes de i (diagonal) y su inversa
		type(Matrix) :: Ji, L, Jn         !< Momentos angulares bf de i y n respectivamente
		type(Matrix) :: RotMu, RotNu, RotN      !< Matrices de rotación
		type(Matrix) :: invBigI
		real(8) :: weight
		real(8) :: maxIval
		
		type(Matrix) :: Ui
		type(Matrix) :: Bi, invBi, tmp
		integer :: fr_sf
		
		real(8) :: Ein
		real(8) :: randNumber
		real(8) :: randDirection
		real(8) :: tmpErot
		
		logical :: debug = .false.
		
		!-----------------------------------------------------------------------
		! Valores para el apaño de 1 átomo + 1 molecula en el momento orbital
		logical :: candidateLCorrection
		real(8) :: Re, muMass, r, m, I_L
		type(RandomSampler) :: rs
		real(8), allocatable :: sample(:,:)
		integer :: trials
		real(8) :: ssum
		real(8) :: l1, l2
		real(8) :: j1, j2
		integer :: fr
		real(8) :: varIn
		real(8) :: varJn
		
		real(8) :: Ilcorr
		!-----------------------------------------------------------------------
		
		if( this.forceInitializing ) then
			call this.initialGuessFragmentsList()
			return
		end if
		
		n = this.nMolecules()
		
		call Jn.columnVector( 3, val=0.0_8 )
		this.L_ = 0.0_8
		this.rotationalEnergy = 0.0_8
		this.LnIm_ = 0.0_8
		
		do i=1,n
			this.clusters(i).J_ = 0.0_8
		end do
		
		maxErot = this.reactorEnergy()-this.vibrationalEnergy_-this.intermolEnergy_
		
		if( maxErot < 0.0_8 .or. n < 2 ) return
		
! 		! Si todos son fragmentos atómicos excepto el último (el último también puede ser atómico)
! 		candidateLCorrection = .true.
! 		do i=1,n-1
! 			if( this.clusters( this.idSorted(i) ).fr() /= 0 ) then
! 				candidateLCorrection = .false.
! 				exit
! 			end if
! 		end do
! 		
! 		if( candidateLCorrection .and. GOptions_useLCorrection ) then
! 				
! 				Re = this.clusters( this.idSorted(n) ).radius()  ! El valor exacto es radius-rcov(1)-rcov(n), el valor que está es para una diatómica
! 				mu = 0.5*this.clusters( this.idSorted(n) ).mass()
! 				
! 				weight = 0.0_8
! 				do i=1,n-1
! 					m = 1.0_8/( 1.0_8/this.clusters( this.idSorted(i) ).mass() + 1.0_8/this.clusters( this.idSorted(n) ).mass() )
! 					r = norm2( this.clusters( this.idSorted(i) ).center()-this.clusters( this.idSorted(n) ).center() )
! 						
! 					if( this.clusters( this.idSorted(n) ).fr() /= 0 ) then
! 						I_L = 1.0_8/( 1.0_8/( mu*Re**2  ) + 1.0_8/( m*r**2 ) )
! 					else
! 						I_L = m*r**2
! 					end if
! 					
! 					weight = weight + GOptions_gammaLCorrection*this.clusters( this.idSorted(n) ).fr()*log(2.0_8*I_L)
! 				end do
! 				
! 				this.LnIm_ = weight
! 					
! 				if( GOptions_printLevel >= 3 ) then
! 					call GOptions_valueReport( "1) LnSqrt2Im", this.LnIm_, indent=2 )
! 				end if
! 				
! 		end if
		
		! Contribución orbital para el caso de 1 atomo y una molecula
		if( n==2 .and. GOptions_useLCorrection ) then
			
			weight = 0.0_8
			
! 			if( this.clusters( this.idSorted(1) ).fr() == 0 .and. this.clusters( this.idSorted(2) ).fr() /= 0 ) then
				! Simplemente he utilizado como momento de inercia el valor asociado al momento orbital del átomo
				! alrededor de la molécula
				varIn = this.clusters( this.idSorted(1) ).mass()*this.clusters( this.idSorted(2) ).mass() &
							*norm2(this.clusters(1).center()-this.clusters(2).center())**2/this.mass()
				r = norm2(this.clusters(1).center()-this.clusters(2).center())
							
				! El factor sqrt, está considerado cuando se convierte LnIm_ a LnWr
				weight = weight + GOptions_gammaLCorrection*maxval( [this.clusters( this.idSorted(1) ).fr(), this.clusters( this.idSorted(2) ).fr()] )*log(2.0_8*varIn)
! 				weight = weight + GOptions_gammaLCorrection*this.clusters( this.idSorted(2) ).fr()*log(varIn/this.clusters( this.idSorted(1) ).mass())
! 				weight = weight + max( 0, this.clusters( this.idSorted(2) ).fr()-1 )*log(r**2)
				
				call random_number(varJn)
				varJn = -sqrt(2.0_8*this.kineticEnergy()*varIn) + 2.0_8*varJn*sqrt(2.0_8*this.kineticEnergy()*varIn)
				this.rotationalEnergy = varJn**2/(2.0_8*varIn)
! 			end if
			
			this.LnIm_ = weight
			
			if( GOptions_printLevel >= 3 ) then
				write(*,"(A)") trim(this.label())
				write(*,"(A,2F10.5)") "L-correction.weight = ", weight, log(2.0_8*varIn)
				write(*,"(A,F10.5)") "L-correction.Erot = ", this.rotationalEnergy
			end if
			
			return
		end if

! 		! Contribución orbital L-correction
! 		! version:  Mon Jul 13 19:23:09 CEST 2015
! 		if( n==2 .and. GOptions_useLCorrection ) then
! 			
! 			if( this.clusters( this.idSorted(1) ).fr() > this.clusters( this.idSorted(2) ).fr() ) then
! 				mu = 1
! 			else
! 				mu = 2
! 			end if
! 			
! 			r = norm2(this.clusters(1).center()-this.clusters(2).center())
! 			Ilcorr = this.clusters( this.idSorted(1) ).mass()*this.clusters( this.idSorted(2) ).mass()*r**2/this.mass()
! 			
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			! Calculo del peso
! 			call invIi.init( 3, 3, 0.0_8 )
! 			do j=3,4-this.clusters( this.idSorted(mu) ).fr(),-1
! 				invIi.data( j, j ) = 1.0_8/this.clusters( this.idSorted(mu) ).diagInertiaTensor.data( j, j )
! 			end do
! 			
! 			call invIi.eigen( eVals=invBi, eVecs=Ui )
! 						
! 			this.LnIm_ = 0.0_8
! 			do j=3,4-this.clusters( this.idSorted(mu) ).fr(),-1
! 				this.LnIm_ = this.LnIm_ - log(2.0_8*invBi.get(j,j))
! 			end do
! 			this.LnIm_ = this.LnIm_ + log(Math_PI**(this.clusters( this.idSorted(mu) ).fr()/2.0)) - log(gamma(this.clusters( this.idSorted(mu) ).fr()/2.0))
! 			
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			! Calculo de la corrección
! 			do j=3,4-this.clusters( this.idSorted(mu) ).fr(),-1
! 				this.LnIm_ = this.LnIm_ - log( 1.0_8 + 1.0_8/( invBi.get(j,j)*Ilcorr ) )
! 			end do
! 			
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			! Calculo de la energía rotacional y de J
! 			this.rotationalEnergy = 0.0_8
! 			do i=1,3
! 				if( ( .not. Math_isInf( 1.0_8/invBi.get(i,i) ) ) .and. 1.0_8/invBi.get(i,i) > maxIval ) then
! 					call random_number( randNumber ) ! [0,1]
! 					call random_number( randDirection ) ! [0,1]
! 					randDirection = 0.0_8+2.0_8*randDirection  ! [0:2]
! 					
! 					this.clusters( this.idSorted(mu) ).J_(i) = &
! 						(-1.0_8)**int(randDirection)*sqrt( 2.0_8*maxErot/abs( invBi.get(i,i) ) )*randNumber
! 					
! 					this.rotationalEnergy = this.rotationalEnergy + 0.5_8*this.clusters( this.idSorted(mu) ).J_(i)**2*invBi.get(i,i)
! 				end if
! 			end do
! 			
! 			if( GOptions_printLevel >= 3 ) then
! 				write(*,"(A)") trim(this.label())
! 				write(*,"(A,2F10.5)") "L-correction.weight = ", weight, log(2.0_8*varIn)
! 				write(*,"(A,F10.5)") "L-correction.Erot = ", this.rotationalEnergy
! 			end if
! 			
! 			return
! 		
! 		end if
		
		! Contribución orbital L-correction
		! version:  Thu Jul 16 12:55:44 CEST 2015
! 		if( GOptions_useLCorrection ) then
! 			
! 			mu = this.nFragments()
! 			
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			! Calculo del peso
! 			call invIi.init( 3, 3, 0.0_8 )
! 			do j=3,4-this.clusters( this.idSorted(mu) ).fr(),-1
! 				invIi.data( j, j ) = 1.0_8/this.clusters( this.idSorted(mu) ).diagInertiaTensor.data( j, j )
! 			end do
! 			
! 			call invIi.eigen( eVals=invBi, eVecs=Ui )
! 						
! 			this.LnIm_ = 0.0_8
! 			do j=3,4-this.clusters( this.idSorted(mu) ).fr(),-1
! 				this.LnIm_ = this.LnIm_ - log(2.0_8*invBi.get(j,j))
! 			end do
! 			this.LnIm_ = this.LnIm_ + log(Math_PI**(this.clusters( this.idSorted(mu) ).fr()/2.0)) - log(gamma(this.clusters( this.idSorted(mu) ).fr()/2.0))
! 			
! 		end if
		
! 		! Contribución orbital para el caso de 1 atomo y una molecula
! 		else if( n==3 .and. GOptions_useLCorrection ) then
! 			
! 			weight = 0.0_8
! 			
! 			if( this.clusters( this.idSorted(1) ).fr() == 0 .and. this.clusters( this.idSorted(2) ).fr() /= 0 ) then
! 				! Simplemente he utilizado como momento de inercia el valor asociado al momento orbital del átomo
! 				! alrededor de la molécula
! 				varIn = this.clusters( this.idSorted(1) ).mass()*this.clusters( this.idSorted(2) ).mass() &
! 							*norm2(this.clusters(1).center()-this.clusters(2).center())**2/this.mass()
! 				r = norm2(this.clusters(1).center()-this.clusters(2).center())
! 							
! 				! El factor sqrt, está considerado cuando se convierte LnIm_ a LnWr
! 				weight = weight + GOptions_gammaLCorrection*maxval( [this.clusters( this.idSorted(1) ).fr(), this.clusters( this.idSorted(2) ).fr()] )*log(2.0_8*varIn)
! ! 				weight = weight + GOptions_gammaLCorrection*this.clusters( this.idSorted(2) ).fr()*log(varIn/this.clusters( this.idSorted(1) ).mass())
! ! 				weight = weight + max( 0, this.clusters( this.idSorted(2) ).fr()-1 )*log(r**2)
! 				
! 				call random_number(varJn)
! 				varJn = -sqrt(2.0_8*this.kineticEnergy()*varIn) + 2.0_8*varJn*sqrt(2.0_8*this.kineticEnergy()*varIn)
! 				this.rotationalEnergy = varJn**2/(2.0_8*varIn)
! 			end if
! 			
! 			this.LnIm_ = weight
! 			
! 			if( GOptions_printLevel >= 3 ) then
! 				write(*,"(A)") trim(this.label())
! 				write(*,"(A,2F10.5)") "L-correction.weight = ", weight, log(2.0_8*varIn)
! 				write(*,"(A,F10.5)") "L-correction.Erot = ", this.rotationalEnergy
! 			end if
! 			
! 			return
! 		end if
		
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
		
		if( debug ) then
			write(*,*) trim(this.label())
		end if
		
		if( effNu > 1 .and. effMu /= effNu ) then
			write(*,"(A,2I8)") "FragmentsList.updateRotationalEnergyJnFull(). effMu /= effNu, ", effMu, effNu
			stop 
		end if
		
! 		if( effMu < 1 .and. GOptions_useLCorrection ) return
		if( effMu < 1 ) return
		
		if( GOptions_useLDOSContrib ) then
			call invBigI.init( (effMu+1)*3, (effNu+1)*3, 0.0_8 )
		else
			call invBigI.init( (effMu)*3, (effNu)*3, 0.0_8 )  !<< Se utiliza este cuando no se mete el L
		end if
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
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
		! Contribución de los N-1 fragmentos a la energía rotacional
		effMu = 1
		do mu=1,n-1
			
			if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
				
				if( debug ) then
					write(*,*) "(mu,effMu) = ", mu, effMu
				end if
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Matriz de rotación que permite transformar a los ejes
				! de N en los ejes de i, pasando por el sistem fix
				! RotMu = R_n^T*R_mu
				RotMu = SpecialMatrix_rotationTransform( &
					this.clusters( this.idSorted(mu) ).inertiaAxes(), &
					this.clusters( this.idSorted(n) ).inertiaAxes() )
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				effNu = 1
				do nu=1,n-1
					
					if( debug ) then
						write(*,*) "(nu,effNu) = ", nu, effNu
					end if
					
					if( this.clusters( this.idSorted(nu) ).fr() /= 0 ) then
					
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! Matriz de rotación que permite transformar a los ejes
						! de N en los ejes de i, pasando por el sistem fix
						! RotNu = R_n^T*R_nu
						RotNu = SpecialMatrix_rotationTransform( &
							this.clusters( this.idSorted(nu) ).inertiaAxes(), &
							this.clusters( this.idSorted(n) ).inertiaAxes() )
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						
						if( mu == nu ) then
						
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							! Se obtiene la matriz inversa del tensor
							! de inercia de i
							call invIi.init( 3, 3, 0.0_8 )
							do j=3,4-this.clusters( this.idSorted(mu) ).fr(),-1
								invIi.data( j, j ) = 1.0_8/this.clusters( this.idSorted(mu) ).diagInertiaTensor.data( j, j )
							end do
							
							if( debug ) then
								write(*,*) "Base(", mu, ",",  nu, ")"
								call invIi.show( formatted=.true., precision=8 )
! 								write(*,*) "Base(", mu, ",",  nu, ") inverse"
! 								Ii = invIi.inverse()
! 								call Ii.show( formatted=.true., precision=8 )
							end if
							
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							! Se obtiene la matriz inversa del tensor
							! de inercia efectivo.
							invIt =  invIi + RotMu*invIn*RotMu.transpose()
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							
							if( debug ) then
								write(*,*) "Effetive to Diag (", effMu, ",",  effNu, ")"
								call invIt.show( formatted=.true., precision=8 )
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
							
							if( debug ) then
								write(*,*) "Extra diag(", effMu, ",",  effNu, ")"
								call invIt.show( formatted=.true., precision=8 )
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
		
		if( GOptions_useLDOSContrib ) then
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Matriz de rotación que permite transformar a los ejes
			! de N en los ejes de i, pasando por el sistem fix
			! Rot = R_n
			RotN = SpecialMatrix_rotationTransform( &
				this.inertiaAxes, &
				this.clusters( this.idSorted(n) ).inertiaAxes() )
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Contribución de L a la energía rotacional
			effMu = 1
			do mu=1,n
				
				if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
				
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					! Matriz de rotación que permite transformar a los ejes
					! de N en los ejes de i, pasando por el sistem fix
					! RotMu = R_n^T*R_mu
					RotMu = SpecialMatrix_rotationTransform( &
						this.clusters( this.idSorted(mu) ).inertiaAxes(), &
						this.clusters( this.idSorted(n) ).inertiaAxes() )
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					
					if( mu == n ) then
					
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! Se obtiene la matriz inversa del tensor
						! de inercia de Il
						call invIl.init( 3, 3, 0.0_8 )
						do j=3,4-this.fl(),-1
							invIl.data( effMu, j ) = 1.0_8/this.diagInertiaTensor( j )
						end do

						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! Se obtiene la matriz inversa del tensor
						! de inercia efectivo.
						invIt =  invIl + RotN*invIn*RotN.transpose()
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						
						if( debug ) then
							write(*,*) "DiagL (", effMu, ")"
							call invIt.show( formatted=.true., precision=8 )
						end if
						
						do i=1,3
							do j=1,3
								call invBigI.set( 3*(effMu-1)+i, 3*(effNu-1)+j, invIt.get(i,j) )
							end do
						end do
						
					else
					
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! Se obtiene la matriz de acoplamieto
						invIt =  RotMu*invIn*RotN.transpose()
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						
						if( debug ) then
							write(*,*) "Extra diagL (", effMu, ")"
							call invIt.show( formatted=.true., precision=8 )
						end if
						
						do i=1,3
							do j=1,3
								call invBigI.set( 3*(effMu-1)+i, 3*(effNu-1)+j, invIt.get(i,j) )
								call invBigI.set( 3*(effNu-1)+i, 3*(effMu-1)+j, invIt.get(i,j) )
							end do
						end do
					
					end if
					
					effMu = effMu + 1
				end if
				
			end do
		end if
		
		if( debug ) then
			write(*,*) "Total"
			call invBigI.show( formatted=.true., precision=8 )
		end if
		
		maxIval = 0.0_8
		do effMu=1,invBigI.nRows
			if( ( .not. Math_isInf( 1.0_8/invBigI.get(effMu,effMu) ) ) .and. 1.0_8/invBigI.get(effMu,effMu) > maxIval ) then
				maxIval = 1.0_8/invBigI.get(effMu,effMu)
			end if
		end do
		
		call invBigI.eigen( eVals=invBi, eVecs=Ui )
		
		if( debug ) then
			write(*,*) "maxIval = ", maxIval
			write(*,*) "Diag Total"
			call invBi.show( formatted=.true., precision=8 )
			Bi = invBi.inverse()
			write(*,*) "Diag Total inverse"
			call Bi.show( formatted=.true., precision=8 )
		end if
		
! 		tol = 0.8*invBi.trace()/invBi.nRows
! 		
! 		if( debug ) then
! 			write(*,*) "Tol = ", tol
! 		end if
		
		fr_sf = 0
		do effMu=1,invBi.nRows
			if( 1.0_8/invBi.get(effMu,effMu) <= maxIval ) then  !< Curiosamente para el caso = funciona
				fr_sf = fr_sf + 1
			end if
		end do
		
		if( debug ) then
			write(*,*) "fr_sf = ", fr_sf
		end if
		
		weight = 0.0_8
		do effMu=invBi.nRows,invBi.nRows-fr_sf+1,-1
			if( debug ) then
				write(*,*) "effMu = ", effMu, invBi.get(effMu,effMu)
			end if
			
			weight = weight - log(invBi.get(effMu,effMu))
		end do
		
		if( GOptions_useLWeightContrib ) then
			do j=3,4-this.fl(),-1
				weight = weight - log(this.diagInertiaTensor(j))
! 				write(*,*) "I", j, " = ", - log(this.diagInertiaTensor(j))
			end do
		end if
		
		if( debug ) then
			write(*,*) "weight = ", weight
		end if
		
		this.LnIm_ = weight
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Genera los momentos angulares en sf
		call Jn.columnVector( 3, val=0.0_8 )
		
		Ein = 0.0_8
		effMu = 1  !< effMu da las coordenadas para localizar un valor en la matriz bigI
		do mu=1,n-1
			if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
			
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Matriz de rotación que permite transformar a los ejes
				! de N en los ejes de i, pasando por el sistem fix
				! RotMu = R_n^T*R_mu
				RotMu = SpecialMatrix_rotationTransform( &
					this.clusters( this.idSorted(effMu) ).inertiaAxes(), &
					this.clusters( this.idSorted(n) ).inertiaAxes() )
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				do while( .true. )
					tmpErot = 0.0_8
					do i=1,3
						if( ( .not. Math_isInf( 1.0_8/invBigI.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) ) &
								.and. 1.0_8/invBigI.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) > maxIval ) then
! 						if( ( invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) > 0.0_8 .and. abs( invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) > tol ) then
							call random_number( randNumber ) ! [0,1]
							call random_number( randDirection ) ! [0,1]
							randDirection = 0.0_8+2.0_8*randDirection  ! [0:2]
							
							this.clusters( this.idSorted(mu) ).J_(i) = &
								(-1.0_8)**int(randDirection)*sqrt( 2.0_8*maxErot/abs( invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) )*randNumber
							
	! 						write(*,*) "mu, i, J = ", mu, i, this.clusters( this.idSorted(effMu) ).J_(i)
							tmpErot = tmpErot + 0.5_8*this.clusters( this.idSorted(mu) ).J_(i)**2*invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i )
							
	! 						write(*,*) "Ei = ", Ein
						end if
					end do
					
					if( tmpErot <= maxErot ) exit
				end do
				Ein = tmpErot
				
				if( debug ) then
					write(*,"(A,A,3F10.1)") trim(this.clusters( this.idSorted(mu) ).label()), ", J = ", this.clusters( this.idSorted(mu) ).J_
				end if
				
				! @todo Hace falta una linea transformadolo J a bf (La energia si se calcula en sf, tal como esta arriba)
				
				! Momento angular para el fragmento mu-esimo
				call Ji.columnVector( 3, values=this.clusters( this.idSorted(mu) ).J_ )
				
				! Contribución a Jn
				Jn = Jn - RotMu*Ji
				
				effMu = effMu + 1
			end if
		end do
		
		this.clusters(this.idSorted(n)).J_ = Jn.data(:,1)
		
		this.rotationalEnergy = Ein
		
		if( debug ) then
			write(*,"(A,A,3F10.1)") trim(this.clusters( this.idSorted(n) ).label()), ", Jcal = ", this.clusters( this.idSorted(n) ).J_
			write(*,"(A,F10.5)") "Energy = ", Ein
		end if
		
		if( debug ) then
			stop
		end if
		
	end subroutine updateRotationalEnergyJnFull
	
	!>
	!! @brief 
	!!
	function buildInertiaRL( r, mass, center ) result( Im )
		real(8), intent(in) :: r(3)
		real(8), intent(in) :: mass
		real(8), optional, intent(in) :: center(3)
		type(Matrix) :: Im
		
		real(8) :: X, Y, Z, m
		
		real(8) :: effCenter(3)
		
		effCenter = 0.0_8
		if( present(center) ) effCenter = center
		
		X = r(1) - effCenter(1)
		Y = r(2) - effCenter(2)
		Z = r(3) - effCenter(3)
		m = mass
		
		call Im.init(3,3)
		
		call Im.set( 1, 1,  m*(Y**2+Z**2) )
		call Im.set( 1, 2, -m*X*Y )
		call Im.set( 1, 3, -m*X*Z )
		call Im.set( 2, 1, -m*Y*X )
		call Im.set( 2, 2,  m*(X**2+Z**2) )
		call Im.set( 2, 3, -m*Y*Z )
		call Im.set( 3, 1, -m*Z*X )
		call Im.set( 3, 2, -m*Z*Y )
		call Im.set( 3, 3,  m*(X**2+Y**2) )
	end function buildInertiaRL
	
	!>
	!! @brief 
	!!
	function buildinvInertiaRL( r, mass, center ) result( Im )
		real(8), intent(in) :: r(3)
		real(8), intent(in) :: mass
		real(8), optional, intent(in) :: center(3)
		type(Matrix) :: Im
		
		real(8) :: X, Y, Z, m
		
		real(8) :: effCenter(3), det
		
		effCenter = 0.0_8
		if( present(center) ) effCenter = center
		
		X = r(1) - effCenter(1)
		Y = r(2) - effCenter(2)
		Z = r(3) - effCenter(3)
		m = mass
		
		det = m*(X**2)*m*(Y**2)-m*X*Y*m*Y*X
		
		call Im.init( 3, 3, 0.0_8 )
		
		call Im.set( 1, 1,  m*(X**2)/det )
		call Im.set( 1, 2,  m*X*Y/det )
		call Im.set( 2, 1,  m*Y*X/det )
		call Im.set( 2, 2,  m*(Y**2)/det )
		call Im.set( 3, 3,  1.0_8/m*(X**2+Y**2) )
	end function buildinvInertiaRL
	
	!>
	!! @brief Actualiza la energía rotacional
	!!
	subroutine updateRotationalEnergyJnFull2( this )
		class(FragmentsList) :: this
		
		real(8) :: maxErot
		integer :: i, j, n
		integer :: mu, nu, effMu, effNu
		type(Matrix) :: In, invIn    !< Tensor de inercia de N proyectado sobre los ejes de i y su inversa
		type(Matrix) :: Ii, invIi    !< Tensor de inercia de i proyectado sobre los ejes de i (diagonal) y su inversa
		type(Matrix) :: Il, invIl
		type(Matrix) :: It, invIt      !< Tensor de inercia de i proyectado sobre los ejes de i (diagonal) y su inversa
		type(Matrix) :: Ji, L, Jn         !< Momentos angulares bf de i y n respectivamente
		type(Matrix) :: RotMu, RotNu, RotN      !< Matrices de rotación
		type(Matrix) :: invBigI
		real(8) :: weight
		real(8) :: maxIval
		
		type(Matrix) :: Ui
		type(Matrix) :: Bi, invBi, tmp
		integer :: fr_sf
		
		real(8) :: Ein
		real(8) :: randNumber
		real(8) :: randDirection
		real(8) :: tmpErot
		
		logical :: debug = .true.
		
		integer :: nEffMu, nEffNu
		type(Matrix) :: Ilmu, invIlmu, Ilnu
		real(8) :: det
		
		if( this.forceInitializing ) then
			call this.initialGuessFragmentsList()
			return
		end if
		
		n = this.nMolecules()
		
		call Jn.columnVector( 3, val=0.0_8 )
		this.L_ = 0.0_8
		this.rotationalEnergy = 0.0_8
		this.LnIm_ = 0.0_8
		
		do i=1,n
			this.clusters(i).J_ = 0.0_8
		end do
		
		maxErot = this.reactorEnergy()-this.vibrationalEnergy_-this.intermolEnergy_
		
		if( maxErot < 0.0_8 .or. n < 2 ) return
		
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
		
		nEffMu = effMu
		nEffNu = effNu
		
		if( debug ) then
			write(*,*) trim(this.label())
		end if
		
		if( effNu > 1 .and. effMu /= effNu ) then
			write(*,"(A,2I8)") "FragmentsList.updateRotationalEnergyJnFull(). effMu /= effNu, ", effMu, effNu
			stop 
		end if
		
! 		if( GOptions_useLDOSContrib ) then
! 			call invBigI.init( (effMu+1)*3, (effNu+1)*3, 0.0_8 )
! 		else
			call invBigI.init( 3*nEffMu+3*(n-1), 3*nEffNu+3*(n-1), 0.0_8 )  !<< Se utiliza este cuando no se mete el L
! 		end if
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
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
		effMu = 1
		do mu=1,n-1
		
			if( debug ) then
				write(*,*) "JJ' (mu,effMu) = ", mu, effMu
			end if
			
			if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Matriz de rotación que permite transformar a los ejes
				! de N en los ejes de i, pasando por el sistem fix
				! RotMu = R_n^T*R_mu
				RotMu = SpecialMatrix_rotationTransform( &
					this.clusters( this.idSorted(mu) ).inertiaAxes(), &
					this.clusters( this.idSorted(n) ).inertiaAxes() )
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				effNu = 1
				do nu=1,n-1
					
					if( debug ) then
						write(*,*) "JJ (nu,effNu) = ", nu, effNu
					end if
					
					if( this.clusters( this.idSorted(nu) ).fr() /= 0 ) then
					
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! Matriz de rotación que permite transformar a los ejes
						! de N en los ejes de i, pasando por el sistem fix
						! RotNu = R_n^T*R_nu
						RotNu = SpecialMatrix_rotationTransform( &
							this.clusters( this.idSorted(nu) ).inertiaAxes(), &
							this.clusters( this.idSorted(n) ).inertiaAxes() )
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						
						if( mu == nu ) then
						
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							! Se obtiene la matriz inversa del tensor
							! de inercia de i
							call invIi.init( 3, 3, 0.0_8 )
							do j=3,4-this.clusters( this.idSorted(mu) ).fr(),-1
								invIi.data( j, j ) = 1.0_8/this.clusters( this.idSorted(mu) ).diagInertiaTensor.data( j, j )
							end do
							
							if( debug ) then
								write(*,*) "Base(", mu, ",",  nu, ")"
								call invIi.show( formatted=.true., precision=8 )
! 								write(*,*) "Base(", mu, ",",  nu, ") inverse"
! 								Ii = invIi.inverse()
! 								call Ii.show( formatted=.true., precision=8 )
							end if
							
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							! Se obtiene la matriz inversa del tensor
							! de inercia efectivo.
							invIt =  invIi + RotMu*invIn*RotMu.transpose()
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							
							if( debug ) then
								write(*,*) "Effetive to Diag (", effMu, ",",  effNu, ")"
								call invIt.show( formatted=.true., precision=8 )
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
							
							if( debug ) then
								write(*,*) "Extra diag(", effMu, ",",  effNu, ")"
								call invIt.show( formatted=.true., precision=8 )
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
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Creación de los Il
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		effMu = 1
		do mu=1,n-1
			
			if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Matriz de rotación que permite transformar a los ejes
				! de N en los ejes de i, pasando por el sistem fix
				! RotMu = R_n^T*R_mu
				RotMu =  SpecialMatrix_rotationTransform( &
					   this.clusters( this.idSorted(mu) ).inertiaAxes(), &
					   this.clusters( this.idSorted(n) ).inertiaAxes() )
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				if( debug ) then
					write(*,*) "JL (mu,effMu) = ", mu, effMu
				end if
				
				invIt = RotMu*invIn
				
				if( debug ) then
					call invIt.show( formatted=.true., precision=8 )
				end if
				
				effNu = 1
				do nu=1,n-1
				
					do i=1,3
						do j=1,3
							call invBigI.set( (effMu-1)+i+3*(effMu-1), 3*(nu-1)+j+3*nEffNu, invIt.get(i,j) )
							call invBigI.set( 3*(nu-1)+j+3*nEffNu, (effMu-1)+i+3*(effMu-1), invIt.get(i,j) )
						end do
					end do
					
				end do
				
				effMu = effMu + 1
			end if
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Tensor de inercia orbital y su inversa
			Ilmu = buildInertiaRL( &
				this.clusters(this.idSorted(mu)).center(), &
				this.clusters( this.idSorted(mu) ).mass(), &
				this.clusters(this.idSorted(n)).center() )
				
			if( &
				abs(Ilmu.get(1,3)) < 1d-10 .and. &
				abs(Ilmu.get(2,3)) < 1d-10 .and. &
				abs(Ilmu.get(3,1)) < 1d-10 .and. &
				abs(Ilmu.get(3,2)) < 1d-10 &
			) then
! 				call invIlmu.init( 3, 3, 0.0_8 )
! 				
! 				det = Ilmu.get(1,1)*Ilmu.get(2,2)-Ilmu.get(1,2)*Ilmu.get(2,1)
! 				write(*,*) "Hola", det
! 				write(*,"(2F40.15)") Ilmu.get(1,1)*Ilmu.get(2,2), Ilmu.get(1,2)*Ilmu.get(2,1)
! 				call invIlmu.set( 1, 1,  Ilmu.get(2,2)/det )
! 				call invIlmu.set( 1, 2, -Ilmu.get(1,2)/det )
! 				call invIlmu.set( 2, 1, -Ilmu.get(2,1)/det )
! 				call invIlmu.set( 2, 2,  Ilmu.get(1,1)/det )
! 				
! 				call invIlmu.set( 3, 3,  1.0_8/Ilmu.get(3,3) )
				write(*,*) "Hola"
				invIlmu = buildinvInertiaRL( &
					this.clusters(this.idSorted(mu)).center(), &
					this.clusters( this.idSorted(mu) ).mass(), &
					this.clusters(this.idSorted(n)).center() )
			else
				invIlmu = Ilmu.inverse()
			end if
			
			if( debug ) then
				write(*,*) "Ilmu(", mu, ")"
				write(*,*) Ilmu.isDiagonal()
				call Ilmu.show( formatted=.true., precision=8 )
				write(*,*) "invIlmu(", mu, ")"
				call invIlmu.show( formatted=.true., precision=8 )
			end if
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			effNu = 1
			do nu=1,n-1
			
				if( debug ) then
					write(*,*) "LL (mu,nu) = ", mu, nu
				end if
				
				if( mu == nu ) then
					
					invIi = invIlmu + invIn
					
					if( debug ) then
						write(*,*) "Base(", mu, ",",  nu, ")"
						call invIi.show( formatted=.true., precision=8 )
					end if
					
				else
				
					invIi = invIn
					
					if( debug ) then
						write(*,*) "Extra diag(", mu, ",",  nu, ")"
						call invIi.show( formatted=.true., precision=8 )
					end if
					
				
				end if
				
				do i=1,3
					do j=1,3
						call invBigI.set( 3*(mu-1)+i+3*nEffMu, 3*(nu-1)+j+3*nEffNu, invIi.get(i,j) )
					end do
				end do
			end do
		end do
		
		if( debug ) then
			write(*,*) "Total"
			call invBigI.show( formatted=.true., precision=8 )
		end if
		
		!! @todo OJJJOOOO apaño
		invBigI = invBigI + 1d-12
		
		maxIval = 0.0_8
		do effMu=1,invBigI.nRows
			if( ( .not. Math_isInf( 1.0_8/invBigI.get(effMu,effMu) ) ) .and. 1.0_8/invBigI.get(effMu,effMu) > maxIval ) then
				maxIval = 1.0_8/invBigI.get(effMu,effMu)
			end if
		end do
		
		call invBigI.eigen( eVals=invBi, eVecs=Ui )
		
		if( debug ) then
			write(*,*) "maxIval = ", maxIval
			write(*,*) "Diag Total"
			call invBi.show( formatted=.true., precision=8 )
			Bi = invBi.inverse()
			write(*,*) "Diag Total inverse"
			call Bi.show( formatted=.true., precision=8 )
		end if
		
! 		tol = 0.8*invBi.trace()/invBi.nRows
! 		
! 		if( debug ) then
! 			write(*,*) "Tol = ", tol
! 		end if
		
		fr_sf = 0
		do effMu=1,invBi.nRows
			if( 1.0_8/invBi.get(effMu,effMu) <= maxIval ) then  !< Curiosamente para el caso = funciona
				fr_sf = fr_sf + 1
			end if
		end do
		
		if( debug ) then
			write(*,*) "fr_sf = ", fr_sf
		end if
		
		weight = 0.0_8
		do effMu=invBi.nRows,invBi.nRows-fr_sf+1,-1
			if( debug ) then
				write(*,*) "effMu = ", effMu, 1.0_8/invBi.get(effMu,effMu), -log(invBi.get(effMu,effMu))
			end if
			
			weight = weight - log(invBi.get(effMu,effMu))
		end do
		
		if( GOptions_useLWeightContrib ) then
			do j=3,4-this.fl(),-1
				weight = weight - log(this.diagInertiaTensor(j))
! 				write(*,*) "I", j, " = ", - log(this.diagInertiaTensor(j))
			end do
		end if
		
		if( debug ) then
			write(*,*) "weight = ", weight
		end if
		
		this.LnIm_ = weight
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Genera los momentos angulares en sf
		call Jn.columnVector( 3, val=0.0_8 )
		
		Ein = 0.0_8
		effMu = 1  !< effMu da las coordenadas para localizar un valor en la matriz bigI
		do mu=1,n-1
			if( this.clusters( this.idSorted(mu) ).fr() /= 0 ) then
			
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Matriz de rotación que permite transformar a los ejes
				! de N en los ejes de i, pasando por el sistem fix
				! RotMu = R_n^T*R_mu
				RotMu = SpecialMatrix_rotationTransform( &
					this.clusters( this.idSorted(effMu) ).inertiaAxes(), &
					this.clusters( this.idSorted(n) ).inertiaAxes() )
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				do while( .true. )
					tmpErot = 0.0_8
					do i=1,3
						if( ( .not. Math_isInf( 1.0_8/invBigI.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) ) &
								.and. 1.0_8/invBigI.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) > maxIval ) then
! 						if( ( invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) > 0.0_8 .and. abs( invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) > tol ) then
							call random_number( randNumber ) ! [0,1]
							call random_number( randDirection ) ! [0,1]
							randDirection = 0.0_8+2.0_8*randDirection  ! [0:2]
							
							this.clusters( this.idSorted(mu) ).J_(i) = &
								(-1.0_8)**int(randDirection)*sqrt( 2.0_8*maxErot/abs( invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i ) ) )*randNumber
							
	! 						write(*,*) "mu, i, J = ", mu, i, this.clusters( this.idSorted(effMu) ).J_(i)
							tmpErot = tmpErot + 0.5_8*this.clusters( this.idSorted(mu) ).J_(i)**2*invBi.get( 3*(effMu-1)+i, 3*(effMu-1)+i )
							
	! 						write(*,*) "Ei = ", Ein
						end if
					end do
					
					if( tmpErot <= maxErot ) exit
				end do
				Ein = tmpErot
				
				if( debug ) then
					write(*,"(A,A,3F10.1)") trim(this.clusters( this.idSorted(mu) ).label()), ", J = ", this.clusters( this.idSorted(mu) ).J_
				end if
				
				! @todo Hace falta una linea transformadolo J a bf (La energia si se calcula en sf, tal como esta arriba)
				
				! Momento angular para el fragmento mu-esimo
				call Ji.columnVector( 3, values=this.clusters( this.idSorted(mu) ).J_ )
				
				! Contribución a Jn
				Jn = Jn - RotMu*Ji
				
				effMu = effMu + 1
			end if
		end do
		
		this.clusters(this.idSorted(n)).J_ = Jn.data(:,1)
		
		this.rotationalEnergy = Ein
		
		if( debug ) then
			write(*,"(A,A,3F10.1)") trim(this.clusters( this.idSorted(n) ).label()), ", Jcal = ", this.clusters( this.idSorted(n) ).J_
			write(*,"(A,F10.5)") "Energy = ", Ein
			stop
		end if
		
		if( GOptions_printLevel >= 3 ) then
			call GOptions_valueReport( "Angular momentum coupling weight", weight, indent=2 )
		end if
		
	end subroutine updateRotationalEnergyJnFull2
	
	!>
	!! @brief Actualiza la energía rotacional
	!!
	subroutine updateRotationalEnergyJn( this )
		class(FragmentsList) :: this
		
		real(8) :: maxErot
		integer :: i, j, n
		type(Matrix) :: Inn, invInn    !< Tensor de inercia de N proyectado sobre los ejes de i y su inversa
		type(Matrix) :: Iii, invIii    !< Tensor de inercia de i proyectado sobre los ejes de i (diagonal) y su inversa
		type(Matrix) :: Il, invIl
		type(Matrix) :: It, invIt      !< Tensor de inercia de i proyectado sobre los ejes de i (diagonal) y su inversa
		type(Matrix) :: Ji, L, Jn         !< Momentos angulares bf de i y n respectivamente
		type(Matrix) :: Rot            !< Matriz de rotación
		
		type(Matrix) :: Xi
! 		type(Matrix) :: Ai, invAi
		type(Matrix) :: Ui
		type(Matrix) :: Bi, invBi, tmp
		real(8) :: scaleFactor
		integer :: fr_sf
		real(8) :: angles(3)
		
		real(8) :: varIn, varJn !< Valores para el apaño de 1 átomo + 1 molecula en el momento orbital
		
		real(8) :: Ein
! 		type(Matrix) :: Ein    !< Matriz auxiliar de 1x1, que representa la contribución a la energía rotacional por la interacción in
		real(8) :: randNumber
		
		type(Matrix) :: tmpMatrix
		
		!-----------------------------------------------------------------------
		! Valores para el apaño de 1 átomo + 1 molecula en el momento orbital
		logical :: candidateLCorrection
		real(8) :: Re, mu, r, m, I_L
		type(RandomSampler) :: rs
		real(8), allocatable :: sample(:,:)
		integer :: trials
		real(8) :: ssum
		real(8) :: l1, l2
		real(8) :: j1, j2
		integer :: fr
		!-----------------------------------------------------------------------
		
		if( this.forceInitializing ) then
			call this.initialGuessFragmentsList()
			return
		end if
		
		n = this.nMolecules()
		
		call Jn.columnVector( 3, val=0.0_8 )
		this.L_ = 0.0_8
		this.rotationalEnergy = 0.0_8
		this.LnIm_ = 0.0_8
		
		do i=1,n
			this.clusters(i).J_ = 0.0_8
		end do
		
		maxErot = this.reactorEnergy()-this.vibrationalEnergy_-this.intermolEnergy_
		
		if( maxErot < 0.0_8 .or. n < 2 ) return
		
		! Si todos son fragmentos atómicos excepto el último (el último también puede ser atómico)
		candidateLCorrection = .true.
		do i=1,n-1
			if( this.clusters( this.idSorted(i) ).fr() /= 0 ) then
				candidateLCorrection = .false.
				exit
			end if
		end do
		
		if( candidateLCorrection .and. GOptions_useLCorrection ) then
				
				Re = this.clusters( this.idSorted(n) ).radius()  ! El valor exacto es radius-rcov(1)-rcov(n), el valor que está es para una diatómica
				mu = 0.5*this.clusters( this.idSorted(n) ).mass()
				
				!---------------------------------------------------
				! Muestreo aleatorio de los l ( l1^2 + l2^2 <= E')
				!---------------------------------------------------
! 				allocate( sample(2,n-1) )
! 				
! 				call rs.init( nDim=2 )
! 				call rs.setRange( 1, [0.0_8,sqrt(this.reactorEnergy_-this.vibrationalEnergy_-this.intermolEnergy_)] )     ! r
! 				call rs.setRange( 2, [0.0_8,2.0_8*MATH_PI] )  ! phi
! 				
! 				trials = 0
! 				do while( .true. )
! 					call rs.uniform( sample )
! 					
! 					ssum = 0.0_8
! 					do i=1,n-1
! 						ssum = ssum + sample(1,i)
! 					end do
! 					
! 					if( ssum <= sqrt(this.reactorEnergy_-this.vibrationalEnergy_-this.intermolEnergy_) ) then
! 						exit
! 					end if
! 					
! 					if( trials > 100 ) then
! 						write(*,*) "blocked1 ...", ssum, sqrt(this.reactorEnergy_-this.vibrationalEnergy_-this.intermolEnergy_)
! 					end if
! 					
! 					trials = trials + 1
! 				end do
				
				!---------------------------------------------------
				
				j1 = 0.0_8
				j2 = 0.0_8
				do i=1,n-1
					m = 1.0_8/( 1.0_8/this.clusters( this.idSorted(i) ).mass() + 1.0_8/this.clusters( this.idSorted(n) ).mass() )
					r = norm2(this.clusters( this.idSorted(i) ).center()-this.clusters( this.idSorted(n) ).center())
						
					if( this.clusters( this.idSorted(n) ).fr() /= 0 ) then
						I_L = 1.0_8/( 1.0_8/( mu*Re**2  ) + 1.0_8/( m*r**2 ) )
					else
						I_L = m*r**2
					end if
					
! 					write(*,*) " m = ", m/amu, m
! 					write(*,*) "mu = ", mu/amu, mu
! 					write(*,*) " r = ", r/angs, r
! 					write(*,*) "Re = ", Re/angs, Re
					
					fr = 0
					select case( this.clusters( this.idSorted(n) ).fr() )
						case(0)
							select case( n )
								case(2)
! 									l1 = 0.0_8
! 									l2 = 0.0_8
									fr = fr + 0
								case(3)
! 									l1 = 0.0_8
! 									l2 = sample(1,i)
									fr = fr + 1
								case default
! 									l1 = sample(1,i)*sin(sample(2,i))
! 									l2 = sample(1,i)*cos(sample(2,i))
									fr = fr + 2
							end select
						case(2)
							select case( n )
								case(2)
! 									l1 = 0.0_8
! 									l2 = sample(1,i)
									fr = fr + 1
								case default
! 									l1 = sample(1,i)*sin(sample(2,i))
! 									l2 = sample(1,i)*cos(sample(2,i))
									fr = fr + 2
							end select
						case(3)
! 							l1 = sample(1,i)*sin(sample(2,i))
! 							l2 = sample(1,i)*cos(sample(2,i))
							fr = fr + 2
					end select
					
					if( this.clusters( this.idSorted(n) ).fr()==2 .and. n==3 ) fr = fr-1
					
! 					write(*,*) "fr(i) = ",	 fr, i, I_L, log(sqrt(2.0_8*I_L)), fr*log(sqrt(2.0_8*I_L))
					
! 					write(*,*) "fr = ", fr
					this.LnIm_ = this.LnIm_ + fr*log(sqrt(2.0_8*I_L))
! 					this.LnSqrt2Im_ = this.LnSqrt2Im_ + log(sqrt(2.0_8*I_L))
					
					
! 					this.rotationalEnergy = this.rotationalEnergy + ( l1**2 + l2**2 )
! 					this.rotationalEnergy = 0.0_8
! 					j1 = j1 - l1
! 					j2 = j2 - l2
				end do
					
! 				deallocate( sample )
				
! 				this.LnSqrt2Im_ = fr*this.LnSqrt2Im_
! 				this.rotationalEnergy = this.rotationalEnergy + ( j1**2 + j2**2 )
				
				if( GOptions_printLevel >= 3 ) then
! 					write(*,*) "fr = ", fr
					call GOptions_valueReport( "1) LnSqrt2Im", this.LnIm_, indent=2 )
! 					call GOptions_valueReport( "Erot", this.rotationalEnergy/eV, "eV", indent=2 )
				end if
		end if
		
! 		! Contribución orbital para el caso de 1 atomo y una molecula
! 		if( n==2 .and. GOptions_useLCorrection ) then
! 			if( this.clusters( this.idSorted(1) ).fr() == 0 .and. this.clusters( this.idSorted(2) ).fr() /= 0 ) then
! 				! Simplemente he utilizado como momento de inercia el valor asociado al momento orbital del átomo
! 				! alrededor de la molécula
! 				varIn = this.clusters( this.idSorted(1) ).mass()*this.clusters( this.idSorted(2) ).mass() &
! 							*norm2(this.clusters(1).center()-this.clusters(2).center())**2/this.mass()
! 							
! ! 				this.LnIm_ = this.LnIm_ + this.clusters( this.idSorted(2) ).fr()*log(sqrt(2.0_8*varIn))
! ! 				this.LnIm_ = this.LnIm_ + 2.0_8*log(sqrt(2.0_8*varIn))
! 				this.LnIm_ = this.LnIm_ + GOptions_frLCorrection*log(sqrt(2.0_8*varIn))
! 				
! ! 				call random_number(varJn)
! ! 				varJn = -sqrt(2.0_8*this.kineticEnergy()*varIn) + 2.0_8*varJn*sqrt(2.0_8*this.kineticEnergy()*varIn)
! ! 				write(*,*) varJn
! ! 				this.rotationalEnergy = varJn**2/(2.0_8*varIn)
! 			end if
! 		end if

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se obtiene la matriz inversa del tensor
		! de inercia de n. Es diagonal en body fix
		call Inn.init( 3, 3, 0.0_8 )
		if( this.clusters( this.idSorted(n) ).fr() /= 0 ) then
			do j=3,4-this.clusters( this.idSorted(n) ).fr(),-1
				Inn.data( j, j ) = this.clusters( this.idSorted(n) ).diagInertiaTensor.data( j, j )
			end do
		end if
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Contribución de los N-1 fragmentos a la energía rotacional
		do i=1,n-1
			if( this.clusters( this.idSorted(i) ).fr() /= 0 ) then
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Se obtiene la matriz inversa del tensor
				! de inercia de i
				Iii = this.clusters( this.idSorted(i) ).diagInertiaTensor
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Matriz de rotación que permite transformar a los ejes
				! de N en los ejes de i, pasando por el sistem fix
				! Rot = R_i^T*R_N
				Rot = SpecialMatrix_rotationTransform( &
					this.clusters( this.idSorted(i) ).inertiaAxes(), &
					this.clusters( this.idSorted(n) ).inertiaAxes() )
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Se obtiene la matriz inversa del tensor
				! de inercia efectivo.
				It =  Iii + Rot.transpose()*Inn*Rot
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				call It.eigen( eVals=Bi, eVecs=Ui )
				invBi = Bi.inverse()
				
				fr_sf = 3
				if( abs( invBi.get(1,1) ) > 10.0_8 ) then
					fr_sf = 2
				end if
! 				fr_sf = this.clusters( this.idSorted(i) ).fr()
				
				write(*,*) ">>> ", trim(this.label()), "  >>  ", trim(this.clusters( this.idSorted(i) ).label())
				write(*,*) "Iii = "
				call Iii.show( formatted=.true. )
				write(*,*) "Bi = "
				call Bi.show( formatted=.true. )
				write(*,*) "fr_sf,fr = ", fr_sf, this.clusters( this.idSorted(i) ).fr()
				
				this.clusters( this.idSorted(i) ).J_(j) = 0.0_8
				do j=3,4-fr_sf,-1
					call random_number( randNumber ) ! [0,1]
					this.clusters( this.idSorted(i) ).J_(j) = sqrt( 2.0_8*maxErot/invBi.get(j,j) )*randNumber
				end do
				
				! Momento angular para el fragmento i-esimo
				call Ji.columnVector( 3, values=this.clusters( this.idSorted(i) ).J_ )
				
				! Contribución a Jn
				Jn = Jn - Rot*Ji
				
				do j=3,4-fr_sf,-1
					Ein = 0.5_8*Ji.get(j,1)**2*invBi.get(j,j)
					
					this.LnIm_ = this.LnIm_ - log(invBi.get(j,j))
				end do
				
				this.rotationalEnergy = this.rotationalEnergy + Ein
			end if
		end do
		
		if( ( Math_isNaN(this.LnIm_) .or. Math_isInf(this.LnIm_) ) .and. GOptions_printLevel >= 2 ) then
			call GOptions_valueReport( "ErotJ", this.rotationalEnergy/eV, "eV", indent=2 )
			
			write(*,*) "-----------------------------"
			do i=1,n-1
				write(*,*) j, this.clusters( this.idSorted(i) ).J_
			end do
			write(*,*) "-----------------------------"
			stop
		end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Contribución de L
		if( this.fl() > 0 .and. this.clusters( this.idSorted(n) ).fr() /= 0 ) then
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Se obtiene la matriz inversa del tensor
			! de inercia de Il
			call Il.init( 3, 3, 0.0_8 )
			do j=3,4-this.fl(),-1
				Il.data( j, j ) = this.diagInertiaTensor( j )
			end do
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Matriz de rotación que permite transformar a los ejes
			! de N en los ejes de i, pasando por el sistem fix
			! Rot = R_i^T*R_N
			Rot = SpecialMatrix_rotationTransform( &
				this.inertiaAxes, &
				this.clusters( this.idSorted(n) ).inertiaAxes() )
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Se obtiene la matriz inversa del tensor
			! de inercia efectivo.
			It =  Il + Rot.transpose()*Inn*Rot
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			call It.eigen( eVals=Bi, eVecs=Ui )
			invBi = Bi.inverse()
			
			fr_sf = 3
			if( abs( invBi.get(1,1) ) > 10.0_8 ) then
				fr_sf = 2
			end if
			
			this.L_ = 0.0_8
			do j=3,4-fr_sf,-1
				call random_number( randNumber ) ! [0,1]
				this.L_(j) = sqrt( 2.0_8*maxErot/Bi.get(j,j) )*randNumber
			end do
			
			! Momento angular para el componente L
			call L.columnVector( 3, values=this.L_ )
			
			! Contribución a Jn
			Jn = Jn - Rot*L
			
			do j=3,4-fr_sf,-1
				Ein = 0.5_8*L.get(j,1)**2*invBi.get(j,j)
				
				this.LnIm_ = this.LnIm_ - log(invBi.get(j,j))
! 				write(*,*) "B", j, " = ", - log(invBi.get(j,j))
			end do
			
			this.rotationalEnergy = this.rotationalEnergy + Ein

! 			if( GOptions_useLWeightContrib ) then
	            do j=3,4-this.fl(),-1
					this.LnIm_ = this.LnIm_ - log(this.diagInertiaTensor(j))
! 					write(*,*) "I", j, " = ", - log(this.diagInertiaTensor(j))
	            end do
! 			end if
		end if
		
		! Puede ser que los componentes tengan problemas de asignación
		! cuando n sea lineal o un átomo, pero la norma de su momento angular
		! es correcta. Bueno, si es un átomo NO
		! @todo Verificar que alcance tiene cuando n es un átomo y se le está asignando un J
		this.clusters(this.idSorted(n)).J_ = Jn.data(:,1)
		
! 		if( GOptions_printLevel >= 3 .or. Math_isNaN(this.LnIm_) .or. Math_isInf(this.LnIm_) ) then
		if( ( Math_isNaN(this.LnIm_) .or. Math_isInf(this.LnIm_) .or. Math_isNaN(this.rotationalEnergy) ) .and. GOptions_printLevel >= 2 ) then
			call GOptions_valueReport( "ErotJL", this.rotationalEnergy/eV, "eV", indent=2 )
			write(*,*) "Im = ", this.diagInertiaTensor
			stop
		end if
	end subroutine updateRotationalEnergyJn
	
	!>
	!! @brief Actualiza la energía rotacional
	!!
	subroutine updateRotationalEnergyL( this )
		class(FragmentsList) :: this
		
		real(8) :: maxErot
		integer :: i, j, n
		type(Matrix) :: Inn, invInn    !< Tensor de inercia de N proyectado sobre los ejes de i y su inversa
		type(Matrix) :: Iii, invIii    !< Tensor de inercia de i proyectado sobre los ejes de i (diagonal) y su inversa
		type(Matrix) :: invIl
		type(Matrix) :: It, invIt      !< Tensor de inercia de i proyectado sobre los ejes de i (diagonal) y su inversa
		type(Matrix) :: Ji, L, Jn         !< Momentos angulares bf de i y n respectivamente
		type(Matrix) :: Rot            !< Matriz de rotación
		
		type(Matrix) :: Xi
! 		type(Matrix) :: Ai, invAi
		type(Matrix) :: Ui
		type(Matrix) :: Bi, invBi, tmp
		real(8) :: scaleFactor
		
		real(8) :: varIn, varJn !< Valores para el apaño de 1 átomo + 1 molecula en el momento orbital
		
		type(Matrix) :: Ein    !< Matriz auxiliar de 1x1, que representa la contribución a la energía rotacional por la interacción in
		real(8) :: randNumber
		
		if( this.forceInitializing ) then
			call this.initialGuessFragmentsList()
			return
		end if
		
		n = this.nMolecules()
		
		call Jn.columnVector( 3, val=0.0_8 )
		this.L_ = 0.0_8
		this.rotationalEnergy = 0.0_8
		this.LnIm_ = 0.0_8
		
		do i=1,n
			this.clusters(i).J_ = 0.0_8
		end do
		
		maxErot = this.reactorEnergy()-this.intermolEnergy_
		
		if( maxErot < 0.0_8 .or. n < 2 ) return
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se obtiene la matriz inversa del tensor
		! de inercia de n
		! Inversa de Inn que es diagonal en body fix
		call invInn.init( 3, 3, 0.0_8 )
		do j=3,4-this.fl(),-1
			invInn.data( j, j ) = 1.0_8/this.diagInertiaTensor( j )
		end do
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Contribución de los N-1 fragmentos a la energía rotational
		do i=1,n
			if( this.clusters( this.idSorted(i) ).fr() /= 0 ) then
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Se obtiene la matriz inversa del tensor
				! de inercia de i
				Iii = this.clusters( this.idSorted(i) ).diagInertiaTensor
				
				! Inversa de Iii que es diagonal en body fix
				call invIii.init( 3, 3, 0.0_8 )
				do j=3,4-this.clusters( this.idSorted(i) ).fr(),-1
					invIii.data( j, j ) = 1.0_8/Iii.data( j, j )
				end do
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Matriz de rotación que permite transformar a los ejes
				! de N en los ejes de i, pasando por el sistem fix
				! Rot = R_i^T*R_N
				Rot = SpecialMatrix_rotationTransform( &
					this.inertiaAxes, &
					this.clusters( this.idSorted(i) ).inertiaAxes() )
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Se obtiene la matriz inversa del tensor
				! de inercia efectivo.
				invIt =  Rot.transpose()*invInn*Rot + invIii
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				call Ui.init( 3, 3, 0.0_8 )
				call Bi.init( 3, 3, 0.0_8 )
				call invIt.eigen( eVals=Bi, eVecs=Ui )
				
				this.clusters( this.idSorted(i) ).J_(j) = 0.0_8
				do j=3,4-this.clusters( this.idSorted(i) ).fr(),-1
					call random_number( randNumber ) ! [0,1]
					this.clusters( this.idSorted(i) ).J_(j) = sqrt( 2.0_8*maxErot/Bi.get(j,j) )*randNumber
				end do
				
				! Momento angular para el fragmento i-esimo
				call Ji.columnVector( 3, values=this.clusters( this.idSorted(i) ).J_ )
				
				! Contribución a Jn
				Jn = Jn - Rot*Ji
				
				! Contribución a la energía rotacional
				Ein = Ji.transpose()*invIt*Ji*0.5_8
				
				this.rotationalEnergy = this.rotationalEnergy + Ein.get(1,1)
				
				do j=3,4-this.clusters( this.idSorted(i) ).fr(),-1
					this.LnIm_ = this.LnIm_ - log(Bi.get(j,j))
				end do
			end if
		end do
		
		if( ( Math_isNaN(this.LnIm_) .or. Math_isInf(this.LnIm_) ) .and. GOptions_printLevel >= 2 ) then
			call GOptions_valueReport( "ErotJ", this.rotationalEnergy/eV, "eV", indent=2 )
			
			write(*,*) "-----------------------------"
			do i=1,n-1
				write(*,*) j, this.clusters( this.idSorted(i) ).J_
			end do
			write(*,*) "-----------------------------"
			stop
		end if
		
		this.L_ = Jn.data(:,1)
		
		if( ( Math_isNaN(this.LnIm_) .or. Math_isInf(this.LnIm_) ) .and. GOptions_printLevel >= 2 ) then
			call GOptions_valueReport( "ErotJL", this.rotationalEnergy/eV, "eV", indent=2 )
			write(*,*) "Im = ", this.diagInertiaTensor
			stop
		end if
	end subroutine updateRotationalEnergyL
	
	!>
	!! @brief Translational weight
	!!
	subroutine updateLnWt( this )
		class(FragmentsList) :: this
		
		integer :: i, n, s
		real(8) :: logMu, Et
		
		integer :: fr
		
		n = this.nMolecules()
		
		Et = this.kineticEnergy()
		
		if( Et < 0.0_8 ) then
			this.LnWt = 0.0_8
! 			write(*,"(A,A,F10.5)") "E < 0", "   "//trim(this.label()), this.LnWt
		else if( this.nMolecules() == 1 ) then
! 			this.LnWt = this.logVfree_ + this.logVtheta_
			this.LnWt = 0.0_8
! 			write(*,"(A,A,F10.5)") "n = 0", "   "//trim(this.label()), this.LnWt
		else
			logMu = 0.0_8
			fr = 0
			do i=1,n
				logMu = logMu + log(this.clusters(i).mass())
				
				if( i /= n ) then
					fr = fr + this.clusters( this.idSorted(i) ).fr()
				end if
			end do
			logMu = logMu - log( this.mass() ) !+ (n-1.0_8)*log(2.0_8)
			
			s = this.ft() + this.fl() + fr
			
			this.LnWt = \
				0.5_8*s*log(2.0_8*Math_PI) - log( Gamma(0.5_8*s) ) &
				+ 1.5*logMu &
				+ this.logVfree_ + this.logVtheta_ &
				+ (0.5_8*s-1.0_8)*log(Et)
				
			this.LnWr =  0.5_8*this.LnIm_
			
! 			write(*,"(A,4I5,A,F10.5)") "s = ", s, this.ft(), this.fl(), fr, "   "//trim(this.label()), this.LnWt
		end if
		
		if( ( Math_isNaN(this.LnWt) .or. Math_isInf(this.LnWt) ) .and. GOptions_printLevel >= 2 ) then
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
			call GOptions_valueReport( "LnWt", this.LnWt, indent=2 )
		end if
	end subroutine updateLnWt
		
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
! 		GOptions_printRigidMoleculeData = .false.
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

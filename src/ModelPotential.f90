!>
!! @brief
!!
module ModelPotential_
	use GOptions_
	use Math_
	use String_
	use UnitsConverter_
	
	use GOptionsM3C_
	use AzizSlamanPotential_
	
	implicit none
	private
	
	public :: &
		ModelPotential_test
		
	!>
	!! Available intermolecular energy functions
	!!
	enum, BIND(c)
		enumerator :: MODEL_NONE = 0
		enumerator :: MODEL_HARDSPHERE
		enumerator :: MODEL_MORSE
		enumerator :: MODEL_LENNARDJONESRE
		enumerator :: MODEL_LENNARDJONESR0
		enumerator :: MODEL_AZIZSLAMAN
		enumerator :: MODEL_HARMONIC
		enumerator :: MODEL_COULOMB
		enumerator :: MODEL_DHARMCOUL
		enumerator :: MODEL_DMORSECOUL
		enumerator :: MODEL_NEXP
		enumerator :: MODEL_FLATR6
		enumerator :: MODEL_FLATEXP
	end enum
	integer, parameter :: MODEL_SIZE = MODEL_FLATEXP
	
	character(*), public, parameter :: MODEL_NAME(0:MODEL_SIZE) = [ &
		"NONE", &
		"HARDSPHERE", &
		"MORSE", &
		"LENNARDJONESRE", &
		"LENNARDJONESR0", &
		"AZIZSLAMAN", &
		"HARMONIC", &
		"COULOMB", &
		"DHARMCOUL", &
		"DMORSECOUL", &
		"NEXP", &
		"FLATR6", &
		"FLATEXP" &
	]
	
	integer, parameter :: MODEL_NUSERPARAMS(0:MODEL_SIZE) = [ &
		0, & ! NONE
		0, & ! HARDSPHERE
		2, & ! MORSE
		4, & ! LENNARDJONESRE
		4, & ! LENNARDJONESR0
		1, & ! AZIZSLAMAN
		3, & ! HARMONIC
		1, & ! COULOMB
		4, & ! DHARMCOUL
		4, & ! DMORSECOUL
		3, & ! NEXP
		3, & ! FLATR6
		4  & ! FLATEXP
	]
	
	type, public :: ModelPotential
		integer, private :: id = 0
		real(8) :: mu = 0.0_8
		real(8), allocatable :: potParams(:)
		
		contains
			generic :: init => initDefault, fromFString
			generic :: assignment(=) => copyModelPotential
			
			procedure :: initDefault
			procedure :: fromFString
			procedure :: copyModelPotential
			final :: destroyModelPotential
			
			procedure :: getId
			procedure :: evaluate
	end type ModelPotential
	
	contains
	
	!>
	!! @brief Constructor
	!!
	subroutine initDefault( this, id, potParams )
		class(ModelPotential) :: this 
		integer, optional, intent(in) :: id
		real(8), optional, intent(in) :: potParams(:)
		
		integer :: idEff
		
		if( present(id) ) then
			this.id = id
		else
			this.id = MODEL_NONE
		end if
		
		if( allocated(this.potParams) ) deallocate(this.potParams)
		
		if( present(potParams) ) then
			allocate( this.potParams(size(potParams)) )
			this.potParams = potParams
		end if
	end subroutine initDefault
	
	!>
	!! @brief Constructor
	!!
	subroutine fromFString( this, fstr, mu )
		class(ModelPotential) :: this 
		character(*), intent(in) :: fstr
		real(8), intent(in) :: mu
		
		integer :: i
		character(100), allocatable :: tokens(:)
		real(8), allocatable :: potUserParams(:)
		logical :: located
		
		this.mu = mu
		
		call FString_split( fstr, tokens, "()" )
		
		if( size(tokens) > 1 ) then
			call FString_toRealArray( tokens(2), potUserParams )
		end if
		
		located = .false.
		call this.initDefault()
		do i=1,size(MODEL_NAME)
			if( trim(adjustl(tokens(1))) == MODEL_NAME(i) ) then
				if( size(potUserParams) < MODEL_NUSERPARAMS(i) ) then
					call GOptions_error( &
						"Incomplete number of parameters for model --"//trim(fstr)//"--", &
						"ModelPotential.fromFString()", &
						"The model --"//trim(fstr)//"-- need "//trim(FString_fromInteger( MODEL_NUSERPARAMS(i) ))//" parameters" &
						)
				end if
				
				! Where is not necessary to use user parameters, for example in COULOMB model
				if( size(potUserParams) /= 0 ) then
					call this.initDefault( i, potUserParams )
				else
					call this.initDefault( i )
				end if
				
				located = .true.
				exit
			end if
		end do
		
		if( allocated(tokens) ) deallocate( tokens )
		if( allocated(potUserParams) ) deallocate( potUserParams )
		
		if( .not. located ) then
			write(*,"(A)") ""
			write(*,"(A)") "### ERROR ### ModelPotential.fromFString()"
			write(*,"(A)") "              Unknown model potential --"//trim(fstr)//"--"
			
			write(*,"(A)") "   Available models:"
			do i=1,size(MODEL_NAME)
				write(*,"(A,A)") "      ", MODEL_NAME(i)
			end do
			
			stop
		end if
	end subroutine fromFString
	
	!*
	! @brief Copy constructor
	!*
	subroutine copyModelPotential( this, other )
		class(ModelPotential), intent(out) :: this
		class(ModelPotential), intent(in) :: other

		this.id = other.id
		
		if( allocated(this.potParams) ) deallocate(this.potParams)
		allocate(this.potParams(size(other.potParams)))
		this.potParams = other.potParams
	end subroutine copyModelPotential
	
	!*
	! @brief Destructor
	!*
	subroutine destroyModelPotential( this )
		type(ModelPotential) :: this
		
		if( allocated(this.potParams) ) deallocate(this.potParams)
	end subroutine destroyModelPotential
	
	!>
	!! @brief Auxiliar function
	!!
	pure function damp( x, s, x0, w, h ) result( output )
		real(8), intent(in) :: x
		integer, intent(in) :: s   ! signo +1 o -1 si s es positivo o no respectivamente
		real(8), intent(in) :: x0  ! valor donde cambia
		real(8), intent(in) :: w   ! ancho
		real(8), intent(in) :: h   ! diferencia max-min
		real(8) :: output
		
		output = h*(0.5_8*s*(1.0_8+tanh(6.0_8*(x-x0)/w))-(s-1.0_8)*0.5_8)
	end function damp
	
	!>
	!! @brief Internal function. Code = HARDSPHERE
	!!
	pure function hardSphere( r, userParams, params ) result( output )
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: userParams(:)
		real(8), optional, intent(in) :: params(:)
		real(8) :: output
		
		real(8) :: E1, Re
		
		Re = params(1) ! Ra+Rb = Re
		
		if( r > Re-GOptionsM3C_overlappingRadius ) then
			output = 0.0_8
		else
			output = (10.0d3)*eV
		end if
	end function hardSphere
	
	!>
	!! @brief Internal function. Code = MORSE
	!!
	pure function morse( r, userParams, params ) result( output )
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: userParams(:)
		real(8), optional, intent(in) :: params(:)
		real(8) :: output
		
		real(8) :: De, alpha, Re
		
		De = userParams(1)
		alpha = userParams(2)
		
		Re = params(1)
		
		! @todo Posiblemente Re debe ser tomado de los params (suma de radios covalentes)
		
		output = De*( exp(-2.0*alpha*(r-Re)) - 2.0*exp(-alpha*(r-Re)) )
	end function morse

	!>
	!! @brief Internal function. Code = LENNARDJONESRE
	!!
	pure function lennardJonesRe( r, userParams, params ) result( output )
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: userParams(:)
		real(8), optional, intent(in) :: params(:)
		real(8) :: output
		
		real(8) :: E0, E1
		real(8) :: De, Re
		integer :: n, m
		real(8) :: Cn, Cm
		
		De = userParams(1)
		Re = userParams(2)
		n  = int(userParams(3))
		m  = int(userParams(4))
		
		! @todo No se puede utilizar GOptions_error, ni write, ni stop en una función pura
		! así que hay que aclarar que siempre se utilizará Cnm=4 hasta que se implementen
		! otros valores
! 		if( n == 12 .and. m == 6 ) then
			Cn = 1.0_8
			Cm = 2.0_8
! 		end if
		
		output = De*( Cn*(Re/r)**n - Cm*(Re/r)**m )
	end function lennardJonesRe
	
	!>
	!! @brief Internal function. Code = LENNARDJONESR0
	!! @todo Esta mal implementada
	!!
	pure function lennardJonesR0( r, userParams, params ) result( output )
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: userParams(:)
		real(8), optional, intent(in) :: params(:)
		real(8) :: output
		
		real(8) :: De, R0, E0
		integer :: n, m
		
		real(8) :: Cnm
		
		De = userParams(1)
		R0 = userParams(2)
		n  = int(userParams(3))
		m  = int(userParams(4))
		
		! @todo No se puede utilizar GOptions_error, ni write, ni stop en una función pura
		! así que hay que aclarar que siempre se utilizará Cnm=4 hasta que se implementen
		! otros valores
! 		if( n == 12 .and. m == 6 ) then
			Cnm = 4.0_8
! 		end if
		
		output = Cnm*De*( (R0/r)**n - (R0/r)**m )
	end function lennardJonesR0
	
	
	!>
	!! @brief Internal function. Code = AZIZSLAMAN
	!! @todo Esta mal implementada
	!!
! 	pure function azizSlaman( r, userParams, params ) result( output )
	function azizSlaman( r, userParams, params ) result( output )
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: userParams(:)
		real(8), optional, intent(in) :: params(:)
		real(8) :: output
		
		integer :: kindAtom
		type(AzizSlamanPotential) :: aziz

		kindAtom  = int(userParams(1))
		
		select case( kindAtom )
			case(1)
				call aziz.init( "He" )
			case(2)
				call aziz.init( "Ar" )
			case default
				stop
		end select
		
		output = aziz.V( r/angs )*cm1
	end function azizSlaman
	
	!>
	!! @brief Internal function. Code = HARMONIC
	!!
	pure function harmonic( r, userParams, params ) result( output )
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: userParams(:)
		real(8), optional, intent(in) :: params(:)
		real(8) :: output
		
		real(8) :: De, alpha, Re
		
		De = userParams(1)
		alpha = userParams(2)
		Re = userParams(3)
		
		output = De*alpha**2*(r-Re)**2-De
		
! 		if( output > 0.0  ) output = 0.0_8
	end function harmonic
	
	!>
	!! @brief Internal function. Code = COULOMB
	!!
	pure function coulomb( r, userParams, params ) result( output )
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: userParams(:)
		real(8), optional, intent(in) :: params(:)
		real(8) :: output
		
		real(8) :: E0, E1
		real(8) :: Re
		real(8) :: q1, q2
		
		real(8) :: Z
		
		E0 = params(1) ! Eab = E0
		E1 = params(2) ! Ea+Eb = E1
		Re = params(3) ! Ra+Rb = Re
		q1 = params(4) ! qa = q1
		q2 = params(5) ! qb = q2
		
		! @todo No estoy seguro de las unidades de Z
		Z = userParams(1)
		
		output = Z*q1*q2/r+E1
	end function coulomb	
	
	!>
	!! @brief Internal function. Code = DHARMCOUL
	!!
	function dHarmCoul( r, userParams, params ) result( output )
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: userParams(:)
		real(8), optional, intent(in) :: params(:)
		real(8) :: output
		
		real(8) :: De, alpha, Re, R0, w, Z
		
		real(8) :: E0, E1
		real(8) :: q1, q2
		
		E0 = params(1) ! Eab = E0
		E1 = params(2) ! Ea+Eb = E1
		Re = params(3) ! Ra+Rb = Re
		q1 = params(4) ! qa = q1
		q2 = params(5) ! qb = q2
		
		alpha = userParams(1)
		R0 = userParams(2)
		w = userParams(3)
		Z = userParams(4)
		
		output = damp( r, -1, R0, w, 1.0_8 )*harmonic( r, userParams=[alpha], params=params ) &
				 + damp( r, 1, R0, w, 1.0_8 )*coulomb( r, userParams=[Z], params=params )
	end function dHarmCoul
	
	!>
	!! @brief Internal function. Code = DMORSECOUL
	!!
	pure function dMorseCoul( r, userParams, params ) result( output )
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: userParams(:)
		real(8), optional, intent(in) :: params(:)
		real(8) :: output
		
		real(8) :: De, alpha, Re, R0, w, Z
		
		real(8) :: E0, E1
		real(8) :: q1, q2
		
		E0 = params(1) ! Eab = E0
		E1 = params(2) ! Ea+Eb = E1
		Re = params(3) ! Ra+Rb = Re
		q1 = params(4) ! qa = q1
		q2 = params(5) ! qb = q2
		
		! La conversión de unidades de hace dentro de cada función que se llama
		alpha = userParams(1)
		De = userParams(2)
		R0 = userParams(3)
		w = userParams(4)
		Z = userParams(5)
		
		output = damp( r, -1, Re+R0, w, 1.0_8 )*morse( r, userParams=[alpha],params=[E0,E0+De,Re,q1,q2] ) &
				 + damp( r, 1, Re+R0, w, 1.0_8 )*coulomb( r, userParams=[Z], params=params )
	end function dMorseCoul
	
	!>
	!! @brief Internal function. Code = NEXP
	!!
	function nexp( r, userParams, params ) result( output )
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: userParams(:)
		real(8), optional, intent(in) :: params(:)
		real(8) :: output
		
		real(8) :: Re
		real(8) :: dE, b, c
		
		real(8) :: A, r0, r1
		
		Re = params(1) ! Ra+Rb = Re
		
		dE = userParams(1)*eV
		c = userParams(2)*angs
		b = userParams(3)/angs
		
		r0 = Re - log(2.0)/(1.2_8/angs)
		r1 = Re + c
! 		write(*,*) Re/angs, r0/angs, r1/angs
! 		stop

		A = dE*exp( b*(r1-Re) )
		
		output = A*exp(-b*(r-Re))
		if( r < r1 ) output = dE
		if( r < r0 ) output = 10.0d3*eV
		
! 		if( r > Re ) then
! ! 			output = A*exp(-b*(r-Re-GOptionsM3C_overlappingRadius)) + E1
! 			output = A*exp(-b*(r-Re-GOptionsM3C_overlappingRadius))
! 		else
! 			output = (10.0d3)*eV
! 		end if
	end function nexp
	
	!>
	!! @brief Internal function. Code = FLATR6
	!!
	function flatR6( r, userParams, params ) result( output )
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: userParams(:)
		real(8), optional, intent(in) :: params(:)
		real(8) :: output
		
		real(8) :: r0, r1, c6
		
		r0 = userParams(1)
		r1 = userParams(2)
		c6 = userParams(3)
		
		if( r < r0 ) then
			output = (10.0d3)*eV
		else if( r < r1 ) then
			output = c6/r1**6
		else
			output = c6/r**6
		end if
	end function flatR6
	
	!>
	!! @brief Internal function. Code = FLATEXP
	!!
	function flatExp( r, userParams, params ) result( output )
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: userParams(:)
		real(8), optional, intent(in) :: params(:)
		real(8) :: output
		
		real(8) :: r0, r1, a, b
		
		r0 = userParams(1)
		r1 = userParams(2)
		a = userParams(3)
		b = userParams(4)
		
		if( r < r0 ) then
			output = (10.0d3)*eV
		else if( r < r1 ) then
			output = a*exp(-b*r1)
		else
			output = a*exp(-b*(r-r0))
		end if
	end function flatExp
	
	!>
	!! @brief
	!!
	function getId( this ) result( id )
		class(ModelPotential), intent(in) :: this
		integer :: id
		
		id = this.id
	end function getId
		
	!>
	!! @brief
	!!
	function evaluate( this, r, params ) result( V )
		class(ModelPotential) :: this
		real(8), intent(in) :: r
		real(8), optional, intent(in) :: params(:)
		real(8) :: V
		
		select case( this.id )
			case( MODEL_HARDSPHERE )
				V = hardSphere( r, this.potParams, params )
			case( MODEL_MORSE )
				V = morse( r, this.potParams, params )
			case( MODEL_LENNARDJONESRE )
				V = lennardJonesRe( r, this.potParams, params )
			case( MODEL_LENNARDJONESR0 )
				V = lennardJonesR0( r, this.potParams, params )
			case( MODEL_AZIZSLAMAN )
				V = azizSlaman( r, this.potParams, params )
			case( MODEL_HARMONIC )
				V = harmonic( r, this.potParams, params )
			case( MODEL_COULOMB )
				V = coulomb( r, this.potParams, params )
			case( MODEL_DHARMCOUL )
				V = dHarmCoul( r, this.potParams, params )
			case( MODEL_DMORSECOUL )
				V = dMorseCoul( r, this.potParams, params )
			case( MODEL_NEXP )
				V = nexp( r, this.potParams, params )
			case( MODEL_FLATR6 )
				V = flatR6( r, this.potParams, params )
			case( MODEL_FLATEXP )
				V = flatExp( r, this.potParams, params )
			case default
				V = 0.0_8
		end select
		
	end function evaluate
		
	!>
	!! @brief Test method
	!!
	subroutine ModelPotential_test()
		
	end subroutine ModelPotential_test
	
end module ModelPotential_

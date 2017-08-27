!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                   !!
!! This file is part of M3C project                                                  !!
!! Copyright (c) 2012-2016 Departamento de Química                                   !!
!!                         Universidad Autónoma de Madrid                            !!
!!                         All rights reserved.                                      !!
!!                                                                                   !!
!!                         * Néstor F. Aguirre (2012-2016)                           !!
!!                           nestor.aguirre@uam.es                                   !!
!!                         * Sergio Díaz-Tendero (2012-2016)                         !!
!!                           sergio.diaztendero@uam.es                               !!
!!                         * M. Paul-Antoine Hervieux (2012-2015)                    !!
!!                           Paul-Antoine.Hervieux@ipcms.unistra.fr                  !!
!!                         * Manuel Alcamí (2012-2016)                               !!
!!                           manuel.alcami@uam.es                                    !!
!!                         * Fernando Martín (2012-2016)                             !!
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

module AzizSlamanPotential_
	use Math_
	use UnitsConverter_
	implicit none
	private
	
	public :: &
		AzizSlamanPotential_test
		
	real(8), private, parameter :: tocm1 = 0.6950356
	
	type, public :: AzizSlamanPotential
		real(8) :: epsilon
		real(8) :: Rm
		real(8) :: A
		real(8) :: alpha
		real(8) :: beta
		real(8) :: c6
		real(8) :: c8
		real(8) :: c10
		real(8) :: D
		real(8) :: Aa
		real(8) :: zeta1
		real(8) :: zeta2
		
		contains
			procedure :: init
			procedure :: V
			procedure :: dV
			procedure :: NdV
			
			procedure, private :: Va
			procedure, private :: dVa
			procedure, private :: Vb
			procedure, private :: dVb
			procedure, private :: F
			procedure, private :: dF
	end type AzizSlamanPotential
	
	contains
	
	subroutine init( this, atom )
		class(AzizSlamanPotential) :: this
		character(*), intent(in) :: atom
		
		select case( trim(atom) )
			case ("He")
				! HFD-B base potential
				! LM2M2
				! R.A. Aziz and M.J. Slaman J.Chem.Phys 94 (1991) 8047
				! Table V.
				this.A = 189635.353_8
				this.alpha = 10.70203539_8
				this.c6 = 1.34687065_8
				this.c8 = 0.41308398_8
				this.c10 = 0.17060159_8
				this.beta = 1.90740649_8
				this.D = 1.4088_8
				this.epsilon = 10.97_8
				this.Rm = 2.9695_8
				this.Aa = 0.0026_8
				this.zeta1 = 1.003535949_8
				this.zeta2 = 1.454790369_8
			case ("Ar")
				! HFD-B base potential
				! LM2M2
				! R.A. Aziz and M.J. Slaman J.Chem.Phys 94 (1991) 8047
				! Table V.
				this.A = 226210.716_8
				this.alpha = 10.77874743_8
				this.c6 = 1.10785136_8
				this.c8 = 0.56072459_8
				this.c10 = 0.34602794_8
				this.beta = 1.8122004_8
				this.D = 1.36_8
				this.epsilon = 143.224_8
				this.Rm = 3.7565_8
				this.Aa = 0.0_8
				this.zeta1 = -1.0_8  !<<< falsos
				this.zeta2 = -2.0_8  !<<< falsos
			case default
				write(*,*) "### ERROR ### AzizSlamanPotential.init(): atom not suppoted ("//trim(atom)//")"
				stop
		end select

	end subroutine
	
	function V( this, R ) result( output )
		class(AzizSlamanPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		output = this.epsilon*( this.Vb( R/this.Rm ) + this.Va( R/this.Rm ) )*tocm1
	end function V
	
	function dV( this, R ) result( output )
		class(AzizSlamanPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		output = (this.epsilon/this.Rm)*( this.dVb( R/this.Rm ) + this.dVa( R/this.Rm ) )*tocm1
	end function dV
	
	function NdV( this, R ) result( output )
		class(AzizSlamanPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		real(8) :: h = 0.00001_8
		
! http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
		output = ( this.V(R-2.0_8*h)-8.0_8*this.V(R-h)+8.0_8*this.V(R+h)-this.V(R+2.0_8*h) )/(12.0_8*h)
! 		output = ( -this.V(R-3.0_8*h)+9.0_8*this.V(R-2.0_8*h)-45.0_8*this.V(R-h) &
! 			   +45.0_8*this.V(R+h)-9.0_8*this.V(R+2.0_8*h)+this.V(R+3.0_8*h) )/(60.0_8*h)
	end function NdV
	
	function Va( this, zeta ) result( output )
		class(AzizSlamanPotential), intent(in) :: this
		real(8), intent(in) :: zeta
		real(8) :: output
		
		if( this.zeta1 <= zeta .and. zeta <= this.zeta2  ) then
			output = this.Aa*( sin( 2.0_8*Math_PI*(zeta-this.zeta1)/(this.zeta2-this.zeta1) - 0.5_8*Math_PI) + 1.0_8 )
		else
			output = 0.0_8
		end if
	end function Va
	
	function dVa( this, zeta ) result( output )
		class(AzizSlamanPotential), intent(in) :: this
		real(8), intent(in) :: zeta
		real(8) :: output
		
		if( this.zeta1 <= zeta .and. zeta <= this.zeta2  ) then
			output = (2.0_8*Math_PI*this.Aa/(this.zeta2-this.zeta1))*sin( 2.0_8*Math_PI*(zeta-this.zeta1)/(this.zeta2-this.zeta1) )
		else
			output = 0.0_8
		end if
	end function dVa
	
	function Vb( this, zeta ) result( output )
		class(AzizSlamanPotential), intent(in) :: this
		real(8), intent(in) :: zeta
		real(8) :: output
		
		output = this.A*exp( -this.alpha*zeta-this.beta*zeta**2.0_8 ) &
			   -( this.c6/zeta**6.0_8 + this.c8/zeta**8.0_8 + this.c10/zeta**10.0_8 )*this.F(zeta)
	end function Vb
	
	function dVb( this, zeta ) result( output )
		class(AzizSlamanPotential), intent(in) :: this
		real(8), intent(in) :: zeta
		real(8) :: output
		
		output = -this.A*(this.alpha+2.0_8*this.beta*zeta)*exp( -this.alpha*zeta-this.beta*zeta**2.0_8 ) &
			   +( 6.0_8*this.c6/zeta**7.0_8 + 8.0_8*this.c8/zeta**9.0_8 + 10.0_8*this.c10/zeta**11.0_8 )*this.F(zeta) &
			   -( this.c6/zeta**6.0_8 + this.c8/zeta**8.0_8 + this.c10/zeta**10.0_8 )*this.dF(zeta)
	end function dVb
	
	function F( this, zeta ) result( output )
		class(AzizSlamanPotential), intent(in) :: this
		real(8), intent(in) :: zeta
		real(8) :: output
		
		if( zeta <= this.D ) then
			output = exp( -(this.D/zeta-1.0_8)**2.0_8 )
		else
			output = 1.0_8
		end if
	end function F
	
	function dF( this, zeta ) result( output )
		class(AzizSlamanPotential), intent(in) :: this
		real(8), intent(in) :: zeta
		real(8) :: output
		
		if( zeta <= this.D ) then
			output = 2.0_8*this.D*(this.D/zeta-1.0_8)*exp( -(this.D/zeta-1.0_8)**2.0_8 )/zeta**2.0_8
		else
			output = 0.0_8
		end if
	end function dF
	
	subroutine AzizSlamanPotential_test()
		real(8) :: r
		type(AzizSlamanPotential) :: potential
		
		call potential.init( "Ar" )
				
		do r = 2.0,10.0,0.001
			write(*,"(3F15.7)") r, potential.V( r ), potential.dV( r )
		end do
		
		write(*,*) ""
		write(*,*) ""
		
	end subroutine AzizSlamanPotential_test
	
end module AzizSlamanPotential_

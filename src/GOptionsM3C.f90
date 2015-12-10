!>
!! @brief
!!
module GOptionsM3C_
	use UnitsConverter_
	use IOStream_
	use Timer_
	
	implicit none
	public
	
	real(8) :: GOptionsM3C_systemRadius = 10.0_8*angs
	real(8) :: GOptionsM3C_randomWalkStepRadius = 2.0_8*angs
	real(8) :: GOptionsM3C_overlappingRadius = 0.0_8*angs
	logical :: GOptionsM3C_useWeightedWalkStep = .false.
	logical :: GOptionsM3C_useRandomWalkers = .true.
	logical :: GOptionsM3C_useZPECorrection = .false.
	logical :: GOptionsM3C_useLCorrection = .false.
	logical :: GOptionsM3C_useLDOSContrib = .false.
	logical :: GOptionsM3C_useLWeightContrib = .false.
	real(8) :: GOptionsM3C_gammaLCorrection = 1.0_8
	logical :: GOptionsM3C_useLReference = .false.
	logical :: GOptionsM3C_useSpinConservationRules = .false.
end module GOptionsM3C_

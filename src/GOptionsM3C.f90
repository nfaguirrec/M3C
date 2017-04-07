!>
!! @brief
!!
module GOptionsM3C_
	use GOptions_
	use UnitsConverter_
	use IOStream_
	use Timer_
	use String_
	
	implicit none
	public
	
	real(8) :: GOptionsM3C_systemRadius = 10.0_8*angs
	real(8) :: GOptionsM3C_randomWalkStepRadius = 2.0_8*angs
	real(8) :: GOptionsM3C_overlappingRadius = 0.0_8*angs
	integer :: GOptionsM3C_radiusType = AtomicElementsDB_COVALENT_RADIUS
	logical :: GOptionsM3C_useWeightedWalkStep = .false.
	logical :: GOptionsM3C_useRandomWalkers = .false.
	logical :: GOptionsM3C_useZPECorrection = .false.
	logical :: GOptionsM3C_useSpinConservationRules = .false.
	type(String) :: GOptionsM3C_angularMomentumCouplingScheme ! = "JJ"
	real(8) :: GOptionsM3C_totalJ(3) = 0.0_8
	type(String) :: GOptionsM3C_structureSamplingMethod ! = "RANDOM"
	logical :: GOptionsM3C_checkAtomicOverlapping = .false.
	real(8) :: GOptionsM3C_atomicOverlappingRadius = 0.0_8*angs
end module GOptionsM3C_

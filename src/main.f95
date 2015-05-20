program QuantumDynamics

use qdprocedures
  
  !Basic configuration
  complex*8, allocatable :: wavefunctions(:,:,:)
  logical                :: reciprocalBase
  logical                :: endAtBoxEdge
  logical                :: storeData

  !Box properties
  !Set nDivZ to 1 to do a 2D simulation, or set nDivY to 1 to do a 1D simulation.
  !For indefinite simulations (please do not use with data storage, that will get
  !    out of hand pretty quickly in terms of storage space if accidentally left
  !    to run for too long) set indefSimulation to .TRUE.
  real*8                 :: subdivisionLength = 1d0
  integer*8              :: nDivX = 100, nDivY = 1, nDivZ = 1
  integer*8              :: nTimeSteps = 100
  real*8                 :: timeStepLength = 1d0
  logical                :: indefSimulation = .FALSE.

  real*8                 :: boxLength

  logical                :: potentialTimeDependent
  real*8, allocatable    :: potential(:,:,:,:)


  boxLength = subdivisionLength*nDivX

!  call potential_insert_block(start(dim),end(dim),depth,sl,potential)
!!  call potential_insert_ellipse(center(dim),radius(dim),depth,sl,potential)
!  call potential_insert_coulomb(center(dim),depth,sl,potential)
!  call potential_insert_slope(start(dim),end(dim),grad(dim+1,dim+1),sl,potential)
!!  call insert_wavefunction_ellipse(center(dim),radius(dim),depth,sl,wavefunctions)  
!  call insert_wavefunction_sine(start(dim),end(dim),momentum(dim),sl,wavefunctions)
!  call insert_wavefunction_gaussian(mean(dim),variance(dim,dim),sl,wavefunctions)
!  call insert_wavefunction_block(start(dim),end(dim),height,sl,wavefunctions)

! The nature of the Coulomb potential depends on the dimensionality: 1D: -k 2D = k*ln(r) <zero at the borders>
! while 3D = k/r. If one wants to use a 

  if (endAtBoxEdge) then
     allocate(wavefunctions(nDivX,nDivY,nDivZ)
  end if
    
  wavefunctions = 0

  allocate(potential(nDivX,nDivY,nDivZ,nTimeSteps))

  potential = 0
  
!  if(reciprocalBase.eqv..FALSE.) then
     call real_space_qd(wavefunctions,potential,potentialTimeDependent,endAtBoxEdge,nTimeSteps)
!  else
!     call reciprocal_space_qd(wavefunctions,potential)
!  end if


! See below for the real/reciprocal space functions

end program QuantumDynamics






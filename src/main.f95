program QuantumDynamics

use qdprocedures
  
  !Basic configuration
  complex*16, allocatable :: wavefunctions(:,:,:)
  logical                :: reciprocalBase
  logical                :: endAtBoxEdge
  logical                :: storeData

  !Box properties
  !Set nDivZ to 1 to do a 2D simulation, or set nDivY to 1 to do a 1D simulation.
  !For indefinite simulations (please do not use with data storage, that will get
  !    out of hand pretty quickly in terms of storage space if accidentally left
  !    to run for too long) set indefSimulation to .TRUE.
  !Note that E in J = 1, m = m_part = 1, and hbar = 1 are sufficient to obtain
  !suitable units for time (hbar/J) and length(

  integer*8              :: nDivX = 1000, nDivY = 1, nDivZ = 1
  integer*8              :: nTimeSteps = 100
  real*8                 :: timeStepLength = 1d0
  logical                :: indefSimulation = .FALSE.
  real*8                 :: deltat = 4d-3, deltax = 1d0
  real*8                 :: accuracy = 1d-60

  real*8                 :: boxLength

  logical                :: potentialTimeDependent = .FALSE.
  real*8, allocatable    :: potential(:,:,:,:)
  integer*8              :: i


  boxLength = subdivisionLength*nDivX


  open(1,file="wavefunction.csv")
  open(2,file="intensities.dat")
  
  allocate(wavefunctions(nDivX,nDivY,nDivZ))
  
    
  wavefunctions = 0

  allocate(potential(nDivX,nDivY,nDivZ,nTimeSteps))

  potential = 0
  do i=400,600
     !wavefunctions(i,1,1) = (0,1d0)*sin((i-4d0)/4d0)+ cos((i-4d0)/4d0)
     wavefunctions(i,1,1) = (0.5d0,0.5d0)
  end do
  wavefunctions = wavefunctions/sum(wavefunctions*conjg(wavefunctions))
  
  
  
  
!  if(reciprocalBase.eqv..FALSE.) then
  call real_space_qd(wavefunctions,potential,potentialTimeDependent, &
       endAtBoxEdge,nTimeSteps,deltat,deltax,accuracy)
!  else
!     call reciprocal_space_qd(wavefunctions,potential)
!  end if


! See below for the real/reciprocal space functions

end program QuantumDynamics






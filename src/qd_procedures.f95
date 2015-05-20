module qdprocedures

  implicit none
  private

  public real_space_qd



contains

  subroutine real_space_qd(wfinout,potential,potentialTimeDependent,endAtBoxEdge,nTimeSteps)

    complex*16, intent(inout) :: wfinout(:,:,:)
    real*8, intent(in)        :: potential(:,:,:,:)
    complex*16                :: kineticMatrix(size(wfinout,1),size(wfinout,2),size(wfinout,3), &
         size(wfinout,1),size(wfinout,2),size(wfinout,3)), &
         initialiser(size(wfinout,1),size(wfinout,2),size(wfinout,3))

    complex*16, allocatable   :: kmDiag(:), kmLinkTo(:,:), kmWeights(:,:), wf(:)
    
    integer*8, intent(in)     :: nTimeSteps
    logical, intent(in)       :: endAtBoxEdge, potentialTimeDependent

    integer*8                 :: i, j, k, timeStep, nDivXEff, nDivYEff, nDivZEff
    integer*8                 :: nDivX, nDivY, nDivZ
    integer*8                 :: maxlinks
    logical                   :: running

    maxlinks = 1
    if(nDivY.gt.1) then
       maxlinks = 2
       if(nDivZ.gt.1) then
          maxlinks = 3
       end if
    end if
    
    nDivXEff = size(wfinout,1)
    nDivYEff = size(wfinout,2)
    nDivZEff = size(wfinout,3)
    nDivX = nDivXEff
    nDivY = nDivYEff
    nDivZ = nDivZEff

    allocate( kmDiag(nDivXEff*nDivYEff*nDivZEff) )
    allocate( kmLinkTo(nDivXEff*nDivYEff*nDivZEff - 2  ,  maxlinks) )
    allocate( kmWeights((nDivXEff*nDivYEff*nDivZEff - 2  ,  maxlinks) )
    allocate( wf(nDivXEff*nDivYEff*nDivZEff) )

    kineticMatrix = 0
    kmDiag = 0
    kmLinkTo = 0
    kmWeights = 0
    

    if(potentialTimeDependent .eqv. .FALSE.) then

       forall(i=1:nDivX-1,j=1:nDivY,k=1:nDivZ)
          kineticMatrix(i,j,k,i,j,k) = 0.5 + potential(i,j,k,1)/2
          kineticMatrix(i+1,j,k,i,j,k) = -0.25d0
          kineticMatrix(i,j,k,i+1,j,k) = -0.25d0
       end forall
       forall(j=1:nDivY,k=1:nDivZ)
          kineticMatrix(nDivX,j,k,nDivX,j,k) = 0.5d0 + potential(nDivX,j,k,1)/2
       end forall
       if (nDivY.gt.1) then
          forall(j=1:nDivX,i=1:nDivY-1,k=1:nDivZ)
             kineticMatrix(j,i,k,j,i,k) = kineticMatrix(j,i,k,j,i,k) + 0.5d0
             kineticMatrix(j,i+1,k,j,i,k) = -0.25d0
             kineticMatrix(j,i,k,j,i+1,k) = -0.25d0
          end forall
          forall(j=1:nDivX,k=1:nDivZ)
             kineticMatrix(j,nDivY,k,j,nDivY,k) = kineticMatrix(j,nDivY,k,j,nDivY,k) + 0.5d0
          end forall
          if (nDivZ.gt.1) then
             forall(j=1:nDivX,k=1:nDivY,i=1:nDivZ-1)
                kineticMatrix(j,k,i,j,k,i) = kineticMatrix(j,k,i,j,k,i) + 0.5d0
                kineticMatrix(j,k,i+1,j,k,i) = -0.25d0
                kineticMatrix(j,k,i,j,k,i+1) = -0.25d0
             end forall
             forall(j=1:nDivX,k=1:nDivY)
                kineticMatrix(j,k,nDivZ,j,k,nDivZ) = kineticMatrix(j,k,nDivZ,j,k,nDivZ) + 0.5d0
          end if
       end if
       
    end if

    reshape(kineticMatrix,(/nDivX*nDivY*nDivZ,nDivX*nDivY*nDivZ/))
    wf = reshape(wfinout,(/nDivX*nDivY*nDivZ/))
    
!    kineticMatrix(nDivXEff,:,:,nDivXEff,:,:) = kineticMatrix(nDivXEff,:,:,nDivXEff,:,:) + 1d0
!    kineticMatrix(:,nDivYEff,:,:,nDivYEff,:) = kineticMatrix(:,nDivYEff,:,:,nDivYEff,:) + 1d0
!    kineticMatrix(:,:,nDivZEff,:,:,nDivZEff) = kineticMatrix(:,:,nDivZEff,:,:,nDivZEff) + 1d0
    
    call scrs_ify(reshape(kineticMatrix,(/nDivX*nDivY*nDivZ,nDivX*nDivY*nDivZ/)), &
         kmDiag,kmLinkTo,kmWeights,maxlinks)


    ! initialiser = A'x - Ax = (A'-A)x = -2ix ; this greatly simplifies the first steps of BiCGstab.
    initialiser = (0,-2d0)*wf
    
    do i=1,nDivX*nDivY*nDivZ
       kmDiag(i) = kmDiag(i) + (0,1d0)
    end do

    call do_bicgstab(wf,initialiser,kmDiag,kmLinkTo,kmWeights,maxlinks)

       
    
      
    
    
    
          
       !       do i=1,(nDivXEff-1)
!          forall(j=1:nDivYEff, k=1:nDivZEff)
!             kineticMatrix(i,j,k,i,j,k) = kineticMatrix(i,j,k,i,j,k) + 1d0
!             kineticMatrix(i,j,k,i,j,k) = kineticMatrix(i,j,k,i,j,k) + potential(i,j,k,1),(/nDivYEff,nDivZEff,nDivYEff,nDivZEff/))
!          end forall
!          kineticMatrix(i,:,:,i+1,:,:) = kineticMatrix(i,:,:,i+1,:,:) - 5d-1
!          kineticMatrix(i+1,:,:,i,:,:) = kineticMatrix(i+1,:,:,i,:,:) - 5d-1
!       end do
!       if (nDivYEff > 3) then
!          do i=1,(nDivYEff-1)
!             forall(j=1:nDivXEff, k=1:nDivZEff)
!                kineticMatrix(j,i,k,j,i,k) = kineticMatrix(j,i,k,j,i,k) + 1d0
!                kineticMatrix(j,i,k,j,i,k) = kineticMatrix(j,i,k,j,i,k) + potential(j,i,k,1),(/nDivYEff,nDivZEff,nDivYEff,nDivZEff/))
!             end forall
!             kineticMatrix(:,i,:,:,i+1,:) = kineticMatrix(:,i,:,:,i+1,:) - 5d-1
!             kineticMatrix(:,i+1,:,:,i,:) = kineticMatrix(:,i+1,:,:,i,:) - 5d-1
!          end do
!          if (nDivZEff > 3) then
!             do i=1,(nDivZEff-1)
!                forall(j=1:nDivXEff, k=1:nDivYEff)
!                   kineticMatrix(j,k,i,j,k,i) = kineticMatrix(i,j,k,i,j,k) + 1d0
!                   kineticMatrix(j,k,i,j,k,i) = kineticMatrix(i,j,k,i,j,k) + potential(i,j,k,1),(/nDivYEff,nDivZEff,nDivYEff,nDivZEff/))
!                end forall
!                kineticMatrix(:,i,:,:,i+1,:) = kineticMatrix(:,i,:,:,i+1,:) - 5d-1
!                kineticMatrix(:,:,i+1,:,:,i) = kineticMatrix(:,:,i+1,:,:,i) - 5d-1
!             end do
!          end if
!       end if

  
  
  kineticMatrix(nDivXEff,:,:,nDivXEff,:,:) = kineticMatrix(nDivXEff,:,:,nDivXEff,:,:) + 1d0
  kineticMatrix(:,nDivYEff,:,:,nDivYEff,:) = kineticMatrix(:,nDivYEff,:,:,nDivYEff,:) + 1d0
  kineticMatrix(:,:,nDivZEff,:,:,nDivZEff) = kineticMatrix(:,:,nDivZEff,:,:,nDivZEff) + 1d0

  kineticMatrix = kineticMatrix*(0d0,1d0)
  ! This made the matrix completely imaginary.
  
  ! After making the kinetic (or total) matrix, do a wavefunction evolution loop:
  print *,kineticMatrix
  running = .TRUE.
  
!  do while (indefSimulation .OR. running)
     ! This might not look too effective, but the potentialTimeDependent condition will be caught
     ! during compilation in gfortran (as seen during debugging)
!     if(potentialTimeDependent == .FALSE.) then
!        call wave_function_evolution()
        

!     if(timeStep == nTimeSteps) then
!        running = .FALSE.
!     end if

!  end do
  
  
  
  
end subroutine real_space_qd
















end module qdprocedures

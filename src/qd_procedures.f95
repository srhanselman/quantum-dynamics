module qdprocedures

use linalgalgorithms
use scrstools

  implicit none
  private

  public real_space_qd



contains

  subroutine real_space_qd(wfinout,potentialin,potentialTimeDependent,endAtBoxEdge,nTimeSteps,deltat, &
       deltax,accuracy)

    complex*16, intent(inout) :: wfinout(:,:,:)
    real*8, intent(in)        :: potentialin(:,:,:,:), deltat, deltax, accuracy
    real*8                    :: deltaxinvsq
    complex*16                :: kineticMatrix(size(wfinout,1),size(wfinout,2),size(wfinout,3), &
         size(wfinout,1),size(wfinout,2),size(wfinout,3))

    complex*16, allocatable   :: kmDiag(:), kmWeights(:,:), wf(:), initialiser(:)
    integer*8, allocatable    :: kmLinkTo(:,:)
    real*8, allocatable       :: potential(:,:)
    
    integer*8, intent(in)     :: nTimeSteps
    logical, intent(in)       :: endAtBoxEdge, potentialTimeDependent

    integer*8                 :: i, j, k, timeStep, nDivXEff, nDivYEff, nDivZEff, t
    integer*8                 :: nDivX, nDivY, nDivZ
    integer*8                 :: maxlinks
    logical                   :: running

    character(len=20)         :: fmt
    
    nDivXEff = size(wfinout,1)
    nDivYEff = size(wfinout,2)
    nDivZEff = size(wfinout,3)
    nDivX = nDivXEff
    nDivY = nDivYEff
    nDivZ = nDivZEff

    write(fmt,'(a, i0, a)') '(I8, ', nDivX*nDivY*nDivZ, '(e20.14,a))'

    if(endAtBoxEdge.eqv..TRUE.) then
       maxlinks = 2
    else
       maxlinks = 1
    end if
    
    if(nDivY.gt.1) then
       if (endAtBoxEdge.eqv..TRUE.) then   
          maxlinks = maxlinks + 2
       else
          maxlinks = maxlinks + 1
       end if
       
       if(nDivZ.gt.1) then
          if (endAtBoxEdge.eqv..TRUE.) then
             maxlinks = maxlinks + 2
          else
             maxlinks = maxlinks + 1
          end if
       end if
    end if
    

    deltaxinvsq = 1/(deltax*deltax)

    allocate( kmDiag(nDivXEff*nDivYEff*nDivZEff) )
    allocate( kmLinkTo(nDivXEff*nDivYEff*nDivZEff - 1 ,  maxlinks) )
    allocate( kmWeights(nDivXEff*nDivYEff*nDivZEff - 1 ,  maxlinks) )
    allocate( wf(nDivXEff*nDivYEff*nDivZEff) )
    allocate( initialiser(nDivXEff*nDivYEff*nDivZEff) )
    allocate( potential(nDivXEff*nDivYEff*nDivZEff , size(potentialin,4) ) )

    kineticMatrix = 0
    kmDiag = 0
    kmLinkTo = 0
    kmWeights = 0

    ! As the potentials are diagonal terms, the initial matrix only contains the kinetic energy;
    ! the time step will be implemented afterwards, while the potential is updated more easily.


    forall(i=1:nDivX-1,j=1:nDivY,k=1:nDivZ)
       kineticMatrix(i,j,k,i,j,k) = 1d0*deltaxinvsq
       kineticMatrix(i+1,j,k,i,j,k) = -0.5d0*deltaxinvsq
       kineticMatrix(i,j,k,i+1,j,k) = -0.5d0
    end forall
    forall(j=1:nDivY,k=1:nDivZ)
       kineticMatrix(nDivX,j,k,nDivX,j,k) = 1d0*deltaxinvsq
    end forall
    if (nDivY.gt.1) then
       forall(j=1:nDivX,i=1:nDivY-1,k=1:nDivZ)
          kineticMatrix(j,i,k,j,i,k) = kineticMatrix(j,i,k,j,i,k) + 1d0*deltaxinvsq
          kineticMatrix(j,i+1,k,j,i,k) = -0.5d0*deltaxinvsq
          kineticMatrix(j,i,k,j,i+1,k) = -0.5d0*deltaxinvsq
       end forall
       forall(j=1:nDivX,k=1:nDivZ)
          kineticMatrix(j,nDivY,k,j,nDivY,k) = kineticMatrix(j,nDivY,k,j,nDivY,k) + 1d0*deltaxinvsq
       end forall
       if (nDivZ.gt.1) then
          forall(j=1:nDivX,k=1:nDivY,i=1:nDivZ-1)
             kineticMatrix(j,k,i,j,k,i) = kineticMatrix(j,k,i,j,k,i) + 1d0*deltaxinvsq
             kineticMatrix(j,k,i+1,j,k,i) = -0.5d0*deltaxinvsq
             kineticMatrix(j,k,i,j,k,i+1) = -0.5d0*deltaxinvsq
          end forall
          forall(j=1:nDivX,k=1:nDivY)
             kineticMatrix(j,k,nDivZ,j,k,nDivZ) = kineticMatrix(j,k,nDivZ,j,k,nDivZ) + 1d0*deltaxinvsq
          end forall
       end if
    end if
       
    if (endAtBoxEdge.eqv..FALSE.) then
       forall(j=1:nDivY,k=1:nDivZ)
          kineticMatrix(1,j,k,nDivX,j,k) = kineticMatrix(1,j,k,nDivX,j,k) - 0.5d0*deltaxinvsq
          kineticMatrix(nDivX,j,k,1,j,k) = kineticMatrix(nDivX,j,k,1,j,k) - 0.5d0*deltaxinvsq
       end forall
       if (nDivY.gt.1) then
          forall(j=1:nDivX,k=1:nDivZ)
             kineticMatrix(j,1,k,j,nDivX,k) = kineticMatrix(j,1,k,j,nDivX,k) - 0.5d0*deltaxinvsq
             kineticMatrix(j,nDivX,k,j,1,k) = kineticMatrix(j,nDivX,k,j,1,k) - 0.5d0*deltaxinvsq
          end forall
          if (nDivZ.gt.1) then
             forall(j=1:nDivX,k=1:nDivY)
                kineticMatrix(j,k,1,j,k,nDivX) = kineticMatrix(j,k,1,j,k,nDivX) - 0.5d0*deltaxinvsq
                kineticMatrix(j,k,nDivX,j,k,1) = kineticMatrix(j,k,nDivX,j,k,1) - 0.5d0*deltaxinvsq
             end forall
          end if
       end if
    end if
    
          


    
    wf = reshape(wfinout,(/nDivX*nDivY*nDivZ/))
    potential = reshape(potentialin,(/nDivX*nDivY*nDivZ , nTimeSteps/))
 
       
    
!    kineticMatrix(nDivXEff,:,:,nDivXEff,:,:) = kineticMatrix(nDivXEff,:,:,nDivXEff,:,:) + 1d0
!    kineticMatrix(:,nDivYEff,:,:,nDivYEff,:) = kineticMatrix(:,nDivYEff,:,:,nDivYEff,:) + 1d0
!    kineticMatrix(:,:,nDivZEff,:,:,nDivZEff) = kineticMatrix(:,:,nDivZEff,:,:,nDivZEff) + 1d0
    
    call scrs_ify(reshape(kineticMatrix,(/nDivX*nDivY*nDivZ,nDivX*nDivY*nDivZ/)), &
         kmDiag,kmLinkTo,kmWeights,maxlinks)

    ! Equation: (1+iTH(t)/2)x_t = (1-iTH(t-1)/2)x_{t-1} -> (H(t) - 2i/T)x_t = (H(t-1) + 2i/T)x_{t-1}
    ! initialiser = (H(t-1) + 2i/T)x - (H(t) - 2i/T)x = (H(t-1)-H(t)+2/T)x = (0.5V(t-1)-0.5V(t)+2/T).
    initialiser = (0,2d0)*wf/deltat

    if(potentialTimeDependent.eqv..FALSE.) then
       do i=1,nDivX*nDivY*nDivZ
          kmDiag(i) = kmDiag(i) + (0,2d0)/deltat
       end do
    else
       do i=1,nDivX*nDivY*nDivZ
          kmDiag(i) = kmDiag(i) + (0,2d0)/deltat + potential(i,1)
       end do
    end if
    
    call do_bicgstab(wf,initialiser,kmDiag,kmLinkTo,kmWeights,accuracy)

    if(potentialTimeDependent.eqv..FALSE.) then
       kmDiag(:) = kmDiag(:) + potential(:,1)
    end if
    
    
    timeloop: do t=2,nTimeSteps
       if(potentialTimeDependent.eqv..TRUE.) then
          kmDiag(:) = kmDiag(:) + 0.25d0*potential(:,(t-1)) + 0.75d0*potential(:,t)
          initialiser = initialiser*(0.5d0*potential(:,(t-1)) - 0.5d0*potential(:,t) + 2d0/deltat)
       else
          initialiser = initialiser*2d0/deltat
       end if

       call do_bicgstab(wf, initialiser,kmDiag,kmLinkTo,kmWeights,accuracy)

       write(1,fmt) t,( REAL(wf(j)),',', j=1,size(wf,1) )
       write(1,fmt) t,( AIMAG(wf(j)),',', j=1,size(wf,1) )
       
       if(potentialTimeDependent.eqv..TRUE.) then
          kmDiag(:) = kmDiag(:) - 0.25d0*potential(:,t-1) - 0.75d0*potential(:,t)
       end if
       
    end do timeloop
    
  
  
  
end subroutine real_space_qd
















end module qdprocedures

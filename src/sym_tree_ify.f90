module SymTreeIfy

  implicit none
  private
  public scrs_ify
  public mult_scrs_vec


contains

  subroutine scrs_ify(matrix,diagonal,linkto,weight,maxlinks)
    ! This subroutine isolates the diagonal from the symmetric matrix, and subsequently 
    complex*16, intent(in)  :: matrix(:,:)
    complex*16, intent(out) :: diagonal(:), weights(:)
    integer*8, intent(in)   :: maxlinks
    integer*8, intent(out)  :: linkto(:)
    integer*8               :: i,j,k,N

    N = size(matrix,1)
    k = 1
    linkto=1
    weight=0

    do i=1,N
       diagonal(i) = matrix(i,i)
    end do

    ! Indeed I will break inner loop 
    do i=1,(N-2)
       innerloop: do j=1,(N-1)
          if (matrix(i,j).ne.0) then
             linkto(i,k)=j
             weights(i,k)=matrix(i,j)
             k = k + 1
             if (k.gt.maxlinks) then
                exit innerloop
             end if
          end if
       end do innerloop
    end do

    weights(N-1,1) = matrix(N-1,N)
  end subroutine scrs_ify

  
  
  subroutine mult_scrs_vec(diagonal,linkto,weights,vec,vecnew)

    ! This subroutine isolates the diagonal from the symmetric matrix, and subsequently

    complex*16, intent(in)  :: diagonal(:), weights(:,:)
    integer*8, intent(in)   :: linkto(:)
    complex*16, intent(in)  :: vec(:)
    complex*16, intent(out) :: vecnew(size(vec,1))
    integer*8               :: i,j,k,N

    vecnew = vec*diagonal
    
    N = size(vec,1)
    

    do i=1,(N-1)
       do j=1,size(weights,2)
          vecnew(linkto(i,j)) = vecnew(linkto(i,j)) + weights(i,j) * vec(i)
          vecnew(i) = vecnew(i) + weights(i,j) * vec(linkto(i,j))
       end do
    end do
    ! Order: m(N-1), prefactor = 2*(add+mult)
       

  end subroutine mult_scrs_vec
  
  





end module SymTreeIfy

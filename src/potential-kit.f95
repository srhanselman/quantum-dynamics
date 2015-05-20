module potentialkit

  public potential_insert_block
!  public potential_insert_ellipse
  public potential_insert_coulomb
  public potential_insert_slope



  subroutine potential_insert_block(start,end,depth,sl,potential)

    real*8, intent(in)    :: start(:), end(:), depth, sl
    integer               :: startelement(:),endelement(:),i,j,k

    real*8, intent(inout) :: potential

    startelement = nint(start/sl)
    endelement = nint(start/sl)
    do i=startelement(1),endelement(1)
       do j=startelement(2),endelement(2)
          do k=startelement(3),endelement(3)
             potential(:,i,j,k) = potential(:,i,j,k) + depth
          end do
       end do
    end do

  end subroutine potential_insert_block


  subroutine potential_insert_coulomb(center(dim),depth,sl,potential)

    real*8, intent(in)   :: center(:), radius(:), depth, sl
    integer 

    





















end module potentialkit

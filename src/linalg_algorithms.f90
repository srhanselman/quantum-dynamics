module linalgalgorithms

  implicit none
  private
  public do_bicgstab



contains

  subroutine do_bicgstab(x,r,diag,linkto,weights,accuracy)

    complex*16, intent(inout) :: x(:)
    complex*16                :: xmult(size(x,1)), beta, rho, rholast, nu, rinit, rlast
    complex*16, intent(inout) :: r(:)

    complex*16, intent(in)    :: diag(:), linkto(:,:), weights(:,:)


    ! First do an initial estimate. For smoothly evolving wavefunctions x_0 should be fine;
    ! this greatly simplifies proceedings as it allows to replace b by an initialiser
    ! (as initialiser = A'x - Ax = (A'-A)x = -2ix).
    rho   = 1
    alpha = 1
    omega = 1
    beta  = 1
    rinit = r
    p = r

    ! Now the big loop is ready for action:
    
    do
       call mult_scrs_vec(diag,linkto,weights,p,nu)
       alpha = rho/(dot_product(rinit,nu))
       s = rlast - alpha*nu
       call mult_src_vec(diag,linkto,weights,s,t)
       omega = dot_product(t,s)/dot_product(t,t)
       
       x = x + alpha*p + omega*s
       
       r = s - omega*t

       rholast = rho
       rho = dot_product(rinit,r)
       beta = rho*alpha/(rholast*omega)

       p = beta*p
       p = p - beta*omega*nu + r
    end do
    

  end subroutine do_bicgstab
  






end module linalgalgorithms

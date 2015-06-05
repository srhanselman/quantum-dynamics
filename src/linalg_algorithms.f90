module linalgalgorithms

use scrstools
  
  implicit none
  private
  public do_bicgstab



contains

  subroutine do_bicgstab(x,r,diag,linkto,weights,accuracy)

   complex*16, intent(inout) :: x(:)
   complex*16              :: xmult(size(x)),beta,rho,rholast,nu(size(x)), &
        rinit(size(x)),rlast(size(x)),alpha,omega,p(size(x)),s(size(x)), &
        t(size(x))
   complex*16, intent(inout) :: r(:)

   complex*16, intent(in)    :: diag(:), weights(:,:)
   integer*8, intent(in)     :: linkto(:,:)
   real*8, intent(in)        :: accuracy

   ! First do an initial estimate. For smoothly evolving wavefunctions x_0 should be fine;
   ! this greatly simplifies proceedings as it allows to replace b by an initialiser
   ! (as initialiser = A'x - Ax = (A'-A)x = -2ix).
   ! Also note that the first steps are performed BEFORE the iteration loop (to make matters
   ! more efficient); 
    rho   = dot_product(r,r)
    alpha = 1
    omega = 1
    beta  = rho
    rinit = r
    p = r

    ! Now the big loop is ready for action:
    
    iteration: do
       call mult_scrs_vec(diag,linkto,weights,p,nu)
       alpha = rho/(dot_product(rinit,nu))
       s = r - alpha*nu
       call mult_scrs_vec(diag,linkto,weights,s,t)
       omega = dot_product(t,s)/dot_product(t,t)

       if (all(abs(alpha*p+omega*s).lt.accuracy)) then
          x=x+alpha*p+omega*s
          print *, x(500)
          exit iteration
       end if
       
      
       x = x + alpha*p + omega*s
       print *,x(500)
       r = s - omega*t

       rholast = rho
       rho = dot_product(rinit,r)
       beta = rho*alpha/(rholast*omega)

       p = beta*p
       p = p - beta*omega*nu + r
    end do iteration

    print *,"done iterating"
    
  end subroutine do_bicgstab
  





end module linalgalgorithms

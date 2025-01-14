module mod_debug
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_param     , only: dims
  use mod_types
  implicit none
  private
  public chk_mean,chk_helmholtz
  contains
  subroutine chk_mean(n,grid_vol_ratio,p,mean)
    !
    ! compute the mean value of an observable over the entire domain
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:) :: grid_vol_ratio
    real(rp), intent(in), dimension(0:,0:,0:) :: p
    real(rp), intent(out) :: mean
    integer :: i,j,k
    mean = 0.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,grid_vol_ratio) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP REDUCTION(+:mean)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          mean = mean + p(i,j,k)*grid_vol_ratio(k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,mean,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine chk_mean
  !
  subroutine chk_helmholtz(lo,hi,dli,dzci,dzfi,alpha,fp,fpp,bc,is_bound,c_or_f,diffmax)
    !
    ! this subroutine checks if the implementation of implicit diffusion is
    ! correct under sanity.f90
    !
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(2) :: dli
    real(rp), intent(in) :: alpha
    real(rp), intent(in), dimension(lo(3)-1:) :: dzfi,dzci
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: fp,fpp
    character(len=1), intent(in), dimension(0:1,3) :: bc
    character(len=1), intent(in), dimension(3) :: c_or_f
    logical         , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(out) :: diffmax
    real(rp) :: val
    integer :: i,j,k
    integer :: idir
    integer, dimension(3) :: q
    q(:) = 0
    do idir = 1,3
      if(bc(1,idir).ne.'P'.and.c_or_f(idir) == 'f'.and.is_bound(1,idir)) q(idir) = 1
    end do
    select case(c_or_f(3))
    case('c')
      diffmax = 0.
      do k=lo(3),hi(3)-q(3)
        do j=lo(2),hi(2)-q(2)
          do i=lo(1),hi(1)-q(1)
            val =  fpp(i,j,k)+(1./alpha)*( &
                  (fpp(i+1,j,k)-2.*fpp(i,j,k)+fpp(i-1,j,k))*(dli(1)**2) + &
                  (fpp(i,j+1,k)-2.*fpp(i,j,k)+fpp(i,j-1,k))*(dli(2)**2) + &
                 ((fpp(i,j,k+1)-fpp(i,j,k  ))*dzci(k ) - &
                  (fpp(i,j,k  )-fpp(i,j,k-1))*dzci(k-1))*dzfi(k) )
            val = val*alpha
            diffmax = max(diffmax,abs(val-fp(i,j,k)))
            !if(abs(val-fp(i,j,k)) > 1.e-8) print*, 'Large difference : ', val-fp(i,j,k),i,j,k
          end do
        end do
      end do
    case('f')
      diffmax = 0.
      do k=lo(3),hi(3)-q(3)
        do j=lo(2),hi(2)-q(2)
          do i=lo(1),hi(1)-q(1)
            val =  fpp(i,j,k)+(1./alpha)*( &
                  (fpp(i+1,j,k)-2.*fpp(i,j,k)+fpp(i-1,j,k))*(dli(1)**2) + &
                  (fpp(i,j+1,k)-2.*fpp(i,j,k)+fpp(i,j-1,k))*(dli(2)**2) + &
                 ((fpp(i,j,k+1)-fpp(i,j,k  ))*dzfi(k+1) - &
                  (fpp(i,j,k  )-fpp(i,j,k-1))*dzfi(k ))*dzci(k) )
            val = val*alpha
            diffmax = max(diffmax,abs(val-fp(i,j,k)))
            !if(abs(val-fp(i,j,k)) > 1.e-8) print*, 'Large difference : ', val,fp(i,j,k),i,j,k
          end do
        end do
      end do
    end select
    call mpi_allreduce(MPI_IN_PLACE,diffmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  end subroutine chk_helmholtz
end module mod_debug

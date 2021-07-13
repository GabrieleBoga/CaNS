module mod_force_time_dependent
  use mod_common_mpi, only: ijk_start
  use mod_types
  use mod_param, only: pi
  implicit none
  private
  public force_time_dependent
  contains
  subroutine force_time_dependent(forcing,f,kf,abc,n,l,dl,zc,dt,up,vp,wp)
    character(len=3), intent(in)       :: forcing
    real(rp), intent(in)               :: f,kf
    real(rp), intent(in), dimension(3) :: abc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in) :: dt
    real(rp), intent(inout), dimension(0:,0:,0:) :: up,vp,wp
    select case(forcing)
    case('abc')
      call force_time_dependent_abc(f,kf,abc,n,l,dl,zc,dt,up,vp,wp)
    case default
      ! do nothing
    end select
    return
  end subroutine force_time_dependent
  !
  subroutine force_time_dependent_abc(f,kf,abc,n,l,dl,zc,dt,up,vp,wp)
    !
    ! simple large-scale ABC forcing
    !
    implicit none
    real(rp), intent(in)               :: f,kf
    real(rp), intent(in), dimension(3) :: abc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in) :: dt
    real(rp), intent(inout), dimension(0:,0:,0:) :: up,vp,wp
    real(rp) :: xcl,ycl,zcl,fx,fy,fz
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,f,abc,kf,l,dl,zc,ijk_start,up,vp,wp,dt) &
    !$OMP PRIVATE(i,j,k,xcl,ycl,zcl,fx,fy,fz)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          zcl =                     zc(k)*2.*pi/l(3)
          ycl = (j+ijk_start(2)-.5)*dl(2)*2.*pi/l(2)
          xcl = (i+ijk_start(1)-.5)*dl(1)*2.*pi/l(1)
          fx  = f*(kf*2.*pi/l(1))**2*(abc(1)*sin(kf*zcl) + abc(3)*cos(kf*ycl))
          fy  = f*(kf*2.*pi/l(2))**2*(abc(2)*sin(kf*xcl) + abc(1)*cos(kf*zcl))
          fz  = f*(kf*2.*pi/l(3))**2*(abc(3)*sin(kf*ycl) + abc(2)*cos(kf*xcl))
          up(i,j,k) = up(i,j,k) + fx*dt
          vp(i,j,k) = vp(i,j,k) + fy*dt
          wp(i,j,k) = wp(i,j,k) + fz*dt
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine force_time_dependent_abc
end module mod_force_time_dependent

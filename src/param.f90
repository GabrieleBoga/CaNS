module mod_param
use mod_types
implicit none
public
!
! parameters
!
real(rp), parameter :: pi = acos(-1._rp)
real(rp), parameter :: small = epsilon(pi)*10**(precision(pi)/2)
character(len=100), parameter :: datadir = 'data/'
real(rp), parameter, dimension(2,3) :: rkcoeff = reshape([32._rp/60._rp,  0._rp        , &
                                                          25._rp/60._rp, -17._rp/60._rp, &
                                                          45._rp/60._rp, -25._rp/60._rp], shape(rkcoeff))
real(rp), parameter, dimension(3)   :: rkcoeff12 = rkcoeff(1,:)+rkcoeff(2,:)
!
! variables to be determined from the input file 'dns.in'
!
integer :: itot,jtot,ktot
real(rp) :: lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,gr
real(rp) :: cfl,dtmin
real(rp) :: uref,lref,rey,visc
!
character(len=100) :: inivel
logical :: is_wallturb
!
integer :: nstep
real(rp) :: time_max,tw_max
logical, dimension(3) :: stop_type
logical :: restart,is_overwrite_save
integer :: nsaves_max
integer :: icheck,iout0d,iout1d,iout2d,iout3d,isave
!
integer, dimension(2) :: dims
integer :: nthreadsmax
!
integer, dimension(0:1,3) :: nb
logical, dimension(0:1,3) :: is_bound
character(len=1), dimension(0:1,3,3) ::  cbcvel
real(rp)        , dimension(0:1,3,3) :: bcvel
character(len=1), dimension(0:1,3)   ::  cbcpre
real(rp)        , dimension(0:1,3)   ::   bcpre
!
real(rp), dimension(3) :: bforce
logical , dimension(3) :: is_forced
real(rp), dimension(3) :: velf
!
integer , dimension(3) :: ng
integer , dimension(3) :: n,n_z,lo,hi
real(rp), dimension(3) :: l
real(rp), dimension(3) :: dl
real(rp), dimension(3) :: dli
!
contains
  subroutine read_input(myid)
  use mpi
  implicit none
  integer, intent(in) :: myid
  integer :: iunit,ierr
    nsaves_max = 0 ! a good default, for backward compatibility
    open(newunit=iunit,file='dns.in',status='old',action='read',iostat=ierr)
      if( ierr == 0 ) then
        read(iunit,*) itot,jtot,ktot
        read(iunit,*) lx,ly,lz
        read(iunit,*) gr
        read(iunit,*) cfl,dtmin
        read(iunit,*) uref,lref,rey
        read(iunit,*) inivel
        read(iunit,*) is_wallturb
        read(iunit,*) nstep, time_max,tw_max
        read(iunit,*) stop_type(1),stop_type(2),stop_type(3)
        read(iunit,*) restart,is_overwrite_save,nsaves_max
        read(iunit,*) icheck,iout0d,iout1d,iout2d,iout3d,isave
        read(iunit,*) cbcvel(0,1,1),cbcvel(1,1,1),cbcvel(0,2,1),cbcvel(1,2,1),cbcvel(0,3,1),cbcvel(1,3,1)
        read(iunit,*) cbcvel(0,1,2),cbcvel(1,1,2),cbcvel(0,2,2),cbcvel(1,2,2),cbcvel(0,3,2),cbcvel(1,3,2)
        read(iunit,*) cbcvel(0,1,3),cbcvel(1,1,3),cbcvel(0,2,3),cbcvel(1,2,3),cbcvel(0,3,3),cbcvel(1,3,3)
        read(iunit,*) cbcpre(0,1  ),cbcpre(1,1  ),cbcpre(0,2  ),cbcpre(1,2  ),cbcpre(0,3  ),cbcpre(1,3  )
        read(iunit,*)  bcvel(0,1,1), bcvel(1,1,1), bcvel(0,2,1), bcvel(1,2,1), bcvel(0,3,1), bcvel(1,3,1)
        read(iunit,*)  bcvel(0,1,2), bcvel(1,1,2), bcvel(0,2,2), bcvel(1,2,2), bcvel(0,3,2), bcvel(1,3,2)
        read(iunit,*)  bcvel(0,1,3), bcvel(1,1,3), bcvel(0,2,3), bcvel(1,2,3), bcvel(0,3,3), bcvel(1,3,3)
        read(iunit,*)  bcpre(0,1  ), bcpre(1,1  ), bcpre(0,2  ), bcpre(1,2  ), bcpre(0,3  ), bcpre(1,3  )
        read(iunit,*)  bforce(1),bforce(2),bforce(3)
        read(iunit,*)  is_forced(1),is_forced(2),is_forced(3)
        read(iunit,*)  velf(1),velf(2),velf(3)
        read(iunit,*) dims(1),dims(2)
        read(iunit,*) nthreadsmax
      else
        if(myid == 0) print*, 'Error reading the input file'
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        error stop
      end if
    close(iunit)
    dx = lx/(1.*itot)
    dy = ly/(1.*jtot)
    dz = lz/(1.*ktot)
    dxi = dx**(-1)
    dyi = dy**(-1)
    dzi = dz**(-1)
    !
    visc = uref*lref/rey
    ng  = [itot,jtot,ktot]
    l   = [lx,ly,lz]
    dl  = [dx,dy,dz]
    dli = [dxi,dyi,dzi]
  end subroutine read_input
end module mod_param

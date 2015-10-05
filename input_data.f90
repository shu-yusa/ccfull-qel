module input_data
  private

  type, public :: inp
    integer :: coup_p, Nch_p
    integer :: coup_t, Nch_t
    integer :: coup,   Nch
    integer :: Nphonon_p,  lambda_p,  nrot_p
    integer :: Nphonon_t,  lambda_t,  nrot_t
    integer :: Nphonon_p2, lambda_p2, nrot_p2
    integer :: Nphonon_t2, lambda_t2, nrot_t2
    integer :: Nphonon,   lambda,   nrot
    integer :: rgrid, Egrid, Jmax, di, nth, rgrid_cut
    integer :: num_stab_pt
    real(8) :: Ap, Zp, r0p, Rp
    real(8) :: At, Zt, r0t, Rt
    real(8) :: omega_p, betaL_p, betaL_pC
    real(8) :: omega_t, betaL_t, betaL_tC
    real(8) :: omega_p2, betaL_p2, betaL_p2C
    real(8) :: omega_t2, betaL_t2, betaL_t2C
    real(8) :: omega,   betaL, betaLC
    real(8) :: E2p, beta2p, beta4p
    real(8) :: E2t, beta2t, beta4t
    real(8) :: E2,  beta2,  beta4
    real(8) :: Emin,  Emax,  dE
    real(8) :: rmin,  rmax,  dr, rcut
    real(8) :: thmin, thmax, dth
    real(8) :: V0, r0, a, Rn
    real(8) :: W0, r0w, aw, Rw
    real(8) :: rmass, Rm, Rc, r0c, del_E

    contains
    
    procedure :: read_input
  end type
  
  contains
!----------------------------------------------------------------------
  subroutine read_input(ip, Fname, dir)
    use global_constant, only : mass
    implicit none
    class(inp), intent(out) :: ip
    real(8), parameter :: e=0.333333333333333d0
    character(len=*), intent(in) :: Fname
    character(len=*), intent(out) :: dir
    integer :: lines

    open(7, file=Fname, status='old', action='read')
    read(7,*) ip%Ap, ip%Zp, ip%At, ip%Zt
    read(7,*) ip%r0p, ip%r0t
    read(7,*) ip%coup_p, ip%coup_t
    read(7,*) ip%Nphonon_t, ip%lambda_t
    read(7,*) ip%omega_t, ip%betaL_t, ip%betaL_tC
    read(7,*) ip%Nphonon_t2, ip%lambda_t2
    read(7,*) ip%omega_t2, ip%betaL_t2, ip%betaL_t2C
    read(7,*) ip%nrot_t
    read(7,*) ip%E2t, ip%beta2t, ip%beta4t
    read(7,*) ip%Nphonon_p, ip%lambda_p
    read(7,*) ip%omega_p, ip%betaL_p, ip%betaL_pC
    read(7,*) ip%Nphonon_p2, ip%lambda_p2
    read(7,*) ip%omega_p2, ip%betaL_p2, ip%betaL_p2C
    read(7,*) ip%nrot_p
    read(7,*) ip%E2p, ip%beta2p, ip%beta4p
    read(7,*) ip%V0, ip%r0, ip%a
    read(7,*) ip%W0, ip%r0w, ip%aw
    read(7,*) ip%r0c
    read(7,*) ip%Emin, ip%Emax, ip%dE
    read(7,*) ip%thmin, ip%thmax, ip%dth
    read(7,*) ip%rmin, ip%rmax, ip%dr, ip%rcut
    read(7,*) ip%Jmax
    read(7,*) ip%num_stab_pt
    read(7,*) ip%del_E
    read(7,*) dir
    close(7)

    if (ip%Ap <= 0.0d0) stop 'Ap must be positive.'
    if (ip%Zp <= 0.0d0) stop 'Zp must be positive.'
    if (ip%At <= 0.0d0) stop 'At must be positive.'
    if (ip%Zt <= 0.0d0) stop 'Zt must be positive.'
    if (ip%r0p <= 0.0d0) stop 'r0p must be positive.'
    if (ip%r0t <= 0.0d0) stop 'r0t must be positive.'
    if (ip%Nphonon_p  < 0) stop 'Nphonon_p cannot be negative.'
    if (ip%Nphonon_t  < 0) stop 'Nphonon_t cannot be negative.'
    if (ip%Nphonon_p2 < 0) stop 'Nphonon_p2 cannot be negative.'
    if (ip%Nphonon_t2 < 0) stop 'Nphonon_t2 cannot be negative.'
    if (ip%lambda_p   < 0) stop 'lambda_p cannot be negative.'
    if (ip%lambda_t   < 0) stop 'lambda_t cannot be negative.'
    if (ip%lambda_p2  < 0) stop 'lambda_p2 cannot be negative.'
    if (ip%lambda_t2  < 0) stop 'lambda_t2 cannot be negative.'
    if (ip%omega_t  < 0) stop 'omega_t cannot be negative.'
    if (ip%omega_t2 < 0) stop 'omega_t2 cannot be negative.'
    if (ip%omega_p  < 0) stop 'omega_p cannot be negative.'
    if (ip%omega_p2 < 0) stop 'omega_p2 cannot be negative.'
    if (ip%nrot_t < 0) stop 'nrot_t cannot be negative.'
    if (ip%nrot_p < 0) stop 'nrot_p cannot be negative.'
    if (ip%E2t < 0) stop 'E2t cannot be negative.'
    if (ip%E2p < 0) stop 'E2p cannot be negative.'
    if (ip%r0 < 0) stop 'r0 cannot be negative.'
    if (ip%r0w < 0) stop 'r0w cannot be negative.'
    if (ip%r0c < 0) stop 'r0c cannot be negative.'
    if (ip%a < 0) stop 'a cannot be negative.'
    if (ip%aw < 0) stop 'aw cannot be negative.'
    if (ip%Jmax < 0) stop 'Jmax cannot be negative.'
    if (ip%Emin < 0) stop 'Emin cannot be negative.'
    if (ip%Emax < 0) stop 'Emax cannot be negative.'
    if (ip%dE < 0) stop 'dE cannot be negative.'
    if (ip%rmin < 0) stop 'rmin cannot be negative.'
    if (ip%rmax < 0) stop 'rmax cannot be negative.'
    if (ip%rcut < 0) stop 'rcut cannot be negative.'
    if (ip%dr < 0) stop 'dr cannot be negative.'
    if (ip%del_E <= 0) stop 'del_E must be positive.'

    select case(ip%coup_p)
      case(-1) ; ip%Nch_p = 1
      case(0)  ; ip%Nch_p = ip%Nphonon_p + ip%Nphonon_p2 + 1
      case(1)  ; ip%Nch_p = ip%nrot_p + ip%Nphonon_p2 + 1
      case default ; stop 'Specify the couping of the projectile correctly.'
    end select
    select case(ip%coup_t)
      case(-1) ; ip%Nch_t = 1
      case(0)  ; ip%Nch_t = ip%Nphonon_t + ip%Nphonon_t2 + 1
      case(1)  ; ip%Nch_t = ip%nrot_t + ip%Nphonon_t2 + 1
      case default ; stop 'Specify the couping of the target correctly.'
    end select
    ip%Nch = ip%Nch_p * ip%Nch_t


    ip%Rn = ip%r0  * (ip%Ap ** e + ip%At ** e)
    ip%Rw = ip%r0w * (ip%Ap ** e + ip%At ** e)
    ip%Rc = ip%r0c * (ip%Ap ** e + ip%At ** e)
    ip%Rp = ip%r0p * ip%Ap ** e
    ip%Rt = ip%r0t * ip%At ** e
    ip%rmass = ip%Ap * ip%At * mass / (ip%Ap + ip%At)
    ip%nth = int((ip%thmax - ip%thmin) / ip%dth) + 1
    ip%rgrid = nint((ip%rmax - ip%rmin) / ip%dr)
! cutting parameter for wn (in RMT)
    if (ip%rcut > ip%rmax) ip%rcut = ip%rmax
    ip%rgrid_cut = min(nint((ip%rcut - ip%rmin) / ip%dr),ip%rgrid)
    ip%rmax =  ip%rmin + dble(ip%rgrid) * ip%dr   ! necessary modification
    if (ip%dE /= 0.0d0) then
      ip%Egrid = nint((ip%Emax - ip%Emin) / ip%dE)
      ip%di = nint(ip%del_E / ip%dE)
    else
      ip%Egrid = 0
      ip%di = 1
    end if
    
    if (ip%coup_t == - 1) then
      if (ip%coup_p == 0) then
        ip%Rm = ip%Rp
        ip%Nphonon = ip%Nphonon_p
        ip%lambda  = ip%lambda_p
        ip%omega   = ip%omega_p
        ip%betaL   = ip%betaL_p
        ip%betaLC  = ip%betaL_pC
        ip%coup    = ip%coup_p
      else if (ip%coup_p == 1) then
        ip%Rm = ip%Rp
        ip%E2    = ip%E2p
        ip%nrot  = ip%nrot_p
        ip%beta2 = ip%beta2p
        ip%beta4 = ip%beta4p
        ip%coup  = ip%coup_p
      end if
    else if (ip%coup_p == -1) then
      if (ip%coup_t == 0) then
        ip%Rm = ip%Rt
        ip%Nphonon = ip%Nphonon_t
        ip%lambda  = ip%lambda_t
        ip%omega   = ip%omega_t
        ip%betaL   = ip%betaL_t
        ip%betaLC  = ip%betaL_tC
        ip%coup    = ip%coup_t
      else if (ip%coup_t == 1) then
        ip%Rm = ip%Rt
        ip%E2    = ip%E2t
        ip%nrot  = ip%nrot_t
        ip%beta2 = ip%beta2t
        ip%beta4 = ip%beta4t
        ip%coup  = ip%coup_t
      end if
    end if

  end subroutine
!----------------------------------------------------------------------!
  function count_lines(Fname) result(N)
  implicit none
  character(len=*), intent(in) :: Fname
  integer :: N, ios
  integer, parameter :: Nfile = 9

  open(Nfile,file=Fname,status="old",action="read",iostat=ios)
    N = 0
    do
      read(Nfile,*,iostat=ios)
      if (ios < 0) then
        exit
      else
        N = N + 1
      end if
    end do
!   write(6,*) "lines =", N
  close(Nfile)

  end function
end module


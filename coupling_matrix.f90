module coupling_matrix
  use input_data, only : inp
  use potentials
  use angular_coup
  private

  type, public, extends(potential) :: coup_mat
    real(8), allocatable, dimension(:) :: e_n
    real(8), private, allocatable, dimension(:,:) :: Vcp_linear
    contains
    procedure :: coup_mat_
    procedure :: make_Vcp
    procedure :: destruct_coup_mat
    procedure :: get_Vcp
    procedure, private :: Vn_rot
    procedure, private :: Vc_rot
    procedure, private :: Vn_vib
    procedure, private :: Vc_vib
    procedure, private :: Vcoup
    procedure, private :: Vcoup_pro_tar
    procedure, private :: eps
    procedure, private :: eps_pro_tar
    procedure, private :: Vn_pro_tar
    procedure, private :: Vc_pro_tar
  end type

  contains
!----------------------------------------------------------------------
  subroutine coup_mat_(this, ip)
    implicit none
    class(coup_mat), intent(out) :: this
    type(inp), target, intent(in) :: ip
    integer :: n, i

    call this%potential%potential_(ip)
    n = ip%Nch * (ip%Nch + 1) / 2
    allocate(this%Vcp_linear(n,ip%rgrid+2))
    allocate(this%e_n(ip%Nch))

! excitation energy
    if (this%ip%coup_t == - 1 .and. this%ip%coup_p == - 1) then
      this%e_n = 0.0d0
    else if (this%ip%coup_t == - 1 .or. this%ip%coup_p == - 1) then
      do i=1, this%ip%Nch
        this%e_n(i) = this%eps(i)
      end do
    else
      call this%eps_pro_tar(this%e_n)
    end if

  end subroutine
!----------------------------------------------------------------------
  elemental subroutine destruct_coup_mat(this)
    implicit none
    class(coup_mat), intent(inout) :: this

    if (allocated(this%Vcp_linear)) deallocate(this%Vcp_linear)
  end subroutine
!----------------------------------------------------------------------
  function get_Vcp(this, ir) result(V)
    implicit none
    class(coup_mat), intent(in) :: this
    integer, intent(in) :: ir
    integer :: n, m, nst, nlen
    real(8) :: V(this%ip%Nch,this%ip%Nch)

    nst = 1
    do n=1, this%ip%Nch
      nlen = this%ip%Nch - n + 1
      V(n:this%ip%Nch,n) = this%Vcp_linear(nst:nst-1+nlen,ir)
      V(n,2:this%ip%Nch) = V(2:this%ip%Nch,n)
      nst = nst + nlen
    end do

    return
  end function
!----------------------------------------------------------------------
  subroutine make_Vcp(this)
    implicit none
    class(coup_mat), intent(inout) :: this
    integer :: i, n, m
    real(8) :: r
    real(8), allocatable :: Vnm(:,:)
    character(len=40), parameter :: FM='(1x,10f8.3)'

! coupling matrix element
    do i=1, this%ip%rgrid+2
      r = this%ip%rmin + dble(i-1) *  this%ip%dr
      if (this%ip%coup_t == - 1 .and. this%ip%coup_p == - 1) then
        this%Vcp_linear(:,i) = 0.0d0
      else if (this%ip%coup_t == - 1 .or. this%ip%coup_p == - 1) then
        call this%Vcoup(r, i)
      else
        call this%Vcoup_pro_tar(r, i)
      end if
    end do

    return
  end subroutine
!----------------------------------------------------------------------!
  subroutine Vcoup(this, r, ir)
    implicit none
    class(coup_mat), intent(inout) :: this
    integer :: n, m, nst, nlen
    integer, intent(in) :: ir
    real(8), intent(in) :: r
    real(8), allocatable, dimension(:,:) :: Vn_cp, Vc_cp
    real(8) :: V

    allocate(Vn_cp(this%ip%Nch,this%ip%Nch))
    allocate(Vc_cp(this%ip%Nch,this%ip%Nch))


    if (this%ip%coup == 0) then
      call this%Vn_vib(ir, r, Vn_cp)
      call this%Vc_vib(r, Vc_cp)
    else if (this%ip%coup == 1) then
      call this%Vn_rot(r, Vn_cp)
      call this%Vc_rot(r, Vc_cp)
    end if

    V = this%Vn(r)
    nst = 1
    do n=1, this%ip%Nch
      nlen = this%ip%Nch - n + 1
      this%Vcp_linear(nst,ir) = (Vn_cp(n,n) - V) + Vc_cp(n,n)
      this%Vcp_linear(nst+1:nst-1+nlen,ir) = Vn_cp(n+1:this%ip%Nch,n) + Vc_cp(n+1:this%ip%Nch,n)
      nst = nst + nlen
    end do

    deallocate(Vn_cp, Vc_cp)
  end subroutine
!----------------------------------------------------------------------
  subroutine Vn_rot(this, r, Vpot)
!   use mkl95_lapack
    use eigen
    implicit none
    class(coup_mat), intent(in) :: this
    integer :: i, n, m, I1, I2
    real(8), intent(in) :: r
    real(8), intent(out) :: Vpot(this%ip%Nch,this%ip%Nch)
    real(8), parameter :: PI = 3.141592653589793d0
    real(8), parameter :: SQRT4PI_INV = 1.0d0/sqrt(4.0d0*PI)
    real(8), parameter :: CONST1 = sqrt(5.0d0/(4.0d0*PI))
    real(8), parameter :: CONST2 = sqrt(9.0d0/(4.0d0*PI))
    real(8), allocatable, save :: O(:,:), Oa(:)
    real(8), dimension(this%ip%Nch) :: Vn
    real(8) :: wg1, wg2, w1, w2, c, f

    if (.not. allocated(O)) then
      allocate(O(this%ip%Nch,this%ip%Nch), Oa(this%ip%Nch))
      w1 = CONST1 * this%ip%beta2 * this%ip%Rm
      w2 = CONST2 * this%ip%beta4 * this%ip%Rm

      if (this%ip%coup_p == 1) then
        do m=1, this%ip%Nrot_p+1
          I2 = 2 * (m - 1)
          do n=m, this%ip%Nrot_p+1
            I1 = 2 * (n - 1)
            call wig3j(I1, 2, I2, 0, 0, 0, wg1)
            call wig3j(I1, 4, I2, 0, 0, 0, wg2)
            O(n, m) = (w1 * wg1 * wg1 + w2 * wg2 * wg2)  &
                    * sqrt(dble((2 * I1 + 1) * (2 * I2 + 1)))
          end do
        end do

!  2nd mode in projectile
        c = this%ip%betaL_p2 * this%ip%Rm * SQRT4PI_INV 
        do m=1, this%ip%Nch
          do n=this%ip%Nrot_p+2, this%ip%Nch
            if (m == 1 .and. n == this%ip%Nrot_p+2) then
              O(n,m) = c
            else if (n == m + 1 .and. m > this%ip%Nrot_p+1) then
              O(n,m) = c * sqrt(dble(m-this%ip%Nrot_p))
            else
              O(n,m) = 0.0d0
            end if
          end do
        end do

      else if (this%ip%coup_t == 1) then
        do m=1, this%ip%nrot_t+1
          I2 = 2 * (m - 1)
          do n=m, this%ip%nrot_t+1
            I1 = 2 * (n - 1)
            call wig3j(I1, 2, I2, 0, 0, 0, wg1)
            call wig3j(I1, 4, I2, 0, 0, 0, wg2)
            O(n, m) = (w1 * wg1 * wg1 + w2 * wg2 * wg2)  &
                    * sqrt(dble((2 * I1 + 1) * (2 * I2 + 1)))
          end do
        end do

!  2nd mode in target
        c = this%ip%betaL_t2 * this%ip%Rm * SQRT4PI_INV 
        do m=1, this%ip%Nch
          do n=this%ip%nrot_t+2, this%ip%Nch
            if (m == 1 .and. n == this%ip%nrot_t+2) then
              O(n,m) = c
            else if (n == m + 1 .and. m > this%ip%nrot_t+1) then
              O(n,m) = c * sqrt(dble(m-this%ip%nrot_t))
            else
              O(n,m) = 0.0d0
            end if
          end do
        end do
      end if
!     call syev(O, Oa, 'V','L')
      call houshld_ql(O, Oa, 'V')
    end if

!   Vpot = 0.0d0
!   f = exp((r - this%ip%Rn)/this%ip%a)
!   Vpot = -O * this%ip%V0/this%ip%a * f / (1.0d0+f)**2
!   forall (n=1:this%ip%Nch)
!     Vpot(n,n) = -this%ip%V0 / (1.0d0 + exp((r - this%ip%Rn)/this%ip%a))
!   end forall


    Vn = - this%ip%V0 / (1.0d0 + exp((r - this%ip%Rn - Oa) / this%ip%a))

    Vpot = 0.0d0
    do m=1, this%ip%Nch
      do n=m, this%ip%Nch
        do i=1, this%ip%Nch
          Vpot(n,m) = Vpot(n,m) + O(n,i) * O(m,i) * Vn(i)
        end do
      end do
    end do


    return
  end subroutine
!----------------------------------------------------------------------
  subroutine Vc_rot(this, r, V)
    use global_constant, only : e2, PI
    implicit none
    class(coup_mat), intent(in) :: this
    integer :: n, m, I1, I2
    real(8), intent(in) :: r
    real(8), intent(out) :: V(this%ip%Nch,this%ip%Nch)
    real(8), parameter :: CONST1 = sqrt(5.0d0/(4.0d0*PI))
    real(8), parameter :: CONST2 = 2.0d0/7.0d0*sqrt(5.0d0/PI)
    real(8), parameter :: CONST3 = 9.0d0/(7.0d0*sqrt(PI))
    real(8), parameter :: SQRT4PI_INV  = 1.0d0/sqrt(4.0d0*PI)
    real(8) :: x, x2, wg1, wg2, w1, w2, ZZe2

    ZZe2 = this%ip%Zp * this%ip%Zt * e2
    w1 = 0.6d0 * ZZe2 * CONST1 * this%ip%beta2 * (1.0d0 + CONST2 * this%ip%beta2)
    w2 = ZZe2 * SQRT4PI_INV * (this%ip%beta4 + CONST3 * this%ip%beta2 * this%ip%beta2)
    if (r > this%ip%Rm) then
      x = this%ip%Rm / r
      x2 = x * x / r
    else
      x = r / this%ip%Rm
      x2 = x * x / this%ip%Rm
    end if

    if (this%ip%coup_p == 1) then
      do m=1, this%ip%Nrot_p+1
        I2 = 2 * (m - 1)
        do n=m, this%ip%Nrot_p+1
          I1 = 2 * (n - 1)
          call wig3j(I1, 2, I2, 0, 0, 0, wg1)
          call wig3j(I1, 4, I2, 0, 0, 0, wg2)
          V(n, m) = (w1 * wg1 * wg1 + w2 * (wg2 * x) ** 2) * x2   &
                        * sqrt(dble((2*I1+1) * (2*I2+1)))
        end do
      end do

!   The second mode of phonon excitation
      w1 = this%ip%betaL_p2C*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_p2+1) * ZZe2
      if (r > this%ip%Rm) then
        x = w1 * (this%ip%Rm / r) ** this%ip%lambda_p2 / r
      else
        x = w1 * (r / this%ip%Rm) ** this%ip%lambda_p2 / this%ip%Rm
      end if

      do m=1, this%ip%Nch
        do n=this%ip%nrot_p+2, this%ip%Nch
          if (m == 1 .and. n == this%ip%nrot_p+2) then
            V(n,m) = x
          else if (n == m + 1 .and. m > this%ip%nrot_p+1) then
            V(n,m) = x * sqrt(dble(m-this%ip%nrot_p))
          else
            V(n,m) = 0.0d0
          end if
        end do
      end do

    else if (this%ip%coup_t == 1) then
      do m=1, this%ip%nrot_t+1
        I2 = 2 * (m - 1)
        do n=m, this%ip%nrot_t+1
          I1 = 2 * (n - 1)
          call wig3j(I1, 2, I2, 0, 0, 0, wg1)
          call wig3j(I1, 4, I2, 0, 0, 0, wg2)
          V(n, m) = (w1 * wg1 * wg1 + w2 * (wg2 * x) ** 2) * x2   &
                        * sqrt(dble((2*I1+1) * (2*I2+1)))
        end do
      end do

!   The second mode of phonon excitation
      w1 = this%ip%betaL_t2C*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_t2+1) * ZZe2
      if (r > this%ip%Rm) then
        x = w1 * (this%ip%Rm / r) ** this%ip%lambda_t2 / r
      else
        x = w1 * (r / this%ip%Rm) ** this%ip%lambda_t2 / this%ip%Rm
      end if

      do m=1, this%ip%Nch
        do n=this%ip%nrot_t+2, this%ip%Nch
          if (m == 1 .and. n == this%ip%nrot_t+2) then
            V(n,m) = x
          else if (n == m + 1 .and. m > this%ip%nrot_t+1) then
            V(n,m) = x * sqrt(dble(m-this%ip%nrot_t))
          else
            V(n,m) = 0.0d0
          end if
        end do
      end do
    end if

    return
  end subroutine
!----------------------------------------------------------------------
  subroutine  Vn_vib(this, ir, r, Vpot)
!   use mkl95_lapack
    use eigen
    implicit none
    class(coup_mat), intent(in) :: this
    integer, intent(in) :: ir
    integer :: n, m, i, n2, idum
    real(8), intent(in) :: r
    real(8), intent(out) :: Vpot(:,:)
    real(8), allocatable, save :: O(:,:), Oa(:)
    real(8), allocatable, save :: O_sub(:,:), Oa_sub(:)
    real(8), dimension(this%ip%Nch) :: Vn
    real(8), parameter :: PI = 3.141592653589793d0
    real(8), parameter :: SQRT4PI_INV = 1.0d0/sqrt(4.0d0*PI)
    real(8) :: c

    if (.not. allocated(O)) then

      c = this%ip%betaL * this%ip%Rm * SQRT4PI_INV 

      allocate(O(this%ip%Nch,this%ip%Nch))
      allocate(Oa(this%ip%Nch))

      if (this%ip%coup_p == 0) then
        do m=1, this%ip%Nphonon_p+1
          do n=m, this%ip%Nphonon_p+1
            if (n == m + 1) then
              O(n,m) = c * sqrt(dble(m))
            else
              O(n,m) = 0.0d0
            end if
          end do
        end do

!  2nd mode in projectile
        c = this%ip%betaL_p2 * this%ip%Rm * SQRT4PI_INV 
        do m=1, this%ip%Nch_p
          do n=max(m,this%ip%Nphonon_p+2), this%ip%Nch_p
            if (m == 1 .and. n == this%ip%Nphonon_p+2) then
              O(n,m) = c
            else if (n == m + 1 .and. m > this%ip%Nphonon_p+1) then
              O(n,m) = c * sqrt(dble(m-this%ip%Nphonon_p))
            else
              O(n,m) = 0.0d0
            end if
          end do
        end do

      else if (this%ip%coup_t == 0) then
        do m=1, this%ip%Nphonon_t+1
          do n=m, this%ip%Nphonon_t+1
            if (n == m + 1) then
              O(n,m) = c * sqrt(dble(m))
            else
              O(n,m) = 0.0d0
            end if
          end do
        end do

! 2nd mode in target
        c = this%ip%betaL_t2 * this%ip%Rm * SQRT4PI_INV 
        do m=1, this%ip%Nch_t
          do n=this%ip%Nphonon_t+2, this%ip%Nch_t
            if (m == 1 .and. n == this%ip%Nphonon_t+2) then
              O(n,m) = c
            else if (n == m + 1 .and. m > this%ip%Nphonon_t+1) then
              O(n,m) = c * sqrt(dble(m-this%ip%Nphonon_t))
            else
              O(n,m) = 0.0d0
            end if
          end do
        end do
      end if

!     call syev(O, Oa, 'V','L')
      call houshld_ql(O, Oa, 'V')

    end if
    

    Vn = - this%ip%V0 / (1.0d0 + exp((r - this%ip%Rn - Oa) / this%ip%a))

    Vpot = 0.0d0
    do m=1, this%ip%Nch
      do n=m, this%ip%Nch
        do i=1, this%ip%Nch
          Vpot(n,m) = Vpot(n,m) + O(n,i) * O(m,i) * Vn(i)
        end do
      end do
    end do

    return
  end subroutine
!----------------------------------------------------------------------
  subroutine Vc_vib(this, r, V)
    use global_constant, only : e2, PI
    implicit none
    class(coup_mat), intent(in) :: this
    integer :: n, m
    real(8), intent(in) :: r
    real(8), intent(out) :: V(:,:)
    real(8), parameter :: CONST1 = sqrt(5.0d0/(4.0d0*PI))
    real(8), parameter :: CONST2 = 2.0d0/7.0d0*sqrt(5.0d0/PI)
    real(8), parameter :: CONST3 = 9.0d0/(7.0d0*sqrt(PI))
    real(8), parameter :: SQRT4PI_INV = 1.0d0/sqrt(4.0d0*PI)
    real(8) :: x, w, ZZe2

    ZZe2 = this%ip%Zp * this%ip%Zt * e2
    w = this%ip%betaLC*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda+1) * ZZe2
    if (r > this%ip%Rm) then
      x = w * (this%ip%Rm / r) ** this%ip%lambda / r
    else
      x = w * (r / this%ip%Rm) ** this%ip%lambda / this%ip%Rm
    end if

    if (this%ip%coup_p == 0) then
      do m=1, this%ip%Nphonon_p+1
        do n=m, this%ip%Nphonon_p+1
          if (n == m + 1) then
            V(n,m) = x * sqrt(dble(m))
          else
            V(n,m) = 0.0d0
          end if
        end do
      end do

! 2nd mode in projectile
      w = this%ip%betaL_p2C*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_p2+1) * ZZe2
      if (r > this%ip%Rm) then
        x = w * (this%ip%Rm / r) ** this%ip%lambda_p2 / r
      else
        x = w * (r / this%ip%Rm) ** this%ip%lambda_p2 / this%ip%Rm
      end if

      do m=1, this%ip%Nch_p
        do n=max(m,this%ip%Nphonon_p+2), this%ip%Nch_p
          if (m == 1 .and. n == this%ip%Nphonon_p+2) then
            V(n,m) = x
          else if (n == m + 1 .and. m > this%ip%Nphonon_p+1) then
            V(n,m) = x * sqrt(dble(m-this%ip%Nphonon_p))
          else
            V(n,m) = 0.0d0
          end if
        end do
      end do

    else if (this%ip%coup_t == 0) then
      do m=1, this%ip%Nphonon_t+1
        do n=m, this%ip%Nphonon_t+1
          if (n == m + 1) then
            V(n,m) = x * sqrt(dble(m))
          else
            V(n,m) = 0.0d0
          end if
        end do
      end do

! 2nd mode in target
      w = this%ip%betaL_t2C*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_t2+1) * ZZe2
      if (r > this%ip%Rm) then
        x = w * (this%ip%Rm / r) ** this%ip%lambda_t2 / r
      else
        x = w * (r / this%ip%Rm) ** this%ip%lambda_t2 / this%ip%Rm
      end if

      do m=1, this%ip%Nch_t
        do n=this%ip%Nphonon_t+2, this%ip%Nch_t
          if (m == 1 .and. n == this%ip%Nphonon_t+2) then
            V(n,m) = x
          else if (n == m + 1 .and. m > this%ip%Nphonon_t+1) then
            V(n,m) = x * sqrt(dble(m-this%ip%Nphonon_t))
          else
            V(n,m) = 0.0d0
          end if
        end do
      end do

    end if

  end subroutine
!----------------------------------------------------------------------
  function eps(this, n) result(e)
    implicit none
    class(coup_mat), intent(in) :: this
    integer, intent(in) :: n
    integer :: I
    real(8) :: e

    if (this%ip%coup_p == 0) then
      if (n < this%ip%Nphonon_p+2) then
        e = dble(n-1) * this%ip%omega
      else
        e = dble(n-this%ip%Nphonon_p-1) * this%ip%omega_p2
      end if
    else if (this%ip%coup_t == 0) then
      if (n < this%ip%Nphonon_t+2) then
        e = dble(n-1) * this%ip%omega
      else
        e = dble(n-this%ip%Nphonon_t-1) * this%ip%omega_t2
      end if
    else if (this%ip%coup_p == 1) then
      if (n < this%ip%nrot_p+2) then 
        I = 2 * (n - 1)
        e = dble(I * (I + 1)) / 6.0d0 * this%ip%E2
      else
        e = dble(n - this%ip%nrot_p - 2) * this%ip%omega_p2
      end if
    else if (this%ip%coup_t == 1) then
      if (n < this%ip%nrot_t+2) then 
        I = 2 * (n - 1)
        e = dble(I * (I + 1)) / 6.0d0 * this%ip%E2
      else
        e = dble(n - this%ip%nrot_t - 2) * this%ip%omega_t2
      end if
    end if

  end function
!----------------------------------------------------------------------!
  subroutine Vcoup_pro_tar(this, r, ir)
    use potentials, only : Vn
    implicit none
    class(coup_mat), intent(inout) :: this
    integer, intent(in) :: ir
    integer :: n, m, nst, nlen
    real(8), intent(in) :: r
    real(8), allocatable, dimension(:,:) :: Vn_cp, Vc_cp
    real(8) :: V

    allocate(Vn_cp(this%ip%Nch,this%ip%Nch), Vc_cp(this%ip%Nch,this%ip%Nch))
    call this%Vn_pro_tar(ir, r, Vn_cp)
    call this%Vc_pro_tar(r, Vc_cp)
    V = this%Vn(r)
    nst = 1
    do n=1, this%ip%Nch
      nlen = this%ip%Nch - n + 1
      this%Vcp_linear(nst,ir) = (Vn_cp(n,n) - V)+ Vc_cp(n,n)
      this%Vcp_linear(nst+1:nst-1+nlen,ir) = Vn_cp(n+1:this%ip%Nch,n) + Vc_cp(n+1:this%ip%Nch,n)
      nst = nst + nlen
    end do

    deallocate(Vn_cp, Vc_cp)

  end subroutine
!----------------------------------------------------------------------
  subroutine Vn_pro_tar(this, ir, r, Vpot)
!   use mkl95_lapack
    use eigen
    use global_constant, only : PI
    implicit none
    class(coup_mat), intent(inout) :: this
    integer :: n, m, i, I1, I2, n1st, idum, st, ed, m_n1st
    integer, intent(in) :: ir
    real(8), intent(in) :: r
    real(8), allocatable, save :: O(:,:), Oa(:)
    real(8), intent(out) :: Vpot(this%ip%Nch,this%ip%Nch)
    real(8), dimension(:,:), allocatable :: Q
    real(8), dimension(this%ip%Nch) :: Vn
    real(8), parameter :: SQRT4PI_INV = 1.0d0/sqrt(4.0d0*PI)
    real(8), parameter :: CONST1 = sqrt(5.0d0/(4.0d0*PI))
    real(8), parameter :: CONST2 = sqrt(9.0d0/(4.0d0*PI))
    real(8) :: ct, wg1, wg2, w1, w2, Qnm

    if (r > this%ip%rcut) then
      Vpot = 0.0d0
      return
    end if

    if (.not. allocated(O)) then
      allocate(O(this%ip%Nch,this%ip%Nch))
      allocate(Oa(this%ip%Nch))
      allocate(Q(this%ip%Nch_t,this%ip%Nch_t))

      do n=1, this%ip%Nch_t
        Q(n:this%ip%Nch_t,n) = 0.0d0
      end do

! 1st mode
      select case(this%ip%coup_t)
        case(0)
          ct = this%ip%betaL_t * this%ip%Rt * SQRT4PI_INV
          do m=1, this%ip%Nphonon_t+1
            do n=m, this%ip%Nphonon_t+1
              if (n == m + 1) then
                Q(n,m) = ct * sqrt(dble(m))
              else
                Q(n,m) = 0.0d0
              end if
            end do
          end do
          n1st = this%ip%Nphonon_t
        case(1)
          w1 = CONST1 * this%ip%beta2t * this%ip%Rt
          w2 = CONST2 * this%ip%beta4t * this%ip%Rt
          do m=1, this%ip%Nrot_t+1
            I2 = 2 * (m - 1)
            do n=m, this%ip%Nrot_t+1
              I1 = 2 * (n - 1)
              call wig3j(I1, 2, I2, 0, 0, 0, wg1)
              call wig3j(I1, 4, I2, 0, 0, 0, wg2)
              Q(n, m) = (w1 * wg1 * wg1 + w2 * wg2 * wg2)  &
                          * sqrt(dble((2 * I1 + 1) * (2 * I2 + 1)))
            end do
          end do
          n1st = this%ip%Nrot_t
      end select

! 2nd mode
      ct = this%ip%betaL_t2 * this%ip%Rt * SQRT4PI_INV 
      do m=1, this%ip%Nch_t
        do n=max(m,n1st+2), this%ip%Nch_t
          if (m == 1 .and. n == n1st+2) then
            Q(n,m) = ct
          else if (n == m + 1 .and. m > n1st+1) then
            Q(n,m) = ct * sqrt(dble(max(0,m-n1st)))
          else
            Q(n,m) = 0.0d0
          end if
        end do
      end do

! projectile excitations
      call sub_Vn(Q, O)

!     call syev(O, Oa, 'V','L')
      call houshld_ql(O, Oa, 'V')
      deallocate(Q)
    end if

    Vn = - this%ip%V0 / (1.0d0 + exp((r - this%ip%Rn - Oa) / this%ip%a))

    Vpot = 0.0d0
    do m=1, this%ip%Nch
      do n=m, this%ip%Nch
        do i=1, this%ip%Nch
          Vpot(n,m) = Vpot(n,m) + O(n,i) * O(m,i) * Vn(i)
        end do
      end do
    end do

    contains
!----------------------------------------------------------------------
    subroutine sub_Vn(Q, O)
      implicit none
      integer :: i, I1, I2, n, m, nt, mt, n1st
      real(8), intent(in) :: Q(:,:)
      real(8), intent(out) :: O(:,:)
      real(8) :: cp, d, w1, w2, wg1, wg2, Nch_t2

      select case(this%ip%coup_p)
! 1st mode in projectile
        case(0)
          cp = this%ip%betaL_p * this%ip%Rp * SQRT4PI_INV
          do m=1, this%ip%Nphonon_p+1
            mt = (m - 1) * this%ip%Nch_t
            do n=m, this%ip%Nphonon_p+1
              nt = (n - 1) * this%ip%Nch_t
              if (n == m + 1) then
                O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
                do i=1, this%ip%Nch_t
                  O(nt+i,mt+i) = cp * sqrt(dble(m))
                end do
              else if (n == m) then
                O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = Q
              else
                O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              end if
            end do
          end do
          n1st = this%ip%Nphonon_p
        case(1)
          w1 = CONST1 * this%ip%beta2p * this%ip%Rp
          w2 = CONST2 * this%ip%beta4p * this%ip%Rp
          do m=1, this%ip%Nrot_p+1
            I2 = 2 * (m - 1)
            mt = (m - 1) * this%ip%Nch_t
            do n=m, this%ip%Nrot_p+1
              nt = (n - 1) * this%ip%Nch_t
              I1 = 2 * (n - 1)
              call wig3j(I1, 2, I2, 0, 0, 0, wg1)
              call wig3j(I1, 4, I2, 0, 0, 0, wg2)
              if (n == m) then
                O(nt+1:nt+this%ip%Nch_t, mt+1:mt+this%ip%Nch_t) = Q
              else
                O(nt+1:nt+this%ip%Nch_t, mt+1:mt+this%ip%Nch_t) = 0.0d0
              end if
              d = (w1 * wg1 * wg1 + w2 * wg2 * wg2)           &
                    * sqrt(dble((2 * I1 + 1) * (2 * I2 + 1)))
              do i=1, this%ip%Nch_t
                O(nt+i,mt+i) = O(nt+i,mt+i) + d
              end do
            end do
          end do
          n1st = this%ip%Nrot_p
      end select

! 2nd mode in projectile
      if (this%ip%Nphonon_p2 > 0) then
        cp = this%ip%betaL_p2 * this%ip%Rp * SQRT4PI_INV
        do m=1, this%ip%Nch_p
          mt = (m - 1) * this%ip%Nch_t
          do n=max(m,n1st+2), this%ip%Nch_p
            nt = (n - 1) * this%ip%Nch_t
            if (m == 1 .and. n == n1st+2) then
              O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              do i=1, this%ip%Nch_t
                O(nt+i,mt+i) = cp
              end do
            else if (n == m + 1 .and. m > n1st+1) then
              O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              do i=1, this%ip%Nch_t
                O(nt+i,mt+i) = cp * sqrt(dble(m-n1st))
              end do
            else if (n == m) then
              O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = Q
            else
              O(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
            end if
          end do
        end do
      end if

      return
    end subroutine
  end subroutine
!----------------------------------------------------------------------!
  subroutine Vc_pro_tar(this, r, V)
    use global_constant, only : e2, PI
    implicit none
    class(coup_mat), intent(in) :: this
    integer :: n, m, I1, I2, n1st
    real(8), intent(in) :: r
    real(8), intent(out) :: V(:,:)
    real(8), dimension(this%ip%Nch_t,this%ip%Nch_t) :: Vt
    real(8), parameter :: CONST1 = sqrt(5.0d0/(4.0d0*PI))
    real(8), parameter :: CONST2 = 2.0d0/7.0d0*sqrt(5.0d0/PI)
    real(8), parameter :: CONST3 = 9.0d0/(7.0d0*sqrt(PI))
    real(8), parameter :: SQRT4PI_INV = 1.0d0/sqrt(4.0d0*PI)
    real(8) :: x, x2, w, wg1, wg2, w1, w2, ZZe2

    ZZe2 = this%ip%Zt * this%ip%Zp * e2

! 1st-mode
    select case (this%ip%coup_t)
      case(0)
        w = this%ip%betaL_tC*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_t+1) * ZZe2
        if (r > this%ip%Rt) then
          x = w * (this%ip%Rt / r) ** this%ip%lambda_t / r
        else
          x = w * (r / this%ip%Rt) ** this%ip%lambda_t / this%ip%Rt
        end if

        do m=1, this%ip%Nphonon_t+1
          do n=m, this%ip%Nphonon_t+1
            if (n == m + 1) then
              Vt(n,m) = x * sqrt(dble(m))
            else
              Vt(n,m) = 0.0d0
            end if
          end do
        end do
        n1st = this%ip%Nphonon_t
      case(1)
        if (r > this%ip%Rt) then
          x = this%ip%Rt / r
          x2 = x * x / r
        else
          x = r / this%ip%Rt
          x2 = x * x / this%ip%Rt
        end if
        w1 = 0.6d0 * ZZe2 * CONST1* this%ip%beta2t * (1.0d0 + CONST2*this%ip%beta2t)
        w2 = ZZe2 * SQRT4PI_INV * (this%ip%beta4t + CONST3*this%ip%beta2t*this%ip%beta2t)
        do m=1, this%ip%Nrot_t+1
          I2 = 2 * (m - 1)
          do n=m, this%ip%Nrot_t+1
            I1 = 2 * (n - 1)
            call wig3j(I1, 2, I2, 0, 0, 0, wg1)
            call wig3j(I1, 4, I2, 0, 0, 0, wg2)
            Vt(n, m) = (w1 * wg1 * wg1 + w2 * (wg2 * x) ** 2) * x2   &
                           * sqrt(dble((2*I1+1) * (2*I2+1)))
          end do
        end do
        n1st = this%ip%Nrot_t
    end select

! 2nd-mode in target
    w = this%ip%betaL_t2C * SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_t2+1) * ZZe2
    if (r > this%ip%Rt) then
      x = w * (this%ip%Rt / r) ** this%ip%lambda_t2 / r
    else
      x = w * (r / this%ip%Rt) ** this%ip%lambda_t2 / this%ip%Rt
    end if

    do m=1, this%ip%Nch_t
      do n=max(m,n1st+2), this%ip%Nch_t
        if (m == 1 .and. n == n1st+2) then
          Vt(n,m) = x
        else if (n == m + 1 .and. m > n1st+1) then
          Vt(n,m) = x * sqrt(dble(max(0,m-n1st)))
        else
          Vt(n,m) = 0.0d0
        end if
      end do
    end do

    call sub_Vc(r, Vt, V)

    contains
!----------------------------------------------------------------------!
    subroutine sub_Vc(r, Vt, V)
      use global_constant, only : e2
      implicit none
      integer :: i, I1, I2, m, n, nt, mt, n1st
      real(8), intent(in) :: r, Vt(:,:)
      real(8), intent(out) :: V(:,:)
      real(8) :: x, x2, w, w1, w2, wg1, wg2, d, ZZe2

      ZZe2 = this%ip%Zt * this%ip%Zp * e2
      select case(this%ip%coup_p)
        case(0)
          w = this%ip%betaL_pC*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_p+1) * ZZe2
          if (r > this%ip%Rp) then
            x = w * (this%ip%Rp / r) ** this%ip%lambda_p / r
          else
            x = w * (r / this%ip%Rp) ** this%ip%lambda_p / this%ip%Rp
          end if
          do m=1, this%ip%Nphonon_p+1
            mt = (m - 1) * this%ip%Nch_t
            do n=m, this%ip%Nphonon_p+1
              nt = (n - 1) * this%ip%Nch_t
              if (n == m + 1) then
                V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
                do i=1, this%ip%Nch_t
                  V(nt+i, mt+i) = x * sqrt(dble(m))
                end do
              else if (n == m) then
                V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = Vt
              else
                V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              end if
            end do
          end do
          n1st = this%ip%Nphonon_p
        case(1)
          if (r > this%ip%Rp) then
            x = this%ip%Rp / r
            x2 = x * x / r
          else
            x = r / this%ip%Rp
            x2 = x * x / this%ip%Rp
          end if
          w1 = 0.6d0 * ZZe2 * CONST1 * this%ip%beta2p *(1.0d0 + CONST2*this%ip%beta2p)
          w2 = ZZe2 * SQRT4PI_INV * (this%ip%beta4p + CONST3*this%ip%beta2p*this%ip%beta2p)
          do m=1, this%ip%Nrot_p+1
            mt = (m - 1) * this%ip%Nch_t
            I2 = 2 * (m - 1)
            do n=m, this%ip%Nrot_p+1
              nt = (n - 1) * this%ip%Nch_t
              I1 = 2 * (n - 1)
              call wig3j(I1, 2, I2, 0, 0, 0, wg1)
              call wig3j(I1, 4, I2, 0, 0, 0, wg2)
              if (n == m) then
                V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = Vt
              else
                V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              end if
              d = (w1 * wg1 * wg1 + w2 * (wg2*x) ** 2) * x2       &
                     * sqrt(dble((2*I1+1) * (2*I2+1)))
              do i=1, this%ip%Nch_t
                V(nt+i,mt+i) = V(nt+i,mt+i) + d
              end do
            end do
          end do
          n1st = this%ip%Nrot_p
      end select

! 2nd mode in projectile
      if (this%ip%Nphonon_p2 > 0) then
        w = this%ip%betaL_p2C*SQRT4PI_INV * 3.0d0/dble(2*this%ip%lambda_p2+1) * ZZe2
        if (r > this%ip%Rp) then
          x = w * (this%ip%Rp / r) ** this%ip%lambda_p2 / r
        else
          x = w * (r / this%ip%Rp) ** this%ip%lambda_p2 / this%ip%Rp
        end if
        do m=1, this%ip%Nch_p
          mt = (m - 1) * this%ip%Nch_t
          do n=max(m,n1st+2), this%ip%Nch_p
            nt = (n - 1) * this%ip%Nch_t
            if (m == 1 .and. n == n1st+2) then
              V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              do i=1, this%ip%Nch_t
                V(nt+i, mt+i) = x
              end do
            else if (n == m + 1 .and. m > n1st+1) then
              V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
              do i=1, this%ip%Nch_t
                V(nt+i, mt+i) = x * sqrt(dble(m-n1st))
              end do
            else if (n == m) then
              V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = Vt
            else
              V(nt+1:nt+this%ip%Nch_t,mt+1:mt+this%ip%Nch_t) = 0.0d0
            end if
          end do
        end do
      end if

    end subroutine
  end subroutine
!----------------------------------------------------------------------!
  subroutine eps_pro_tar(this, eps)
    implicit none
    class(coup_mat), intent(in) :: this
    integer :: n, m, I, n1st
    real(8), intent(out) :: eps(this%ip%Nch)
    real(8), dimension(this%ip%Nch_t) :: eps_t

    select case(this%ip%coup_t)
      case(0)
        forall (n=1:this%ip%Nphonon_t+1) eps_t(n) = dble(n-1)*this%ip%omega_t
        n1st = this%ip%Nphonon_t
      case(1)
        do n=1, this%ip%Nrot_t + 1
          I = 2 * (n - 1)
          eps_t(n) = dble(I*(I+1)) / 6.0d0 * this%ip%E2t
        end do
        n1st = this%ip%Nrot_t
    end select

!$  2nd phonon in the target
    forall (n=n1st+2:this%ip%Nch_t) 
      eps_t(n) = dble(n-n1st-1) * this%ip%omega_t2
    end forall

      
    select case(this%ip%coup_p)
      case(0)
        do n=1, this%ip%Nphonon_p+1
          m = (n - 1) * this%ip%Nch_t
          eps(m+1:m+this%ip%Nch_t) = dble(n - 1) * this%ip%omega_p + eps_t
        end do
        n1st = this%ip%Nphonon_p
      case(1)
        do n=1, this%ip%Nrot_p+1
          m = (n - 1) * this%ip%Nch_t
          I = 2 * (n - 1)
          eps(m+1:m+this%ip%Nch_t) = dble(I*(I+1))/6.0d0*this%ip%E2p + eps_t
        end do
        n1st = this%ip%Nrot_p
    end select
! 2nd mode in the projectile
    do n=n1st+2, this%ip%Nch_p
      m = (n - 1) * this%ip%Nch_t
      eps(m+1:m+this%ip%Nch_t) = dble(n-n1st-1)*this%ip%omega_p2 + eps_t
    end do

    return
  end subroutine
!----------------------------------------------------------------------!
end module


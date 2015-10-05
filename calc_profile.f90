module calc_profile
  use input_data
  use potentials
  use coupling_matrix
  use relative_potential

  type, public :: profiler
    type(inp), private, pointer :: ip
    type(coup_mat), private, pointer :: cm
    type(rel_pot), private, pointer :: pot
    contains
    procedure :: profiler_
    procedure :: pot_prof
    procedure :: Reaction_prof
  end type

  contains

!----------------------------------------------------------------------!
  subroutine profiler_(this, ip, pot, cm)
    implicit none
    class(profiler), intent(out) :: this
    type(inp), intent(in), target :: ip
    type(coup_mat), intent(in), target :: cm
    type(rel_pot), intent(in), target :: pot

    this%ip => ip
    this%cm => cm
    this%pot => pot
  end subroutine
!-----------------------------------------------------------------------
  subroutine pot_prof(this)
    implicit none
    class(profiler), intent(in) :: this
    integer :: i
    real(8) :: r

    open(8,file='potential.dat');
    open(9,file='pot_nucl.dat');
    do i=1, this%ip%rgrid_cut
      r = this%ip%rmin + dble(i-1) * this%ip%dr
      write(8,'(1x,3f8.2)') r, dble(this%pot%get_Vrel(0,i)),aimag(this%pot%get_Vrel(0,i))
      write(9,'(1x,2f8.2)') r, this%pot%Vn(r)
    end do
    close(8)
    close(9)
  end subroutine
!-----------------------------------------------------------------------
  subroutine check_directory(dir)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer :: ios, chdir, getcwd
    character(len=*), intent(in) :: dir
    character(len=300) :: cwd, command
    character(len=1) :: q = 'y'

    
    ios = getcwd(cwd)
    ios = chdir(dir)
    if (ios /= 0) then
      write(output_unit,*) '*************************************************'
      write(output_unit,*) 'Directory '
      write(output_unit,*)
      write(output_unit,*) trim(dir(len(trim(cwd))+2:))
      write(output_unit,*)
      write(output_unit,*) 'does not exist.'
      write(output_unit,*)
      write(output_unit,'(a,$)')' Do you make the directory ? (y/n) : '
      do
        read(5,*) q
        if (trim(q) == 'y' .or. trim(q) == 'Y') then
          command = 'mkdir '//trim(dir)
          call system(command)
          write(output_unit,*)
          write(output_unit,*) 'Directory ',trim(dir(len(trim(cwd))+2:)), &
                     ' is created.'
          exit
        else if (trim(q) == 'n' .or. trim(q) == 'N') then
          write(output_unit,*)
          write(output_unit,*) 'Make the directory ', &
                      trim(dir(len(trim(cwd))+2:))
          write(output_unit,*)
          write(output_unit,*) '*********************************************'
          exit
        else 
          write(output_unit,'(a,$)') " Answer 'y' or 'n' : "
        end if
      end do
    end if
    
    ios = chdir(cwd)

    if (trim(q) == 'n' .or. trim(q) == 'N') stop

    return
  end subroutine
!----------------------------------------------------------------------!
  subroutine Reaction_prof(this, Fname)
    use, intrinsic :: iso_fortran_env
    implicit none
    class(profiler), intent(in) :: this
    integer, parameter :: N = 7
    integer :: m, k, l
    real(8) :: rp, rb
    character(len=*), intent(in) :: Fname
    character(len=20), parameter :: FM='(1x,a,i3," ",a)'
    character(len=20), parameter :: FM2='(1x,a,f7.4,a)'
    character(len=20), parameter :: FM3='(1x,a,i2)'
    character(len=20), parameter :: FM4='(1x,a,f7.4)'
    character(len=20), parameter :: FM5='(1x,a,f9.4,a)'
    character(len=50), parameter :: FM6='(1x,i3,a,i3,a,f7.3," MeV")'
    character(len=20), parameter :: FM8='(1x,a,f11.4,a)'
    character(len=40),parameter::FM7='(1x,i3,a,i3,a,i3,a,f7.3," MeV")'
    character(len=2), save :: X(111)
    data X /'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na',  &
            'Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti',  &
            'V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As',  &
            'Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru',  &
            'Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs',  &
            'Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',  &
            'Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir',  &
            'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',  &
            'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es',  &
            'Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds',  &
            'Rg'/
   
    rb = this%cm%rb
    rp = this%cm%rp
    
    open(N, file=Fname)
    write(output_unit,*)  '------------------------------------------------'
    write(output_unit,FM) 'Target     nucleus : ',nint(this%ip%At), X(nint(this%ip%Zt))
    write(output_unit,FM) 'Projectile nucleus : ',nint(this%ip%Ap), X(nint(this%ip%Zp))
    write(output_unit,*)
    write(N,*)  '------------------------------------------------'
    write(N,FM) ' Target     nucleus : ',nint(this%ip%At), X(nint(this%ip%Zt))
    write(N,FM) ' Projectile nucleus : ',nint(this%ip%Ap), X(nint(this%ip%Zp))
    write(N,*)
    select case(this%ip%coup_t)
      case(-1)
        write(output_unit,*) 'No excitation in the target'
        write(N,*) ' No excitation in the target'
      case(0)
        write(output_unit,*) 'Phonon excitation in the target.'
        write(N,*) ' Phonon excitation in the target :'
        write(N,FM2) '  r0      =', this%ip%r0t     ,' fm'
        write(N,FM2) '  Rt      =', this%ip%Rt    ,' fm'
        write(N,FM4) '  betaN   =', this%ip%betaL_t
        write(N,FM4) '  betaC   =', this%ip%betaL_tC
        write(N,FM2) '  omega   =', this%ip%omega_t , ' MeV'
        write(N,FM3) '  Lambda  =', this%ip%Lambda_t
        write(N,FM3) '  Nphonon =', this%ip%Nphonon_t
      case(1)
        write(output_unit,*) 'Rotational excitation in the target.'
        write(N,*) ' Rotational excitation in the target :'
        write(N,FM2) '  r0      =', this%ip%r0t    ,' fm'
        write(N,FM2) '  Rt      =', this%ip%Rt    ,' fm'
        write(N,FM4) '  beta2   =', this%ip%beta2t
        write(N,FM4) '  beta4   =', this%ip%beta4t
        write(N,FM2) '  E2      =', this%ip%E2t    ,' MeV'
        write(N,FM3) '  Nrot    =', this%ip%Nrot_t
    end select
    write(N,*)
    if (this%ip%coup_t /= -1 .and. this%ip%Nphonon_t2 > 0) then
        write(output_unit,*) 'Phonon excitation in the target.'
        write(N,*) ' Phonon excitation in the target :'
        write(N,FM4) '  betaN   =', this%ip%betaL_t2
        write(N,FM4) '  betaC   =', this%ip%betaL_t2C
        write(N,FM2) '  omega   =', this%ip%omega_t2 , ' MeV'
        write(N,FM3) '  Lambda  =', this%ip%Lambda_t2
        write(N,FM3) '  Nphonon =', this%ip%Nphonon_t2
    end if
    write(N,*)
    select case(this%ip%coup_p)
      case(-1)
        write(output_unit,*) 'No excitation in the projectile'
        write(N,*) ' No excitation in the projectile'
      case(0)
        write(output_unit,*) 'Phonon excitation in the projectile.'
        write(N,*) ' Phonon excitation in the projectile :'
        write(N,FM2) '  r0      =', this%ip%r0p     ,' fm'
        write(N,FM2) '  Rp      =', this%ip%Rp     ,' fm'
        write(N,FM4) '  betaL   =', this%ip%betaL_p
        write(N,FM4) '  betaC   =', this%ip%betaL_pC
        write(N,FM2) '  omega   =', this%ip%omega_p , ' MeV'
        write(N,FM3) '  Lambda  =', this%ip%Lambda_p
        write(N,FM3) '  Nphonon =', this%ip%Nphonon_p
      case(1)
        write(output_unit,*) 'Rotational excitation in the projectile :'
        write(N,*) ' Rotational excitation in the projectile :'
        write(N,FM2) '  r0      =', this%ip%r0p    ,' fm'
        write(N,FM2) '  Rp      =', this%ip%Rp    ,' fm'
        write(N,FM4) '  beta2   =', this%ip%beta2p
        write(N,FM4) '  beta4   =', this%ip%beta4p
        write(N,FM2) '  E2      =', this%ip%E2p    ,' MeV'
        write(N,FM3) '  Nrot    =', this%ip%Nrot_p
    end select
    write(N,*)
    if (this%ip%coup_p /= -1 .and. this%ip%Nphonon_p2 > 0) then
        write(output_unit,*) 'Phonon excitation in the projectile.'
        write(N,*) ' Phonon excitation in the projectile :'
        write(N,FM4) '  betaL   =', this%ip%betaL_p2
        write(N,FM4) '  betaC   =', this%ip%betaL_p2C
        write(N,FM2) '  omega   =', this%ip%omega_p2 , ' MeV'
        write(N,FM3) '  Lambda  =', this%ip%Lambda_p2
        write(N,FM3) '  Nphonon =', this%ip%Nphonon_p2
    end if
    write(N,*)
    write(output_unit,*)
    write(output_unit,'(1x,a,i5)') 'The number of channels = ', this%ip%Nch
    write(output_unit,*) '----------------------------------------------'
    write(N,*)
    write(N,'(1x,a,i5)') ' The number of channels = ', this%ip%Nch
    write(N,*) '----------------------------------------------'
    write(N,*) ' ** Potential Profile (Woods-Saxon) **'
    write(N,FM5) '  V0 =', this%ip%V0, ' MeV'
    write(N,FM5) '  r0 =', this%ip%r0, ' fm'
    write(N,FM5) '  Rn =', this%ip%Rn, ' fm'
    write(N,FM5) '  a  =', this%ip%a, ' fm'
    write(N,*) ' ** Coulomb radius **'
    write(N,FM5) '  r0c=', this%ip%r0c, ' fm'
    write(N,FM5) '  Rc =', this%ip%Rc, ' fm'
    write(N,*)
    write(N,FM5) '  W0 =', this%ip%W0,  ' MeV'
    write(N,FM5) '  r0w=', this%ip%r0w, ' fm'
    write(N,FM5) '  Rw =', this%ip%Rw, ' fm'
    write(N,FM5) '  aw =', this%ip%aw,  ' fm'
    write(N,*)
    write(N,*) ' Minimum of the Coulomb pocket'
    write(N,FM5) '  r    =', rp, ' fm'
    write(N,FM5) '  V(r) =', this%pot%Vn(rp) + this%pot%Vc(rp), ' MeV'
    write(N,*)
    write(N,*) ' Maximum of the Coulomb barrier'
    write(N,FM5) '  r    =', rb, ' fm'
    write(N,FM5) '  V(r) =', this%pot%Vn(rb) + this%pot%Vc(rb), ' MeV'
    write(N,*) '------------------------------------------------'
    write(N,*)
    write(N,*) ' ** The numer of stabilization points **'
    write(N,'(2x,i4,a)') this%ip%num_stab_pt, ' points'
    write(N,*)
    write(N,*) '------------------------------------------------'
    write(N,*) ' ** other information **                      '
    write(N,FM5) ' rmin =', this%ip%rmin, ' fm'
    write(N,FM5) ' rmax =', this%ip%rmax, ' fm'
    write(N,FM5) ' dr   =', this%ip%dr,   ' fm'
    write(N,FM5) ' rcut =', this%ip%rcut,   ' fm'
    write(N,FM5)
    write(N,FM5) ' Emin =', this%ip%Emin, ' MeV'
    write(N,FM5) ' Emax =', this%ip%Emax, ' MeV'
    write(N,FM5) ' dE   =', this%ip%dE,   ' MeV'
    write(N,FM5)
    write(N,*)
    write(N,'(1x,a,i5)') ' Jmax =', this%ip%Jmax
    write(N,*) '----------------------------------------------'
    write(N,*) ' Channels :'
    if (this%ip%coup_p == - 1 .and. this%ip%coup_t == - 1) then
      write(output_unit,*) '1 ch.: n_t = 0, n_p = 0, E* = 0.000 MeV'
      write(N,*) '1 ch.: n_t = 0, n_p = 0, E* = 0.000 MeV'
    else if (this%ip%coup_p == - 1 .and. this%ip%coup_t /= - 1) then
      do m=1, this%ip%Nch
        write(output_unit,FM6) m,'ch.: n_t =',m-1,' n_p = 0, E* =',this%cm%e_n(m)
        write(N,FM6) m,'ch.: n_t =',m-1,' n_p = 0, E* =',this%cm%e_n(m)
      end do
    else if (this%ip%coup_p /= - 1 .and. this%ip%coup_t == - 1) then
      do m=1, this%ip%Nch
        write(output_unit,FM6) m,'ch.: n_t = 0, n_p = ',m-1, ' E* =',this%cm%e_n(m)
        write(N,FM6) m,'ch.: n_t = 0, n_p = ',m-1, ' E* =',this%cm%e_n(m)
      end do
    else
      write(output_unit,*) 'Nch_p =', this%ip%Nch_p, 'Nch_t =', this%ip%Nch_t
      do k=1, this%ip%Nch_p
        L = (k - 1) * this%ip%Nch_t
        do m=L+1, L+this%ip%Nch_t
          write(output_unit,FM7)m,' ch: nt =',m-L-1,', np =',k-1,', E* =',this%cm%e_n(m)
          write(N,FM7)m,' ch: nt =',m-L-1,', np =',k-1,', E* =',this%cm%e_n(m)
        end do
      end do
    end if

    close(N)

    return
  end subroutine
end module


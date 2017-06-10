  program main
    use, intrinsic :: iso_fortran_env
    use global_constant, only : PI, mass
    use input_data
    use relative_potential
    use coupling_matrix
    use coupled_channels
    use calc_profile
    implicit none
    logical :: b
    integer :: i, n, m, ios, getcwd, chdir
    real(8), allocatable, dimension(:) :: Ea, Esig, angl
    real(8), allocatable, dimension(:) :: sig_iel, sig_fus, spin
    real(8), allocatable, dimension(:,:) :: sig_iel_n
    real(8), allocatable, dimension(:,:) :: dsig_iel, dsig_qel, dsig_R
    real(8), allocatable, dimension(:,:,:) :: dsig_qel_n
    complex(8), allocatable, dimension(:,:) :: fc
    complex(8), allocatable, dimension(:,:,:) :: fN
    character(len=500) :: dir, cwd, wd
    character(len=500) :: dir2='angular'
    character(len=500) :: dir3='angular_dist'
    character(len=500) :: dir4="Q_val_dist"
    character(len=100), parameter :: INPF="input_scat"
    type(inp) :: ip
    type(cc_scat) :: cc_calc
    type(rel_pot) :: vrp
    type(coup_mat) :: cm
    type(profiler) :: prof

    call ip%read_input(trim(INPF), dir)

    ios = getcwd(cwd)
    call check_directory(trim(cwd)//'/'//dir)
    call check_directory(trim(cwd)//'/'//trim(dir)//'/'//dir2)
    call check_directory(trim(cwd)//'/'//trim(dir)//'/'//dir3)
    call check_directory(trim(cwd)//'/'//trim(dir)//'/'//dir4)

    call vrp%rel_pot_(ip)
    call cm%coup_mat_(ip)
    call cc_calc%cc_scat_(ip, vrp, cm)
    call prof%profiler_(ip, vrp, cm)

    allocate(Ea(ip%Egrid+1))
    allocate(Esig(ip%Egrid+1))
    allocate(angl(ip%nth))
    allocate(fN(ip%Nch,ip%nth,ip%Egrid+1))
    allocate(fc(ip%nth,ip%Egrid+1))
    allocate(dsig_iel(ip%nth,ip%Egrid+1))
    allocate(dsig_qel(ip%nth,ip%Egrid+1))
    allocate(dsig_qel_n(ip%Nch,ip%nth,ip%Egrid+1))
    allocate(dsig_R(ip%nth,ip%Egrid+1))
    allocate(spin(ip%Egrid+1))
    allocate(sig_fus(ip%Egrid+1))
    allocate(sig_iel_n(ip%Nch,ip%Egrid+1))
    allocate(sig_iel(ip%Egrid+1))

    forall (i=1:ip%Egrid+1) Ea(i) = ip%Emin + dble(i-1) * ip%dE
    forall (i=1:ip%nth)   angl(i) =ip%thmin + dble(i-1) * ip%dth
    angl = angl * PI / 180.0d0

    b = cm%find_rmin(0)
    call vrp%make_Vrel()
    ios = chdir(trim(cwd)//'/'//trim(dir))
    call prof%pot_prof()
    call prof%Reaction_prof("calc_info")
    ios = chdir(trim(cwd))
    write(output_unit,*)
    write(output_unit,*) 'Energy(MeV)  Fusion cross section(mb)'

    call cm%make_Vcp()

    do i=1, ip%Egrid+1
      call cc_calc%cc_scattering(Ea(i), spin(i), sig_fus(i), sig_iel_n(:,i))
      fc(:,i) = cc_calc%scat_amp_Coul(Ea(i), angl)
      call cc_calc%scat_amp_nucl(Ea(i), angl, fN(:,:,i))
      write(output_unit,'(1x f8.4, 4x,es18.6)') Ea(i), sig_fus(i)
    end do

    Esig = Ea * sig_fus
    sig_iel = sum(sig_iel_n(2:ip%Nch,:),dim=1)
    dsig_R = abs(fc) ** 2 * 10.0d0
    dsig_qel_n(1,:,:) = abs(fc + fN(1,:,:)) ** 2 * 10.0d0
    do n=2, ip%Nch
      dsig_qel_n(n,:,:) = abs(fN(n,:,:)) ** 2 * 10.0d0
    end do
    dsig_iel = sum(dsig_qel_n(2:ip%Nch,:,:),dim=1) 
    dsig_qel = dsig_iel + dsig_qel_n(1,:,:)

    ios = chdir(trim(cwd)//'/'//trim(dir))
    call output()
    call qel_corrected_angle()

    ios = chdir(trim(cwd))

    call cc_calc%destruct_cc_scat()
    call vrp%destruct_rel_pot()
    call cm%destruct_coup_mat()
    write(output_unit,*) "calculation finished."
    write(output_unit,*) 'directory : '//trim(dir)

    deallocate(Ea,  Esig, angl)
    deallocate(sig_fus, sig_iel, sig_iel_n, spin)
    deallocate(fc, fN)
    deallocate(dsig_iel, dsig_qel, dsig_qel_n, dsig_R)

    contains
!********************************************************************!
    subroutine output()
      implicit none
      integer :: t, di_u, di_l
      real(8) :: E, Dfus, Dqel
      character(len=45), parameter :: FM1='(1x,a,f8.3,a)'
      character(len=45), parameter :: FM2='(1x,a,es13.4,a)'
      character(len=45), parameter :: FM3='(1x,a,i3,a,es13.4,a)'
      character(len=50), parameter :: FM4='(1x,f8.3,3es13.4)'
      character(len=50), parameter :: FM5='(1x,f7.3,f8.2,3x,a,es14.4)'
      character(len=50), parameter :: FM6='(19x,i3,1x,2es14.4)'
      character(len=50), parameter :: FM7='(8x,f8.2,3x,a,es14.4)'
      character(len=50), parameter :: FM8='(1x,2f10.3,es14.4)'
      character(len=50), parameter :: FM9='(1x,2f10.3,es11.3)'
      character(len=50), parameter :: FM10='(1x,f10.3,es11.3)'
      character(len=500) :: c, c1, c2, c3, c4

!-- Fusion cross section -------------------------------
      open(7,file='fusion.dat')
      do i=1, ip%Egrid+1
        write(7,*) Ea(i), sig_fus(i)
      end do
      close(7)

!-- Fusion barrier distribution -----------------------------
      open(7,file='fus_bar_dist.dat')
      do i=ip%di+1, ip%Egrid-ip%di+1
        Dfus = 2.0d0 * ((Esig(i+ip%di)-Esig(i)) / (Ea(i+ip%di)-Ea(i))  &
                     -  (Esig(i)-Esig(i-ip%di)) / (Ea(i)-Ea(i-ip%di))) &
                     / (Ea(i+ip%di) - Ea(i-ip%di))
        write(7,FM4) Ea(i), Dfus
      end do
      close(7)

!-- Energy dependence of differential scattering cross sections -----
      do n=1, ip%nth
        t = nint(angl(n) * 1800.0d0 / PI)   ! per 10 degree
        if (mod(t,100) == 0) then
          write(c,'(f6.2)') angl(n) * 180.0d0 / PI
          c1 = trim(dir2)//'/qel_'//trim(adjustl(c))//'deg.dat'
!         c2 = trim(dir2)//'/el_'//trim(adjustl(c))//'deg.dat'
!         c3 = trim(dir2)//'/iel_'//trim(adjustl(c))//'deg.dat'
          open(7,file=trim(c1))
!         open(8,file=trim(c2))
!         open(9,file=trim(c3))
          do i=1, ip%Egrid+1
            E = 2.0d0 * Ea(i) * sin(0.5d0 * angl(n)) &
                     / (1.0d0 + sin(0.5d0 * angl(n)))
            write(7,FM9) Ea(i), E, dsig_qel(n,i)/dsig_R(n,i)
!           write(8,FM10) Ea(i), dsig_qel_n(1,n,i)/dsig_R(n,i)
!           write(9,FM10) Ea(i), dsig_iel(n,i)/dsig_R(n,i)
          end do
          close(7)
!         close(8)
!         close(9)
        end if
      end do

!-- qel. barrier distribution -----------------------------
      if (mod(ip%di,2) == 0) then
        di_u = ip%di / 2
        di_l = di_u
      else
        di_u = ceiling(dble(ip%di)/2.0d0)
        di_l = ip%di / 2
      end if
      do n=1, ip%nth
        t = nint(angl(n) * 1800.0d0 / PI)   ! per 10 degree
        if (t > 1100 .and. mod(t,100) == 0) then
          write(c,'(f6.2)') angl(n) * 180.0d0 / PI
         c1 = trim(dir2)//'/qel_bar'//trim(adjustl(c))//'deg.dat'
          open(7,file=trim(c1))
          do i=di_l+1, ip%Egrid-di_u+1
            E = 0.5d0 * (Ea(i+di_u) + Ea(i-di_l))
            E = 2.0d0 * E * sin(0.5d0 * angl(n))            &
                     / (1.0d0 + sin(0.5d0 * angl(n)))
            Dqel = - (dsig_qel(n,i+di_u) / dsig_R(n,i+di_u)   &
                   -  dsig_qel(n,i-di_l) / dsig_R(n,i-di_l))  &
                   / (Ea(i+di_u) - Ea(i-di_l))
            write(7,FM8) Ea(i), E, Dqel
          end do
          close(7)
        end if
      end do


!-- angular distribution of elastic differential scattering cross section
      do i=1, ip%Egrid+1
        if (mod(nint(Ea(i)*10.0d0),50) == 0) then     ! per 5 MeV
          write(c,'(f6.2)') Ea(i)
          open(7,file=trim(dir3)//'/diff_el_'//trim(adjustl(c))//'MeV.dat')
          do n=1, ip%nth
            write(7,*) angl(n)*180.0d0/PI, dsig_qel_n(1,n,i)/dsig_R(n,i)
          end do
          close(7)
        end if
      end do


!-- Q-value distribution for scattering
      do i=1, ip%Egrid+1
        if (mod(nint(Ea(i)*10.0d0),50) == 0) then   ! per 5 MeV
          do n=1, ip%nth
            t = nint(angl(n) * 1800.0d0 / PI)
!           if (t > 1100 .and. mod(t,50) == 0) then    ! per 5 degree
            if (t > 1100 .and. mod(t,100) == 0) then   ! per 10 degree
              write(c1,'(f6.2)') Ea(i)
              write(c2,'(f6.2)') angl(n) * 180.0d0 / PI
              c = trim(dir4)//'/Qdist_'//trim(adjustl(c1))//'MeV_'// &
                  trim(adjustl(c2))//'deg.dat'
              open(7,file=trim(c))
                do m=1, ip%Nch
                  write(7,*) cm%e_n(m), dsig_qel_n(m,n,i), &
                             dsig_qel_n(m,n,i)/dsig_R(n,i)
                end do
              close(7)
            end if
          end do
        end if
      end do


    end subroutine
!**********************************************************************!
    subroutine qel_corrected_angle()
      implicit none
      integer :: i, n, m
      integer :: s, di_u, di_l
      integer :: t(ip%Nch)
      integer, dimension(ip%Nch) :: k
      integer, dimension(ip%nth,ip%Egrid+1) :: kcm2
      real(8), allocatable, dimension(:,:,:) :: dsig_qel_n_lab
      real(8), allocatable, dimension(:,:) :: dsig_R_lab
      real(8), dimension(ip%Nch) :: Kcm, vcm, tcm, tlab
      real(8) :: E, Dqel, Elab, V, E2
      character(len=50), parameter :: FM8='(1x,2f10.3,2es14.4)'
      character(len=50), parameter :: FM9='(1x,2f10.3,2es11.3)'
      character(len=500) :: c, c1, c2, c3, c4

      allocate(dsig_qel_n_lab(ip%Nch,ip%nth,ip%Egrid+1))
      allocate(dsig_R_lab(ip%nth,ip%Egrid+1))
      dsig_qel_n_lab = 0.0d0
      do i=1, ip%Egrid+1
        Kcm = Ea(i) - cm%e_n
        Elab = (ip%Ap + ip%At) / ip%At * Ea(i)
        vcm = sqrt(2.0d0*ip%rmass*Kcm) / (ip%Ap*mass)
        V = sqrt(2.0d0*(ip%Ap*mass)*Elab) / ((ip%Ap+ip%At)*mass)
        do n=1, ip%nth
          tcm = -V/vcm*sin(angl(n))**2 + cos(angl(n))*sqrt(1.0d0-(V/vcm*sin(angl(n)))**2)
          tcm = acos(tcm) * 180.0d0 / PI
          tcm = nint(tcm * 10.0d0) / 10.0d0
          k = nint((tcm - ip%thmin) / ip%dth) + 1
          kcm2(n,i) = k(1)
          do m=1, ip%Nch
            if (k(m) > 0 .and. k(m) <= ip%nth) then
              dsig_qel_n_lab(m,n,i) = dsig_qel_n(m,k(m),i)
            end if
          end do
          dsig_R_lab(n,i) = dsig_R(k(1),i)
        end do
      end do
      dsig_qel = sum(dsig_qel_n_lab(:,:,:),dim=1) 

!-- Energy dependence of differential scattering cross sections -----
      do n=1, ip%nth
        s = nint(angl(n) * 1800.0d0 / PI)   ! per 10 degree
        if (mod(s,100) == 0) then
          write(c,'(f6.2)') angl(n) * 180.0d0 / PI
          c1 = trim(dir2)//'/qel_'//trim(adjustl(c))//'deg_lab.dat'
          open(7,file=trim(c1))
          do i=1, ip%Egrid+1
            E2= 2.0d0 * Ea(i) * sin(0.5d0 * angl(kcm2(n,i))) &
                     / (1.0d0 + sin(0.5d0 * angl(kcm2(n,i))))
            write(7,FM9) Ea(i), E2, dsig_qel(n,i)/dsig_R_lab(n,i)
          end do
          close(7)
        end if
      end do

!-- qel. barrier distribution -----------------------------
      if (mod(ip%di,2) == 0) then
        di_u = ip%di / 2
        di_l = di_u
      else
        di_u = ceiling(dble(ip%di)/2.0d0)
        di_l = ip%di / 2
      end if
      do n=1, ip%nth
        s = nint(angl(n) * 1800.0d0 / PI)   ! per 10 degree
        if (s > 1100 .and. mod(s,100) == 0) then
          write(c,'(f6.2)') angl(n) * 180.0d0 / PI
         c1 = trim(dir2)//'/qel_bar'//trim(adjustl(c))//'deg_lab.dat'
          open(7,file=trim(c1))
          do i=di_l+1, ip%Egrid-di_u+1
            E = 0.5d0 * (Ea(i+di_u) + Ea(i-di_l))
            E2= 2.0d0 * E * sin(0.5d0 * angl(kcm2(n,i))) &
                 / (1.0d0 + sin(0.5d0 * angl(kcm2(n,i))))
            Dqel = - (dsig_qel(n,i+di_u) / dsig_R_lab(n,i+di_u)   &
                   -  dsig_qel(n,i-di_l) / dsig_R_lab(n,i-di_l))  &
                   / (Ea(i+di_u) - Ea(i-di_l))
            write(7,FM8) Ea(i), E2, Dqel
          end do
          close(7)
        end if
      end do

      deallocate(dsig_qel_n_lab)
      deallocate(dsig_R_lab)

    end subroutine
!**********************************************************************!
  end program


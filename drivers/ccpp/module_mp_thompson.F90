!>\file module_mp_thompson.F90
!! This file contains the entity of Thompson MP scheme.

!>\ingroup aathompson

module module_mp_thompson

   use module_mp_thompson_params
   use module_mp_thompson_main
   use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
   use module_mp_radar

#ifdef MPI
   use mpi
#endif

   contains
!>\ingroup aathompson
!! This subroutine calculates simplified cloud species equations and create
!! lookup tables in Thomspson scheme.
!>\section gen_thompson_init thompson_init General Algorithm
!> @{
      subroutine thompson_init(is_aerosol_aware_in,       &
                               merra2_aerosol_aware_in,   &
                               mpicomm, mpirank, mpiroot, &
                               threads, errmsg, errflg)

         implicit none

         logical, intent(in) :: is_aerosol_aware_in
         logical, intent(in) :: merra2_aerosol_aware_in
         integer, intent(in) :: mpicomm, mpirank, mpiroot
         integer, intent(In) :: threads
         character(len=*), intent(inout) :: errmsg
         integer,          intent(inout) :: errflg

         integer :: i, j, k, l, m, n
         logical :: micro_init
         real(wp) :: stime, etime
         logical, parameter :: precomputed_tables = .FALSE.

! Set module variable is_aerosol_aware/merra2_aerosol_aware
         is_aerosol_aware = is_aerosol_aware_in
         merra2_aerosol_aware = merra2_aerosol_aware_in
         if (is_aerosol_aware .and. merra2_aerosol_aware) then
            errmsg = 'Logic error in thompson_init: only one of the two options can be true, ' // &
                     'not both: is_aerosol_aware or merra2_aerosol_aware'
            errflg = 1
            return
         end if
         if (mpirank==mpiroot) then
            if (is_aerosol_aware) then
               write (*,'(a)') 'Using aerosol-aware version of Thompson microphysics'
            else if(merra2_aerosol_aware) then
               write (*,'(a)') 'Using merra2 aerosol-aware version of Thompson microphysics'
            else
               write (*,'(a)') 'Using non-aerosol-aware version of Thompson microphysics'
            end if
         end if

         av_g(idx_bg1) = av_g_old
         bv_g(idx_bg1) = bv_g_old
         dimNRHG = NRHG1
         
         micro_init = .FALSE.

!> - Allocate space for lookup tables (J. Michalakes 2009Jun08).

         if (.NOT. ALLOCATED(tcg_racg) ) then
            ALLOCATE(tcg_racg(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
            micro_init = .TRUE.
         endif

         if (.NOT. ALLOCATED(tmr_racg)) ALLOCATE(tmr_racg(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tcr_gacr)) ALLOCATE(tcr_gacr(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
         ! if (.NOT. ALLOCATED(tmg_gacr)) ALLOCATE(tmg_gacr(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tnr_racg)) ALLOCATE(tnr_racg(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tnr_gacr)) ALLOCATE(tnr_gacr(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))

         if (.NOT. ALLOCATED(tcs_racs1)) ALLOCATE(tcs_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tmr_racs1)) ALLOCATE(tmr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tcs_racs2)) ALLOCATE(tcs_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tmr_racs2)) ALLOCATE(tmr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tcr_sacr1)) ALLOCATE(tcr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tms_sacr1)) ALLOCATE(tms_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tcr_sacr2)) ALLOCATE(tcr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tms_sacr2)) ALLOCATE(tms_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tnr_racs1)) ALLOCATE(tnr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tnr_racs2)) ALLOCATE(tnr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tnr_sacr1)) ALLOCATE(tnr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
         if (.NOT. ALLOCATED(tnr_sacr2)) ALLOCATE(tnr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))

         if (.NOT. ALLOCATED(tpi_qcfz)) ALLOCATE(tpi_qcfz(ntb_c,nbc,45,ntb_IN))
         if (.NOT. ALLOCATED(tni_qcfz)) ALLOCATE(tni_qcfz(ntb_c,nbc,45,ntb_IN))

         if (.NOT. ALLOCATED(tpi_qrfz)) ALLOCATE(tpi_qrfz(ntb_r,ntb_r1,45,ntb_IN))
         if (.NOT. ALLOCATED(tpg_qrfz)) ALLOCATE(tpg_qrfz(ntb_r,ntb_r1,45,ntb_IN))
         if (.NOT. ALLOCATED(tni_qrfz)) ALLOCATE(tni_qrfz(ntb_r,ntb_r1,45,ntb_IN))
         if (.NOT. ALLOCATED(tnr_qrfz)) ALLOCATE(tnr_qrfz(ntb_r,ntb_r1,45,ntb_IN))

         if (.NOT. ALLOCATED(tps_iaus)) ALLOCATE(tps_iaus(ntb_i,ntb_i1))
         if (.NOT. ALLOCATED(tni_iaus)) ALLOCATE(tni_iaus(ntb_i,ntb_i1))
         if (.NOT. ALLOCATED(tpi_ide)) ALLOCATE(tpi_ide(ntb_i,ntb_i1))

         if (.NOT. ALLOCATED(t_Efrw)) ALLOCATE(t_Efrw(nbr,nbc))
         if (.NOT. ALLOCATED(t_Efsw)) ALLOCATE(t_Efsw(nbs,nbc))

         if (.NOT. ALLOCATED(tnr_rev)) ALLOCATE(tnr_rev(nbr, ntb_r1, ntb_r))
         if (.NOT. ALLOCATED(tpc_wev)) ALLOCATE(tpc_wev(nbc,ntb_c,nbc))
         if (.NOT. ALLOCATED(tnc_wev)) ALLOCATE(tnc_wev(nbc,ntb_c,nbc))

         if (.NOT. ALLOCATED(tnccn_act))                                   &
               ALLOCATE(tnccn_act(ntb_arc,ntb_arw,ntb_art,ntb_arr,ntb_ark))

         if_micro_init: if (micro_init) then

!> - From Martin et al. (1994), assign gamma shape parameter mu for cloud
!! drops according to general dispersion characteristics (disp=~0.25
!! for maritime and 0.45 for continental)
!.. disp=SQRT((mu+2)/(mu+1) - 1) so mu varies from 15 for Maritime
!.. to 2 for really dirty air.  This not used in 2-moment cloud water
!.. scheme and nu_c used instead and varies from 2 to 15 (integer-only).
         mu_c_l = min(15.0_wp, (1000.e6/Nt_c_l + 2.))
         mu_c_o = min(15.0_wp, (1000.e6/Nt_c_o + 2.))

!> - Compute Schmidt number to one-third used numerous times
         Sc3 = Sc**(1./3.)

!> - Compute minimum ice diam from mass, min snow/graupel mass from diam
         D0i = (xm0i/am_i)**(1./bm_i)
         xm0s = am_s * D0s**bm_s
         xm0g = am_g(NRHG) * D0g**bm_g

!> - Compute constants various exponents and gamma() associated with cloud,
!! rain, snow, and graupel
         do n = 1, 15
            cce(1,n) = n + 1.
            cce(2,n) = bm_r + n + 1.
            cce(3,n) = bm_r + n + 4.
            cce(4,n) = n + bv_c + 1.
            cce(5,n) = bm_r + n + bv_c + 1.
            ccg(1,n) = WGAMMA(cce(1,n))
            ccg(2,n) = WGAMMA(cce(2,n))
            ccg(3,n) = WGAMMA(cce(3,n))
            ccg(4,n) = WGAMMA(cce(4,n))
            ccg(5,n) = WGAMMA(cce(5,n))
            ocg1(n) = 1./ccg(1,n)
            ocg2(n) = 1./ccg(2,n)
         enddo

         cie(1) = mu_i + 1.
         cie(2) = bm_i + mu_i + 1.
         cie(3) = bm_i + mu_i + bv_i + 1.
         cie(4) = mu_i + bv_i + 1.
         cie(5) = mu_i + 2.
         cie(6) = bm_i*0.5 + mu_i + bv_i + 1.
         cie(7) = bm_i*0.5 + mu_i + 1.
         cig(1) = WGAMMA(cie(1))
         cig(2) = WGAMMA(cie(2))
         cig(3) = WGAMMA(cie(3))
         cig(4) = WGAMMA(cie(4))
         cig(5) = WGAMMA(cie(5))
         cig(6) = WGAMMA(cie(6))
         cig(7) = WGAMMA(cie(7))
         oig1 = 1./cig(1)
         oig2 = 1./cig(2)
         obmi = 1./bm_i

         cre(1) = bm_r + 1.
         cre(2) = mu_r + 1.
         cre(3) = bm_r + mu_r + 1.
         cre(4) = bm_r*2. + mu_r + 1.
         cre(5) = mu_r + bv_r + 1.
         cre(6) = bm_r + mu_r + bv_r + 1.
         cre(7) = bm_r*0.5 + mu_r + bv_r + 1.
         cre(8) = bm_r + mu_r + bv_r + 3.
         cre(9) = mu_r + bv_r + 3.
         cre(10) = mu_r + 2.
         cre(11) = 0.5*(bv_r + 5. + 2.*mu_r)
         cre(12) = bm_r*0.5 + mu_r + 1.
         cre(13) = bm_r*2. + mu_r + bv_r + 1.
         do n = 1, 13
            crg(n) = WGAMMA(cre(n))
         enddo
         obmr = 1./bm_r
         ore1 = 1./cre(1)
         org1 = 1./crg(1)
         org2 = 1./crg(2)
         org3 = 1./crg(3)

         cse(1) = bm_s + 1.
         cse(2) = bm_s + 2.
         cse(3) = bm_s*2.
         cse(4) = bm_s + bv_s + 1.
         cse(5) = bm_s*2. + bv_s + 1.
         cse(6) = bm_s*2. + 1.
         cse(7) = bm_s + mu_s + 1.
         cse(8) = bm_s + mu_s + 2.
         cse(9) = bm_s + mu_s + 3.
         cse(10) = bm_s + mu_s + bv_s + 1.
         cse(11) = bm_s*2. + mu_s + bv_s + 1.
         cse(12) = bm_s*2. + mu_s + 1.
         cse(13) = bv_s + 2.
         cse(14) = bm_s + bv_s
         cse(15) = mu_s + 1.
         cse(16) = 1.0 + (1.0 + bv_s)/2.
         cse(17) = bm_s + bv_s + 2.

         do n = 1, 17
            csg(n) = WGAMMA(cse(n))
         enddo
         oams = 1./am_s
         obms = 1./bm_s
         ocms = oams**obms

         cge(1,:) = bm_g + 1.
         cge(2,:) = mu_g + 1.
         cge(3,:) = bm_g + mu_g + 1.
         cge(4,:) = bm_g*2. + mu_g + 1.
         cge(10,:) = mu_g + 2.
         cge(12,:) = bm_g*0.5 + mu_g + 1.
         do m = 1, NRHG
             cge(5,m) = bm_g*2. + mu_g + bv_g(m) + 1.
             cge(6,m) = bm_g + mu_g + bv_g(m) + 1.
             cge(7,m) = bm_g*0.5 + mu_g + bv_g(m) + 1.
             cge(8,m) = mu_g + bv_g(m) + 1.                                  ! not used
             cge(9,m) = mu_g + bv_g(m) + 3.
             cge(11,m) = 0.5*(bv_g(m) + 5. + 2.*mu_g)
         enddo
         do m = 1, NRHG
             do n = 1, 12
                 cgg(n,m) = WGAMMA(cge(n,m))
             enddo
         enddo
         oamg = 1./am_g
         obmg = 1./bm_g
         do m = 1, NRHG
             oamg(m) = 1./am_g(m)
             ocmg(m) = oamg(m)**obmg
         enddo
         oge1 = 1./cge(1,1)
         ogg1 = 1./cgg(1,1)
         ogg2 = 1./cgg(2,1)
         ogg3 = 1./cgg(3,1)

!+---+-----------------------------------------------------------------+
!> - Simplify various rate equations
!+---+-----------------------------------------------------------------+

!>  - Compute rain collecting cloud water and cloud ice
         t1_qr_qc = PI*.25*av_r * crg(9)
         t1_qr_qi = PI*.25*av_r * crg(9)
         t2_qr_qi = PI*.25*am_r*av_r * crg(8)

!>  - Compute graupel collecting cloud water
!         t1_qg_qc = PI*.25*av_g * cgg(9)

!>  - Compute snow collecting cloud water
         t1_qs_qc = PI*.25*av_s

!>  - Compute snow collecting cloud ice
         t1_qs_qi = PI*.25*av_s

!>  - Compute evaporation of rain; ignore depositional growth of rain
         t1_qr_ev = 0.78 * crg(10)
         t2_qr_ev = 0.308*Sc3*SQRT(av_r) * crg(11)

!>  - Compute sublimation/depositional growth of snow
         t1_qs_sd = 0.86
         t2_qs_sd = 0.28*Sc3*SQRT(av_s)

!>  - Compute melting of snow
         t1_qs_me = PI*4.*C_sqrd*olfus * 0.86
         t2_qs_me = PI*4.*C_sqrd*olfus * 0.28*Sc3*SQRT(av_s)

!>  - Compute sublimation/depositional growth of graupel
         t1_qg_sd = 0.86 * cgg(10,1)
!         t2_qg_sd = 0.28*Sc3*SQRT(av_g) * cgg(11)

!>  - Compute melting of graupel
         t1_qg_me = PI*4.*C_cube*olfus * 0.86 * cgg(10,1)
!         t2_qg_me = PI*4.*C_cube*olfus * 0.28*Sc3*SQRT(av_g) * cgg(11)

!>  - Compute constants for helping find lookup table indexes
         nic2 = nint(log10(r_c(1)))
         nii2 = nint(log10(r_i(1)))
         nii3 = nint(log10(Nt_i(1)))
         nir2 = nint(log10(r_r(1)))
         nir3 = nint(log10(N0r_exp(1)))
         nis2 = nint(log10(r_s(1)))
         nig2 = nint(log10(r_g(1)))
         nig3 = nint(log10(N0g_exp(1)))
         niIN2 = nint(log10(Nt_IN(1)))

!>  - Create bins of cloud water (from min diameter up to 100 microns)
         Dc(1) = D0c*1.0_dp
         dtc(1) = D0c*1.0_dp
         do n = 2, nbc
            Dc(n) = Dc(n-1) + 1.e-6_dp
            dtc(n) = (Dc(n) - Dc(n-1))
         enddo

!>  - Create bins of cloud ice (from min diameter up to 2x min snow size)
         xDx(1) = D0i*1.0_dp
         xDx(nbi+1) = D0s*2.0_dp
         do n = 2, nbi
            xDx(n) = exp(real(n-1, kind=dp)/real(nbi, kind=dp) &
                     *log(real(xDx(nbi+1)/xDx(1), kind=dp)) + log(real(xDx(1), kind=dp)))
         enddo
         do n = 1, nbi
            Di(n) = sqrt(real(xDx(n)*xDx(n+1), kind=dp))
            dti(n) = xDx(n+1) - xDx(n)
         enddo

!>  - Create bins of rain (from min diameter up to 5 mm)
         xDx(1) = D0r*1.0_dp
         xDx(nbr+1) = 0.005_dp
         do n = 2, nbr
            xDx(n) = exp(real(n-1, kind=dp)/real(nbr, kind=dp) &
                     *log(real(xDx(nbr+1)/xDx(1), kind=dp)) + log(real(xDx(1), kind=dp)))
         enddo
         do n = 1, nbr
            Dr(n) = sqrt(real(xDx(n)*xDx(n+1), kind=dp))
            dtr(n) = xDx(n+1) - xDx(n)
         enddo

!>  - Create bins of snow (from min diameter up to 2 cm)
         xDx(1) = D0s*1.0_dp
         xDx(nbs+1) = 0.02_dp
         do n = 2, nbs
            xDx(n) = exp(real(n-1, kind=dp)/real(nbs, kind=dp) &
                     *log(real(xDx(nbs+1)/xDx(1), kind=dp)) + log(real(xDx(1), kind=dp)))
         enddo
         do n = 1, nbs
            Ds(n) = sqrt(real(xDx(n)*xDx(n+1), kind=dp))
            dts(n) = xDx(n+1) - xDx(n)
         enddo

!>  - Create bins of graupel (from min diameter up to 5 cm)
         xDx(1) = D0g*1.0_dp
         xDx(nbg+1) = 0.05_dp
         do n = 2, nbg
            xDx(n) = exp(real(n-1, kind=dp)/real(nbg, kind=dp) &
                     *log(real(xDx(nbg+1)/xDx(1), kind=dp)) + log(real(xDx(1), kind=dp)))
         enddo
         do n = 1, nbg
            Dg(n) = sqrt(real(xDx(n)*xDx(n+1), kind=dp))
            dtg(n) = xDx(n+1) - xDx(n)
         enddo

!>  - Create bins of cloud droplet number concentration (1 to 3000 per cc)
         xDx(1) = 1.0_dp
         xDx(nbc+1) = 3000.0_dp
         do n = 2, nbc
            xDx(n) = exp(real(n-1, kind=dp)/real(nbc, kind=dp)                          &
                     *log(real(xDx(nbc+1)/xDx(1), kind=dp)) + log(real(xDx(1), kind=dp)))
         enddo
         do n = 1, nbc
            t_Nc(n) = sqrt(real(xDx(n)*xDx(n+1), kind=dp)) * 1.e6_dp
         enddo
         nic1 = log(real(t_Nc(nbc)/t_Nc(1), kind=dp))

!+---+-----------------------------------------------------------------+
!> - Create lookup tables for most costly calculations
!+---+-----------------------------------------------------------------+

! Assign mpicomm to module variable
         mpi_communicator = mpicomm

! Standard tables are only written by master MPI task;
! (physics init cannot be called by multiple threads,
!  hence no need to test for a specific thread number)
         if (mpirank==mpiroot) then
            thompson_table_writer = .true.
         else
            thompson_table_writer = .false.
         end if

         precomputed_tables_1: if (.not.precomputed_tables) then

         call cpu_time(stime)

         do m = 1, ntb_r
            do k = 1, ntb_r1
               do n = 1, dimNRHG

               do j = 1, ntb_g
                  do i = 1, ntb_g1
                     tcg_racg(i,j,n,k,m) = 0.0_dp
                     tmr_racg(i,j,n,k,m) = 0.0_dp
                     tcr_gacr(i,j,n,k,m) = 0.0_dp
                     !tmg_gacr(i,j,k,m) = 0.0_dp
                     tnr_racg(i,j,n,k,m) = 0.0_dp
                     tnr_gacr(i,j,n,k,m) = 0.0_dp
                  enddo
                  enddo
               enddo
            enddo
         enddo

         do m = 1, ntb_r
            do k = 1, ntb_r1
               do j = 1, ntb_t
                  do i = 1, ntb_s
                     tcs_racs1(i,j,k,m) = 0.0_dp
                     tmr_racs1(i,j,k,m) = 0.0_dp
                     tcs_racs2(i,j,k,m) = 0.0_dp
                     tmr_racs2(i,j,k,m) = 0.0_dp
                     tcr_sacr1(i,j,k,m) = 0.0_dp
                     tms_sacr1(i,j,k,m) = 0.0_dp
                     tcr_sacr2(i,j,k,m) = 0.0_dp
                     tms_sacr2(i,j,k,m) = 0.0_dp
                     tnr_racs1(i,j,k,m) = 0.0_dp
                     tnr_racs2(i,j,k,m) = 0.0_dp
                     tnr_sacr1(i,j,k,m) = 0.0_dp
                     tnr_sacr2(i,j,k,m) = 0.0_dp
                  enddo
               enddo
            enddo
         enddo

         do m = 1, ntb_IN
            do k = 1, 45
               do j = 1, ntb_r1
                  do i = 1, ntb_r
                     tpi_qrfz(i,j,k,m) = 0.0_dp
                     tni_qrfz(i,j,k,m) = 0.0_dp
                     tpg_qrfz(i,j,k,m) = 0.0_dp
                     tnr_qrfz(i,j,k,m) = 0.0_dp
                  enddo
               enddo
               do j = 1, nbc
                  do i = 1, ntb_c
                     tpi_qcfz(i,j,k,m) = 0.0_dp
                     tni_qcfz(i,j,k,m) = 0.0_dp
                  enddo
               enddo
            enddo
         enddo

         do j = 1, ntb_i1
            do i = 1, ntb_i
               tps_iaus(i,j) = 0.0_dp
               tni_iaus(i,j) = 0.0_dp
               tpi_ide(i,j) = 0.0_dp
            enddo
         enddo

         do j = 1, nbc
            do i = 1, nbr
               t_Efrw(i,j) = 0.0
            enddo
            do i = 1, nbs
               t_Efsw(i,j) = 0.0
            enddo
         enddo

         do k = 1, ntb_r
            do j = 1, ntb_r1
               do i = 1, nbr
                  tnr_rev(i,j,k) = 0.0_dp
               enddo
            enddo
         enddo

         do k = 1, nbc
            do j = 1, ntb_c
               do i = 1, nbc
                  tpc_wev(i,j,k) = 0.0_dp
                  tnc_wev(i,j,k) = 0.0_dp
               enddo
            enddo
         enddo

         do m = 1, ntb_ark
            do l = 1, ntb_arr
               do k = 1, ntb_art
                  do j = 1, ntb_arw
                     do i = 1, ntb_arc
                        tnccn_act(i,j,k,l,m) = 1.0
                     enddo
                  enddo
               enddo
            enddo
         enddo

         if (mpirank==mpiroot) write (*,*)'creating microphysics lookup tables ... '
         if (mpirank==mpiroot) write (*,'(a, f5.2, a, f5.2, a, f5.2, a, f5.2)') &
            ' using: mu_c_o=',mu_c_o,' mu_i=',mu_i,' mu_r=',mu_r,' mu_g=',mu_g

!>  - Call table_ccnact() to read a static file containing CCN activation of aerosols. The
!! data were created from a parcel model by Feingold & Heymsfield with
!! further changes by Eidhammer and Kriedenweis
         if (mpirank==mpiroot) write(*,*) '  calling table_ccnAct routine'
         call table_ccnAct(errmsg,errflg)
         if (.not. errflg==0) return

!>  - Call table_efrw() and table_efsw() to creat collision efficiency table 
!! between rain/snow and cloud water
         if (mpirank==mpiroot) write(*,*) '  creating qc collision eff tables'
         call table_Efrw
         call table_Efsw

!>  - Call table_dropevap() to creat rain drop evaporation table
         if (mpirank==mpiroot) write(*,*) '  creating rain evap table'
         call table_dropEvap

!>  - Call qi_aut_qs() to create conversion of some ice mass into snow category
         if (mpirank==mpiroot) write(*,*) '  creating ice converting to snow table'
         call qi_aut_qs

         call cpu_time(etime)
         if (mpirank==mpiroot) print '("Calculating Thompson tables part 1 took ",f10.3," seconds.")', etime-stime

         end if precomputed_tables_1

!>  - Call radar_init() to initialize various constants for computing radar reflectivity
         call cpu_time(stime)
         xam_r = am_r
         xbm_r = bm_r
         xmu_r = mu_r
         xam_s = am_s
         xbm_s = bm_s
         xmu_s = mu_s
         xam_g = am_g(idx_bg1)
         xbm_g = bm_g
         xmu_g = mu_g
         call radar_init
         call cpu_time(etime)
         if (mpirank==mpiroot) print '("Calling radar_init took ",f10.3," seconds.")', etime-stime


         if_not_iiwarm: if (.not. iiwarm) then

         precomputed_tables_2: if (.not.precomputed_tables) then

         call cpu_time(stime)

!>  - Call qr_acr_qg() to create rain collecting graupel & graupel collecting rain table
         if (mpirank==mpiroot) write(*,*) '  creating rain collecting graupel table'
         call cpu_time(stime)
         call qr_acr_qg(1)
         call cpu_time(etime)
         if (mpirank==mpiroot) print '("Computing rain collecting graupel table took ",f10.3," seconds.")', etime-stime

!>  - Call qr_acr_qs() to create rain collecting snow & snow collecting rain table
         if (mpirank==mpiroot) write (*,*) '  creating rain collecting snow table'
         call cpu_time(stime)
         call qr_acr_qs
         call cpu_time(etime)
         if (mpirank==mpiroot) print '("Computing rain collecting snow table took ",f10.3," seconds.")', etime-stime

!>  - Call freezeh2o() to create cloud water and rain freezing (Bigg, 1953) table
         if (mpirank==mpiroot) write(*,*) '  creating freezing of water drops table'
         call cpu_time(stime)
         call freezeH2O(threads)
         call cpu_time(etime)
         if (mpirank==mpiroot) print '("Computing freezing of water drops table took ",f10.3," seconds.")', etime-stime

         call cpu_time(etime)
         if (mpirank==mpiroot) print '("Calculating Thompson tables part 2 took ",f10.3," seconds.")', etime-stime

         end if precomputed_tables_2

         endif if_not_iiwarm

         if (mpirank==mpiroot) write(*,*) ' ... DONE microphysical lookup tables'

         endif if_micro_init

      end subroutine thompson_init
!> @}

!>\ingroup aathompson
!!This is a wrapper routine designed to transfer values from 3D to 1D.
!!\section gen_mpgtdriver Thompson mp_gt_driver General Algorithm
!> @{
      subroutine mp_gt_driver(qv, qc, qr, qi, qs, qg, qb, ni, nr, nc, ng,&
                              nwfa, nifa, nwfa2d, nifa2d,             &
                              tt, th, pii,                            &
                              p, w, dz, dt_in, dt_inner,              &
                              sedi_semi, decfl, lsm,                  &
                              RAINNC, RAINNCV,                        &
                              SNOWNC, SNOWNCV,                        &
                              ICENC, ICENCV,                          &
                              GRAUPELNC, GRAUPELNCV, SR,              &
#if ( WRF_CHEM == 1 )
                              rainprod, evapprod,                     &
#endif
                              refl_10cm, diagflag, do_radar_ref,      &
                              max_hail_diam_sfc,                      &
                              vt_dbz_wt, first_time_step,             &
                              re_cloud, re_ice, re_snow,              &
                              has_reqc, has_reqi, has_reqs,           &
                              aero_ind_fdb, rand_perturb_on,          &
                              kme_stoch,                              &
                              rand_pert, spp_prt_list, spp_var_list,  &
                              spp_stddev_cutoff, n_var_spp,           &
                              ids,ide, jds,jde, kds,kde,              &  ! domain dims
                              ims,ime, jms,jme, kms,kme,              &  ! memory dims
                              its,ite, jts,jte, kts,kte,              &  ! tile dims
                              fullradar_diag, istep, nsteps,          &
                              errmsg, errflg,                         &
                              ! Extended diagnostics, array pointers
                              ! only associated if ext_diag flag is .true.
                              ext_diag,                               &
                              !vts1, txri, txrc,                       &
                              prw_vcdc,                               &
                              prw_vcde, tpri_inu, tpri_ide_d,         &
                              tpri_ide_s, tprs_ide, tprs_sde_d,       &
                              tprs_sde_s, tprg_gde_d,                 &
                              tprg_gde_s, tpri_iha, tpri_wfz,         &
                              tpri_rfz, tprg_rfz, tprs_scw, tprg_scw, &
                              tprg_rcs, tprs_rcs,                     &
                              tprr_rci, tprg_rcg,                     &
                              tprw_vcd_c, tprw_vcd_e, tprr_sml,       &
                              tprr_gml, tprr_rcg,                     &
                              tprr_rcs, tprv_rev, tten3, qvten3,      &
                              qrten3, qsten3, qgten3, qiten3, niten3, &
                              nrten3, ncten3, qcten3,                 &
                              pfils, pflls)

         implicit none

!..Subroutine arguments
         integer, intent(in):: ids,ide, jds,jde, kds,kde, &
                              ims,ime, jms,jme, kms,kme, &
                              its,ite, jts,jte, kts,kte
         real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout):: &
                           qv, qc, qr, qi, qs, qg, ni, nr
         real(wp), dimension(ims:ime, kms:kme, jms:jme), optional, intent(inout):: &
                           tt, th
         real(wp), dimension(ims:ime, kms:kme, jms:jme), optional, intent(in):: &
                           pii
         real(wp), dimension(ims:ime, kms:kme, jms:jme), optional, intent(inout):: &
                           nc, nwfa, nifa, qb, ng
         real(wp), dimension(ims:ime, jms:jme), optional, intent(in):: nwfa2d, nifa2d
         integer, dimension(ims:ime, jms:jme), intent(in):: lsm
         real(wp), dimension(ims:ime, kms:kme, jms:jme), optional, intent(inout):: &
                           re_cloud, re_ice, re_snow
         real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout):: pfils, pflls
         integer, intent(in) :: rand_perturb_on, kme_stoch, n_var_spp
         real(wp), dimension(:,:), intent(in) :: rand_pert
         real(wp), dimension(:), intent(in) :: spp_prt_list, spp_stddev_cutoff
         character(len=10), dimension(:), intent(in) :: spp_var_list
         integer, intent(in):: has_reqc, has_reqi, has_reqs
#if ( WRF_CHEM == 1 )
         real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout):: &
                           rainprod, evapprod
#endif
         real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in):: &
                           p, w, dz
         real(wp), dimension(ims:ime, jms:jme), intent(inout):: &
                           RAINNC, RAINNCV, SR
         real(wp), dimension(ims:ime, jms:jme), optional, intent(inout)::      &
                           SNOWNC, SNOWNCV,                              &
                           ICENC, ICENCV,                                &
                           GRAUPELNC, GRAUPELNCV
         real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout)::       &
                           refl_10cm
         real(wp), dimension(ims:ime, jms:jme), intent(inout)::       &
                           max_hail_diam_sfc
         real(wp), dimension(ims:ime, kms:kme, jms:jme), optional, intent(inout):: &
                           vt_dbz_wt
         logical, intent(in) :: first_time_step
         real(wp), intent(in):: dt_in, dt_inner
         logical, intent(in) :: sedi_semi
         integer, intent(in) :: decfl
         ! To support subcycling: current step and maximum number of steps
         integer, intent (in) :: istep, nsteps
         logical, intent (in) :: fullradar_diag 
         ! Extended diagnostics, array pointers only associated if ext_diag flag is .true.
         logical, intent (in) :: ext_diag
         logical, optional, intent(in):: aero_ind_fdb
         real(wp), dimension(:,:,:), intent(inout)::                     &
                           !vts1, txri, txrc,                       &
                           prw_vcdc,                               &
                           prw_vcde, tpri_inu, tpri_ide_d,         &
                           tpri_ide_s, tprs_ide,                   &
                           tprs_sde_d, tprs_sde_s, tprg_gde_d,     &
                           tprg_gde_s, tpri_iha, tpri_wfz,         &
                           tpri_rfz, tprg_rfz, tprs_scw, tprg_scw, &
                           tprg_rcs, tprs_rcs,                     &
                           tprr_rci, tprg_rcg,                     &
                           tprw_vcd_c, tprw_vcd_e, tprr_sml,       &
                           tprr_gml, tprr_rcg,                     &
                           tprr_rcs, tprv_rev, tten3, qvten3,      &
                           qrten3, qsten3, qgten3, qiten3, niten3, &
                           nrten3, ncten3, qcten3

   !..Local variables
         real(wp), dimension(kts:kte):: &
                           qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, qb1d, &
                           ni1d, nr1d, nc1d, ng1d, nwfa1d, nifa1d, &
                           t1d, p1d, w1d, dz1d, rho, dBZ, pfil1, pfll1
   !..Extended diagnostics, single column arrays
         real(wp), dimension(:), allocatable::                              &
                           !vtsk1, txri1, txrc1,                       &
                           prw_vcdc1,                                 &
                           prw_vcde1, tpri_inu1, tpri_ide1_d,         &
                           tpri_ide1_s, tprs_ide1,                    &
                           tprs_sde1_d, tprs_sde1_s, tprg_gde1_d,     &
                           tprg_gde1_s, tpri_iha1, tpri_wfz1,         &
                           tpri_rfz1, tprg_rfz1, tprs_scw1, tprg_scw1,&
                           tprg_rcs1, tprs_rcs1,                      &
                           tprr_rci1, tprg_rcg1,                      &
                           tprw_vcd1_c, tprw_vcd1_e, tprr_sml1,       &
                           tprr_gml1, tprr_rcg1,                      &
                           tprr_rcs1, tprv_rev1,  tten1, qvten1,      &
                           qrten1, qsten1, qgten1, qiten1, niten1,    &
                           nrten1, ncten1, qcten1

         real(wp), dimension(kts:kte):: re_qc1d, re_qi1d, re_qs1d
#if ( WRF_CHEM == 1 )
      real(wp), dimension(kts:kte):: &
                        rainprod1d, evapprod1d
#endif
         real(wp), dimension(its:ite, jts:jte):: pcp_ra, pcp_sn, pcp_gr, pcp_ic
         real(wp) :: dt, pptrain, pptsnow, pptgraul, pptice
         real(wp) :: qc_max, qr_max, qs_max, qi_max, qg_max, ni_max, nr_max
         real(wp) :: ygra1, zans1
         real(dp) :: lamg, lam_exp, lamr, N0_min, N0_exp
         integer:: lsml
         real(wp) :: rand1, rand2, rand3, rand_pert_max
         integer:: i, j, k, m
         integer:: imax_qc,imax_qr,imax_qi,imax_qs,imax_qg,imax_ni,imax_nr
         integer:: jmax_qc,jmax_qr,jmax_qi,jmax_qs,jmax_qg,jmax_ni,jmax_nr
         integer:: kmax_qc,kmax_qr,kmax_qi,kmax_qs,kmax_qg,kmax_ni,kmax_nr
         integer:: i_start, j_start, i_end, j_end
         logical, optional, intent(in) :: diagflag
         integer, optional, intent(in) :: do_radar_ref
         logical :: melti = .false.
         integer :: ndt, it

         ! CCPP error handling
         character(len=*), optional, intent(  out) :: errmsg
         integer,          optional, intent(  out) :: errflg

         ! CCPP
         if (present(errmsg)) errmsg = ''
         if (present(errflg)) errflg = 0

         ! No need to test for every subcycling step
         test_only_once: if (first_time_step .and. istep==1) then
            ! Activate this code when removing the guard above
      
            if ( (present(tt) .and. (present(th) .or. present(pii))) .or. &
               (.not.present(tt) .and. .not.(present(th) .and. present(pii))) ) then
               if (present(errmsg) .and. present(errflg)) then
                  write(errmsg, '(a)') 'Logic error in mp_gt_driver: provide either tt or th+pii'
                  errflg = 1
                  return
               else
                  write(*,'(a)') 'Logic error in mp_gt_driver: provide either tt or th+pii'
                  stop
               end if
            end if
   
            if (is_aerosol_aware .and. (.not.present(nc)     .or. &
                                       .not.present(nwfa)   .or. &
                                       .not.present(nifa)   .or. &
                                       .not.present(nwfa2d) .or. &
                                       .not.present(nifa2d)      )) then
               if (present(errmsg) .and. present(errflg)) then
                  write(errmsg, '(*(a))') 'Logic error in mp_gt_driver: provide nc, nwfa, nifa, nwfa2d', &
                                          ' and nifa2d for aerosol-aware version of Thompson microphysics'
                  errflg = 1
                  return
               else
                  write(*, '(*(a))') 'Logic error in mp_gt_driver: provide nc, nwfa, nifa, nwfa2d', &
                                    ' and nifa2d for aerosol-aware version of Thompson microphysics'
                  stop
               end if
            else if (merra2_aerosol_aware .and. (.not.present(nc)   .or. &
                                                .not.present(nwfa) .or. &
                                                .not.present(nifa)      )) then
               if (present(errmsg) .and. present(errflg)) then
                  write(errmsg, '(*(a))') 'Logic error in mp_gt_driver: provide nc, nwfa, and nifa', &
                                          ' for merra2 aerosol-aware version of Thompson microphysics'
                  errflg = 1
                  return
               else
                  write(*, '(*(a))') 'Logic error in mp_gt_driver: provide nc, nwfa, and nifa', &
                                    ' for merra2 aerosol-aware version of Thompson microphysics'
                  stop
               end if
            else if (.not.is_aerosol_aware .and. .not.merra2_aerosol_aware .and. &
                     (present(nwfa) .or. present(nifa) .or. present(nwfa2d) .or. present(nifa2d))) then
               write(*,*) 'WARNING, nc/nwfa/nifa/nwfa2d/nifa2d present but is_aerosol_aware/merra2_aerosol_aware are FALSE'
            end if
         end if test_only_once

         ! These must be alwyas allocated
         !allocate (vtsk1(kts:kte))
         !allocate (txri1(kts:kte))
         !allocate (txrc1(kts:kte))
         allocate_extended_diagnostics: if (ext_diag) then
            allocate (prw_vcdc1(kts:kte))
            allocate (prw_vcde1(kts:kte))
            allocate (tpri_inu1(kts:kte))
            allocate (tpri_ide1_d(kts:kte))
            allocate (tpri_ide1_s(kts:kte))
            allocate (tprs_ide1(kts:kte))
            allocate (tprs_sde1_d(kts:kte))
            allocate (tprs_sde1_s(kts:kte))
            allocate (tprg_gde1_d(kts:kte))
            allocate (tprg_gde1_s(kts:kte))
            allocate (tpri_iha1(kts:kte))
            allocate (tpri_wfz1(kts:kte))
            allocate (tpri_rfz1(kts:kte))
            allocate (tprg_rfz1(kts:kte))
            allocate (tprs_scw1(kts:kte))
            allocate (tprg_scw1(kts:kte))
            allocate (tprg_rcs1(kts:kte))
            allocate (tprs_rcs1(kts:kte))
            allocate (tprr_rci1(kts:kte))
            allocate (tprg_rcg1(kts:kte))
            allocate (tprw_vcd1_c(kts:kte))
            allocate (tprw_vcd1_e(kts:kte))
            allocate (tprr_sml1(kts:kte))
            allocate (tprr_gml1(kts:kte))
            allocate (tprr_rcg1(kts:kte))
            allocate (tprr_rcs1(kts:kte))
            allocate (tprv_rev1(kts:kte))
            allocate (tten1(kts:kte))
            allocate (qvten1(kts:kte))
            allocate (qrten1(kts:kte))
            allocate (qsten1(kts:kte))
            allocate (qgten1(kts:kte))
            allocate (qiten1(kts:kte))
            allocate (niten1(kts:kte))
            allocate (nrten1(kts:kte))
            allocate (ncten1(kts:kte))
            allocate (qcten1(kts:kte))
         end if allocate_extended_diagnostics

!+---+
         i_start = its
         j_start = jts
         i_end   = ite
         j_end   = jte

!..For idealized testing by developer.
!     if ( (ide-ids+1).gt.4 .and. (jde-jds+1).lt.4 .and.                &
!          ids.eq.its.and.ide.eq.ite.and.jds.eq.jts.and.jde.eq.jte) then
!        i_start = its + 2
!        i_end   = ite
!        j_start = jts
!        j_end   = jte
!     endif

!     dt = dt_in
         RAINNC(:,:) = 0.0
         SNOWNC(:,:) = 0.0
         ICENC(:,:) = 0.0
         GRAUPELNC(:,:) = 0.0
         pcp_ra(:,:) = 0.0
         pcp_sn(:,:) = 0.0
         pcp_gr(:,:) = 0.0
         pcp_ic(:,:) = 0.0
         pfils(:,:,:) = 0.0
         pflls(:,:,:) = 0.0
         rand_pert_max = 0.0
         ndt = max(nint(dt_in/dt_inner),1)
         dt = dt_in/ndt
         if(dt_in .le. dt_inner) dt= dt_in

      !Get the Thompson MP SPP magnitude and standard deviation cutoff,
      !then compute rand_pert_max

         if (rand_perturb_on .ne. 0) then
         do k =1,n_var_spp
            select case (spp_var_list(k))
            case('mp')
               rand_pert_max = spp_prt_list(k)*spp_stddev_cutoff(k)
            end select
         enddo
         endif

      do it = 1, ndt

         qc_max = 0.
         qr_max = 0.
         qs_max = 0.
         qi_max = 0.
         qg_max = 0
         ni_max = 0.
         nr_max = 0.
         imax_qc = 0
         imax_qr = 0
         imax_qi = 0
         imax_qs = 0
         imax_qg = 0
         imax_ni = 0
         imax_nr = 0
         jmax_qc = 0
         jmax_qr = 0
         jmax_qi = 0
         jmax_qs = 0
         jmax_qg = 0
         jmax_ni = 0
         jmax_nr = 0
         kmax_qc = 0
         kmax_qr = 0
         kmax_qi = 0
         kmax_qs = 0
         kmax_qg = 0
         kmax_ni = 0
         kmax_nr = 0

         j_loop:  do j = j_start, j_end
            i_loop:  do i = i_start, i_end

!+---+-----------------------------------------------------------------+
!..Introduce stochastic parameter perturbations by creating as many scalar rand1, rand2, ...
!.. variables as needed to perturb different pieces of microphysics. gthompsn  21Mar2018
! Setting spp_mp_opt to 1 gives graupel Y-intercept pertubations (2^0)
!                   2 gives cloud water distribution gamma shape parameter perturbations (2^1)
!                   4 gives CCN & IN activation perturbations (2^2)
!                   3 gives both 1+2
!                   5 gives both 1+4
!                   6 gives both 2+4
!                   7 gives all 1+2+4
! For now (22Mar2018), standard deviation should be up to 0.75 and cut-off at 3.0
! stddev in order to constrain the various perturbations from being too extreme.
!+---+-----------------------------------------------------------------+
               rand1 = 0.0
               rand2 = 0.0
               rand3 = 0.0
               if (rand_perturb_on .ne. 0) then
                  if (MOD(rand_perturb_on,2) .ne. 0) rand1 = rand_pert(i,1)
                  m = RSHIFT(ABS(rand_perturb_on),1)
                  if (MOD(m,2) .ne. 0) rand2 = rand_pert(i,1)*2.
                  m = RSHIFT(ABS(rand_perturb_on),2)
                  if (MOD(m,2) .ne. 0) rand3 = 0.25*(rand_pert(i,1)+rand_pert_max)
                  m = RSHIFT(ABS(rand_perturb_on),3)
               endif
      !+---+-----------------------------------------------------------------+

               pptrain = 0.
               pptsnow = 0.
               pptgraul = 0.
               pptice = 0.
               RAINNCV(i,j) = 0.
               IF ( PRESENT (snowncv) ) THEN
                  SNOWNCV(i,j) = 0.
               ENDIF
               IF ( PRESENT (icencv) ) THEN
                  ICENCV(i,j) = 0.
               ENDIF
               IF ( PRESENT (graupelncv) ) THEN
                  GRAUPELNCV(i,j) = 0.
               ENDIF
               SR(i,j) = 0.

               do k = kts, kte
                  if (present(tt)) then
                     t1d(k) = tt(i,k,j)
                  else
                     t1d(k) = th(i,k,j)*pii(i,k,j)
                  end if
                  p1d(k) = p(i,k,j)
                  w1d(k) = w(i,k,j)
                  dz1d(k) = dz(i,k,j)
                  qv1d(k) = qv(i,k,j)
                  qc1d(k) = qc(i,k,j)
                  qi1d(k) = qi(i,k,j)
                  qr1d(k) = qr(i,k,j)
                  qs1d(k) = qs(i,k,j)
                  qg1d(k) = qg(i,k,j)
                  ni1d(k) = ni(i,k,j)
                  nr1d(k) = nr(i,k,j)
                  rho(k) = RoverRv*p1d(k) / (R*t1d(k)*(qv1d(k)+RoverRv))

            ! These arrays are always allocated and must be initialized
            !vtsk1(k) = 0.
            !txrc1(k) = 0.
            !txri1(k) = 0.
                  initialize_extended_diagnostics: if (ext_diag) then
                     prw_vcdc1(k) = 0.
                     prw_vcde1(k) = 0.
                     tpri_inu1(k) = 0.
                     tpri_ide1_d(k) = 0.
                     tpri_ide1_s(k) = 0.
                     tprs_ide1(k) = 0.
                     tprs_sde1_d(k) = 0.
                     tprs_sde1_s(k) = 0.
                     tprg_gde1_d(k) = 0.
                     tprg_gde1_s(k) = 0.
                     tpri_iha1(k) = 0.
                     tpri_wfz1(k) = 0.
                     tpri_rfz1(k) = 0.
                     tprg_rfz1(k) = 0.
                     tprs_scw1(k) = 0.
                     tprg_scw1(k) = 0.
                     tprg_rcs1(k) = 0.
                     tprs_rcs1(k) = 0.
                     tprr_rci1(k) = 0.
                     tprg_rcg1(k) = 0.
                     tprw_vcd1_c(k) = 0.
                     tprw_vcd1_e(k) = 0.
                     tprr_sml1(k) = 0.
                     tprr_gml1(k) = 0.
                     tprr_rcg1(k) = 0.
                     tprr_rcs1(k) = 0.
                     tprv_rev1(k) = 0.
                     tten1(k) = 0.
                     qvten1(k) = 0.
                     qrten1(k) = 0.
                     qsten1(k) = 0.
                     qgten1(k) = 0.
                     qiten1(k) = 0.
                     niten1(k) = 0.
                     nrten1(k) = 0.
                     ncten1(k) = 0.
                     qcten1(k) = 0.
                  endif initialize_extended_diagnostics
               enddo

               lsml = lsm(i,j)
               if (is_aerosol_aware .or. merra2_aerosol_aware) then
                  do k = kts, kte
                     nc1d(k) = nc(i,k,j)
                     nwfa1d(k) = nwfa(i,k,j)
                     nifa1d(k) = nifa(i,k,j)
                  enddo
               else
                  do k = kts, kte
                     if(lsml == 1) then
                        nc1d(k) = Nt_c_l/rho(k)
                     else
                        nc1d(k) = Nt_c_o/rho(k)
                     endif
                     nwfa1d(k) = 11.1E6
                     nifa1d(k) = naIN1*0.01
                  enddo
               endif

!..If not the variable-density graupel-hail hybrid, then set the vol mixing
!.. ratio to mass mixing ratio divided by constant density (500kg/m3) value.

               if (is_hail_aware) then
                  do k = kts, kte
                      ng1d(k) = ng(i,k,j)
                      qb1d(k) = qb(i,k,j)
                  enddo
              else

                  do k = kte, kts, -1
                      if (qg1d(k).gt.R1) then
                          ygra1 = alog10(max(1.E-9, qg1d(k)*rho(k)))
                          zans1 = 3.0 + 2./7.*(ygra1+8.)
                          zans1 = MAX(2., MIN(zans1, 6.))
                          N0_exp = 10.**(zans1)
                          lam_exp = (N0_exp*am_g(idx_bg1)*cgg(1,1)/(rho(k)*qg1d(k)))**oge1
                          lamg = lam_exp * (cgg(3,1)*ogg2*ogg1)**obmg
                          ng1d(k) = cgg(2,1)*ogg3*rho(k)*qg1d(k)*lamg**bm_g / am_g(idx_bg1)
                          ng1d(k) = MAX(R2, ng1d(k)/rho(k))
                          qb1d(k) = qg1d(k)/rho_g(idx_bg1)
                      else
                          ng1d(k) = 0
                          qb1d(k) = 0
                      endif
                  enddo
              endif

!> - Call mp_thompson()
               call mp_thompson(qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, qb1d, ni1d,     &
                           nr1d, nc1d, ng1d, nwfa1d, nifa1d, t1d, p1d, w1d, dz1d,   &
                           pptrain, pptsnow, pptgraul, pptice, &
#if ( WRF_CHEM == 1 )
                     rainprod1d, evapprod1d, &
#endif
                           rand1, rand2, rand3, &
                           ext_diag,                                        & 
                           sedi_semi, decfl,                                &
                           !vtsk1, txri1, txrc1,                            &
                           prw_vcdc1, prw_vcde1,                            &
                           tpri_inu1, tpri_ide1_d, tpri_ide1_s, tprs_ide1,  &
                           tprs_sde1_d, tprs_sde1_s,                        &
                           tprg_gde1_d, tprg_gde1_s, tpri_iha1, tpri_wfz1,  &
                           tpri_rfz1, tprg_rfz1, tprs_scw1, tprg_scw1,      &
                           tprg_rcs1, tprs_rcs1, tprr_rci1,                 &
                           tprg_rcg1, tprw_vcd1_c,                          &
                           tprw_vcd1_e, tprr_sml1, tprr_gml1, tprr_rcg1,    &
                           tprr_rcs1, tprv_rev1,                            &
                           tten1, qvten1, qrten1, qsten1,                   &
                           qgten1, qiten1, niten1, nrten1, ncten1, qcten1,  &
                           pfil1, pfll1, lsml,                              &
                           kts, kte, dt, i, j)

               pcp_ra(i,j) = pcp_ra(i,j) + pptrain
               pcp_sn(i,j) = pcp_sn(i,j) + pptsnow
               pcp_gr(i,j) = pcp_gr(i,j) + pptgraul
               pcp_ic(i,j) = pcp_ic(i,j) + pptice
               RAINNCV(i,j) = pptrain + pptsnow + pptgraul + pptice
               RAINNC(i,j) = RAINNC(i,j) + pptrain + pptsnow + pptgraul + pptice
               IF ( PRESENT(snowncv) .AND. PRESENT(snownc) ) THEN
                  ! Add ice to snow if separate ice not present
                  IF ( .NOT.PRESENT(icencv) .OR. .NOT.PRESENT(icenc) ) THEN
                     SNOWNCV(i,j) = pptsnow + pptice
                     SNOWNC(i,j) = SNOWNC(i,j) + pptsnow + pptice
                  ELSE
                     SNOWNCV(i,j) = pptsnow
                     SNOWNC(i,j) = SNOWNC(i,j) + pptsnow
                  ENDIF
               ENDIF
               ! Use separate ice if present (as in FV3)
               IF ( PRESENT(icencv) .AND. PRESENT(icenc) ) THEN
                  ICENCV(i,j) = pptice
                  ICENC(i,j) = ICENC(i,j) + pptice
               ENDIF
               IF ( PRESENT(graupelncv) .AND. PRESENT(graupelnc) ) THEN
                  GRAUPELNCV(i,j) = pptgraul
                  GRAUPELNC(i,j) = GRAUPELNC(i,j) + pptgraul
               ENDIF
               SR(i,j) = (pptsnow + pptgraul + pptice) / (RAINNCV(i,j)+R1)

!..Reset lowest model level to initial state aerosols (fake sfc source).
!.. Changed 13 May 2013 to fake emissions in which nwfa2d is aerosol
!.. number tendency (number per kg per second).
               if (is_aerosol_aware) then
                  if ( PRESENT (aero_ind_fdb) ) then
                  if ( .not. aero_ind_fdb) then
                     nwfa1d(kts) = nwfa1d(kts) + nwfa2d(i,j)*dt
                     nifa1d(kts) = nifa1d(kts) + nifa2d(i,j)*dt
                  endif
                  else
                  nwfa1d(kts) = nwfa1d(kts) + nwfa2d(i,j)*dt
                  nifa1d(kts) = nifa1d(kts) + nifa2d(i,j)*dt
                  end if

                  do k = kts, kte
                     nc(i,k,j) = nc1d(k)
                     nwfa(i,k,j) = nwfa1d(k)
                     nifa(i,k,j) = nifa1d(k)
                  enddo
               endif

               if (merra2_aerosol_aware) then
                  do k = kts, kte
                     nc(i,k,j) = nc1d(k)
                     nwfa(i,k,j) = nwfa1d(k)
                     nifa(i,k,j) = nifa1d(k)
                  enddo
               endif
              
               if (is_hail_aware) then
                  do k = kts, kte
                      ng(i,k,j) = ng1d(k)
                      qb(i,k,j) = qb1d(k)
                  enddo
              endif

               do k = kts, kte
                  qv(i,k,j) = qv1d(k)
                  qc(i,k,j) = qc1d(k)
                  qi(i,k,j) = qi1d(k)
                  qr(i,k,j) = qr1d(k)
                  qs(i,k,j) = qs1d(k)
                  qg(i,k,j) = qg1d(k)
                  ni(i,k,j) = ni1d(k)
                  nr(i,k,j) = nr1d(k)
                  pfils(i,k,j) = pfils(i,k,j) + pfil1(k)
                  pflls(i,k,j) = pflls(i,k,j) + pfll1(k)
                  if (present(tt)) then
                     tt(i,k,j) = t1d(k)
                  else
                     th(i,k,j) = t1d(k)/pii(i,k,j)
                  endif
#if ( WRF_CHEM == 1 )
            rainprod(i,k,j) = rainprod1d(k)
            evapprod(i,k,j) = evapprod1d(k)
#endif
                  if (qc1d(k) .gt. qc_max) then
                     imax_qc = i
                     jmax_qc = j
                     kmax_qc = k
                     qc_max = qc1d(k)
                  elseif (qc1d(k) .lt. 0.0) then
                     write(*,'(a,e16.7,a,3i8)') 'WARNING, negative qc ', qc1d(k),        &
                                 ' at i,j,k=', i,j,k
                  endif
                  if (qr1d(k) .gt. qr_max) then
                     imax_qr = i
                     jmax_qr = j
                     kmax_qr = k
                     qr_max = qr1d(k)
                  elseif (qr1d(k) .lt. 0.0) then
                     write(*,'(a,e16.7,a,3i8)') 'WARNING, negative qr ', qr1d(k),        &
                                 ' at i,j,k=', i,j,k
                  endif
                  if (nr1d(k) .gt. nr_max) then
                     imax_nr = i
                     jmax_nr = j
                     kmax_nr = k
                     nr_max = nr1d(k)
                  elseif (nr1d(k) .lt. 0.0) then
                     write(*,'(a,e16.7,a,3i8)') 'WARNING, negative nr ', nr1d(k),        &
                                 ' at i,j,k=', i,j,k
                  endif
                  if (qs1d(k) .gt. qs_max) then
                     imax_qs = i
                     jmax_qs = j
                     kmax_qs = k
                     qs_max = qs1d(k)
                  elseif (qs1d(k) .lt. 0.0) then
                     write(*,'(a,e16.7,a,3i8)') 'WARNING, negative qs ', qs1d(k),        &
                                 ' at i,j,k=', i,j,k
                  endif
                  if (qi1d(k) .gt. qi_max) then
                     imax_qi = i
                     jmax_qi = j
                     kmax_qi = k
                     qi_max = qi1d(k)
                  elseif (qi1d(k) .lt. 0.0) then
                     write(*,'(a,e16.7,a,3i8)') 'WARNING, negative qi ', qi1d(k),        &
                                 ' at i,j,k=', i,j,k
                  endif
                  if (qg1d(k) .gt. qg_max) then
                     imax_qg = i
                     jmax_qg = j
                     kmax_qg = k
                     qg_max = qg1d(k)
                  elseif (qg1d(k) .lt. 0.0) then
                     write(*,'(a,e16.7,a,3i8)') 'WARNING, negative qg ', qg1d(k),        &
                                 ' at i,j,k=', i,j,k
                  endif
                  if (ni1d(k) .gt. ni_max) then
                     imax_ni = i
                     jmax_ni = j
                     kmax_ni = k
                     ni_max = ni1d(k)
                  elseif (ni1d(k) .lt. 0.0) then
                     write(*,'(a,e16.7,a,3i8)') 'WARNING, negative ni ', ni1d(k),        &
                                 ' at i,j,k=', i,j,k
                  endif
                  if (qv1d(k) .lt. 0.0) then
                     write(*,'(a,e16.7,a,3i8)') 'WARNING, negative qv ', qv1d(k),        &
                                 ' at i,j,k=', i,j,k
                     if (k.lt.kte-2 .and. k.gt.kts+1) then
                        write(*,*) '   below and above are: ', qv(i,k-1,j), qv(i,k+1,j)
                        qv(i,k,j) = max(1.E-7, 0.5*(qv(i,k-1,j) + qv(i,k+1,j)))
                     else
                        qv(i,k,j) = 1.E-7
                     endif
                  endif
               enddo

               assign_extended_diagnostics: if (ext_diag) then
                  do k=kts,kte
                     !vts1(i,k,j)       = vtsk1(k)
                     !txri(i,k,j)       = txri(i,k,j)       + txri1(k)
                     !txrc(i,k,j)       = txrc(i,k,j)       + txrc1(k)
                     prw_vcdc(i,k,j)   = prw_vcdc(i,k,j)   + prw_vcdc1(k)
                     prw_vcde(i,k,j)   = prw_vcde(i,k,j)   + prw_vcde1(k)
                     tpri_inu(i,k,j)   = tpri_inu(i,k,j)   + tpri_inu1(k) 
                     tpri_ide_d(i,k,j) = tpri_ide_d(i,k,j) + tpri_ide1_d(k)
                     tpri_ide_s(i,k,j) = tpri_ide_s(i,k,j) + tpri_ide1_s(k)
                     tprs_ide(i,k,j)   = tprs_ide(i,k,j)   + tprs_ide1(k)
                     tprs_sde_s(i,k,j) = tprs_sde_s(i,k,j) + tprs_sde1_s(k)
                     tprs_sde_d(i,k,j) = tprs_sde_d(i,k,j) + tprs_sde1_d(k)
                     tprg_gde_d(i,k,j) = tprg_gde_d(i,k,j) + tprg_gde1_d(k)
                     tprg_gde_s(i,k,j) = tprg_gde_s(i,k,j) + tprg_gde1_s(k)
                     tpri_iha(i,k,j)   = tpri_iha(i,k,j)   + tpri_iha1(k)
                     tpri_wfz(i,k,j)   = tpri_wfz(i,k,j)   + tpri_wfz1(k)
                     tpri_rfz(i,k,j)   = tpri_rfz(i,k,j)   + tpri_rfz1(k)
                     tprg_rfz(i,k,j)   = tprg_rfz(i,k,j)   + tprg_rfz1(k)
                     tprs_scw(i,k,j)   = tprs_scw(i,k,j)   + tprs_scw1(k)
                     tprg_scw(i,k,j)   = tprg_scw(i,k,j)   + tprg_scw1(k)
                     tprg_rcs(i,k,j)   = tprg_rcs(i,k,j)   + tprg_rcs1(k)
                     tprs_rcs(i,k,j)   = tprs_rcs(i,k,j)   + tprs_rcs1(k)
                     tprr_rci(i,k,j)   = tprr_rci(i,k,j)   + tprr_rci1(k)
                     tprg_rcg(i,k,j)   = tprg_rcg(i,k,j)   + tprg_rcg1(k)
                     tprw_vcd_c(i,k,j) = tprw_vcd_c(i,k,j) + tprw_vcd1_c(k)
                     tprw_vcd_e(i,k,j) = tprw_vcd_e(i,k,j) + tprw_vcd1_e(k)
                     tprr_sml(i,k,j)   = tprr_sml(i,k,j)   + tprr_sml1(k)
                     tprr_gml(i,k,j)   = tprr_gml(i,k,j)   + tprr_gml1(k)
                     tprr_rcg(i,k,j)   = tprr_rcg(i,k,j)   + tprr_rcg1(k)
                     tprr_rcs(i,k,j)   = tprr_rcs(i,k,j)   + tprr_rcs1(k)
                     tprv_rev(i,k,j)   = tprv_rev(i,k,j)   + tprv_rev1(k)
                     tten3(i,k,j)      = tten3(i,k,j)      + tten1(k) 
                     qvten3(i,k,j)     = qvten3(i,k,j)     + qvten1(k)
                     qrten3(i,k,j)     = qrten3(i,k,j)     + qrten1(k)
                     qsten3(i,k,j)     = qsten3(i,k,j)     + qsten1(k)
                     qgten3(i,k,j)     = qgten3(i,k,j)     + qgten1(k)
                     qiten3(i,k,j)     = qiten3(i,k,j)     + qiten1(k) 
                     niten3(i,k,j)     = niten3(i,k,j)     + niten1(k)
                     nrten3(i,k,j)     = nrten3(i,k,j)     + nrten1(k)
                     ncten3(i,k,j)     = ncten3(i,k,j)     + ncten1(k)
                     qcten3(i,k,j)     = qcten3(i,k,j)     + qcten1(k)
                  enddo
               endif assign_extended_diagnostics

               if (ndt>1 .and. it==ndt) then
                  SR(i,j) = (pcp_sn(i,j) + pcp_gr(i,j) + pcp_ic(i,j)) / (RAINNC(i,j)+R1)
                  RAINNCV(i,j) = RAINNC(i,j)
                  IF ( PRESENT (snowncv) ) THEN
                     SNOWNCV(i,j) = SNOWNC(i,j)
                  ENDIF
                  IF ( PRESENT (icencv) ) THEN
                     ICENCV(i,j) = ICENC(i,j)
                  ENDIF
                  IF ( PRESENT (graupelncv) ) THEN
                     GRAUPELNCV(i,j) = GRAUPELNC(i,j)
                  ENDIF
               endif 

         ! Diagnostic calculations only for last step
         ! if Thompson MP is called multiple times
               last_step_only: IF ((ndt>1 .and. it==ndt) .or. &
                                 (nsteps>1 .and. istep==nsteps) .or. &
                                 (nsteps==1 .and. ndt==1)) THEN

                  ! max_hail_diam_sfc(i,j) = hail_mass_99th_percentile(kts, kte, qg1d, t1d, p1d, qv1d)

!> - Call calc_refl10cm()

                  diagflag_present: IF ( PRESENT (diagflag) ) THEN
                     if (diagflag .and. do_radar_ref == 1) then
            !
                        ! Only set melti to true at the output times
                        if (fullradar_diag) then
                           melti=.true.
                        else
                           melti=.false.
                        endif
            !
!                        if (present(vt_dbz_wt)) then
!                           call calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d,   &
!                                             t1d, p1d, dBZ, rand1, kts, kte, i, j, &
!                                             melti, vt_dbz_wt(i,:,j),              &
!                                             first_time_step)
!                        else
!AAJ                           call calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d,   &
!AAJ                                             t1d, p1d, dBZ, rand1, kts, kte, i, j, &
!AAJ                                             melti)
                           call calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d,   &
                                             t1d, p1d, dBZ, kts, kte, i, j)
                                         
!                        endif
                        do k = kts, kte
                           refl_10cm(i,k,j) = max(-35., dBZ(k))
                        enddo
                     endif
                  ENDIF diagflag_present

                  IF (has_reqc.ne.0 .and. has_reqi.ne.0 .and. has_reqs.ne.0) THEN
                     do k = kts, kte
                        re_qc1d(k) = re_qc_min
                        re_qi1d(k) = re_qi_min
                        re_qs1d(k) = re_qs_min
                     enddo
         !> - Call calc_effectrad()
!                     call calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d,  &
!                          re_qc1d, re_qi1d, re_qs1d, lsml, kts, kte)
                     call calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d,  &
                                          re_qc1d, re_qi1d, re_qs1d, kts, kte)
                     do k = kts, kte
                        re_cloud(i,k,j) = max(re_qc_min, min(re_qc1d(k), re_qc_max))
                        re_ice(i,k,j)   = max(re_qi_min, min(re_qi1d(k), re_qi_max))
                        re_snow(i,k,j)  = max(re_qs_min, min(re_qs1d(k), re_qs_max))
                     enddo
                  ENDIF
               ENDIF last_step_only
            enddo i_loop
         enddo j_loop

! DEBUG - GT
!      write(*,'(a,7(a,e13.6,1x,a,i3,a,i3,a,i3,a,1x))') 'MP-GT:', &
!         'qc: ', qc_max, '(', imax_qc, ',', jmax_qc, ',', kmax_qc, ')', &
!         'qr: ', qr_max, '(', imax_qr, ',', jmax_qr, ',', kmax_qr, ')', &
!         'qi: ', qi_max, '(', imax_qi, ',', jmax_qi, ',', kmax_qi, ')', &
!         'qs: ', qs_max, '(', imax_qs, ',', jmax_qs, ',', kmax_qs, ')', &
!         'qg: ', qg_max, '(', imax_qg, ',', jmax_qg, ',', kmax_qg, ')', &
!         'ni: ', ni_max, '(', imax_ni, ',', jmax_ni, ',', kmax_ni, ')', &
!         'nr: ', nr_max, '(', imax_nr, ',', jmax_nr, ',', kmax_nr, ')'
! END DEBUG - GT
      enddo ! end of nt loop

      do j = j_start, j_end
        do k = kts, kte
          do i = i_start, i_end
            pfils(i,k,j) = pfils(i,k,j)/dt_in
            pflls(i,k,j) = pflls(i,k,j)/dt_in
          enddo
        enddo 
      enddo

      ! These are always allocated
      !deallocate (vtsk1)
      !deallocate (txri1)
      !deallocate (txrc1)
      deallocate_extended_diagnostics: if (ext_diag) then
         deallocate (prw_vcdc1)
         deallocate (prw_vcde1)
         deallocate (tpri_inu1)
         deallocate (tpri_ide1_d)
         deallocate (tpri_ide1_s)
         deallocate (tprs_ide1)
         deallocate (tprs_sde1_d)
         deallocate (tprs_sde1_s)
         deallocate (tprg_gde1_d)
         deallocate (tprg_gde1_s)
         deallocate (tpri_iha1)
         deallocate (tpri_wfz1)
         deallocate (tpri_rfz1)
         deallocate (tprg_rfz1)
         deallocate (tprs_scw1)
         deallocate (tprg_scw1)
         deallocate (tprg_rcs1)
         deallocate (tprs_rcs1)
         deallocate (tprr_rci1)
         deallocate (tprg_rcg1)
         deallocate (tprw_vcd1_c)
         deallocate (tprw_vcd1_e)
         deallocate (tprr_sml1)
         deallocate (tprr_gml1)
         deallocate (tprr_rcg1)
         deallocate (tprr_rcs1)
         deallocate (tprv_rev1)
         deallocate (tten1)
         deallocate (qvten1)
         deallocate (qrten1)
         deallocate (qsten1)
         deallocate (qgten1)
         deallocate (qiten1)
         deallocate (niten1)
         deallocate (nrten1)
         deallocate (ncten1)
         deallocate (qcten1)
      end if deallocate_extended_diagnostics

   end subroutine mp_gt_driver
!> @}

!>\ingroup aathompson
   subroutine thompson_finalize()

      implicit none

      if (ALLOCATED(tcg_racg)) DEALLOCATE(tcg_racg)
      if (ALLOCATED(tmr_racg)) DEALLOCATE(tmr_racg)
      if (ALLOCATED(tcr_gacr)) DEALLOCATE(tcr_gacr)
      ! if (ALLOCATED(tmg_gacr)) DEALLOCATE(tmg_gacr)
      if (ALLOCATED(tnr_racg)) DEALLOCATE(tnr_racg)
      if (ALLOCATED(tnr_gacr)) DEALLOCATE(tnr_gacr)

      if (ALLOCATED(tcs_racs1)) DEALLOCATE(tcs_racs1)
      if (ALLOCATED(tmr_racs1)) DEALLOCATE(tmr_racs1)
      if (ALLOCATED(tcs_racs2)) DEALLOCATE(tcs_racs2)
      if (ALLOCATED(tmr_racs2)) DEALLOCATE(tmr_racs2)
      if (ALLOCATED(tcr_sacr1)) DEALLOCATE(tcr_sacr1)
      if (ALLOCATED(tms_sacr1)) DEALLOCATE(tms_sacr1)
      if (ALLOCATED(tcr_sacr2)) DEALLOCATE(tcr_sacr2)
      if (ALLOCATED(tms_sacr2)) DEALLOCATE(tms_sacr2)
      if (ALLOCATED(tnr_racs1)) DEALLOCATE(tnr_racs1)
      if (ALLOCATED(tnr_racs2)) DEALLOCATE(tnr_racs2)
      if (ALLOCATED(tnr_sacr1)) DEALLOCATE(tnr_sacr1)
      if (ALLOCATED(tnr_sacr2)) DEALLOCATE(tnr_sacr2)

      if (ALLOCATED(tpi_qcfz)) DEALLOCATE(tpi_qcfz)
      if (ALLOCATED(tni_qcfz)) DEALLOCATE(tni_qcfz)

      if (ALLOCATED(tpi_qrfz)) DEALLOCATE(tpi_qrfz)
      if (ALLOCATED(tpg_qrfz)) DEALLOCATE(tpg_qrfz)
      if (ALLOCATED(tni_qrfz)) DEALLOCATE(tni_qrfz)
      if (ALLOCATED(tnr_qrfz)) DEALLOCATE(tnr_qrfz)

      if (ALLOCATED(tps_iaus)) DEALLOCATE(tps_iaus)
      if (ALLOCATED(tni_iaus)) DEALLOCATE(tni_iaus)
      if (ALLOCATED(tpi_ide))  DEALLOCATE(tpi_ide)

      if (ALLOCATED(t_Efrw)) DEALLOCATE(t_Efrw)
      if (ALLOCATED(t_Efsw)) DEALLOCATE(t_Efsw)

      if (ALLOCATED(tnr_rev)) DEALLOCATE(tnr_rev)
      if (ALLOCATED(tpc_wev)) DEALLOCATE(tpc_wev)
      if (ALLOCATED(tnc_wev)) DEALLOCATE(tnc_wev)

      if (ALLOCATED(tnccn_act)) DEALLOCATE(tnccn_act)

   end subroutine thompson_finalize

!+---+-----------------------------------------------------------------+

! !+---+-----------------------------------------------------------------+
! !ctrlL
! !+---+-----------------------------------------------------------------+
! !..Creation of the lookup tables and support functions found below here.
! !+---+-----------------------------------------------------------------+
! !>\ingroup aathompson
!! Rain collecting graupel (and inverse).  Explicit CE integration.
   subroutine qr_acr_qg(NRHGtable)

      implicit none

      INTEGER, INTENT(IN) ::NRHGtable

!..Local variables
      integer:: i, j, k, m, n, n2, n3, idx_bg
      integer:: km, km_s, km_e
      real(dp), dimension(nbg):: N_g
      real(dp), dimension(nbg,NRHGtable):: vg
      real(dp), dimension(nbr):: vr, N_r
      real(dp) :: N0_r, N0_g, lam_exp, lamg, lamr
      real(dp) :: massg, massr, dvg, dvr, t1, t2, z1, z2, y1, y2
      logical :: force_read_thompson, write_thompson_tables
      logical :: lexist,lopen
      integer :: good,ierr

      force_read_thompson = .false.
      write_thompson_tables = .false.
!+---+


      good = 0
        INQUIRE(FILE=qr_acr_qg_file, EXIST=lexist)
#ifdef MPI
        call MPI_BARRIER(mpi_communicator,ierr)
#endif
        IF ( lexist ) THEN
          OPEN(63,file=qr_acr_qg_file,form="unformatted",err=1234)
!sms$serial begin
          READ(63,err=1234) tcg_racg
          READ(63,err=1234) tmr_racg
          READ(63,err=1234) tcr_gacr
          ! READ(63,err=1234) tmg_gacr
          READ(63,err=1234) tnr_racg
          READ(63,err=1234) tnr_gacr
!sms$serial end
          good = 1
 1234     CONTINUE
          IF ( good .NE. 1 ) THEN
            INQUIRE(63,opened=lopen)
            IF (lopen) THEN
              IF( force_read_thompson ) THEN
                write(0,*) "Error reading "//qr_acr_qg_file//" Aborting because force_read_thompson is .true."
                return
              ENDIF
              CLOSE(63)
            ELSE
              IF( force_read_thompson ) THEN
                write(0,*) "Error opening "//qr_acr_qg_file//" Aborting because force_read_thompson is .true."
                return
              ENDIF
            ENDIF
         ELSE
            INQUIRE(63,opened=lopen)
            IF (lopen) THEN
              CLOSE(63)
            ENDIF
          ENDIF
        ELSE
          IF( force_read_thompson ) THEN
            write(0,*) "Non-existent "//qr_acr_qg_file//" Aborting because force_read_thompson is .true."
            return
          ENDIF
        ENDIF

      IF (.NOT. good .EQ. 1 ) THEN
        if (thompson_table_writer) then
          write_thompson_tables = .true.
          write(0,*) "ThompMP: computing qr_acr_qg"
        endif
        do n2 = 1, nbr
!        vr(n2) = av_r*Dr(n2)**bv_r * exp(real(-fv_r*Dr(n2), kind=dp))
         vr(n2) = -0.1021 + 4.932E3*Dr(n2) - 0.9551E6*Dr(n2)*Dr(n2)     &
              + 0.07934E9*Dr(n2)*Dr(n2)*Dr(n2)                          &
              - 0.002362E12*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
        enddo
        do n3 = 1, NRHGtable
           do n = 1, nbg
             ! idx_bg indexes module coefficients, not vg
             !!         if (.not. is_hail_aware) idx_bg = idx_bg1
             idx_bg = 6
             !!         if (is_hail_aware) idx_bg = n3
             vg(n,n3) = av_g(idx_bg)*Dg(n)**bv_g(idx_bg)
           enddo
        enddo

!..Note values returned from wrf_dm_decomp1d are zero-based, add 1 for
!.. fortran indices.  J. Michalakes, 2009Oct30.

#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
        CALL wrf_dm_decomp1d ( ntb_r*ntb_r1, km_s, km_e )
#else
        km_s = 0
        km_e = ntb_r*ntb_r1 - 1
#endif

        do km = km_s, km_e
         m = km / ntb_r1 + 1
         k = mod( km , ntb_r1 ) + 1

         lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
         lamr = lam_exp * (crg(3)*org2*org1)**obmr
         N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
         do n2 = 1, nbr
            N_r(n2) = N0_r*Dr(n2)**mu_r *exp(real(-lamr*Dr(n2), kind=dp))*dtr(n2)
         enddo
         do n3 = 1, NRHGtable
            idx_bg = 6
            !!          if (.not. is_hail_aware) idx_bg = idx_bg1
            !!          if (is_hail_aware) idx_bg = n3

            do j = 1, ntb_g
               do i = 1, ntb_g1
                   lam_exp = (N0g_exp(i)*am_g(idx_bg)*cgg(1,1)/r_g(j))**oge1
                   lamg = lam_exp * (cgg(3,1)*ogg2*ogg1)**obmg
                   N0_g = N0g_exp(i)/(cgg(2,1)*lam_exp) * lamg**cge(2,1)
                   do n = 1, nbg
                       N_g(n) = N0_g*Dg(n)**mu_g * DEXP(-lamg*Dg(n))*dtg(n)
                   enddo

            t1 = 0.0_dp
            t2 = 0.0_dp
            z1 = 0.0_dp
            z2 = 0.0_dp
            y1 = 0.0_dp
            y2 = 0.0_dp
            do n2 = 1, nbr
               massr = am_r * Dr(n2)**bm_r
               do n = 1, nbg
                  massg = am_g(idx_bg) * Dg(n)**bm_g

                  dvg = 0.5d0*((vr(n2) - vg(n,n3)) + abs(real(vr(n2)-vg(n,n3), kind=dp)))
                  dvr = 0.5d0*((vg(n,n3) - vr(n2)) + abs(real(vg(n,n3)-vr(n2), kind=dp)))

                  t1 = t1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                      *dvg*massg * N_g(n)* N_r(n2)
                  z1 = z1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                      *dvg*massr * N_g(n)* N_r(n2)
                  y1 = y1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                      *dvg       * N_g(n)* N_r(n2)

                  t2 = t2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                      *dvr*massr * N_g(n)* N_r(n2)
                  y2 = y2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                      *dvr       * N_g(n)* N_r(n2)
                  z2 = z2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                      *dvr*massg * N_g(n)* N_r(n2)
               enddo
 97            continue
            enddo
            tcg_racg(i,j,n3,k,m) = t1
            tmr_racg(i,j,n3,k,m) = min(z1, r_r(m)*1.0_dp)
            tcr_gacr(i,j,n3,k,m) = t2
            ! tmg_gacr(i,j,k,m) = min(z2, r_g(j)*1.0_dp)
            tnr_racg(i,j,n3,k,m) = y1
            tnr_gacr(i,j,n3,k,m) = y2
         enddo
         enddo
        enddo
      enddo 

        IF ( write_thompson_tables ) THEN
          write(0,*) "Writing "//qr_acr_qg_file//" in Thompson MP init"
          OPEN(63,file=qr_acr_qg_file,form="unformatted",err=9234)
          WRITE(63,err=9234) tcg_racg
          WRITE(63,err=9234) tmr_racg
          WRITE(63,err=9234) tcr_gacr
          ! WRITE(63,err=9234) tmg_gacr
          WRITE(63,err=9234) tnr_racg
          WRITE(63,err=9234) tnr_gacr
          CLOSE(63)
          RETURN    ! ----- RETURN
 9234     CONTINUE
          write(0,*) "Error writing "//qr_acr_qg_file
          return
        ENDIF
      ENDIF

   end subroutine qr_acr_qg
!+---+-----------------------------------------------------------------+
!ctrlL
!+---+-----------------------------------------------------------------+
!>\ingroup aathompson
!!Rain collecting snow (and inverse).  Explicit CE integration.
   subroutine qr_acr_qs

      implicit none

!..Local variables
      integer:: i, j, k, m, n, n2
      integer:: km, km_s, km_e
      real(dp), dimension(nbr):: vr, D1, N_r
      real(dp), dimension(nbs):: vs, N_s
      real(dp) :: loga_, a_, b_, second, M0, M2, M3, Mrat, oM3
      real(dp) :: N0_r, lam_exp, lamr, slam1, slam2
      real(dp) :: dvs, dvr, masss, massr
      real(dp) :: t1, t2, t3, t4, z1, z2, z3, z4
      real(dp) :: y1, y2, y3, y4
      logical force_read_thompson, write_thompson_tables
      logical lexist,lopen
      integer good,ierr

!+---+

      force_read_thompson = .false.
      write_thompson_tables = .false.

      good = 0
        INQUIRE(FILE=qr_acr_qs_file, EXIST=lexist)
#ifdef MPI
        call MPI_BARRIER(mpi_communicator,ierr)
#endif
        IF ( lexist ) THEN
          !write(0,*) "ThompMP: read "//qr_acr_qs_file//" instead of computing"
          OPEN(63,file=qr_acr_qs_file,form="unformatted",err=1234)
!sms$serial begin
          READ(63,err=1234)tcs_racs1
          READ(63,err=1234)tmr_racs1
          READ(63,err=1234)tcs_racs2
          READ(63,err=1234)tmr_racs2
          READ(63,err=1234)tcr_sacr1
          READ(63,err=1234)tms_sacr1
          READ(63,err=1234)tcr_sacr2
          READ(63,err=1234)tms_sacr2
          READ(63,err=1234)tnr_racs1
          READ(63,err=1234)tnr_racs2
          READ(63,err=1234)tnr_sacr1
          READ(63,err=1234)tnr_sacr2
!sms$serial end
          good = 1
 1234     CONTINUE
          IF ( good .NE. 1 ) THEN
            INQUIRE(63,opened=lopen)
            IF (lopen) THEN
              IF( force_read_thompson ) THEN
                write(0,*) "Error reading "//qr_acr_qs_file//" Aborting because force_read_thompson is .true."
                return
              ENDIF
              CLOSE(63)
            ELSE
              IF( force_read_thompson ) THEN
                write(0,*) "Error opening "//qr_acr_qs_file//" Aborting because force_read_thompson is .true."
                return
              ENDIF
            ENDIF
          ELSE
            INQUIRE(63,opened=lopen)
            IF (lopen) THEN
              CLOSE(63)
            ENDIF
          ENDIF
        ELSE
          IF( force_read_thompson ) THEN
            write(0,*) "Non-existent "//qr_acr_qs_file//" Aborting because force_read_thompson is .true."
            return
          ENDIF
        ENDIF

      IF (.NOT. good .EQ. 1 ) THEN
        if (thompson_table_writer) then
          write_thompson_tables = .true.
          write(0,*) "ThompMP: computing qr_acr_qs"
        endif
        do n2 = 1, nbr
!        vr(n2) = av_r*Dr(n2)**bv_r * exp(real(-fv_r*Dr(n2), kind=dp))
         vr(n2) = -0.1021 + 4.932E3*Dr(n2) - 0.9551E6*Dr(n2)*Dr(n2)     &
              + 0.07934E9*Dr(n2)*Dr(n2)*Dr(n2)                          &
              - 0.002362E12*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
         D1(n2) = (vr(n2)/av_s)**(1./bv_s)
        enddo
        do n = 1, nbs
         vs(n) = 1.5*av_s*Ds(n)**bv_s * exp(real(-fv_s*Ds(n), kind=dp))
        enddo

!..Note values returned from wrf_dm_decomp1d are zero-based, add 1 for
!.. fortran indices.  J. Michalakes, 2009Oct30.

#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
        CALL wrf_dm_decomp1d ( ntb_r*ntb_r1, km_s, km_e )
#else
        km_s = 0
        km_e = ntb_r*ntb_r1 - 1
#endif

        do km = km_s, km_e
         m = km / ntb_r1 + 1
         k = mod( km , ntb_r1 ) + 1

         lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
         lamr = lam_exp * (crg(3)*org2*org1)**obmr
         N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
         do n2 = 1, nbr
            N_r(n2) = N0_r*Dr(n2)**mu_r * exp(real(-lamr*Dr(n2), kind=dp))*dtr(n2)
         enddo

         do j = 1, ntb_t
            do i = 1, ntb_s

!..From the bm_s moment, compute plus one moment.  If we are not
!.. using bm_s=2, then we must transform to the pure 2nd moment
!.. (variable called "second") and then to the bm_s+1 moment.

               M2 = r_s(i)*oams*1.0_dp
               if (bm_s.gt.2.0-1.E-3 .and. bm_s.lt.2.0+1.E-3) then
                  loga_ = sa(1) + sa(2)*Tc(j) + sa(3)*bm_s &
                     + sa(4)*Tc(j)*bm_s + sa(5)*Tc(j)*Tc(j) &
                     + sa(6)*bm_s*bm_s + sa(7)*Tc(j)*Tc(j)*bm_s &
                     + sa(8)*Tc(j)*bm_s*bm_s + sa(9)*Tc(j)*Tc(j)*Tc(j) &
                     + sa(10)*bm_s*bm_s*bm_s
                  a_ = 10.0**loga_
                  b_ = sb(1) + sb(2)*Tc(j) + sb(3)*bm_s &
                     + sb(4)*Tc(j)*bm_s + sb(5)*Tc(j)*Tc(j) &
                     + sb(6)*bm_s*bm_s + sb(7)*Tc(j)*Tc(j)*bm_s &
                     + sb(8)*Tc(j)*bm_s*bm_s + sb(9)*Tc(j)*Tc(j)*Tc(j) &
                     + sb(10)*bm_s*bm_s*bm_s
                  second = (M2/a_)**(1./b_)
               else
                  second = M2
               endif

               loga_ = sa(1) + sa(2)*Tc(j) + sa(3)*cse(1) &
                  + sa(4)*Tc(j)*cse(1) + sa(5)*Tc(j)*Tc(j) &
                  + sa(6)*cse(1)*cse(1) + sa(7)*Tc(j)*Tc(j)*cse(1) &
                  + sa(8)*Tc(j)*cse(1)*cse(1) + sa(9)*Tc(j)*Tc(j)*Tc(j) &
                  + sa(10)*cse(1)*cse(1)*cse(1)
               a_ = 10.0**loga_
               b_ = sb(1)+sb(2)*Tc(j)+sb(3)*cse(1) + sb(4)*Tc(j)*cse(1) &
                  + sb(5)*Tc(j)*Tc(j) + sb(6)*cse(1)*cse(1) &
                  + sb(7)*Tc(j)*Tc(j)*cse(1) + sb(8)*Tc(j)*cse(1)*cse(1) &
                  + sb(9)*Tc(j)*Tc(j)*Tc(j)+sb(10)*cse(1)*cse(1)*cse(1)
               M3 = a_ * second**b_

               oM3 = 1./M3
               Mrat = M2*(M2*oM3)*(M2*oM3)*(M2*oM3)
               M0   = (M2*oM3)**mu_s
               slam1 = M2 * oM3 * Lam0
               slam2 = M2 * oM3 * Lam1

               do n = 1, nbs
                  N_s(n) = Mrat*(Kap0*exp(real(-slam1*Ds(n), kind=dp)) &
                      + Kap1*M0*Ds(n)**mu_s * exp(real(-slam2*Ds(n), kind=dp)))*dts(n)
               enddo

               t1 = 0.0_dp
               t2 = 0.0_dp
               t3 = 0.0_dp
               t4 = 0.0_dp
               z1 = 0.0_dp
               z2 = 0.0_dp
               z3 = 0.0_dp
               z4 = 0.0_dp
               y1 = 0.0_dp
               y2 = 0.0_dp
               y3 = 0.0_dp
               y4 = 0.0_dp
               do n2 = 1, nbr
                  massr = am_r * Dr(n2)**bm_r
                  do n = 1, nbs
                     masss = am_s * Ds(n)**bm_s

                     dvs = 0.5d0*((vr(n2) - vs(n)) + DABS(vr(n2)-vs(n)))
                     dvr = 0.5d0*((vs(n) - vr(n2)) + DABS(vs(n)-vr(n2)))

                     if (massr .gt. 1.5*masss) then
                     t1 = t1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvs*masss * N_s(n)* N_r(n2)
                     z1 = z1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvs*massr * N_s(n)* N_r(n2)
                     y1 = y1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvs       * N_s(n)* N_r(n2)
                     else
                     t3 = t3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvs*masss * N_s(n)* N_r(n2)
                     z3 = z3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvs*massr * N_s(n)* N_r(n2)
                     y3 = y3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvs       * N_s(n)* N_r(n2)
                     endif

                     if (massr .gt. 1.5*masss) then
                     t2 = t2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvr*massr * N_s(n)* N_r(n2)
                     y2 = y2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvr       * N_s(n)* N_r(n2)
                     z2 = z2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvr*masss * N_s(n)* N_r(n2)
                     else
                     t4 = t4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvr*massr * N_s(n)* N_r(n2)
                     y4 = y4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvr       * N_s(n)* N_r(n2)
                     z4 = z4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvr*masss * N_s(n)* N_r(n2)
                     endif

                  enddo
               enddo
               tcs_racs1(i,j,k,m) = t1
               tmr_racs1(i,j,k,m) = min(z1, r_r(m)*1.0_dp)
               tcs_racs2(i,j,k,m) = t3
               tmr_racs2(i,j,k,m) = z3
               tcr_sacr1(i,j,k,m) = t2
               tms_sacr1(i,j,k,m) = z2
               tcr_sacr2(i,j,k,m) = t4
               tms_sacr2(i,j,k,m) = z4
               tnr_racs1(i,j,k,m) = y1
               tnr_racs2(i,j,k,m) = y3
               tnr_sacr1(i,j,k,m) = y2
               tnr_sacr2(i,j,k,m) = y4
            enddo
         enddo
        enddo

        IF ( write_thompson_tables ) THEN
          write(0,*) "Writing "//qr_acr_qs_file//" in Thompson MP init"
          OPEN(63,file=qr_acr_qs_file,form="unformatted",err=9234)
          WRITE(63,err=9234)tcs_racs1
          WRITE(63,err=9234)tmr_racs1
          WRITE(63,err=9234)tcs_racs2
          WRITE(63,err=9234)tmr_racs2
          WRITE(63,err=9234)tcr_sacr1
          WRITE(63,err=9234)tms_sacr1
          WRITE(63,err=9234)tcr_sacr2
          WRITE(63,err=9234)tms_sacr2
          WRITE(63,err=9234)tnr_racs1
          WRITE(63,err=9234)tnr_racs2
          WRITE(63,err=9234)tnr_sacr1
          WRITE(63,err=9234)tnr_sacr2
          CLOSE(63)
          RETURN    ! ----- RETURN
 9234     CONTINUE
          write(0,*) "Error writing "//qr_acr_qs_file
        ENDIF
      ENDIF

   end subroutine qr_acr_qs
!+---+-----------------------------------------------------------------+
!ctrlL
!+---+-----------------------------------------------------------------+
!>\ingroup aathompson
!! This is a literal adaptation of Bigg (1954) probability of drops of
!! a particular volume freezing.  Given this probability, simply freeze
!! the proportion of drops summing their masses.
   subroutine freezeH2O(threads)

      implicit none

!..Interface variables
      integer, intent(in):: threads

!..Local variables
      integer:: i, j, k, m, n, n2
      INTEGER:: km, km_s, km_e
      real(dp) :: N_r, N_c
      real(dp), dimension(nbr):: massr
      real(dp), dimension(nbc):: massc
      real(dp) :: sum1, sum2, sumn1, sumn2, &
                         prob, vol, Texp, orho_w, &
                         lam_exp, lamr, N0_r, lamc, N0_c, y
      integer :: nu_c
      real(wp) :: T_adjust
      logical force_read_thompson, write_thompson_tables
      logical lexist,lopen
      integer good,ierr

!+---+
      force_read_thompson = .false.
      write_thompson_tables = .false.

      good = 0
        INQUIRE(FILE=freeze_h2o_file,EXIST=lexist)
#ifdef MPI
        call MPI_BARRIER(mpi_communicator,ierr)
#endif
        IF ( lexist ) THEN
          !write(0,*) "ThompMP: read "//freeze_h2o_file//" instead of computing"
          OPEN(63,file=freeze_h2o_file,form="unformatted",err=1234)
!sms$serial begin
          READ(63,err=1234)tpi_qrfz
          READ(63,err=1234)tni_qrfz
          READ(63,err=1234)tpg_qrfz
          READ(63,err=1234)tnr_qrfz
          READ(63,err=1234)tpi_qcfz
          READ(63,err=1234)tni_qcfz
!sms$serial end
          good = 1
 1234     CONTINUE
          IF ( good .NE. 1 ) THEN
            INQUIRE(63,opened=lopen)
            IF (lopen) THEN
              IF( force_read_thompson ) THEN
                write(0,*) "Error reading "//freeze_h2o_file//" Aborting because force_read_thompson is .true."
                return
              ENDIF
              CLOSE(63)
            ELSE
              IF( force_read_thompson ) THEN
                write(0,*) "Error opening "//freeze_h2o_file//" Aborting because force_read_thompson is .true."
                return
              ENDIF
            ENDIF
          ELSE
            INQUIRE(63,opened=lopen)
            IF (lopen) THEN
              CLOSE(63)
            ENDIF
          ENDIF
        ELSE
          IF( force_read_thompson ) THEN
            write(0,*) "Non-existent "//freeze_h2o_file//" Aborting because force_read_thompson is .true."
            return
          ENDIF
        ENDIF

      IF (.NOT. good .EQ. 1 ) THEN
        if (thompson_table_writer) then
          write_thompson_tables = .true.
          write(0,*) "ThompMP: computing freezeH2O"
        endif

        orho_w = 1./rho_w2

        do n2 = 1, nbr
         massr(n2) = am_r*Dr(n2)**bm_r
        enddo
        do n = 1, nbc
         massc(n) = am_r*Dc(n)**bm_r
        enddo

!..Freeze water (smallest drops become cloud ice, otherwise graupel).
        do km = km_s, km_e
         m = km / 45 + 1
         k = mod( km , 45 ) + 1
         T_adjust = MAX(-3.0, MIN(3.0 - ALOG10(Nt_IN(m)), 3.0))
!         print*, ' Freezing water for temp = ', -k
         Texp = exp( real(k, kind=dp) - T_adjust*1.0_dp ) - 1.0_dp
!$OMP PARALLEL DO SCHEDULE(dynamic) num_threads(threads) &
!$OMP PRIVATE(j,i,lam_exp,lamr,N0_r,sum1,sum2,sumn1,sumn2,n2,N_r,vol,prob)
         do j = 1, ntb_r1
            do i = 1, ntb_r
               lam_exp = (N0r_exp(j)*am_r*crg(1)/r_r(i))**ore1
               lamr = lam_exp * (crg(3)*org2*org1)**obmr
               N0_r = N0r_exp(j)/(crg(2)*lam_exp) * lamr**cre(2)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sumn1 = 0.0_dp
               sumn2 = 0.0_dp
               do n2 = nbr, 1, -1
                  N_r = N0_r*Dr(n2)**mu_r*exp(real(-lamr*Dr(n2), kind=dp))*dtr(n2)
                  vol = massr(n2)*orho_w
                  prob = max(0.0_dp, 1.0_dp - exp(-120.0_dp*vol*5.2e-4_dp * Texp))
                  if (massr(n2) .lt. xm0g) then
                     sumn1 = sumn1 + prob*N_r
                     sum1 = sum1 + prob*N_r*massr(n2)
                  else
                     sumn2 = sumn2 + prob*N_r
                     sum2 = sum2 + prob*N_r*massr(n2)
                  endif
                  if ((sum1+sum2).ge.r_r(i)) EXIT
               enddo
               tpi_qrfz(i,j,k,m) = sum1
               tni_qrfz(i,j,k,m) = sumn1
               tpg_qrfz(i,j,k,m) = sum2
               tnr_qrfz(i,j,k,m) = sumn2
            enddo
         enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(dynamic) num_threads(threads) &
!$OMP PRIVATE(j,i,nu_c,lamc,N0_c,sum1,sumn2,vol,prob,N_c)
         do j = 1, nbc
            nu_c = min(15, nint(1000.E6/t_Nc(j)) + 2)
            do i = 1, ntb_c
               lamc = (t_Nc(j)*am_r* ccg(2,nu_c) * ocg1(nu_c) / r_c(i))**obmr
               N0_c = t_Nc(j)*ocg1(nu_c) * lamc**cce(1,nu_c)
               sum1 = 0.0_dp
               sumn2 = 0.0_dp
               do n = nbc, 1, -1
                  vol = massc(n)*orho_w
                  prob = max(0.0_dp, 1.0_dp - exp(-120.0_dp*vol*5.2e-4_dp * Texp))
                  N_c = N0_c*Dc(n)**nu_c*EXP(-lamc*Dc(n))*dtc(n)
                  sumn2 = min(t_Nc(j), sumn2 + prob*N_c)
                  sum1 = sum1 + prob*N_c*massc(n)
                  if (sum1 .ge. r_c(i)) EXIT
               enddo
               tpi_qcfz(i,j,k,m) = sum1
               tni_qcfz(i,j,k,m) = sumn2
            enddo
         enddo
!$OMP END PARALLEL DO
        enddo
        
        IF ( write_thompson_tables ) THEN
          write(0,*) "Writing "//freeze_h2o_file//" in Thompson MP init"
          OPEN(63,file=freeze_h2o_file,form="unformatted",err=9234)
          WRITE(63,err=9234)tpi_qrfz
          WRITE(63,err=9234)tni_qrfz
          WRITE(63,err=9234)tpg_qrfz
          WRITE(63,err=9234)tnr_qrfz
          WRITE(63,err=9234)tpi_qcfz
          WRITE(63,err=9234)tni_qcfz
          CLOSE(63)
          RETURN    ! ----- RETURN
 9234     CONTINUE
          write(0,*) "Error writing "//freeze_h2o_file
          return
        ENDIF
      ENDIF

   end subroutine freezeH2O

!+---+-----------------------------------------------------------------+
!ctrlL
!+---+-----------------------------------------------------------------+
!>\ingroup aathompson
!! Cloud ice converting to snow since portion greater than min snow
!! size.  Given cloud ice content (kg/m**3), number concentration
!! (#/m**3) and gamma shape parameter, mu_i, break the distrib into
!! bins and figure out the mass/number of ice with sizes larger than
!! D0s.  Also, compute incomplete gamma function for the integration
!! of ice depositional growth from diameter=0 to D0s.  Amount of
!! ice depositional growth is this portion of distrib while larger
!! diameters contribute to snow growth (as in Harrington et al. 1995).
   subroutine qi_aut_qs

      implicit none

!..Local variables
      integer:: i, j, n2
      real(dp), dimension(nbi):: N_i
      real(dp) :: N0_i, lami, Di_mean, t1, t2
      real(wp) :: xlimit_intg

!+---+

      do j = 1, ntb_i1
         do i = 1, ntb_i
            lami = (am_i*cig(2)*oig1*Nt_i(j)/r_i(i))**obmi
            Di_mean = (bm_i + mu_i + 1.) / lami
            N0_i = Nt_i(j)*oig1 * lami**cie(1)
            t1 = 0.0_dp
            t2 = 0.0_dp
            if (SNGL(Di_mean) .gt. 5.*D0s) then
             t1 = r_i(i)
             t2 = Nt_i(j)
             tpi_ide(i,j) = 0.0_dp
            elseif (SNGL(Di_mean) .lt. D0i) then
             t1 = 0.0_dp
             t2 = 0.0_dp
             tpi_ide(i,j) = 1.0_dp
            else
             xlimit_intg = lami*D0s
             tpi_ide(i,j) = GAMMP(mu_i+2.0, xlimit_intg) * 1.0_dp
             do n2 = 1, nbi
               N_i(n2) = N0_i*Di(n2)**mu_i * exp(real(-lami*Di(n2), kind=dp))*dti(n2)
               if (Di(n2).ge.D0s) then
                  t1 = t1 + N_i(n2) * am_i*Di(n2)**bm_i
                  t2 = t2 + N_i(n2)
               endif
             enddo
            endif
            tps_iaus(i,j) = t1
            tni_iaus(i,j) = t2
         enddo
      enddo

   end subroutine qi_aut_qs
!ctrlL
!+---+-----------------------------------------------------------------+
!>\ingroup aathompson
!! Variable collision efficiency for rain collecting cloud water using
!! method of Beard and Grover, 1974 if a/A less than 0.25; otherwise
!! uses polynomials to get close match of Pruppacher & Klett Fig 14-9.
   subroutine table_Efrw

      implicit none

!..Local variables
      real(dp) :: vtr, stokes, reynolds, Ef_rw
      real(dp) :: p, yc0, F, G, H, z, K0, X
      integer:: i, j

      do j = 1, nbc
      do i = 1, nbr
         Ef_rw = 0.0
         p = Dc(j)/Dr(i)
         if (Dr(i).lt.50.E-6 .or. Dc(j).lt.3.E-6) then
          t_Efrw(i,j) = 0.0
         elseif (p.gt.0.25) then
          X = Dc(j)*1.e6_dp
          if (Dr(i) .lt. 75.e-6) then
             Ef_rw = 0.026794*X - 0.20604
          elseif (Dr(i) .lt. 125.e-6) then
             Ef_rw = -0.00066842*X*X + 0.061542*X - 0.37089
          elseif (Dr(i) .lt. 175.e-6) then
             Ef_rw = 4.091e-06*X*X*X*X - 0.00030908*X*X*X               &
                   + 0.0066237*X*X - 0.0013687*X - 0.073022
          elseif (Dr(i) .lt. 250.e-6) then
             Ef_rw = 9.6719e-5*X*X*X - 0.0068901*X*X + 0.17305*X        &
                   - 0.65988
          elseif (Dr(i) .lt. 350.e-6) then
             Ef_rw = 9.0488e-5*X*X*X - 0.006585*X*X + 0.16606*X         &
                   - 0.56125
          else
             Ef_rw = 0.00010721*X*X*X - 0.0072962*X*X + 0.1704*X        &
                   - 0.46929
          endif
         else
          vtr = -0.1021 + 4.932E3*Dr(i) - 0.9551E6*Dr(i)*Dr(i) &
              + 0.07934E9*Dr(i)*Dr(i)*Dr(i) &
              - 0.002362E12*Dr(i)*Dr(i)*Dr(i)*Dr(i)
          stokes = Dc(j)*Dc(j)*vtr*rho_w2/(9.*1.718E-5*Dr(i))
          reynolds = 9.*stokes/(p*p*rho_w2)

          F = log(real(reynolds, kind=dp))
          G = -0.1007_dp - 0.358_dp*F + 0.0261_dp*F*F
          K0 = exp(G)
          z = log(stokes/(K0+1.e-15_dp))
          H = 0.1465D0 + 1.302D0*z - 0.607D0*z*z + 0.293D0*z*z*z
          yc0 = 2.0_dp/PI * ATAN(H)
          Ef_rw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))

         endif

         t_Efrw(i,j) = max(0.0, min(SNGL(Ef_rw), 0.95))

      enddo
      enddo

   end subroutine table_Efrw
!ctrlL
!+---+-----------------------------------------------------------------+
!>\ingroup aathompson
!! Variable collision efficiency for snow collecting cloud water using
!! method of Wang and Ji, 2000 except equate melted snow diameter to
!! their "effective collision cross-section."
   subroutine table_Efsw

      implicit none

!..Local variables
      real(dp) :: Ds_m, vts, vtc, stokes, reynolds, Ef_sw
      real(dp) :: p, yc0, F, G, H, z, K0
      integer:: i, j

      do j = 1, nbc
      vtc = 1.19e4_dp * (1.0e4_dp*Dc(j)*Dc(j)*0.25_dp)
      do i = 1, nbs
         vts = av_s*Ds(i)**bv_s * exp(real(-fv_s*Ds(i), kind=dp)) - vtc
         Ds_m = (am_s*Ds(i)**bm_s / am_r)**obmr
         p = Dc(j)/Ds_m
         if (p.gt.0.25 .or. Ds(i).lt.D0s .or. Dc(j).lt.6.E-6 &
               .or. vts.lt.1.E-3) then
          t_Efsw(i,j) = 0.0
         else
          stokes = Dc(j)*Dc(j)*vts*rho_w2/(9.*1.718E-5*Ds_m)
          reynolds = 9.*stokes/(p*p*rho_w2)

          F = log(real(reynolds, kind=dp))
          G = -0.1007_dp - 0.358_dp*F + 0.0261_dp*F*F
          K0 = exp(G)
          z = log(stokes/(K0+1.e-15_dp))
          H = 0.1465D0 + 1.302D0*z - 0.607D0*z*z + 0.293D0*z*z*z
          yc0 = 2.0_dp/PI * ATAN(H)
          Ef_sw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))

          t_Efsw(i,j) = max(0.0, min(SNGL(Ef_sw), 0.95))
         endif

      enddo
      enddo

   end subroutine table_Efsw
!ctrlL
!+---+-----------------------------------------------------------------+
!ctrlL
!+---+-----------------------------------------------------------------+
!>\ingroup aathompson
!! Integrate rain size distribution from zero to D-star to compute the
!! number of drops smaller than D-star that evaporate in a single
!! timestep.  Drops larger than D-star dont evaporate entirely so do
!! not affect number concentration.
   subroutine table_dropEvap

      implicit none

!..Local variables
      integer:: i, j, k, n
      real(dp), dimension(nbc):: N_c, massc
      real(dp) :: summ, summ2, lamc, N0_c
      integer:: nu_c
!      real(dp) :: Nt_r, N0, lam_exp, lam
!      real(wp) :: xlimit_intg

      do n = 1, nbc
         massc(n) = am_r*Dc(n)**bm_r
      enddo

      do k = 1, nbc
         nu_c = min(15, nint(1000.E6/t_Nc(k)) + 2)
         do j = 1, ntb_c
            lamc = (t_Nc(k)*am_r* ccg(2,nu_c)*ocg1(nu_c) / r_c(j))**obmr
            N0_c = t_Nc(k)*ocg1(nu_c) * lamc**cce(1,nu_c)
            do i = 1, nbc
!-GT           tnc_wev(i,j,k) = GAMMP(nu_c+1., SNGL(Dc(i)*lamc))*t_Nc(k)
               N_c(i) = N0_c* Dc(i)**nu_c*EXP(-lamc*Dc(i))*dtc(i)
!     if(j.eq.18 .and. k.eq.50) print*, ' N_c = ', N_c(i)
               summ = 0.
               summ2 = 0.
               do n = 1, i
                  summ = summ + massc(n)*N_c(n)
                  summ2 = summ2 + N_c(n)
               enddo
!      if(j.eq.18 .and. k.eq.50) print*, '  DEBUG-TABLE: ', r_c(j), t_Nc(k), summ2, summ
               tpc_wev(i,j,k) = summ
               tnc_wev(i,j,k) = summ2
            enddo
         enddo
      enddo

!
!..To do the same thing for rain.
!
!     do k = 1, ntb_r
!     do j = 1, ntb_r1
!        lam_exp = (N0r_exp(j)*am_r*crg(1)/r_r(k))**ore1
!        lam = lam_exp * (crg(3)*org2*org1)**obmr
!        N0 = N0r_exp(j)/(crg(2)*lam_exp) * lam**cre(2)
!        Nt_r = N0 * crg(2) / lam**cre(2)
!        do i = 1, nbr
!           xlimit_intg = lam*Dr(i)
!           tnr_rev(i,j,k) = GAMMP(mu_r+1.0, xlimit_intg) * Nt_r
!        enddo
!     enddo
!     enddo

! TO APPLY TABLE ABOVE
!..Rain lookup table indexes.
!         Dr_star = sqrt(-2.0_dp*DT * t1_evap/(2.*PI) &
!                 * 0.78*4.*diffu(k)*xsat*rvs/rho_w)
!         idx_d = nint(1.0 + real(nbr, kind=wp) * log(real(Dr_star/D0r, kind=dp))             &
!               / log(real(Dr(nbr)/D0r, kind=dp)))
!         idx_d = max(1, min(idx_d, nbr))
!
!         nir = nint(log10(real(rr(k), kind=wp)))
!         do nn = nir-1, nir+1
!            n = nn
!            if ( (rr(k)/10.**nn).ge.1.0 .and. &
!                 (rr(k)/10.**nn).lt.10.0) goto 154
!         enddo
!154      continue
!         idx_r = int(rr(k)/10.**n) + 10*(n-nir2) - (n-nir2)
!         idx_r = max(1, min(idx_r, ntb_r))
!
!         lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
!         lam_exp = lamr * (crg(3)*org2*org1)**bm_r
!         N0_exp = org1*rr(k)/am_r * lam_exp**cre(1)
!         nir = nint(log10(real(N0_exp, kind=dp))
!         do nn = nir-1, nir+1
!            n = nn
!            if ( (N0_exp/10.**nn).ge.1.0 .and. &
!                 (N0_exp/10.**nn).lt.10.0) goto 155
!         enddo
!155      continue
!         idx_r1 = int(N0_exp/10.**n) + 10*(n-nir3) - (n-nir3)
!         idx_r1 = max(1, min(idx_r1, ntb_r1))
!
!         pnr_rev(k) = min(nr(k)*odts, SNGL(tnr_rev(idx_d,idx_r1,idx_r) &   ! RAIN2M
!                    * odts))

   end subroutine table_dropEvap
!
!ctrlL
!+---+-----------------------------------------------------------------+
!>\ingroup aathompson
!! Fill the table of CCN activation data created from parcel model run
!! by Trude Eidhammer with inputs of aerosol number concentration,
!! vertical velocity, temperature, lognormal mean aerosol radius, and
!! hygroscopicity, kappa.  The data are read from external file and
!! contain activated fraction of CCN for given conditions.
   subroutine table_ccnAct(errmess,errflag)

      implicit none

!..Error handling variables
      character(len=*), intent(inout) :: errmess
      integer,          intent(inout) :: errflag

!..Local variables
      integer:: iunit_mp_th1, i
      logical:: opened

      iunit_mp_th1 = -1
      do_loop_ccn : do i = 20, 99
         INQUIRE (i, OPENED=opened)
         if (.not. opened) then
            iunit_mp_th1 = i
            exit do_loop_ccn
         endif
      enddo do_loop_ccn

      if (iunit_mp_th1 < 0) then
         write(0,*) 'module_mp_thompson: table_ccnAct: '//   &
                   'Can not find unused fortran unit to read in lookup table.'
         return
      endif

      !WRITE(*, '(A,I2)') 'module_mp_thompson: opening CCN_ACTIVATE.BIN on unit ',iunit_mp_th1
      OPEN(iunit_mp_th1, FILE='CCN_ACTIVATE.BIN',                      &
            FORM='UNFORMATTED', STATUS='OLD', CONVERT='BIG_ENDIAN', ERR=9009)

!sms$serial begin
      READ(iunit_mp_th1, ERR=9010) tnccn_act
!sms$serial end

      return
 9009 CONTINUE
      WRITE(errmess , '(A,I2)') 'module_mp_thompson: error opening CCN_ACTIVATE.BIN on unit ',iunit_mp_th1
      errflag = 1
      return
 9010 CONTINUE
      WRITE(errmess , '(A,I2)') 'module_mp_thompson: error reading CCN_ACTIVATE.BIN on unit ',iunit_mp_th1
      errflag = 1
      return

   end subroutine table_ccnAct

! !
! !-------------------------------------------------------------------
!    SUBROUTINE semi_lagrange_sedim(km,dzl,wwl,rql,precip,pfsan,dt,R1)
! !-------------------------------------------------------------------
! !
! ! This routine is a semi-Lagrangain forward advection for hydrometeors
! ! with mass conservation and positive definite advection
! ! 2nd order interpolation with monotonic piecewise parabolic method is used.
! ! This routine is under assumption of decfl < 1 for semi_Lagrangian
! !
! ! dzl    depth of model layer in meter
! ! wwl    terminal velocity at model layer m/s
! ! rql    dry air density*mixing ratio
! ! precip precipitation at surface 
! ! dt     time step
! !
! ! author: hann-ming henry juang <henry.juang@noaa.gov>
! !         implemented by song-you hong
! ! reference: Juang, H.-M., and S.-Y. Hong, 2010: Forward semi-Lagrangian advection
! !         with mass conservation and positive definiteness for falling
! !         hydrometeors. *Mon.  Wea. Rev.*, *138*, 1778-1791
! !
!       implicit none

!       integer, intent(in) :: km
!       real(wp), intent(in) ::  dt, R1
!       real(wp), intent(in) :: dzl(km),wwl(km)
!       real(wp), intent(out) :: precip
!       real(wp), intent(inout) :: rql(km)
!       real(wp), intent(out)  :: pfsan(km)
!       integer ::  k,m,kk,kb,kt
!       real(wp) :: tl,tl2,qql,dql,qqd
!       real(wp) :: th,th2,qqh,dqh
!       real(wp) :: zsum,qsum,dim,dip,con1,fa1,fa2
!       real(wp) :: allold, decfl
!       real(wp) :: dz(km), ww(km), qq(km)
!       real(wp) :: wi(km+1), zi(km+1), za(km+2)
!       real(wp) :: qn(km)
!       real(wp) :: dza(km+1), qa(km+1), qmi(km+1), qpi(km+1)
!       real(wp) :: net_flx(km)
! !
!       precip = 0.0
!       qa(:) = 0.0
!       qq(:) = 0.0
!       dz(:) = dzl(:)
!       ww(:) = wwl(:)
!       do k = 1,km
!         if(rql(k).gt.R1) then 
!           qq(k) = rql(k) 
!         else 
!           ww(k) = 0.0 
!         endif
!         pfsan(k) = 0.0
!         net_flx(k) = 0.0
!       enddo
! ! skip for no precipitation for all layers
!       allold = 0.0
!       do k=1,km
!         allold = allold + qq(k)
!       enddo
!       if(allold.le.0.0) then
!          return 
!       endif
! !
! ! compute interface values
!       zi(1)=0.0
!       do k=1,km
!         zi(k+1) = zi(k)+dz(k)
!       enddo
! !     n=1
! ! plm is 2nd order, we can use 2nd order wi or 3rd order wi
! ! 2nd order interpolation to get wi
!       wi(1) = ww(1)
!       wi(km+1) = ww(km)
!       do k=2,km
!         wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
!       enddo
! ! 3rd order interpolation to get wi
!       fa1 = 9./16.
!       fa2 = 1./16.
!       wi(1) = ww(1)
!       wi(2) = 0.5*(ww(2)+ww(1))
!       do k=3,km-1
!         wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
!       enddo
!       wi(km) = 0.5*(ww(km)+ww(km-1))
!       wi(km+1) = ww(km)

! ! terminate of top of raingroup
!       do k=2,km
!         if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
!       enddo

! ! diffusivity of wi
!       con1 = 0.05
!       do k=km,1,-1
!         decfl = (wi(k+1)-wi(k))*dt/dz(k)
!         if( decfl .gt. con1 ) then
!           wi(k) = wi(k+1) - con1*dz(k)/dt
!         endif
!       enddo
! ! compute arrival point
!       do k=1,km+1
!         za(k) = zi(k) - wi(k)*dt
!       enddo
!       za(km+2) = zi(km+1)

!       do k=1,km+1
!         dza(k) = za(k+1)-za(k)
!       enddo

! ! computer deformation at arrival point
!       do k=1,km
!         qa(k) = qq(k)*dz(k)/dza(k)
!       enddo
!       qa(km+1) = 0.0

! ! estimate values at arrival cell interface with monotone
!       do k=2,km
!         dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
!         dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
!         if( dip*dim.le.0.0 ) then
!           qmi(k)=qa(k)
!           qpi(k)=qa(k)
!         else
!           qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
!           qmi(k)=2.0*qa(k)-qpi(k)
!           if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
!             qpi(k) = qa(k)
!             qmi(k) = qa(k)
!           endif
!         endif
!       enddo
!       qpi(1)=qa(1)
!       qmi(1)=qa(1)
!       qmi(km+1)=qa(km+1)
!       qpi(km+1)=qa(km+1)

! ! interpolation to regular point
!       qn = 0.0
!       kb=1
!       kt=1
!       intp : do k=1,km
!              kb=max(kb-1,1)
!              kt=max(kt-1,1)
! ! find kb and kt
!              if( zi(k).ge.za(km+1) ) then
!                exit intp
!              else
!                find_kb : do kk=kb,km
!                          if( zi(k).le.za(kk+1) ) then
!                            kb = kk
!                            exit find_kb
!                          else
!                            cycle find_kb
!                          endif
!                enddo find_kb
!                find_kt : do kk=kt,km+2
!                          if( zi(k+1).le.za(kk) ) then
!                            kt = kk
!                            exit find_kt
!                          else
!                            cycle find_kt
!                          endif
!                enddo find_kt
!                kt = kt - 1
! ! compute q with piecewise constant method
!                if( kt.eq.kb ) then
!                  tl=(zi(k)-za(kb))/dza(kb)
!                  th=(zi(k+1)-za(kb))/dza(kb)
!                  tl2=tl*tl
!                  th2=th*th
!                  qqd=0.5*(qpi(kb)-qmi(kb))
!                  qqh=qqd*th2+qmi(kb)*th
!                  qql=qqd*tl2+qmi(kb)*tl
!                  qn(k) = (qqh-qql)/(th-tl)
!                else if( kt.gt.kb ) then
!                  tl=(zi(k)-za(kb))/dza(kb)
!                  tl2=tl*tl
!                  qqd=0.5*(qpi(kb)-qmi(kb))
!                  qql=qqd*tl2+qmi(kb)*tl
!                  dql = qa(kb)-qql
!                  zsum  = (1.-tl)*dza(kb)
!                  qsum  = dql*dza(kb)
!                  if( kt-kb.gt.1 ) then
!                  do m=kb+1,kt-1
!                    zsum = zsum + dza(m)
!                    qsum = qsum + qa(m) * dza(m)
!                  enddo
!                  endif
!                  th=(zi(k+1)-za(kt))/dza(kt)
!                  th2=th*th
!                  qqd=0.5*(qpi(kt)-qmi(kt))
!                  dqh=qqd*th2+qmi(kt)*th
!                  zsum  = zsum + th*dza(kt)
!                  qsum  = qsum + dqh*dza(kt)
!                  qn(k) = qsum/zsum
!                endif
!                cycle intp
!              endif

!        enddo intp

! ! rain out
!       sum_precip: do k=1,km
!                     if( za(k).lt.0.0 .and. za(k+1).le.0.0 ) then
!                       precip = precip + qa(k)*dza(k)
!                       net_flx(k) =  qa(k)*dza(k)
!                       cycle sum_precip
!                     else if ( za(k).lt.0.0 .and. za(k+1).gt.0.0 ) then
!                       th = (0.0-za(k))/dza(k)
!                       th2 = th*th
!                       qqd = 0.5*(qpi(k)-qmi(k))
!                       qqh = qqd*th2+qmi(k)*th
!                       precip = precip + qqh*dza(k)
!                       net_flx(k) = qqh*dza(k)
!                       exit sum_precip
!                     endif
!                     exit sum_precip
!       enddo sum_precip

! ! calculating precipitation fluxes
!       do k=km,1,-1
!          if(k == km) then
!            pfsan(k) = net_flx(k)
!          else
!            pfsan(k) = pfsan(k+1) + net_flx(k)
!          end if
!       enddo
! !
! ! replace the new values
!       rql(:) = max(qn(:),R1)

!    END SUBROUTINE semi_lagrange_sedim

! !>\ingroup aathompson
! !! @brief Calculates graupel size distribution parameters
! !!
! !! Calculates graupel intercept and slope parameters for
! !! for a vertical column 
! !!  
! !! @param[in]    kts     integer start index for vertical column
! !! @param[in]    kte     integer end index for vertical column
! !! @param[in]    rand1   real random number for stochastic physics
! !! @param[in]    rg      real array, size(kts:kte) for graupel mass concentration [kg m^3]
! !! @param[out]   ilamg   double array, size(kts:kte) for inverse graupel slope parameter [m]
! !! @param[out]   N0_g    double array, size(kts:kte) for graupel intercept paramter [m-4]
!    subroutine graupel_psd_parameters(kts, kte, rand1, rg, ilamg, N0_g)

!       implicit none

!       integer, intent(in) :: kts, kte
!       real(wp), intent(in) :: rand1
!       real(wp), intent(in) :: rg(:)
!       real(dp), intent(out) :: ilamg(:), N0_g(:)

!       integer :: k
!       real(wp) :: ygra1, zans1
!       real(dp) :: N0_exp, lam_exp, lamg

!       do k = kte, kts, -1
!          ygra1 = alog10(max(1.e-9, rg(k)))
!          zans1 = 3.4 + 2./7.*(ygra1+8.) + rand1
!          N0_exp = 10.**(zans1)
!          N0_exp = max(real(gonv_min, kind=dp), min(N0_exp, real(gonv_max, kind=dp)))
!          lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
!          lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
!          ilamg(k) = 1./lamg
!          N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
!       enddo

!    end subroutine graupel_psd_parameters

! !>\ingroup aathompson
! !! @brief Calculates graupel/hail maximum diameter
! !!
! !! Calculates graupel/hail maximum diameter (currently the 99th percentile of mass distribtuion)
! !! for a vertical column 
! !!  
! !! @param[in]    kts             integer start index for vertical column
! !! @param[in]    kte             integer end index for vertical column
! !! @param[in]    qg              real array, size(kts:kte) for graupel mass mixing ratio [kg kg^-1]
! !! @param[in]    temperature     double array, size(kts:kte) temperature [K]
! !! @param[in]    pressure        double array, size(kts:kte) pressure [Pa]
! !! @param[in]    qv              real array, size(kts:kte) water vapor mixing ratio [kg kg^-1]
! !! @param[out]   max_hail_diam   real maximum hail diameter [m]
!    function hail_mass_99th_percentile(kts, kte, qg, temperature, pressure, qv) result(max_hail_diam)

!       implicit none
      
!       integer, intent(in) :: kts, kte
!       real(wp), intent(in) :: qg(:), temperature(:), pressure(:), qv(:)
!       real(wp) :: max_hail_diam

!       integer :: k
!       real(wp) :: rho(kts:kte), rg(kts:kte), max_hail_column(kts:kte)
!       real(dp) :: ilamg(kts:kte), N0_g(kts:kte)
!       real(wp), parameter :: random_number = 0.

!       max_hail_column = 0.
!       rg = 0.
!       do k = kts, kte
!          rho(k) = RoverRv*pressure(k) / (R*temperature(k)*(max(1.e-10, qv(k))+RoverRv))
!          if (qg(k) .gt. R1) then
!             rg(k) = qg(k)*rho(k)
!          else
!             rg(k) = R1
!          endif 
!       enddo 

!       call graupel_psd_parameters(kts, kte, random_number, rg, ilamg, N0_g)

!       where(rg .gt. 1.e-9) max_hail_column = 10.05 * ilamg
!       max_hail_diam = max_hail_column(kts)
      
!    end function hail_mass_99th_percentile

!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
end module module_mp_thompson
!+---+-----------------------------------------------------------------+

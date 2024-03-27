!=================================================================================================================
!
module module_mp_thompson

    use mpas_kind_types, only: wp => RKIND, sp => R4KIND, dp => R8KIND

    use module_mp_thompson_params
    use module_mp_thompson_main
    use module_mp_thompson_utils
    use mpas_atmphys_functions, only: gammp, wgamma, rslf, rsif
    use mpas_atmphys_utilities
    use mpas_io_units, only : mpas_new_unit, mpas_release_unit
    use mp_radar

contains

!=================================================================================================================
    subroutine thompson_init(l_mp_tables)
        implicit none
!=================================================================================================================

!input arguments:
        logical,intent(in) :: l_mp_tables

        integer, parameter :: open_OK = 0
        integer :: i, j, k, l, m, n
        integer :: istat
        logical :: micro_init
        integer :: mp_unit

!..Allocate space for lookup tables (J. Michalakes 2009Jun08).
        ! AAJ hail aware is off for now
        ! if (PRESENT(ng)) then
        !     is_hail_aware = .TRUE.
        !     dimNRHG = NRHG
        ! else
        ! av_g(idx_bg1) = av_g_old
        ! bv_g(idx_bg1) = bv_g_old
        ! dimNRHG = NRHG1
        ! endif

        dimNRHG = NRHG
        micro_init = .false.

!=================================================================================================================

!..Allocate space for lookup tables (J. Michalakes 2009Jun08).

        if (.not. allocated(tcg_racg) ) then
            allocate(tcg_racg(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
            micro_init = .true.
        endif

        if (.not. allocated(tmr_racg)) allocate(tmr_racg(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
        if (.not. allocated(tcr_gacr)) allocate(tcr_gacr(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
        ! if (.not. allocated(tmg_gacr)) allocate(tmg_gacr(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
        if (.not. allocated(tnr_racg)) allocate(tnr_racg(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))
        if (.not. allocated(tnr_gacr)) allocate(tnr_gacr(ntb_g1,ntb_g,dimNRHG,ntb_r1,ntb_r))

        if (.not. allocated(tcs_racs1)) allocate(tcs_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tmr_racs1)) allocate(tmr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tcs_racs2)) allocate(tcs_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tmr_racs2)) allocate(tmr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tcr_sacr1)) allocate(tcr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tms_sacr1)) allocate(tms_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tcr_sacr2)) allocate(tcr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tms_sacr2)) allocate(tms_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tnr_racs1)) allocate(tnr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tnr_racs2)) allocate(tnr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tnr_sacr1)) allocate(tnr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
        if (.not. allocated(tnr_sacr2)) allocate(tnr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))

        if (.not. allocated(tpi_qcfz)) allocate(tpi_qcfz(ntb_c,nbc,45,ntb_in))
        if (.not. allocated(tni_qcfz)) allocate(tni_qcfz(ntb_c,nbc,45,ntb_in))

        if (.not. allocated(tpi_qrfz)) allocate(tpi_qrfz(ntb_r,ntb_r1,45,ntb_in))
        if (.not. allocated(tpg_qrfz)) allocate(tpg_qrfz(ntb_r,ntb_r1,45,ntb_in))
        if (.not. allocated(tni_qrfz)) allocate(tni_qrfz(ntb_r,ntb_r1,45,ntb_in))
        if (.not. allocated(tnr_qrfz)) allocate(tnr_qrfz(ntb_r,ntb_r1,45,ntb_in))

        if (.not. allocated(tps_iaus)) allocate(tps_iaus(ntb_i,ntb_i1))
        if (.not. allocated(tni_iaus)) allocate(tni_iaus(ntb_i,ntb_i1))
        if (.not. allocated(tpi_ide)) allocate(tpi_ide(ntb_i,ntb_i1))

        if (.not. allocated(t_efrw)) allocate(t_efrw(nbr,nbc))
        if (.not. allocated(t_efsw)) allocate(t_efsw(nbs,nbc))

        if (.not. allocated(tnr_rev)) allocate(tnr_rev(nbr, ntb_r1, ntb_r))
        if (.not. allocated(tpc_wev)) allocate(tpc_wev(nbc,ntb_c,nbc))
        if (.not. allocated(tnc_wev)) allocate(tnc_wev(nbc,ntb_c,nbc))

        if (.not. allocated(tnccn_act)) allocate(tnccn_act(ntb_arc,ntb_arw,ntb_art,ntb_arr,ntb_ark))

        if (micro_init) then

!..From Martin et al. (1994), assign gamma shape parameter mu for cloud
!.. drops according to general dispersion characteristics (disp=~0.25
!.. for Maritime and 0.45 for Continental).
!.. disp=SQRT((mu+2)/(mu+1) - 1) so mu varies from 15 for Maritime
!.. to 2 for really dirty air.  This not used in 2-moment cloud water
!.. scheme and nu_c used instead and varies from 2 to 15 (integer-only).
!     mu_c = MIN(15., (1000.E6/Nt_c + 2.))

!..Schmidt number to one-third used numerous times.
            Sc3 = Sc**(1./3.)

!..Compute min ice diam from mass, min snow/graupel mass from diam.
            D0i = (xm0i/am_i)**(1./bm_i)
            xm0s = am_s * D0s**bm_s
            xm0g = am_g(NRHG) * D0g**bm_g

!..These constants various exponents and gamma() assoc with cloud,
!.. rain, snow, and graupel.
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
!..Simplify various rate eqns the best we can now.
!+---+-----------------------------------------------------------------+

!..Rain collecting cloud water and cloud ice
            t1_qr_qc = PI*.25*av_r * crg(9)
            t1_qr_qi = PI*.25*av_r * crg(9)
            t2_qr_qi = PI*.25*am_r*av_r * crg(8)

!..Graupel collecting cloud water
!     t1_qg_qc = PI*.25*av_g * cgg(9)

!..Snow collecting cloud water
            t1_qs_qc = PI*.25*av_s

!..Snow collecting cloud ice
            t1_qs_qi = PI*.25*av_s

!..Evaporation of rain; ignore depositional growth of rain.
            t1_qr_ev = 0.78 * crg(10)
            t2_qr_ev = 0.308*Sc3*SQRT(av_r) * crg(11)

!..Sublimation/depositional growth of snow
            t1_qs_sd = 0.86
            t2_qs_sd = 0.28*Sc3*SQRT(av_s)

!..Melting of snow
            t1_qs_me = PI*4.*C_sqrd*olfus * 0.86
            t2_qs_me = PI*4.*C_sqrd*olfus * 0.28*Sc3*SQRT(av_s)

!..Sublimation/depositional growth of graupel
            t1_qg_sd = 0.86 * cgg(10,1)
!     t2_qg_sd = 0.28*Sc3*SQRT(av_g) * cgg(11)

!..Melting of graupel
            t1_qg_me = PI*4.*C_cube*olfus * 0.86 * cgg(10,1)
!     t2_qg_me = PI*4.*C_cube*olfus * 0.28*Sc3*SQRT(av_g) * cgg(11)

!..Constants for helping find lookup table indexes.
            nic2 = NINT(ALOG10(r_c(1)))
            nii2 = NINT(ALOG10(r_i(1)))
            nii3 = NINT(ALOG10(Nt_i(1)))
            nir2 = NINT(ALOG10(r_r(1)))
            nir3 = NINT(ALOG10(N0r_exp(1)))
            nis2 = NINT(ALOG10(r_s(1)))
            nig2 = NINT(ALOG10(r_g(1)))
            nig3 = NINT(ALOG10(N0g_exp(1)))
            niIN2 = NINT(ALOG10(Nt_IN(1)))

!..Create bins of cloud water (from min diameter up to 100 microns).
            Dc(1) = D0c*1.0d0
            dtc(1) = D0c*1.0d0
            do n = 2, nbc
                Dc(n) = Dc(n-1) + 1.0D-6
                dtc(n) = (Dc(n) - Dc(n-1))
            enddo

!..Create bins of cloud ice (from min diameter up to 5x min snow size).
            xDx(1) = D0i*1.0d0
            xDx(nbi+1) = 2.0d0*D0s
            do n = 2, nbi
                xDx(n) = DEXP(REAL(n-1,KIND=dp)/REAL(nbi,KIND=dp) & 
                        *DLOG(xDx(nbi+1)/xDx(1)) +DLOG(xDx(1)))
            enddo
            do n = 1, nbi
                Di(n) = DSQRT(xDx(n)*xDx(n+1))
                dti(n) = xDx(n+1) - xDx(n)
            enddo

!..Create bins of rain (from min diameter up to 5 mm).
            xDx(1) = D0r*1.0d0
            xDx(nbr+1) = 0.005d0
            do n = 2, nbr
                xDx(n) = DEXP(REAL(n-1,KIND=dp)/REAL(nbr,KIND=dp) &
                    *DLOG(xDx(nbr+1)/xDx(1)) +DLOG(xDx(1)))
            enddo
            do n = 1, nbr
                Dr(n) = DSQRT(xDx(n)*xDx(n+1))
                dtr(n) = xDx(n+1) - xDx(n)
            enddo

!..Create bins of snow (from min diameter up to 2 cm).
            xDx(1) = D0s*1.0d0
            xDx(nbs+1) = 0.02d0
            do n = 2, nbs
                xDx(n) = DEXP(REAL(n-1,KIND=dp)/REAL(nbs,KIND=dp) &
                    *DLOG(xDx(nbs+1)/xDx(1)) +DLOG(xDx(1)))
            enddo
            do n = 1, nbs
                Ds(n) = DSQRT(xDx(n)*xDx(n+1))
                dts(n) = xDx(n+1) - xDx(n)
            enddo

!..Create bins of graupel (from min diameter up to 5 cm).
            xDx(1) = D0g*1.0d0
            xDx(nbg+1) = 0.05d0
            do n = 2, nbg
                xDx(n) = DEXP(REAL(n-1,KIND=dp)/REAL(nbg,KIND=dp) &
                    *DLOG(xDx(nbg+1)/xDx(1)) +DLOG(xDx(1)))
            enddo
            do n = 1, nbg
                Dg(n) = DSQRT(xDx(n)*xDx(n+1))
                dtg(n) = xDx(n+1) - xDx(n)
            enddo

!..Create bins of cloud droplet number concentration (1 to 3000 per cc).
            xDx(1) = 1.0d0
            xDx(nbc+1) = 3000.0d0
            do n = 2, nbc
                xDx(n) = DEXP(REAL(n-1,KIND=dp)/REAL(nbc,KIND=dp) &
                    *DLOG(xDx(nbc+1)/xDx(1)) +DLOG(xDx(1)))
            enddo
            do n = 1, nbc
                t_Nc(n) = DSQRT(xDx(n)*xDx(n+1)) * 1.D6
            enddo
            nic1 = DLOG(t_Nc(nbc)/t_Nc(1))

!+---+-----------------------------------------------------------------+
!..Create lookup tables for most costly calculations.
!+---+-----------------------------------------------------------------+

            do m = 1, ntb_r
                do k = 1, ntb_r1
                    do n = 1, dimNRHG
                        do j = 1, ntb_g
                            do i = 1, ntb_g1
                                tcg_racg(i,j,n,k,m) = 0.0d0
                                tmr_racg(i,j,n,k,m) = 0.0d0
                                tcr_gacr(i,j,n,k,m) = 0.0d0
                                !tmg_gacr(i,j,n,k,m) = 0.0d0
                                tnr_racg(i,j,n,k,m) = 0.0d0
                                tnr_gacr(i,j,n,k,m) = 0.0d0
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            do m = 1, ntb_r
                do k = 1, ntb_r1
                    do j = 1, ntb_t
                        do i = 1, ntb_s
                            tcs_racs1(i,j,k,m) = 0.0d0
                            tmr_racs1(i,j,k,m) = 0.0d0
                            tcs_racs2(i,j,k,m) = 0.0d0
                            tmr_racs2(i,j,k,m) = 0.0d0
                            tcr_sacr1(i,j,k,m) = 0.0d0
                            tms_sacr1(i,j,k,m) = 0.0d0
                            tcr_sacr2(i,j,k,m) = 0.0d0
                            tms_sacr2(i,j,k,m) = 0.0d0
                            tnr_racs1(i,j,k,m) = 0.0d0
                            tnr_racs2(i,j,k,m) = 0.0d0
                            tnr_sacr1(i,j,k,m) = 0.0d0
                            tnr_sacr2(i,j,k,m) = 0.0d0
                        enddo
                    enddo
                enddo
            enddo

            do m = 1, ntb_IN
                do k = 1, 45
                    do j = 1, ntb_r1
                        do i = 1, ntb_r
                            tpi_qrfz(i,j,k,m) = 0.0d0
                            tni_qrfz(i,j,k,m) = 0.0d0
                            tpg_qrfz(i,j,k,m) = 0.0d0
                            tnr_qrfz(i,j,k,m) = 0.0d0
                        enddo
                    enddo
                    do j = 1, nbc
                        do i = 1, ntb_c
                            tpi_qcfz(i,j,k,m) = 0.0d0
                            tni_qcfz(i,j,k,m) = 0.0d0
                        enddo
                    enddo
                enddo
            enddo

            do j = 1, ntb_i1
                do i = 1, ntb_i
                    tps_iaus(i,j) = 0.0d0
                    tni_iaus(i,j) = 0.0d0
                    tpi_ide(i,j) = 0.0d0
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
                        tnr_rev(i,j,k) = 0.0d0
                    enddo
                enddo
            enddo

            do k = 1, nbc
                do j = 1, ntb_c
                    do i = 1, nbc
                        tpc_wev(i,j,k) = 0.0d0
                        tnc_wev(i,j,k) = 0.0d0
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

!..Check that the look-up tables are available.
            if(.not. l_mp_tables) return

!..Read a static file containing CCN activation of aerosols. The
!.. data were created from a parcel model by Feingold & Heymsfield with
!.. further changes by Eidhammer and Kriedenweis.
            
! AAJ aerosol aware is off. Add this to the build table functionality
! for MPAS
            ! if (is_aerosol_aware) then
            !     call table_ccnAct
            !  endif

!..Collision efficiency between rain/snow and cloud water.
!     call physics_message('--- creating qc collision eff tables')
            call table_Efrw
            call table_Efsw

!..Drop evaporation.
!     call physics_message('--- creating rain evap table')
            call table_dropEvap

!..Rain collecting graupel & graupel collecting rain.
#if defined(mpas)
            call mpas_new_unit(mp_unit, unformatted = .true.)
#else
            mp_unit = 11
#endif
            open(unit=mp_unit,file='CCN_ACTIVATE.BIN',form='UNFORMATTED',status='OLD',action='READ', &
            iostat = istat)
            if(istat /= open_OK) &
            call physics_error_fatal('subroutine thompson_init: ' // &
            'failure opening CCN_ACTIVATE.Bin')
            read(mp_unit) tnccn_act
            close(unit=mp_unit)

            open(unit=mp_unit,file='MP_THOMPSON_QRacrQG_DATA.DBL',form='UNFORMATTED',status='OLD',action='READ', &
                iostat = istat)
            if(istat /= open_OK) &
                call physics_error_fatal('subroutine thompson_init: ' // &
                'failure opening MP_THOMPSON_QRacrQG.DBL')
            read(mp_unit) tcg_racg
            read(mp_unit) tmr_racg
            read(mp_unit) tcr_gacr
            ! read(mp_unit) tmg_gacr
            read(mp_unit) tnr_racg
            read(mp_unit) tnr_gacr
            close(unit=mp_unit)
!     write(0,*) '--- end read MP_THOMPSON_QRacrQG.DBL'
!     write(0,*) 'max tcg_racg =',maxval(tcg_racg)
!     write(0,*) 'min tcg_racg =',minval(tcg_racg)
!     write(0,*) 'max tmr_racg =',maxval(tmr_racg)
!     write(0,*) 'min tmr_racg =',minval(tmr_racg)
!     write(0,*) 'max tcr_gacr =',maxval(tcr_gacr)
!     write(0,*) 'min tcr_gacr =',minval(tcr_gacr)
!     write(0,*) 'max tmg_gacr =',maxval(tmg_gacr)
!     write(0,*) 'min tmg_gacr =',minval(tmg_gacr)
!     write(0,*) 'max tnr_racg =',maxval(tnr_racg)
!     write(0,*) 'min tnr_racg =',minval(tnr_racg)
!     write(0,*) 'max tnr_gacr =',maxval(tnr_gacr)
!     write(0,*) 'min tnr_gacr =',minval(tnr_gacr)

!..Rain collecting snow & snow collecting rain.
            open(unit=mp_unit,file='MP_THOMPSON_QRacrQS_DATA.DBL',form='UNFORMATTED',status='OLD',action='READ', &
                iostat=istat)
            if(istat /= open_OK) &
                call physics_error_fatal('subroutine thompson_init: ' // &
                'failure opening MP_THOMPSON_QRacrQS.DBL')
            read(mp_unit) tcs_racs1
            read(mp_unit) tmr_racs1
            read(mp_unit) tcs_racs2
            read(mp_unit) tmr_racs2
            read(mp_unit) tcr_sacr1
            read(mp_unit) tms_sacr1
            read(mp_unit) tcr_sacr2
            read(mp_unit) tms_sacr2
            read(mp_unit) tnr_racs1
            read(mp_unit) tnr_racs2
            read(mp_unit) tnr_sacr1
            read(mp_unit) tnr_sacr2
            close(unit=mp_unit)
!     write(0,*) '--- end read MP_THOMPSON_QRacrQS.DBL'
!     write(0,*) 'max tcs_racs1 =',maxval(tcs_racs1)
!     write(0,*) 'min tcs_racs1 =',minval(tcs_racs1)
!     write(0,*) 'max tmr_racs1 =',maxval(tmr_racs1)
!     write(0,*) 'min tmr_racs1 =',minval(tmr_racs1)
!     write(0,*) 'max tcs_racs2 =',maxval(tcs_racs2)
!     write(0,*) 'min tcs_racs2 =',minval(tcs_racs2)
!     write(0,*) 'max tmr_racs2 =',maxval(tmr_racs2)
!     write(0,*) 'min tmr_racs2 =',minval(tmr_racs2)
!     write(0,*) 'max tcr_sacr1 =',maxval(tcr_sacr1)
!     write(0,*) 'min tcr_sacr1 =',minval(tcr_sacr1)
!     write(0,*) 'max tms_sacr1 =',maxval(tms_sacr1)
!     write(0,*) 'min tms_sacr1 =',minval(tms_sacr1)
!     write(0,*) 'max tcr_sacr2 =',maxval(tcr_sacr2)
!     write(0,*) 'min tcr_sacr2 =',minval(tcr_sacr2)
!     write(0,*) 'max tms_sacr2 =',maxval(tms_sacr2)
!     write(0,*) 'min tms_sacr2 =',minval(tms_sacr2)
!     write(0,*) 'max tnr_racs1 =',maxval(tnr_racs1)
!     write(0,*) 'min tnr_racs1 =',minval(tnr_racs1)
!     write(0,*) 'max tnr_racs2 =',maxval(tnr_racs2)
!     write(0,*) 'min tnr_racs2 =',minval(tnr_racs2)
!     write(0,*) 'max tnr_sacr1 =',maxval(tnr_sacr1)
!     write(0,*) 'min tnr_sacr1 =',minval(tnr_sacr1)
!     write(0,*) 'max tnr_sacr2 =',maxval(tnr_sacr2)
!     write(0,*) 'min tnr_sacr2 =',minval(tnr_sacr2)

!..Cloud water and rain freezing (Bigg, 1953).
            open(unit=mp_unit,file='MP_THOMPSON_freezeH2O_DATA.DBL',form='UNFORMATTED',status='OLD',action='READ', &
                iostat=istat)
            if(istat /= open_OK) &
                call physics_error_fatal('subroutine thompson_init: ' // &
                'failure opening MP_THOMPSON_freezeH2O.DBL')
            read(mp_unit) tpi_qrfz
            read(mp_unit) tni_qrfz
            read(mp_unit) tpg_qrfz
            read(mp_unit) tnr_qrfz
            read(mp_unit) tpi_qcfz
            read(mp_unit) tni_qcfz
            close(unit=mp_unit)
!     write(0,*) '--- end read MP_THOMPSON_freezeH2O.DBL:'
!     write(0,*) 'max tpi_qrfz =',maxval(tpi_qrfz)
!     write(0,*) 'min tpi_qrfz =',minval(tpi_qrfz)
!     write(0,*) 'max tni_qrfz =',maxval(tni_qrfz)
!     write(0,*) 'min tni_qrfz =',minval(tni_qrfz)
!     write(0,*) 'max tpg_qrfz =',maxval(tpg_qrfz)
!     write(0,*) 'min tpg_qrfz =',minval(tpg_qrfz)
!     write(0,*) 'max tnr_qrfz =',maxval(tnr_qrfz)
!     write(0,*) 'min tnr_qrfz =',minval(tnr_qrfz)
!     write(0,*) 'max tpi_qcfz =',maxval(tpi_qcfz)
!     write(0,*) 'min tpi_qcfz =',minval(tpi_qcfz)
!     write(0,*) 'max tni_qcfz =',maxval(tni_qcfz)
!     write(0,*) 'min tni_qcfz =',minval(tni_qcfz)

!..Conversion of some ice mass into snow category.
            open(unit=mp_unit,file='MP_THOMPSON_QIautQS_DATA.DBL',form='UNFORMATTED',status='OLD',action='READ', &
                iostat=istat)
            if(istat /= open_OK) &
                call physics_error_fatal('subroutine thompson_init: ' // &
                'failure opening MP_THOMPSON_QIautQS.DBL')
            read(mp_unit) tpi_ide
            read(mp_unit) tps_iaus
            read(mp_unit) tni_iaus
            close(unit=mp_unit)
#if defined(mpas)
            call mpas_release_unit(mp_unit)
#endif
!     write(0,*) '--- end read MP_THOMPSON_QIautQS.DBL '
!     write(0,*) 'max tps_iaus =',maxval(tps_iaus)
!     write(0,*) 'min tps_iaus =',minval(tps_iaus)
!     write(0,*) 'max tni_iaus =',maxval(tni_iaus)
!     write(0,*) 'min tni_iaus =',minval(tni_iaus)

!..Initialize various constants for computing radar reflectivity.
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

        endif

    END SUBROUTINE thompson_init
!
!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
!ctrlL
!+---+-----------------------------------------------------------------+
!..This is a wrapper routine designed to transfer values from 3D to 1D.
!+---+-----------------------------------------------------------------+
    SUBROUTINE mp_gt_driver(qv, qc, qr, qi, qs, qg, qb, ni, nr, nc, ng,&
        nwfa, nifa, nwfa2d, nifa2d,                       &
        th, pii, p, w, dz, dt_in, itimestep,       &
        RAINNC, RAINNCV, &
        SNOWNC, SNOWNCV, &
        GRAUPELNC, GRAUPELNCV, SR, &
        refl_10cm, diagflag, do_radar_ref,      &
        re_cloud, re_ice, re_snow,              &
        has_reqc, has_reqi, has_reqs,           &
#if defined(mpas)
        ntc,muc,rainprod,evapprod, &
#endif
        ids,ide, jds,jde, kds,kde, &             ! domain dims
        ims,ime, jms,jme, kms,kme, &             ! memory dims
        its,ite, jts,jte, kts,kte)               ! tile dims

        implicit none

!..Subroutine arguments
        INTEGER, INTENT(IN):: ids,ide, jds,jde, kds,kde, &
            ims,ime, jms,jme, kms,kme, &
            its,ite, jts,jte, kts,kte
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
            qv, qc, qr, qi, qs, qg, ni, nr, th
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), OPTIONAL, INTENT(INOUT):: &
            nc, nwfa, nifa, qb, ng
        REAL, DIMENSION(ims:ime, jms:jme), OPTIONAL, INTENT(IN):: nwfa2d, nifa2d
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
            re_cloud, re_ice, re_snow
        INTEGER, INTENT(IN):: has_reqc, has_reqi, has_reqs
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
            pii, p, w, dz
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT):: &
            RAINNC, RAINNCV, SR
        REAL, DIMENSION(ims:ime, jms:jme), OPTIONAL, INTENT(INOUT)::      &
            SNOWNC, SNOWNCV, GRAUPELNC, GRAUPELNCV
#if defined(mpas)
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN):: &
            ntc,muc
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
            rainprod,evapprod
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT), OPTIONAL:: &
            refl_10cm
#else
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT)::       &
            refl_10cm
#endif
        REAL, INTENT(IN):: dt_in
        INTEGER, INTENT(IN):: itimestep

!..Local variables
        REAL, DIMENSION(kts:kte):: &
            qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, qb1d,     &
            ni1d, nr1d, nc1d, ng1d, nwfa1d, nifa1d,       &
            t1d, p1d, w1d, dz1d, rho, dBZ
        REAL, DIMENSION(kts:kte):: re_qc1d, re_qi1d, re_qs1d
#if defined(mpas)
        REAL, DIMENSION(kts:kte):: &
            rainprod1d, evapprod1d
#endif
        REAL, DIMENSION(its:ite, jts:jte):: pcp_ra, pcp_sn, pcp_gr, pcp_ic
        REAL:: dt, pptrain, pptsnow, pptgraul, pptice
        REAL:: qc_max, qr_max, qs_max, qi_max, qg_max, ni_max, nr_max
        REAL:: nwfa1
        REAL:: ygra1, zans1
        DOUBLE PRECISION:: lamg, lam_exp, lamr, N0_min, N0_exp
        INTEGER:: i, j, k
        INTEGER:: imax_qc,imax_qr,imax_qi,imax_qs,imax_qg,imax_ni,imax_nr
        INTEGER:: jmax_qc,jmax_qr,jmax_qi,jmax_qs,jmax_qg,jmax_ni,jmax_nr
        INTEGER:: kmax_qc,kmax_qr,kmax_qi,kmax_qs,kmax_qg,kmax_ni,kmax_nr
        INTEGER:: i_start, j_start, i_end, j_end
        LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
        INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref
        CHARACTER*256:: mp_debug

!+---+

        i_start = its
        j_start = jts
        i_end   = MIN(ite, ide-1)
        j_end   = MIN(jte, jde-1)

!..For idealized testing by developer.
!     if ( (ide-ids+1).gt.4 .and. (jde-jds+1).lt.4 .and.                &
!          ids.eq.its.and.ide.eq.ite.and.jds.eq.jts.and.jde.eq.jte) then
!        i_start = its + 2
!        i_end   = ite
!        j_start = jts
!        j_end   = jte
!     endif

        dt = dt_in

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
        do i = 1, 256
            mp_debug(i:i) = char(0)
        enddo

!     if (.NOT. is_aerosol_aware .AND. PRESENT(nc) .AND. PRESENT(nwfa)  &
!               .AND. PRESENT(nifa) .AND. PRESENT(nwfa2d)) then
!        write(mp_debug,*) 'WARNING, nc-nwfa-nifa-nwfa2d present but is_aerosol_aware is FALSE'
!        CALL wrf_debug(0, mp_debug)
!     endif

        j_loop:  do j = j_start, j_end
            i_loop:  do i = i_start, i_end

                pptrain = 0.
                pptsnow = 0.
                pptgraul = 0.
                pptice = 0.
                RAINNCV(i,j) = 0.
                IF ( PRESENT (snowncv) ) THEN
                    SNOWNCV(i,j) = 0.
                ENDIF
                IF ( PRESENT (graupelncv) ) THEN
                    GRAUPELNCV(i,j) = 0.
                ENDIF
                SR(i,j) = 0.

#if defined(mpas)
                Nt_c = ntc(i,j)
                mu_c = muc(i,j)
#endif
                do k = kts, kte
                    t1d(k) = th(i,k,j)*pii(i,k,j)
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
                enddo
                if (is_aerosol_aware) then
                    do k = kts, kte
                        nc1d(k) = nc(i,k,j)
                        nwfa1d(k) = nwfa(i,k,j)
                        nifa1d(k) = nifa(i,k,j)
                    enddo
                    nwfa1 = nwfa2d(i,j)
                else
                    do k = kts, kte
                        rho(k) = 0.622*p1d(k)/(R*t1d(k)*(qv1d(k)+0.622))
!AAJ TEST                        nc1d(k) = Nt_c/rho(k)
                        nc1d(k) = 100.e6/rho(k)
                        nwfa1d(k) = 11.1E6/rho(k)
                        nifa1d(k) = naIN1*0.01/rho(k)
                    enddo
                    nwfa1 = 11.1E6
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

                call mp_thompson(qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, qb1d, ni1d,     &
                    nr1d, nc1d, ng1d, nwfa1d, nifa1d, t1d, p1d, w1d, dz1d,  &
                    pptrain, pptsnow, pptgraul, pptice, &
#if defined(mpas)
                    rainprod1d, evapprod1d, &
#endif
                    kts, kte, dt, i, j)

                pcp_ra(i,j) = pptrain
                pcp_sn(i,j) = pptsnow
                pcp_gr(i,j) = pptgraul
                pcp_ic(i,j) = pptice
                RAINNCV(i,j) = pptrain + pptsnow + pptgraul + pptice
                RAINNC(i,j) = RAINNC(i,j) + pptrain + pptsnow + pptgraul + pptice
                IF ( PRESENT(snowncv) .AND. PRESENT(snownc) ) THEN
                    SNOWNCV(i,j) = pptsnow + pptice
                    SNOWNC(i,j) = SNOWNC(i,j) + pptsnow + pptice
                ENDIF
                IF ( PRESENT(graupelncv) .AND. PRESENT(graupelnc) ) THEN
                    GRAUPELNCV(i,j) = pptgraul
                    GRAUPELNC(i,j) = GRAUPELNC(i,j) + pptgraul
                ENDIF
                SR(i,j) = (pptsnow + pptgraul + pptice)/(RAINNCV(i,j)+1.e-12)



!..Reset lowest model level to initial state aerosols (fake sfc source).
!.. Changed 13 May 2013 to fake emissions in which nwfa2d is aerosol
!.. number tendency (number per kg per second).
                if (is_aerosol_aware) then
!-GT        nwfa1d(kts) = nwfa1
                    nwfa1d(kts) = nwfa1d(kts) + nwfa2d(i,j)*dt_in
                    nifa1d(kts) = nifa1d(kts) + nifa2d(i,j)*dt_in
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
                    nc(i,k,j) = nc1d(k) !AAJ testing nc output
                    qv(i,k,j) = qv1d(k)
                    qc(i,k,j) = qc1d(k)
                    qi(i,k,j) = qi1d(k)
                    qr(i,k,j) = qr1d(k)
                    qs(i,k,j) = qs1d(k)
                    qg(i,k,j) = qg1d(k)
                    ni(i,k,j) = ni1d(k)
                    nr(i,k,j) = nr1d(k)
                    th(i,k,j) = t1d(k)/pii(i,k,j)
#if defined(mpas)
                    rainprod(i,k,j) = rainprod1d(k)
                    evapprod(i,k,j) = evapprod1d(k)
#endif
                    if (qc1d(k) .gt. qc_max) then
                        imax_qc = i
                        jmax_qc = j
                        kmax_qc = k
                        qc_max = qc1d(k)
                    elseif (qc1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative qc ', qc1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (qr1d(k) .gt. qr_max) then
                        imax_qr = i
                        jmax_qr = j
                        kmax_qr = k
                        qr_max = qr1d(k)
                    elseif (qr1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative qr ', qr1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (nr1d(k) .gt. nr_max) then
                        imax_nr = i
                        jmax_nr = j
                        kmax_nr = k
                        nr_max = nr1d(k)
                    elseif (nr1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative nr ', nr1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (qs1d(k) .gt. qs_max) then
                        imax_qs = i
                        jmax_qs = j
                        kmax_qs = k
                        qs_max = qs1d(k)
                    elseif (qs1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative qs ', qs1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (qi1d(k) .gt. qi_max) then
                        imax_qi = i
                        jmax_qi = j
                        kmax_qi = k
                        qi_max = qi1d(k)
                    elseif (qi1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative qi ', qi1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (qg1d(k) .gt. qg_max) then
                        imax_qg = i
                        jmax_qg = j
                        kmax_qg = k
                        qg_max = qg1d(k)
                    elseif (qg1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative qg ', qg1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (ni1d(k) .gt. ni_max) then
                        imax_ni = i
                        jmax_ni = j
                        kmax_ni = k
                        ni_max = ni1d(k)
                    elseif (ni1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative ni ', ni1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                    endif
                    if (qv1d(k) .lt. 0.0) then
                        write(mp_debug,*) 'WARNING, negative qv ', qv1d(k),        &
                            ' at i,j,k=', i,j,k
!            CALL wrf_debug(150, mp_debug)
                        if (k.lt.kte-2 .and. k.gt.kts+1) then
                            write(mp_debug,*) '   below and above are: ', qv(i,k-1,j), qv(i,k+1,j)
!               CALL wrf_debug(150, mp_debug)
                            qv(i,k,j) = MAX(1.E-7, 0.5*(qv(i,k-1,j) + qv(i,k+1,j)))
                        else
                            qv(i,k,j) = 1.E-7
                        endif
                    endif
                enddo

!        IF ( PRESENT (diagflag) ) THEN
!        if (diagflag .and. do_radar_ref == 1) then
!         call calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d,       &
!                     t1d, p1d, dBZ, kts, kte, i, j)
!         do k = kts, kte
!            refl_10cm(i,k,j) = MAX(-35., dBZ(k))
!         enddo
!        endif
!        ENDIF

                IF (has_reqc.ne.0 .and. has_reqi.ne.0 .and. has_reqs.ne.0) THEN
                    do k = kts, kte
                        re_qc1d(k) = 2.49E-6
                        re_qi1d(k) = 4.99E-6
                        re_qs1d(k) = 9.99E-6
                    enddo
                    call calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d,  &
                        re_qc1d, re_qi1d, re_qs1d, kts, kte)
                    do k = kts, kte
                        re_cloud(i,k,j) = MAX(2.49E-6, MIN(re_qc1d(k), 50.E-6))
                        re_ice(i,k,j)   = MAX(4.99E-6, MIN(re_qi1d(k), 125.E-6))
                        re_snow(i,k,j)  = MAX(9.99E-6, MIN(re_qs1d(k), 999.E-6))
                    enddo
                ENDIF

            enddo i_loop
        enddo j_loop

! DEBUG - GT
        write(mp_debug,'(a,7(a,e13.6,1x,a,i3,a,i3,a,i3,a,1x))') 'MP-GT:', &
            'qc: ', qc_max, '(', imax_qc, ',', jmax_qc, ',', kmax_qc, ')', &
            'qr: ', qr_max, '(', imax_qr, ',', jmax_qr, ',', kmax_qr, ')', &
            'qi: ', qi_max, '(', imax_qi, ',', jmax_qi, ',', kmax_qi, ')', &
            'qs: ', qs_max, '(', imax_qs, ',', jmax_qs, ',', kmax_qs, ')', &
            'qg: ', qg_max, '(', imax_qg, ',', jmax_qg, ',', kmax_qg, ')', &
            'ni: ', ni_max, '(', imax_ni, ',', jmax_ni, ',', kmax_ni, ')', &
            'nr: ', nr_max, '(', imax_nr, ',', jmax_nr, ',', kmax_nr, ')'
!     CALL wrf_debug(150, mp_debug)
! END DEBUG - GT

        do i = 1, 256
            mp_debug(i:i) = char(0)
        enddo

    END SUBROUTINE mp_gt_driver

!+---+-----------------------------------------------------------------+
! !ctrlL
!+---+-----------------------------------------------------------------+
!..Creation of the lookup tables and support functions found below here.
!+---+-----------------------------------------------------------------+
!..Rain collecting graupel (and inverse).  Explicit CE integration.
!+---+-----------------------------------------------------------------+

    subroutine qr_acr_qg(NRHGtable)
        implicit none

        INTEGER, INTENT(IN) ::NRHGtable

        !..Local variables
        INTEGER:: i, j, k, m, n, n2, n3, idx_bg
        INTEGER:: km, km_s, km_e
        DOUBLE PRECISION, DIMENSION(nbg):: N_g
        DOUBLE PRECISION, DIMENSION(nbg,NRHGtable):: vg
        DOUBLE PRECISION, DIMENSION(nbr):: vr, N_r
        DOUBLE PRECISION:: N0_r, N0_g, lam_exp, lamg, lamr
        DOUBLE PRECISION:: massg, massr, dvg, dvr, t1, t2, z1, z2, y1, y2

        !+---+

        do n2 = 1, nbr
            !        vr(n2) = av_r*Dr(n2)**bv_r * DEXP(-fv_r*Dr(n2))
            vr(n2) = -0.1021 + 4.932E3*Dr(n2) - 0.9551E6*Dr(n2)*Dr(n2)     &
                + 0.07934E9*Dr(n2)*Dr(n2)*Dr(n2)                           &
                - 0.002362E12*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
        enddo

        do n3 = 1, NRHGtable
            do n = 1, nbg
                idx_bg = n3
                vg(n,n3) = av_g(idx_bg)*Dg(n)**bv_g(idx_bg)
            enddo
        enddo

        km_s = 0
        km_e = ntb_r*ntb_r1 - 1

        do km = km_s, km_e
            m = km / ntb_r1 + 1
            k = mod( km , ntb_r1 ) + 1

            lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
            lamr = lam_exp * (crg(3)*org2*org1)**obmr
            N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
            do n2 = 1, nbr
                N_r(n2) = N0_r*Dr(n2)**mu_r *DEXP(-lamr*Dr(n2))*dtr(n2)
            enddo

            do n3 = 1, NRHGtable
                idx_bg = n3

                do j = 1, ntb_g
                    do i = 1, ntb_g1
                        lam_exp = (N0g_exp(i)*am_g(idx_bg)*cgg(1,1)/r_g(j))**oge1
                        lamg = lam_exp * (cgg(3,1)*ogg2*ogg1)**obmg
                        N0_g = N0g_exp(i)/(cgg(2,1)*lam_exp) * lamg**cge(2,1)
                        do n = 1, nbg
                            N_g(n) = N0_g*Dg(n)**mu_g * DEXP(-lamg*Dg(n))*dtg(n)
                        enddo

                        t1 = 0.0d0
                        t2 = 0.0d0
                        z1 = 0.0d0
                        z2 = 0.0d0
                        y1 = 0.0d0
                        y2 = 0.0d0
                        do n2 = 1, nbr
                            massr = am_r * Dr(n2)**bm_r
                            do n = 1, nbg
                                massg = am_g(idx_bg) * Dg(n)**bm_g

                                dvg = 0.5d0*((vr(n2) - vg(n,n3)) + DABS(vr(n2)-vg(n,n3)))
                                dvr = 0.5d0*((vg(n,n3) - vr(n2)) + DABS(vg(n,n3)-vr(n2)))

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
97                          continue
                        enddo
                        tcg_racg(i,j,n3,k,m) = t1
                        tmr_racg(i,j,n3,k,m) = DMIN1(z1, r_r(m)*1.0d0)
                        tcr_gacr(i,j,n3,k,m) = t2
                        ! tmg_gacr(i,j,n3,k,m) = DMIN1(z2, r_g(j)*1.0d0)
                        tnr_racg(i,j,n3,k,m) = y1
                        tnr_gacr(i,j,n3,k,m) = y2
                    enddo
                enddo
            enddo
        enddo

    end subroutine qr_acr_qg
!+---+-----------------------------------------------------------------+
!ctrlL
!+---+-----------------------------------------------------------------+
!..Rain collecting snow (and inverse).  Explicit CE integration.
!+---+-----------------------------------------------------------------+

    subroutine qr_acr_qs

        implicit none

!..Local variables
        INTEGER:: i, j, k, m, n, n2
        INTEGER:: km, km_s, km_e
        DOUBLE PRECISION, DIMENSION(nbr):: vr, D1, N_r
        DOUBLE PRECISION, DIMENSION(nbs):: vs, N_s
        DOUBLE PRECISION:: loga_, a_, b_, second, M0, M2, M3, Mrat, oM3
        DOUBLE PRECISION:: N0_r, lam_exp, lamr, slam1, slam2
        DOUBLE PRECISION:: dvs, dvr, masss, massr
        DOUBLE PRECISION:: t1, t2, t3, t4, z1, z2, z3, z4
        DOUBLE PRECISION:: y1, y2, y3, y4

!+---+

        do n2 = 1, nbr
!        vr(n2) = av_r*Dr(n2)**bv_r * DEXP(-fv_r*Dr(n2))
            vr(n2) = -0.1021 + 4.932E3*Dr(n2) - 0.9551E6*Dr(n2)*Dr(n2)     &
                + 0.07934E9*Dr(n2)*Dr(n2)*Dr(n2)                           &
                - 0.002362E12*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
            D1(n2) = (vr(n2)/av_s)**(1./bv_s)
        enddo
        do n = 1, nbs
            vs(n) = 1.5*av_s*Ds(n)**bv_s * DEXP(-fv_s*Ds(n))
        enddo

        km_s = 0
        km_e = ntb_r*ntb_r1 - 1

        do km = km_s, km_e
            m = km / ntb_r1 + 1
            k = mod( km , ntb_r1 ) + 1

            lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
            lamr = lam_exp * (crg(3)*org2*org1)**obmr
            N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
            do n2 = 1, nbr
                N_r(n2) = N0_r*Dr(n2)**mu_r * DEXP(-lamr*Dr(n2))*dtr(n2)
            enddo

            do j = 1, ntb_t
                do i = 1, ntb_s

!..From the bm_s moment, compute plus one moment.  If we are not
!.. using bm_s=2, then we must transform to the pure 2nd moment
!.. (variable called "second") and then to the bm_s+1 moment.

                    M2 = r_s(i)*oams *1.0d0
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
                        N_s(n) = Mrat*(Kap0*DEXP(-slam1*Ds(n)) &
                            + Kap1*M0*Ds(n)**mu_s * DEXP(-slam2*Ds(n)))*dts(n)
                    enddo

                    t1 = 0.0d0
                    t2 = 0.0d0
                    t3 = 0.0d0
                    t4 = 0.0d0
                    z1 = 0.0d0
                    z2 = 0.0d0
                    z3 = 0.0d0
                    z4 = 0.0d0
                    y1 = 0.0d0
                    y2 = 0.0d0
                    y3 = 0.0d0
                    y4 = 0.0d0
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
                    tmr_racs1(i,j,k,m) = DMIN1(z1, r_r(m)*1.0d0)
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

    end subroutine qr_acr_qs
!+---+-----------------------------------------------------------------+
!ctrlL
!+---+-----------------------------------------------------------------+
!..This is a literal adaptation of Bigg (1954) probability of drops of
!..a particular volume freezing.  Given this probability, simply freeze
!..the proportion of drops summing their masses.
!+---+-----------------------------------------------------------------+

    subroutine freezeH2O

        implicit none

        !..Local variables
        INTEGER:: i, j, k, m, n, n2
        INTEGER:: km, km_s, km_e
        DOUBLE PRECISION :: N_r, N_c
        DOUBLE PRECISION, DIMENSION(nbr):: massr
        DOUBLE PRECISION, DIMENSION(nbc):: massc
        DOUBLE PRECISION:: sum1, sum2, sumn1, sumn2, &
            prob, vol, Texp, orho_w, &
            lam_exp, lamr, N0_r, lamc, N0_c, y
        INTEGER:: nu_c
        REAL:: T_adjust

        !+---+

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
            Texp = DEXP( REAL(k,KIND=dp) - T_adjust*1.0D0 ) - 1.0D0
            do j = 1, ntb_r1
                do i = 1, ntb_r
                    lam_exp = (N0r_exp(j)*am_r*crg(1)/r_r(i))**ore1
                    lamr = lam_exp * (crg(3)*org2*org1)**obmr
                    N0_r = N0r_exp(j)/(crg(2)*lam_exp) * lamr**cre(2)
                    sum1 = 0.0d0
                    sum2 = 0.0d0
                    sumn1 = 0.0d0
                    sumn2 = 0.0d0
                    do n2 = nbr, 1, -1
                        N_r = N0_r*Dr(n2)**mu_r*DEXP(-lamr*Dr(n2))*dtr(n2)
                        vol = massr(n2)*orho_w
                        prob = MAX(0.0D0, 1.0D0 - DEXP(-120.0D0*vol*5.2D-4 * Texp))
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

            do j = 1, nbc
                nu_c = MIN(15, NINT(1000.E6/t_Nc(j)) + 2)
                do i = 1, ntb_c
                    lamc = (t_Nc(j)*am_r* ccg(2,nu_c) * ocg1(nu_c) / r_c(i))**obmr
                    N0_c = t_Nc(j)*ocg1(nu_c) * lamc**cce(1,nu_c)
                    sum1 = 0.0d0
                    sumn2 = 0.0d0
                    do n = nbc, 1, -1
                        vol = massc(n)*orho_w
                        prob = MAX(0.0D0, 1.0D0 - DEXP(-120.0D0*vol*5.2D-4 * Texp))
                        N_c = N0_c*Dc(n)**nu_c*EXP(-lamc*Dc(n))*dtc(n)
                        sumn2 = MIN(t_Nc(j), sumn2 + prob*N_c)
                        sum1 = sum1 + prob*N_c*massc(n)
                        if (sum1 .ge. r_c(i)) EXIT
                    enddo
                    tpi_qcfz(i,j,k,m) = sum1
                    tni_qcfz(i,j,k,m) = sumn2
                enddo
            enddo
        enddo

    end subroutine freezeH2O
!+---+-----------------------------------------------------------------+
!ctrlL
!+---+-----------------------------------------------------------------+
!..Cloud ice converting to snow since portion greater than min snow
!.. size.  Given cloud ice content (kg/m**3), number concentration
!.. (#/m**3) and gamma shape parameter, mu_i, break the distrib into
!.. bins and figure out the mass/number of ice with sizes larger than
!.. D0s.  Also, compute incomplete gamma function for the integration
!.. of ice depositional growth from diameter=0 to D0s.  Amount of
!.. ice depositional growth is this portion of distrib while larger
!.. diameters contribute to snow growth (as in Harrington et al. 1995).
!+---+-----------------------------------------------------------------+

    subroutine qi_aut_qs

        implicit none

!..Local variables
        INTEGER:: i, j, n2
        DOUBLE PRECISION, DIMENSION(nbi):: N_i
        DOUBLE PRECISION:: N0_i, lami, Di_mean, t1, t2
        REAL:: xlimit_intg

!+---+

        do j = 1, ntb_i1
            do i = 1, ntb_i
                lami = (am_i*cig(2)*oig1*Nt_i(j)/r_i(i))**obmi
                Di_mean = (bm_i + mu_i + 1.) / lami
                N0_i = Nt_i(j)*oig1 * lami**cie(1)
                t1 = 0.0d0
                t2 = 0.0d0
                if (SNGL(Di_mean) .gt. 5.*D0s) then
                    t1 = r_i(i)
                    t2 = Nt_i(j)
                    tpi_ide(i,j) = 0.0D0
                elseif (SNGL(Di_mean) .lt. D0i) then
                    t1 = 0.0D0
                    t2 = 0.0D0
                    tpi_ide(i,j) = 1.0D0
                else
                    xlimit_intg = lami*D0s
                    tpi_ide(i,j) = GAMMP(mu_i+2.0, xlimit_intg) * 1.0D0
                    do n2 = 1, nbi
                        N_i(n2) = N0_i*Di(n2)**mu_i * DEXP(-lami*Di(n2))*dti(n2)
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
!..Variable collision efficiency for rain collecting cloud water using
!.. method of Beard and Grover, 1974 if a/A less than 0.25; otherwise
!.. uses polynomials to get close match of Pruppacher & Klett Fig 14-9.
!+---+-----------------------------------------------------------------+

    subroutine table_Efrw

        implicit none

!..Local variables
        DOUBLE PRECISION:: vtr, stokes, reynolds, Ef_rw
        DOUBLE PRECISION:: p, yc0, F, G, H, z, K0, X
        INTEGER:: i, j

        do j = 1, nbc
            do i = 1, nbr
                Ef_rw = 0.0
                p = Dc(j)/Dr(i)
                if (Dr(i).lt.50.E-6 .or. Dc(j).lt.3.E-6) then
                    t_Efrw(i,j) = 0.0
                elseif (p.gt.0.25) then
                    X = Dc(j)*1.D6
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

                    F = DLOG(reynolds)
                    G = -0.1007D0 - 0.358D0*F + 0.0261D0*F*F
                    K0 = DEXP(G)
                    z = DLOG(stokes/(K0+1.D-15))
                    H = 0.1465D0 + 1.302D0*z - 0.607D0*z*z + 0.293D0*z*z*z
                    yc0 = 2.0D0/PI * ATAN(H)
                    Ef_rw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))

                endif

                t_Efrw(i,j) = MAX(0.0, MIN(SNGL(Ef_rw), 0.95))

            enddo
        enddo

    end subroutine table_Efrw
!ctrlL
!+---+-----------------------------------------------------------------+
!..Variable collision efficiency for snow collecting cloud water using
!.. method of Wang and Ji, 2000 except equate melted snow diameter to
!.. their "effective collision cross-section."
!+---+-----------------------------------------------------------------+

    subroutine table_Efsw

        implicit none

!..Local variables
        DOUBLE PRECISION:: Ds_m, vts, vtc, stokes, reynolds, Ef_sw
        DOUBLE PRECISION:: p, yc0, F, G, H, z, K0
        INTEGER:: i, j

        do j = 1, nbc
            vtc = 1.19D4 * (1.0D4*Dc(j)*Dc(j)*0.25D0)
            do i = 1, nbs
                vts = av_s*Ds(i)**bv_s * DEXP(-fv_s*Ds(i)) - vtc
                Ds_m = (am_s*Ds(i)**bm_s / am_r)**obmr
                p = Dc(j)/Ds_m
                if (p.gt.0.25 .or. Ds(i).lt.D0s .or. Dc(j).lt.6.E-6 &
                    .or. vts.lt.1.E-3) then
                    t_Efsw(i,j) = 0.0
                else
                    stokes = Dc(j)*Dc(j)*vts*rho_w2/(9.*1.718E-5*Ds_m)
                    reynolds = 9.*stokes/(p*p*rho_w2)

                    F = DLOG(reynolds)
                    G = -0.1007D0 - 0.358D0*F + 0.0261D0*F*F
                    K0 = DEXP(G)
                    z = DLOG(stokes/(K0+1.D-15))
                    H = 0.1465D0 + 1.302D0*z - 0.607D0*z*z + 0.293D0*z*z*z
                    yc0 = 2.0D0/PI * ATAN(H)
                    Ef_sw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))

                    t_Efsw(i,j) = MAX(0.0, MIN(SNGL(Ef_sw), 0.95))
                endif

            enddo
        enddo

    end subroutine table_Efsw
!ctrlL

!+---+-----------------------------------------------------------------+
!..Integrate rain size distribution from zero to D-star to compute the
!.. number of drops smaller than D-star that evaporate in a single
!.. timestep.  Drops larger than D-star dont evaporate entirely so do
!.. not affect number concentration.
!+---+-----------------------------------------------------------------+

    subroutine table_dropEvap

        implicit none

!..Local variables
        INTEGER:: i, j, k, n
        DOUBLE PRECISION, DIMENSION(nbc):: N_c, massc
        DOUBLE PRECISION:: summ, summ2, lamc, N0_c
        INTEGER:: nu_c
!      DOUBLE PRECISION:: Nt_r, N0, lam_exp, lam
!      REAL:: xlimit_intg

        do n = 1, nbc
            massc(n) = am_r*Dc(n)**bm_r
        enddo

        do k = 1, nbc
            nu_c = MIN(15, NINT(1000.E6/t_Nc(k)) + 2)
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
!         Dr_star = DSQRT(-2.D0*DT * t1_evap/(2.*PI) &
!                 * 0.78*4.*diffu(k)*xsat*rvs/rho_w)
!         idx_d = NINT(1.0 + FLOAT(nbr) * DLOG(Dr_star/D0r)             &
!               / DLOG(Dr(nbr)/D0r))
!         idx_d = MAX(1, MIN(idx_d, nbr))
!
!         nir = NINT(ALOG10(rr(k)))
!         do nn = nir-1, nir+1
!            n = nn
!            if ( (rr(k)/10.**nn).ge.1.0 .and. &
!                 (rr(k)/10.**nn).lt.10.0) goto 154
!         enddo
!154      continue
!         idx_r = INT(rr(k)/10.**n) + 10*(n-nir2) - (n-nir2)
!         idx_r = MAX(1, MIN(idx_r, ntb_r))
!
!         lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
!         lam_exp = lamr * (crg(3)*org2*org1)**bm_r
!         N0_exp = org1*rr(k)/am_r * lam_exp**cre(1)
!         nir = NINT(DLOG10(N0_exp))
!         do nn = nir-1, nir+1
!            n = nn
!            if ( (N0_exp/10.**nn).ge.1.0 .and. &
!                 (N0_exp/10.**nn).lt.10.0) goto 155
!         enddo
!155      continue
!         idx_r1 = INT(N0_exp/10.**n) + 10*(n-nir3) - (n-nir3)
!         idx_r1 = MAX(1, MIN(idx_r1, ntb_r1))
!
!         pnr_rev(k) = MIN(nr(k)*odts, SNGL(tnr_rev(idx_d,idx_r1,idx_r) &   ! RAIN2M
!                    * odts))

    end subroutine table_dropEvap
!
!ctrlL
#if !defined (mpas)
!+---+-----------------------------------------------------------------+
!..Fill the table of CCN activation data created from parcel model run
!.. by Trude Eidhammer with inputs of aerosol number concentration,
!.. vertical velocity, temperature, lognormal mean aerosol radius, and
!.. hygroscopicity, kappa.  The data are read from external file and
!.. contain activated fraction of CCN for given conditions.
!+---+-----------------------------------------------------------------+

    subroutine table_ccnAct

        USE module_domain
        USE module_dm
        implicit none

        LOGICAL, EXTERNAL:: wrf_dm_on_monitor

!..Local variables
        INTEGER:: iunit_mp_th1, i
        LOGICAL:: opened
        CHARACTER*64 errmess

        iunit_mp_th1 = -1
        IF ( wrf_dm_on_monitor() ) THEN
            DO i = 20,99
                INQUIRE ( i , OPENED = opened )
                IF ( .NOT. opened ) THEN
                    iunit_mp_th1 = i
                    GOTO 2010
                ENDIF
            ENDDO
2010        CONTINUE
        ENDIF
#if defined(DM_PARALLEL) && !defined(STUBMPI)
        CALL wrf_dm_bcast_bytes ( iunit_mp_th1 , IWORDSIZE )
#endif
        IF ( iunit_mp_th1 < 0 ) THEN
            CALL wrf_error_fatal ( 'module_mp_thompson: table_ccnAct: '//   &
                'Can not find unused fortran unit to read in lookup table.')
        ENDIF

        IF ( wrf_dm_on_monitor() ) THEN
            WRITE(errmess, '(A,I2)') 'module_mp_thompson: opening CCN_ACTIVATE.BIN on unit ',iunit_mp_th1
!       CALL wrf_debug(150, errmess)
            OPEN(iunit_mp_th1,FILE='CCN_ACTIVATE.BIN',                      &
                FORM='UNFORMATTED',STATUS='OLD',ERR=9009)
        ENDIF

#define DM_BCAST_MACRO(A) CALL wrf_dm_bcast_bytes(A, size(A)*R4SIZE)

        IF ( wrf_dm_on_monitor() ) READ(iunit_mp_th1,ERR=9010) tnccn_act
#if defined(DM_PARALLEL) && !defined(STUBMPI)
        DM_BCAST_MACRO(tnccn_act)
#endif


        RETURN
9009    CONTINUE
        WRITE( errmess , '(A,I2)' ) 'module_mp_thompson: error opening CCN_ACTIVATE.BIN on unit ',iunit_mp_th1
        CALL wrf_error_fatal(errmess)
        RETURN
9010    CONTINUE
        WRITE( errmess , '(A,I2)' ) 'module_mp_thompson: error reading CCN_ACTIVATE.BIN on unit ',iunit_mp_th1
        CALL wrf_error_fatal(errmess)

    end subroutine table_ccnAct
#endif

!+---+-----------------------------------------------------------------+

!+---+-----------------------------------------------------------------+
END MODULE module_mp_thompson
!+---+-----------------------------------------------------------------+

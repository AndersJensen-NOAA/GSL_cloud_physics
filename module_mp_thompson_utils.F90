module module_mp_thompson_utils

#if defined(mpas)
    use mpas_kind_types, only: wp => RKIND, sp => R4KIND, dp => R8KIND
#else
#ifndef CCPP
#define CCPP
#endif
#endif

#ifdef(CCPP)
    use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
#endif

    use module_mp_thompson_params
    use module_mp_radar
contains
!+---+-----------------------------------------------------------------+
!..Compute _radiation_ effective radii of cloud water, ice, and snow.
!.. These are entirely consistent with microphysics assumptions, not
!.. constant or otherwise ad hoc as is internal to most radiation
!.. schemes.  Since only the smallest snowflakes should impact
!.. radiation, compute from first portion of complicated Field number
!.. distribution, not the second part, which is the larger sizes.
!+---+-----------------------------------------------------------------+
    subroutine calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d,   &
    &                re_qc1d, re_qi1d, re_qs1d, lsm, kts, kte)

        IMPLICIT NONE

        !..Sub arguments
        INTEGER, INTENT(IN):: kts, kte
        INTEGER, OPTIONAL, INTENT(IN):: lsm
        REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
        &                    t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d
        REAL, DIMENSION(kts:kte), INTENT(INOUT):: re_qc1d, re_qi1d, re_qs1d
        !..Local variables
        INTEGER:: k
        REAL, DIMENSION(kts:kte):: rho, rc, nc, ri, ni, rs
        REAL:: smo2, smob, smoc
        REAL:: tc0, loga_, a_, b_
        DOUBLE PRECISION:: lamc, lami
        LOGICAL:: has_qc, has_qi, has_qs
        INTEGER:: inu_c
        real, dimension(15), parameter:: g_ratio = (/24,60,120,210,336,   &
        &                504,720,990,1320,1716,2184,2730,3360,4080,4896/)

        has_qc = .false.
        has_qi = .false.
        has_qs = .false.

        do k = kts, kte
            rho(k) = 0.622*p1d(k)/(R*t1d(k)*(qv1d(k)+0.622))
            rc(k) = MAX(R1, qc1d(k)*rho(k))
            nc(k) = MAX(2., MIN(nc1d(k)*rho(k), Nt_c_max))

            if (.NOT. (is_aerosol_aware .or. merra2_aerosol_aware)) then
                nc(k) = Nt_c_l
                if(present(lsml)) then
                    if( lsml == 1) then
                        nc(k) = Nt_c_l
                    else
                        nc(k) = Nt_c_o
                    endif
                endif
            endif
            ! if (.NOT. is_aerosol_aware) nc(k) = Nt_c
            if (rc(k).gt.R1 .and. nc(k).gt.R2) has_qc = .true.
            ri(k) = MAX(R1, qi1d(k)*rho(k))
            ni(k) = MAX(R2, ni1d(k)*rho(k))
            if (ri(k).gt.R1 .and. ni(k).gt.R2) has_qi = .true.
            rs(k) = MAX(R1, qs1d(k)*rho(k))
            if (rs(k).gt.R1) has_qs = .true.
        enddo

        if (has_qc) then
            do k = kts, kte
                if (rc(k).le.R1 .or. nc(k).le.R2) CYCLE
                if (nc(k).lt.100) then
                    inu_c = 15
                elseif (nc(k).gt.1.E10) then
                    inu_c = 2
                else
                    inu_c = MIN(15, NINT(1000.E6/nc(k)) + 2)
                endif
                lamc = (nc(k)*am_r*g_ratio(inu_c)/rc(k))**obmr
                re_qc1d(k) = MAX(2.51E-6, MIN(SNGL(0.5D0 * DBLE(3.+inu_c)/lamc), 50.E-6))
            enddo
        endif

        if (has_qi) then
            do k = kts, kte
                if (ri(k).le.R1 .or. ni(k).le.R2) CYCLE
                lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
                re_qi1d(k) = MAX(2.51E-6, MIN(SNGL(0.5D0 * DBLE(3.+mu_i)/lami), 125.E-6))
            enddo
        endif

        if (has_qs) then
            do k = kts, kte
                if (rs(k).le.R1) CYCLE
                tc0 = MIN(-0.1, t1d(k)-273.15)
                smob = rs(k)*oams

                !..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
                !.. then we must compute actual 2nd moment and use as reference.
                if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
                    smo2 = smob
                else
                    loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
                    &         + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
                    &         + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
                    &         + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
                    &         + sa(10)*bm_s*bm_s*bm_s
                    a_ = 10.0**loga_
                    b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
                    &         + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
                    &         + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
                    &         + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
                    &         + sb(10)*bm_s*bm_s*bm_s
                    smo2 = (smob/a_)**(1./b_)
                endif
                !..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
                loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
                &         + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
                &         + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
                &         + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
                &         + sa(10)*cse(1)*cse(1)*cse(1)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
                &        + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
                &        + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
                &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
                smoc = a_ * smo2**b_
                re_qs1d(k) = MAX(5.01E-6, MIN(0.5*(smoc/smob), 999.E-6))
            enddo
        endif

    end subroutine calc_effectRad

!>\ingroup aathompson
!! Compute radar reflectivity assuming 10 cm wavelength radar and using
!! Rayleigh approximation.  Only complication is melted snow/graupel
!! which we treat as water-coated ice spheres and use Uli Blahak's
!! library of routines.  The meltwater fraction is simply the amount
!! of frozen species remaining from what initially existed at the
!! melting level interface.

    subroutine calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d, &
        t1d, p1d, dBZ, kts, kte, ii, jj, melti,       &
        vt_dBZ, first_time_step)

        IMPLICIT NONE

!..Sub arguments
        INTEGER, INTENT(IN):: kts, kte, ii, jj
!   REAL, INTENT(IN):: rand1
        REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
            qv1d, qc1d, qr1d, nr1d, qs1d, qg1d, t1d, p1d
        REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ
        REAL, DIMENSION(kts:kte), OPTIONAL, INTENT(INOUT):: vt_dBZ
        LOGICAL, OPTIONAL, INTENT(IN) :: first_time_step

!..Local variables
        LOGICAL :: do_vt_dBZ
        LOGICAL :: allow_wet_graupel
        LOGICAL :: allow_wet_snow
        REAL, DIMENSION(kts:kte):: temp, pres, qv, rho, rhof
        REAL, DIMENSION(kts:kte):: rc, rr, nr, rs, rg

        DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilamg, N0_r, N0_g
        REAL, DIMENSION(kts:kte):: mvd_r
        REAL, DIMENSION(kts:kte):: smob, smo2, smoc, smoz
        REAL:: oM3, M0, Mrat, slam1, slam2, xDs
        REAL:: ils1, ils2, t1_vts, t2_vts, t3_vts, t4_vts
        REAL:: vtr_dbz_wt, vts_dbz_wt, vtg_dbz_wt

        REAL, DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel

        DOUBLE PRECISION:: N0_exp, N0_min, lam_exp, lamr, lamg
        REAL:: a_, b_, loga_, tc0, SR
        DOUBLE PRECISION:: fmelt_s, fmelt_g

        INTEGER:: i, k, k_0, kbot, n
        LOGICAL, OPTIONAL, INTENT(IN):: melti
        LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg

        DOUBLE PRECISION:: cback, x, eta, f_d
        REAL:: xslw1, ygra1, zans1

!+---+
        if (present(vt_dBZ) .and. present(first_time_step)) then
            do_vt_dBZ = .true.
            if (first_time_step) then
!           no bright banding, to be consistent with hydrometeor retrieval in GSI
                allow_wet_snow = .false.
            else
                allow_wet_snow = .true.
            endif
            allow_wet_graupel = .false.
        else
            do_vt_dBZ = .false.
            allow_wet_snow = .true.
            allow_wet_graupel = .false.
        endif

        do k = kts, kte
            dBZ(k) = -35.0
        enddo

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
        do k = kts, kte
            temp(k) = t1d(k)
            qv(k) = MAX(1.E-10, qv1d(k))
            pres(k) = p1d(k)
            rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
            rhof(k) = SQRT(RHO_NOT/rho(k))
            rc(k) = MAX(R1, qc1d(k)*rho(k))
            if (qr1d(k) .gt. R1) then
                rr(k) = qr1d(k)*rho(k)
                nr(k) = MAX(R2, nr1d(k)*rho(k))
                lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
                ilamr(k) = 1./lamr
                N0_r(k) = nr(k)*org2*lamr**cre(2)
                mvd_r(k) = (3.0 + mu_r + 0.672) * ilamr(k)
                L_qr(k) = .true.
            else
                rr(k) = R1
                nr(k) = R1
                mvd_r(k) = 50.E-6
                L_qr(k) = .false.
            endif
            if (qs1d(k) .gt. R2) then
                rs(k) = qs1d(k)*rho(k)
                L_qs(k) = .true.
            else
                rs(k) = R1
                L_qs(k) = .false.
            endif
            if (qg1d(k) .gt. R2) then
                rg(k) = qg1d(k)*rho(k)
                L_qg(k) = .true.
            else
                rg(k) = R1
                L_qg(k) = .false.
            endif
        enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope, and useful moments for snow.
!+---+-----------------------------------------------------------------+
        do k = kts, kte
            smo2(k) = 0.
            smob(k) = 0.
            smoc(k) = 0.
            smoz(k) = 0.
        enddo
        if (ANY(L_qs .eqv. .true.)) then
            do k = kts, kte
                if (.not. L_qs(k)) CYCLE
                tc0 = MIN(-0.1, temp(k)-273.15)
                smob(k) = rs(k)*oams

!..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
!.. then we must compute actual 2nd moment and use as reference.
                if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
                    smo2(k) = smob(k)
                else
                    loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
                    &         + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
                    &         + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
                    &         + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
                    &         + sa(10)*bm_s*bm_s*bm_s
                    a_ = 10.0**loga_
                    b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
                    &         + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
                    &         + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
                    &         + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
                    &         + sb(10)*bm_s*bm_s*bm_s
                    smo2(k) = (smob(k)/a_)**(1./b_)
                endif

!..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
                loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
                &         + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
                &         + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
                &         + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
                &         + sa(10)*cse(1)*cse(1)*cse(1)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
                &        + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
                &        + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
                &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
                smoc(k) = a_ * smo2(k)**b_

!..Calculate bm_s*2 (th) moment.  Useful for reflectivity.
                loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(3) &
                &         + sa(4)*tc0*cse(3) + sa(5)*tc0*tc0 &
                &         + sa(6)*cse(3)*cse(3) + sa(7)*tc0*tc0*cse(3) &
                &         + sa(8)*tc0*cse(3)*cse(3) + sa(9)*tc0*tc0*tc0 &
                &         + sa(10)*cse(3)*cse(3)*cse(3)
                a_ = 10.0**loga_
                b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(3) + sb(4)*tc0*cse(3) &
                &        + sb(5)*tc0*tc0 + sb(6)*cse(3)*cse(3) &
                &        + sb(7)*tc0*tc0*cse(3) + sb(8)*tc0*cse(3)*cse(3) &
                &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(3)*cse(3)*cse(3)
                smoz(k) = a_ * smo2(k)**b_
            enddo
        endif

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for graupel.
!+---+-----------------------------------------------------------------+

        if (ANY(L_qg .eqv. .true.)) then
            do k = kte, kts, -1
                ygra1 = alog10(max(1.e-9, rg(k)))
                zans1 = 3.0 + 2./7.*(ygra1+8.) + rand1
                N0_exp = 10.**(zans1)
                N0_exp = max(dble(gonv_min), min(N0_exp, dble(gonv_max)))
                lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
                lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
                ilamg(k) = 1./lamg
                N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
            enddo
            ! call graupel_psd_parameters(kts, kte, rand1, rg, ilamg, N0_g)
        endif

!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
        k_0 = kts
        if ( melti ) then
            K_LOOP:do k = kte-1, kts, -1
                if ((temp(k).gt.273.15) .and. L_qr(k)                         &
                &                            .and. (L_qs(k+1).or.L_qg(k+1)) ) then
                    k_0 = MAX(k+1, k_0)
                    EXIT K_LOOP
                endif
            enddo K_LOOP
        endif
!+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+

        do k = kts, kte
            ze_rain(k) = 1.e-22
            ze_snow(k) = 1.e-22
            ze_graupel(k) = 1.e-22
            if (L_qr(k)) ze_rain(k) = N0_r(k)*crg(4)*ilamr(k)**cre(4)
            if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
            &                           * (am_s/900.0)*(am_s/900.0)*smoz(k)
            if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
            &                              * (am_g/900.0)*(am_g/900.0)         &
            &                              * N0_g(k)*cgg(4)*ilamg(k)**cge(4)
        enddo

!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

        if (.not. iiwarm .and. melti .and. k_0.ge.2) then
            do k = k_0-1, kts, -1

!..Reflectivity contributed by melting snow
                if (allow_wet_snow .and. L_qs(k) .and. L_qs(k_0) ) then
                    SR = MAX(0.01, MIN(1.0 - rs(k)/(rs(k) + rr(k)), 0.99))
                    fmelt_s = DBLE(SR*SR)
                    eta = 0.d0
                    oM3 = 1./smoc(k)
                    M0 = (smob(k)*oM3)
                    Mrat = smob(k)*M0*M0*M0
                    slam1 = M0 * Lam0
                    slam2 = M0 * Lam1
                    do n = 1, nrbins
                        x = am_s * xxDs(n)**bm_s
                        call rayleigh_soak_wetgraupel (x, DBLE(ocms), DBLE(obms), &
                        &              fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                        &              CBACK, mixingrulestring_s, matrixstring_s,          &
                        &              inclusionstring_s, hoststring_s,                    &
                        &              hostmatrixstring_s, hostinclusionstring_s)
                        f_d = Mrat*(Kap0*DEXP(-slam1*xxDs(n))                     &
                        &              + Kap1*(M0*xxDs(n))**mu_s * DEXP(-slam2*xxDs(n)))
                        eta = eta + f_d * CBACK * simpson(n) * xdts(n)
                    enddo
                    ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
                endif

!..Reflectivity contributed by melting graupel
                if (allow_wet_graupel .and. L_qg(k) .and. L_qg(k_0) ) then
                    SR = MAX(0.01, MIN(1.0 - rg(k)/(rg(k) + rr(k)), 0.99))
                    fmelt_g = DBLE(SR*SR)
                    eta = 0.d0
                    lamg = 1./ilamg(k)
                    do n = 1, nrbins
                        x = am_g * xxDg(n)**bm_g
                        call rayleigh_soak_wetgraupel (x, DBLE(ocmg), DBLE(obmg), &
                        &              fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                        &              CBACK, mixingrulestring_g, matrixstring_g,          &
                        &              inclusionstring_g, hoststring_g,                    &
                        &              hostmatrixstring_g, hostinclusionstring_g)
                        f_d = N0_g(k)*xxDg(n)**mu_g * DEXP(-lamg*xxDg(n))
                        eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
                    enddo
                    ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
                endif

            enddo
        endif

        do k = kte, kts, -1
            dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
        enddo

!..Reflectivity-weighted terminal velocity (snow, rain, graupel, mix).
        if (do_vt_dBZ) then
            do k = kte, kts, -1
                vt_dBZ(k) = 1.E-3
                if (rs(k).gt.R2) then
                    Mrat = smob(k) / smoc(k)
                    ils1 = 1./(Mrat*Lam0 + fv_s)
                    ils2 = 1./(Mrat*Lam1 + fv_s)
                    t1_vts = Kap0*csg(5)*ils1**cse(5)
                    t2_vts = Kap1*Mrat**mu_s*csg(11)*ils2**cse(11)
                    ils1 = 1./(Mrat*Lam0)
                    ils2 = 1./(Mrat*Lam1)
                    t3_vts = Kap0*csg(6)*ils1**cse(6)
                    t4_vts = Kap1*Mrat**mu_s*csg(12)*ils2**cse(12)
                    vts_dbz_wt = rhof(k)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)
                    if (temp(k).ge.273.15 .and. temp(k).lt.275.15) then
                        vts_dbz_wt = vts_dbz_wt*1.5
                    elseif (temp(k).ge.275.15) then
                        vts_dbz_wt = vts_dbz_wt*2.0
                    endif
                else
                    vts_dbz_wt = 1.E-3
                endif

                if (rr(k).gt.R1) then
                    lamr = 1./ilamr(k)
                    vtr_dbz_wt = rhof(k)*av_r*crg(13)*(lamr+fv_r)**(-cre(13))      &
                        / (crg(4)*lamr**(-cre(4)))
                else
                    vtr_dbz_wt = 1.E-3
                endif

                if (rg(k).gt.R2) then
                    lamg = 1./ilamg(k)
                    vtg_dbz_wt = rhof(k)*av_g*cgg(5)*lamg**(-cge(5))               &
                        / (cgg(4)*lamg**(-cge(4)))
                else
                    vtg_dbz_wt = 1.E-3
                endif

                vt_dBZ(k) = (vts_dbz_wt*ze_snow(k) + vtr_dbz_wt*ze_rain(k)      &
                    + vtg_dbz_wt*ze_graupel(k))                        &
                    / (ze_rain(k)+ze_snow(k)+ze_graupel(k))
            enddo
        endif

    end subroutine calc_refl10cm
!


end module module_mp_thompson_utils

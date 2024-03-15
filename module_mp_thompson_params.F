module module_mp_thompson_params

    use mpas_kind_types

    implicit none

    logical, parameter :: iiwarm = .false.
    logical, public :: is_aerosol_aware = .false.

    LOGICAL, PARAMETER :: dustyIce = .true.
    LOGICAL, PARAMETER :: homogIce = .true.
    LOGICAL, PARAMETER:: is_hail_aware = .false.

    INTEGER, PARAMETER :: IFDRY = 0
    REAL, PARAMETER :: T_0 = 273.15
    REAL, PARAMETER :: PI = 3.1415926536

!..Densities of rain, snow, graupel, and cloud ice.
    REAL, PARAMETER:: rho_w2 = 1000.0 !AAJ change to rho_w2 to solve MPAS same var name conflict
    REAL :: rho_s2 = 100.0 !AAJ change to rho_s2 to solve MPAS same var name conflict
    ! REAL, PARAMETER:: rho_g = 500.0
    REAL, PARAMETER:: rho_i = 890.0
    INTEGER, PARAMETER:: NRHG = 9
    INTEGER, PARAMETER:: NRHG1 = 1
    INTEGER :: dimNRHG

    REAL, DIMENSION(NRHG), PARAMETER:: &
        rho_g = (/50., 100., 200., 300., 400., 500., 600., 700., 800./)
    INTEGER, PARAMETER :: idx_bg1 = 6 ! index for rhog when mp=8 or 28

!..Prescribed number of cloud droplets.  Set according to known data or
!.. roughly 100 per cc (100.E6 m^-3) for Maritime cases and
!.. 300 per cc (300.E6 m^-3) for Continental.  Gamma shape parameter,
!.. mu_c, calculated based on Nt_c is important in autoconversion
!.. scheme.  In 2-moment cloud water, Nt_c represents a maximum of
!.. droplet concentration and nu_c is also variable depending on local
!.. droplet number concentration.
!     REAL, PARAMETER:: Nt_c = 100.E6
    REAL, PARAMETER:: Nt_c_max = 1999.E6
    REAL:: Nt_c

!..Declaration of constants for assumed CCN/IN aerosols when none in
!.. the input data.  Look inside the init routine for modifications
!.. due to surface land-sea points or vegetation characteristics.
    REAL, PARAMETER:: naIN0 = 1.5E6
    REAL, PARAMETER:: naIN1 = 0.5E6
    REAL, PARAMETER:: naCCN0 = 300.0E6
    REAL, PARAMETER:: naCCN1 = 50.0E6

!..Generalized gamma distributions for rain, graupel and cloud ice.
!.. N(D) = N_0 * D**mu * exp(-lamda*D);  mu=0 is exponential.
    REAL, PARAMETER:: mu_r = 0.0
    REAL, PARAMETER:: mu_g = 0.0
    REAL, PARAMETER:: mu_i = 0.0
    REAL:: mu_c

!..Sum of two gamma distrib for snow (Field et al. 2005).
!.. N(D) = M2**4/M3**3 * [Kap0*exp(-M2*Lam0*D/M3)
!..    + Kap1*(M2/M3)**mu_s * D**mu_s * exp(-M2*Lam1*D/M3)]
!.. M2 and M3 are the (bm_s)th and (bm_s+1)th moments respectively
!.. calculated as function of ice water content and temperature.
    REAL, PARAMETER:: mu_s = 0.6357
    REAL, PARAMETER:: Kap0 = 490.6
    REAL, PARAMETER:: Kap1 = 17.46
    REAL, PARAMETER:: Lam0 = 20.78
    REAL, PARAMETER:: Lam1 = 3.29

!..Y-intercept parameter for graupel is not constant and depends on
!.. mixing ratio.  Also, when mu_g is non-zero, these become equiv
!.. y-intercept for an exponential distrib and proper values are
!.. computed based on same mixing ratio and total number concentration.
    REAL, PARAMETER:: gonv_min = 1.E2
    REAL, PARAMETER:: gonv_max = 1.E6

!..Mass power law relations:  mass = am*D**bm
!.. Snow from Field et al. (2005), others assume spherical form.
    REAL, PARAMETER:: am_r = PI*rho_w2/6.0
    REAL, PARAMETER:: bm_r = 3.0
    REAL, PARAMETER:: am_s = 0.069
    REAL, PARAMETER:: bm_s = 2.0
!    REAL, PARAMETER:: am_g = PI*rho_g/6.0
    REAL, DIMENSION(NRHG), PARAMETER:: am_g = (/              &
    &            PI*rho_g(1)/6.0, PI*rho_g(2)/6.0, PI*rho_g(3)/6.0,    &
    &            PI*rho_g(4)/6.0, PI*rho_g(5)/6.0, PI*rho_g(6)/6.0,    &
    &            PI*rho_g(7)/6.0, PI*rho_g(8)/6.0, PI*rho_g(9)/6.0 /)
    REAL, PARAMETER:: bm_g = 3.0
    REAL, PARAMETER:: am_i = PI*rho_i/6.0
    REAL, PARAMETER:: bm_i = 3.0

!..Fallspeed power laws relations:  v = (av*D**bv)*exp(-fv*D)
!.. Rain from Ferrier (1994), ice, snow, and graupel from
!.. Thompson et al (2008). Coefficient fv is zero for graupel/ice.
    REAL, PARAMETER:: av_r = 4854.0
    REAL, PARAMETER:: bv_r = 1.0
    REAL, PARAMETER:: fv_r = 195.0
    REAL, PARAMETER:: av_s = 40.0
    REAL, PARAMETER:: bv_s = 0.55
    REAL, PARAMETER:: fv_s = 100.0
!    REAL, PARAMETER:: av_g = 442.0
!    REAL, PARAMETER:: bv_g = 0.89
    REAL, PARAMETER:: av_g_old = 442.0
    REAL, PARAMETER:: bv_g_old = 0.89
    REAL, DIMENSION(NRHG):: & ! Computed from A. Heymsfield: Best - Reynolds relation
    &    av_g = (/ 45.9173813, 67.0867386, 98.0158463, 122.353378,     &
    &              143.204224, 161.794724, 178.762115, 194.488785,     &
    &              209.225876/)
    REAL, DIMENSION(NRHG):: & ! Computed from A. Heymsfield: Best - Reynolds relation
    &    bv_g = (/0.640961647, 0.640961647, 0.640961647, 0.640961647,  &
    &             0.640961647, 0.640961647, 0.640961647, 0.640961647,  &
    &             0.640961647/)
    REAL, PARAMETER:: a_coeff = 0.47244157
    REAL, PARAMETER:: b_coeff = 0.54698726
    REAL, PARAMETER:: av_i = 1493.9
    REAL, PARAMETER:: bv_i = 1.0
    REAL, PARAMETER:: av_c = 0.316946E8
    REAL, PARAMETER:: bv_c = 2.0

!..Capacitance of sphere and plates/aggregates: D**3, D**2
    REAL, PARAMETER:: C_cube = 0.5
    REAL, PARAMETER:: C_sqrd = 0.15

!..Collection efficiencies.  Rain/snow/graupel collection of cloud
!.. droplets use variables (Ef_rw, Ef_sw, Ef_gw respectively) and
!.. get computed elsewhere because they are dependent on stokes
!.. number.
    REAL, PARAMETER:: Ef_si = 0.05
    REAL, PARAMETER:: Ef_rs = 0.95
    REAL, PARAMETER:: Ef_rg = 0.75
    REAL, PARAMETER:: Ef_ri = 0.95

!..Minimum microphys values
!.. R1 value, 1.E-12, cannot be set lower because of numerical
!.. problems with Paul Field's moments and should not be set larger
!.. because of truncation problems in snow/ice growth.
    REAL, PARAMETER:: R1 = 1.E-12
    REAL, PARAMETER:: R2 = 1.E-6
    REAL, PARAMETER:: eps = 1.E-15

!..Constants in Cooper curve relation for cloud ice number.
    REAL, PARAMETER:: TNO = 5.0
    REAL, PARAMETER:: ATO = 0.304

!..Rho_not used in fallspeed relations (rho_not/rho)**.5 adjustment.
    REAL, PARAMETER:: rho_not = 101325.0/(287.05*298.0)

!..Schmidt number
    REAL, PARAMETER:: Sc = 0.632
    REAL:: Sc3

!..Homogeneous freezing temperature
    REAL, PARAMETER:: HGFR = 235.16

!..Water vapor and air gas constants at constant pressure
    REAL, PARAMETER:: Rv = 461.5
    REAL, PARAMETER:: oRv = 1./Rv
    REAL, PARAMETER:: R = 287.04
    REAL, PARAMETER:: Cp2 = 1004.0 !AAJ change to Cp2 to solve MPAS same var name conflict
    REAL, PARAMETER:: R_uni = 8.314                           ! J (mol K)-1

    DOUBLE PRECISION, PARAMETER:: k_b = 1.38065E-23           ! Boltzmann constant [J/K]
    DOUBLE PRECISION, PARAMETER:: M_w = 18.01528E-3           ! molecular mass of water [kg/mol]
    DOUBLE PRECISION, PARAMETER:: M_a = 28.96E-3              ! molecular mass of air [kg/mol]
    DOUBLE PRECISION, PARAMETER:: N_avo = 6.022E23            ! Avogadro number [1/mol]
    DOUBLE PRECISION, PARAMETER:: ma_w = M_w / N_avo          ! mass of water molecule [kg]
    REAL, PARAMETER:: ar_volume = 4./3.*PI*(2.5e-6)**3        ! assume radius of 0.025 micrometer, 2.5e-6 cm

!..Enthalpy of sublimation, vaporization, and fusion at 0C.
    REAL, PARAMETER:: lsub = 2.834E6
    REAL, PARAMETER:: lvap0 = 2.5E6
    REAL, PARAMETER:: lfus = lsub - lvap0
    REAL, PARAMETER:: olfus = 1./lfus

!..Ice initiates with this mass (kg), corresponding diameter calc.
!..Min diameters and mass of cloud, rain, snow, and graupel (m, kg).
    REAL, PARAMETER:: xm0i = 1.E-12
    REAL, PARAMETER:: D0c = 1.E-6
    REAL, PARAMETER:: D0r = 50.E-6
    REAL, PARAMETER:: D0s = 300.E-6
    REAL, PARAMETER:: D0g = 350.E-6
    REAL:: D0i, xm0s, xm0g

!..Lookup table dimensions
    INTEGER, PARAMETER:: nbins = 100
    INTEGER, PARAMETER:: nbc = nbins
    INTEGER, PARAMETER:: nbi = nbins
    INTEGER, PARAMETER:: nbr = nbins
    INTEGER, PARAMETER:: nbs = nbins
    INTEGER, PARAMETER:: nbg = nbins
    INTEGER, PARAMETER:: ntb_c = 37
    INTEGER, PARAMETER:: ntb_i = 64
    INTEGER, PARAMETER:: ntb_r = 37
    INTEGER, PARAMETER:: ntb_s = 37
    INTEGER, PARAMETER:: ntb_g = 37
    INTEGER, PARAMETER:: ntb_g1 = 37
    INTEGER, PARAMETER:: ntb_r1 = 37
    INTEGER, PARAMETER:: ntb_i1 = 55
    INTEGER, PARAMETER:: ntb_t = 9
    INTEGER:: nic1, nic2, nii2, nii3, nir2, nir3, nis2, nig2, nig3
    INTEGER, PARAMETER:: ntb_arc = 7
    INTEGER, PARAMETER:: ntb_arw = 9
    INTEGER, PARAMETER:: ntb_art = 7
    INTEGER, PARAMETER:: ntb_arr = 5
    INTEGER, PARAMETER:: ntb_ark = 4
    INTEGER, PARAMETER:: ntb_IN = 55
    INTEGER:: niIN2

    DOUBLE PRECISION, DIMENSION(nbins+1):: xDx
    DOUBLE PRECISION, DIMENSION(nbc):: Dc, dtc
    DOUBLE PRECISION, DIMENSION(nbi):: Di, dti
    DOUBLE PRECISION, DIMENSION(nbr):: Dr, dtr
    DOUBLE PRECISION, DIMENSION(nbs):: Ds, dts
    DOUBLE PRECISION, DIMENSION(nbg):: Dg, dtg
    DOUBLE PRECISION, DIMENSION(nbc):: t_Nc

!..Lookup tables for cloud water content (kg/m**3).
    REAL, DIMENSION(ntb_c), PARAMETER:: &
        r_c = (/1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
        1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
        1.e-2/)

!..Lookup tables for cloud ice content (kg/m**3).
    REAL, DIMENSION(ntb_i), PARAMETER:: &
        r_i = (/1.e-10,2.e-10,3.e-10,4.e-10, &
        5.e-10,6.e-10,7.e-10,8.e-10,9.e-10, &
        1.e-9,2.e-9,3.e-9,4.e-9,5.e-9,6.e-9,7.e-9,8.e-9,9.e-9, &
        1.e-8,2.e-8,3.e-8,4.e-8,5.e-8,6.e-8,7.e-8,8.e-8,9.e-8, &
        1.e-7,2.e-7,3.e-7,4.e-7,5.e-7,6.e-7,7.e-7,8.e-7,9.e-7, &
        1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
        1.e-3/)

!..Lookup tables for rain content (kg/m**3).
    REAL, DIMENSION(ntb_r), PARAMETER:: &
        r_r = (/1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
        1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
        1.e-2/)

!..Lookup tables for graupel content (kg/m**3).
    REAL, DIMENSION(ntb_g), PARAMETER:: &
        r_g = (/1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
        1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
        1.e-2/)

!..Lookup tables for snow content (kg/m**3).
    REAL, DIMENSION(ntb_s), PARAMETER:: &
        r_s = (/1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
        1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
        1.e-2/)

!..Lookup tables for rain y-intercept parameter (/m**4).
    REAL, DIMENSION(ntb_r1), PARAMETER:: &
        N0r_exp = (/1.e6,2.e6,3.e6,4.e6,5.e6,6.e6,7.e6,8.e6,9.e6, &
        1.e7,2.e7,3.e7,4.e7,5.e7,6.e7,7.e7,8.e7,9.e7, &
        1.e8,2.e8,3.e8,4.e8,5.e8,6.e8,7.e8,8.e8,9.e8, &
        1.e9,2.e9,3.e9,4.e9,5.e9,6.e9,7.e9,8.e9,9.e9, &
        1.e10/)

!..Lookup tables for graupel y-intercept parameter (/m**4).
    REAL, DIMENSION(ntb_g1), PARAMETER:: &
        N0g_exp = (/1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2, &
        1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3, &
        1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4, &
        1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5, &
        1.e6/)

!..Lookup tables for ice number concentration (/m**3).
    REAL, DIMENSION(ntb_i1), PARAMETER:: &
        Nt_i = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0, &
        1.e1,2.e1,3.e1,4.e1,5.e1,6.e1,7.e1,8.e1,9.e1, &
        1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2, &
        1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3, &
        1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4, &
        1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5, &
        1.e6/)

!..Aerosol table parameter: Number of available aerosols, vertical
!.. velocity, temperature, aerosol mean radius, and hygroscopicity.
    REAL, DIMENSION(ntb_arc), PARAMETER:: &
        ta_Na = (/10.0, 31.6, 100.0, 316.0, 1000.0, 3160.0, 10000.0/)
    REAL, DIMENSION(ntb_arw), PARAMETER:: &
        ta_Ww = (/0.01, 0.0316, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0/)
    REAL, DIMENSION(ntb_art), PARAMETER:: &
        ta_Tk = (/243.15, 253.15, 263.15, 273.15, 283.15, 293.15, 303.15/)
    REAL, DIMENSION(ntb_arr), PARAMETER:: &
        ta_Ra = (/0.01, 0.02, 0.04, 0.08, 0.16/)
    REAL, DIMENSION(ntb_ark), PARAMETER:: &
        ta_Ka = (/0.2, 0.4, 0.6, 0.8/)

!..Lookup tables for IN concentration (/m**3) from 0.001 to 1000/Liter.
    REAL, DIMENSION(ntb_IN), PARAMETER:: &
        Nt_IN = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0, &
        1.e1,2.e1,3.e1,4.e1,5.e1,6.e1,7.e1,8.e1,9.e1, &
        1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2, &
        1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3, &
        1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4, &
        1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5, &
        1.e6/)

!..For snow moments conversions (from Field et al. 2005)
    REAL, DIMENSION(10), PARAMETER:: &
        sa = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
        0.31255,   0.000204,  0.003199, 0.0,      -0.015952/)
    REAL, DIMENSION(10), PARAMETER:: &
        sb = (/ 0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
        0.060366,  0.000079,  0.000594, 0.0,      -0.003577/)

!..Temperatures (5 C interval 0 to -40) used in lookup tables.
    REAL, DIMENSION(ntb_t), PARAMETER:: &
        Tc = (/-0.01, -5., -10., -15., -20., -25., -30., -35., -40./)

!..Lookup tables for various accretion/collection terms.
!.. ntb_x refers to the number of elements for rain, snow, graupel,
!.. and temperature array indices.  Variables beginning with t-p/c/m/n
!.. represent lookup tables.  Save compile-time memory by making
!.. allocatable (2009Jun12, J. Michalakes).
!     INTEGER, PARAMETER:: R8SIZE = 8
!     INTEGER, PARAMETER:: R4SIZE = 4
    INTEGER, PARAMETER:: R8SIZE = R8KIND !MPAS kind
    INTEGER, PARAMETER:: R4SIZE = R4KIND !MPAS kind
    REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:,:,:)::             &
        tcg_racg, tmr_racg, tcr_gacr,                           & ! tmg_gacr
        tnr_racg, tnr_gacr
    REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:,:)::             &
        tcs_racs1, tmr_racs1, tcs_racs2, tmr_racs2,             &
        tcr_sacr1, tms_sacr1, tcr_sacr2, tms_sacr2,             &
        tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2
    REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:,:)::             &
        tpi_qcfz, tni_qcfz
    REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:,:)::             &
        tpi_qrfz, tpg_qrfz, tni_qrfz, tnr_qrfz
    REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:)::                 &
        tps_iaus, tni_iaus, tpi_ide
    REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:):: t_Efrw
    REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:):: t_Efsw
    REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:):: tnr_rev
    REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:)::               &
        tpc_wev, tnc_wev
    REAL (KIND=R4SIZE), ALLOCATABLE, DIMENSION(:,:,:,:,:):: tnccn_act

!..Variables holding a bunch of exponents and gamma values (cloud water,
!.. cloud ice, rain, snow, then graupel).
    REAL, DIMENSION(5,15):: cce, ccg
    REAL, DIMENSION(15)::  ocg1, ocg2
    REAL, DIMENSION(7):: cie, cig
    REAL:: oig1, oig2, obmi
    REAL, DIMENSION(13):: cre, crg
    REAL:: ore1, org1, org2, org3, obmr
    REAL, DIMENSION(17):: cse, csg
    REAL:: oams, obms, ocms
    REAL, DIMENSION(12,NRHG):: cge, cgg
    REAL:: oge1, ogg1, ogg2, ogg3, obmg
    REAL, DIMENSION(NRHG):: oamg, ocmg

!..Declaration of precomputed constants in various rate eqns.
    REAL:: t1_qr_qc, t1_qr_qi, t2_qr_qi, t1_qg_qc, t1_qs_qc, t1_qs_qi
    REAL:: t1_qr_ev, t2_qr_ev
    REAL:: t1_qs_sd, t2_qs_sd, t1_qg_sd, t2_qg_sd
    REAL:: t1_qs_me, t2_qs_me, t1_qg_me, t2_qg_me

end module module_mp_thompson_params

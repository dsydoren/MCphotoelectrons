
!-------------------------------------------------
!
SUBROUTINE Initiate_en_collisions

  USE ParallelOperationValues
  USE photoelectrons
  USE rng_wrapper
  USE PhysicalConstants


  USE MPI

  IMPLICIT NONE

  integer i, itmp

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)

  INTEGER init_random_seed
  INTEGER j
  REAL(8) myran

! prepare collsion frequencies
  CALL Prepare_N2_elastic           !  1
  CALL Prepare_N2_ionN2             !  2
  CALL Prepare_N2_ionN              !  3
  CALL Prepare_N2_vibtot            !  4
  CALL Prepare_N2_triplet_A3Zup     !  5
  CALL Prepare_N2_triplet_B3Pg      !  6
  CALL Prepare_N2_triplet_C3Pg      !  7
  CALL Prepare_N2_triplet_W3Du      !  8
  CALL Prepare_N2_triplet_Bp3Zum    !  9
  CALL Prepare_N2_singlet_a1Pg      ! 10
  CALL Prepare_N2_singlet_b1Pu      ! 11
  CALL Prepare_N2_singlet_bp1Zup    ! 12
  CALL Prepare_N2_singlet_c4p1Zu    ! 13
  CALL Prepare_N2_singlet_w1Du      ! 14
  CALL Prepare_N2_singlet_c1Pu      ! 15
  CALL Prepare_N2_singlet_ap1Zum    ! 16
  CALL Prepare_N2_singlet_app1Zgp   ! 17
  CALL Prepare_N2_singlet_other1Pu  ! 18
  CALL Prepare_N2_hls_15p8eV        ! 19
  CALL Prepare_N2_hls_VUV           ! 20
  CALL Prepare_N2_hls_17p3eV        ! 21
  CALL Prepare_N2_hls_Rydberg       ! 22
  CALL Prepare_N2_hls_tripman       ! 23

  CALL Prepare_O2_elastic           ! 24
  CALL Prepare_O2_ionO2             ! 25
  CALL Prepare_O2_ionO              ! 26
  CALL Prepare_O2_vibtot            ! 27
  CALL Prepare_O2_Rydberg           ! 28
  CALL Prepare_O2_AApc              ! 29
  CALL Prepare_O2_a1Du              ! 30
  CALL Prepare_O2_b1Zgp             ! 31
  CALL Prepare_O2_longband          ! 32
  CALL Prepare_O2_secondband        ! 33
  CALL Prepare_O2_13Pg              ! 34
  CALL Prepare_O2_8p9eV             ! 35
  CALL Prepare_O2_B3Zu              ! 36

  CALL Prepare_O_elastic            ! 37
  CALL Prepare_O_iontot             ! 38
  CALL Prepare_O_Rydberg            ! 39
  CALL Prepare_O_1D                 ! 40
  CALL Prepare_O_1S                 ! 41
  CALL Prepare_O_5P                 ! 42
  CALL Prepare_O_5S                 ! 43
  CALL Prepare_O_3s3S               ! 44
  CALL Prepare_O_3d3D               ! 45
  CALL Prepare_O_3sp3D              ! 46
  CALL Prepare_O_3p3P               ! 47
  CALL Prepare_O_5d3D               ! 48
  CALL Prepare_O_4d3D               ! 49
  CALL Prepare_O_2p53P              ! 50
  CALL Prepare_O_3spp3P             ! 51
  CALL Prepare_O_4dp3P              ! 52

  CALL Prepare_He_elastic           ! 53
  CALL Prepare_He_ionHe             ! 54
  CALL Prepare_He_21P               ! 55
  CALL Prepare_He_31P               ! 56
  CALL Prepare_He_41P               ! 57
  CALL Prepare_He_21S               ! 58
  CALL Prepare_He_23S               ! 59
  CALL Prepare_He_23P               ! 60

  CALL Prepare_Max_frequency_of_en_collisions

  init_random_seed = 912345678

! initialize random numbers generators
  IF (Rank_of_process.EQ.0) THEN 
     PRINT '(/2x,"Process ",i4," : Seed for well_random_seed: ",i12)', Rank_of_process, init_random_seed
     CALL well_random_seed(init_random_seed)
     DO j = 1, 1000000
        myran = well_random_number()
     END DO

     DO i = 1, N_of_processes - 1
        itmp = INT(2000000000.0_8 * well_random_number())                       
        CALL MPI_SEND(itmp, 1, MPI_INTEGER, i, 101, MPI_COMM_WORLD, ierr)
     END DO
  ELSE
     CALL MPI_RECV(init_random_seed, 1, MPI_INTEGER, 0, 101, MPI_COMM_WORLD, stattus, ierr)

     PRINT '(/2x,"Process ",i4," : Seed for well_random_seed: ",i12)', Rank_of_process, init_random_seed

     CALL well_random_seed(init_random_seed)

     DO j = 1, 1000000
        myran = well_random_number()
     END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! set reaction thresholds for collisions between electrons and neutrals
! N2:                         threshold [eV] {from Table 1 of [Majeed and Strickland]}
!
  threshold_eV(1)  =  0.0_8  !  1 CS_N2_elastic_cm2               n/a
  threshold_eV(2)  = 15.6_8  !  2 CS_N2_ionN2_cm2                15.6
  threshold_eV(3)  = 29.5_8  !  3 CS_N2_ionN_cm2                 29.5 ??? my assumption based on the low-energy end of the cross-section curve 
  threshold_eV(4)  =  0.9_8  !  4 CS_N2_vibtot_cm2                0.9
  threshold_eV(5)  =  6.2_8  !  5 CS_N2_triplet_A3Zup_cm2         6.2
  threshold_eV(6)  =  7.4_8  !  6 CS_N2_triplet_B3Pg_cm2          7.4
  threshold_eV(7)  = 11.0_8  !  7 CS_N2_triplet_C3Pg_cm2         11.0
  threshold_eV(8)  =  7.5_8  !  8 CS_N2_triplet_W3Du_cm2          7.5
  threshold_eV(9)  =  8.0_8  !  9 CS_N2_triplet_Bp3Zum_cm2        8.0
  threshold_eV(10) =  9.1_8  ! 10 CS_N2_singlet_a1Pg_cm2          9.1
  threshold_eV(11) = 12.6_8  ! 11 CS_N2_singlet_b1Pu_cm2         12.6
  threshold_eV(12) = 14.2_8  ! 12 CS_N2_singlet_bp1Zup_cm2       14.2
  threshold_eV(13) = 12.9_8  ! 13 CS_N2_singlet_c4p1Zu_cm2       12.9
  threshold_eV(14) =  8.89_8 ! 14 CS_N2_singlet_w1Du_cm2          8.89
  threshold_eV(15) = 12.9_8  ! 15 CS_N2_singlet_c1Pu_cm2         12.9
  threshold_eV(16) =  8.4_8  ! 16 CS_N2_singlet_ap1Zum_cm2        8.4
  threshold_eV(17) = 12.3_8  ! 17 CS_N2_singlet_app1Zgp_cm2      12.3
  threshold_eV(18) = 12.6_8  ! 18 CS_N2_singlet_other1Pu_cm2     12.6
  threshold_eV(19) = 16.4_8  ! 19 CS_N2_hls_15p8eV_cm2           16.4
  threshold_eV(20) = 23.7_8  ! 20 CS_N2_hls_VUV_cm2              23.7
  threshold_eV(21) = 17.4_8  ! 21 CS_N2_hls_17p3eV_cm2           17.4
  threshold_eV(22) = 40.0_8  ! 22 CS_N2_hls_Rydberg_cm2          40.0
  threshold_eV(23) = 11.0_8  ! 23 CS_N2_hls_tripman_cm2          11.0
!
! O2:                         threshold [eV] {from Table 2 of [Majeed and Strickland]}
!
  threshold_eV(24) =  0.0_8  ! 24 CS_O2_elastic_cm2               n/a
  threshold_eV(25) = 12.08_8 ! 25 CS_O2_ionO2_cm2                12.1   ! changed to 12.08 to make sure that ionization by 1025.7A line occurs
  threshold_eV(26) = 22.5_8  ! 26 CS_O2_ionO_cm2                 22.5  ??? my assumption based on the low-energy end of the cross-section curve
  threshold_eV(27) =  0.5_8  ! 27 CS_O2_vibtot_cm2                0.5 = average of 0.3(v=1) 0.4(v=2) 0.6(v=3) 0.8(v=4) ???
  threshold_eV(28) = 16.0_8  ! 28 CS_O2_Rydberg_cm2              16.0
  threshold_eV(29) =  4.5_8  ! 29 CS_O2_AApc_cm2                  4.5
  threshold_eV(30) =  1.0_8  ! 30 CS_O2_a1Du_cm2                  1.0
  threshold_eV(31) =  1.6_8  ! 31 CS_O2_b1Zgp_cm2                 1.6
  threshold_eV(32) = 10.0_8  ! 32 CS_O2_longband_cm2             10.0
  threshold_eV(33) = 10.3_8  ! 33 CS_O2_secondband_cm2           10.3
  threshold_eV(34) =  7.6_8  ! 34 CS_O2_13Pg_cm2                  7.6
  threshold_eV(35) =  8.9_8  ! 35 CS_O2_8p9eV_cm2                 8.9
  threshold_eV(36) =  8.3_8  ! 36 CS_O2_B3Zu_cm2                  8.3
!
! O:                          threshold [eV] {from Table 3 of [Majeed and Strickland]}
!
  threshold_eV(37) =  0.0_8  ! 37 CS_O_elastic_cm2                n/a
  threshold_eV(38) = 13.6_8  ! 38 CS_O_iontot_cm2                13.6
  threshold_eV(39) = 14.0_8  ! 39 CS_O_Rydberg_cm2               14.0
  threshold_eV(40) =  2.0_8  ! 40 CS_O_1D_cm2                     2.0
  threshold_eV(41) =  4.2_8  ! 41 CS_O_1S_cm2                     4.2
  threshold_eV(42) = 10.7_8  ! 42 CS_O_5P_cm2                    10.7
  threshold_eV(43) =  9.3_8  ! 43 CS_O_5S_cm2                     9.3
  threshold_eV(44) =  9.5_8  ! 44 CS_O_3s3S_cm2                   9.5
  threshold_eV(45) = 12.1_8  ! 45 CS_O_3d3D_cm2                  12.1  ??? cr-sect starts at 12.0, we take 12.0 instead of 12.1 ???
  threshold_eV(46) = 12.5_8  ! 46 CS_O_3sp3D_cm2                 12.5
  threshold_eV(47) = 11.0_8  ! 47 CS_O_3p3P_cm2                  11.0
  threshold_eV(48) = 13.0_8  ! 48 CS_O_5d3D_cm2                  13.0
  threshold_eV(49) = 12.8_8  ! 49 CS_O_4d3D_cm2                  12.8
  threshold_eV(50) = 15.0_8  ! 50 CS_O_2p53P_cm2                 15.0
  threshold_eV(51) = 14.0_8  ! 51 CS_O_3spp3P_cm2                14.0
  threshold_eV(52) = 16.0_8  ! 52 CS_O_4dp3P_cm2                 16.0
!
! He
!
  threshold_eV(53) =  0.0_8   ! CS_He_elastic_cm2                  n/a
  threshold_eV(54) = 24.6_8   ! CS_He_ionHe_cm2                   24.6
  threshold_eV(55) = 21.220_8 ! CS_He_21P_cm2                     21.22  ! 
  threshold_eV(56) = 23.087_8 ! CS_He_31P_cm2                     23.087 !
  threshold_eV(57) = 23.724_8 ! CS_He_41P_cm2                     23.724 ! these thresholds are from Biagi database, www.lxcat.net
  threshold_eV(58) = 20.620_8 ! CS_He_21S_cm2                     20.62  !
  threshold_eV(59) = 19.820_8 ! CS_He_23S_cm2                     19.82  !
  threshold_eV(60) = 20.960_8 ! CS_He_23P_cm2                     20.96  !

! set reaction thresholds for photoionization 
! according to tables in the book mentioned below, thresholds for some photoionization processes may be 
! different from those for ionization by electrons, for example threshold for process N2+e -> N+ + N + 2e
! is 29.5 eV (threshold_eV(3)) while the lowest energy spectrum bin with nonzero cross section for 
! N2 + hnu -> N+ + N + e is 500-550A where 550 A gives 22.54 eV and 500 A gives 24.79 eV
! for consistency with the photoionization tables, 
! we use the longest wave in the highest energy spectrum bin with nonzero cross section to calculate the threshold 

  threshold_photoion_He_eV       = h_Planck_Js * c_ms * 1.0d10 / (e_Cl * 550.0_8)  ! 22.54
 
  threshold_photoion_O_4S_eV     = h_Planck_Js * c_ms * 1.0d10 / (e_Cl * 950.0_8)
  threshold_photoion_O_2D_eV     = h_Planck_Js * c_ms * 1.0d10 / (e_Cl * 750.0_8)
  threshold_photoion_O_2P_eV     = h_Planck_Js * c_ms * 1.0d10 / (e_Cl * 700.0_8)
  threshold_photoion_O_4P_eV     = h_Planck_Js * c_ms * 1.0d10 / (e_Cl * 450.0_8)
  threshold_photoion_O_2Pst_eV   = h_Planck_Js * c_ms * 1.0d10 / (e_Cl * 304.0_8)   ! because there is non-zero cross-section for line 9 303.78A

  threshold_photoion_N2_to_N2_eV = h_Planck_Js * c_ms * 1.0d10 / (e_Cl * 800.0_8)
  threshold_photoion_N2_to_N_eV  = h_Planck_Js * c_ms * 1.0d10 / (e_Cl * 550.0_8)

  threshold_photoion_O2_to_O2_eV = h_Planck_Js * c_ms * 1.0d10 / (e_Cl * 1050.0_8)
  threshold_photoion_O2_to_O_eV  = h_Planck_Js * c_ms * 1.0d10 / (e_Cl * 700.0_8)

END SUBROUTINE Initiate_en_collisions

!--------------------------------------
!
SUBROUTINE Prepare_Max_frequency_of_en_collisions

  use combined_collisions

  implicit none

  integer alloc_err
  integer n
  real(8) energy_eV

! functions
  REAL(8) CSV_N2_elastic_m3s               ! 1
  REAL(8) CSV_N2_ionN2_m3s                 ! 2
  REAL(8) CSV_N2_ionN_m3s                  ! 3
  REAL(8) CSV_N2_vibtot_m3s                ! 4
  REAL(8) CSV_N2_triplet_A3Zup_m3s         ! 5
  REAL(8) CSV_N2_triplet_B3Pg_m3s          ! 6
  REAL(8) CSV_N2_triplet_C3Pg_m3s          ! 7
  REAL(8) CSV_N2_triplet_W3Du_m3s          ! 8
  REAL(8) CSV_N2_triplet_Bp3Zum_m3s        ! 9
  REAL(8) CSV_N2_singlet_a1Pg_m3s          ! 10
  REAL(8) CSV_N2_singlet_b1Pu_m3s          ! 11
  REAL(8) CSV_N2_singlet_bp1Zup_m3s        ! 12
  REAL(8) CSV_N2_singlet_c4p1Zu_m3s        ! 13
  REAL(8) CSV_N2_singlet_w1Du_m3s          ! 14
  REAL(8) CSV_N2_singlet_c1Pu_m3s          ! 15
  REAL(8) CSV_N2_singlet_ap1Zum_m3s        ! 16
  REAL(8) CSV_N2_singlet_app1Zgp_m3s       ! 17
  REAL(8) CSV_N2_singlet_other1Pu_m3s      ! 18
  REAL(8) CSV_N2_hls_15p8eV_m3s            ! 19
  REAL(8) CSV_N2_hls_VUV_m3s               ! 20
  REAL(8) CSV_N2_hls_17p3eV_m3s            ! 21
  REAL(8) CSV_N2_hls_Rydberg_m3s           ! 22
  REAL(8) CSV_N2_hls_tripman_m3s           ! 23

  REAL(8) CSV_O2_elastic_m3s               ! 24
  REAL(8) CSV_O2_ionO2_m3s                 ! 25
  REAL(8) CSV_O2_ionO_m3s                  ! 26
  REAL(8) CSV_O2_vibtot_m3s                ! 27
  REAL(8) CSV_O2_Rydberg_m3s               ! 28
  REAL(8) CSV_O2_AApc_m3s                  ! 29
  REAL(8) CSV_O2_a1Du_m3s                  ! 30
  REAL(8) CSV_O2_b1Zgp_m3s                 ! 31
  REAL(8) CSV_O2_longband_m3s              ! 32
  REAL(8) CSV_O2_secondband_m3s            ! 33
  REAL(8) CSV_O2_13Pg_m3s                  ! 34
  REAL(8) CSV_O2_8p9eV_m3s                 ! 35
  REAL(8) CSV_O2_B3Zu_m3s                  ! 36

  REAL(8) CSV_O_elastic_m3s                ! 37
  REAL(8) CSV_O_iontot_m3s                 ! 38
  REAL(8) CSV_O_Rydberg_m3s                ! 39
  REAL(8) CSV_O_1D_m3s                     ! 40
  REAL(8) CSV_O_1S_m3s                     ! 41
  REAL(8) CSV_O_5P_m3s                     ! 42
  REAL(8) CSV_O_5S_m3s                     ! 43
  REAL(8) CSV_O_3s3S_m3s                   ! 44
  REAL(8) CSV_O_3d3D_m3s                   ! 45
  REAL(8) CSV_O_3sp3D_m3s                  ! 46
  REAL(8) CSV_O_3p3P_m3s                   ! 47
  REAL(8) CSV_O_5d3D_m3s                   ! 48
  REAL(8) CSV_O_4d3D_m3s                   ! 49
  REAL(8) CSV_O_2p53P_m3s                  ! 50
  REAL(8) CSV_O_3spp3P_m3s                 ! 51
  REAL(8) CSV_O_4dp3P_m3s                  ! 52
  
  REAL(8) CSV_He_elastic_m3s               ! 53
  REAL(8) CSV_He_ionHe_m3s                 ! 54
  REAL(8) CSV_He_21P_m3s                   ! 55
  REAL(8) CSV_He_31P_m3s                   ! 56
  REAL(8) CSV_He_41P_m3s                   ! 57
  REAL(8) CSV_He_21S_m3s                   ! 58
  REAL(8) CSV_He_23S_m3s                   ! 59
  REAL(8) CSV_He_23P_m3s                   ! 60

  N_lowen = int((max_lowen_w_eV+1.0d-6*lowen_dw_eV) / lowen_dw_eV)
  allocate(lowen_CSV_N2_total_m3s(0:N_lowen), stat = alloc_err)
  allocate(lowen_CSV_O2_total_m3s(0:N_lowen), stat = alloc_err)
  allocate(lowen_CSV_O_total_m3s(0:N_lowen), stat = alloc_err)
  allocate(lowen_CSV_He_total_m3s(0:N_lowen), stat = alloc_err)

  N_highen     = int(max_highen_w_eV+1.0d-6)
  min_N_highen = int(max_lowen_w_eV +1.0d-6)
  allocate(highen_CSV_N2_total_m3s(min_N_highen:N_highen), stat = alloc_err)
  allocate(highen_CSV_O2_total_m3s(min_N_highen:N_highen), stat = alloc_err)
  allocate(highen_CSV_O_total_m3s( min_N_highen:N_highen), stat = alloc_err)
  allocate(highen_CSV_He_total_m3s(min_N_highen:N_highen), stat = alloc_err)

  do n = 0, N_lowen
     energy_eV = n * lowen_dw_eV

     lowen_CSV_N2_total_m3s(n) = CSV_N2_elastic_m3s(energy_eV) + &               ! 1
                               & CSV_N2_ionN2_m3s(energy_eV) + &                 ! 2
                               & CSV_N2_ionN_m3s(energy_eV) + &                  ! 3
                               & CSV_N2_vibtot_m3s(energy_eV) + &                ! 4
                               & CSV_N2_triplet_A3Zup_m3s(energy_eV) + &         ! 5
                               & CSV_N2_triplet_B3Pg_m3s(energy_eV) + &          ! 6
                               & CSV_N2_triplet_C3Pg_m3s(energy_eV) + &          ! 7
                               & CSV_N2_triplet_W3Du_m3s(energy_eV) + &          ! 8
                               & CSV_N2_triplet_Bp3Zum_m3s(energy_eV) + &        ! 9
                               & CSV_N2_singlet_a1Pg_m3s(energy_eV) + &          ! 10
                               & CSV_N2_singlet_b1Pu_m3s(energy_eV) + &          ! 11
                               & CSV_N2_singlet_bp1Zup_m3s(energy_eV) + &        ! 12
                               & CSV_N2_singlet_c4p1Zu_m3s(energy_eV) + &        ! 13
                               & CSV_N2_singlet_w1Du_m3s(energy_eV) + &          ! 14
                               & CSV_N2_singlet_c1Pu_m3s(energy_eV) + &          ! 15
                               & CSV_N2_singlet_ap1Zum_m3s(energy_eV) + &        ! 16
                               & CSV_N2_singlet_app1Zgp_m3s(energy_eV) + &       ! 17
                               & CSV_N2_singlet_other1Pu_m3s(energy_eV) + &      ! 18
                               & CSV_N2_hls_15p8eV_m3s(energy_eV) + &            ! 19
                               & CSV_N2_hls_VUV_m3s(energy_eV) + &               ! 20
                               & CSV_N2_hls_17p3eV_m3s(energy_eV) + &            ! 21
                               & CSV_N2_hls_Rydberg_m3s(energy_eV) + &           ! 22
                               & CSV_N2_hls_tripman_m3s(energy_eV)               ! 23

     lowen_CSV_O2_total_m3s(n) = CSV_O2_elastic_m3s(energy_eV) + &      ! 24
                               & CSV_O2_ionO2_m3s(energy_eV) + &        ! 25
                               & CSV_O2_ionO_m3s(energy_eV) + &         ! 26
                               & CSV_O2_vibtot_m3s(energy_eV) + &       ! 27
                               & CSV_O2_Rydberg_m3s(energy_eV) + &      ! 28
                               & CSV_O2_AApc_m3s(energy_eV) + &         ! 29
                               & CSV_O2_a1Du_m3s(energy_eV) + &         ! 30
                               & CSV_O2_b1Zgp_m3s(energy_eV) + &        ! 31
                               & CSV_O2_longband_m3s(energy_eV) + &     ! 32
                               & CSV_O2_secondband_m3s(energy_eV) + &   ! 33
                               & CSV_O2_13Pg_m3s(energy_eV) + &         ! 34
                               & CSV_O2_8p9eV_m3s(energy_eV) + &        ! 35
                               & CSV_O2_B3Zu_m3s(energy_eV)             ! 36


     lowen_CSV_O_total_m3s(n) = CSV_O_elastic_m3s(energy_eV) + &    ! 37
                              & CSV_O_iontot_m3s(energy_eV) + &     ! 38
                              & CSV_O_Rydberg_m3s(energy_eV) + &    ! 39
                              & CSV_O_1D_m3s(energy_eV) + &         ! 40
                              & CSV_O_1S_m3s(energy_eV) + &         ! 41
                              & CSV_O_5P_m3s(energy_eV) + &         ! 42
                              & CSV_O_5S_m3s(energy_eV) + &         ! 43
                              & CSV_O_3s3S_m3s(energy_eV) + &       ! 44
                              & CSV_O_3d3D_m3s(energy_eV) + &       ! 45
                              & CSV_O_3sp3D_m3s(energy_eV) + &      ! 46
                              & CSV_O_3p3P_m3s(energy_eV) + &       ! 47
                              & CSV_O_5d3D_m3s(energy_eV) + &       ! 48
                              & CSV_O_4d3D_m3s(energy_eV) + &       ! 49
                              & CSV_O_2p53P_m3s(energy_eV) + &      ! 50
                              & CSV_O_3spp3P_m3s(energy_eV) + &     ! 51
                              & CSV_O_4dp3P_m3s(energy_eV)          ! 52
 
     lowen_CSV_He_total_m3s(n) = CSV_He_elastic_m3s(energy_eV) + &    ! 53
                               & CSV_He_ionHe_m3s(energy_eV) + &      ! 54
                               & CSV_He_21P_m3s(energy_eV) + &        ! 55
                               & CSV_He_31P_m3s(energy_eV) + &        ! 56
                               & CSV_He_41P_m3s(energy_eV) + &        ! 57
                               & CSV_He_21S_m3s(energy_eV) + &        ! 58
                               & CSV_He_23S_m3s(energy_eV) + &        ! 59
                               & CSV_He_23P_m3s(energy_eV)            ! 60
 end do

  do n = min_N_highen, N_highen
     energy_eV = dble(n)

     highen_CSV_N2_total_m3s(n) = CSV_N2_elastic_m3s(energy_eV) + &               ! 1
                                & CSV_N2_ionN2_m3s(energy_eV) + &                 ! 2
                                & CSV_N2_ionN_m3s(energy_eV) + &                  ! 3
                                & CSV_N2_vibtot_m3s(energy_eV) + &                ! 4
                                & CSV_N2_triplet_A3Zup_m3s(energy_eV) + &         ! 5
                                & CSV_N2_triplet_B3Pg_m3s(energy_eV) + &          ! 6
                                & CSV_N2_triplet_C3Pg_m3s(energy_eV) + &          ! 7
                                & CSV_N2_triplet_W3Du_m3s(energy_eV) + &          ! 8
                                & CSV_N2_triplet_Bp3Zum_m3s(energy_eV) + &        ! 9
                                & CSV_N2_singlet_a1Pg_m3s(energy_eV) + &          ! 10
                                & CSV_N2_singlet_b1Pu_m3s(energy_eV) + &          ! 11
                                & CSV_N2_singlet_bp1Zup_m3s(energy_eV) + &        ! 12
                                & CSV_N2_singlet_c4p1Zu_m3s(energy_eV) + &        ! 13
                                & CSV_N2_singlet_w1Du_m3s(energy_eV) + &          ! 14
                                & CSV_N2_singlet_c1Pu_m3s(energy_eV) + &          ! 15
                                & CSV_N2_singlet_ap1Zum_m3s(energy_eV) + &        ! 16
                                & CSV_N2_singlet_app1Zgp_m3s(energy_eV) + &       ! 17
                                & CSV_N2_singlet_other1Pu_m3s(energy_eV) + &      ! 18
                                & CSV_N2_hls_15p8eV_m3s(energy_eV) + &            ! 19
                                & CSV_N2_hls_VUV_m3s(energy_eV) + &               ! 20
                                & CSV_N2_hls_17p3eV_m3s(energy_eV) + &            ! 21
                                & CSV_N2_hls_Rydberg_m3s(energy_eV) + &           ! 22
                                & CSV_N2_hls_tripman_m3s(energy_eV)               ! 23

     highen_CSV_O2_total_m3s(n) = CSV_O2_elastic_m3s(energy_eV) + &      ! 24
                                & CSV_O2_ionO2_m3s(energy_eV) + &        ! 25
                                & CSV_O2_ionO_m3s(energy_eV) + &         ! 26
                                & CSV_O2_vibtot_m3s(energy_eV) + &       ! 27
                                & CSV_O2_Rydberg_m3s(energy_eV) + &      ! 28
                                & CSV_O2_AApc_m3s(energy_eV) + &         ! 29
                                & CSV_O2_a1Du_m3s(energy_eV) + &         ! 30
                                & CSV_O2_b1Zgp_m3s(energy_eV) + &        ! 31
                                & CSV_O2_longband_m3s(energy_eV) + &     ! 32
                                & CSV_O2_secondband_m3s(energy_eV) + &   ! 33
                                & CSV_O2_13Pg_m3s(energy_eV) + &         ! 34
                                & CSV_O2_8p9eV_m3s(energy_eV) + &        ! 35
                                & CSV_O2_B3Zu_m3s(energy_eV)             ! 36


     highen_CSV_O_total_m3s(n) = CSV_O_elastic_m3s(energy_eV) + &    ! 37
                               & CSV_O_iontot_m3s(energy_eV) + &     ! 38
                               & CSV_O_Rydberg_m3s(energy_eV) + &    ! 39
                               & CSV_O_1D_m3s(energy_eV) + &         ! 40
                               & CSV_O_1S_m3s(energy_eV) + &         ! 41
                               & CSV_O_5P_m3s(energy_eV) + &         ! 42
                               & CSV_O_5S_m3s(energy_eV) + &         ! 43
                               & CSV_O_3s3S_m3s(energy_eV) + &       ! 44
                               & CSV_O_3d3D_m3s(energy_eV) + &       ! 45
                               & CSV_O_3sp3D_m3s(energy_eV) + &      ! 46
                               & CSV_O_3p3P_m3s(energy_eV) + &       ! 47
                               & CSV_O_5d3D_m3s(energy_eV) + &       ! 48
                               & CSV_O_4d3D_m3s(energy_eV) + &       ! 49
                               & CSV_O_2p53P_m3s(energy_eV) + &      ! 50
                               & CSV_O_3spp3P_m3s(energy_eV) + &     ! 51
                               & CSV_O_4dp3P_m3s(energy_eV)          ! 52

     highen_CSV_He_total_m3s(n) = CSV_He_elastic_m3s(energy_eV) + &    ! 53
                                & CSV_He_ionHe_m3s(energy_eV) + &      ! 54
                                & CSV_He_21P_m3s(energy_eV) + &        ! 55
                                & CSV_He_31P_m3s(energy_eV) + &        ! 56
                                & CSV_He_41P_m3s(energy_eV) + &        ! 57
                                & CSV_He_21S_m3s(energy_eV) + &        ! 58
                                & CSV_He_23S_m3s(energy_eV) + &        ! 59
                                & CSV_He_23P_m3s(energy_eV)            ! 60
  end do

END SUBROUTINE Prepare_Max_frequency_of_en_collisions

!--------------------------------------
!
REAL(8) FUNCTION fl_Max_frequency_of_en_collisions_s1(energy_eV, Nn_He_m3, Nn_O_m3, Nn_N2_m3, Nn_O2_m3)

  USE combined_collisions

  IMPLICIT NONE

  REAL(8) energy_eV
  REAL(8) Nn_He_m3, Nn_O_m3, Nn_N2_m3, Nn_O2_m3

  REAL(8) temp, alfa
  INTEGER n

  REAL(8) my_CSV_N2_total_m3s
  REAL(8) my_CSV_O2_total_m3s
  REAL(8) my_CSV_O_total_m3s
  REAL(8) my_CSV_He_total_m3s

  IF (energy_eV.LE.max_lowen_w_eV) THEN
     temp = energy_eV / lowen_dw_eV
     n = MIN(INT(temp), N_lowen-1)
     alfa = MAX(0.0_8, MIN(temp - n, 1.0_8))
     my_CSV_N2_total_m3s = lowen_CSV_N2_total_m3s(n) * (1.0_8 - alfa) + lowen_CSV_N2_total_m3s(n+1) * alfa
     my_CSV_O2_total_m3s = lowen_CSV_O2_total_m3s(n) * (1.0_8 - alfa) + lowen_CSV_O2_total_m3s(n+1) * alfa
     my_CSV_O_total_m3s  = lowen_CSV_O_total_m3s(n)  * (1.0_8 - alfa) + lowen_CSV_O_total_m3s(n+1)  * alfa
     my_CSV_He_total_m3s = lowen_CSV_He_total_m3s(n) * (1.0_8 - alfa) + lowen_CSV_He_total_m3s(n+1)  * alfa
  ELSE
     n = MAX(min_N_highen, INT(energy_eV))
     alfa = MAX(0.0_8, MIN(energy_eV - n, 1.0_8))
     my_CSV_N2_total_m3s = highen_CSV_N2_total_m3s(n) * (1.0_8 - alfa) + highen_CSV_N2_total_m3s(n+1) * alfa
     my_CSV_O2_total_m3s = highen_CSV_O2_total_m3s(n) * (1.0_8 - alfa) + highen_CSV_O2_total_m3s(n+1) * alfa
     my_CSV_O_total_m3s  = highen_CSV_O_total_m3s(n)  * (1.0_8 - alfa) + highen_CSV_O_total_m3s(n+1)  * alfa
     my_CSV_He_total_m3s = highen_CSV_O_total_m3s(n)  * (1.0_8 - alfa) + highen_CSV_O_total_m3s(n+1)  * alfa
  END IF

  fl_Max_frequency_of_en_collisions_s1 = my_CSV_N2_total_m3s * Nn_N2_m3 + &
                                       & my_CSV_O2_total_m3s * Nn_O2_m3 + &
                                       & my_CSV_O_total_m3s  * Nn_O_m3 + &
                                       & my_CSV_He_total_m3s * Nn_He_m3

if (fl_Max_frequency_of_en_collisions_s1.lt.0.0_8) then
   print '("## neg-coll-error ## :: ",9(2x,e12.5))', &
        & energy_eV, &
        & my_CSV_N2_total_m3s, my_CSV_O2_total_m3s, my_CSV_O_total_m3s, my_CSV_He_total_m3s, &
        & Nn_N2_m3, Nn_O2_m3, Nn_O_m3, Nn_He_m3
   stop
end if

END FUNCTION fl_Max_frequency_of_en_collisions_s1

!--------------------------------------
!
REAL(8) FUNCTION fl_Frequency_of_Collision_s1(kind, energy_eV, Nn_He_m3, Nn_O_m3, Nn_N2_m3, Nn_O2_m3)

  USE Photoelectrons

  implicit none

  INTEGER kind
  REAL(8) energy_eV
  REAL(8) Nn_He_m3, Nn_O_m3, Nn_N2_m3, Nn_O2_m3

! functions
  REAL(8) CSV_N2_elastic_m3s               ! 1
  REAL(8) CSV_N2_ionN2_m3s                 ! 2
  REAL(8) CSV_N2_ionN_m3s                  ! 3
  REAL(8) CSV_N2_vibtot_m3s                ! 4
  REAL(8) CSV_N2_triplet_A3Zup_m3s         ! 5
  REAL(8) CSV_N2_triplet_B3Pg_m3s          ! 6
  REAL(8) CSV_N2_triplet_C3Pg_m3s          ! 7
  REAL(8) CSV_N2_triplet_W3Du_m3s          ! 8
  REAL(8) CSV_N2_triplet_Bp3Zum_m3s        ! 9
  REAL(8) CSV_N2_singlet_a1Pg_m3s          ! 10
  REAL(8) CSV_N2_singlet_b1Pu_m3s          ! 11
  REAL(8) CSV_N2_singlet_bp1Zup_m3s        ! 12
  REAL(8) CSV_N2_singlet_c4p1Zu_m3s        ! 13
  REAL(8) CSV_N2_singlet_w1Du_m3s          ! 14
  REAL(8) CSV_N2_singlet_c1Pu_m3s          ! 15
  REAL(8) CSV_N2_singlet_ap1Zum_m3s        ! 16
  REAL(8) CSV_N2_singlet_app1Zgp_m3s       ! 17
  REAL(8) CSV_N2_singlet_other1Pu_m3s      ! 18
  REAL(8) CSV_N2_hls_15p8eV_m3s            ! 19
  REAL(8) CSV_N2_hls_VUV_m3s               ! 20
  REAL(8) CSV_N2_hls_17p3eV_m3s            ! 21
  REAL(8) CSV_N2_hls_Rydberg_m3s           ! 22
  REAL(8) CSV_N2_hls_tripman_m3s           ! 23

  REAL(8) CSV_O2_elastic_m3s               ! 24
  REAL(8) CSV_O2_ionO2_m3s                 ! 25
  REAL(8) CSV_O2_ionO_m3s                  ! 26
  REAL(8) CSV_O2_vibtot_m3s                ! 27
  REAL(8) CSV_O2_Rydberg_m3s               ! 28
  REAL(8) CSV_O2_AApc_m3s                  ! 29
  REAL(8) CSV_O2_a1Du_m3s                  ! 30
  REAL(8) CSV_O2_b1Zgp_m3s                 ! 31
  REAL(8) CSV_O2_longband_m3s              ! 32
  REAL(8) CSV_O2_secondband_m3s            ! 33
  REAL(8) CSV_O2_13Pg_m3s                  ! 34
  REAL(8) CSV_O2_8p9eV_m3s                 ! 35
  REAL(8) CSV_O2_B3Zu_m3s                  ! 36

  REAL(8) CSV_O_elastic_m3s                ! 37
  REAL(8) CSV_O_iontot_m3s                 ! 38
  REAL(8) CSV_O_Rydberg_m3s                ! 39
  REAL(8) CSV_O_1D_m3s                     ! 40
  REAL(8) CSV_O_1S_m3s                     ! 41
  REAL(8) CSV_O_5P_m3s                     ! 42
  REAL(8) CSV_O_5S_m3s                     ! 43
  REAL(8) CSV_O_3s3S_m3s                   ! 44
  REAL(8) CSV_O_3d3D_m3s                   ! 45
  REAL(8) CSV_O_3sp3D_m3s                  ! 46
  REAL(8) CSV_O_3p3P_m3s                   ! 47
  REAL(8) CSV_O_5d3D_m3s                   ! 48
  REAL(8) CSV_O_4d3D_m3s                   ! 49
  REAL(8) CSV_O_2p53P_m3s                  ! 50
  REAL(8) CSV_O_3spp3P_m3s                 ! 51
  REAL(8) CSV_O_4dp3P_m3s                  ! 52

  REAL(8) CSV_He_elastic_m3s               ! 53
  REAL(8) CSV_He_ionHe_m3s                 ! 54
  REAL(8) CSV_He_21P_m3s                   ! 55
  REAL(8) CSV_He_31P_m3s                   ! 56
  REAL(8) CSV_He_41P_m3s                   ! 57
  REAL(8) CSV_He_21S_m3s                   ! 58
  REAL(8) CSV_He_23S_m3s                   ! 59
  REAL(8) CSV_He_23P_m3s                   ! 60

! fool proof check
  if ((kind.lt.1).or.(kind.gt.N_colkind)) then
     print '("Error-1 in fl_Frequency_of_Collision_s1, illegal kind of collision :: kind/min/max = ",3(2x,i3))', kind, 1, N_colkind
     stop
  end if

  fl_Frequency_of_Collision_s1 = 0.0_8

  SELECT CASE (kind)

  CASE (1)
     fl_Frequency_of_Collision_s1 = CSV_N2_elastic_m3s(energy_eV) * Nn_N2_m3
  CASE (2) 
     fl_Frequency_of_Collision_s1 = CSV_N2_ionN2_m3s(energy_eV) * Nn_N2_m3
  CASE (3)
     fl_Frequency_of_Collision_s1 = CSV_N2_ionN_m3s(energy_eV) * Nn_N2_m3
  CASE (4)
     fl_Frequency_of_Collision_s1 = CSV_N2_vibtot_m3s(energy_eV) * Nn_N2_m3
  CASE (5)
     fl_Frequency_of_Collision_s1 = CSV_N2_triplet_A3Zup_m3s(energy_eV) * Nn_N2_m3
  CASE (6)
     fl_Frequency_of_Collision_s1 = CSV_N2_triplet_B3Pg_m3s(energy_eV) * Nn_N2_m3
  CASE (7)
     fl_Frequency_of_Collision_s1 = CSV_N2_triplet_C3Pg_m3s(energy_eV) * Nn_N2_m3
  CASE (8)
     fl_Frequency_of_Collision_s1 = CSV_N2_triplet_W3Du_m3s(energy_eV) * Nn_N2_m3
  CASE (9)
     fl_Frequency_of_Collision_s1 = CSV_N2_triplet_Bp3Zum_m3s(energy_eV) * Nn_N2_m3
  CASE (10)
     fl_Frequency_of_Collision_s1 = CSV_N2_singlet_a1Pg_m3s(energy_eV) * Nn_N2_m3
  CASE (11)
     fl_Frequency_of_Collision_s1 = CSV_N2_singlet_b1Pu_m3s(energy_eV) * Nn_N2_m3
  CASE (12)
     fl_Frequency_of_Collision_s1 = CSV_N2_singlet_bp1Zup_m3s(energy_eV) * Nn_N2_m3
  CASE (13)
     fl_Frequency_of_Collision_s1 = CSV_N2_singlet_c4p1Zu_m3s(energy_eV) * Nn_N2_m3
  CASE (14)
     fl_Frequency_of_Collision_s1 = CSV_N2_singlet_w1Du_m3s(energy_eV) * Nn_N2_m3
  CASE (15)
     fl_Frequency_of_Collision_s1 = CSV_N2_singlet_c1Pu_m3s(energy_eV) * Nn_N2_m3
  CASE (16)
     fl_Frequency_of_Collision_s1 = CSV_N2_singlet_ap1Zum_m3s(energy_eV) * Nn_N2_m3
  CASE (17)
     fl_Frequency_of_Collision_s1 = CSV_N2_singlet_app1Zgp_m3s(energy_eV) * Nn_N2_m3
  CASE (18)
     fl_Frequency_of_Collision_s1 = CSV_N2_singlet_other1Pu_m3s(energy_eV) * Nn_N2_m3
  CASE (19)
     fl_Frequency_of_Collision_s1 = CSV_N2_hls_15p8eV_m3s(energy_eV) * Nn_N2_m3
  CASE (20)
     fl_Frequency_of_Collision_s1 = CSV_N2_hls_VUV_m3s(energy_eV) * Nn_N2_m3
  CASE (21)
     fl_Frequency_of_Collision_s1 = CSV_N2_hls_17p3eV_m3s(energy_eV) * Nn_N2_m3
  CASE (22)
     fl_Frequency_of_Collision_s1 = CSV_N2_hls_Rydberg_m3s(energy_eV) * Nn_N2_m3
  CASE (23)
     fl_Frequency_of_Collision_s1 = CSV_N2_hls_tripman_m3s(energy_eV) * Nn_N2_m3

  CASE (24)
     fl_Frequency_of_Collision_s1 = CSV_O2_elastic_m3s(energy_eV) * Nn_O2_m3
  CASE (25)
     fl_Frequency_of_Collision_s1 = CSV_O2_ionO2_m3s(energy_eV) * Nn_O2_m3
  CASE (26)
     fl_Frequency_of_Collision_s1 = CSV_O2_ionO_m3s(energy_eV) * Nn_O2_m3
  CASE (27)
     fl_Frequency_of_Collision_s1 = CSV_O2_vibtot_m3s(energy_eV) * Nn_O2_m3
  CASE (28)
     fl_Frequency_of_Collision_s1 = CSV_O2_Rydberg_m3s(energy_eV) * Nn_O2_m3
  CASE (29)
     fl_Frequency_of_Collision_s1 = CSV_O2_AApc_m3s(energy_eV) * Nn_O2_m3
  CASE (30)
     fl_Frequency_of_Collision_s1 = CSV_O2_a1Du_m3s(energy_eV) * Nn_O2_m3
  CASE (31)
     fl_Frequency_of_Collision_s1 = CSV_O2_b1Zgp_m3s(energy_eV) * Nn_O2_m3
  CASE (32)
     fl_Frequency_of_Collision_s1 = CSV_O2_longband_m3s(energy_eV) * Nn_O2_m3
  CASE (33)
     fl_Frequency_of_Collision_s1 = CSV_O2_secondband_m3s(energy_eV) * Nn_O2_m3
  CASE (34)
     fl_Frequency_of_Collision_s1 = CSV_O2_13Pg_m3s(energy_eV) * Nn_O2_m3
  CASE (35)
     fl_Frequency_of_Collision_s1 = CSV_O2_8p9eV_m3s(energy_eV) * Nn_O2_m3
  CASE (36)
     fl_Frequency_of_Collision_s1 = CSV_O2_B3Zu_m3s(energy_eV) * Nn_O2_m3

  CASE (37)
     fl_Frequency_of_Collision_s1 = CSV_O_elastic_m3s(energy_eV) * Nn_O_m3
  CASE (38)
     fl_Frequency_of_Collision_s1 = CSV_O_iontot_m3s(energy_eV) * Nn_O_m3
  CASE (39)
     fl_Frequency_of_Collision_s1 = CSV_O_Rydberg_m3s(energy_eV) * Nn_O_m3
  CASE (40)
     fl_Frequency_of_Collision_s1 = CSV_O_1D_m3s(energy_eV) * Nn_O_m3
  CASE (41)
     fl_Frequency_of_Collision_s1 = CSV_O_1S_m3s(energy_eV) * Nn_O_m3
  CASE (42)
     fl_Frequency_of_Collision_s1 = CSV_O_5P_m3s(energy_eV) * Nn_O_m3
  CASE (43)
     fl_Frequency_of_Collision_s1 = CSV_O_5S_m3s(energy_eV) * Nn_O_m3
  CASE (44)
     fl_Frequency_of_Collision_s1 = CSV_O_3s3S_m3s(energy_eV) * Nn_O_m3
  CASE (45)
     fl_Frequency_of_Collision_s1 = CSV_O_3d3D_m3s(energy_eV) * Nn_O_m3
  CASE (46)
     fl_Frequency_of_Collision_s1 = CSV_O_3sp3D_m3s(energy_eV) * Nn_O_m3
  CASE (47)
     fl_Frequency_of_Collision_s1 = CSV_O_3p3P_m3s(energy_eV) * Nn_O_m3
  CASE (48)
     fl_Frequency_of_Collision_s1 = CSV_O_5d3D_m3s(energy_eV) * Nn_O_m3
  CASE (49)
     fl_Frequency_of_Collision_s1 = CSV_O_4d3D_m3s(energy_eV) * Nn_O_m3
  CASE (50)
     fl_Frequency_of_Collision_s1 = CSV_O_2p53P_m3s(energy_eV) * Nn_O_m3
  CASE (51)
     fl_Frequency_of_Collision_s1 = CSV_O_3spp3P_m3s(energy_eV) * Nn_O_m3
  CASE (52)
     fl_Frequency_of_Collision_s1 = CSV_O_4dp3P_m3s(energy_eV) * Nn_O_m3

  CASE (53)
     fl_Frequency_of_Collision_s1 = CSV_He_elastic_m3s(energy_eV) * Nn_He_m3
  CASE (54)
     fl_Frequency_of_Collision_s1 = CSV_He_ionHe_m3s(energy_eV) * Nn_He_m3
  CASE (55)
     fl_Frequency_of_Collision_s1 = CSV_He_21P_m3s(energy_eV) * Nn_He_m3
  CASE (56)
     fl_Frequency_of_Collision_s1 = CSV_He_31P_m3s(energy_eV) * Nn_He_m3
  CASE (57)
     fl_Frequency_of_Collision_s1 = CSV_He_41P_m3s(energy_eV) * Nn_He_m3
  CASE (58)
     fl_Frequency_of_Collision_s1 = CSV_He_21S_m3s(energy_eV) * Nn_He_m3
  CASE (59)
     fl_Frequency_of_Collision_s1 = CSV_He_23S_m3s(energy_eV) * Nn_He_m3
  CASE (60)
     fl_Frequency_of_Collision_s1 = CSV_He_23P_m3s(energy_eV) * Nn_He_m3

  END SELECT

END FUNCTION fl_Frequency_of_Collision_s1

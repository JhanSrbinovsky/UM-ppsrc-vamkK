! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for McClatchy atymospheric profiles.

MODULE MCC_data

  USE scm_utils, ONLY:                                                        &
      rmdi, imdi, z_th, z_rh, nmod_lv                                         &
    , zhook_in, zhook_out, jprb, lhook, dr_hook  

  USE s_interp_mod, ONLY: interp1d

  USE atmos_constants_mod,  ONLY: cp, r, kappa, p_zero, recip_kappa

  IMPLICIT NONE

!=============================================================================
!
! Description:
! Contains Standard Mcclatchey Atmospheric profiles to use in place of
! missing initial profile data when the SCM top of model height is higher
! the inital profile data provided by the forcings.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!=============================================================================
!
! Declarations:

  ! Mcclatchy Tropical Profile
  !===========================
  REAL :: MCC_trp_z(63)
  REAL :: MCC_trp_p(63)
  REAL :: MCC_trp_t(63)
  REAL :: MCC_trp_q(63)
  REAL :: MCC_trp_th(63)
  REAL :: MCC_trp_o3(63)
  REAL :: MCC_trp_rho_air(63)
  REAL :: MCC_trp_rho_h2o(63)
  REAL :: MCC_trp_rho_o3(63)
  REAL :: MCC_trp_exner(63)

  ! Mcclatchy Mid-latitude winter profile
  !======================================
  REAL :: MCC_mlw_z(33)
  REAL :: MCC_mlw_p(33)
  REAL :: MCC_mlw_t(33)
  REAL :: MCC_mlw_q(33)
  REAL :: MCC_mlw_th(33)
  REAL :: MCC_mlw_o3(33)
  REAL :: MCC_mlw_exner(33)

  ! Mcclatchy Mid-latitude summer profile
  !======================================
  REAL :: MCC_mls_z(33)
  REAL :: MCC_mls_p(33)
  REAL :: MCC_mls_t(33)
  REAL :: MCC_mls_q(33)
  REAL :: MCC_mls_th(33)
  REAL :: MCC_mls_o3(33)
  REAL :: MCC_mls_exner(33)


  ! Arrays to be allocated and used by the SCM
  !===========================================
  REAL, ALLOCATABLE :: &
    MCC_rh_p(:)        &
  , MCC_rh_exner(:)    &
  , MCC_rh_t(:)        &

  , MCC_th_p(:)        &
  , MCC_th_exner(:)    &
  , MCC_th_t(:)        &

  , MCC_th(:)          &
  , MCC_q(:)           &
  , MCC_o3(:)



  ! McClatchey standard atmospheric profiles
  !=========================================

  !===================
  ! TROPICAL LATITUDE
  !===================
  ! z (m)
  DATA MCC_trp_z /         0.0,                                               &
         1.0e+3,        2.0e+3,        3.0e+3,        4.0e+3,        5.0e+3,  &
         6.0e+3,        7.0e+3,        8.0e+3,        9.0e+3,       10.0e+3,  &
        11.0e+3,       12.0e+3,       13.0e+3,       14.0e+3,       15.0e+3,  &
        16.0e+3,       17.0e+3,       18.0e+3,       19.0e+3,       20.0e+3,  &
        21.0e+3,       22.0e+3,       23.0e+3,       24.0e+3,       25.0e+3,  &
        26.0e+3,       27.0e+3,       28.0e+3,       29.0e+3,       30.0e+3,  &
        31.0e+3,       32.0e+3,       33.0e+3,       34.0e+3,       35.0e+3,  &
        36.0e+3,       37.0e+3,       38.0e+3,       39.0e+3,       40.0e+3,  &
        41.0e+3,       42.0e+3,       43.0e+3,       44.0e+3,       45.0e+3,  &
        46.0e+3,       47.0e+3,       48.0e+3,       49.0e+3,       50.0e+3,  &
        51.0e+3,       52.0e+3,       53.0e+3,       54.0e+3,       55.0e+3,  &
        56.0e+3,       57.0e+3,       58.0e+3,       59.0e+3,       60.0e+3,  &
        70.0e+3,      100.0e+3 /


  ! Pressure (Pa)
  DATA MCC_trp_p /    101300.0,                                               &
     90400.0000,    80500.0000,    71500.0000,    63300.0000,    55900.0000,  &
     49200.0000,    43200.0000,    37800.0000,    32900.0000,    28600.0000,  &
     24700.0000,    21300.0000,    18200.0000,    15600.0000,    13200.0000,  &
     11100.0000,     9370.0000,     7890.0000,     6660.0000,     5650.0000,  &
      4800.0000,     4090.0000,     3500.0000,     3000.0000,     2570.0000,  &
      2208.2500,     1902.3100,     1640.0600,     1413.9100,     1220.0000,  &
      1056.1100,      916.3230,      795.5350,      690.5970,      600.0000,  &
       522.9200,      456.7370,      399.1650,      348.8040,      305.0000,  &
       267.4090,      234.8470,      206.2800,      181.0880,      159.0000,  &
       139.8950,      123.2810,      108.7530,       96.1237,       85.4000,  &
        76.6229,       69.3784,       63.2475,       57.9262,       53.1941,  &
        48.8917,       44.9037,       41.1476,       37.5647,       34.1139,  &
         5.7900,        0.0300 /


  ! Temperature (K)
  DATA MCC_trp_t /     300.000,                                               &
        294.000,       288.000,       284.000,       277.000,       270.000,  &
        264.000,       257.000,       250.000,       244.000,       237.000,  &
        230.000,       224.000,       217.000,       210.000,       204.000,  &
        197.000,       195.000,       199.000,       203.000,       207.000,  &
        211.000,       215.000,       217.000,       219.000,       221.000,  &
        223.071,       225.244,       227.478,       229.740,       232.000,  &
        234.235,       236.446,       238.641,       240.824,       243.000,  &
        245.172,       247.349,       249.538,       251.751,       254.000,  &
        256.291,       258.586,       260.839,       262.997,       265.000,  &
        266.777,       268.268,       269.390,       270.025,       270.000,  &
        269.172,       267.693,       265.742,       263.450,       260.910,  &
        258.193,       255.349,       252.417,       249.426,       246.399,  &
        219.000,       210.000 /


  ! Air density (kg/m^3)
  DATA MCC_trp_rho_air / 1.167e+0,                                            &
       1.064e+0,      9.689e-1,      8.756e-1,      7.951e-1,      7.199e-1,  &
       6.501e-1,      5.855e-1,      5.258e-1,      4.708e-1,      4.202e-1,  &
       3.740e-1,      3.316e-1,      2.929e-1,      2.578e-1,      2.260e-1,  &
       1.972e-1,      1.676e-1,      1.382e-1,      1.145e-1,      9.515e-2,  &
       7.938e-2,      6.645e-2,      5.618e-2,      4.763e-2,      4.045e-2,  &
     3.44680e-2,    2.94284e-2,    2.51289e-2,    2.14436e-2,      1.831e-2,  &
     1.56925e-2,    1.34877e-2,    1.16046e-2,    9.98552e-3,      8.600e-3,  &
     7.43004e-3,    6.43291e-3,    5.57205e-3,    4.82499e-3,      4.181e-3,  &
     3.63425e-3,    3.16588e-3,    2.75939e-3,    2.40469e-3,      2.097e-3,  &
     1.83394e-3,    1.60773e-3,    1.41191e-3,    1.24318e-3,      1.101e-3,  &
     9.85472e-4,    8.90852e-4,    8.11422e-4,    7.43033e-4,    6.82678e-4,  &
     6.28179e-4,    5.77963e-4,    5.30898e-4,    4.86178e-4,    4.43241e-4,  &
           RMDI,          RMDI /


  ! Water content (kg/m^3)
  DATA MCC_trp_rho_h2o / 1.90e-2,                                             &
        1.30e-2,       9.30e-3,       4.70e-3,       2.20e-3,       1.50e-3,  &
        8.50e-4,       4.70e-4,       2.50e-4,       1.20e-4,       5.00e-5,  &
        1.70e-5,       6.00e-6,       1.80e-6,       1.00e-6,       7.60e-7,  &
        6.40e-7,       5.60e-7,       5.00e-7,       4.90e-7,       4.50e-7,  &
        5.10e-7,       5.10e-7,       5.40e-7,       6.00e-7,       6.70e-7,  &
     6.79450e-7,    6.31733e-7,    5.49854e-7,    4.53289e-7,       3.60e-7,  &
     2.84044e-7,    2.24225e-7,    1.76770e-7,    1.39187e-7,       1.10e-7,  &
     8.82139e-8,    7.20326e-8,    5.98502e-8,    5.04737e-8,       4.30e-8,  &
     3.67682e-8,    3.14295e-8,    2.67728e-8,    2.26579e-8,       1.90e-8,  &
     1.57512e-8,    1.28621e-8,    1.03099e-8,    8.10658e-9,       6.30e-9,  &
     4.93467e-9,    3.92203e-9,    3.16550e-9,    2.59498e-9,    2.15954e-9,  &
     1.82223e-9,    1.55622e-9,    1.34205e-9,    1.16563e-9,    1.01672e-9,  &
           rmdi,          rmdi /


  ! Ozone density (kg/m^3)
  DATA MCC_trp_rho_o3 / 5.60e-08,                                             &
        5.60e-8,       5.40e-8,       5.10e-8,       4.70e-8,       4.50e-8,  &
        4.30e-8,       4.10e-8,       3.90e-8,       3.90e-8,       3.90e-8,  &
        4.10e-8,       4.30e-8,       4.50e-8,       4.50e-8,       4.70e-8,  &
        4.70e-8,       6.90e-8,       9.00e-8,       1.40e-7,       1.90e-7,  &
        2.40e-7,       2.80e-7,       3.20e-7,       3.40e-7,       3.40e-7,  &
     3.30075e-7,    3.13946e-7,    2.93019e-7,    2.68201e-7,       2.40e-7,  &
     2.08931e-7,    1.76652e-7,    1.45052e-7,    1.16090e-7,       9.20e-8,  &
     7.45581e-8,    6.23765e-8,    5.36522e-8,    4.69395e-8,       4.10e-8,  &
     3.49216e-8,    2.87668e-8,    2.28363e-8,    1.74539e-8,       1.30e-8,  &
     9.78389e-9,    7.57047e-9,    6.06109e-9,    5.02979e-9,       4.30e-9,  &
     3.74043e-9,    3.29423e-9,    2.93037e-9,    2.62635e-9,    2.36586e-9,  &
     2.13705e-9,    1.93136e-9,    1.74263e-9,    1.56644e-9,    1.39966e-9,  &
           rmdi,          rmdi /


  ! Specific humidity (kg/kg)
  ! Calculated from MCC_trp_rho_h20, MCC_trp_rho_air, upper two values at
  ! 70km,100km obtained from separate MMC file
  DATA MCC_trp_q /  1.60202e-2,                                               &
     1.20706e-2,    9.50726e-3,    5.33909e-3,    2.75931e-3,    2.07929e-3,  &
     1.30578e-3,    8.02089e-4,    4.75240e-4,    2.54820e-4,    1.18977e-4,  &
     4.54525e-5,    1.80938e-5,    6.14540e-6,    3.87896e-6,    3.36282e-6,  &
     3.24543e-6,    3.34128e-6,    3.61793e-6,    4.27946e-6,    4.72935e-6,  &
     6.42475e-6,    7.67488e-6,    9.61187e-6,    1.25969e-5,    1.65634e-5,  &
     1.97121e-5,    2.14663e-5,    2.18809e-5,    2.11382e-5,    1.96610e-5,  &
     1.81003e-5,    1.66241e-5,    1.52325e-5,    1.39387e-5,    1.27905e-5,  &
     1.18725e-5,    1.11974e-5,    1.07410e-5,    1.04608e-5,    1.02845e-5,  &
     1.01170e-5,    9.92747e-6,    9.70234e-6,    9.42229e-6,    9.06048e-6,  &
     8.58865e-6,    8.00010e-6,    7.30204e-6,    6.52080e-6,    5.72204e-6,  &
     5.00739e-6,    4.40254e-6,    3.90116e-6,    3.49240e-6,    3.16333e-6,  &
     2.90080e-6,    2.69259e-6,    2.52788e-6,    2.39753e-6,    2.29382e-6,  &
     3.24647e-6,    3.26000e-6 /


  ! Specific Ozone (kg/kg)
  ! Calculated from MCC_trp_rho_o3, MCC_trp_rho_air, upper two values at
  ! 70km,100km obtained from separate MMC file
  DATA MCC_trp_o3 / 4.79863e-8,                                               &
     5.26316e-8,    5.57333e-8,    5.82458e-8,    5.91121e-8,    6.25087e-8,  &
     6.61437e-8,    7.00256e-8,    7.41727e-8,    8.28377e-8,    9.28129e-8,  &
     1.09626e-7,    1.29674e-7,    1.53636e-7,    1.74554e-7,    2.07965e-7,  &
     2.38337e-7,    4.11694e-7,    6.51230e-7,    1.22271e-6,    1.99684e-6,  &
     3.02342e-6,    4.21368e-6,    5.69594e-6,    7.13831e-6,    8.40537e-6,  &
     9.57618e-6,    1.06680e-5,    1.16605e-5,    1.25071e-5,    1.31074e-5,  &
     1.33139e-5,    1.30971e-5,    1.24994e-5,    1.16257e-5,    1.06976e-5,  &
     1.00346e-5,    9.69637e-6,    9.62872e-6,    9.72832e-6,    9.80617e-6,  &
     9.60893e-6,    9.08643e-6,    8.27578e-6,    7.25822e-6,    6.19929e-6,  &
     5.33487e-6,    4.70877e-6,    4.29281e-6,    4.04589e-6,    3.90553e-6,  &
     3.79556e-6,    3.69783e-6,    3.61139e-6,    3.53462e-6,    3.46554e-6,  &
     3.40196e-6,    3.34165e-6,    3.28241e-6,    3.22194e-6,    3.15778e-6,  &
     9.33767e-7,    8.60000e-8 /


  ! Theta (K)
  ! Calculated from MCC_trp_t and MCC_trp_p
  DATA MCC_trp_th /     298.895,                                              &
        302.598,        306.407,      312.559,       315.648,       318.793,  &
        323.285,        326.623,      330.079,       335.189,       338.861,  &
        342.915,        348.399,      353.020,       357.010,       363.759,  &
        369.099,        383.467,      411.029,       440.086,       470.341,  &
        502.283,        535.748,      565.335,       596.227,       628.855,  &
        662.856,        698.438,      735.893,       775.384,       816.706,  &
        859.258,        903.266,      949.211,       997.389,      1047.647,  &
       1099.349,       1152.818,     1208.649,      1267.254,      1328.534,  &
       1391.836,       1457.358,     1525.536,      1596.458,      1669.507,  &
       1743.290,       1817.492,     1891.641,      1964.145,      2031.452,  &
       2088.937,       2137.236,     2178.473,      2214.583,      2247.273,  &
       2278.093,       2308.424,     2339.569,      2372.791,      2409.404,  &
       3553.987,      15321.322 /


  !============================================================================
  ! Mid Latitude Winter
  !============================================================================
  ! z (m)
  DATA MCC_mlw_z /                                                            &
         0.0e+0,         1.0e+3,       2.0e+3,        3.0e+3,        4.0e+3,  &
         5.0e+3,         6.0e+3,       7.0e+3,        8.0e+3,        9.0e+3,  &
         1.0e+4,         1.1e+4,       1.2e+4,        1.3e+4,        1.4e+4,  &
         1.5e+4,         1.6e+4,       1.7e+4,        1.8e+4,        1.9e+4,  &
         2.0e+4,         2.1e+4,       2.2e+4,        2.3e+4,        2.4e+4,  &
         2.5e+4,         3.0e+4,       3.5e+4,        4.0e+4,        4.5e+4,  &
         5.0e+4,         7.0e+4,       1.0e+5/


  ! Pressure (Pa)
  DATA MCC_mlw_p /                                                            &
       1.018e+5,       8.973e+4,     7.897e+4,      6.938e+4,      6.081e+4,  &
       5.313e+4,       4.627e+4,     4.016e+4,      3.473e+4,      2.992e+4,  &
       2.568e+4,       2.199e+4,     1.882e+4,      1.610e+4,      1.378e+4,  &
       1.178e+4,       1.007e+4,     8.610e+3,      7.350e+3,      6.280e+3,  &
       5.370e+3,       4.580e+3,     3.910e+3,      3.340e+3,      2.860e+3,  &
       2.430e+3,       1.110e+3,     5.180e+2,      2.530e+2,      1.290e+2,  &
       6.820e+1,       4.670e+0,     3.000e-2/


  ! Temperature (K)
  DATA MCC_mlw_t /                                                            &
       2.722e+2,       2.687e+2,     2.652e+2,      2.617e+2,      2.557e+2,  &
       2.497e+2,       2.437e+2,     2.377e+2,      2.317e+2,      2.257e+2,  &
       2.197e+2,       2.192e+2,     2.187e+2,      2.182e+2,      2.177e+2,  &
       2.172e+2,       2.167e+2,     2.162e+2,      2.157e+2,      2.152e+2,  &
       2.152e+2,       2.152e+2,     2.152e+2,      2.152e+2,      2.152e+2,  &
       2.152e+2,       2.174e+2,     2.278e+2,      2.432e+2,      2.585e+2,  &
       2.657e+2,       2.307e+2,     2.102e+2/


  ! Water mass Fraction - Assuming kg/(kg of air) i.e. specific humidity
  DATA MCC_mlw_q /                                                            &
     2.69023e-3,     2.15146e-3,   1.73577e-3,    1.25677e-3,    8.33132e-4,  &
     5.10052e-4,     2.85757e-4,   1.45599e-4,    6.70241e-5,    3.46395e-5,  &
     1.84184e-5,     1.27002e-5,   9.06968e-6,    6.68740e-6,    5.12239e-6,  &
     4.04232e-6,     4.00000e-6,   3.99855e-6,    3.99831e-6,    4.00196e-6,  &
     4.00460e-6,     4.00215e-6,   4.00757e-6,    4.00738e-6,    4.00086e-6,  &
     4.00000e-6,     3.99887e-6,   4.00050e-6,    4.00000e-6,    3.99770e-6,  &
     3.99821e-6,     3.99943e-6,   4.00000e-6/


  ! Ozone mass Fraction - Assuming kg/(kg of air) i.e. specific ozone
  DATA MCC_mlw_o3 /                                                           &
     4.61183e-8,     4.64716e-8,   4.72516e-8,    5.30877e-8,    5.91644e-8,  &
     7.82620e-8,     9.67644e-8,   1.30818e-7,    1.72347e-7,    2.59796e-7,  &
     3.92927e-7,     6.00686e-7,   8.66955e-7,    1.16640e-6,    1.45058e-6,  &
     1.79894e-6,     2.22222e-6,   2.80979e-6,    3.45117e-6,    4.22812e-6,  &
     5.17836e-6,     5.79436e-6,   6.78447e-6,    7.20221e-6,    7.78546e-6,  &
     8.60759e-6,     1.06561e-5,   1.16102e-5,    1.13103e-5,    7.46697e-6,  &
     4.80232e-6,     1.21968e-6,   8.60000e-8/

  ! Potential temperature K, calculated from mcc_mlw_t, mcc_mlw_p
  DATA MCC_mlw_th /                                                           &
     2.70817e+2,     2.77147e+2,   2.83701e+2,    2.90503e+2,    2.94735e+2,  &
     2.99135e+2,     3.03706e+2,   3.08457e+2,    3.13409e+2,    3.18573e+2,  &
     3.23939e+2,     3.37843e+2,   3.52398e+2,    3.67623e+2,    3.83449e+2,  &
     4.00093e+2,     4.17461e+2,   4.35554e+2,    4.54635e+2,    4.74429e+2,  &
     4.96123e+2,     5.19192e+2,   5.43185e+2,    5.68189e+2,    5.93934e+2,  &
     6.22226e+2,     7.86245e+2,   1.02422e+3,    1.34180e+3,    1.72878e+3,  &
     2.13173e+3,     3.98094e+3,   1.53359e+4/



  !===========================================================================
  ! Mid Latitude Summer
  !===========================================================================
  ! z (m)
  DATA MCC_mls_z /                                                            &
        0.00e+0,        1.00e+3,      2.00e+3,       3.00e+3,       4.00e+3,  &
        5.00e+3,        6.00e+3,      7.00e+3,       8.00e+3,       9.00e+3,  &
        1.00e+4,        1.10e+4,      1.20e+4,       1.30e+4,       1.40e+4,  &
        1.50e+4,        1.60e+4,      1.70e+4,       1.80e+4,       1.90e+4,  &
        2.00e+4,        2.10e+4,      2.20e+4,       2.30e+4,       2.40e+4,  &
        2.50e+4,        3.00e+4,      3.50e+4,       4.00e+4,       4.50e+4,  &
        5.00e+4,        7.00e+4,      1.04e+5/


  ! Pressure (Pa)
  DATA MCC_mls_p /                                                            &
       1.013e+5,       9.020e+4,     8.020e+4,      7.100e+4,      6.280e+4,  &
       5.540e+4,       4.870e+4,     4.260e+4,      3.720e+4,      3.240e+4,  &
       2.810e+4,       2.430e+4,     2.090e+4,      1.790e+4,      1.530e+4,  &
       1.300e+4,       1.100e+4,     9.500e+3,      8.120e+3,      6.950e+3,  &
       5.950e+3,       5.100e+3,     4.370e+3,      3.760e+3,      3.220e+3,  &
       2.770e+3,       1.320e+3,     6.520e+2,      3.330e+2,      1.760e+2,  &
       9.510e+1,       6.710e+0,     3.000e-2/

  ! Temperature (K)
  DATA MCC_mls_t /                                                            &
       2.940e+2,       2.900e+2,     2.850e+2,      2.790e+2,      2.730e+2,  &
       2.671e+2,       2.610e+2,     2.547e+2,      2.482e+2,      2.417e+2,  &
       2.352e+2,       2.288e+2,     2.223e+2,      2.169e+2,      2.158e+2,  &
       2.158e+2,       2.158e+2,     2.158e+2,      2.160e+2,      2.170e+2,  &
       2.182e+2,       2.194e+2,     2.206e+2,      2.218e+2,      2.230e+2,  &
       2.242e+2,       2.342e+2,     2.453e+2,      2.575e+2,      2.697e+2,  &
       2.762e+2,       2.191e+2,     2.099e+2/


  ! Water mass Fraction - Assuming kg/(kg of air) i.e. specific humidity
  DATA MCC_mls_q /                                                            &
     1.15285e-2,     8.50953e-3,   5.93178e-3,    3.85394e-3,    2.35280e-3,  &
     1.38199e-3,     9.35976e-4,   6.36296e-4,    4.02018e-4,    2.52608e-4,  &
     1.54460e-4,     5.91849e-5,   1.97224e-5,    5.77370e-6,    4.02832e-6,  &
     4.00247e-6,     3.99814e-6,   4.00347e-6,    4.00101e-6,    3.99715e-6,  &
     4.00002e-6,     4.00085e-6,   3.99921e-6,    3.99601e-6,    3.99563e-6,  &
     3.99600e-6,     3.99782e-6,   3.99570e-6,    3.99527e-6,    3.99828e-6,  &
     4.00152e-6,     4.00210e-6,   3.99655e-6/


  ! Ozone mass Fraction - Assuming kg/(kg of air) i.e. specific ozone
  DATA MCC_mls_o3 /                                                           &
     4.94078e-8,     5.49002e-8,   6.08388e-8,    6.96631e-8,    7.96714e-8,  &
     9.12116e-8,     1.06047e-7,   1.28631e-7,    1.51235e-7,    1.84104e-7,  &
     2.16197e-7,     2.97276e-7,   3.66360e-7,    5.21720e-7,    7.28741e-7,  &
     9.05321e-7,     1.18255e-6,   1.56488e-6,    2.13794e-6,    2.86791e-6,  &
     3.57896e-6,     4.44539e-6,   5.21636e-6,    5.75696e-6,    6.36121e-6,  &
     6.96976e-6,     1.01855e-5,   9.93524e-6,    9.10034e-6,    5.71811e-6,  &
     3.58470e-6,     8.06044e-7,   8.63576e-8/


  ! Potential temperature K, calculated from mcc_mls_t, mcc_mls_p
  DATA MCC_mls_th /                                                          &
     2.92917e+2,     2.98670e+2,   3.03539e+2,    3.07672e+2,    3.11795e+2, &
     3.16179e+2,     3.20545e+2,   3.24996e+2,    3.29203e+2,    3.33485e+2, &
     3.37986e+2,     3.42721e+2,   3.47633e+2,    3.54537e+2,    3.68910e+2, &
     3.86481e+2,     4.05369e+2,   4.22704e+2,    4.42495e+2,    4.64744e+2, &
     4.88517e+2,     5.13314e+2,   5.39404e+2,    5.66134e+2,    5.94968e+2, &
     6.24450e+2,     8.06106e+2,   1.03275e+3,    1.31347e+3,    1.65053e+3, &
     2.01522e+3,     3.40896e+3,   1.53140e+4/


CONTAINS
!=============================================================================
  SUBROUTINE get_mcc(latitude, month)

    IMPLICIT NONE

    REAL,    INTENT(In) :: latitude(:,:)
    INTEGER, INTENT(In) :: month
    INTEGER :: k

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('GET_MCC',zhook_in,zhook_handle)

    IF (ABS(latitude(1,1)) <= 30.0) THEN

      ! Tropical latitude
      CALL interp1d(mcc_trp_z, mcc_trp_q,  z_th, mcc_q)
      CALL interp1d(mcc_trp_z, mcc_trp_o3, z_th, mcc_o3)

      ! Get Temperature profile on rho/theta levels
      CALL interp1d(mcc_trp_z, mcc_trp_t,  z_rh, mcc_rh_t)
      CALL interp1d(mcc_trp_z, mcc_trp_t,  z_th, mcc_th_t)

      ! Exner MCC
      mcc_trp_exner(:) = (mcc_trp_p(:)/p_zero)**(kappa)

      ! linearly interpolate for exner on UM theta levels
      CALL interp1d(mcc_trp_z(:), mcc_trp_exner(:),           &
                    z_th(:), mcc_th_exner(:))
      CALL interp1d(mcc_trp_z(:), mcc_trp_exner(:),           &
                    z_rh(:), mcc_rh_exner(:))

    ELSE IF ( (ABS(latitude(1,1)) >  30.0) .AND.              &
              (ABS(latitude(1,1)) <= 60.0)) THEN

      ! Mid-latitude
      ! At some stage, should do linear tranisiton from one profile to
      ! another, rather than a step change.

      SELECT CASE(month)

      CASE (9:12,1:2) ! Winter (includes autumn)

        ! Specific Humidity/Ozone on theta levels
        CALL interp1d(mcc_mlw_z, mcc_mlw_q,  z_th, mcc_q)
        CALL interp1d(mcc_mlw_z, mcc_mlw_o3, z_th, mcc_o3)

        ! Get Temperature profile on rho/theta levels
        CALL interp1d(mcc_mlw_z, mcc_mlw_t,  z_rh, mcc_rh_t)
        CALL interp1d(mcc_mlw_z, mcc_mlw_t,  z_th, mcc_th_t)

        ! Exner MCC
        mcc_mlw_exner(:) = (mcc_mlw_p(:)/p_zero)**(kappa)

        ! linearly interpolate for exner on UM theta levels
        CALL interp1d(mcc_mlw_z(:), mcc_mlw_exner(:),           &
                      z_th(:),      mcc_th_exner(:))
        CALL interp1d(mcc_mlw_z(:), mcc_mlw_exner(:),           &
                      z_rh(:),      mcc_rh_exner(:))

      CASE (3:8) ! Summer (includes spring)

        ! Specific Humidity/Ozone on theta levels
        CALL interp1d(mcc_mls_z, mcc_mls_q,  z_th, mcc_q)
        CALL interp1d(mcc_mls_z, mcc_mls_o3, z_th, mcc_o3)

        ! Get Temperature profile on rho/theta levels
        CALL interp1d(mcc_mls_z, mcc_mls_t,  z_rh, mcc_rh_t)
        CALL interp1d(mcc_mls_z, mcc_mls_t,  z_th, mcc_th_t)

        ! Exner mcc
        mcc_mls_exner(:) = (mcc_mls_p(:)/p_zero)**(kappa)

        ! linearly interpolate for exner on UM theta levels
        CALL interp1d(mcc_mls_z(:), mcc_mls_exner(:),           &
                      z_th(:),      mcc_th_exner(:))
        CALL interp1d(mcc_mls_z(:), mcc_mls_exner(:),           &
                      z_rh(:),      mcc_rh_exner(:))

!!$          ! Get Pressure profile on rho/theta levels
!!$          Call interp1d(mcc_mls_z, mcc_mls_p,  z_rh, mcc_rh_p, .true.,    &
!!$                        old_t = mcc_mls_t,                                &
!!$                        new_t = mcc_rh_t)
!!$
!!$          Call interp1d(mcc_mls_z, mcc_mls_p,  z_th, mcc_th_p, .true.,    &
!!$                        old_t = mcc_mls_t,                                &
!!$                        new_t = mcc_th_t)


      END SELECT

    ELSE
      ! Sub-arctic
    END IF



    ! Calculate MCC theta
    WHERE (mcc_rh_exner /= rmdi)                                              &
      mcc_rh_p = p_zero*mcc_rh_exner**recip_kappa

    WHERE ((mcc_th_exner /= rmdi).AND.(mcc_th_t /= rmdi))                     &
      mcc_th = mcc_th_t / mcc_th_exner

!!$    Do k=1, nmod_lv
!!$      mcc_th(k) = mcc_th_t(k)                                             &
!!$                * ((1.0e+5/mcc_th_p(k))**(R/cp))
!!$    End Do

    IF (lhook) CALL dr_hook('GET_MCC',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE get_mcc

!=============================================================================

  SUBROUTINE alloc_mcc(n_new_lv)

    IMPLICIT NONE

    INTEGER :: n_new_lv

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('ALLOC_MCC',zhook_in,zhook_handle)

    ALLOCATE( MCC_rh_p     (n_new_lv+1), &
              MCC_rh_exner (n_new_lv+1), &
              MCC_rh_t     (n_new_lv+1))

    ALLOCATE( MCC_th_p     (n_new_lv), &
              MCC_th_exner (n_new_lv), &
              MCC_th_t     (n_new_lv) )

    ALLOCATE( MCC_th (n_new_lv), &
              MCC_q  (n_new_lv), &
              MCC_o3 (n_new_lv) )

    IF (lhook) CALL dr_hook('ALLOC_MCC',zhook_out,zhook_handle)

  END SUBROUTINE alloc_mcc

!=============================================================================

  SUBROUTINE dealloc_mcc

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('DEALLOC_MCC',zhook_in,zhook_handle)

    DEALLOCATE( MCC_rh_p,  MCC_rh_exner, MCC_rh_t, &
                MCC_th_p,  MCC_th_exner, MCC_th_t, &
                MCC_th,    MCC_q,        MCC_o3)

    IF (lhook) CALL dr_hook('DEALLOC_MCC',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE dealloc_mcc

!=============================================================================

END MODULE MCC_data


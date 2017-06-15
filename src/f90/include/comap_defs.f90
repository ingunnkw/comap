module comap_defs
  use healpix_types
  implicit none

  integer(i4b), parameter :: SM_CES = 1, SM_CIRC = 2, SM_LISSA=3

  integer(i4b), parameter :: ST_NUM = 32, ST_MJD = 1, ST_AZ = 2, ST_EL = 3, &
   & ST_DK = 4, ST_LST = 5, ST_PWV = 6, ST_PWV_CHANGE = 7, &
   & ST_HUMIDITY = 8, ST_HUMIDITY_CHANGE = 9, ST_WIND = 10, &
   & ST_TENC = 12, ST_TENC_CHANGE = 13, ST_CRYO = 14, ST_CRYO_CHANGE = 15, &
   & ST_T_AMBIENT = 16, ST_MAX = 16

  integer(i4b), parameter :: NUM_DET_STATS = 32, DET_STAT_TYPEB = 1, &
       & DET_STAT_WEATHER1 = 2, DET_STAT_WEATHER2 = 3, DET_STAT_SIGMA0 = 4, &
       & DET_STAT_ALPHA = 5, DET_STAT_FKNEE = 6, DET_STAT_JUMP = 7, &
       & DET_STAT_BIAS = 8, DET_STAT_SSS = 9, DET_STAT_10HZ = 10, &
       & DET_STAT_GAIN = 11, DET_STAT_1_2HZ = 12, DET_STAT_SIDELOBE_HIT = 13, &
       & DET_STAT_SIDELOBE_AZ = 14, DET_STAT_SIDELOBE_EL = 15, &
       & DET_STAT_NU_LOW = 16, DET_STAT_NU_HIGH = 17, DET_STAT_AZORDER = 18, &
       & DET_STAT_MAX = 18

  integer(i4b), parameter :: NUM_FILTR_PAR = 5, FILTR_HIGH_NU_SCAN = 1, &
       & FILTR_HIGH_ALPHA = 2, FILTR_LOW_NU = 3, FILTR_LOW_ALPHA = 4, &
       & FILTR_AZ_ORDER = 5

  real(dp), parameter :: fwhm2sigma = 1d0/sqrt(8d0*log(2d0))

end module

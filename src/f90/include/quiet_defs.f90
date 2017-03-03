module quiet_defs
  use healpix_types
  implicit none

  integer(i4b), parameter :: STAT_NUM = 32, STAT_MJD = 1, STAT_AZ = 2, STAT_EL = 3, &
   & STAT_DK = 4, STAT_LST = 5, STAT_PWV = 6, STAT_PWV_CHANGE = 7, &
   & STAT_HUMIDITY = 8, STAT_HUMIDITY_CHANGE = 9, STAT_WIND = 10, &
   & STAT_TENC = 12, STAT_TENC_CHANGE = 13, STAT_CRYO = 14, STAT_CRYO_CHANGE = 15, &
   & STAT_T_AMBIENT = 16, STAT_MAX = 16

  integer(i4b), parameter :: NUM_DIODE_STATS = 32, DIODE_STAT_TYPEB = 1, &
       & DIODE_STAT_WEATHER1 = 2, DIODE_STAT_WEATHER2 = 3, DIODE_STAT_SIGMA0 = 4, &
       & DIODE_STAT_ALPHA = 5, DIODE_STAT_FKNEE = 6, DIODE_STAT_JUMP = 7, &
       & DIODE_STAT_BIAS = 8, DIODE_STAT_SSS = 9, DIODE_STAT_10HZ = 10, &
       & DIODE_STAT_GAIN = 11, DIODE_STAT_1_2HZ = 12, DIODE_STAT_SIDELOBE_HIT = 13, &
       & DIODE_STAT_SIDELOBE_AZ = 14, DIODE_STAT_SIDELOBE_EL = 15, &
       & DIODE_STAT_NU_LOW = 16, DIODE_STAT_NU_HIGH = 17, DIODE_STAT_AZORDER = 18, &
       & DIODE_STAT_MAX = 18

  integer(i4b), parameter :: NUM_FILTER_PAR = 5, FILTER_HIGH_NU_SCAN = 1, &
       & FILTER_HIGH_ALPHA = 2, FILTER_LOW_NU = 3, FILTER_LOW_ALPHA = 4, &
       & FILTER_AZ_ORDER = 5

  real(dp), parameter :: fwhm2sigma = 1d0/sqrt(8d0*log(2d0))

end module

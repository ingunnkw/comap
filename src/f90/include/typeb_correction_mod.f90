module typeb_correction_mod

   use healpix_types    ! dp, i2b, i4b
   use quiet_fileutils  ! for getlun()
   use quiet_module_mod ! for get_num_modules(),get_num_diodes()
   use math_tools
   use l1_read_mod

   implicit none

   integer, private :: nmodules, ndiodes

   real(dp), parameter, private :: volts2bits = 2.0d0**16.0d0
   real(dp), parameter, private :: bits2volts = 1.0d0/volts2bits

! Noise tolerance
   real(dp) typeb_thresh

! Version of Type-B correction to apply
   integer typeb_vers

! Were Type-B parameters read OK from file
   logical(lgt), private :: typeb_init_ok, debug
   character(len=1024), private :: typeb_params_dir

! Specification of keywords in anomaly.txt Type-B parameter files from Buder/Chicago
   integer, parameter :: nkeywords = 7
   character*24, parameter :: keywords(nkeywords) = &
        (/'__DEMOD_SAMPLING_RATE__ ', & ! Sampling rate for noise calculation (Hz; usually 100)
          '__N_ITERATIONS__        ', & ! Iteration count niter for Type-B fitting (usually 3)
          '__HEIGHT_FACTOR__       ', & ! Heights to be multiplied by this factor (usually 2.0)
          '__PREAMP_RMS_FACTOR__   ', & ! Preamp low-pass filter RMS reduction factor
          '__GLITCHING_CHANNELS__  ', & ! Is channel glitching (1) or not (0)?
          '__DETECTOR_RMS__        ', & ! Detector RMS (XXX UNITS XXX)
          '__CENTERS_HEIGHTS__     '/)        ! Glitch parameters ([]=not used):
                                       !       1. Location y0 / bits
                                       !       2. [Error on location dy0 / bits]
                                       !       3. Amplitude (height) A /bits
                                       !       4. [Error on amplitude dA / bits]

   integer, parameter :: nknownglitches = 62 ! Number of glitches in anomaly.txt
   integer, parameter :: nglitches1 = 31     ! Number of glitches before line break
   integer, parameter :: nblanklines = 2     ! Number of lines to jump in glitch list

! Specification of keywords in anomaly.txt Type-B parameter files from Buder/Chicago, for files of version 2
   !integer, parameter :: nkey2words = 8
   !character*24, parameter :: key2words(nkey2words) = &
   !     (/'__DEMOD_SAMPLING_RATE__ ', & ! Sampling rate for noise calculation (Hz; usually 100)
   !       '__N_ITERATIONS__        ', & ! Iteration count niter for Type-B fitting (usually 3)
   !       '__HEIGHT_FACTOR__       ', & ! Heights to be multiplied by this factor (usually 2.0)
   !       '__PREAMP_RMS_FACTOR__   ', & ! Preamp low-pass filter RMS reduction factor
   !       '__GLITCHING_CHANNELS__  ', & ! Is channel glitching (0+) or not (-1)?
                                        ! Integer >= 0 indicates glitch pattern
   !       '__DETECTOR_RMS__        ', & ! Detector RMS (XXX UNITS XXX)
   !       '__CENTERS_HEIGHTS__     ', & ! Glitch parameters ([]=not used):
   !                                    !       1. Location y0 / bits
   !                                    !       2. [Error on location dy0 / bits]
   !                                    !       3. Amplitude (height) A /bits
   !                                    !       4. [Error on amplitude dA / bits]
    !                                   !__C_H__ followed by glitch pattern int
    !      '__FORMAT_VERSION2__     ' /) ! Indicate version 2 format

   ! This parameter applies to version 2 format
   !integer, parameter :: n2knownglitches = 38 ! Number of glitches in anomaly.txt

   ! Number of glitching patterns to look for
   !integer nglitchpatt
   ! v1 - only 1 pattern
   ! v2 - currently 2 patterns (Q)
   ! v2 - there may be 4+ patterns (W) - IMPLEMENT PROPERLY LATER
   !integer, parameter :: n1glitchpattern = 1
   !integer, parameter :: n2glitchpattern = 2

   ! For v2, number of known glitches in each pattern
!  PATTERN    0  1  2...
   ! Have to assume these don't change:
   integer, parameter :: nglitch(2) = (/ 42,38 /)
   integer, parameter :: nglitch_v3(7) = (/ 102,102,102,102,28,102,102 /)

   type typeb_param_struct
      integer   demod_sampling_rate                    ! Hz
      integer   n_iterations                           ! Usually 3
      real(dp)  height_factor                          ! Usually 2.0d0
      real(dp)  multi_glitch_coverage_sigma            ! Usually 5.0
      real(dp), allocatable :: preamp_rms_factor(:,:)  ! Usually < 1
      ! Next array is called __CENTER_HEIGHT_PATTERN__ after v2
      integer,  allocatable :: glitching_channels(:,:) ! Glitching (1) or not (0)
      real(dp), allocatable :: detector_rms_init(:,:)  ! Volts
      logical(lgt), allocatable :: detector_rms_converged(:,:) ! logical
      real(dp), allocatable :: detector_rms(:,:)       ! Volts
      ! Next array will eventually be removed - for back compatability
      real(dp), allocatable :: centers_heights(:,:)    ! Volts
      ! Next array and refs to will eventually be removed
!      real(dp), allocatable :: centers_heights3(:,:,:) ! Volts (v2 addition)
      ! And that one replaced by this one:
      real(dp), allocatable :: centers_heights22(:,:,:) ! Volts (v2,v3)
   end type typeb_param_struct

   type(typeb_param_struct) typeb_params

   logical(lgt) multi

contains

!  function gget_phase_offset(phase) result(offset)
!    implicit none
!    integer(i2b) :: phase(:)
!    integer(i4b) :: offset
!    offset = 0
!    if(size(phase) < 1) return
!    if(phase(1) == 1)   return
!    offset = 1
!  end function


!============================================================

   subroutine correct_typeb(myid,data_in)
      implicit none

      integer(i4b),      intent(in)    :: myid
      type(data_struct), intent(inout) :: data_in

      real(dp), allocatable :: two_five_Hz_data(:)

      integer(i4b) imod,idiode,phase_offset,iter,maxiter,i
      real(dp) pre_rms_two_five_Hz, rms_two_five_Hz, current_delta,x,y
      logical(lgt) verbose,init
 
      verbose=.true.

      if (verbose) write(*,*) 'Starting the B correction...'

!   So to be completely clear:
!
!   (i) We have up, down at 100 Hz
!   (ii) Demodulate so we have (for each pair of
!   singly-modulated points) x = (up - down) / 2.0 at 50
!   Hz. The rms is therefore estimated from
!   *double*-demodulated data.
!   (iii) For adjacent DD points, calculate y = x1 - x2
!   (NO factor of 2), giving data at 25 Hz.
!   __DETECTOR_RMS__ = the sigma (NO factor of sqrt(n))
!   of the distribution of y's.

!The __DETECTOR_RMS__ is in units of bit*sqrt(s) so there are some factors:
!I divide by 2 when calculating y: so y = (x1 - x2) / 2
!The sigma is multiplied by 0.2 to convert the noise from 25 Hz to 1
!Hz.  You might need an additional conversion factor if you start with
!data not in bits.  It doesn't really matter where you put the factors
!of 2 as long as you end up with the result in bit*sqrt(s).


!double_demod(data)                                                     
!rms1 = calc_rms(data)                                                  
!do
!   data2 = data
!   deglitch(data2,rms1)
!   rms2 = calc_rms(data2)
!   if(abs(rms2-rms1) > tol) exit
!   rms1 = rms2
!end do

!      ! Demodulate 100 Hz -> 50 Hz to return array data_out
!      ! Has to be done to the whole array
!      !write(*,*) 'e.g.',data_in%RQ(1)%demod(1,1)
!      phase_offset = gget_phase_offset(data_in%phase)
!      call double_demodulate(data_in=data_in,data_out=data_out,offset=phase_offset)
!      write(*,*) 'e.g.',data_in%RQ(1)%demod(1,1)
!      write(*,*) 'e.g.',data_out%RQ(1)%demod(1,1)
!      write(*,*) 'DDed OK'

      ! Make first guess at noise on uncorrected data
      rms_two_five_Hz=0.0d0
      do imod = 1,nmodules
         do idiode = 1,ndiodes
            if  (typeb_params%glitching_channels(imod,idiode) < 0) cycle
            ! Input is in volts
            call calculate_quad_bit_noise(data_in%RQ(imod)% &
                    demod(idiode,:),rms_two_five_Hz,.false.)
            !write(*,*) '(volts)',imod,idiode,rms_two_five_Hz
            typeb_params%detector_rms(imod,idiode)=rms_two_five_Hz
       !     if (rms_two_five_Hz<1.0e-8) write(*,*) imod,idiode,rms_two_five_Hz
       !     if (rms_two_five_Hz-0.2<1.0e-8) write(*,*) imod,idiode,rms_two_five_Hz
         enddo
      enddo
      write(*,*) 'First guesses at noise made..'

      !write(*,*) typeb_params%detector_rms(:,:)
      !stop

      ! Initialize the noise
      typeb_thresh=0.0002d0
      typeb_thresh=typeb_thresh / (0.2*65536.0)
      write(*,*) 'Thresh is (volts)', typeb_thresh

!double_demod(data)                                                     
!rms_pre = calc_rms(data)                                                  
!do
!   deglitch(data,rms_pre)
!   rms = calc_rms(data)
!   if(abs(rms-rms_pre) < tol) exit
!   rms_pre = rms
!enddo

      maxiter=50
      pre_rms_two_five_Hz=0.0d0
      modules: do imod = 18, 18
         diodes: do idiode = 1, ndiodes
            iter=1
            noise_iter: do while ( &
                 .not. typeb_params%detector_rms_converged(imod,idiode))
               pre_rms_two_five_Hz=typeb_params%detector_rms(imod,idiode)
               write(*,*) 't',typeb_params%detector_rms(imod,idiode)
               ! Skip the diode if not a glitcher...
               if  (typeb_params%glitching_channels(imod,idiode) < 0) then
                  typeb_params%detector_rms_converged(imod,idiode)=.true.
                  write(*,*) 'Non-glitcher m x d',imod,idiode
                  exit
               endif
               ! ...or if already converged - NOT NEEDED??
               if (typeb_params%detector_rms_converged(imod,idiode)) then
           !       write(*,*) 'Skipping converged m x d', imod, idiode
                  exit
               endif

               write(*,*) 'Iteration m x d', imod, idiode, ':', iter
               ! Correct the data for Type B
               ! Input is in volts
               x=sum(data_in%RQ(imod)%demod(idiode-1,:))
               ! typeb_params%detector_rms(imod,idiode) passed implicitly
               multi=.false.
               call apply_typeb_correction_v3(data_in%RQ(imod)% &
                    demod(idiode-1,:),data_in%RQ(imod)% &
                    avg(idiode-1,:),6880_i2b, &
                    myid,imod,idiode,verbose,multi)
               y=sum(data_in%RQ(imod)%demod(idiode-1,:))
               !call assert(y/=x,'Bad - B corr did nothing!')
               write(*,*) 'mod x di B-corrected', imod,idiode
               ! Calculate the 25-Hz rms
               call calculate_quad_bit_noise(data_in%RQ(imod)% &
                    demod(idiode,:),rms_two_five_Hz,.false.)
               typeb_params%detector_rms(imod,idiode)=rms_two_five_Hz
               write(*,*) 't2',typeb_params%detector_rms(imod,idiode)
               current_delta=pre_rms_two_five_Hz-rms_two_five_Hz
               write(*,*) 'pre noise is   ',pre_rms_two_five_Hz
               write(*,*) 'white noise is ',rms_two_five_Hz
               write(*,*) 'delta is       ',current_delta
               write(*,*) 'thresh noise is',typeb_thresh
               !stop

               test_noise: if (abs(current_delta) < typeb_thresh) then
                  write(*,*) 'converged: ',imod,idiode
                  typeb_params%detector_rms_converged(imod,idiode) = .true.
               endif test_noise
               pre_rms_two_five_Hz=rms_two_five_Hz

               iter=iter+1
               if (iter > maxiter) then
                  write(*,*) '***maxiter reached!!'
                  exit
               endif

            enddo noise_iter
         enddo diodes
      enddo modules

      write(*,*) 'Type-B correction finished..'

      return

   end subroutine correct_typeb

!============================================================

   subroutine calculate_quad_bit_noise(upHz_array,rms,verbose)
      implicit none

      real(dp),              intent(in)  :: upHz_array(:) ! Volts
      logical(lgt),          intent(in)  :: verbose
      real(dp),              intent(out) :: rms ! Volts
      real(dp), allocatable :: downHz_array(:)

      ! Difference 50 Hz -> 25 Hz and return 25-Hz data                
      !write(*,*) 'going into second demodulation..'                   
      call demodulate2(upHz_array,downHz_array)

      ! Check the white noise for this diode                           
      call measure_white_noise(downHz_array,rms)
      if (verbose) write(*,*) 'white noise is',rms

      ! Clean up
      deallocate(downHz_array)

      ! We're still in Volts

   end subroutine calculate_quad_bit_noise

!============================================================

   subroutine measure_white_noise(array,rms)
      implicit none

      real(dp), allocatable, intent(in)  :: array(:)
      real(dp),              intent(out) :: rms
      real(dp) n

      n   = real(size(array),dp)
      rms = sqrt(sum(array**2) / n)

      return

   end subroutine measure_white_noise

!============================================================

   subroutine convert_white_noise_to_1Hz(std_two_five_Hz,std_one_Hz)
      implicit none

      real(dp), intent(in)  :: std_two_five_Hz
      real(dp), intent(out) :: std_one_Hz

      std_one_Hz = 0.2 * std_two_five_Hz

      return

   end subroutine convert_white_noise_to_1Hz

!============================================================

   subroutine demodulate2(a,b)
      implicit none

      ! Demodulate e.g. 50-Hz data to 25-Hz data by differencing

      real(dp),              intent(in)   :: a(:)
      real(dp), allocatable, intent(out)  :: b(:)
      integer n

      ! Demodulate down to 25 Hz
      n = size(a)/2

      allocate(b(n))

      !write(*,*) n
      !write(*,*) size(a)
      !write(*,*) size(b)

      b = (a(1:2*n:2) - a(2:2*n:2))/2

      !write(*,*) 'Demodulated to 25 Hz (a=',size(a),'b=',size(b),')'
      !write(*,*) 'a',size(a)
      !write(*,*) 'b',size(b)

      !write(*,*) a%RQ(1)%demod(1,1)
      !write(*,*) 'here2'

      return

   end subroutine demodulate2

!============================================================

   subroutine initialize_typeb_correction_mod(unit,paramfile)
      ! NB this is for v3 only
      implicit none
      integer(i4b)     :: unit
      character(len=*) :: paramfile
      logical(lgt) havetypeb_vers,verbose
      integer(i4b) status

      typeb_init_ok = .false.

      ! Fetch module, diode information
      call initialize_module_mod(paramfile)
      ndiodes = get_num_diodes()
      nmodules = get_num_modules()

      ! And the Type-B parameters
      call get_parameter(unit, paramfile, 'DEBUG_TYPEB',      par_lgt=debug)
      call get_parameter(unit, paramfile, 'TYPEB_PARAMS_DIR', par_string=typeb_params_dir)
      havetypeb_vers = .false.
      call get_parameter(unit, paramfile, 'TYPEB_VERS',       par_int=typeb_vers,par_present=havetypeb_vers)
      !if (.not. havetypeb_vers) typeb_vers = 1 ! Dangerous!
      call get_parameter(unit, paramfile, 'TYPEB_THRESH',       par_dp=typeb_thresh)

      !debug=.true.

      if (debug) write(*,*) 'B info:', debug,typeb_vers

      ! Initialize the expected number of glitching patterns

      ! Number of glitching patterns to look for
      ! And expected length of each glitching pattern
      ! v1 - only 1 pattern
      ! v2 - currently 2 patterns (Q)
      ! v3 - currently 7 patterns (W) - IMPLEMENT PROPERLY LATER
      !nglitchpatt=0
      !if (typeb_vers == 1) then
      !   nglitchpatt=n1glitchpattern
      !elseif (typeb_vers == 2) then
      !   ! I think this is now redundant - stored in nglitch(:)
      !   !nglitchpatt=n2glitchpattern
      !   !nglitch(1)=38 ! IB has switched patterns round..
      !   !nglitch(2)=42
      !else
      !   write(*,*) 'Disallowed or mangled Type-B version', typeb_vers
      !   return
      !endif

      ! W-band initialization is different:
      ! all the parameters are read once from file
      ! and are assumed to be the same for all CESs
      ! The noise is later re-estimated
      !write(*,*) typeb_vers
      !write(*,*) typeb_thresh
      write(*,*) trim(typeb_params_dir)

      if (typeb_vers /= 3) stop

! Initialize typeb_params structure

      allocate(typeb_params%preamp_rms_factor(nmodules,ndiodes))
      allocate(typeb_params%glitching_channels(nmodules,ndiodes))
      allocate(typeb_params%detector_rms_init(nmodules,ndiodes))
      allocate(typeb_params%detector_rms(nmodules,ndiodes))
      allocate(typeb_params%detector_rms_converged(nmodules,ndiodes))
      allocate(typeb_params%centers_heights22(maxval(nglitch_v3),&
           4,0:size(nglitch_v3)-1))

      typeb_params%demod_sampling_rate = 0
      typeb_params%n_iterations = 0
      typeb_params%height_factor = 0.0d0
      typeb_params%preamp_rms_factor = 0.0d0
      typeb_params%glitching_channels = 0
      typeb_params%detector_rms_init = 0.0d0
      typeb_params%detector_rms = 0.0d0
      typeb_params%detector_rms_converged = .false.
      typeb_params%centers_heights22 = 0.0d0

! Now read the params -> struct

      typeb_params_dir=trim(typeb_params_dir) // 'wband/'
      if (debug) write(*,*) 'Reading',trim(typeb_params_dir)

      call read_typeb_params_v3(typeb_params_dir,status,verbose)
      if (debug) write(*,*) 'Type-B params initialized'

      typeb_init_ok = .true.

      return

   end subroutine initialize_typeb_correction_mod

!============================================================

   subroutine read_typeb_params_v3(typeb_params_dir,status,verbose)

      implicit none

      ! Arguments
      character(len=512),       intent(inout) :: typeb_params_dir
      integer(i4b),             intent(out)   :: status

      character(len=512) filename

      logical, optional :: verbose
      logical           :: verbose_

      ! Outputs are returned via globals

      ! Local working variables
      integer(i4b) iunit,imod,iglitch,ipatt,icol,idiode,id,ifile,ik
      integer(i4b), parameter :: lmax = 1000 ! Max number of lines to read
      character(len=80) key
      character(len=32) dummy
      character(len=16), dimension(5) :: files

      ! Aim is to read a set of parameter.txt data files into
      ! one structure of parameters for ALL CESs

      verbose_ = .false.; if(present(verbose)) verbose_ = verbose

      !verbose_=.true.

      ! Check typeb_vers agrees
      if (typeb_vers /= 3) stop

      ! Set up files to be read
      files = (/ 'params.txt  ', 'preamp.txt  ', 'mask.txt    ', &
                 'rms_init.txt', 'patterns.txt' /)

      if (verbose_) write(*,*) 'files ', files

      ! Loop over files
      paramfiles: do ifile = 1, size(files)

         ! Operations common to all files
         filename=''
         if (verbose_) write(*,*) 'dir ', trim(typeb_params_dir)
         if (verbose_) write(*,*) 'file ', trim(files(ifile))
         filename=trim(typeb_params_dir) // trim(files(ifile))
         if (verbose_) write(*,*) 'Reading file ', trim(filename)
         ! Fetch a file handle
         iunit=getlun()
         open(iunit, file=filename, status='OLD', iostat=status)
         if (verbose_) write(*,*) 'status of open ', status
         if (status.ne.0) then
            status = status
            return
         endif
         if (verbose_) write(*,*) 'status of open ', status

         ! Read the single, probably fixed parameters from params.txt
         solo_params: if (ifile == 1) then
            ! Check version (first line)
            read(iunit,'(A)',end=10,err=20) key
            if (verbose_) write(*,*) ifile, key
            if (trim(key)/='__FORMAT_VERSION2__') then
               write(*,*) 'Type-B version problem', trim(key)
               stop
            endif

            lines: do while (status == 0)

               ! Read each line
               ! The order of the data does not in theory matter
               ! This is also the loop's exit - 10 or 20
               read(iunit,'(A)',end=10,err=20) key

               ! Ignore comments and blank lines
               if (key(:1)=='#') cycle
               if (trim(key)=='') cycle

               ! Extract sampling rate
               sampling_rate: if (key(:23)=='__DEMOD_SAMPLING_RATE__') then
                  read(key(24:),*) dummy
                  read(dummy,'(I3)') typeb_params%demod_sampling_rate
                  if (verbose_)write(*,*)'rate  ',typeb_params%demod_sampling_rate
                  cycle
               endif sampling_rate

               ! Extract n_iterations
               niter: if (key(:16)=='__N_ITERATIONS__') then
                  read(key(17:),*) dummy
                  read(dummy,'(I3)') typeb_params%n_iterations
                  if (verbose_) write(*,*) 'niter ', typeb_params%n_iterations
                  cycle
               endif niter

               ! Extract height factor
               height: if (key(1:17)=='__HEIGHT_FACTOR__') then
                  read(key(18:),*) dummy
                  read(dummy,'(F4.2)') typeb_params%height_factor
                  if (verbose_) write(*,*) 'height', typeb_params%height_factor
                  cycle
               endif height

               ! Extract multi glitch coverage sigma [new for W band]
               sigma_multi: if (key(1:31)=='__MULTI_GLITCH_COVERAGE_SIGMA__') then
                  read(key(32:),*) dummy
                  read(dummy,'(F4.2)') typeb_params%multi_glitch_coverage_sigma
                  if (verbose_) write(*,*) 'multi-glitch sigma', &
                       typeb_params%multi_glitch_coverage_sigma
                  cycle
               endif sigma_multi
            enddo lines
         endif solo_params

         ! Read the preamp rms factors
         preamp_rms: if (ifile == 2) then
            ! Find the keyword
            do ik = 1,lmax
               read(iunit,'(A)',end=10,err=20) key
               if (trim(key)=='__PREAMP_RMS_FACTOR__') exit
            enddo

            ! Read the next nmodules lines
            do imod=1,nmodules
               !write(*,*) imod!, shape(typeb_params%preamp_rms_factor)
               read(iunit,*,iostat=status) id, &
                    (typeb_params%preamp_rms_factor(imod,idiode),idiode=1,ndiodes)
               !write(*,'(i2,1X,4(1X,f5.3))') id, (typeb_params%preamp_rms_factor(imod,idiode),idiode=1,ndiodes)
               call assert(id==imod-1,'imod mismatch '//itoa(imod)//' '//itoa(id))
            enddo
            if (verbose_) write(*,*) 'e.g. ',typeb_params%preamp_rms_factor(9,3)
         endif preamp_rms

         ! Extract center height patterns (formerly, glitching_channels)
         mask: if (ifile == 3) then
            ! Find the keyword
            do ik = 1,lmax
               read(iunit,'(A)',end=10,err=20) key
               if (trim(key)=='__CENTER_HEIGHT_PATTERN__') exit ! the loop..
            enddo

            ! Read the next nmodules lines
            do imod=1,nmodules
               read(iunit,*,iostat=status) id, &
               (typeb_params%glitching_channels(imod,idiode),idiode=1,ndiodes)
               !write(*,'(i2,1X,4(1X,i2))') id, (typeb_params%glitching_channels(imod,idiode),idiode=1,ndiodes)
               call assert(id==imod-1,'imod mismatch '//itoa(imod)//' '//itoa(id))
            enddo
            if (verbose_) write(*,*) 'e.g. ',typeb_params%glitching_channels(9,3)
         endif mask

         ! Read detector rms first guesses
         detector_rms: if (ifile == 4) then
            ! Find the keyword
            do ik = 1,lmax
               read(iunit,'(A)',end=10,err=20) key
               if (trim(key)=='__DETECTOR_RMS__') exit
            enddo

            ! Read the next nmodules lines
            do imod = 1, nmodules
               read(iunit,*,iostat=status) id, &
                    (typeb_params%detector_rms_init(imod,idiode),idiode=1,ndiodes)
               !write(*,'(i2,1X,4(1X,f11.10))') id, (typeb_params%detector_rms(imod,idiode),idiode=1,ndiodes)
               call assert(id==imod-1,'imod mismatch '//itoa(imod)//' '//itoa(id))
            enddo
            if (verbose_) write(*,*) 'e.g. ', typeb_params%detector_rms_init(3,3)
            typeb_params%detector_rms=typeb_params%detector_rms_init

         endif detector_rms

         ! Read centers heights (x nglitchpatt)
         patterns:if (ifile == 5) then

            pattlines: do while (status == 0) ! Loop over lines/patterns
               if (verbose_) write(*,*) 'status ',status
               read(iunit,'(A)',iostat=status) key ! Read line
               if (key(:1)=='#') cycle pattlines
               if (trim(key)=='') cycle pattlines

               ! Find the keyword
               !header: do ik = 1,lmax
               !read(iunit,'(A)',end=10,err=20) key
               if (key(1:19)=='__CENTERS_HEIGHTS__') then
                  if (verbose_) write(*,*) 'key ',trim(key)
               endif

               ipatt=-99
               ! Catch the ID of this pattern
               read(key(20:),*) dummy
               read(dummy,'(I2)') ipatt
               !write(*,*) 'ipatt ',ipatt

               ! Check ipatt is acceptable
               if (.not. is_pattern_ok(ipatt)) stop
               ! Read the next nglitch lines for this ipatt
               if (verbose_) write(*,*) 'Read pattern ', ipatt
               glitches: do iglitch = 1, nglitch_v3(ipatt+1)
                  !write(*,*) iglitch,ipatt,nglitch_v3(ipatt+1)
                  read(iunit,'(A)',iostat=status) key
                  !write(*,*) 'k3 ', trim(key)
                  read(key,*,iostat=status) &
                  (typeb_params%centers_heights22(iglitch,icol,ipatt),icol=1,4)
                  if (verbose_) write(*,'(i3,1X,4(1X,f18.10))') iglitch, &
                       (typeb_params%centers_heights22(iglitch,icol,ipatt),icol=1,4)
               enddo glitches
            enddo pattlines

         endif patterns

         ! Close the file handle since the end of the file has been reached
10       close(iunit)

      ! End loop over files
      enddo paramfiles


!      ! There has been an error
20    if (status < 0) write(*,*) &
           'Type-B parameters read OK (last EOF reached)'
      !write(*,*) 'Error reading the type-B file'
      !write(*,*) 'bad file', files(ifile)
      !write(*,*) 'status was',status
!      stop

      if (verbose_) write(*,*) typeb_params%demod_sampling_rate
      if (verbose_) write(*,*) typeb_params%n_iterations
      if (verbose_) write(*,*) typeb_params%height_factor
      if (verbose_) write(*,*) typeb_params%multi_glitch_coverage_sigma

      return

   end subroutine read_typeb_params_v3

!============================================================

   subroutine read_typeb_params(filename,keywords,status,verbose)

      ! Read Type-B parameters from file

      implicit none

      ! Arguments
      character*128,            intent(in)    :: filename
      character*24,             intent(in)    :: keywords(:)
      integer,                  intent(out)   :: status

      logical, optional :: verbose
      logical           :: verbose_

      ! Local working variables
      integer :: keywords_length(size(keywords))
      character*80 line
      character*32 dummy
      integer nkeywords
      integer(i4b) iunit

      ! Very local working variables
      integer ikey, imod, idiode, iglitch, rglitch, icol, iblank
      integer id ! Read in but not used

      verbose_ = .false.; if(present(verbose)) verbose_ = verbose

      ! Check typeb_vers agrees
      if (typeb_vers /= 1) stop

      nkeywords = size(keywords)

      ! Determine length of each keyword
      do ikey = 1, nkeywords
         keywords_length(ikey) = len(trim(keywords(ikey)))
      enddo

      dummy=''

      iunit=getlun()
      open (iunit, file=filename, status='OLD', iostat=status)
      if (status.ne.0) then
         !write(*,*) 'read error'
         status = status
         return
      else
         do while (status.eq.0)
            read(iunit,'(A)',iostat=status) line
            if (line(1:1).ne.'#' .and. trim(line).ne.'') then

                if (line(:keywords_length(1))==keywords(1)) then
                  if (verbose_) write(*,*) 'MATCH 1 ', keywords(1)
                  read(line(keywords_length(1)+1:),*) dummy
                  read(dummy,'(I3)') typeb_params%demod_sampling_rate
                  if (verbose_) write(*,*) 'rate  ', typeb_params%demod_sampling_rate
               endif

               if (line(:keywords_length(2))==keywords(2)) then
                  if (verbose_) write(*,*) 'MATCH 2 ', keywords(2)
                  read(line(keywords_length(2)+1:),*) dummy
                  read(dummy,'(I3)') typeb_params%n_iterations
                  if (verbose_) write(*,*) 'niter ', typeb_params%n_iterations
               endif

               if (line(:keywords_length(3))==keywords(3)) then
                  if (verbose_) write(*,*) 'MATCH 3 ', keywords(3)
                  read(line(keywords_length(3)+1:),*) dummy
                  read(dummy,'(F4.2)') typeb_params%height_factor
                  if (verbose_) write(*,*) 'height', typeb_params%height_factor
               endif

               if (line(:keywords_length(4))==keywords(4)) then
                  if (verbose_) write(*,*) 'MATCH 4 ', keywords(4)
                  do imod = 1, nmodules
                     read(iunit,*,iostat=status) &
                          id, (typeb_params%preamp_rms_factor(imod,idiode),idiode=1,ndiodes)
                  enddo
                  if (verbose_) write(*,*) typeb_params%preamp_rms_factor(9,3)
               endif

               if (line(:keywords_length(5))==keywords(5)) then
                  if (verbose_) write(*,*) 'MATCH 5 ', keywords(5)
                  do imod = 1, nmodules
                     read(iunit,*,iostat=status) &
                          id, (typeb_params%glitching_channels(imod,idiode),idiode=1,ndiodes)
                  enddo
                  if (verbose_) write(*,*) typeb_params%glitching_channels(2,2)
               endif

               if (line(:keywords_length(6))==keywords(6)) then
                  if (verbose_) write(*,*) 'MATCH 6 ', keywords(6)
                  do imod = 1, nmodules
                     read(iunit,*,iostat=status) &
                          id, (typeb_params%detector_rms(imod,idiode),idiode=1,ndiodes)
                  enddo
                  if (verbose_) write(*,*) typeb_params%detector_rms(3,3)
               endif

               if (line(:keywords_length(7))==keywords(7)) then
                  if (verbose_) write(*,*) 'MATCH 7 ', keywords(7)
                  ! This formatting is horrible!
                  rglitch=0

                  do iglitch = 1,nglitches1+nblanklines
                     read(iunit,'(A)',iostat=status) line
                     if (trim(line)=='') then
                        ! If one blank line, read more blank lines (set by nblanklines)
                        do iblank = 1, nblanklines-1
                           read(iunit,'(A)',iostat=status) line
                           if (trim(line).ne.'') then ! Exit if subsequent lines not blank
                              status = 10
                              return
                           endif
                        enddo
                        rglitch=iglitch
                        exit
                     endif
                     read(line,*) (typeb_params%centers_heights(iglitch,icol),icol=1,4)
                  enddo

                  do iglitch = rglitch,nknownglitches
                     read(iunit,'(A)',iostat=status) line
                     read(line,*) (typeb_params%centers_heights(iglitch,icol),icol=1,4)
                  enddo

               endif

            endif
         enddo ! End status loop
         close(iunit)
      endif ! If file opened OK

      ! Handle status cases
      if (status.eq.-1) then
         status = 0 ! -1 is EOF, so indicate read OK
      else
         status = status ! File in unexpected format, so indicate failure
      endif

      return

   end subroutine read_typeb_params


!============================================================

   function is_pattern_ok(ipatt)

      ! Tests measured __CENTER_HEIGHTS__ pattern ID
      ! Check it is either 0, 1 (or 2?)

      implicit none

      integer,      intent(in)  :: ipatt
      logical(lgt) is_pattern_ok
      integer, dimension(7), parameter :: ok_patterns = (/ 0,1,2,3,4,5,6 /)
      !integer, dimension(2), parameter :: ok_patterns = (/ 0,1,2 /)
      integer iavail
      is_pattern_ok=.false.

      do iavail = 1,size(ok_patterns(:))
         if (ipatt == ok_patterns(iavail)) then
            is_pattern_ok = .true.
            exit ! Exit loop
         endif
      enddo

      if (.not. is_pattern_ok) write(*,*) '***ipatt is not OK', ipatt, ok_patterns

   end function is_pattern_ok

!============================================================

   subroutine read_typeb_params_v22(filename,status,verbose)

      implicit none

      ! Arguments
      character*128,            intent(in)    :: filename
      integer,                  intent(out)   :: status

      logical, optional :: verbose
      logical           :: verbose_

      ! Outputs are returned via globals

      ! Local working variables
      integer(i4b) iunit,imod,iglitch,ipatt,icol,idiode,id
      character(len=80) key
      character(len=32) dummy

      ! Aim is to read one anomaly.txt data file into
      ! one structure of parameters for this CES

      verbose_ = .false.; if(present(verbose)) verbose_ = verbose

      ! Check typeb_vers agrees
      if (typeb_vers /= 2) stop

      ! Fetch a file handle
      iunit=getlun()
      open (iunit, file=filename, status='OLD', iostat=status)
      if (status.ne.0) then
         status = status
         return
      endif

      ! Read first line
      ! This is the only line treated uniquely like this
      read(iunit,'(A)',end=10,err=20) key
      ! Check the file version
      if (trim(key)/='__FORMAT_VERSION2__') then
         write(*,*) 'Type-B version problem', trim(key)
         stop
      endif

      ! Now the version is checked, continue reading remaining lines
      do
         ! Read each line
         ! The order of the data does not in theory matter
         ! Usually the file format is the same for all
         read(iunit,'(A)',end=10,err=20) key

         ! Ignore comments and blank lines
         if (key(:1)=='#') cycle
         if (trim(key)=='') cycle

         ! Extract sampling rate
         if (key(:23)=='__DEMOD_SAMPLING_RATE__') then
            read(key(24:),*) dummy
            read(dummy,'(I3)') typeb_params%demod_sampling_rate
            if (verbose_) write(*,*) 'rate  ', typeb_params%demod_sampling_rate
            cycle
         endif

         ! Extract n_iterations
         if (key(:16)=='__N_ITERATIONS__') then
            read(key(17:),*) dummy
            read(dummy,'(I3)') typeb_params%n_iterations
            if (verbose_) write(*,*) 'niter ', typeb_params%n_iterations
            cycle
         endif

         ! Extract height factor
         if (key(1:17)=='__HEIGHT_FACTOR__') then
            read(key(18:),*) dummy
            read(dummy,'(F4.2)') typeb_params%height_factor
            if (verbose_) write(*,*) 'height', typeb_params%height_factor
            cycle
         endif

         ! Extract preamp rms factors
         if (trim(key)=='__PREAMP_RMS_FACTOR__') then
            ! Read the next nmodules lines
            do imod=1,nmodules
               read(iunit,*,iostat=status) id, &
               (typeb_params%preamp_rms_factor(imod,idiode),idiode=1,ndiodes)
               !write(*,'(i2,1X,4(1X,f5.3))') id, (typeb_params%preamp_rms_factor(imod,idiode),idiode=1,ndiodes)
               call assert(id==imod-1,'imod mismatch '//itoa(imod)//' '//itoa(id))
            enddo
            if (verbose_) write(*,*) typeb_params%preamp_rms_factor(9,3)
            cycle
         endif

         ! Extract center height patterns (formerly, glitching_channels)
         if (trim(key)=='__CENTER_HEIGHT_PATTERN__') then
            ! Read the next nmodules lines
            do imod=1,nmodules
               read(iunit,*,iostat=status) id, &
               (typeb_params%glitching_channels(imod,idiode),idiode=1,ndiodes)
               !write(*,'(i2,1X,4(1X,i2))') id, (typeb_params%glitching_channels(imod,idiode),idiode=1,ndiodes)
               call assert(id==imod-1,'imod mismatch '//itoa(imod)//' '//itoa(id))
            enddo
            if (verbose_) write(*,*) typeb_params%glitching_channels(9,3)
            cycle
         endif

         ! Read detector rms
         if (trim(key)=='__DETECTOR_RMS__') then
            ! Read the next nmodules lines
            do imod = 1, nmodules
               read(iunit,*,iostat=status) id, &
                    (typeb_params%detector_rms(imod,idiode),idiode=1,ndiodes)
               !write(*,'(i2,1X,4(1X,f11.10))') id, (typeb_params%detector_rms(imod,idiode),idiode=1,ndiodes)
               call assert(id==imod-1,'imod mismatch '//itoa(imod)//' '//itoa(id))
            enddo
            if (verbose_) write(*,*) typeb_params%detector_rms(3,3)
            cycle
         endif

         ! Read centers heights (x nglitchpatt)
         if (key(1:19)=='__CENTERS_HEIGHTS__') then
            ipatt=-99
            ! Catch the ID of this pattern
            read(key(20:),*) dummy
            read(dummy,'(I2)') ipatt
            ! Check ipatt is acceptable
            if (.not. is_pattern_ok(ipatt)) stop
            ! Read the next nglitch lines for this ipatt
            if (verbose_) write(*,*) 'Read pattern ', ipatt
            do iglitch = 1, nglitch(ipatt+1)
               !write(*,*) iglitch,ipatt,nglitch(ipatt+1)
               read(iunit,*,iostat=status) &
               (typeb_params%centers_heights22(iglitch,icol,ipatt),icol=1,4)
               if (verbose_) write(*,'(i2,1X,4(1X,f18.10))') iglitch, (typeb_params%centers_heights22(iglitch,icol,ipatt),icol=1,4)
            enddo
            cycle
         endif

      ! End loop over lines
      enddo

      ! There has been an error
20    write(*,*) 'There has been an error in reading the type-B file!'
      stop

      ! Close the file handle since the end of the file has been reached
10    close(iunit)
      return

   end subroutine read_typeb_params_v22

!============================================================

   function dG_NewtonIter(dG_prev,y,A,sigma)
      ! Apply one iteration of Type-B correction
      ! Based on Kusaka's C++ implementation
      real(dp) dG_NewtonIter
      real(dp), intent(in) :: dG_prev,y,A,sigma
      ! Local working variables
      real(dp) P_G,erf_y,denom,numer
      real(dp), parameter :: pi = 3.1415926535

      P_G = exp(-((dG_prev+y)**2.0d0)/(2.0d0*sigma*sigma))/(sqrt(2.0d0*pi*sigma*sigma))

      erf_y = erf((dG_prev+y)/(sqrt(2.0d0)*sigma))

      denom = 1.0d0 + (A*P_G)
      numer = (A/2.0d0) * erf_y

      dG_NewtonIter = -(dG_prev + numer) / denom

      return

   end function dG_NewtonIter

!============================================================

   function dG_Newton(y,A,sigma,niter)
      ! Loop over niter iterations of Type-B correction
      ! Based on Kusaka's C++ implementation
      real(dp) dG_Newton
      real(dp), intent(in) :: y,A,sigma
      integer, intent(in) :: niter
      ! Local working variables
      real(dp) dG
      integer iter

      dG = 0.0d0
      do iter = 1,niter
         dG = dG + dG_NewtonIter(dG,y,A,sigma)
      enddo

      dG_Newton = dG

      return

   end function dG_Newton

!============================================================

   subroutine bcorrect(ya,yd,A,y0,sigma,niter,decorrect,avcorrect)

      ! Apply Type-B correction to a single data point
      ! All input/output quantities are in volts

      implicit none

      ! Arguments
      real(dp), intent(inout) :: yd                   ! Demod sample value (volts)
      real(dp), intent(in)    :: ya                   ! Average sample value (volts)
      real(dp), intent(in)    :: A, y0                ! Glitch amplitude (volts) and location (volts)
      real(dp), intent(in)    :: sigma                ! White noise at 800 kHz (volts)
      integer,  intent(in)    :: niter                ! Number of iterations for Newton-Raphson
      logical,  intent(in)    :: decorrect,avcorrect  ! Correct DE, AV or not

      ! Local working variables
      real(dp) dG_1, dG_2
      real(dp) yd_bits,ya_bits
      real(dp) sigma_bits,A_bits,y0_bits

      ! First convert inputs from volts to bits
      yd_bits    = volts2bits * yd
      ya_bits    = volts2bits * ya
      y0_bits    = volts2bits * y0
      A_bits     = volts2bits * A
      sigma_bits = volts2bits * sigma

      if (decorrect) then
         ! Apply Type-B correction with all quantities in bits
         dG_1 = dG_Newton(ya_bits+yd_bits-y0_bits,A_bits,sigma_bits,niter)
         dG_2 = dG_Newton(ya_bits-yd_bits-y0_bits,A_bits,sigma_bits,niter)
         yd_bits = yd_bits + (dG_1 - dG_2)/2.0d0
      else
         yd_bits = yd_bits
      endif

      ! Convert corrected demod value to volts
      yd = bits2volts*yd_bits

      if (avcorrect) write(*,*) 'AV correction not implemented..'

      return

   end subroutine bcorrect

!============================================================

   subroutine init_typeb_correction(run,seg,myid,init_ok,verbose)
      implicit none

      ! Arguments
      integer,                  intent(in)  :: run
      integer,                  intent(in)  :: seg
      integer(i4b),             intent(in)  :: myid
      logical(lgt),             intent(out) :: init_ok

      logical, optional :: verbose
      logical           :: verbose_

      integer sseg,status

      ! Local working variables
      character*128 fullpath
      character*5 crun        ! For CES runids up to 99999
      character*3 cseg        ! For CES segments up to 999
      character*3 ctypeb_vers ! For Type-B versions up to 999
      integer imod, idiode, iglitch, icol, ipatt

      verbose_ = .false.; if(present(verbose)) verbose_ = verbose

      ! Translate the CES numbering from Oslo to Chicago scheme, if needed
      if (typeb_vers == 1) then
         ! This was an imprecise conversion
         sseg = seg - 1
      elseif (typeb_vers == 2) then
         ! Correction was already made to the files, at v2_translated
         sseg = seg
      else
         write(*,*) 'Disallowed or mangled Type-B version', typeb_vers
         return
      endif

      ! Identify Type-B parameter file from run, seg
      ! Convert chars to integers
      write(crun,'(I5)') run
      write(cseg,'(I3)') sseg
      write(ctypeb_vers,'(I3)') typeb_vers

      ! Generate full path to Type-B parameters for this CES
      if (typeb_vers==1) then
         fullpath = trim(typeb_params_dir)//"/v"//trim(adjustl(ctypeb_vers))//"/"//trim(adjustl(crun))//'.'//trim(adjustl(cseg))//'/anomaly.txt'
      elseif (typeb_vers==2) then
         fullpath = trim(typeb_params_dir)//"/v"//trim(adjustl(ctypeb_vers))//"_translated/"//trim(adjustl(crun))//'.'//trim(adjustl(cseg))//'/anomaly.txt'
      else
         write(*,*) 'Mangled Type-B version', typeb_vers
         write(*,*) 'Attempting to continue...'
         fullpath = trim(typeb_params_dir)//"/"//trim(adjustl(crun))//'.'//trim(adjustl(cseg))//'/anomaly.txt'
      endif

! Initialize typeb_params structure

      allocate(typeb_params%preamp_rms_factor(nmodules,ndiodes))
      allocate(typeb_params%glitching_channels(nmodules,ndiodes))
      allocate(typeb_params%detector_rms(nmodules,ndiodes))
      if (typeb_vers == 1) then
         allocate(typeb_params%centers_heights(nknownglitches,4))
      elseif (typeb_vers == 2) then
         allocate(typeb_params%centers_heights22(maxval(nglitch),4,0:size(nglitch)-1))
      endif

      typeb_params%demod_sampling_rate = 0
      typeb_params%n_iterations = 0
      typeb_params%height_factor = 0.0d0
      typeb_params%preamp_rms_factor = 0.0d0
      typeb_params%glitching_channels = 0
      typeb_params%detector_rms = 0.0d0
      if (typeb_vers == 1) then
         typeb_params%centers_heights = 0.0d0
      elseif (typeb_vers == 2) then
         typeb_params%centers_heights22 = 0.0d0
      endif

! Read in Type-B parameters for this CES from file (directly) to
! typeb_params structure

      typeb_init_ok = .false.

      if (debug) write(*,*) 'id = ',myid,': Reading Type-B parameters from', trim(fullpath)

      status = 0

! A bit ugly, but have separated these two reads out
      if (typeb_vers == 1) then
          call read_typeb_params(fullpath,keywords,status)
      elseif (typeb_vers == 2) then
         call read_typeb_params_v22(fullpath,status,verbose=debug)
      else
         write(*,*) 'Unsupported or mangled Type-B parameter version', typeb_vers
         return
      endif

      if (status==0) then
          if(debug) write(*,*) 'id = ',myid,': Type-B parameters read OK'
         typeb_init_ok = .true.
      else
          if(debug) write(*,*) 'id = ',myid,': Unable to read Type-B parameters from file (not there/wrong format)'
          if(debug) write(*,*) 'id = ',myid,': No Type-B correction will be applied for CES ',crun//'.'//trim(adjustl(cseg))
         typeb_init_ok = .false.
      endif

! Now do some housekeeping

      ! Convert detector_rms and centers_heights from bits to volts
      typeb_params%detector_rms    = bits2volts * typeb_params%detector_rms
      if (typeb_vers==1) then
         typeb_params%centers_heights = bits2volts * typeb_params%centers_heights
      elseif (typeb_vers==2) then
         typeb_params%centers_heights22 = bits2volts * typeb_params%centers_heights22
      endif

      ! Apply height factor to glitch amplitudes and their errors
      if (typeb_vers==1) then
         typeb_params%centers_heights(:,3) = typeb_params%height_factor &
           * typeb_params%centers_heights(:,3)
         typeb_params%centers_heights(:,4) = typeb_params%height_factor &
           * typeb_params%centers_heights(:,4)
      elseif (typeb_vers==2) then
         typeb_params%centers_heights22(:,3,:) = typeb_params%height_factor &
              * typeb_params%centers_heights22(:,3,:)
         typeb_params%centers_heights22(:,4,:) = typeb_params%height_factor &
              * typeb_params%centers_heights22(:,4,:)
      endif

! Inspect read-in parameters if required

      if (verbose_) then
         write(*,*) 'status ', status
         write(*,*) 'rate   ', typeb_params%demod_sampling_rate
         write(*,*) 'niter  ', typeb_params%n_iterations
         write(*,*) 'height ', typeb_params%height_factor
         write(*,*) ''

         write(*,*) 'preamp factors'
         do imod = 1, nmodules
            write(*,*) (typeb_params%preamp_rms_factor(imod,idiode),idiode=1,ndiodes)
         enddo

         write(*,*) 'glitching mask'
         do imod = 1, nmodules
            write(*,*) (typeb_params%glitching_channels(imod,idiode),idiode=1,ndiodes)
         enddo

         write(*,*) 'detector rms'
         do imod = 1, nmodules
            write(*,*) (typeb_params%detector_rms(imod,idiode),idiode=1,ndiodes)
         enddo

         if (typeb_vers==1) then
            write(*,*) 'glitch parameters'
            do iglitch = 1, nknownglitches
               write(*,*) (typeb_params%centers_heights(iglitch,icol),icol=1,4)
            enddo
         elseif (typeb_vers==2) then
            write(*,*) 'more glitch (type-2, RQ08Q2) parameters'
            do iglitch = 1,maxval(nglitch)
               write(*,*) (typeb_params%centers_heights22(iglitch,icol,:),icol=1,4)
            enddo
         endif
         do ipatt = 0, size(nglitch)-1
            write(*,*) 'glitch pattern', ipatt
            do iglitch = 1, maxval(nglitch)
               write(*,*) (typeb_params%centers_heights22(iglitch,icol,ipatt),icol=1,4)
            enddo
         enddo
      endif

      status = 0
      init_ok = typeb_init_ok

      return

   end subroutine init_typeb_correction

!============================================================

   subroutine nearest(a,b,d)

      ! Returns location d of value b(d) nearest to a in a 1-D array of numbers b(:)

      implicit none

      ! Arguments
      real(dp), intent(in)  :: a    ! The search value
      real(dp), intent(in)  :: b(:) ! The array to be searched
      integer,  intent(out) :: d    ! The location of b(d) nearest to a

      ! Local working variables
      real(dp) db(size(b))
      integer i,c(1)

      c = 0
      d = 0.0d0

      do i = 1, size(b)
         db(i) = abs(a-b(i))
      end do

      c = minloc(db(:))
      d = c(1)

      return

   end subroutine nearest

!============================================================

   function nearest_func(a,b)

      ! Returns location d of value b(d) nearest to a in a 1-D array of numbers b(:)

      implicit none

      ! Arguments
      real(dp), intent(in)  :: a    ! The search value
      real(dp), intent(in)  :: b(:) ! The array to be searched
      !integer,  intent(out) :: d    ! The location of b(d) nearest to a
      integer nearest_func

      ! Local working variables
      real(dp) db(size(b))
      integer i,c(1)

      c = 0
      nearest_func = -99

      do i = 1, size(b)
         db(i) = abs(a-b(i))
      end do

      c = minloc(db(:))
      nearest_func=c(1)

   end function nearest_func

!============================================================

   subroutine apply_typeb_correction_v3(demod,average,scale,&
        myid,imod,idiode,verbose,multi)

      ! Apply Type-B correction to demod data for one *diode*

      use l1_read_mod

      implicit none

      ! Arguments
      integer(i4b),       intent(in)    :: myid       ! MPI id
      integer(i4b),       intent(in)    :: imod       ! MPI id
      integer(i4b),       intent(in)    :: idiode     ! MPI id
!      integer(i2b),       intent(in) :: scale(:)   ! Data scale
      integer(i2b),       intent(in) :: scale   ! Data scale
      real(dp),           intent(inout) :: demod(:)   ! Demod data   (volts)
      real(dp),           intent(inout) :: average(:) ! Average data (volts)
      !logical(lgt), optional,    intent(in) :: converged     ! Already converged?

      logical, parameter :: decorrect = .true.  ! DO correct demod data
      logical, parameter :: avcorrect = .false. ! DO NOT correct average data

      integer status

      ! Local working variables
      real(dp) A, y0, sigma
      integer nsamples, niter
      integer isamp, iglitch, ipatt

      logical, optional :: verbose
      logical verbose_
      logical(lgt) multi

      if (typeb_vers /= 3) return

      verbose_ = .false.; if(present(verbose)) verbose_ = verbose

      if (.not. typeb_init_ok) then
          if(debug) write(*,*) "id = ", myid,": No Type-B correction applied"
         return
      endif

      !iglitch=0
      ! Ignore an already-converged diode
      !if (converged) then
      !   if (verbose) write(*,*) 'Skipping (already converged) ',idiode
      !   cycle
      !endif

      ! If this is a glitching channel
      !if (typeb_params%glitching_channels(imod,idiode)<0) then
      !   if (verbose) write(*,*) 'Skipping (non-glitching) ',idiode
      !   return
      !endif

      ! Select the glitching pattern
      ipatt=typeb_params%glitching_channels(imod,idiode)
      if (verbose_) write(*,*) 'Using pattern ', ipatt, ' for diode ', imod,idiode

      ! Calculate sigma using, inter alia, scale
      sigma = sqrt(dble(scale) & 
      !sigma = sqrt(dble(scale(1)) & 
           *dble(typeb_params%demod_sampling_rate)) &
           * typeb_params%detector_rms(imod,idiode) &
           * typeb_params%preamp_rms_factor(imod,idiode)

      niter = typeb_params%n_iterations

      write(*,*) 'i',typeb_params%detector_rms(imod,idiode)

      ! Jump out if correcting multiple glitches
      if (multi) then
         call bcorrect_multi(ipatt,imod,idiode,average(:),demod(:),sigma,niter)
         return
      endif

      ! Loop over samples for this channel
      nsamples = size(demod(:))
      iglitch = 0

      over_samples: do isamp = 1, nsamples
         ! Find location iglitch of nearest glitch to this sample value
         call nearest(average(isamp), &
              typeb_params%centers_heights22(:,1,ipatt),iglitch)

         ! Apply Type-B correction for nearest glitch to this sample
         ! NB average, demod are defined on [0,3], so subtract 1
         ! A is already corrected for height factor
         A = typeb_params%centers_heights22(iglitch,3,ipatt)
         y0 = typeb_params%centers_heights22(iglitch,1,ipatt)
         call bcorrect(average(isamp),demod(isamp),&
              A,y0,sigma,niter,decorrect,avcorrect)
      enddo over_samples

      status = 0

      return

   end subroutine apply_typeb_correction_v3

!============================================================

   subroutine bcorrect_multi(ipatt,imod,idiode,avg,demod,sigma,niter)
      implicit none

      integer,  intent(in)    :: ipatt,imod,idiode,niter
      real(dp), intent(in)    :: sigma
      real(dp), intent(in)    :: avg(:)
      real(dp), intent(inout) :: demod(:)

      integer isamp,nsamp

      !call bcorrect(average(idiode-1,isamp),demod(idiode-1,isamp), &
      !     A,y0,sigma,niter,decorrect,avcorrect)

      nsamp=size(demod)

      samples: do isamp = 1,nsamp
       write(*,*) isamp  
      enddo samples

   end subroutine bcorrect_multi

!============================================================

   subroutine apply_typeb_correction(demod,average,scale,&
        imod,myid,verbose,converged)

      ! Apply Type-B correction to demod data for one module

      use l1_read_mod

      implicit none

      ! Arguments
      integer(i4b),              intent(in)    :: myid         ! MPI id
      integer,                   intent(in)    :: imod         ! Module number [1,19]
      integer(i2b),  intent(in)    :: scale(:)     ! Data scale
      real(dp), allocatable,     intent(inout) :: demod(:,:)   ! Demod data   (volts)
      real(dp), allocatable,     intent(inout) :: average(:,:) ! Average data (volts)
      logical(lgt), optional,    intent(in) :: converged(:) ! Which diodes already converged?

      logical, parameter :: decorrect = .true.  ! DO correct demod data
      logical, parameter :: avcorrect = .false. ! DO NOT correct average data

      integer status

      ! Local working variables
      real(dp) A, y0, sigma
      integer nsamples, niter
      integer isamp, idiode, iglitch, ipatt

      logical, optional :: verbose
      logical verbose_

      verbose_ = .false.; if(present(verbose)) verbose_ = verbose

      if (.not. typeb_init_ok) then
          if(debug) write(*,*) "id = ", myid,": No Type-B correction applied"
         return
      endif

      niter = typeb_params%n_iterations
      nsamples = size(demod(0,:))

      iglitch=0
      ! Select Type-B version
      if (typeb_vers==1) then
         ! Loop over diodes
         do idiode = 1, ndiodes
            ! If this is a glitching channel
            if (typeb_params%glitching_channels(imod,idiode)>0) then
               ! Loop over samples for this channel
               do isamp = 1, nsamples
                  ! Calculate sigma for this CES-module-diode
                  ! > sigma_800 = sqrt(6880*100) * detector_rms [in anomaly.txt] *
                  ! > preamp_rms_factor [in anomaly.txt] (eqn 3 in Buder -> Zwart 0907)
                  ! Since these are usually all the same (6880), only set once
                  if (isamp == 1) then
                     sigma = sqrt(dble(scale(isamp))*dble(typeb_params%demod_sampling_rate)) &
                          * typeb_params%detector_rms(imod,idiode) &
                          * typeb_params%preamp_rms_factor(imod,idiode)
                  endif

                  ! Find location iglitch of nearest glitch to this sample value
                  call nearest(average(idiode-1,isamp),typeb_params%centers_heights(:,1),iglitch)
                  iglitch=nearest_func(average(idiode-1,isamp),typeb_params%centers_heights(:,1))
                  ! Apply Type-B correction for nearest glitch to this sample
                  ! NB average, demod are defined on [0,3], so subtract 1
                  A = typeb_params%centers_heights(iglitch,3) ! Already corrected for height factor
                  y0 = typeb_params%centers_heights(iglitch,1)
                  call bcorrect(average(idiode-1,isamp),demod(idiode-1,isamp), &
                       A,y0,sigma,niter,decorrect,avcorrect)

               enddo ! End loop over samples
            endif ! If glitching
         enddo ! End loop over diodes

      elseif (typeb_vers==2 .or. typeb_vers==3) then
         ! Loop over diodes
         over_diodes: do idiode = 1, ndiodes
            ! Ignore any already-converged diodes
            if (converged(idiode) .and. typeb_vers==3) then
               if (verbose) write(*,*) 'Skipping (already converged) ',&
                    imod,idiode
               cycle
            endif
            ! If this is a glitching channel
            if (typeb_params%glitching_channels(imod,idiode)>-1) then
               ! Select the glitching pattern
               ipatt=typeb_params%glitching_channels(imod,idiode)
               if (verbose_) write(*,*) 'Using pattern ', ipatt, ' for diode ', imod,idiode
               ! Loop over samples for this channel
               iglitch = 0
               over_samples: do isamp = 1, nsamples
                  ! Calculate sigma for this CES-module-diode
                  ! > sigma_800 = sqrt(6880*100) * detector_rms [in anomaly.txt] *
                  ! > preamp_rms_factor [in anomaly.txt] (eqn 3 in Buder -> Zwart 0907)
                  ! Since these are usually all the same (6880), only set once
                  if (isamp == 1) then
                     sigma = sqrt(dble(scale(isamp))*dble(typeb_params%demod_sampling_rate)) &
                          * typeb_params%detector_rms(imod,idiode) &
                          * typeb_params%preamp_rms_factor(imod,idiode)
                  endif

                  ! Find location iglitch of nearest glitch to this sample value
                  call nearest(average(idiode-1,isamp),typeb_params%centers_heights22(:,1,ipatt),iglitch)

                  ! Apply Type-B correction for nearest glitch to this sample
                  ! NB average, demod are defined on [0,3], so subtract 1
                  ! A is already corrected for height factor
                  A = typeb_params%centers_heights22(iglitch,3,ipatt)
                  y0 = typeb_params%centers_heights22(iglitch,1,ipatt)
                  call bcorrect(average(idiode-1,isamp),demod(idiode-1,isamp),&
                       A,y0,sigma,niter,decorrect,avcorrect)

               enddo over_samples ! End loop over samples
            endif ! If glitching
         enddo over_diodes ! End loop over diodes

      else
         write(*,*) 'Disallowed or mangled Type-B version', typeb_vers
         return
      endif ! If typeb_vers

      status = 0

      return

    end subroutine apply_typeb_correction

!============================================================

   subroutine deallocate_typeb_param_struct()
      ! Deallocate typeb_params
      implicit none
      integer status

      typeb_params%demod_sampling_rate = 0
      typeb_params%n_iterations = 0
      typeb_params%height_factor = 0.0d0

      if (allocated(typeb_params%preamp_rms_factor))  deallocate(typeb_params%preamp_rms_factor)
      if (allocated(typeb_params%glitching_channels)) deallocate(typeb_params%glitching_channels)
      if (allocated(typeb_params%detector_rms_init))  deallocate(typeb_params%detector_rms_init)
      if (allocated(typeb_params%detector_rms))       deallocate(typeb_params%detector_rms)
      if (allocated(typeb_params%detector_rms_converged)) deallocate(typeb_params%detector_rms_converged)
      if (allocated(typeb_params%centers_heights))    deallocate(typeb_params%centers_heights)
      if (allocated(typeb_params%centers_heights22))  deallocate(typeb_params%centers_heights22)
      status = 0

      write(*,*) 'Deallocated Type-B structure..'

      return
   end subroutine deallocate_typeb_param_struct

!============================================================


 
end module typeb_correction_mod

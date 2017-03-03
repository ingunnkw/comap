program map_editor
  use healpix_types
  use map_editor_utils
  use map_editor_simple_ops_mod
  use map_editor_complex_ops_mod
  implicit none

  ! *******************************************************
  ! *                                                     *
  ! *          Utility for editing Healpix maps           *
  ! *                                                     *
  ! *  written by Hans Kristian Eriksen, November, 2004   *
  ! *                                                     *
  ! *         Copyright 2004. All rights reserved.        *
  ! *                                                     *
  ! *******************************************************

  !
  ! ChangeLog:
  ! 
  ! November 2, 2004  -- First version
  !

  integer(i4b)       :: iargc, lmin, lmax, lcut, nside_out, unit, i, l, m, n, seed
  integer(i4b)       :: nside, ordering, nmaps, component, s_max, ncol
  character(len=3)   :: suffix
  character(len=128) :: mapname_in1, mapname_in2, infofile, mapname_out, maskfile, rmsfile
  character(len=128) :: beamfile_in, beamfile_out, option
  character(len=128) :: string_real, string_int, operation, beaminfo, covartype, map2mask_file, outprefix
  real(dp)           :: value, sigma_0, rms, cl0, r_fill
  real(dp)           :: fwhm_in, fwhm_out, md(4), fact, f(3)

  real(dp),     pointer, dimension(:,:) :: map, map2, resmap
  logical(lgt), pointer, dimension(:)   :: mask
  character(len=80), dimension(180)  :: header

  if (iargc() == 0) then
     write(*,*) ' '
     write(*,*) 'Usage: map_editor [operation] [arg1] [arg2] ...'
     write(*,*) ' '
     write(*,*) 'where [operation] = {help, scale, add_offset, log, ln, exp, abs, inv, sqrt, '
     write(*,*) '                     add, subtract, multiply, divide, half_sum, half_diff, '
     write(*,*) '                     weighted_sum, smooth, add_gaussian_noise, subtract_low_l, '
     write(*,*) '                     rms2mask, QU2P, shift_columns, ptsrc, compute_mean_stddev,'
     write(*,*) '                     median_filter_misspix, scale_TQU}'
     write(*,*) ''
     write(*,*) '   For more information on any given operation, execute the '
     write(*,*) '   program with options "help + operation", e.g.:'
     write(*,*) ''
     write(*,*) '       map_editor help smooth'
     write(*,*) ''
     stop
  end if

  unit = 25

  call getarg(1,operation)
  n = len(trim(operation))
  suffix = operation(n-2:n)

  if (trim(operation) == 'scale' .or. trim(operation) == 'add_offset' .or. &
       & trim(operation) == 'log' .or. trim(operation) == 'ln' .or. &
       & trim(operation) == 'exp' .or. trim(operation) == 'sqrt' .or. &
       & trim(operation) == 'abs' .or. trim(operation) == 'inv' .or. &
       & trim(operation) == 'hitcount2rms' .or. trim(operation) == 'rms2mask' .or. &
       & trim(operation) == 'amp2mask' .or. trim(operation) == 'asinh' .or. &
       & trim(operation) == 't2a' .or. trim(operation) == 'a2t' .or. &
       & trim(operation) == 'max_scalar' .or. trim(operation) == 'min_scalar' .or. &
       & trim(operation) == 'hitcount2mask' .or. trim(operation) == 'QU2P' .or. &
       & trim(operation) == 'missing2mask') then

     call getarg(2,mapname_in1)
     call getarg(3,mapname_out)
     if (iargc() == 4) then
        call getarg(4,string_real)
        read(string_real,*) value
     else
        value = 0.
     end if

     call initialize_single_map(mapname_in1, nside, ordering, nmaps, header, map)

     call operate_on_single_map(map, operation, value)
     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'ring2nest' .or. trim(operation) == 'nest2ring') then

     call getarg(2,mapname_in1)
     call getarg(3,mapname_out)

     call initialize_single_map(mapname_in1, nside, ordering, nmaps, header, map)

     if (trim(operation) == 'ring2nest') then

        if (ordering == 2) then
           write(*,*) 'Error: Map is already in NESTED format'
           stop
        end if

        do i = 1, nmaps
           call convert_ring2nest(nside, map(:,i))
        end do

        ordering = 2

     else

        if (ordering == 1) then
           write(*,*) 'Error: Map is already in RING format'
           stop
        end if

        do i = 1, nmaps
           call convert_nest2ring(nside, map(:,i))
        end do

        ordering = 1

     end if

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')


  else if (trim(operation) == 'scale_TQU') then

     if (iargc() /= 6) then
        write(*,*) 'Usage: map_editor scale_TQU [infile] [outfile] [T scale] [Q scale] [U scale]'
        stop
     end if

     call getarg(2,mapname_in1)
     call getarg(3,mapname_out)
     call getarg(4,string_real)
     read(string_real,*) f(1)
     call getarg(5,string_real)
     read(string_real,*) f(2)
     call getarg(6,string_real)
     read(string_real,*) f(3)

     call initialize_single_map(mapname_in1, nside, ordering, nmaps, header, map)
     do i = 1, 3
        map(:,i) = map(:,i) * f(i)
     end do
     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'apply_mask') then

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)
     call getarg(4,mapname_out)
     if (iargc() == 5) then
        call getarg(5,string_real)
        read(string_real,*) fact
     end if
     call initialize_single_map(mapname_in1, nside, ordering, nmaps, header, map)
     if (iargc() == 4) then
        call apply_mask(maskfile, nside, ordering, map)
     else if (iargc() == 5) then
        call apply_mask(maskfile, nside, ordering, map, fact)
     else
        write(*,*) "Usage: map_editor apply_mask [input map] [mask] [output map]"
     end if
     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'compute_mean_stddev') then

     call getarg(2,outprefix)
     call getarg(3,maskfile)

     call output_mean_and_stddev(nside, ordering, nmaps, header, map, map2)
     call write_result_map(trim(outprefix)//'_mean.fits', nside, ordering, header, map, suffix=='_dp')
     call write_result_map(trim(outprefix)//'_stddev.fits', nside, ordering, header, map2, suffix=='_dp')

  else if (trim(operation) == 'add' .or. trim(operation) == 'subtract' .or. &
       & trim(operation) == 'multiply' .or. trim(operation) == 'divide' .or. &
       & trim(operation) == 'max' .or. trim(operation) == 'min' .or. &
       & trim(operation) == 'half_sum' .or. trim(operation) == 'half_diff' .or. &
       & trim(operation) == 'splice') then

     ! Simple two map operations
     call getarg(2,mapname_in1)
     call getarg(3,mapname_in2)
     call getarg(4,mapname_out)

     call initialize_single_map(mapname_in1, nside, ordering, nmaps, header, map)
     call initialize_second_map(mapname_in2, nside, ordering, nmaps, map2)

     call operate_on_two_maps(map, map2, operation)
     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'weighted_sum') then

     call getarg(2,infofile)
     call getarg(3,mapname_out)
     call compute_weighted_sum(infofile, nside, ordering, nmaps, resmap, header)
     call write_result_map(mapname_out, nside, ordering, header, resmap, suffix=='_dp')

  else if (trim(operation) == 'ptsrc') then

     call output_pointsource

  else if (trim(operation) == 'smooth') then

     call getarg(2,beaminfo)
     call getarg(3,mapname_in1)
     call getarg(4,string_int)
     read(string_int,*) lmin
     call getarg(5,string_int)
     read(string_int,*) lmax
     call getarg(6,string_int)
     read(string_int,*) nside

     if (iargc() == 10) then
        call getarg(10,string_int)
        read(string_int,*) r_fill
     else
        r_fill = -1.d0
     end if

     if (trim(beaminfo) == 'f2f') then
        ! Input beam from file, output beam from file
        call getarg(7,beamfile_in)
        call getarg(8,beamfile_out)
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .false., &
             & beamfile_in=beamfile_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'f2g') then
        ! Input beam from file, output beam from gaussian
        call getarg(7,beamfile_in)
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .false., &
             & beamfile_in=beamfile_in, fwhm_out=fwhm_out)
     else if (trim(beaminfo) == 'g2f') then
        ! Input beam from gaussian, output beam from file
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,beamfile_out)
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .false., &
             & fwhm_in=fwhm_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'g2g') then
        ! Input beam from gaussian, output beam from gaussian
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .false., &
             & fwhm_in=fwhm_in, fwhm_out=fwhm_out)
     else if (trim(beaminfo) == 'f2f_EB') then
        ! Input beam from file, output beam from file
        call getarg(7,beamfile_in)
        call getarg(8,beamfile_out)
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .true., &
             & beamfile_in=beamfile_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'f2g_EB') then
        ! Input beam from file, output beam from gaussian
        call getarg(7,beamfile_in)
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .true., &
             & beamfile_in=beamfile_in, fwhm_out=fwhm_out)
     else if (trim(beaminfo) == 'g2f_EB') then
        ! Input beam from gaussian, output beam from file
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,beamfile_out)
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .true., &
             & fwhm_in=fwhm_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'g2g_EB') then
        ! Input beam from gaussian, output beam from gaussian
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .true., &
             & fwhm_in=fwhm_in, fwhm_out=fwhm_out)
     else 
        write(*,*) 'Invalid beam option. Exiting.'
        stop
     end if

     call getarg(9,mapname_out)

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'shift_columns') then

     call getarg(2,mapname_in1)
     call getarg(3,mapname_out)
     call getarg(4,string_int)
     read(string_int,*) ncol

     call shift_columns(mapname_in1, ncol, nside, ordering, header, map)

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'add_gaussian_noise') then

     call getarg(2,mapname_in1)
     call getarg(3,covartype)
     call getarg(4,rmsfile)
     call getarg(5,string_int)
     read(string_int,*) seed
     call getarg(6,mapname_out)
     
     if (trim(covartype) == 'sqrt_N') then

        call getarg(7,map2mask_file)
        
        ! Assume that second file contains sqrt_N of Gaussian noise
        call add_gaussian_noise_sqrt_N(mapname_in1, rmsfile, map2mask_file, seed, nside, ordering, nmaps, map, header)

     else

        if (iargc() == 6) then

           ! Assume that second file contains RMS of Gaussian noise
           call add_gaussian_noise(mapname_in1, rmsfile, seed, nside, ordering, nmaps, map, header)

        else if (iargc() == 7) then
           ! Assume that second file contains Nobs of Gaussian noise
           call getarg(7,string_real)
           read(string_real,*) sigma_0
           call add_gaussian_noise(mapname_in1, rmsfile, seed, nside, ordering, nmaps, map, header, sigma_0)
        end if

     end if

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'subtract_mono_dipole') then

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)
     call getarg(4,mapname_out)
     if (iargc() == 8) then
        call getarg(5,string_real)
        read(string_real,*) md(1)
        call getarg(6,string_real)
        read(string_real,*) md(2)
        call getarg(7,string_real)
        read(string_real,*) md(3)
        call getarg(8,string_real)
        read(string_real,*) md(4)
        call subtract_mono_dipole(mapname_in1, maskfile, nside, ordering, map, header, md)
     else
        call subtract_mono_dipole(mapname_in1, maskfile, nside, ordering, map, header)
     end if

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'subtract_mono_dipole_highl') then

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)
     call getarg(4,string_int)
     read(string_int,*) lmax
     call getarg(5,string_int)
     read(string_int,*) lcut
     call getarg(6,mapname_out)
     call subtract_mono_dipole_highl(mapname_in1, maskfile, lmax, lcut, nside, ordering, map, header)

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'partrans') then

     call getarg(2,mapname_in1)
     call getarg(3,string_int)
     read(string_int,*) nside_out
     call getarg(4,mapname_out)
     call initialize_single_map(mapname_in1, nside, ordering, nmaps, header, map)
     write(*,*) nside, ordering, nside_out
     call qu_transport_map(nside_out, map)
     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'make_co_region_map') then

     call getarg(2,mapname_in1)
     call getarg(3,mapname_out)

     call make_co_region_map(mapname_in1, nside, ordering, map)

     call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=size(map,2)==3)
     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'extract_multipole_range') then

     call getarg(2,mapname_in1)
     call getarg(3,string_int)
     read(string_int,*) lmin
     call getarg(4,string_int)
     read(string_int,*) lmax
     call getarg(5,mapname_out)

     call extract_multipole_range(mapname_in1, lmin, lmax, nside, ordering, map, header)

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'ud_grade') then

     if (iargc() /= 4) then
        write(*,*) 'Usage: map_editor ud_grade [input map] [nside_out] [output map]'
        stop
     end if
     call getarg(2,mapname_in1)
     call getarg(3,string_int)
     read(string_int,*) nside
     call getarg(4,mapname_out)
     call ud_grade_map_editor(mapname_in1, nside, mapname_out, suffix=='_dp')

  else if (trim(operation) == 'median_filter_source_holes') then

     call median_filter_holes

  else if (trim(operation) == 'median_filter_misspix') then

     call median_filter_misspix

  else if (trim(operation) == 'print_map_to_ascii') then

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)
     call getarg(4,mapname_out)

     call print_map_to_ascii(mapname_in1, maskfile, mapname_out)

  else if (trim(operation) == 'print_two_maps_to_ascii') then

     call getarg(2,mapname_in1)
     call getarg(3,mapname_in2)
     call getarg(4,maskfile)
     call getarg(5,mapname_out)

     call print_two_maps_to_ascii(mapname_in1, mapname_in2, maskfile, mapname_out)

  else if (trim(operation) == 'print_isolatitude') then

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)
     call getarg(4,mapname_out)

     call print_isolatitude(mapname_in1, maskfile, mapname_out)

  else if (trim(operation) == 'print_isolatitude_var') then

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)
     call getarg(4,mapname_out)

     call print_isolatitude_var(mapname_in1, maskfile, mapname_out)

  else if (trim(operation) == 'fit_ideal_dust') then

     call fit_ideal_dust

  else if (trim(operation) == 'help') then

     call getarg(2,option)        
     call print_help(option)

  else if (trim(operation) == 'print_maximum') then

     call getarg(2,mapname_in1)
     call getarg(3,string_int)
     read(string_int,*) component
     call getarg(4,string_int)
     read(string_int,*) fwhm_in

     call print_maximum(mapname_in1, component, fwhm_in)

  else if (trim(operation) == 'print_stats') then

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)

     call print_stats(mapname_in1, maskfile)

  else if (trim(operation) == 'print_stats_col') then

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)

     call print_stats_col(mapname_in1, maskfile)

  else if (trim(operation) == 'summarize_detector_angles') then

     call getarg(2,mapname_in1)
     call getarg(3,string_int)
     read(string_int,*) s_max
     call getarg(4,string_int)
     read(string_int,*) nside
     call getarg(5,mapname_out)

     call summarize_detector_angles(mapname_in1, s_max, nside, mapname_out)

  else if (trim(operation) == 'compute_spectral_index_map') then

     call compute_index_map

  else if (trim(operation) == 'maskcount') then

     call getarg(2,mapname_in1)
     call maskcount(mapname_in1)

  else if (trim(operation) == 'badcount') then

     call getarg(2,mapname_in1)
     call maskcount(mapname_in1)

  else if (trim(operation) == 'merge_maps') then

     call merge_maps(suffix=='_dp')

  else if (trim(operation) == 'make_ptsrc_map') then

     call make_ptsrc_map

  else if (trim(operation) == 'print_scaled_gal_avg') then

     call print_scaled_gal_avg

  else if (trim(operation) == 'fit_line_to_ASCII_data') then
     
     call fit_line

  else if (trim(operation) == 'convert_beam') then

     call getarg(2,string_real)
     read(string_real,*) fwhm_in
     call getarg(3,string_int)
     read(string_int,*) lmax
     call getarg(4,string_int)
     read(string_int,*) nmaps
     call getarg(5,beamfile_out)

     if (iargc() == 6) then
        call getarg(6,beamfile_in)
        call convert_beam(fwhm_in, lmax, nmaps, beamfile_out, beamfile_in)
     else
        call convert_beam(fwhm_in, lmax, nmaps, beamfile_out)
     end if

  else
     write(*,*) 'Unknown operation. Exiting.'
     stop 
  end if


  ! Clean up arrays and exit 
  ! Compaq's F90 compiler doesn't seem to allow allocated() calls on pointers
!  if (allocated(map))    deallocate(map)
!  if (allocated(map2))   deallocate(map2)
!  if (allocated(resmap)) deallocate(resmap)
!  if (allocated(mask))   deallocate(mask)


contains

  subroutine print_help(option)
    implicit none

    character(len=*), intent(in) :: option

    if (trim(option) == 'scale' .or. trim(option) == 'add_offet' .or. &
         & trim(option) == 'log' .or. trim(option) == 'ln' .or. &
         & trim(option) == 'exp' .or. &
         & trim(option) == 'abs') then
       
       write(*,*) ''
       write(*,*) '   The following operations act on a single map, and require '
       write(*,*) '   only an input filename and an output filename:'
       write(*,*) ''
       write(*,*) '        log, ln, exp, abs '
       write(*,*) ''
       write(*,*) '   The following operations require an additional floating'
       write(*,*) '   point argument:'
       write(*,*) ''
       write(*,*) '        scale, add_offset'
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '        map_editor log inmap.fits outmap.fits'
       write(*,*) '        map_editor scale inmap.fits outmap.fits 2.0'
       write(*,*) ''

  else if (trim(option) == 'add' .or. trim(option) == 'subtract' .or. &
       & trim(option) == 'multiply' .or. trim(option) == 'divide' .or. &
       & trim(option) == 'half_sum' .or. trim(option) == 'half_diff') then
    
       write(*,*) ''
       write(*,*) '   The following operations act on two maps, and require '
       write(*,*) '   two input filenames and one output filename:'
       write(*,*) ''
       write(*,*) '        add, subtract, multiply, divide, half_sum, half_diff '
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '        map_editor subtract inmap1.fits inmap2.fits outmap.fits'
       write(*,*) ''

    else if (trim(option) == 'smooth') then

       write(*,*) ''
       write(*,*) '   The smoothing operation reads a map, computes its spherical '
       write(*,*) '   harmonics transform, deconvolve a given beam (and pixel window),'
       write(*,*) '   convolve with a new given beam (and pixel window), and finally'
       write(*,*) '   outputs resulting the inverse spherical harmonics transform.'
       write(*,*) ''
       write(*,*) '   Both the input and output beams can be given on two formats,'
       write(*,*) '   either in the form of the FWHM of a Gaussian beam, or as a'
       write(*,*) '   FITS ascii table.'
       write(*,*) ''
       write(*,*) '   The beam id must be one of the following strings: '
       write(*,*) ''
       write(*,*) '        f2f -> Input and output beams are read from files'
       write(*,*) '        f2g -> Input beam is read from file, output beam is Gaussian'
       write(*,*) '        g2f -> Input beam is Gaussian, output beam is read from file'
       write(*,*) '        g2g -> Input and output beams are Gaussian'
       write(*,*) ''
       write(*,*) '   If a file is requested, then [input/output beam] must be a filename.'
       write(*,*) '   If a Gaussian is requested, then it must be a floating point number,'
       write(*,*) '   giving the FWHM of the Gaussian beam in arcmin.'
       write(*,*) ''
       write(*,*) '   The variables lmin, lmax and Nside may formally be choosen freely '
       write(*,*) '   (independent of the input map), but a good rule of thumb is'
       write(*,*) '   lmax < 3*Nside_in, and preferrably even lmax < 2*Nside_in.'
       write(*,*) ''
       write(*,*) '   The command line format is as follows:'
       write(*,*) ''
       write(*,*) '        map_editor smooth [beam id] [input filename] [lmin] [lmax]'
       write(*,*) '               [nside_out] [input beam] [output beam]'
       write(*,*) '               [output filename] [radius fill in pixels; optional]'
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '        map_editor smooth f2g map.fits 2 1024 512 beam.fits 60. smoothed.fits'
       write(*,*) ''

    else if (trim(option) == 'add_gaussian_noise') then

       write(*,*) ''
       write(*,*) '   Gaussian noise can be added to an existing map by calling the program'
       write(*,*) '   with the "add_gaussian_noise" option. This mode of operation assumes'
       write(*,*) '   an input filename, an RMS/Nobs filename, a seed and an output filename.'
       write(*,*) '   If only four arguments are provided, the second file is assumed to '
       write(*,*) '   contain the noise RMS. If five is provided (the fifth being '
       write(*,*) '   the standard deviation of a single detector hit, sigma_0), it is'
       write(*,*) '   assumed to contain the number of observations for each pixel, Nobs.'
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '        map_editor add_gaussian_noise map.fits type rms.fits -128341 noisymap.fits '
       write(*,*) '        map_editor add_gaussian_noise map.fits type nobs.fits -128341 noisymap.fits 100.'
       write(*,*) ''

    else if (trim(option) == 'convert_beam') then

       write(*,*) ''
       write(*,*) '   Convert beamfiles between formats'
       write(*,*) ''
       write(*,*) '         map_editor convert_beam [fwhm] [lmax] [nmaps] [outfile] [infile]'
       write(*,*) ''
       write(*,*) '   Note that infile is optional.'
       write(*,*) ''

    else if (trim(option) == 'help') then

       write(*,*) ''
       write(*,*) '   Execute the program as follows to get more information about an operation:'
       write(*,*) ''
       write(*,*) '         map_editor help [operation]'
       write(*,*) ''
       write(*,*) '   For a list of all available operation, run the program without any'
       write(*,*) '   arguments.'
       write(*,*) ''

    else

       write(*,*) ''
       write(*,*) '   Unknown option.'
       write(*,*) ''

    end if

  end subroutine print_help

end program map_editor

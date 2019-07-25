# comap

TODO list:

Highest priority:
- Automate the pipeline from level 1 to P(k) (fix all the bugs that make it crash at random times) (Marie/Håvard)

l2gen:
- Implement new calibration procedure (Håvard)
- Include weather data from yr in l2 format (Maren)
- Add option to do polyfilter on whole feed (Håvard)
- Test David's calibration suggestion (Håvard)
- Make code to find/deal with spikes (someone) (store info for diagnostics)

- Simplify the masking procedure and implement the "reason" array properly (Håvard)
- Try tricks with pca-filter before polyfilter (Håvard)

tod2comap: 
- Fix rms of the "all feeds"-map (Marie)
- Add possibility to make map from simulations present in l2 (Maren and Marie) (probably not final soulution)

Simulations:
- Discuss how to deal with simulated tods in a way where we will not have to store tods or maps in memory (everyone)

Scan_detect:
- Implement jack-knife data in runlist (Håvard)

Diagnostics:
- Make code for combining correlation plots from multiple scans (Håvard and/or Maren)
    - make it possible to compare data with different flags to find the corresponding features in the corr-plot (Håvard and/or Maren)
- Make tsys-diagnostic plot (Håvard)
- make diagnotics, maps, powerspectra available to everyone automatically (Håvard + ovro/Tim)

map2ps:
- implement paralell x perpendicular power spectrum (Håvard)

Other:
- Make "dir" instead of protodir for a clean start with the new science data?

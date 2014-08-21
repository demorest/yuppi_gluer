yuppi_gluer
============

This code takes SIGPROC filterbank files that record a single sub-band of a single
scan of an observation and stitches them together in frequency and time to create
an observation-long data file in PSRFITS format.  Gaps in the data (in frequency or time)
will be filled with zeros in the output file.

Usage
------ 

<pre>
$ ./fil2fits 
   -o outfile [-g outlenGB] [-src src_name] [-raj raj] [-dej dej] [--] infiles ...
      Converts sub-band SIGPROC .fil data into full-band PSRFITS
         -o: Name of the output psrfits file
             1 char* value
         -g: Approx length in GB of output data files
             1 int value between 1 and 1000
             default: `10'
       -src: Source Name
             1 char* value
             default: `BLANK'
       -raj: Source RA  (J2000) in hhmmss.ss
             1 double value between 0 and 240000.00
             default: `0.00'
       -dej: Source Dec (J2000) in ddmmss.ss
             1 double value between -900000.00 and 900000.00
             default: `0.00'
    infiles: Data files
             1...16384 values
  version: 21Apr14
</pre>



yuppi_gluer
============

This code takes SIGPROC filterbank files that record a single sub-band of a single
scan of an observation and stitches them together in frequency and time to create
an observation-long data file in PSRFITS format.  Gaps in the data (in frequency or time)
will be filled with zeros in the output file.

Install
---------

Everything is included, so you should just be able to download and run `make`

The code does require the environment variable `YUPPI_DIR` to be set to the 
directory containing the PSRFITS template (in this case, `yuppi_PSRFITS_v4.3_search_template.txt`).


Usage
------ 

Running the executable `fil2fits` without any arguments gives all the options:

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

The only required arguments here are the base name for the output PSRFITS file
and the input SIGPROC filterbank files that you want to combine.  

The optional arguments `src`, `raj`, and `dej` will by default be filled in 
with the appropriate values from the input filterbank file header string.  They
are provided here only for convenience in replace missing or incorrect values
in the filterbank headers.

The `-g` option sets the maximum size of individual output files.  If there is more
data beyond this value, a new PSRFITS file will be created.  

The input files can be in any order and wildcards may be used.


A Few Warnings
----------------

1. Currently, gaps between the filterbank data files (in frequency or time) will
be filled with zeros regardless of how large the gap is.  This will be fixed (see below),
but for now be careful of accidentally included scans from different days, as the full
time between them will be filled with zeros.

2. In order to determine which sub-bands are in the same scan, we bin the start times
of the filterbank files in 10 second chunks to assign scan numbers.  We do this because
each sub-band file in each scan does not necessarily have exactly the same start time. 
This means that scans shorter than 10 seconds will not be handled properly.

3. Depending on the number and size of the sub-band files, combining them into one FITS
file can be a bit slow, but is read/write limited.  Possibly can improve this by chunking 
the reads and writes a bit better.


Still To Do
---------------

1. Set maximum gap size for frequency and time.

2. Make the bin time for determining scan number an optional parameter.

3. Look for ways to improve processing speed.
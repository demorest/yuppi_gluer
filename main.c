#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "presto2.h"
#include "vectors2.h"
#include "sigproc_fb2.h"
#include "main_cmd.h"
#include "psrfits.h"

#define DEBUG 0

int comp_fb(const void *a, const void *b);
int comp_bin_index(const void *a, const void *b);
int comp_bin_sbin(const void *a, const void *b);
int comp_bin_scan(const void *a, const void *b);
int get_sign(float x);
int check_fbs(sigprocfb *fb, int numfiles);
void init_sigproc(sigprocfb *fb);
void coord_double2str(double x, char *outstr, int is_ra);

typedef struct BIN_INDEX {
  double fch1;  /* Channel 1 Frequency */
  int sbin;     /* Scan bin number */
  int scan;     /* Scan number */
  int index;    /* Original Index */
} bin_index;

int main(int argc, char *argv[])
{
  FILE **infiles;
  unsigned char **rawdata;
  int headerlen;
  int numfiles = 1, nchans, nbits;
  int min_mjd, min_sec, imjd, smjd;
  int tmp_mjd, tmp_sec;
  double tmp_offs;
  int spec_per_row, bytes_per_subband, nchans_subband;
  int rownum, breakout, vals_read, specnum, ichan, datidx;
  int fch1idx, fsign,  datidx_start, rawidx_start;
  double max_freq, min_freq, tmp_min, tmp_max, df, dt, offs;
  char source_name[80], str_ra[20], str_dec[20];
  char basename[80];
  int ii, status, tmpval;
  sigprocfb *fbs;
  struct psrfits pf;
  int *idx;
  Cmdline *cmd;

  /* Call usage() if we have no command line arguments */
  if (argc == 1) {
    Program = argv[0];
    usage();
    exit(0);
  }

  /* Parse the command line using the excellent program Clig */
  cmd = parseCmdline(argc, argv);

  strcpy(basename, cmd->outfile);
  numfiles = cmd->argc;
  infiles = (FILE **) malloc(numfiles * sizeof(FILE *));
  idx = gen_ivect(numfiles);

  fbs = (sigprocfb *) malloc(numfiles * sizeof(sigprocfb));

  /* Initialize filterbank structures and fill them with header info */
  for(ii=0; ii<numfiles; ii++)
    {
      init_sigproc(&fbs[ii]);
      //printf("Working on file %s\n", cmd->argv[ii]);
      infiles[ii] = chkfopen(cmd->argv[ii], "rb");
      headerlen = read_filterbank_header(&fbs[ii], infiles[ii]);
      if ( !strcmp("ERROR", fbs[ii].inpfile) )
	strcpy(fbs[ii].inpfile, cmd->argv[ii]);
      //print_filterbank_header(&fbs[ii]);
    }

  /* Check headers to make sure files are consistent */
  if((tmpval = check_fbs(fbs, numfiles)) != 0)
    return 1;

  /* If necessary, make command line changes to header info */
  for(ii=0; ii<numfiles; ii++)
    {
      if((tmpval = strcmp(cmd->src_name, "BLANK")) == 0)
	strcpy(fbs[ii].source_name, cmd->src_name);
      if(fabs(cmd->raj) > 0.0)
	fbs[ii].src_raj = cmd->raj;
      if(fabs(cmd->dej) > 0.0)
	fbs[ii].src_dej = cmd->dej;
    }

  /* Calculate time vals and find minima */
  /* Also find min/max frequencies */
  min_mjd = -1;
  min_sec = -1;
  tmp_mjd = 0;
  tmp_sec = 0;
  tmp_offs = 0.0;
  min_freq = -1.0;
  max_freq = -1.0;
  tmp_min = 0.0;
  tmp_max = 0.0;
  for(ii=0; ii<numfiles; ii++)
    {
      tmp_mjd = (int) fbs[ii].tstart;
      tmp_sec = (int) ((fbs[ii].tstart - tmp_mjd) * 24 * 3600 + 0.5);
      tmp_offs = (fbs[ii].tstart - tmp_mjd) * 24 * 3600 - tmp_sec;
      if( tmp_mjd < min_mjd || (min_mjd < 0 && min_sec < 0)){
	min_mjd = tmp_mjd;
	min_sec = tmp_sec;
      }
      if( tmp_sec < min_sec && tmp_mjd == min_mjd)
	min_sec = tmp_sec;
      fbs[ii].imjd = tmp_mjd;
      fbs[ii].smjd = tmp_sec;
      fbs[ii].offs = tmp_offs;
      fbs[ii].offs_samp = (int) roundf(tmp_offs / fbs[ii].tsamp); 

      tmp_min = (fbs[ii].foff > 0) ? fbs[ii].fch1 : \
	fbs[ii].fch1 + (fbs[ii].nchans - 1) * fbs[ii].foff;

      tmp_max = (fbs[ii].foff < 0) ? fbs[ii].fch1 : \
        fbs[ii].fch1 + (fbs[ii].nchans - 1) * fbs[ii].foff;

      if(tmp_min < min_freq || min_freq < 0)
	min_freq = tmp_min;
      if(tmp_max > max_freq || max_freq < 0)
	max_freq = tmp_max;
    }
  
  /* Get the channel width from one of the headers */
  df = fabs(fbs[0].foff);
  nchans = (int) (roundf( (max_freq - min_freq) / df ) + 1);

#if (DEBUG)
  printf("nchans = %d of width %.1f making a total bw of %.1f\n", nchans, df, nchans * df);
  printf("total bandwidth = %.1f - %.1f + %.1f = %.1f\n", max_freq, min_freq, df, df + max_freq - min_freq);
#endif
  
  printf("Max freq = %f\n", max_freq);
  printf("Min freq = %f\n", min_freq);
  
  /* Determine Scan Number */
  int last_bin, nscans, bin_dt;
  int *sidx, *sblocks, *refsec;
  bin_index *bidx;
  bidx = (bin_index *) malloc(numfiles * sizeof(bin_index));
  
  /* Calculate the seconds from start of obs and the bin values */
  bin_dt = 10;  // NEEDS TO BE AN INPUT PARAMETER!
  for(ii=0; ii<numfiles; ii++)
    {
      fbs[ii].obs_sec = (fbs[ii].imjd - min_mjd) * 24 * 3600 + (fbs[ii].smjd - min_sec);
      bidx[ii].fch1 = fbs[ii].fch1;
      bidx[ii].sbin = fbs[ii].obs_sec / bin_dt;
      bidx[ii].scan = 0;
      bidx[ii].index = ii;
    }
  
  /* Sort the bin_index structure by sbin */
  qsort(bidx, numfiles, sizeof(bin_index), comp_bin_sbin);
  
  /* Find the scan numbers */
  last_bin = -1;
  nscans = -1;  
  for(ii=0; ii<numfiles; ii++)
    {
      if(bidx[ii].sbin > last_bin){
	nscans++;
	last_bin = bidx[ii].sbin;
      }
      
      bidx[ii].scan = nscans;
      fbs[bidx[ii].index].scan = nscans;
    }
  nscans += 1;

  /* Sort by scan, then fch1 */
  qsort(bidx, numfiles, sizeof(bin_index), comp_bin_scan);
  
  /* Find number of blocks per scan and where they start */
  /* Also find the smallest obs_sec per scan */
  sidx = (int *) (malloc)(nscans * sizeof(int));
  sblocks = (int *) (malloc)(nscans * sizeof(int));
  refsec = (int *) (malloc)(nscans * sizeof(int));
  for(ii=0; ii<nscans; ii++){
    sidx[ii] = sblocks[ii] = 0;
    refsec[ii] = -1;
  }
  int maxblocks = 0;
  int iscan = 0;  
  for(ii=0; ii<numfiles; ii++)
    {
      iscan = bidx[ii].scan;
      if(sblocks[iscan] == 0)
        sidx[iscan] = ii;
      sblocks[iscan] += 1;
      if(sblocks[iscan] > maxblocks)
	maxblocks = sblocks[iscan];
      if(fbs[bidx[ii].index].obs_sec < refsec[iscan] || refsec[iscan] < 0)
	refsec[iscan] = fbs[bidx[ii].index].obs_sec;
    }

  /* Find offset in seconds between start of scan and start of block */
  for(ii=0; ii<numfiles; ii++)
    {
      fbs[ii].scan_sec = fbs[ii].obs_sec - refsec[fbs[ii].scan];
    }

#if (DEBUG)
  printf("\nThere are %d scans\n", nscans);
  for(ii=0; ii<nscans; ii++)
    printf("Scan %d has %d blocks\n", ii, sblocks[ii]);

  printf("\nFilename, Scan number:\n");
  for(ii=0; ii<numfiles; ii++)
    printf("File: %s,  Scan: %d\n", fbs[ii].inpfile, fbs[ii].scan);

  printf("\n Scan, fch1:\n");
  for(ii=0; ii<numfiles; ii++)
    printf("Scan: %d, fch1: %.2f, obs_sec: %d, scan_sec: %d, offs_samp: %d, offs = %.2f\n", fbs[bidx[ii].index].scan, fbs[bidx[ii].index].fch1, fbs[bidx[ii].index].obs_sec, fbs[bidx[ii].index].scan_sec, fbs[bidx[ii].index].offs_samp, fbs[bidx[ii].index].offs * 1000000);
#endif
  
  /* Get sorted indices (by scan, fch1) */
  for(ii=0; ii<numfiles; ii++){
    idx[ii] = bidx[ii].index;
#if (DEBUG)
    printf("ii = %d, idx = %d\n", ii, idx[ii]);
#endif
  }

  /* Get header parameters and initialize psrfits struct */
  nbits = fbs[idx[0]].nbits;
  dt = fbs[idx[0]].tsamp;
  strcpy(source_name, fbs[idx[0]].source_name);
  coord_double2str(fbs[idx[0]].src_raj, str_ra, 1);   // is_ra = 1
  coord_double2str(fbs[idx[0]].src_dej, str_dec, 0);  // is_ra = 0

  fill_psrfits_struct(nbits, &pf, dt, source_name, nchans, df, \
		      min_mjd, min_sec, 0.0, str_ra, str_dec, min_freq, basename);

  spec_per_row = pf.hdr.nsblk;
  pf.multifile = 1;
  pf.rows_per_file = (int) ((cmd->outlenGB * 1073741824.0)
			    / (double) pf.sub.bytes_per_subint);
  nchans_subband = fbs[idx[0]].nchans;
  bytes_per_subband = spec_per_row * nchans_subband;

  status = psrfits_create(&pf);
  if (status) {
    printf("\nError (%d) creating PSRFITS...\n\n", status);
    status = 0;
  }

#if (DEBUG)
  printf("Created the psrfits file\n");
#endif

  int nblocks, jj, istart,  srownum;
  int *pad_sec, *pad_samp; // array of off_sec and offs_samp
  int skips, pads;
  int *read_dat;
  unsigned char *tmpdata;

  /* Allocate space for arrays we'll need */
  pad_sec = gen_ivect(maxblocks);
  pad_samp = gen_ivect(maxblocks);
  read_dat = gen_ivect(maxblocks);
  tmpdata = gen_bvect(bytes_per_subband);

  /* Allocate memory for raw data buffer */
  rawdata = (unsigned char **) malloc(maxblocks * sizeof(unsigned char*));
  for (ii = 0; ii < maxblocks; ii++)
    rawdata[ii] = gen_bvect(bytes_per_subband);
  
  /* Initialize subint data rows */
  for(ii=0; ii<nchans; ii++){
    pf.sub.dat_weights[ii] = 0.0;
    for(specnum=0; specnum<spec_per_row; specnum++)
      pf.sub.data[specnum * nchans + ii] = (unsigned char) 0;
  }
  
  rownum = 0;
  datidx = 0;
  nblocks = 0;
  istart = 0;

  /* Loop over the Scans */
  for(iscan=0; iscan<nscans; iscan++){
    /* Get number of scans per block and start index */
    nblocks = sblocks[iscan];
    istart = sidx[iscan];
    srownum = 0;
    breakout = 0;
    
    /* Initialize some useful arrays and clear rawdata buffer */
    for(ii=0; ii<nblocks; ii++){
      pad_sec[ii] = fbs[idx[ii + istart]].scan_sec;
      pad_samp[ii] = fbs[idx[ii + istart]].offs_samp;
      read_dat[ii] = 0;
      for(jj=0; jj<bytes_per_subband; jj++){
	rawdata[ii][jj] = (unsigned char) 0;
	tmpdata[jj] = (unsigned char) 0;
      }
    }

    /* Read in data for one scan by row */
    while(1){
      printf("\rWorking on row %d, (Scan %d, Scan Row %d)\r", ++rownum, iscan, ++srownum);
      fflush(stdout);

      for(ii=0; ii<nblocks; ii++){
	/* If rownum < refsec[ii], then we need to skip some rows */
	/* Since this will only happen at start of scan, we can   */
	/* just skip reading in the data, since everything is     */
	/* zero initialiazed.                                     */
	if(rownum < refsec[iscan] + 1){
	  printf("break: rownum = %d, refsec = %d, ii=%d\n", rownum, refsec[iscan], ii);
	  break;
	}

	/* If integer second offsets, fill data buffer with zeros */
	if(pad_sec[ii] > 0){
	  vals_read = 0;
	  for(jj=0; jj<bytes_per_subband; jj++){
	    rawdata[ii][jj] = (unsigned char) 0;
	    vals_read++;
	  }
	  pad_sec[ii] -= 1;
	}
	
	/* If this is the first time reading data (not just padding), */
	/* then we have to check for and apply sample offsets.        */
	else if(read_dat[ii] == 0){
	  /* If block starts earlier than the start time, we have to */
	  /* ignore some samples.                                    */
	  skips = 0;
	  pads = 0;
	  if(pad_samp[ii] < 0){
	    printf("pad_samp[%d] = %d\n", ii, pad_samp[ii]);
	    skips = -1 * pad_samp[ii] * nchans_subband;
	    vals_read = 0;
	    vals_read = fread(tmpdata, sizeof(unsigned char), skips, 
			      infiles[idx[ii + istart]]);
	    if(vals_read < skips){
	      printf("Read %d padding vals (expected %d) for scan %d, block %d\n",
		     vals_read, skips, iscan, ii);
	    }
	    vals_read = 0;
	    vals_read = fread(rawdata[ii], sizeof(unsigned char),
			      bytes_per_subband, infiles[idx[ii + istart]]);
	    read_dat[ii] += 1;
	  }

	  /* If block starts later than the start time, we need to pad */
	  /* some samples.                                             */
	  else if(pad_samp[ii] > 0){
	    printf("pad_samp[%d] = %d\n", ii, pad_samp[ii]);
	    pads = pad_samp[ii] * nchans_subband;
	    vals_read = 0;
	    vals_read = fread(tmpdata, sizeof(unsigned char), bytes_per_subband - pads, 
			      infiles[idx[ii + istart]]);
	    if(vals_read < (bytes_per_subband - pads)){
              printf("Read %d padding vals (expected %d) for scan %d, block %d\n",
                     vals_read, bytes_per_subband - pads, iscan, ii);
            }
	    for(jj=0; jj<bytes_per_subband; jj++){
	      if(jj<pads){
		rawdata[ii][jj] = (unsigned char) 0;
		vals_read++;
	      }
	      else
		rawdata[ii][jj] = tmpdata[jj-pads];
	    }
	    read_dat[ii] += 1;
	  }
	  
	  /* Otherwise, no padding is necessary and we read in normally */
	  else{
	    vals_read = 0;
            vals_read = fread(rawdata[ii], sizeof(unsigned char),
                              bytes_per_subband, infiles[idx[ii + istart]]);
	    read_dat[ii] += 1;
	  }
	}
	
	/* If there are no second offsets, and it's not the first time */
	/* reading in data, just read in normally.                     */
	else{
	  vals_read = 0;
	  vals_read = fread(rawdata[ii], sizeof(unsigned char),
			    bytes_per_subband, infiles[idx[ii + istart]]);
	}
	
	if(vals_read < bytes_per_subband){
	  printf("\n\nRead %d vals (expected %d) for scan %d, block %d\n",
		 vals_read, bytes_per_subband, iscan, ii);
	  printf("inp file: %s, istart = %d\n", fbs[idx[ii + istart]].inpfile, istart);
	  breakout = 1;
	  if(iscan < nscans-1){
	    printf("Done with scan %d, on to scan %d\n", iscan, iscan+1);
	    rownum--; // do this since rownum incremented but not written earlier
	  }
	  if(iscan == nscans-1)
	    printf("Done with last scan, finishing up...\n");
	  break;
	}
      }

      if(breakout == 1)
	break;
      
      /* Write out the data from rawdata to psrfits structure */
      for (ii = 0; ii < nblocks; ii++) {
	fch1idx = (int) roundf((fbs[idx[ii + istart]].fch1 - min_freq)/df);
	fsign = get_sign(fbs[idx[ii + istart]].foff);

	/* Loop over all the spectra per row */
	for (specnum = 0; specnum < spec_per_row; specnum++) {
	  datidx_start = specnum * nchans + fch1idx;
	  rawidx_start = specnum * nchans_subband;
	  for (ichan = 0; ichan < nchans_subband; ichan++) {
	    pf.sub.data[datidx_start + fsign * ichan] = rawdata[ii][rawidx_start + ichan];
	    if(srownum==1 && specnum==0)
	      pf.sub.dat_weights[fch1idx + ichan] = 1.0;
	  }
	}
      }
      
      /* Now write the subint row */
      pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
      status = psrfits_write_subint(&pf);
      if (status) {
	printf("\nError (%d) writing PSRFITS...\n\n", status);
	break;
      }
    }
  }      

  /* Free up everything */
  free(rawdata);
  free(tmpdata);
  free(read_dat);
  free(pad_samp);
  free(pad_sec);
  free(refsec);
  free(sblocks);
  free(sidx);
  free(bidx);
  free(fbs);
  free(idx);
  free(infiles);

  return 0;
}


/* Comparison function to sort sigproc header structs */
int comp_fb(const void *a, const void *b)
{
  sigprocfb *v1, *v2;
  v1 = (sigprocfb *) a;
  v2 = (sigprocfb *) b;
  if (v1->fch1 == v2->fch1)
    return 0;
  else
    {
      if (v1->fch1 < v2->fch1)
        return -1;
      else
        return 1;
    }
}

/* Comparison function to sort bin_index by index */
int comp_bin_index(const void *a, const void *b)
{
  bin_index *v1, *v2;
  v1 = (bin_index *) a;
  v2 = (bin_index *) b;
  if (v1->index == v2->index)
    return 0;
  else
    {
      if (v1->index < v2->index)
        return -1;
      else
        return 1;
    }
}

/* Comparison function to sort bin_index by bin */
int comp_bin_sbin(const void *a, const void *b)
{
  bin_index *v1, *v2;
  v1 = (bin_index *) a;
  v2 = (bin_index *) b;
  if (v1->sbin == v2->sbin)
    return 0;
  else
    {
      if (v1->sbin < v2->sbin)
	return -1;
      else
        return 1;
    }
}

/* Comparison function to sort bin_index by scan, then fch1 */
int comp_bin_scan(const void *a, const void *b)
{
  bin_index *v1, *v2;
  v1 = (bin_index *) a;
  v2 = (bin_index *) b;
  if (v1->scan == v2->scan)
    {
      if(v1->fch1 < v2->fch1)
	return -1;
      else
	return 1;
    }
  else
    {
      if (v1->scan < v2->scan)
	return -1;
      else
        return 1;
    }
}



/* Return the sign of a float */
int get_sign(float x)
{
  return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

/* Initialize a sigprocfb variable */
void init_sigproc(sigprocfb *fb)
{
  strcpy(fb->inpfile, "");     /* Input filename */
  strcpy(fb->source_name, ""); /* Source name */
  strcpy(fb->ifstream, "");    /* Y=IF present, X=IF not present (4 possible values) */
  fb->N = 0;                   /* Number of points (in time) in the file */
  fb->tstart = 0.0;            /* MJD start time */
  fb->tsamp = 0.0;             /* Sampling time in sec */
  fb->src_raj = 0.0;           /* Source RA  (J2000) in hhmmss.ss */
  fb->src_dej = 0.0;           /* Source DEC (J2000) in ddmmss.ss */
  fb->az_start = 0.0;          /* Starting azimuth in deg */
  fb->za_start = 0.0;          /* Starting zenith angle in deg */
  fb->fch1 = 0.0;              /* Highest channel frequency (MHz) */
  fb->foff = 0.0;              /* Channel stepsize (MHz) */
  fb->machine_id = -1;         /* Instrument ID (see backend_name() */
  fb->telescope_id = -1;       /* Telescope ID (see telescope_name() */
  fb->nchans = 0;              /* Number of finterbank channels */
  fb->nbits = 0;               /* Number of bits in the filterbank samples */
  fb->nifs = 0;                /* Number if IFs present */
  fb->nbeams = 0;              /* Number of beams in the observing system */
  fb->ibeam = 0;               /* Beam number used for this data */
  fb->sumifs = 1;              /* Whether the IFs are summed or not */
  fb->headerlen = 0;           /* Header length in bytes */
  fb->imjd = 0;                /* Integer MJD */
  fb->smjd = 0;                /* Seconds past 0h MJD (rounded) */
  fb->offs = 0.0;              /* Offset from integer second */
  fb->offs_samp = 0;           /* Offset in time samples */
  fb->obs_sec = 0;             /* Seconds since start of observation */
  fb->scan_sec = 0;            /* Seconds since start of scan */
  fb->scan = 0;                /* Scan number */
}

/* Check filterbank headers for consistency */
int check_fbs(sigprocfb *fbs, int numfiles)
{
  int ii, retval, tmp;
  char source_name[80];
  char inpfile[80];
  double foff, tsamp;
  int nbits;
  
  retval = 0;
  tmp = 0;
  strcpy(inpfile, fbs[0].inpfile);
  strcpy(source_name, fbs[0].source_name);
  tsamp = fbs[0].tsamp;
  foff = fabs(fbs[0].foff);
  nbits = fbs[0].nbits;
  
  for(ii=1; ii<numfiles; ii++){
    if((tmp = strcmp(source_name, fbs[ii].source_name)) != 0){
      printf("%s and %s have different source_name\n", inpfile, fbs[ii].inpfile);
      printf("%s , %s\n", source_name, fbs[ii].source_name);
      retval++;
    }
    if(foff != fabs(fbs[ii].foff)){
      printf("%s and %s have different |foff|\n", inpfile, fbs[ii].inpfile);
      printf("%f , %f\n", foff, fabs(fbs[ii].foff));
      retval++;
    }
    if(tsamp != fbs[ii].tsamp){
      printf("%s and %s have different tsamp\n", inpfile, fbs[ii].inpfile);
      printf("%f , %f\n", tsamp, fbs[ii].tsamp);
      retval++;
    }
    if(nbits != fbs[ii].nbits){
      printf("%s and %s have different nbits\n", inpfile, fbs[ii].inpfile);
      printf("%d , %d\n", nbits, fbs[ii].nbits);
      retval++;
    }
  }

  return retval;
}

/* Convert raj/dej double to string */
void coord_double2str(double x, char *outstr, int is_ra)
{
  int hh, mm;
  float ss;
  
  hh = (int) (x / 10000.0);
  mm = (int) (fabs(x - hh * 10000) / 100.0);
  ss = fabs(x - hh * 10000) - mm * 100;
  if(is_ra) // Print RA string if is_ra
    sprintf(outstr, "%.2i:%.2i:%07.4f", hh, mm, ss);
  else     // Else pring DEC string
    sprintf(outstr, "%+.2i:%.2i:%07.4f", hh, mm, ss);

}

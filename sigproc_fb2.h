typedef struct SIGPROCFB {
  char inpfile[80];      /* Input filename */
  char source_name[80];  /* Source name */
  char ifstream[8];      /* Y=IF present, X=IF not present (4 possible values) */
  long long N;           /* Number of points (in time) in the file */
  double tstart;         /* MJD start time */
  double tsamp;          /* Sampling time in sec */
  double src_raj;        /* Source RA  (J2000) in hhmmss.ss */
  double src_dej;        /* Source DEC (J2000) in ddmmss.ss */
  double az_start;       /* Starting azimuth in deg */
  double za_start;       /* Starting zenith angle in deg */
  double fch1;           /* Highest channel frequency (MHz) */
  double foff;           /* Channel stepsize (MHz) */
  int machine_id;        /* Instrument ID (see backend_name() */
  int telescope_id;      /* Telescope ID (see telescope_name() */
  int nchans;            /* Number of finterbank channels */
  int nbits;             /* Number of bits in the filterbank samples */
  int nifs;              /* Number if IFs present */
  int nbeams;            /* Number of beams in the observing system */
  int ibeam;             /* Beam number used for this data */
  int sumifs;            /* Whether the IFs are summed or not */
  int headerlen;         /* Header length in bytes */
  int imjd;              /* Integer MJD */
  int smjd;              /* Seconds past 0h MJD (rounded) */
  double offs;           /* Offset from integer second */
  int offs_samp;         /* Offset in time samples */
  int obs_sec;           /* Seconds since start of observation */
  int scan_sec;          /* Seconds since the start of the scan */
  int scan;              /* Scan number */
} sigprocfb;

/* sigproc_fb.c */
int read_filterbank_header(sigprocfb *fb, FILE *inputfile);
void write_filterbank_header(sigprocfb *fb, FILE *outfile);
void print_filterbank_header(sigprocfb *fb);
void get_filterbank_static(int *ptsperbyte, int *bytesperpt, 
                           int *bytesperblk, float *clip_sigma);
void set_filterbank_static(int ptsperbyte, int ptsperblk, 
                           int bytesperpt, int bytesperblk, 
                           int numchan, float clip_sigma, double dt);
void set_filterbank_padvals(float *fpadvals, int good_padvals);
int skip_to_filterbank_rec(FILE *infiles[], int numfiles, int rec);
int read_filterbank_rawblock(FILE *infiles[], int numfiles, unsigned char *data, 
			     int *padding);
int read_filterbank_rawblocks(FILE *infiles[], int numfiles, unsigned char rawdata[], 
			      int numblocks, int *padding);
void get_filterbank_channel(int channum, float chandat[], unsigned char rawdata[], 
			    int numblocks);
void convert_filterbank_block(int *indata, unsigned char *outdata);

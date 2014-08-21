#include <sys/times.h>
#include "chkio.h"
#include <math.h>
#include <string.h>
#include <time.h>
#include "clk_tck.h"
#include "meminfo.h"


#ifndef SQRT2
#define SQRT2         1.4142135623730950488016887242096980785696718753769
#endif
#ifndef DBLCORRECT
#define DBLCORRECT    1e-14
#endif
#ifndef PI
#define PI            3.1415926535897932384626433832795028841971693993751
#endif
#ifndef TWOPI
#define TWOPI         6.2831853071795864769252867665590057683943387987502
#endif
#ifndef DEGTORAD
#define DEGTORAD      0.017453292519943295769236907684886127134428718885417
#endif
#ifndef RADTODEG
#define RADTODEG      57.29577951308232087679815481410517033240547246656
#endif
#ifndef PIBYTWO
#define PIBYTWO       1.5707963267948966192313216916397514420985846996876
#endif
#ifndef SOL
#define SOL           299792458.0
#endif
#ifndef SECPERJULYR
#define SECPERJULYR   31557600.0
#endif
#ifndef SECPERDAY
#define SECPERDAY     86400.0
#endif
#ifndef ARCSEC2RAD
#define ARCSEC2RAD    4.8481368110953599358991410235794797595635330237270e-6
#endif
#ifndef SEC2RAD
#define SEC2RAD       7.2722052166430399038487115353692196393452995355905e-5
#endif
#ifndef __GNUC__
#define __inline__
#endif
/* Maximum number of input files to try and patch together */
#define MAXPATCHFILES 100
/* Blocksize to use when reading datafiles or subbands */
#define SUBSBLOCKLEN 1024

/* various function-like macros */

#ifndef SWAP
/* Swaps two variables of undetermined type */
#define SWAP(a,b) tempzz=(a);(a)=(b);(b)=tempzz;
#endif

#ifndef POWER
/* Returns unnormalized Fourier power  */
/*   Requires the following variables in calling function */
/*   double powargr, powargi; */
#define POWER(r,i) (powargr=(r),powargi=(i),\
		    powargr*powargr+powargi*powargi)
#endif

#ifndef PHASE
/* Returns Fourier phase (degrees)  */
/*   Requires the following variables in calling function */
/*   double phsargr, phsargi, phstmp; */
#define PHASE(r,i) (phsargr=(r),phsargi=(i),\
		    ((phstmp=RADTODEG*atan2(phsargi,phsargr)) > 0.0) ? \
		    phstmp : phstmp+360.0)
#endif
  
#ifndef RADIAN_PHASE
/* Returns Fourier phase (radians)  */
/*   Requires the following variables in calling function */
/*   double radargr, radargi, radtmp; */
#define RADIAN_PHASE(r,i) (radargr=(r),radargi=(i),\
		    ((radtmp=atan2(radargi,radargr)) > 0.0) ? \
		    radtmp : radtmp+TWOPI)
#endif

#define GET_BIT(c, n) (*(c+(n>>3)) >> (7-(n&7)) & 1)
#define SET_BIT(c, n) (*(c+(n>>3)) |= 1 << (7-(n&7)))
#define UNSET_BIT(c, n) (*(c+(n>>3)) &= ~(1 << (7-(n&7))))

/*   Number of bins (total) to average for local power:   */
/*     Must be an even number (1/2 on each side).         */
#define NUMLOCPOWAVG  20

/*  Number of bins next to freq in question to            */
/*    ignore (on each side) when determining local power. */
#define DELTAAVGBINS  5

/*  Number of bins on each side of the central frequency  */
/*    to sum for Fourier interpolation (low accuracy)     */
#define NUMFINTBINS   16

/* Used for raw-data handling */
typedef enum {
  IF0, IF1, SUMIFS
} IFs;

/*  Constants used in the interpolation routines */
typedef enum {
  LOWACC, HIGHACC
} presto_interp_acc;

/*  Constants used in the binary search routines */
typedef enum {
  INTERBIN, INTERPOLATE
} presto_interptype;

typedef enum {
  NO_CHECK_ALIASED, CHECK_ALIASED
} presto_checkaliased;

/*  Constants used in the correlation/convolution routines */
typedef enum {
  CONV, CORR, INPLACE_CONV, INPLACE_CORR
} presto_optype;

typedef enum {
  FFTDK, FFTD, FFTK, NOFFTS
} presto_ffts;

typedef enum {
  RAW, PREPPED, FFT, SAME
} presto_datainf;


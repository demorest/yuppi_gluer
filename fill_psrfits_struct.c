#include <stdio.h>
#include <math.h>
#include <string.h>
#include "chkio.h"
#include "vectors2.h"
#include "psrfits.h"
//#include "slalib.h"

#ifndef DEGTORAD
#define DEGTORAD 0.017453292519943295769236907684886127134428718885417
#endif
#ifndef RADTODEG
#define RADTODEG 57.29577951308232087679815481410517033240547246656
#endif
#ifndef SOL
#define SOL 299792458.0
#endif
#ifndef SWAP
/* Swaps two variables of undetermined type */
#define SWAP(a,b) tempzz=(a);(a)=(b);(b)=tempzz;
#endif

// Return the beam FWHM in degrees for obs_freq in MHz
// and dish_diam in m
static double beam_FWHM(double obs_freq, double dish_diam)
{
    double lambda = SOL / (obs_freq * 1e6);
    return 1.2 * lambda / dish_diam * RADTODEG;
}


void fill_psrfits_struct(int numbits, struct psrfits *pf, float dt, char* source, int nchantot, float df, int imjd, int smjd, float offs, char* ra, char* dec, float lofreq, char *basename)
{
    int ii;

    pf->filenum = 0;            // Crucial for initialization
    strcpy(pf->basefilename, basename);
    pf->hdr.nsblk = (int) (roundf(1.0 / dt));

    // Now set values for our hdrinfo structure
    strcpy(pf->hdr.telescope, "VLA");
    strcpy(pf->hdr.obs_mode, "SEARCH");
    strcpy(pf->hdr.backend, "YUPPI");
    strcpy(pf->hdr.frontend, "HFF");
    strcpy(pf->hdr.observer, " ");
    strcpy(pf->hdr.project_id, " ");
    strcpy(pf->hdr.source, source);
    strcpy(pf->hdr.date_obs, " ");
    //pf->hdr.scanlen = get_hdr_double(h, "obs_time");

    strcpy(pf->hdr.poln_type, "LIN");   // set based on known receivers
    strcpy(pf->hdr.poln_order, "AA");
    pf->hdr.summed_polns = 0;
    
    strcpy(pf->hdr.track_mode, "TRACK");  // Potentially not-true?
    strcpy(pf->hdr.cal_mode, "OFF");      // Potentially not-true?
    strcpy(pf->hdr.feed_mode, "FA");      // check this...

    pf->hdr.dt = dt;
    pf->hdr.fctr = lofreq + (nchantot/2.0 - 0.5)*df;
    pf->hdr.BW = nchantot*df;
    pf->hdr.beam_FWHM = beam_FWHM(pf->hdr.fctr, 34.0);
    pf->hdr.nchan = nchantot;
    pf->hdr.orig_nchan = nchantot;
    pf->hdr.df = df;
    pf->hdr.orig_df = df;
    pf->hdr.nbits = numbits;
    pf->hdr.npol = 1;
    pf->hdr.start_day = imjd;
    pf->hdr.start_sec = smjd;
    pf->hdr.MJD_epoch = (long double)imjd + ((long double)smjd + (long double)offs)/((long double)(24*3600));
    pf->hdr.scan_number = 0;
    strcpy(pf->hdr.ra_str, ra);
    strcpy(pf->hdr.dec_str, dec);
    
    pf->hdr.ra2000 = 0.0;
    pf->hdr.dec2000 = 0.0;

    pf->hdr.azimuth = 0.0;
    pf->hdr.zenith_ang = 0.0;
    pf->hdr.rcvr_polns = 1;
    pf->hdr.offset_subint = 0;
    pf->hdr.onlyI = 0;
    pf->hdr.ds_time_fact = 1;
    pf->hdr.ds_freq_fact = 1;
    pf->hdr.chan_dm = 0.0;
    pf->hdr.fd_hand = 1.0;
    pf->hdr.be_phase = 0;     // This is almost certainly not correct
    pf->hdr.fd_sang = pf->hdr.fd_xyph = 0.0;    // This is almost certainly not correct
    pf->hdr.feed_angle = 0.0;   // This is almost certainly not correct
    pf->hdr.cal_freq = pf->hdr.cal_dcyc = pf->hdr.cal_phs = 0.0;        //  ditto

    // Now set values for our subint structure
    pf->sub.tel_az = 0.0;
    pf->sub.tel_zen = 0.0;
    pf->sub.lst = 0.0;
    pf->sub.tsubint = pf->hdr.nsblk * pf->hdr.dt;
    pf->sub.ra = pf->hdr.ra2000;
    pf->sub.dec = pf->hdr.dec2000;
    pf->sub.offs = 0.5 * pf->sub.tsubint;
    //slaEqgal(pf->hdr.ra2000 * DEGTORAD, pf->hdr.dec2000 * DEGTORAD,
    //         &pf->sub.glon, &pf->sub.glat);
    //pf->sub.glon *= RADTODEG;
    //pf->sub.glat *= RADTODEG;
    // The following three are unknown or hard to get, I think (SMR)
    pf->sub.feed_ang = 0.0;
    pf->sub.pos_ang = 0.0;
    pf->sub.par_ang = 0.0;

    pf->sub.bytes_per_subint = (pf->hdr.nbits/8) * pf->hdr.nchan * pf->hdr.npol * pf->hdr.nsblk;

    pf->sub.FITS_typecode = TBYTE;        // newer, back to bytes
    //pf->sub.FITS_typecode = TFLOAT;     // new
    //pf->sub.FITS_typecode = TBYTE;      // 11 = byte

    // Create and initialize the subint arrays
    pf->sub.dat_freqs = gen_fvect(pf->hdr.nchan);
    pf->sub.dat_weights = gen_fvect(pf->hdr.nchan);
    for (ii = 0; ii < pf->hdr.nchan; ii++) {
        pf->sub.dat_freqs[ii] = lofreq + ii * pf->hdr.df;
        //pf->sub.dat_weights[ii] = 1.0;
    }

    // The following are re-set to try to preserve the band shape later
    pf->sub.dat_offsets = gen_fvect(pf->hdr.nchan * pf->hdr.npol);
    pf->sub.dat_scales = gen_fvect(pf->hdr.nchan * pf->hdr.npol);
    for (ii = 0; ii < pf->hdr.nchan * pf->hdr.npol; ii++) {
        pf->sub.dat_offsets[ii] = 0.0;
        pf->sub.dat_scales[ii] = 1.0;
    }

    // This is the raw data block that will be updated 
    // for each row of the PSRFITS file

    pf->sub.data = gen_bvect( pf->sub.bytes_per_subint );   // newest 
    //pf->sub.data = gen_bvect( (pf->sub.bytes_per_subint)/4 );   // newer
    //pf->sub.data = gen_fvect( (pf->sub.bytes_per_subint)/4 ); // new
    //pf->sub.data = gen_bvect(pf->sub.bytes_per_subint);
}


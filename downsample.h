/* downsample.c */

void convert_4bit_to_8bit(unsigned char *indata, unsigned char *outdata, int N);
void pf_4bit_to_8bit(struct psrfits *pf);
void convert_8bit_to_4bit(unsigned char *indata, unsigned char *outdata, int N);
void pf_8bit_to_4bit(struct psrfits *pf);
void get_stokes_I(struct psrfits *pf);
void downsample_time(struct psrfits *pf);
void guppi_update_ds_params(struct psrfits *pf);

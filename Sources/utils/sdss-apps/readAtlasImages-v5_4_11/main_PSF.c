/*
 * Read a psField file produced by photo
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dervish.h"
#include "phConsts.h"
#include "phVariablePsf.h"
#include "phFits.h"
/*
 * some symbols to prevent dervish.o from being loaded from libatlas.a
 * if we are also linking against real dervish
 */

int verbose = 0;

static void usage(void);

int
main(int ac, char *av[])
{
   FITS *fits;				/* the table in question */
   char *infile, *outfile;		/* input and output filenames */
   int hdu;				/* desired table in file */
   PSF_BASIS *basis = NULL;		/* PSF_BASIS derived from kl_comps */
   PSF_KL_COMP *klc;			/* PSF_KL_COMP read from file */
   PSF_REG * psf_reg;			/*    desired PSF */
   PIXDATATYPE regType;			/* type of desired REGION */
   int row;				/* desired row */
   float rowc, colc;			/* desired position in frame */

   while(ac > 1 && (av[1][0] == '-' || av[1][0] == '+')) {
      switch (av[1][1]) {
       case '?':
       case 'h':
	 usage();
	 exit(0);
	 break;
       case 'i':
	 fprintf(stderr,"SDSS read_PSF. Id: %s\n", phPhotoVersion());
	 exit(0);
	 break;
       case 'v':
	 verbose++;
	 break;
       default:
	 shError("Unknown option %s\n",av[1]);
	 break;
      }
      ac--;
      av++;
   }
   if(ac <= 5) {
      shError("\
You must specify an input file, an HDU, a row and column position\n\
and an output file\n");
      exit(1);
   }
   infile = av[1]; hdu = atoi(av[2]);
   rowc = atof(av[3]); colc = atof(av[4]);
   outfile = av[5];
/*
 * dummy calls to pull .o files out of the real libdervish.a if we are
 * linking against it
 */
   (void)shTypeGetFromName("RHL");
/*
 * open file
 */
   if((fits = open_fits_table(infile, hdu)) == NULL) {
      exit(1);
   }
/*
 * Read KL components
 */
   for(row = 1; row <= fits->naxis2; row++) {
      if((klc = read_KLC(fits, row)) == NULL) {
	 exit(1);
      }

/*
 * Create the PSF_BASIS; we need the first PSF_KL_BASIS to do this
 */
      if(basis == NULL) {
	 PSF_KERNEL *kern = phPsfKernelNew(fits->naxis2 - 1, 1, 1,
					   klc->nrow_b, klc->ncol_b);
	 basis = phPsfBasisNew(kern, KL_BASIS, -1, -1, klc->reg->type);
	 phPsfKernelDel(kern);		/* decrement reference counter */
      }
      
      phPsfBasisSetFromKLComp(basis, row - 1, klc, 0);

      klc->reg = NULL; phPsfKLCompDel(klc);
   }
/*
 * Reconstruct PSF at specified point.
 */
#if 1
   regType = TYPE_U16;			/* create an unsigned short region */
#else
   regType = TYPE_FL32;			/* create a float region */
#endif
   psf_reg = phPsfKLReconstruct(basis, rowc, colc, regType);
/*
 * and write it to a file
 */
   write_fits_file(outfile, psf_reg->reg, (regType == TYPE_U16) ? 16 : -32);

   phPsfRegDel(psf_reg);
   phPsfBasisDel(basis);

   phFitsDel(fits);
   
   return(0);
}

/*****************************************************************************/

static void
usage(void)
{
   char **line;

   static char *msg[] = {
      "Usage: read_PSF [options] input-file hdu output-file",
      "Your options are:",
      "       -?      This message",
      "       -h      This message",
      "       -i      Print an ID string and exit",
      "       -v      Turn up verbosity (repeat flag for more chatter)",
      NULL,
   };

   for(line = msg;*line != NULL;line++) {
      fprintf(stderr,"%s\n",*line);
   }
}

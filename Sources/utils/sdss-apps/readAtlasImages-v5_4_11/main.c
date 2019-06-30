/*
 * Read an atlas image table produced by photo
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dervish.h"
#include "phFits.h"
#include "phConsts.h"
/*
 * some symbols to prevent dervish.o from being loaded from libatlas.a
 * if we are also linking against real dervish
 */

int verbose = 0;

static void set_background(REGION *reg, int bkgd);
static void usage(void);

int
main(int ac, char *av[])
{
   ATLAS_IMAGE *ai;			/* the desired atlas image */
   int bkgd = SOFT_BIAS;		/* desired background level */
   int color = 0;			/* desired color */
   FITS *fits;				/* the table in question */
   char *infile, *outfile;		/* input and output filenames */
   int row0, col0;			/* origin of image in region */
   REGION *reg;				/* region to write */
   int row;				/* desired row */

   while(ac > 1 && (av[1][0] == '-' || av[1][0] == '+')) {
      switch (av[1][1]) {
       case '?':
       case 'h':
	 usage();
	 exit(0);
	 break;
       case 'b':
	 if(ac == 2) {
	    shError("-b requires a number");
	 } else {
	    ac--; av++;
	    bkgd = atoi(av[1]);
	 }
	 break;
       case 'c':
	 if(ac == 2) {
	    shError("-c requires a number");
	 } else {
	    ac--; av++;
	    color = atoi(av[1]);
	 }
	 break;
       case 'i':
	 fprintf(stderr,"SDSS read_atlas_images. Id: %s\n", phPhotoVersion());
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
   if(ac <= 3) {
      shError("You must specify an input file, a row, and an output file\n");
      exit(1);
   }
   infile = av[1]; row = atoi(av[2]); outfile = av[3];
/*
 * dummy calls to pull .o files out of the real libdervish.a if we are
 * linking against it
 */
   (void)shTypeGetFromName("RHL");
/*
 * open file
 */
   if((fits = open_fits_table(infile, 1)) == NULL) {
      exit(1);
   }
/*
 * read atlas image
 */
   if((ai = read_atlas_image(fits,row)) == NULL) {
      exit(1);
   }
   if(ai->id < 0) {			/* no atlas image for this object */
      shError("Object %d has no atlas image", row);
      exit(1);
   }
   shAssert(ai->master_mask != NULL);

   if(color < 0 || color >= ai->ncolor) {
      shError("Invalid color; please choose a number in 0..%d", ai->ncolor-1);
      exit(1);
   }
/*
 * convert it to a region
 */
   reg = shRegNew("atlas image",
		  ai->master_mask->rmax - ai->master_mask->rmin + 1,
		  ai->master_mask->cmax - ai->master_mask->cmin + 1, TYPE_U16);
   set_background(reg, bkgd);
   row0 = ai->master_mask->rmin + ai->drow[color];
   col0 = ai->master_mask->cmin + ai->dcol[color];
   phRegionSetFromAtlasImage(ai, color, reg, row0, col0, 0.0);
/*
 * and write it to a file
 */
   write_fits_file(outfile, reg, 16);

   phAtlasImageDel(ai, 1);
   shRegDel(reg);
   phFitsDel(fits);
   
   return(0);
}

/*****************************************************************************/

static void
usage(void)
{
   char **line;

   static char *msg[] = {
      "Usage: read_atlas_image [options] input-file row output-file",
      "Your options are:",
      "       -?      This message",
      "       -b #    Set background level to #",
      "       -c #    Use colour # (0..ncolor-1; default 0)",
      "       -h      This message",
      "       -i      Print an ID string and exit",
      "       -v      Turn up verbosity (repeat flag for more chatter)",
      NULL,
   };

   for(line = msg;*line != NULL;line++) {
      fprintf(stderr,"%s\n",*line);
   }
}

/*****************************************************************************/

static void
set_background(REGION *reg,
	       int bkgd)
{
   U16 *row0;
   int i;

   row0 = reg->rows[0];
   for(i = 0; i < reg->ncol; i++) {
      row0[i] = bkgd;
   }
   for(i = 1; i < reg->nrow; i++) {
      memcpy(reg->rows[i], row0, reg->ncol*sizeof(U16));
   }
}

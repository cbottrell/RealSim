/*
 * Read an fpM file produced by photo
 */
#define CALC_POLYGONS 1			/* define this if you want to use
					   the -p flag to calculate polgonal
					   approximations to masks */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dervish.h"
#include "phFits.h"
#include "phConsts.h"
#include "phSpanUtil.h"
#if CALC_POLYGONS
#  include "phGeometry.h"
#endif
/*
 * some symbols to prevent dervish.o from being loaded from libatlas.a
 * if we are also linking against real dervish
 */

int verbose = 0;

static void usage(void);
static int get_hdu(const char *str);

int
main(int ac, char *av[])
{
   FITS *fits;				/* the table in question */
   char *infile, *outfile;		/* input and output filenames */
   int hdu;				/* desired table in file */
   MASK *mask;				/* mask to write */
   OBJMASK *om;				/* objmask read from file */
   FILE *out = NULL;			/* output file for polygons */
   POLYGONS *poly;			/* polygonal approximation to om */
   int polygons = 0;			/* find polygonal approximation? */
   int row;				/* desired row */

   while(ac > 1 && (av[1][0] == '-' || av[1][0] == '+')) {
      switch (av[1][1]) {
       case '?':
       case 'h':
	 usage();
	 exit(0);
	 break;
       case 'i':
	 fprintf(stderr,"SDSS read_mask. Id: %s\n", phPhotoVersion());
	 exit(0);
	 break;
#if CALC_POLYGONS
       case 'p':
	 polygons = 1;
	 break;
#endif
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
      shError("You must specify an input file, a table, and an output file\n");
      exit(1);
   }
   infile = av[1]; hdu = get_hdu(av[2]); outfile = av[3];
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
 * Create MASK
 */
   if(polygons) {
      out = fopen(outfile, "w");
      if(out == NULL) {
	 fprintf(stderr,"Cannot open %s for write\n", outfile);
	 exit(1);
      }
      mask = NULL;
   } else {
      mask = shMaskNew("from file", fits->maskrows, fits->maskcols);
      shMaskClear(mask);
   }
/*
 * read OBJMASKs and set bits in mask
 */
   if(verbose) {
      if(polygons) {
	 printf("Writing polygons to %s\n", outfile);
      } else {
	 printf("reading %4d OBJMASKs from HDU %2d into a %dx%d MASK\n",
		fits->naxis2, hdu, mask->nrow, mask->ncol);
      }
   }
   for(row = 1; row <= fits->naxis2; row++) {
      if((om = read_objmask(fits, row)) == NULL) {
	 exit(1);
      }

      if(polygons) {
#if CALC_POLYGONS
	 char hdr[100];
	 sprintf(hdr, "%d (npix %d)", om->id, om->npix);
	 poly = phPolygonsFromObjmask(om);
	 phPolygonsPrint(out, hdr, poly);
	 phPolygonsDel(poly);
	 if(row < fits->naxis2) {
	    fprintf(out, "\n");
	 }
#endif
      } else {
	 phMaskSetFromObjmask(om, mask, 1);
      }

      phObjmaskDel(om);
   }
/*
 * and write it to a file
 */
   if(!polygons) {
      write_fits_file(outfile, mask, 8);
      shMaskDel(mask);
   }
/*
 * cleanup
 */
   if(out != NULL) {
      fclose(out);
   }
   phFitsDel(fits);
   
   return(0);
}

/*****************************************************************************/
/*
 * Return the desired HDU
 */
static int
get_hdu(const char *str)
{
   char *end = NULL;
   int hdu;

   hdu = strtol(str, &end, 10);
   if(end != str) {			/* parsing succeeded */
      return(hdu);
   }

   if(strcmp(str, "INTERP") == 0) {	/* pixel's value interpolated */
      return((int)S_MASK_INTERP + 1);
   } else if(strcmp(str, "SATUR") == 0) { /* is/was saturated  */
      return((int)S_MASK_SATUR + 1);
   } else if(strcmp(str, "NOTCHECKED") == 0) { /* NOT examined for an object */
      return((int)S_MASK_NOTCHECKED + 1);
   } else if(strcmp(str, "OBJECT") == 0) { /* part of some object */
      return((int)S_MASK_OBJECT + 1);
   } else if(strcmp(str, "BRIGHTOBJECT") == 0) { /* part of bright object */
      return((int)S_MASK_BRIGHTOBJECT + 1);
   } else if(strcmp(str, "BINOBJECT") == 0) { /* part of binned object */
      return((int)S_MASK_BINOBJECT + 1);
   } else if(strcmp(str, "CATOBJECT") == 0) { /* part of a catalogued object */
      return((int)S_MASK_CATOBJECT + 1);
   } else if(strcmp(str, "SUBTRACTED") == 0) { /* vals subtracted from pixel */
      return((int)S_MASK_SUBTRACTED + 1);
   } else if(strcmp(str, "GHOST") == 0) { /* part of a ghost */
      return((int)S_MASK_GHOST + 1);
   } else if(strcmp(str, "CR") == 0) { /* part of a cosmic ray */
      return((int)S_MASK_CR + 1);
   } else {
      fprintf(stderr, "Invalid bitplane: %s\n", str);
      exit(1);
   }

   return(0);				/* NOTREACHED */
}

/*****************************************************************************/

static void
usage(void)
{
   char **line;

   static char *msg[] = {
      "Usage: read_mask [options] input-file hdu output-file",
      "The HDU may be a number (>= 1), or a name such as CR or INTERP",
      "If -p is absent, the ourput-file will be a 8-bit FITS image;",
      "If -p is specified, the ourput-file will be a list of vertices.",
      "",
      "Your options are:",
      "       -?      This message",
      "       -h      This message",
      "       -p      Calculate polygonal approximations (if so compiled)",
      "       -i      Print an ID string and exit",
      "       -v      Turn up verbosity (repeat flag for more chatter)",
      NULL,
   };

   for(line = msg;*line != NULL;line++) {
      fprintf(stderr,"%s\n",*line);
   }
}

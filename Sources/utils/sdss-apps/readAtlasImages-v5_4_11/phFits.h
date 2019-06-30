#if !defined(PHFITS_H)
#define PHFITS_H

#include <fcntl.h>
#include "phObjc.h"

typedef struct {
   int fd, hfd;				/* file descriptors for main table
					   and heap */
   int bitpix;				/* value of bitpix */
   int naxis;				/* number of axes */
   int naxis1, naxis2;			/* size of table */
   int theap;				/* offset of heap */
   int pcount;				/* size of heap */
   off_t data_start;			/* start of data */
   off_t heap_start;			/* start of heap */

   int maskrows, maskcols;		/* read from PDU of fpM files */
} FITS;

FITS *phFitsNew(const char *file);
void phFitsDel(FITS *fits);

FITS *open_fits_table(const char *name, int hdu);
ATLAS_IMAGE *read_atlas_image(const FITS *fits, int row);
OBJMASK *read_objmask(const FITS *fits, int row);
#if defined(PH_VARIABLE_PSF_H)
   PSF_KL_COMP *read_KLC(const FITS *fits, int row);
#endif

void write_fits_file(const char *name, const void *reg, int bitpix);
#endif

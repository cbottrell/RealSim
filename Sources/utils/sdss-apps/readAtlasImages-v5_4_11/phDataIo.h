#if !defined PHDATAIO_H
#define PHDATAIO_H

#include "phObjc.h"
#include "phRice.h"

/*****************************************************************************/
/*
 * flatten/inflate OBJMASKs to/from strings
 */
int
phObjmaskFlatten(const OBJMASK *om,	/* OBJMASK to flatten */
		 unsigned char *buff,	/* buffer to flatten into; or NULL */
		 int len);		/* dimension of buff */
int
phSpanmaskFlatten(const SPANMASK *sm,	/* SPANMASK to flatten */
		  unsigned char *buff,	/* buffer to flatten into; or NULL */
		  int len);		/* dimension of buff */
int
phAtlasImageFlatten(const ATLAS_IMAGE *ai, /* ATLAS_IMAGE to flatten */
		    unsigned char *buff, /* buffer to flatten into; or NULL */
		    int len);		/* dimension of buff */

OBJMASK *
phObjmaskInflate(OBJMASK *om,		/* OBJMASK to fill, or NULL */
		 unsigned char *buff,	/* buffer to inflate from */
		 int *plen);		/* if non-NULL:
					   on input, length of buff or 0
					   on output, number of bytes read */
SPANMASK *
phSpanmaskInflate(SPANMASK *sm,		/* SPANMASK to fill, or NULL */
		  unsigned char *buff,	/* buffer to inflate from */
		  int *plen);		/* if non-NULL:
					   on input, length of buff or 0
					   on output, number of bytes read */
ATLAS_IMAGE *
phAtlasImageInflate(ATLAS_IMAGE *ai,	/* ATLAS_IMAGE to fill, or NULL */
		  unsigned char *buff,	/* buffer to inflate from */
		  int *plen);		/* if non-NULL:
					   on input, length of buff or 0
					   on output, number of bytes read */


#endif

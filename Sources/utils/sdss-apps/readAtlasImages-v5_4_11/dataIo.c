/*
 * The code in this file is used to flatten/inflate complex photo data types
 * into network-byteorder strings for output to FITS binary tables
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dervish.h"
#include "phSpanUtil.h"
#include "phObjc.h"
#include "phDataIo.h"

/*****************************************************************************/
/*
 * utility functions to byte swap data
 */
#if defined(SDSS_LITTLE_ENDIAN)
static void
swap2(unsigned char *buff,
      int n)
{
   int i;
   char tmp;

   for(i = 0; i < n; i += 2) {
      tmp = buff[i]; buff[i] = buff[i + 1]; buff[i + 1] = tmp;
   }
}

static void
swap24(unsigned char *buff,
      int n)
{
   int i;
   char tmp;

   for(i = 0; i < n; i += 4) {
      tmp = buff[i];     buff[i] =     buff[i + 3]; buff[i + 3] = tmp;
      tmp = buff[i + 1]; buff[i + 1] = buff[i + 2]; buff[i + 2] = tmp;
   }
}
#else
static void
swap2(unsigned char *buff,		/* NOTUSED */
     int n)				/* NOTUSED */
{}

static void
swap24(unsigned char *buff,		/* NOTUSED */
      int n)				/* NOTUSED */
{}
#endif

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Flatten/inflate an SPANMASK into/from a network-byte-order string of
 * bytes.
 *
 * The buffer must be long enough; if it isn't an assertion will fail. You
 * can find out the needed size by first calling the routine with buff == NULL,
 */
int
phSpanmaskFlatten(const SPANMASK *sm,	/* SPANMASK to flatten */
		  unsigned char *buff,	/* buffer to flatten into; or NULL */
		  int len)		/* dimension of buff */
{
   unsigned char *base;			/* saved value of buff at misc times */
   unsigned char *const buff0 = buff;	/* initial value of buff */
   int j,k;
   int nel;				/* number of elements on a chain */
   int nrow, ncol;			/* == sm->n{row,col} */
   OBJMASK *om;				/* an OBJMASK on one of the chains */
   int size;				/* size required */
   
   size = 2*sizeof(int);
   if(sm != NULL) {
      for(j = 0; j < NMASK_TYPES; j++) {
	 size += sizeof(int);
	 nel = (sm->masks[j] == NULL) ? 0 : sm->masks[j]->nElements;
	 for(k = 0; k < nel; k++) {
	    om = shChainElementGetByPos(sm->masks[j], k);
	    size += phObjmaskFlatten(om, NULL, 0);
	 }
      }
   }
   
   if(buff == NULL) {
      return(size);
   }
   shAssert(len >= size);
/*
 * first the ints, describing the size
 */
   base = buff;
   nrow = (sm == NULL) ? 0 : sm->nrow;
   ncol = (sm == NULL) ? 0 : sm->ncol;
   memcpy(buff, (void *)&nrow, sizeof(int)); buff += sizeof(int);
   memcpy(buff, (void *)&ncol, sizeof(int)); buff += sizeof(int);
   swap24(base, buff - base);
   if(sm == NULL) {
      shAssert(buff <= buff0 + len);
      return(size);
   }
/*
 * then the OBJMASKs themselves
 */
   for(j = 0; j < NMASK_TYPES; j++) {
      nel = (sm == NULL || sm->masks[j] == NULL) ? 0 : sm->masks[j]->nElements;
      base = buff;
      memcpy(buff, (void *)&nel, sizeof(int)); buff += sizeof(int);
      swap24(base, buff - base);

      for(k = 0; k < nel; k++) {
	 om = shChainElementGetByPos(sm->masks[j], k);
	 buff += phObjmaskFlatten(om, buff, buff0 + len - buff);
      }
   }
   shAssert(buff <= buff0 + len);
   
   return(size);
}

/*
 * <AUTO EXTRACT>
 *
 * Inflate an SPANMASK from a network-byte-order string of bytes.
 */
SPANMASK *
phSpanmaskInflate(SPANMASK *sm,		/* SPANMASK to fill, or NULL */
		  unsigned char *buff,	/* buffer to inflate from */
		  int *plen)		/* if non-NULL:
					   on input, length of buff or 0
					   on output, number of bytes read */
{
   unsigned char *const buff0 = buff;	/* initial value of buff */
   int j,k;
   int len;				/* number of bytes read */
   int nel = 0;				/* number of elements on a chain */
   int nrow = 0, ncol = 0;		/* spanmask's size */
   static char *objmask_type = (char *)UNKNOWN;
   OBJMASK *om;				/* an OBJMASK on one of the chains */
   
   if(objmask_type == (char *)UNKNOWN) {
      objmask_type = (char *)shTypeGetFromName("OBJMASK");
   }

   shAssert(buff != NULL);
/*
 * first the ints, describing the size
 */
   swap24(buff,2*sizeof(int));
   memcpy((void *)&nrow, buff, sizeof(int)); buff += sizeof(int);
   memcpy((void *)&ncol, buff, sizeof(int)); buff += sizeof(int);

   if(nrow == 0 && ncol == 0) {
      phSpanmaskDel(sm);
      sm = NULL;
   } else {
      if(sm == NULL) {
	 sm = phSpanmaskNew(nrow, ncol);
      }
/*
 * then the OBJMASKs themselves
 */
      for(j = 0; j < NMASK_TYPES; j++) {
	 swap24(buff,sizeof(int));
	 memcpy((void *)&nel, buff, sizeof(int)); buff += sizeof(int);
	 sm->masks[j] = (nel == 0) ? NULL : shChainNew(objmask_type);

	 for(k = 0; k < nel; k++) {
	    len = (plen == NULL || *plen == 0) ? 0 : (buff0 + *plen) - buff;
	    
	    om = phObjmaskInflate(NULL, buff, &len); buff += len;
	    shChainElementAddByPos(sm->masks[j], om, objmask_type, k, AFTER);
	 }
      }
   }

   if(plen != NULL) {
      if(*plen > 0) {
	 shAssert(buff - buff0 <= *plen);
      }
      *plen = buff - buff0;
   }
   
   return(sm);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Flatten/inflate an OBJMASK into/from a network-byte-order string of
 * bytes.
 *
 * The buffer must be long enough; if it isn't an assertion will fail. You
 * can find out the needed size by first calling the routine with buff == NULL,
 */
int
phObjmaskFlatten(const OBJMASK *om,	/* OBJMASK to flatten */
		 unsigned char *buff,	/* buffer to flatten into; or NULL */
		 int len)		/* dimension of buff */
{
   unsigned char *base;			/* saved value of buff at misc times */
   unsigned char *const buff0 = buff;	/* initial value of buff */
   int i;
   int npix, nspan;			/* == om->n{pix,span} */
   int size;				/* size required */
   
   size = sizeof(int);			/* we always write nspan */
   if(om == NULL || om->nspan == 0) {
      nspan = npix = 0;
   } else {
      nspan = om->nspan; npix = om->npix;

      size += 8*sizeof(int) + 3*nspan*sizeof(short);
      if(om->data != NULL) {
	 size += npix*sizeof(PIX);
      }
   }
   
   if(buff == NULL) {
      return(size);
   }
   shAssert(len >= size);
/*
 * first the ints, describing e.g. the size and bounding box
 */
   base = buff;
   memcpy(buff, (void *)&nspan, sizeof(int)); buff += sizeof(int);
   swap24(base, buff - base);

   if(om != NULL && om->nspan > 0) {
      base = buff;
      memcpy(buff, (void *)&om->row0, sizeof(int)); buff += sizeof(int);
      memcpy(buff, (void *)&om->col0, sizeof(int)); buff += sizeof(int);
      memcpy(buff, (void *)&om->rmin, sizeof(int)); buff += sizeof(int);
      memcpy(buff, (void *)&om->rmax, sizeof(int)); buff += sizeof(int);
      memcpy(buff, (void *)&om->cmin, sizeof(int)); buff += sizeof(int);
      memcpy(buff, (void *)&om->cmax, sizeof(int)); buff += sizeof(int);
      memcpy(buff, (void *)&om->npix, sizeof(int)); buff += sizeof(int);
      i = (om->data != NULL) ? 1 : 0;
      memcpy(buff, (void *)&i, sizeof(int)); buff += sizeof(int);
      
      swap24(base, buff - base);
/*
 * then the spans themselves
 */
      base = buff;
      for(i = 0; i < nspan; i++) {
	 memcpy(buff, (void *)&om->s[i], sizeof(SPAN)); buff += sizeof(SPAN);
      }
      
      swap2(base, buff - base);
/*
 * and then the data, if present
 */
      if(om->data != NULL) {
	 base = buff;
	 memcpy(buff, (void *)om->data, om->npix*sizeof(PIX));
	 buff += npix*sizeof(PIX); 
	 
	 swap2(base, buff - base);
      }
   }

   shAssert(buff <= buff0 + len);
   
   return(size);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Inflate an OBJMASK from a network-byte-order string of bytes.
 */
OBJMASK *
phObjmaskInflate(OBJMASK *om,		/* OBJMASK to fill, or NULL */
		 unsigned char *buff,	/* buffer to inflate from */
		 int *plen)		/* if non-NULL:
					   on input, length of buff or 0
					   on output, number of bytes read */
{
   unsigned char *const buff0 = buff;	/* initial value of buff */
   int i;
   int has_data = 0;			/* flattened OBJMASK had data != NULL*/
   int npix = 0, nspan = 0;		/* == om->n{pix,span} */
   
   shAssert(buff != NULL);
/*
 * first the ints, describing e.g. the size and bounding box
 */
   swap24(buff,sizeof(int));
   memcpy((void *)&nspan, buff, sizeof(int)); buff += sizeof(int);

   if(nspan == 0) {
      phObjmaskDel(om);
      om = phObjmaskNew(0);
   } else {
      if(om == NULL) {
	 om = phObjmaskNew(nspan); om->nspan = nspan;
      } else {
	 shAssert(om->size >= nspan);
      }
      
      swap24(buff,8*sizeof(int));
      memcpy((void *)&om->row0, buff, sizeof(int)); buff += sizeof(int);
      memcpy((void *)&om->col0, buff, sizeof(int)); buff += sizeof(int);
      memcpy((void *)&om->rmin, buff, sizeof(int)); buff += sizeof(int);
      memcpy((void *)&om->rmax, buff, sizeof(int)); buff += sizeof(int);
      memcpy((void *)&om->cmin, buff, sizeof(int)); buff += sizeof(int);
      memcpy((void *)&om->cmax, buff, sizeof(int)); buff += sizeof(int);
      memcpy((void *)&om->npix, buff, sizeof(int)); buff += sizeof(int);
      memcpy((void *)&has_data, buff, sizeof(int)); buff += sizeof(int);
/*
 * then the spans themselves
 */
      swap2(buff,nspan*3*sizeof(short));
      
      for(i = 0; i < nspan; i++) {
	 memcpy((void *)&om->s[i], buff, sizeof(SPAN)); buff += sizeof(SPAN);
      }
/*
 * and then the data, if present
 */
      if(has_data) {
	 npix = om->npix;
	 
	 swap2(buff,npix*sizeof(PIX));
	 
	 om->data = shRealloc(om->data, npix*sizeof(PIX));
	 memcpy((void *)om->data, buff, npix*sizeof(PIX));
	 buff += npix*sizeof(PIX);
      }
   }

   if(plen != NULL) {
      if(*plen > 0) {
	 shAssert(buff - buff0 <= *plen);
      }
      *plen = buff - buff0;
   }

   return(om);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Flatten/inflate an ATLAS_IMAGE into/from a network-byte-order string of
 * bytes.
 *
 * The buffer must be long enough; if it isn't an assertion will fail. You
 * can find out the needed size by first calling the routine with buff == NULL,
 */
int
phAtlasImageFlatten(const ATLAS_IMAGE *ai, /* ATLAS_IMAGE to flatten */
		    unsigned char *buff, /* buffer to flatten into; or NULL */
		    int len)		/* dimension of buff */
{
   unsigned char *base;			/* saved value of buff at misc times */
   unsigned char *const buff0 = buff;	/* initial value of buff */
   int compress_protocol;		/* how was data compressed? */
   int i;
   int ncolor, npix;			/* == ai->n{color,pix} */
   RICE_SEG *rseg = NULL;		/* Rice-coded pixel data */
   int size;				/* size required */

#if FLOATING_PHOTO
   compress_protocol = NOT_COMPRESSED;
#else
   compress_protocol = RICE_COMPRESSION_V2;
#endif
   
   shAssert(ai != NULL);
   ncolor = ai->ncolor; npix = ai->npix;
   size = 4*sizeof(int);

   if(ai->ncolor == 0) {
      shAssert(ai->master_mask == NULL);
   } else {
      shAssert(ai->master_mask != NULL);
      size += phObjmaskFlatten(ai->master_mask, NULL, 0);

      for(i = 0; i < ncolor; i++) {
	 size += 2*sizeof(int);
	 size += phSpanmaskFlatten(ai->regmask[i], NULL, 0);
	 size += phObjmaskFlatten(ai->mask[i], NULL, 0);
	 size += npix*sizeof(PIX);
	 if(compress_protocol != NOT_COMPRESSED) {
	    size += npix/RICE_J + 2 + 4*sizeof(int); /* max expansion */
	 }
      }
   }

   if(buff == NULL) {
      return(size);
   }
   shAssert(len >= size);
/*
 * first the ints, describing e.g. the number of colours present
 */
   base = buff;
   memcpy(buff, (void *)&ai->ncolor, sizeof(int)); buff += sizeof(int);
   memcpy(buff, (void *)&ai->id, sizeof(int)); buff += sizeof(int);
   memcpy(buff, (void *)&compress_protocol, sizeof(int)); buff += sizeof(int);
   if(compress_protocol >= RICE_COMPRESSION_V2) {
      memcpy(buff, (void *)&ai->npix, sizeof(int)); buff += sizeof(int);
   }

   swap24(base, buff - base);

   if(ai->ncolor == 0) {
      return(buff - buff0);
   }
/*
 * next the master_mask
 */
   buff += phObjmaskFlatten(ai->master_mask, buff, buff0 + len - buff);
/*
 * and then the data in each band
 */
   if(compress_protocol != NOT_COMPRESSED) {
      int nblock = npix/RICE_J;		/* number of Rice coded blocks */
      if(nblock*RICE_J != npix) nblock++;

      rseg = phRiceSegNew(nblock);
   }

   for(i = 0; i < ncolor; i++) {
      base = buff;
      memcpy(buff, (void *)&ai->drow[i], sizeof(int)); buff += sizeof(int);
      memcpy(buff, (void *)&ai->dcol[i], sizeof(int)); buff += sizeof(int);

      swap24(base, buff - base);

      buff += phSpanmaskFlatten(ai->regmask[i], buff, buff0 + len - buff);
      buff += phObjmaskFlatten(ai->mask[i], buff, buff0 + len - buff);
/*
 * compress data, if so desired
 */
      shAssert(ai->pix[i] != NULL);
      base = buff;
      if(compress_protocol != NOT_COMPRESSED) {
#if FLOATING_PHOTO
	 shFatal("You cannot get here: "
		 "(compress_protocol == %d) != NOT_COMPRESSED",
		 compress_protocol);
#else
	 phRiceSegCompress(ai->pix[i], npix, rseg, RICE_J);
	 buff += phRiceSegFlatten(rseg, buff, buff0 + len - buff);
#endif
      } else {
	 memcpy(buff, (void *)ai->pix[i], npix*sizeof(PIX));
	 buff += npix*sizeof(PIX); 
      
	 swap2(base, buff - base);
      }
   }
   phRiceSegDel(rseg);

   shAssert(buff <= buff0 + len && buff - buff0 <= size);
   
   return(buff - buff0);
}

/*
 * <AUTO EXTRACT>
 *
 * Inflate an ATLAS_IMAGE from a network-byte-order string of bytes.
 */
ATLAS_IMAGE *
phAtlasImageInflate(ATLAS_IMAGE *ai,	/* ATLAS_IMAGE to fill, or NULL */
		  unsigned char *buff,	/* buffer to inflate from */
		  int *plen)		/* if non-NULL:
					   on input, length of buff or 0
					   on output, number of bytes read */
{
   unsigned char *const buff0 = buff;	/* initial value of buff */
   int compress_protocol = 0;		/* how was data compressed? */
   int i;
   int len;				/* number of bytes read */
   int ncolor = 0, npix;		/* == ai->n{color,pix} */
   RICE_SEG *rseg = NULL;		/* Rice-coded pixel data */

   shAssert(buff != NULL);
/*
 * first the ints, describing e.g. the number of colours present
 */
   swap24(buff,3*sizeof(int));
   memcpy((void *)&ncolor, buff, sizeof(int)); buff += sizeof(int);

   if(ai == NULL) {
      ai = phAtlasImageNew(ncolor);
   } else {
      *(int *)&ai->ncolor = ncolor;
   }
   ai->shallow_copy = 0;

   memcpy((void *)&ai->id, buff, sizeof(int)); buff += sizeof(int);
   memcpy((void *)&compress_protocol, buff, sizeof(int)); buff += sizeof(int);
/*
 * Allow for changes in the protocol
 */
   if(compress_protocol >= RICE_COMPRESSION_V2) {
      swap24(buff,sizeof(int));
      memcpy((void *)&npix, buff, sizeof(int)); buff += sizeof(int);
   }

   if(ncolor == 0) {
      return(ai);
   }
/*
 * next the master_mask
 */
   len = (plen == NULL || *plen == 0) ? 0 : (buff0 + *plen) - buff;

   ai->master_mask = phObjmaskInflate(NULL, buff, &len); buff += len;
   if(compress_protocol < RICE_COMPRESSION_V2) {
      npix = ai->master_mask->npix;
   }
   ai->npix = npix;
/*
 * and then the data in each band
 */
   for(i = 0; i < ncolor; i++) {
      swap24(buff, 2*sizeof(int));
      memcpy((void *)&ai->drow[i], buff, sizeof(int)); buff += sizeof(int);
      memcpy((void *)&ai->dcol[i], buff, sizeof(int)); buff += sizeof(int);

      len = (plen == NULL || *plen == 0) ? 0 : (buff0 + *plen) - buff;
      ai->regmask[i] = phSpanmaskInflate(NULL, buff, &len); buff += len;

      len = (plen == NULL || *plen == 0) ? 0 : (buff0 + *plen) - buff;
      ai->mask[i] = phObjmaskInflate(NULL, buff, &len); buff += len;
/*
 * get pixel data, maybe uncompressing as we go
 */
      ai->pix[i] = shMalloc(npix*sizeof(PIX));
      if(compress_protocol == NOT_COMPRESSED) {
	 swap2(buff, npix*sizeof(PIX));
	 memcpy((void *)ai->pix[i], buff, npix*sizeof(PIX));
	 buff += npix*sizeof(PIX);
      } else {
#if FLOATING_PHOTO
	 shFatal("You cannot get here: "
		 "(compress_protocol == %d) != NOT_COMPRESSED",
		 compress_protocol);
#else
	 len = (plen == NULL || *plen == 0) ? 0 : (buff0 + *plen) - buff;
	 rseg = phRiceSegInflate(rseg, buff, &len); buff += len;
	 phRiceSegUncompress(ai->pix[i], rseg);
#endif
      }
   }
   phRiceSegDel(rseg);

   if(plen != NULL) {
      if(*plen > 0) {
	 shAssert(buff - buff0 <= *plen);
      }
      *plen = buff - buff0;
   }
   
   return(ai);
}

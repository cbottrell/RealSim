/*
 * Create ATLAS_IMAGEs from REGIONs and OBJC lists
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <alloca.h>
#include "dervish.h"
#include "phRandom.h"
#include "phMeasureObj.h"
#include "phOffset.h"

static PIX *ai_pix_set(PIX *pix, const OBJMASK *om, const REGION *reg,
		       int drow, int dcol);
static void ai_reg_set(const OBJMASK *om, const PIX *pix, REGION *reg,
		       int drow, int dcol, float sky);
#if !defined(STAND_ALONE)
static void ai_reg_set_val(const OBJMASK *om, REGION *reg,
			   int val, float sigma, RANDOM *rand,
			   int drow, int dcol);
#endif

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Create an ATLAS_IMAGE 
 */
ATLAS_IMAGE *
phAtlasImageNew(int ncolor)		/* desired number of colours */
{
   int c;
   ATLAS_IMAGE *new_ai = shMalloc(sizeof(ATLAS_IMAGE));

   shAssert(ncolor <= NCOLOR);
   *(int *)&new_ai->ncolor = ncolor;
   new_ai->id = new_ai->parent = -1;
   new_ai->run = new_ai->rerun = new_ai->camCol = new_ai->field = -1;
   new_ai->shallow_copy = 0;
   new_ai->master_mask = NULL;
   
   new_ai->npix = 0;
   for(c = 0;c < NCOLOR;c++) {
      new_ai->drow[c] = new_ai->dcol[c] = 0;
      new_ai->regmask[c] = NULL;
      new_ai->mask[c] = NULL;
      new_ai->pix[c] = NULL;
   }

   return(new_ai);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Make an atlas image from an OBJC
 */
ATLAS_IMAGE *
phAtlasImageNewFromObjc(const OBJC *objc)
{
   int c;
   ATLAS_IMAGE *new;
   OBJECT1 *obj1;

   shAssert(objc != NULL);

   new = phAtlasImageNew(objc->ncolor);
   new->shallow_copy = 1;

   new->id = objc->id;

   for(c = 0;c < objc->ncolor;c++) {
      obj1 = objc->color[c];
      shAssert(obj1 != NULL);

      if(obj1->region != NULL) {
	 new->regmask[c] = (SPANMASK *)obj1->region->mask;
	 if(new->regmask[c] != NULL) {
	    shAssert(new->regmask[c]->cookie == SPAN_COOKIE);
	 }
      }
      new->mask[c] = obj1->mask;
   }

   return(new);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Destroy an ATLAS_IMAGE. If deep is true, all structures allocated within the
 * ATLAS_IMAGE will be freed; if it's false they will not. The pix[] arrays
 * are deleted whether or not deep is true
 */
void
phAtlasImageDel(ATLAS_IMAGE *ai, int deep)
{
   int c;
   
   if(ai == NULL) return;

   if(p_shMemRefCntrGet(ai) > 0) {	/* still referenced somewhere */
      p_shMemRefCntrDecr(ai);
      return;
   }

   shAssert(!(deep && ai->shallow_copy));

   for(c = 0;c < ai->ncolor;c++) {
      shFree(ai->pix[c]);

      if(deep) {
	 phSpanmaskDel(ai->regmask[c]);
	 phObjmaskDel(ai->mask[c]);
      }
   }
   if(deep) {
      phObjmaskDel(ai->master_mask);
   }

   shFree(ai);
}

/***************************************************************************
 * <AUTO EXTRACT>
 *
 * delete the ATLAS_IMAGEs associated with the given OBJC and its family
 */
void
phAtlasImageDelFromObjc(OBJC *objc,	/* OBJC whose AIs are to be deleted */
			int deep)	/* should we destroy siblings
					   and children? */
{
   if(objc == NULL) return;
   
   if(deep) {
      phAtlasImageDelFromObjc(objc->children, 1); /* n.b. will recurse */
      phAtlasImageDelFromObjc(objc->sibbs, 1); /*      down lists */
   }

   phAtlasImageDel(objc->aimage, deep);
   objc->aimage = NULL;
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Return a copy of an ATLAS_IMAGE. If deep is true, all structures allocated
 * within the ATLAS_IMAGE will be copied; if it's false they will not (they'll
 * be set to point to the ones in the original).
 *
 * The dcol[], drow[], and pix[] arrays are copied whether or not deep is true
 */
ATLAS_IMAGE *
phAtlasImageCopy(const ATLAS_IMAGE *old, /* atlas image to copy */
		 int deep)		/* do a deep copy? */
{
   int c;
   ATLAS_IMAGE *ai;			/* atlas image to return */
   int ncolor;				/* number of colours in atlas images */

   shAssert(old != NULL);
   shAssert(old->master_mask != NULL && old->master_mask->npix >= 0);

   ncolor = old->ncolor;
   ai = phAtlasImageNew(ncolor);
   ai->shallow_copy = deep ? 0 : 1;

   ai->npix = old->npix;
   for(c = 0;c < ncolor;c++) {
      ai->drow[c] = old->drow[c];
      ai->dcol[c] = old->dcol[c];

      ai->pix[c] = shMalloc(old->npix*sizeof(PIX));
      memcpy(ai->pix[c], old->pix[c], old->npix*sizeof(PIX));
      
      if(deep) {
	 ai->regmask[c] = phSpanmaskCopy(old->regmask[c], 0, 0);
	 ai->mask[c] = phObjmaskCopy(old->mask[c], 0, 0);
      } else {
	 ai->regmask[c] = old->regmask[c];
	 ai->mask[c] = old->mask[c];
      }
   }

   if(deep) {
      ai->master_mask = phObjmaskCopy(old->master_mask, 0, 0);
   } else {
      ai->master_mask = old->master_mask;
   }
   

   return(ai);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Given an OBJC with its canonical centre set and a FIELDPARAMS,
 * set the pixel values pix[] in the OBJC's ATLAS_IMAGE
 *
 * If val is >= 0, the data regions in fparams will have their values
 * replaced by val after the atlas images are extracted. If sigma is
 * greater than zero, the values will have N(0,sigma^2) noise added
 * (why bother? because the deblender assumes that the sky is noisy
 * when correcting for biases)
 */
#if !defined(STAND_ALONE)
void
phAtlasImageCut(OBJC *objc,		/* set the ATLAS_IMAGE in here */
		int color,		/* which colour? (-1 => all) */
		const FIELDPARAMS *fparams, /* all about the frame */
		int val,		/* value to set */
		float sigma,		/* s.d. of val */
		RANDOM *rand)		/* random numbers; may be NULL if
					   sigma <= 0 */
{
   ATLAS_IMAGE *ai;			/* the atlas image to set */
   int c0, c1;				/* range of colours to cut */
   int c;
   const REGION *data;			/* == fparams->frame[c].data */
   float drow, dcol;			/* offsets from reference colour */

   shAssert(objc != NULL);
   shAssert(objc->aimage != NULL && objc->aimage->master_mask != NULL);
   shAssert(objc->ncolor == objc->aimage->ncolor);
   shAssert(color < 0 || color < objc->ncolor);
   shAssert(fparams != NULL);
   shAssert(objc->flags3 & OBJECT3_HAS_CENTER);
   shAssert(sigma <= 0 || rand != NULL);

   ai = objc->aimage;
   if(color >= 0) {
      c0 = c1 = color;
   } else {
      c0 = 0; c1 = objc->ncolor - 1;
   }

   ai->npix = ai->master_mask->npix;
   for(c = c0;c <= c1;c++) {
      data = fparams->frame[c].data;
      shAssert(data != NULL);
/*
 * calculate the offsets from the reference colour
 */
      phOffsetDo(fparams, data->nrow/2, data->ncol/2,
		 fparams->ref_band_index, c,
		 0, NULL, NULL, &drow, NULL, &dcol, NULL);
      ai->drow[c] = (drow > 0) ? drow + 0.5 : -(-drow + 0.5);
      ai->dcol[c] = (dcol > 0) ? dcol + 0.5 : -(-dcol + 0.5);
/*
 * copy over mask fields from objc->color, if they exist
 */
      if(objc->color[c] != NULL) {
	 ai->mask[c] = objc->color[c]->mask;
	 if(objc->color[c]->mask != NULL) {
	    objc->color[c]->mask->refcntr++;
	 }
      }
/*
 * and set the ATLAS_IMAGE->pix fields
 */
      shFree(ai->pix[c]);		/* it may be the wrong size */
      ai->pix[c] = ai_pix_set(NULL, ai->master_mask, data,
						     ai->drow[c], ai->dcol[c]);
/*
 * if so desired, set the cut-out pixels to some value
 */
      if(val >= 0) {
	 phRegionSetValFromAtlasImage(ai, c, (REGION *)fparams->frame[c].data,
				      val, sigma, rand, 0, 0);
	 
      }
   }
}
#endif

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Set the pix[] fields in an atlas image from a region
 */
void
phAtlasImageSetFromRegion(ATLAS_IMAGE *ai, /* the atlas image with pix to set*/
			  int c,	/* which colour? */
			  const REGION *data)	/* region with desired data */
{
   shAssert(ai != NULL && ai->master_mask != NULL);
   shAssert(ai->pix[0] == NULL || ai->npix == ai->master_mask->npix);
   shAssert(c >= 0 && c < ai->ncolor);
   shAssert(data != NULL);

   ai->npix = ai->master_mask->npix;
   ai->pix[c] = ai_pix_set(ai->pix[c], ai->master_mask, data,
						     ai->drow[c], ai->dcol[c]);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Set the pix[] fields in an atlas image to a certain value within a
 * specified OBJMASK
 */
void
phAtlasImageSetInObjmask(ATLAS_IMAGE *ai, /* the atlas image with pix to set*/
			 int c,		/* which colour? */
			 const OBJMASK *om, /* pixels to set */
			 const PIX val)	/* desired value */
{
   int drow, dcol;			/* how much to offset om's coords */
   int i, j, k;
   const OBJMASK *mm;			/* == ai->master_mask */
   SPAN *mp;				/* == mm->s[i] */
   PIX *ptr;				/* pointer to pix */
   SPAN *sp;				/* == om->s[i] */
   int x1, x2, y;			/* unpacked from sp */

   shAssert(c >= 0 && c < ai->ncolor);
   shAssert(ai != NULL && ai->master_mask != NULL && ai->pix[c] != NULL);
   shAssert(ai->master_mask->row0 == 0 && ai->master_mask->col0 == 0);
   shAssert(om != NULL && om->npix >= 0);
/*
 * Go through the objmask om looking for overlaps with the master_mask.
 *
 * Note that the master_mask's must be offset by (drow, dcol) to bring
 * into the c-band coordinate system, so we have to shift om the other way
 */
   drow = om->row0 - ai->drow[c]; dcol = om->col0 - ai->dcol[c];

   mm = ai->master_mask;
   ptr = ai->pix[c];
   j = 0;				/* counter in mm->s */
   mp = &mm->s[j]; j++;

   for(i = 0;i < om->nspan;i++) {
      sp = &om->s[i];
      y = sp->y; x1 = sp->x1; x2 = sp->x2;

      y += drow; x1 += dcol; x2 += dcol;
/*
 * advance mm's counter j until mm's current span is on same row as om's,
 * and om's span doesn't lie entirely to the left of mm's
 */
      if(y < mp->y) {			/* below mm's current span */
	 continue;
      } else {				/* above mm's current span */
	 while(mp->y < y || mp->x2 < x1) {
	    ptr += mp->x2 - mp->x1 + 1;
	    if(j == mm->nspan) {
	       j++;			/* signal that we ran out of spans */
	       break;
	    }
	    mp = &mm->s[j]; j++;
	 }

	 if(j > mm->nspan) {		/* we ran out of spans in mm */
	    break;
	 }
      }

      if(mp->y > y || mp->x1 > x2) {	/* no overlap in this row */
	 continue;
      }
/*
 * If we get here, there's an overlap. Set the requisite values in pix
 */
      shAssert(y == mp->y && x1 <= mp->x2 && x2 >= mp->x1);
      if(x1 < mp->x1) {
	 x1 = mp->x1;
      }
      if(x2 > mp->x2) {
	 x2 = mp->x2;
      }
      for(k = x1; k <= x2; k++) {
	 ptr[k - mp->x1] = val;
      }
   }
#if !defined(NDEBUG)
   if(j <= mm->nspan) {			/* we haven't finished with mm */
      for(;;) {
	 ptr += mp->x2 - mp->x1 + 1;
	 if(j == mm->nspan) {
	    break;
	 }
	 mp = &mm->s[j]; j++;
      }
   }
   shAssert(ptr == ai->pix[c] + ai->npix);
#endif
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Given an ATLAS_IMAGE containing pixel data, set the corresponding pixels
 * in a region
 *
 * Because we may want to restore the atlas image data into a smaller
 * region than the one that it was created from, the region's physical
 * origin is taken to be at [row0,col0] in the atlas image's coordinate
 * system.
 *
 * Note that the atlas image is also shifted by the appropriate drow, dcol
 * to take it to the proper coordinate system for its colour
 */
void
phRegionSetFromAtlasImage(const ATLAS_IMAGE *ai, /* here's the data */
			  int c,	/* colour to set */
			  REGION *reg,	/* the region in question */
			  int row0, int col0, /* offset ai's origin to here */
			  float sky)	/* add this to all ai's pixels */
{
   shAssert(ai != NULL);
   shAssert(c >= 0 && c < ai->ncolor);
   shAssert(ai->master_mask != NULL && ai->master_mask->npix >= 0);

   if(ai->npix == 0) {			/* don't set region */
      return;
   }
   
#if !defined(NDEBUG)
   if(ai->master_mask->npix == 0) {
      shAssert(ai->pix[c] == NULL);
   } else {
      shAssert(ai->pix[c] != NULL && ai->npix == ai->master_mask->npix);
   }
#endif
   shAssert(reg != NULL && reg->type == TYPE_PIX);

   ai_reg_set(ai->master_mask, ai->pix[c], reg,
	      ai->drow[c] - row0, ai->dcol[c] - col0, sky);
}

/*****************************************************************************/
/*
 * Set the value of a part of a REGION wherever an atlas image has been cut.
 *
 * If sigma is greater than zero, the values will have N(0,sigma^2) noise
 * added (why bother? because the deblender assumes that the sky is noisy
 * when correcting for biases)
 */
#if !defined(STAND_ALONE)
void
phRegionSetValFromAtlasImage(const ATLAS_IMAGE *ai, /* atlas image */
			     int c,	/* which band? */
			     REGION *reg, /* region to set */
			     int val,	/* value to set */
			     float sigma, /* s.d. of val */
			     RANDOM *rand, /* random numbers; may be NULL if
					      sigma <= 0 */
			     int row0,	/* offset of ai in reg */
			     int col0)
{
   shAssert(ai != NULL);
   shAssert(c >= 0 && c < ai->ncolor);
   shAssert(ai->master_mask != NULL && ai->master_mask->npix >= 0);
   shAssert(reg != NULL && reg->type == TYPE_PIX);
   if(sigma > 0 && rand == NULL) {
      sigma = 0;
   }

   ai_reg_set_val(ai->master_mask, reg, val, sigma, rand,
		  ai->drow[c] - row0, ai->dcol[c] - col0);
}
#endif

/*****************************************************************************/
/*
 * Return an array of all the pixel values corresponding to pixels present
 * in the OBJMASK, the values being those in the REGION.
 *
 * We allow the origin of the region to be different from that in the mask;
 * the convention is that pixel [drow, dcol] in the region is [0,0] in the mask
 */
static PIX *
ai_pix_set(PIX *pix,			/* space to hold the pixels, or NULL */
	   const OBJMASK *om,		/* the OBJMASK specifying the pixels */
	   const REGION *reg,		/* the region wherein they dwell */
	   int drow, int dcol)		/* how much to offset om's coords */
{
   int i;
   int npix;				/* number of pixels in OBJMASK */
   PIX **rows;				/* == reg->ROWS */
   PIX *rptr;				/* == &reg->ROWS[y][x] */
   PIX *ptr;				/* pointer to pix */
   SPAN *sp;				/* == om->s[i] */
   int nrow, ncol;			/* == reg->n{row,col} */
   int x1, x2, y;			/* unpacked from sp */

   shAssert(om != NULL && om->npix >= 0);
   shAssert(reg != NULL && reg->type == TYPE_PIX);
/*
 * allocate pix if needs be
 */
   if(pix == NULL) {
      pix = shMalloc(om->npix*sizeof(PIX));
   }
/*
 * Set pix from the region. We must be careful as the
 * offset OBJMASK may lie outside the REGION
 */
   drow -= reg->row0; dcol -= reg->col0;

   rows = reg->ROWS;
   nrow = reg->nrow; ncol = reg->ncol;
   ptr = pix;

   for(i = 0;i < om->nspan;i++) {
      sp = &om->s[i];
      y = sp->y; x1 = sp->x1; x2 = sp->x2;
      npix = x2 - x1 + 1;

      y += drow; x1 += dcol; x2 += dcol;
      if(y >= 0 && y < nrow) {
	 int ncopy = npix;		/* number of pixels to copy */
	 int col0 = x1;			/* starting column in region */
	 int ptr0 = 0;			/* starting index in ptr */

	 if(x1 < 0) {			/* span starts off left of region */
	    col0 = 0;
	    ptr0 = -x1;
	    if(x2 < 0) {
	       ptr0 = npix;
	    }

	    memset(ptr,'\0',ptr0*sizeof(PIX));
	    ncopy -= ptr0;
	 }
	 if(x2 >= ncol) {		/* span extends to right of region */
	    x2 = x2 - ncol;
	    
	    ncopy -= x2;
	    if(ncopy > 0) {		/* x1 is < ncol */
	       memset(&ptr[ptr0 + ncopy],'\0',x2*sizeof(PIX));
	    }
	 }
	 
	 if(ncopy > 0) {
	    rptr = &rows[y][col0];
	    memcpy(&ptr[ptr0],rptr,ncopy*sizeof(PIX));
	 }
      } else {
	 memset(ptr,'\0',npix*sizeof(PIX));
      }
      ptr += npix;
   }
   shAssert(ptr == pix + om->npix);

   return(pix);
}

/*****************************************************************************/
/*
 * Return an array of all the pixel values corresponding to pixels present
 * in the OBJMASK, the values being those in the REGION.
 *
 * We allow the origin of the region to be different from that in the mask;
 * the convention is that pixel [drow, dcol] in the region is [0,0] in the mask
 */
static void
ai_reg_set(const OBJMASK *om,		/* the OBJMASK specifying the pixels */
	   const PIX *pix,		/* the values to use */
	   REGION *reg,			/* the region wherein they dwell */
	   int drow, int dcol,		/* how much to offset om's coords */
	   float sky)			/* add this value to all pixels */
{
   int i;
   int npix;				/* number of pixels in OBJMASK */
   PIX **rows;				/* == reg->ROWS */
   const PIX *ptr;			/* pointer to pix */
   SPAN *sp;				/* == om->s[i] */
   int nrow, ncol;			/* == reg->n{row,col} */
   int x1, x2, y;			/* unpacked from sp */
/*
 * Set the pixels in the region.
 *
 * We must be careful as the offset OBJMASK may lie outside the REGION
 */
   ptr = pix;
   rows = reg->ROWS;
   nrow = reg->nrow; ncol = reg->ncol;

   for(i = 0;i < om->nspan;i++) {
      sp = &om->s[i];
      y = sp->y; x1 = sp->x1; x2 = sp->x2;
      npix = x2 - x1 + 1;

      y += drow; x1 += dcol; x2 += dcol;
      if(y >= 0 && y < nrow) {
	 int col0 = x1;			/* starting column in image */
	 int ncopy = npix;		/* number of pixels to copy */
	 int ptr0 = 0;			/* starting index in ptr[] */
	    
	 if(x1 < 0) {
	    col0 = 0;
	    ptr0 = -x1;
	    
	    ncopy -= ptr0;
	 }
	 if(x2 >= ncol) {
	    ncopy -= (x2 - ncol + 1);
	 }

	 if(ncopy > 0) {
	    if(sky == 0) {
	       memcpy(&rows[y][col0],&ptr[ptr0],ncopy*sizeof(PIX));
	    } else {
	       int j;
	       for(j = 0; j < ncopy; j++) {
		  rows[y][col0 + j] = ptr[ptr0 + j] + sky;
	       }
	    }
	 }
      }
      ptr += npix;
   }
   shAssert(ptr == pix + om->npix);
}

/*****************************************************************************/
/*
 * Set all the pixels in the OBJMASK to val
 *
 * We allow the origin of the region to be different from that in the mask;
 * the convention is that pixel [drow, dcol] in the region is [0,0] in the mask
 */
#if !defined(STAND_ALONE)
static void
ai_reg_set_val(const OBJMASK *om,	/* the OBJMASK specifying the pixels */
	       REGION *reg,		/* the region to set */
	       int val,			/* the desired region */
	       float sigma,		/* s.d. of val */
	       RANDOM *rand,		/* random numbers; may be NULL if
					   sigma <= 0 */
	       int drow, int dcol)	/* how much to offset om's coords */
{
   float fac,r = 0,v1,v2;		/* for Gaussian noise; initialise r to
					   appease the IRIX 5.3 cc */
   int i, j;
   const float inorm = 1.0/(float)((1U<<(8*sizeof(int)-1)) - 1);
   int npix;				/* number of pixels in OBJMASK */
   PIX **rows;				/* == reg->ROWS */
   PIX *rptr;				/* == &reg->ROWS[y][x] */
   SPAN *sp;				/* == om->s[i] */
   int nrow, ncol;			/* == reg->n{row,col} */
   int x1, x2, y;			/* unpacked from sp */
   float var2 = 0;			/* twice the variance corrected
					   for dither noise */
   shAssert(om->npix >= 0);
   shAssert(reg != NULL && reg->type == TYPE_PIX);
/*
 * set up the random number generator
 */
   if(sigma > 0) {
      var2 = 2*(sigma*sigma - 1/12.0);	/* allow for dither noise */

      if(var2 <= 0) {
	 sigma = -1;			/* desired sigma's less than dither
					   noise */
      }
      shAssert(rand != NULL);
   }
/*
 * Set pixels to val; we must be careful as the offset OBJMASK may lie
 * outside the REGION
 */
   rows = reg->ROWS;
   nrow = reg->nrow; ncol = reg->ncol;

   DECLARE_PHRANDOM(rand);

   for(i = 0;i < om->nspan;i++) {
      sp = &om->s[i];
      y = sp->y; x1 = sp->x1; x2 = sp->x2;
      npix = x2 - x1 + 1;

      y += drow; x1 += dcol; x2 += dcol;
      if(y >= 0 && y < nrow) {
	 int ncopy = npix;		/* number of pixels to copy */
	 int col0 = x1;			/* starting column in region */

	 if(x1 < 0) {			/* span starts off left of region */
	    col0 = 0;
	    ncopy -= (-x1);
	 }
	 if(x2 >= ncol) {		/* span extends to right of region */
	    x2 = x2 - ncol;
	    
	    ncopy -= x2;
	 }
	 
	 if(ncopy > 0) {
	    rptr = &rows[y][col0];
	    
	    if(sigma > 0) {
	       j = 0;
	       if(ncopy & 01) {
		  rptr[j++] = val +
		    sqrt(var2/2)*phGaussdev() + (PHRANDOM & 0x1);
	       }
	       
	       while(j < ncopy) {
		  do {
		     v1 = PHRANDOM*inorm;
		     v2 = PHRANDOM*inorm;
		     r = v1*v1+v2*v2;
		  } while (r >= 1.0);
		  fac = sqrt(-var2*log(r)/r);
		  
		  rptr[j++] = val + fac*v1 + (PHRANDOM & 0x1);
		  rptr[j++] = val + fac*v2 + (PHRANDOM & 0x1);
	       }
	    } else {
	       for(j = 0;j < ncopy;j++) {
		  rptr[j] = val;
	       }
	    }
	 }
      }
   }

   END_PHRANDOM(rand);
}
#endif

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Set all the pixels in an atlas image that are not an objmask to val
 *
 * The objmask is supposed to be in c's coordinate system
 */
void
phAtlasImageSetIfNotInObjmask(ATLAS_IMAGE *ai, /* the atlas image */
			      int c,	/* use this band */
			      const OBJMASK *omask, /* not these pixels */
			      const PIX val) /* value to use */
{
   int ai_drow, ai_dcol;		/* == ai->{row,dcol}[c] */
   int i;
   int i0 = -1;				/* start of settable section of pptr */
   int mmi, omi;			/* counters in omask and mmask */
   int omi0;				/* starting value of omi */
   const OBJMASK *mmask;		/* == ai->master_mask */
   int mx1, mx2, my;			/* unpacked from mmask->s[] */
   PIX *pptr;				/* pointer to ai->pix[c] */
   int x1, x2, y;			/* unpacked from omask->s[] */

   shAssert(ai != NULL && ai->master_mask != NULL);
   shAssert(c >= 0 && c < ai->ncolor && ai->pix[c] != NULL);
   
   mmask = ai->master_mask;
   ai_drow = ai->drow[c];
   ai_dcol = ai->dcol[c];

   if(mmask->nspan == 0) {
      return;
   }
   
   pptr = ai->pix[c];
   mx1 = mx2 = 0;
   for(omi0 = mmi = 0; mmi < mmask->nspan; mmi++) {
      if(mmi > 0) {
	 pptr += mx2 - mx1 + 1;		/* skip over previous span */
      }

      my = mmask->s[mmi].y + ai_drow;
      mx1 = mmask->s[mmi].x1 + ai_dcol;
      mx2 = mmask->s[mmi].x2 + ai_dcol;

      i0 = mx1;
/*
 * skip all spans in omask below and to the left of this mmask span.
 * We may need to process several spans in omask, and they may be on
 * the same row as the next span in mmask --- hence the need to keep
 * the starting value of omi, omi0
 */
      if(omi0 >= omask->nspan) {	/* finished omask */
	 break;
      }
      omi = omi0;
      y = omask->s[omi0].y; x1 = omask->s[omi0].x1; x2 = omask->s[omi0].x2;
      while(y < my || (my == y && x2 < mx1)) {
	 if(++omi == omask->nspan) {
	    y = omask->rmax + 1;
	    break;
	 }
	 y = omask->s[omi].y; x1 = omask->s[omi].x1; x2 = omask->s[omi].x2;
      }
      if(y > omask->rmax) {		/* finished omask */
	 break;
      }
      omi0 = omi;
/*
 * process all omask spans that intersect this mmask span
 */
      while(my == y && x1 <= mx2) {
	 if(x1 < mx1) x1 = mx1;
	 
	 for(i = i0;i < x1;i++) {
	    pptr[i - mx1] = val;
	 }
	 i0 = x2 + 1;

	 if(++omi == omask->nspan) {
	    y = omask->rmax + 1;
	    break;
	 }
	 y = omask->s[omi].y; x1 = omask->s[omi].x1; x2 = omask->s[omi].x2;
      }
      if(y > omask->rmax) {		/* finished omask */
	 break;
      }

      for(i = i0;i <= mx2;i++) {
	 pptr[i - mx1] = val;
      }
   }
/*
 * now set the rest of mmask; both the rest of this span and all of each of
 * the rest of the spans
 */
   if(i0 > mx2) i0 = mx2;		/* not beyond the end of the span */
   i0 = &pptr[i0 - mx1] - ai->pix[c];	/* relative to start of pix array */
   shAssert(i0 >= 0);

   pptr = ai->pix[c];
   for(i = i0; i < ai->npix; i++) {
      pptr[i] = val;
   }
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Add the pixels in one atlas image to another:
 *    "ai->pix[c] += cai->pix[c] + bias;"
 *
 * Both atlas images must be the same size.
 * If c is in range, just add that band; if it is -1, add all bands
 *
 * See also phAtlasImagesMinusEquals()
 */
void
phAtlasImagesPlusEquals(ATLAS_IMAGE *ai, /* atlas image to add to */
			const ATLAS_IMAGE *cai,	/* atlas image to add */
			int bias,	/* value to add to cai->pix[] */
			int c)		/* which color to use, or -1 */
{
   PIX *ai_pix;				/* == ai->pix[c] */
   const PIX *cai_pix;			/* == cai->pix[c] */
   int c0, c1;				/* range of colours to add */
   int i;
   int npix;				/* == ai->npix */

   shAssert(ai != NULL && cai != NULL);
   shAssert(ai->ncolor == cai->ncolor && ai->npix == cai->npix);
   shAssert(c >= -1 && c < ai->ncolor);
   npix = ai->npix;

   if(c == -1) {
      c0 = 0; c1 = ai->ncolor;
   } else {
      c0 = c; c1 = c + 1;
   }

   for(c = c0; c < c1; c++) {
      ai_pix = ai->pix[c];
      cai_pix = cai->pix[c];
      if(bias == 0) {
	 for(i = 0; i < npix; i++) {
	    ai_pix[i] += cai_pix[i];
	 }
      } else {
	 for(i = 0; i < npix; i++) {
	    ai_pix[i] += cai_pix[i] + bias;
	 }
      }
   }
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Subtract the pixels in one atlas image from another:
 *    "ai->pix[c] -= cai->pix[c] + bias;"
 *
 * Both atlas images must be the same size.
 * If c is in range, just add that band; if it is -1, add all bands
 *
 * See also phAtlasImagesPlusEquals()
 */
void
phAtlasImagesMinusEquals(ATLAS_IMAGE *ai, /* atlas image to add to */
			 const ATLAS_IMAGE *cai, /* atlas image to add */
			 int bias,	/* value to add to cai->pix[] */
			 int c)		/* which color to use, or -1 */
{
   PIX *ai_pix;				/* == ai->pix[c] */
   const PIX *cai_pix;			/* == cai->pix[c] */
   int c0, c1;				/* range of colours to add */
   int i;
   int npix;				/* == ai->npix */

   shAssert(ai != NULL && cai != NULL);
   shAssert(ai->ncolor == cai->ncolor && ai->npix == cai->npix);
   shAssert(c >= -1 && c < ai->ncolor);
   npix = ai->npix;

   if(c == -1) {
      c0 = 0; c1 = ai->ncolor;
   } else {
      c0 = c; c1 = c + 1;
   }

   for(c = c0; c < c1; c++) {
      ai_pix = ai->pix[c];
      cai_pix = cai->pix[c];
      if(bias == 0) {
	 for(i = 0; i < npix; i++) {
	    ai_pix[i] -= cai_pix[i];
	 }
      } else {
	 for(i = 0; i < npix; i++) {
	    ai_pix[i] -= cai_pix[i] + bias;
	 }
      }
   }
}

#if !defined(STAND_ALONE)
/*****************************************************************************/
/*
 * Consolidate the good pixels in an ATLAS_IMAGE's pix[] arrays; the
 * number of good pixels is returned
 */
struct sect {
   int i0, i1;				/* starting/ending indices */
};					/* bad sections of pix[] */

static void
consolidate_pix(ATLAS_IMAGE *aimage,	/* the atlas image in question */
		struct sect *bad,	/* bad sections of pix[] */
		int nbad)		/* number of bad sections */
{
   int c;				/* colour counter */
   int idst;				/* index of destination section */
   struct sect *good;			/* the good sections */
   int ncolor = aimage->ncolor;		/* number of bands */
   int ngood;				/* number of good sections */
   const int npix = aimage->npix;	/* number of pixels in pix[] */
   int npix_good;			/* number of good pixels */
   int i;
   int n;				/* number of pixels to copy */
   PIX **pix = aimage->pix;		/* data arrays */

   if(nbad == 0 || npix == 0) {
      return;
   }
/*
 * find the _good_ sections, merging contiguous bad ones
 */
   good = alloca((nbad + 1)*sizeof(struct sect));
   ngood = 0;

   i = 0;
   if(bad[0].i0 == 0) {
      good[ngood].i0 = bad[i].i1 + 1;
      i++;
   } else {
      good[ngood].i0 = 0;
   }

   for(; i < nbad; i++) {
      good[ngood].i1 = bad[i].i0 - 1;
      if(good[ngood].i1 >= good[ngood].i0) { /* not a null good section */
	 ngood++;
      }
      good[ngood].i0 = bad[i].i1 + 1;
   }

   if(good[ngood].i0 != npix) {
      good[ngood++].i1 = npix - 1;
   }
/*
 * Copy the good sections to make them contiguous.
 *
 * We could realloc, but shRealloc is too dumb to make this worth doing;
 * the memory involved is (probably) tiny anyway
 */
   shAssert(ngood <= nbad + 1);
   idst = npix_good = 0;
   for(i = 0; i < ngood; i++) {
      int isrc = good[i].i0;
      n = good[i].i1 - good[i].i0 + 1;
      npix_good += n;
      shAssert(idst + n <= npix && isrc + n <= npix);
      if(isrc != idst) {
	 for(c = 0; c < ncolor; c++) {
	    memmove(&pix[c][idst], &pix[c][isrc], n*sizeof(PIX));
	 }
      }
      idst += n;
   }

   aimage->npix = npix_good;
}

/*
 * Trim an OBJC's atlas image to lie within the specified rectangle
 */
void
phAtlasImageTrimToRect(OBJC *objc,	/* the object in question */
		       int rmin, int cmin, /* LLC of rectangle */
		       int rmax, int cmax) /* URC of rectangle */
{
   ATLAS_IMAGE *aimage;			/* == objc->aimage */
   struct sect *bad_sections;		/* bad sections of pix[] */
   int ipix;				/* current index into pix[] */
   int i;
   OBJMASK *mmask;			/* == objc->aimage->master_mask */
   int nbad;				/* number of bad sections */
   int npix_span;			/* number of pixels in a span */
   SPAN *s;				/* a span in the mmask */

   shAssert(objc != NULL && objc->aimage != NULL);
   aimage = objc->aimage;
   mmask = aimage->master_mask;
   shAssert(mmask != NULL);
/*
 * Allocate the book-keeping for bad sections
 */
   bad_sections = alloca(mmask->nspan*sizeof(struct sect));
   nbad = 0;

   ipix = 0;				/* where we've got to aimage->pix */
   for(i = 0; i < mmask->nspan; i++, ipix += npix_span) {
      s = &mmask->s[i];
      npix_span = s->x2 - s->x1 + 1;

      if(s->y < rmin || s->y > rmax || s->x2 < cmin || s->x1 > cmax) {
	 bad_sections[nbad].i0 = ipix;
	 bad_sections[nbad].i1 = ipix + npix_span - 1;
	 nbad++;
	 
	 s->x2 = s->x1 - 1;		/* delete entire span */
	 continue;
      }
/*
 * we now know that at least one pixel in the middle of the span's good
 */
      if(s->x1 < cmin) {
	 bad_sections[nbad].i0 = ipix;
	 bad_sections[nbad].i1 = ipix + (cmin - 1) - s->x1;
	 nbad++;

	 s->x1 = cmin;
	 continue;
      }

      if(s->x2 > cmax) {
	 bad_sections[nbad].i0 = ipix + (cmax + 1) - s->x1;
	 bad_sections[nbad].i1 = ipix + s->x2 - s->x1;
	 nbad++;
	 
	 s->x2 = cmax;
	 continue;
      }
/*
 * An entirely good span; that's easy
 */
      ;
   }
/*
 * Actually move the pixels around.
 */
   shAssert(nbad <= mmask->nspan);
   consolidate_pix(aimage, bad_sections, nbad);

   phCanonizeObjmask(mmask, 1);
   shAssert(aimage->pix[0] == NULL || aimage->npix == mmask->npix);
}

/*****************************************************************************/
/*
 * Return the index in the ->pix arrays of the pixel (row, col)
 */
int
phFindAtlasImagePixel(const ATLAS_IMAGE *aimage, /* the atlas image */
			int row, int col) /* the desired pixel */
{
   int i;
   const OBJMASK *mmask;		/* == aimage->master_mask */
   int n;				/* number of pixels before (row,col) */
   
   if(aimage == NULL || aimage->master_mask == NULL) {
      shError("atlas image or master mask is NULL");
      return(-1);
   }
   mmask = aimage->master_mask;

   for(i = n = 0; i < mmask->npix; i++) {
      if(mmask->s[i].y < row || mmask->s[i].x1 > col || mmask->s[i].x2 < col) {
	 n += mmask->s[i].x2 - mmask->s[i].x1 + 1;
	 continue;
      } else if(mmask->s[i].y > row) {
	 shError("Pixel (%d, %d) isn't in the atlas image", row, col);
	 return(-1);
      }

      if(!(mmask->s[i].y == row &&
	   mmask->s[i].x1 <= col && mmask->s[i].x2 >= col)) {
	 shError("Bug: you cannot get here");
	 return(-1);
      }

      n += col - mmask->s[i].x1;

      return(n);
   }

   shError("Pixel (%d, %d) isn't in the atlas image", row, col);
   return(-1);
}
#endif

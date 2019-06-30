/*
 * <INTRO>
 * Utilities to handle masks made from SPANs
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <alloca.h>
#include "dervish.h"
#include "phConsts.h"
#include "phSpanUtil.h"
#if !defined(STAND_ALONE)
#  include "phUtils.h"
#endif

#ifndef MAXSHORT
#  define MAXSHORT MAX_S16
#endif

#if !defined(STAND_ALONE)
/*
 * Because SPAN.y and SPAN.x1 are both less than 16 bits, we can generate
 * a 32-bit sort key SORTKEY(x,y) by putting y in the upper 16 bits, and
 * x1 in the lower; if the spans are sorted on this they are sorted first
 * by y and then by x1. This fails if x < 0, so add an x-bias of 10000
 */
#define SORTKEY(X,Y) (((Y) << 16) | (10000 + X))

static void heap_sort(SPAN *s, unsigned int n);
static void insertion_sort(SPAN *s, unsigned int n);
static int is_canonical_objmaskChain(const CHAIN *chain);
#endif

static TYPE objmask_type = UNKNOWN;	/* == shTypeGetFromName("OBJMASK") */

/*****************************************************************************/
/*
 * MUST be called before calling any SPAN/OBJMASK routines
 */
void
initSpanObjmask(void)
{
   objmask_type = shTypeGetFromName("OBJMASK");
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Make a new OBJMASK
 */
OBJMASK *
phObjmaskNew(
	    int nspans			/* number of spans */
	    )
{
   static int id = 0;
   OBJMASK *new;
   
   new = shMalloc(sizeof(OBJMASK));
   *(int *)&new->id = id++;
   new->refcntr = 1;
   new->s = (nspans == 0) ? NULL : shMalloc(nspans*sizeof(SPAN));
   new->size = nspans;
   new->nspan = 0;
   new->row0 = new->col0 = 0;
   new->rmin = new->rmax = 0;
   new->cmin = new->cmax = 0;
   new->npix = -1; new->sum = 0.0;
   new->data = NULL;
   new->user = NULL;

   return(new);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Grab more memory for a OBJMASK
 *
 * Return sv
 */
OBJMASK *
phObjmaskRealloc(OBJMASK *sv,		/* Old Spanvec; may be NULL */
		 int nspans)		/* new number of spans */
{
   if(sv == NULL) {
      sv = phObjmaskNew(nspans);
   } else {
      sv->s = shRealloc(sv->s, nspans*sizeof(SPAN));
      sv->size = nspans;
      if(sv->nspan > nspans) {
	 sv->nspan = nspans;
      }
   }

   return(sv);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Delete a OBJMASK
 */
void
phObjmaskDel(OBJMASK *sv)
{
   if(sv == NULL) return;

   if(--sv->refcntr > 0) {		/* still referenced somewhere */
      return;
   }
   if(p_shMemRefCntrGet(sv) > 0) {	/* still referenced somewhere */
      p_shMemRefCntrDecr(sv);
      return;
   }
   
   shFree(sv->s);
   shFree(sv->data);
   shFree(sv);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Delete a OBJMASK chain
 */
void
phObjmaskChainDel(CHAIN *mask)
{
   OBJMASK *sv;

   if(mask == NULL) return;

   shAssert(mask->type == objmask_type);
   
   while((sv = shChainElementRemByPos(mask, 0)) != NULL) {
      phObjmaskDel(sv);
   }
   shChainDel(mask);
}    

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Make a copy of a OBJMASK (doesn't copy OBJMASK->data)
 */
OBJMASK *
phObjmaskCopy(const OBJMASK *sv,	/* OBJMASK to be copied */
	      float fdr, float fdc	/* shift by this much */
	      )
{
   OBJMASK *new;
   int dr, dc;

   if(sv == NULL) {
      return(NULL);
   }
   new = phObjmaskNew(sv->nspan);

   dr = (fdr < 0) ?  -(-fdr + 0.5) : fdr + 0.5;
   dc = (fdc < 0) ?  -(-fdc + 0.5) : fdc + 0.5;

   new->nspan = sv->nspan;
   new->row0 = sv->row0;
   new->col0 = sv->col0;
   new->npix = sv->npix;
   
   if(dr == 0 && dc == 0) {
      memcpy(new->s,sv->s,sv->nspan*sizeof(SPAN));
      new->cmin = sv->cmin;
      new->rmin = sv->rmin;
      new->cmax = sv->cmax;
      new->rmax = sv->rmax;
   } else {
      int i, n = sv->nspan;
      SPAN *const snew = new->s, *s = sv->s;

      for(i = 0;i < n;i++) {
	 snew[i].y = s[i].y + dr;
	 snew[i].x1 = s[i].x1 + dc;
	 snew[i].x2 = s[i].x2 + dc;
      }
      new->cmin = sv->cmin + dc;
      new->rmin = sv->rmin + dr;
      new->cmax = sv->cmax + dc;
      new->rmax = sv->rmax + dr;
   }

   return new;
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Make a copy of a complete OBJMASK CHAIN
 */
CHAIN *
phObjmaskChainCopy(const CHAIN *chain,	/* chain to copy */
		   float dr, float dc)	/* shifted by this amount */
{
    int i;
    int len;				/* length of chain */
    CHAIN *om;
    char *cobjmask_type = (char *)objmask_type;
    OBJMASK *obj, *nobj;

    shAssert(chain != NULL && chain->type == objmask_type);

    om = shChainNew(cobjmask_type);
    len = chain->nElements;
    for(i = 0; i < len; i++) {
       obj = shChainElementGetByPos(chain, i);
       nobj = phObjmaskCopy(obj, dr, dc);
       shChainElementAddByPos(om, nobj, cobjmask_type, TAIL, AFTER);
    }

    return(om);
}    

#if !defined(STAND_ALONE)

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Grow an OBJMASK out by <n> pixels in all directions.
 *
 * Note that this grow routine isn't very good; it grows correctly along
 * the rows, and grows the very topmost and bottommost spans properly in
 * the row direction, but it doesn't do a full isotropic grow.
 *
 * This may or may not matter; let us leave it for now (it's fine for objects
 * with only a few pixels).
 */
OBJMASK *
phObjmaskGrow(const OBJMASK *om,	/* OBJMASK to grow */
	      const REGION *reg,	/* REGION that om lives in */
	      int n)			/* by how many pixels to grow */
{
   int i,j,k;
   int nrmin, nrmax;			/* number of spans at the top and
					   bottom of obj respectively */
   int ngmin, ngmax;			/* number of spans to grow obj by
					   for each top or bottom input span */
   int nspan;				/* == om->nspan */
   OBJMASK *onew;			/* == onew->sv[0] */
   int r_rmin, r_rmax, r_cmin, r_cmax;	/* bounding box for reg */
   int x1, x2, y;

   shAssert(om != NULL && om->nspan > 0);
   shAssert(phObjmaskIsCanonical(om));
   shAssert(n >= 0);
   shAssert(reg != NULL);

   nspan = om->nspan;
   r_rmin = reg->row0;
   r_rmax = r_rmin + reg->nrow - 1;
   r_cmin = reg->col0;
   r_cmax = r_cmin + reg->ncol - 1;
   shAssert(om->rmin >= r_rmin && om->rmax <= r_rmax);
/*
 * we expect to have to add n levels at the top and bottom of the mask,
 * so let us face up to the memcpy involved. This is one reason why this
 * routine doesn't grow the OBJECT in place -- it wouldn't be appreciably more
 * efficient
 *
 * See how many new SPANs we'll need (if merging occurs this is an
 * overestimate, but not by much).
 */
   y = om->rmin;
   for(i = 0;i < nspan;i++) {
      if(om->s[i].y != y) {
	 break;
      }
   }
   nrmin = i;				/* number of spans with y == rmin */
   ngmin = (y - r_rmin > n ? n : y - r_rmin);

   y = om->rmax;
   for(i = nspan - 1;i >= 0;i--) {
      if(om->s[i].y != y) {
	 break;
      }
   }
   nrmax = nspan - 1 - i;		/* number of spans with y == rmax */
   ngmax = (r_rmax - y > n) ? n : r_rmax - y; /* number of spans to add */
   shAssert(ngmin >= 0 && ngmax >= 0);
/*
 * OK, time to make the output mask
 */
   onew = phObjmaskNew(nspan + nrmin*ngmin + nrmax*ngmax);
   onew->row0 = om->row0; onew->col0 = om->col0;
/*
 * Grow the objmask down below rmin
 */
   y = om->rmin;
   for(i = j = 0;i < nrmin;i++) {
      for(k = ngmin;k > 0;k--) {
	 x1 = om->s[i].x1 - (n - k + 1);
	 if(x1 < r_cmin) x1 = r_cmin;
	 x2 = om->s[i].x2 + (n - k + 1);
	 if(x2 > r_cmax) x2 = r_cmax;
	 
	 onew->s[j].y = y - k;
	 onew->s[j].x1 = x1;
	 onew->s[j].x2 = x2;
	 j++;
      }
   }
/*
 * Grow the easy spans that already exist
 */
   for(i = 0;i < nspan;i++) {
      x1 = om->s[i].x1 - n;
      if(x1 < r_cmin) x1 = r_cmin;
      x2 = om->s[i].x2 + n;
      if(x2 > r_cmax) x2 = r_cmax;

      onew->s[j].y = om->s[i].y;
      onew->s[j].x1 = x1;
      onew->s[j].x2 = x2;
      j++;
   }
/*
 * and grow up above rmax
 */
   y = om->rmax;
   for(i = nspan - 1;i > nspan - nrmax - 1;i--) {
      for(k = 1;k <= ngmax;k++) {
	 x1 = om->s[i].x1 - (n - k + 1);
	 if(x1 < r_cmin) x1 = r_cmin;
	 x2 = om->s[i].x2 + (n - k + 1);
	 if(x2 > r_cmax) x2 = r_cmax;
	 
	 onew->s[j].y = y + k;
	 onew->s[j].x1 = x1;
	 onew->s[j].x2 = x2;
	 j++;
      }
   }

   onew->nspan = j;
/*
 * some of those new spans may overlap
 */
   phCanonizeObjmask(onew, 1);
   
   return(onew);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Grow a chain of OBJMASKs, merging those that grow into each other.
 */
void
phObjmaskChainGrow(CHAIN *chain,	/* the chain of OBJMASKs to grow */
		  const REGION *reg,	/* REGION that objects live in */
		  int n)		/* by how many pixels to grow */
{
   CURSOR_T curs1, curs2;		/* cursors for chain */
   int i;
   OBJMASK *om1;			/* OBJMASK to grow */
   OBJMASK *om2;			/* another OBJMASK on the chain */
   OBJMASK *uni;			/* union of om1 and om2 */
   char *objmask_type = (char *)shTypeGetFromName("OBJMASK");

   shAssert(chain != NULL && chain->type == (TYPE)objmask_type);
   shAssert(reg != NULL && reg->type == TYPE_PIX);
   
   for(i = shChainSize(chain) - 1; i >= 0; i--) {
      om1 = shChainElementGetByPos(chain, i);
      om2 = phObjmaskGrow(om1, reg, n);

      (void)shChainElementChangeByPos(chain, i, om2);
      phObjmaskDel(om1);
   }
/*
 * We've grown all of the objmasks, now try merging them
 */
   curs1 = shChainCursorNew(chain);
   curs2 = shChainCursorNew(chain);
   while((om1 = shChainWalk(chain,curs1,NEXT)) != NULL) {
      shChainCursorSet(chain, curs2, HEAD);
      while((om2 = shChainWalk(chain,curs2,NEXT)) != om1) {
	 if(phObjmaskIntersect(om1, om2, 0, 0)) {
/*
 * The two overlap. Replace om1 by union of objmasks, and delete om2
 */
	    uni = phObjmaskUnion(om1, om2);
	    (void)shChainElementChangeByCursor(chain, curs1, uni);
	    phObjmaskDel(om1);
	    (void)shChainElementRemByCursor(chain, curs2);
	    phObjmaskDel(om2);

	    om1 = uni;
	 }
      }
   }
   
   shChainCursorDel(chain, curs1);
   shChainCursorDel(chain, curs2);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Set the bounding box for an OBJMASK; we update npix while we are about it
 */
void
phObjmaskBBSet(OBJMASK *om)
{
   int i;
   int npix;
   int nspans = om->nspan;
   int xmin, xmax;
   SPAN *spans = om->s;
   int x1, x2;				/* unpacked from spans */

   npix = 0;
   xmin = MAXSHORT;
   xmax = -MAXSHORT;
   for(i = 0; i < nspans; i++) {
      x1 = spans[i].x1; 
      x2 = spans[i].x2; 
      if(x1 < xmin) xmin = x1;
      if(x2 > xmax) xmax = x2;
      npix += x2 - x1 + 1;
   }
   om->npix = npix;
   
   if(nspans > 0) {
      om->cmin = xmin;
      om->rmin = spans[0].y;
      om->cmax = xmax;
      om->rmax = spans[nspans-1].y;
   } else {
      om->rmin = 0;
      om->cmin = 0;
      om->rmax = 0;
      om->cmax = 0;
   }
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Make an Objmask given the LL and UR corners (x == col, y == row)
 */
OBJMASK *
phObjmaskFromRect(int x1, int y1, int x2, int y2)
{
   int nspan = y2 - y1 + 1;
   OBJMASK *om;
   int i;

   shAssert(x1 <= x2);
   shAssert(y1 <= y2);
   
   om = phObjmaskNew(nspan);

   for(i = 0;i < nspan;i++) {
      om->s[i].y = y1 + i;
      om->s[i].x1 = x1;
      om->s[i].x2 = x2;
   }
   om->nspan = nspan;

   om->cmin = x1; om->cmax = x2;
   om->rmin = y1; om->rmax = y2;
   om->npix = (y2 - y1 + 1)*(x2 - x1 + 1);
   
   return(om);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Make an circular Objmask given the centre and radius (x == col, y == row)
 */
OBJMASK *
phObjmaskFromCircle(int yc, int xc, int r)
{
   int npix = 0;			/* number of pixels in OBJMASK */
   int nspan = 2*r + 1;
   OBJMASK *om;
   const int r2 = r*r;
   int i;
   int hlen;				/* half-length of chord */
   
   om = phObjmaskNew(nspan);

   npix = 0;
   for(i = 0;i <= r;i++) {
      hlen = sqrt(r2 - i*i);
      om->s[r + i].y = yc + i;
      om->s[r + i].x1 = xc - hlen;
      om->s[r + i].x2 = xc + hlen;

      om->s[r - i].y = yc - i;
      om->s[r - i].x1 = xc - hlen;
      om->s[r - i].x2 = xc + hlen;

      npix += 2*(2*hlen + 1);
   }
   om->nspan = nspan;

   om->cmin = xc - r; om->cmax = xc + r;
   om->rmin = yc - r; om->rmax = yc + r;
   om->npix = npix;
   
   return(om);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Fill out the OBJMASK->data field of an OBJMASK from a region; the
 * OBJMASK may be NULL in which case nothing is set
 */
void
phObjmaskSetFromRegion(OBJMASK *om,	/* the OBJMASK to set */
		       const REGION *reg) /* the region wherein it dwells */
{
   int i;
   int npix;				/* number of pixels in OBJMASK */
   PIX *ptr;				/* pointer to allocated space */
   SPAN *sp;				/* == om->s[i] */
#if !defined(NDEBUG)
   int nrow, ncol;			/* == reg->n{row,col} */
#endif

   if(om == NULL) return;
   
   shAssert(om->npix >= 0);
   shAssert(reg != NULL && reg->type == TYPE_PIX);
#if !defined(NDEBUG)
   nrow = reg->nrow; ncol = reg->ncol;
#endif
/*
 * allocate space for data
 */
   if(om->data != NULL) {
      shFree(om->data);
   }

   om->data = shMalloc(om->npix*sizeof(PIX));
/*
 * and fill that space from the region
 */
   ptr = om->data;
   for(i = 0;i < om->nspan;i++) {
      sp = &om->s[i];
      npix = sp->x2 - sp->x1 + 1;
      shAssert(sp->y >= 0 && sp->y < nrow && sp->x1 >= 0 && sp->x2 < ncol);
      memcpy(ptr,&reg->ROWS[sp->y][sp->x1],npix*sizeof(PIX));

      ptr += npix;
   }
   shAssert(ptr == om->data + om->npix);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Set pixels in a REGION from an OBJMASK. The OBJMASK may be NULL
 */
void
phRegionSetFromObjmask(REGION *reg,	/* the region to set */
		       const OBJMASK *om) /* and the OBJMASK */
{
   int i;
   int npix;				/* number of pixels in SPAN */
   PIX *ptr;				/* pointer to OBJMASK's data */
   SPAN *sp;				/* == om->s[i] */
#if !defined(NDEBUG)
   int nrow, ncol;			/* == reg->n{row,col} */
#endif

   if(om == NULL) return;

   shAssert(om->nspan == 0 || om->data != NULL);
   shAssert(reg != NULL && reg->type == TYPE_PIX);
#if !defined(NDEBUG)
   nrow = reg->nrow; ncol = reg->ncol;
#endif
/*
 * Set values in the region
 */
   ptr = om->data;
   for(i = 0;i < om->nspan;i++) {
      sp = &om->s[i];
      npix = sp->x2 - sp->x1 + 1;
      shAssert(sp->y >= 0 && sp->y < nrow && sp->x1 >= 0 && sp->x2 < ncol);
      memcpy(&reg->ROWS[sp->y][sp->x1],ptr,npix*sizeof(PIX));

      ptr += npix;
   }
   shAssert(ptr == om->data + om->npix);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Set pixels in a REGION to a scalar <val> wherever an OBJMASK has pixels;
 * the OBJMASK may be NULL
 */
void
phRegionSetValFromObjmask(REGION *reg,	/* the region to set */
			  const OBJMASK *om, /* the OBJMASK defining pixels */
			  int val)	/* to set to this value */
{
   int i,j;
   PIX *row;				/* pointer to row of REGION */
   SPAN *sp;				/* == om->s[i] */
   int nrow, ncol;			/* == reg->n{row,col} */
   int row0, col0;			/* == reg->{row,col}0 */
   int y, x1, x2;			/* == sp->{y, x1, x2} */

   if(om == NULL) return;
   
   shAssert(reg != NULL && reg->type == TYPE_PIX);
   nrow = reg->nrow; ncol = reg->ncol;
   row0 = reg->row0; col0 = reg->col0;
/*
 * Set values in the region
 */
   for(i = 0;i < om->nspan;i++) {
      sp = &om->s[i];
      y = sp->y;
      if(y < row0 || y >= row0 + nrow) {
	 continue;
      }
      x1 = sp->x1; x2 = sp->x2;
      if(x1 < col0) { x1 = col0; }
      if(x2 >= col0 + ncol) { x2 = col0 + ncol - 1; }
      
      row = &reg->ROWS[y - row0][x1 - col0];
      for(j = x2 - x1;j >= 0;j--) {
	 row[j] = val;
      }
   }
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Set pixels in a REGION to a scalar <val> wherever an OBJMASK on the
 * CHAIN objmasks has pixels
 */
void
phRegionSetValFromObjmaskChain(REGION *reg, /* the region to set */
			       const CHAIN *chain, /* the CHAIN of OBJMASKs
						      specifying pixels */
			       int val)	/* to set to this value */
{
   int i;
   const OBJMASK *om;			/* an OBJMASK from chain */
   int rmin, rmax, cmin, cmax;		/* bounding box for region */
   int size;				/* number of elements in chain */
   
   shAssert(reg != NULL && reg->type == TYPE_PIX);
   shAssert(chain != NULL && chain->type == objmask_type);

   rmin = reg->row0; rmax = reg->row0 + reg->nrow - 1;
   cmin = reg->col0; cmax = reg->col0 + reg->ncol - 1;

   size = shChainSize(chain);
   for(i = 0; i < size; i++) {
      om = shChainElementGetByPos(chain, i);

      if(om->rmax < rmin || om->cmax < cmin ||
	 om->rmin > rmax || om->cmin > cmax) {
	 continue;			/* no overlap */
      }

      phRegionSetValFromObjmask(reg, om, val);
   }
}
#endif /* !defined(STAND_ALONE) */

/*****************************************************************************/
/*
 * Make a SPANMASK but don't allocate chains.
 */
static SPANMASK *
make_spanmask(void)
{
   SPANMASK *sm = shMalloc(sizeof(*sm));
   sm->cookie = SPAN_COOKIE;
   *(SPANMASK **)(&sm->parent) = NULL;
   sm->nchild = 0;

   return(sm);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Create a complete SPANMASK, corresponding to an nrow*ncol bitmask
 */
SPANMASK *
phSpanmaskNew(int nrow, int ncol)
{
    SPANMASK *sm;
    int im;
    char *cobjmask_type = (char *)objmask_type;

    shAssert(NMASK_TYPES == (int)S_NMASK_TYPES);

    sm = make_spanmask();
    sm->nrow = nrow; sm->ncol = ncol;
    
    for(im = 0; im < NMASK_TYPES; im++)
	sm->masks[im] = shChainNew(cobjmask_type);
    return sm;
}    

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Make a sub-SPANMASK (currently, simply return the whole mask)
 */
SPANMASK *
phSubSpanmaskNew(
		 const SPANMASK *ism
		 )
{
   int i;
   SPANMASK *sm = make_spanmask();

   ((SPANMASK *)ism)->nchild++;
   sm->parent = ism;

   for(i = 0;i < NMASK_TYPES;i++) {
      sm->masks[i] = ism->masks[i];
   }
   
   return(sm);
}    

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Make a copy of a complete SPANMASK
 */
SPANMASK *
phSpanmaskCopy(
	       const SPANMASK *ism,	/* original mask */
	       float drow, float dcol	/* shift of new wrt old mask */
	       )
{
    SPANMASK *sm;
    int im;

    if(ism == NULL) {
       return(NULL);
    }

    sm = make_spanmask();
    for(im = 0; im < NMASK_TYPES; im++) {
       sm->masks[im] = phObjmaskChainCopy(ism->masks[im], drow, dcol);
    }

    return sm;
}    

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Delete a complete SPANMASK
 */
void
phSpanmaskDel(
	      SPANMASK *mask
	      )
{
    int im;
    
    if(mask == NULL) return;
      
    shAssert(mask->cookie == SPAN_COOKIE);

    if(mask->parent != NULL) {
       ((SPANMASK *)mask->parent)->nchild--; /* we don't own the OBJMASKs,
						   so don't free them */
    } else {
       if(mask->nchild != 0) {
	  shError("Attempt to delete SPANMASK with live children\n");
       } else {
	  for(im = 0; im < NMASK_TYPES; im++) {
	     phObjmaskChainDel(mask->masks[im]);
	  }
       }
    }
    shFree(mask);
}    

#if !defined(STAND_ALONE)
/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Add an OBJMASK to a specified plane of a SPANMASK, ensuring that the
 * SPANMASK remains canonical
 */
void
phObjmaskAddToSpanmask(const OBJMASK *om, /* OBJMASK to be added */
		       SPANMASK *smask,	/* SPANMASK to be augmented */
		       S_MASKTYPE which) /* desired plane of smask */
{
   CURSOR_T curs;			/* cursor for chain */
   CHAIN *mask;				/* a plane of smask */
   OBJMASK *sv;				/* an OBJMASK on the mask CHAIN */
   int rmin, cmin;			/* unpacked from om */

   shAssert(om != NULL && smask != NULL && which >=0 && which < S_NMASK_TYPES);

   mask = smask->masks[which];
   curs = shChainCursorNew(mask);

   rmin = om->rmin; cmin = om->cmin;
   while((sv = shChainWalk(mask,curs,NEXT)) != NULL) {
      if(sv->rmin >= rmin) {
	 if(sv->rmin > rmin || sv->cmin > cmin) {
	    break;
	 }
      }
   }
   if(sv == NULL) {
      shChainElementAddByPos(mask, (void *)om,
			      (char *)objmask_type, TAIL, AFTER);
   } else {
      shChainElementAddByCursor(mask, (void *)om,
			      (char *)objmask_type, curs, BEFORE);
   }
   shChainCursorDel(mask, curs);
}
   
/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Create a chain for a SPANMASK from a rectangle.
 */
CHAIN *
phSpanmaskFromRect(
		   int x1, 
		   int y1, 
		   int x2, 
		   int y2
		   )
{
    CHAIN *mask;
    char *cobjmask_type = (char *)objmask_type;

    mask = shChainNew(cobjmask_type);
    
    (void)shChainElementAddByPos(mask,
				 phObjmaskFromRect(x1,y1,x2,y2),
				 cobjmask_type,TAIL,AFTER);
    return mask;
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Change row0, col0, bounding box, and spans so that SPANMASK
 * looks more like a submask.
 */
SPANMASK *
phSpanmaskResetCorner(
		      SPANMASK *mask,	/* mask to change */
		      int row0,		/* desired value of row0 */
		      int col0		/*                  and of col0 */
		      )
{
   int i;
   CURSOR_T curs;			/* cursor for chain */
   CHAIN *chain;			/* a chain of OBJMASKs */
   OBJMASK *obj;			/* one of the OBJMASKs */

   for(i = 0;i < NMASK_TYPES;i++) {
      chain = mask->masks[i];
      curs = shChainCursorNew(chain);
      while((obj = shChainWalk(chain,curs,NEXT)) != NULL) {
	 phObjmaskResetCorner(obj, row0, col0);
      }
      shChainCursorDel(chain, curs);
   }

   return(mask);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Change row0, col0, bounding box, and spans so that OBJMASK
 * looks more like a submask.
 */
OBJMASK *
phObjmaskResetCorner(
		     OBJMASK *mask,
		     int row0,
		     int col0
		     )
{
    int i,drow, dcol;
    SPAN *spans;
    
    shAssert(row0 >= 0);
    shAssert(col0 >= 0);

    if ((row0 != mask->row0) || (col0 != mask->col0)) {
	
	drow = row0 - mask->row0;
	dcol = col0 - mask->col0;

	spans = mask->s;
	for (i=0; i<mask->nspan; i++) {
	    spans[i].y -= drow;
	    spans[i].x1 -= dcol;
	    spans[i].x2 -= dcol;
	}
	mask->rmin -= drow;
	mask->rmax -= drow;
	mask->cmin -= dcol;
	mask->cmax -= dcol;
	mask->row0 = row0;
	mask->col0 = col0;
    }

    return(mask);
}

/*****************************************************************************/
/*
 * Find the span in an OBJMASK that contains pixel (rc, cc)
 */
int
phObjmaskFindSpan(const OBJMASK *om,	/* the mask */
		  int rc, int cc,	/* desired pixel */
		  int i0)		/* initial guess for index, or -1 */
{
   int i;
   const int cmin = (om == NULL) ? -1 : om->cmin;
   const int cmax = (om == NULL) ? -1 : om->cmax;
   const int rmin = (om == NULL) ? -1 : om->rmin;
   const int rmax = (om == NULL) ? -1 : om->rmax;
   const int nspan = (om == NULL) ? -1 : om->nspan;
   const SPAN *s = (om == NULL) ? NULL : om->s;
/*
 * Is position above/below top/bottom boundary of mask?
 */
   if(om == NULL || nspan == 0 ||
      rc < rmin || rc > rmax || cc < cmin || cc > cmax) {
      return(-1);
   }
/*
 * find the first span corresponding to the position
 */
   if(i0 < 0) {
      i0 = (rc - rmin)*nspan/(float)(rmax - rmin + 1);
   }

   if(i0 < 0) i0 = 0;
   if(i0 >= nspan) i0 = nspan - 1;

   i = i0;				/* initial guess */

   if(s[i].y >= rc) {
      for(; i >= 0 && s[i].y >= rc; i--) continue;
      if(i < 0 || s[i].y < rc) i++;
   } else if(s[i].y < rc) {
      for(; i < nspan && s[i].y < rc; i++) continue;
      if(i >= nspan) i--;
   }
/*
 * index i should now give the first span with y == rc
 *
 * Find a span that includes cc
 */
   shAssert(i >= 0 && i < nspan);

   for(; i < nspan && s[i].y == rc; i++) {
      if(cc >= s[i].x1 && cc <= s[i].x2) {
	 break;
      }
   }
   if(i == nspan || s[i].y != rc) {
      return(-1);			/* no overlap */
   }

   return(i);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * return an OBJMASK's transpose
 */
OBJMASK *
phObjmaskTranspose(const OBJMASK *om)
{
   int cmin;				/* == om->cmin */
   int i;
   int id;				/* a span's ID == (index in st) + 1 */
   int *id0;				/* transposed span ids for this */
   int *idm;				/*             and previous row */
   int onspan;				/* == om->nspan */
   int tnspan;				/* == trans->nspan */
   SPAN *so, *st;			/* == {om,trans}->s */
   OBJMASK *trans;			/* the required transpose */
   int width;				/* == om->cmax - om->cmin + 1 */
   int x;				/* counter in x */
   int x1, x2;				/* == so[].x[12] - cmin */
   int y, ym;				/* y for this and previous row */

   if(om == NULL) {
      return(NULL);
   } else if(om->nspan == 0) {
      return(phObjmaskNew(0));
   }

   cmin = om->cmin;
   width = om->cmax - om->cmin + 1;	/* width of bounding box */
   trans = phObjmaskNew(width + 3);
   so = om->s; st = trans->s;
   onspan = om->nspan; tnspan = trans->nspan;

   id0 = alloca(2*width*sizeof(int));
   idm = id0 + width;
   memset(id0, '\0', width*sizeof(int));
/*
 * Now go through om looking at each span creating transposed spans.
 * If a pixel is in the same column as one in the previous span simply
 * increment the proper (transposed) spans's x2 value, otherwise start
 * a new (transposed) span. The arrays id0 and idm contain the ids for
 * the transposed spans for this and the previous row, which are one
 * more than the spans' indices in trans->s
 */
   ym = so[0].y - 1;
   for(i = 0; i < onspan; i++) {
      y = so[i].y;
      x1 = so[i].x1 - cmin; x2 = so[i].x2 - cmin;
      if(y != ym) {			/* switch id0 and idm */
	 int *tmp = id0; id0 = idm; idm = tmp;
	 memset(id0, '\0', width*sizeof(int));
      }
      for(x = x1; x <= x2; x++) {
	 if((id = idm[x]) == 0) {
	    id0[x] = ++tnspan;
	    if(tnspan >= trans->size) {
	       phObjmaskRealloc(trans, 1.5*trans->size + 1);
	       st = trans->s;
	    }
	    st[tnspan - 1].y = x + cmin;
	    st[tnspan - 1].x1 = st[tnspan - 1].x2 = y;
	 } else {
	    st[id - 1].x2++;
	    id0[x] = id;
	 }
      }
      
   }
   trans->nspan = tnspan;

   phCanonizeObjmask(trans, 1);

   return(trans);
}

/*****************************************************************************/
/*
 * comparison function for qsort
 */
static int
compar(const void *pa, const void *pb)
{
   const OBJMASK *a = *(OBJMASK **)pa;
   const OBJMASK *b = *(OBJMASK **)pb;

   if(a->rmin < b->rmin) {
      return(-1);
   } else if(a->rmin == b->rmin) {
      return(a->cmin - b->cmin);
   } else {
      return(1);
   }
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Put a SPANMASK into canonical form, namely sorted by LL corner of
 * bounding box.
 *
 * If canonize_objmasks is true, first ensure that each OBJMASK is canonical
 * (sort in y, sort in x, find overlaps, find the bounding box).
 *
 * If nearly_sorted is true, use an algorithm that assumes that the
 * SPANMASK is almost in proper order
 */
void
phCanonizeObjmaskChain(CHAIN *chain,	/* CHAIN of OBJMASKs to work on */
		       int canonize_objmasks, /* canonize each OBJMASK too? */
		       int nearly_sorted) /* are OBJMASKs almost sorted? */
{
   CURSOR_T crsr;			/* cursor for crChain */
   OBJMASK *om;				/* an OBJMASK from the chain */
   
   shAssert(chain != NULL && chain->type == objmask_type);

   if(canonize_objmasks) {
      crsr = shChainCursorNew(chain);
      while((om = shChainWalk(chain, crsr, NEXT)) != NULL) {
	 phCanonizeObjmask(om, nearly_sorted);
      }
      shChainCursorDel(chain, crsr);
   }
/*
 * Sort by the lower-left corner of the bounding box, in Y then X
 */
   shChainQsort(chain, compar);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Put a OBJMASK into canonical form: sort in y, sort in x, find
 * overlaps, find the bounding box.
 *
 * (Yes, we know that it should be "canonicalize", not canonise)
 */
void
phCanonizeObjmask(
		  OBJMASK *sv,		/* OBJMASK to work on */
		  int nearly_sorted	/* is sv almost sorted already? */
		  )
{
    SPAN *spans;		/* Array of spans */
    int total_spans;
    int i;
    int npix;				/* number of pixels in mask */
    int y, x1, x2;			/* unpacked from spans[i].{y,x1,x2} */
    int ymin, ymax, xmin, xmax;
    int nspans;
    SPAN *s1;

    shAssert(sv != NULL);
    shAssert(sv->nspan <= sv->size);
    spans = sv->s;
/*
 * Sort in Y and X
 *
 * Don't use the system qsort, as naive quick sort (which qsort often is)
 * has an N^2 worst-case running time; also heap_sort is simple enough that
 * the comparison can be coded inline. In addition, if the list is nearly
 * sorted, it's still better to use an insertion sort (with linear time and
 * a small coefficient for sorted lists).
 */
    if(nearly_sorted) {
       insertion_sort(spans, sv->nspan);
    } else {
       heap_sort(spans, sv->nspan);
    }
				/* Uniquify spans */
    i = 0;
    while(i < sv->nspan) {
	s1 = &spans[i++];
	if(s1->x1 > s1->x2) {		/* empty span */
	   s1->y = -MAXSHORT + 1;	/* tag as invalid */
	   continue;
	}
	while(i < sv->nspan &&
	      s1->y == spans[i].y && s1->x2 >= spans[i].x1 - 1) {
	    spans[i].y = -MAXSHORT + 1; /* tag as invalid */
	    if(spans[i].x2 > s1->x2)
		s1->x2 = spans[i].x2;
	    i++;
	}
    }
				/* Compress; and find bounding box and npix */
    ymin = MAXSHORT;
    ymax = -MAXSHORT;
    xmin = MAXSHORT;
    xmax = -MAXSHORT;
    nspans = sv->nspan;
    npix = 0;
    for(total_spans = 0, i = 0; i < nspans; i++) {
	if(spans[i].y == -MAXSHORT + 1)
	    continue;
	spans[total_spans++] = spans[i];
	y = spans[i].y;
	x1 = spans[i].x1;
	x2 = spans[i].x2;
	if(y < ymin) ymin = y;
	if(y > ymax) ymax = y;
	if(x1 < xmin) xmin = x1;
	if(x2 > xmax) xmax = x2;
	npix += (x2 - x1 + 1);
    }
    if(nspans > 0) {
        sv->npix = npix;
	sv->rmin = ymin;
	sv->rmax = ymax;
	sv->cmin = xmin;
	sv->cmax = xmax;
    } else {
        sv->npix = 0;
	sv->rmin = 0;
	sv->rmax = 0;
	sv->cmin = 0;
	sv->cmax = 0;
    }
    sv->nspan = total_spans;	/* Should I realloc?  Probably a waste */
				/* of time. */
    /* We now have St. OBJMASK */
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Determine if two OBJMASKs overlap.
 */
int
phObjmaskIntersect(const OBJMASK *sv1, 
		   const OBJMASK *sv2,
		   float fdy,		/* offset to be applied to sv2 */
		   float fdx)
{
   const int dx = (fdx < 0) ?  -(-fdx + 0.5) : fdx + 0.5;
   const int dy = (fdy < 0) ?  -(-fdy + 0.5) : fdy + 0.5;
   int irow1, irow2;
   int xoff1, xoff2, yoff1, yoff2;
   SPAN *s1, *s2;

   if(sv2->cmax + dx < sv1->cmin /* check bounding box */
      || sv1->cmax < sv2->cmin + dx
      || sv2->rmax + dy < sv1->rmin
      || sv1->rmax < sv2->rmin + dy)
       return 0;

   xoff1 = sv1->col0;
   xoff2 = sv2->col0 + dx;
   yoff1 = sv1->row0;
   yoff2 = sv2->row0 + dy;
   s1 = sv1->s;
   s2 = sv2->s;

   irow1 = 0;
   irow2 = 0;
   while(irow1 < sv1->nspan && irow2 < sv2->nspan) {
       int icol1, icol2;

       /* find a span with a common y */
       if(s1[irow1].y + yoff1 < s2[irow2].y + yoff2) {
	   irow1++;
	   continue;
       }
       if(s2[irow2].y + yoff2 < s1[irow1].y + yoff1) {
	   irow2++;
	   continue;
       }
			/* loop over spans with this y */
       for(icol1 = irow1, icol2 = irow2;
	   icol1 < sv1->nspan && s1[icol1].y == s1[irow1].y
	   && s2[icol2].y == s2[irow2].y; icol1++) {

	   while(icol2 < sv2->nspan && s2[icol2].y == s2[irow2].y
		 && s2[icol2].x1 + xoff2 <= s1[icol1].x2 + xoff1) {

	       if(s2[icol2].x2 + xoff2 >= s1[icol1].x1 + xoff1) {
		   return 1;
	       }
	       icol2++;
	   }
       }
/*
 * No intersection at this y value; proceed to next spans
 */
       irow1 = icol1; irow2 = icol2;
   }
   return 0;
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Merge the second OBJMASK into the first; if the first is NULL, it
 * will be allocated
 *
 * It is asserted that the OBJMASKs have the same col0 and row0
 *
 * The return value is sv1, possible as allocated
 */
OBJMASK *
phObjmaskMerge(OBJMASK *sv1,		/* target OBJMASK */
	       const OBJMASK *sv2,	/* input OBJMASK */
	       float fdy,		/* Shift sv2 by this amount */
	       float fdx)		/* before merging */
{
   const int dx = (fdx < 0) ?  -(-fdx + 0.5) : fdx + 0.5;
   const int dy = (fdy < 0) ?  -(-fdy + 0.5) : fdy + 0.5;
    int i;
    int i1, i2;				/* counters for sv1, sv2 */
    long k1, k2;			/* keys in sort; >= 32bits */
    SPAN *new_spans;			/* the output SPANs */
    int nspan;
    const int nspan1 = (sv1 == NULL) ? 0 : sv1->nspan; /* unalias them */
    const int nspan2 = (sv2 == NULL) ? 0 : sv2->nspan;

    shAssert(sv2 != NULL);
    if(sv1 != NULL) {
       shAssert(sv1->col0 == sv2->col0 && sv1->row0 == sv2->row0);
    }

    nspan = nspan1 + nspan2;
/*
 * special case one OBJMASK being empty
 */
    if(nspan2 == 0) {
       return(sv1);
    } else if(nspan1 == 0) {
       sv1 = phObjmaskRealloc(sv1, nspan2);
       sv1->row0 = sv2->row0;
       sv1->col0 = sv2->col0;
       sv1->rmin = sv2->rmin + dy;
       sv1->rmax = sv2->rmax + dy;
       sv1->cmin = sv2->cmin + dx;
       sv1->cmax = sv2->cmax + dx;
       sv1->npix = sv2->npix;
       for(i = 0; i < nspan2; i++) {
	  SPAN *s1 = &sv1->s[i];	/* unalias them */
	  SPAN *s2 = &sv2->s[i];
	  s1->y = s2->y + dy;
	  s1->x1 = s2->x1 + dx;
	  s1->x2 = s2->x2 + dx;
       }
       sv1->nspan = nspan;
       return(sv1);
    }
/*
 * neither is empty, we have real work to do. Go through both OBJMASKs,
 * inserting them into the output in sorted order
 */
    new_spans = shMalloc(nspan*sizeof(SPAN));

    i1 = i2 = 0;
    k1 = SORTKEY(sv1->s[i1].x1, sv1->s[i1].y);
    k2 = SORTKEY(sv2->s[i2].x1 + dx, sv2->s[i2].y + dy);
    
    i = 0;
    for(;;) {
       if(k1 < k2) {
	  new_spans[i++] = sv1->s[i1++];
	  if(i1 == nspan1) {
	     while(i2 < nspan2) {
		SPAN *s1 = &new_spans[i++]; /* unalias them */
		SPAN *s2 = &sv2->s[i2++];
		s1->y = s2->y + dy;
		s1->x1 = s2->x1 + dx;
		s1->x2 = s2->x2 + dx;
	     }
	     break;
	  }
	  k1 = SORTKEY(sv1->s[i1].x1, sv1->s[i1].y);
       } else {
	  SPAN *s1 = &new_spans[i++]; /* unalias them */
	  SPAN *s2 = &sv2->s[i2++];

	  s1->y = s2->y + dy;
	  s1->x1 = s2->x1 + dx;
	  s1->x2 = s2->x2 + dx;
	  if(i2 == nspan2) {
	     while(i1 < nspan1) {
		new_spans[i++] = sv1->s[i1++];
	     }
	     break;
	  }
	  k2 = SORTKEY(sv2->s[i2].x1 + dx, sv2->s[i2].y + dy);
       }
    }

    shFree(sv1->s);
    sv1->s = new_spans;
    sv1->nspan = sv1->size = nspan;
    phCanonizeObjmask(sv1,1);

    return(sv1);
}

/*****************************************************************************/
/*
 * Set a circular patch in an OBJMASK, which is then returned.
 *
 * If the OBJMASK is NULL it will be created. 
 */
OBJMASK *
phObjmaskMergeWithCircle(OBJMASK *om,	/* mask to set (or create if NULL) */
			 float rowc, float colc, /* centre of circle */
			 float r)	/* radius */
{
   OBJMASK *disk = phObjmaskFromCircle(rowc, colc, r);
   
   if(om == NULL) {
      return(disk);
   } else {
      phObjmaskMerge(om, disk, 0, 0);
      phObjmaskDel(disk);
      
      return(om);
   }
}

/*****************************************************************************/
/*
 * Clear a circular patch in an OBJMASK
 */
void
phObjmaskClearWithCircle(OBJMASK *om,	/* mask to clear part of */
			 float rowc, float colc, /* centre of circle */
			 float r)	/* radius */
{
   OBJMASK *disk = phObjmaskFromCircle(rowc, colc, r);
   
   if(om != NULL) {
      phObjmaskAndNotObjmask(om, disk);
      phObjmaskDel(disk);
   }
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Merge a chain of OBJMASKs into a single OBJMASK.
 * We assert that the chain is canonical (i.e. sorted in y-then-x, and with
 * row0,col0 the same for OBJMASKs)
 */
OBJMASK *
phMergeObjmaskChain(
		    const CHAIN *chain /* array of OBJMASK */
		    )
{
   CURSOR_T curs;			/* cursor for chain */
   OBJMASK *sv;				/* OBJMASKs within chain */
   OBJMASK *merge;			/* Merged OBJMASK */
   SPAN *mspans;
   int ispan;
   int total_spans;		/* total spans in merged spanvec */

   shAssert(chain != NULL && chain->type == objmask_type);
   shAssert(is_canonical_objmaskChain(chain));

   if(chain->nElements == 0) {
      merge = phObjmaskNew(1);
      merge->npix = merge->nspan = 0;

      return(merge);
   }

   curs = shChainCursorNew(chain);
   total_spans = 0;
   while((sv = shChainWalk(chain,curs,NEXT)) != NULL) {
       total_spans += sv->nspan;
   }
   merge = phObjmaskNew(total_spans);
   merge->nspan = total_spans;
   sv = (OBJMASK *)shChainElementGetByPos((CHAIN *)chain,HEAD);
   merge->row0 = sv->row0; merge->col0 = sv->col0;

   mspans = merge->s;
   ispan = 0;
   (void)shChainCursorSet(chain, curs, HEAD);
   while((sv = shChainWalk(chain,curs,NEXT)) != NULL) {
       int i;
       const int sv_nspan = sv->nspan;	/* unalias them */
       SPAN *sv_s = sv->s;
       for(i = 0; i < sv_nspan; i++) {
	  mspans[ispan++] = sv_s[i];
       }
   }
   shChainCursorDel(chain, curs);

   phCanonizeObjmask(merge,0);

   return merge;
}

/*
 * <AUTO EXTRACT>
 * Merge two chains of OBJMASKs into a single chain. Note that both
 * of the chains are destroyed
 *
 * In the resulting CHAIN, it is guaranteed that the OBJMASKs are sorted
 */
CHAIN *
phObjmaskChainMerge(
		   CHAIN *c1,		/* array of OBJMASK */
		   CHAIN *c2		/* array of OBJMASK */
		   )
{
   CURSOR_T cursj;			/* cursor for join */
   CURSOR_T curs1, curs2;		/* cursor for c1 and c2 */
   OBJMASK *el1, *el2;			/* elements of c1 and c2 */
   CHAIN *join;				/* result */
   long k1, k2;				/* keys in sort; >= 32bits */
   const char *cobjmask_type = (char *)objmask_type;
   
   shAssert(c1 != NULL);
   shAssert(c2 != NULL);
   shAssert(c1->type == objmask_type);
   shAssert(shChainTypeGet(c1) == shChainTypeGet(c2));
   
   if(shChainSize(c1) == 0) {
      shChainDel(c1);
      return(c2);
   } else if(shChainSize(c2) == 0) {
      shChainDel(c2);
      return(c1);
   }
/*
 * OK, we have to do some work as we have two non-empty chains
 */
   shAssert(is_canonical_objmaskChain(c1));
   shAssert(is_canonical_objmaskChain(c2));
/*
 * Build a sorted list
 */
   join = shChainNew(cobjmask_type);
   cursj = shChainCursorNew(join);
   curs1 = shChainCursorNew(c1);
   curs2 = shChainCursorNew(c2);
   
   el1 = shChainWalk(c1,curs1,NEXT);	/* we know that size >= 1 */
   el2 = shChainWalk(c2,curs2,NEXT);	/* "    "    "   "   "  " */
   k1 = SORTKEY(el1->cmin, el1->rmin);
   k2 = SORTKEY(el2->cmin, el2->rmin);

   for(;;) {
      if(k1 < k2) {
	 shChainElementAddByPos(join, el1, cobjmask_type, TAIL, AFTER);
	 if((el1 = shChainWalk(c1,curs1,NEXT)) == NULL) {
	    do {
	       shChainElementAddByPos(join, el2, cobjmask_type, TAIL, AFTER);
	       shChainElementRemByCursor(c2, curs2);
	    } while((el2 = shChainWalk(c2,curs2,NEXT)) != NULL);
	    break;
	 }
	 k1 = SORTKEY(el1->cmin, el1->rmin);
      } else {
	 shChainElementAddByPos(join, el2, cobjmask_type, TAIL, AFTER);
	 if((el2 = shChainWalk(c2,curs2,NEXT)) == NULL) {
	    do {
	       shChainElementAddByPos(join, el1, cobjmask_type, TAIL, AFTER);
	       shChainElementRemByCursor(c1, curs1);
	    } while((el1 = shChainWalk(c1,curs1,NEXT)) != NULL);
	    break;
	 }
	 k2 = SORTKEY(el2->cmin, el2->rmin);
      }
   }

   shChainCursorDel(c1,curs1); shChainDel(c1);
   shChainCursorDel(c2,curs2); shChainDel(c2);
   shChainCursorDel(join,cursj);

   shAssert(is_canonical_objmaskChain(join));

   return(join);
}

/*****************************************************************************/
/*
 * Find the inclusive range of chain elements that could contain a row
 * in the range y1 -- y2 inclusive
 */
static void
find_span_range(const CHAIN *chain,	/* chain to search */
		int y1,			/* lower value of y */
		int y2,			/* upper value of y  */
		int *p0,		/* lower limit to search */
		int *p1)		/* upper limit to search */
{
   const int chainsize = chain->nElements;
   int p;				/* index for chain */
#if 0
   unsigned int p0_l, p0_u;		/* range for p0 while searching */
#endif
   unsigned int p1_l, p1_u;		/* range for p1 while searching */
   OBJMASK *sv;

   shAssert(chainsize > 0);
/*
 * Ensure that the chain is indexed
 */
   if(chainsize > 3 && !shChainIsIndexed(chain)) {
      shChainIndexMake(chain, SH_CHAIN_INDEX_BY_POS);
   }
/*
 * Find the part of the chain that could intersect; we desire that
 * indices p0 and p1
 */
   p1_l = 0; p1_u = chainsize - 1;
#if 0
   p0_l = 0; p0_u = chainsize - 1;
/*
 * find p0 first, remembering anything we learn about p1 as we proceed
 */
   for(;;) {
      p = (p0_l + p0_u)/2;
      if(p == p0_l) {
	 break;
      }
      sv = shChainElementGetByPos(chain, p);

      if(sv->rmax < y1) {		/* too small */
	 p0_l = p;
      } else {
	 p0_u = p;
	 if(sv->rmin > y2) {		/* too large even for p1 */
	    p1_u = p;
	 } else {
	    if(p > p1_l) {
	       p1_l = p;
	    }
	 }
      }      
   }
   *p0 = p0_l;
#else
   *p0 = 0;
#endif
/*
 * got p0; look for p1
 */
   if(p1_l < *p0) {
      p1_l = *p0;
   }
   for(;;) {
      p = (p1_l + p1_u + 1)/2;
      if(p == p1_u) {
	 break;
      }
      sv = shChainElementGetByPos(chain, p);

      if(sv->rmin > y2) {		/* too large */
	 p1_u = p;
      } else {
	 p1_l = p;
      }      
   }
   *p1 = p1_u;
/*
 * Deal with complications such as elements with the same rmin,
 * but rmax out of order
 */
   for(p = *p0; p > 0; p--) {
      sv = shChainElementGetByPos(chain, p);
      if(sv->rmax < y1) {
	 break;
      }
   }
   *p0 = p;

   for(p = *p1; p < chainsize - 1; p++) {
      sv = shChainElementGetByPos(chain, p);
      if(y2 < sv->rmin) {
	 break;
      }
   }
   *p1 = p;
}
	
/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Determine if a rectangle intersects the spans in mask.
 */
int
phRectIntersectMask(const CHAIN *chain, /* array of OBJMASK */
		    int x1,		/* boundaries of rectangle */
		    int y1,
		    int x2,
		    int y2)
{
   int i, j;
   int i0, i1;				/* range of i */
   SPAN *s;
   OBJMASK *sv;
   int row0;				/* row0 for first element of chain */
   int rx1, rx2, ry1, ry2;		/* "reduced" x1 etc.; row0 subtracted*/

   shAssert(chain != NULL && chain->type == objmask_type);
   shAssert(x1 <= x2);
   shAssert(y1 <= y2);

   if(shChainSize(chain) < 1) {
      return(0);
   }
/*
 * find the range of possible OBJMASKs
 */
   i0 = i1 = 0;				/* make gcc happy */
   find_span_range(chain, y1, y2, &i0, &i1);

   sv = shChainElementGetByPos(chain, 0); row0 = sv->row0;
   for(i = i0; i <= i1; i++) {
      sv = shChainElementGetByPos(chain, i);
      shAssert(sv->row0 == row0);
      
      ry1 = y1 - sv->row0; ry2 = y2 - sv->row0;
      rx1 = x1 - sv->col0; rx2 = x2 - sv->col0;
      if(rx2 < sv->cmin || rx1 > sv->cmax || ry2 < sv->rmin || ry1 > sv->rmax){
	 continue;			/* no match on bounding box */
      }
      
      s = sv->s;
      for(j = 0; j < sv->nspan; j++) {
	 if(ry2 < s[j].y) {		/* can't intersect, as s is sorted */
	    break;
	 }
	 if(ry1 > s[j].y || rx1 > s[j].x2 || rx2 < s[j].x1) {
	    continue;
	 }
	 return 1;
      }
   }

   return 0;
}
    
/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Determine if a pixel intersects the spans in a OBJMASK.
 */
int
phPixIntersectObjmask(
		      const OBJMASK *sv, /* input SPANVEC */
		      int x,		/* pixel position */
		      int y
		      )
{
   SPAN *s;
   int offy;
   int offx;
   int i;

   shAssert(sv != NULL);

   if(x < sv->cmin || sv->cmax < x /* check bounding boxes */
      || y < sv->rmin || sv->rmax < y)
       return 0;
   s = sv->s;
   offy = sv->row0;
   offx = sv->col0;
   for(i = 0; i < sv->nspan; i++) {
       if(s[i].y + offy < y || y < s[i].y + offy
	  || s[i].x2 + offx < x || x < s[i].x1 + offx)
	   continue;
       return 1;
   }
   return 0;
}
/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Determine if a pixel intersects the spans in mask.
 */
int
phPixIntersectMask(const CHAIN *chain,	/* array of OBJMASKs */
		   int x,		/* pixel position */
		   int y)
{
   int chainsize;			/* number of elements on chain */
   OBJMASK *sv;
   int i;
   int i0, i1;				/* range of possible values of i */
   int p;				/* index for chain */
   int p0, p1;				/* p0 <= p <= p1 could match */
   int row0;				/* row0 for first element of chain */

   shAssert(chain != NULL && chain->type == objmask_type);
   chainsize = chain->nElements;

   if(chainsize < 1) {
      return(0);
   }
/*
 * find the range of possible OBJMASKs
 */
   p0 = p1 = 0;				/* make gcc happy */
   find_span_range(chain, y, y, &p0, &p1);
/*
 * OK, we have to look more carefully at the remaining OBJMASKs
 */
   sv = shChainElementGetByPos(chain, 0); row0 = sv->row0;
   for(p = p0; p <= p1; p++) {
      sv = shChainElementGetByPos(chain, p);
      shAssert(sv->row0 == row0);
      
      if(sv->rmin > y) {                /*no further OBJMASKs can include x,y*/
         return(0);
      }
      if(y > sv->rmax) {
         continue;
      }
      
      if(x >= sv->cmin && x <= sv->cmax) { /* point is within bounding box */
	 const int nspan = sv->nspan;
	 const SPAN *s = sv->s;
	 const int y0 = y - sv->row0;
	 const int x0 = x - sv->col0;
/*
 * Most OBJMASKs are simple, with one SPAN for each row of the image. We
 * can use this to avoid a full binary search
 */
	 i = y - sv->rmin;		/* initial guess */

	 if(i >= nspan) {		/* hmm, not a good guess */
	    i0 = -1; i1 = nspan;
	 } else {
	    unsigned int step = 1;	/* how much to step up/down */
	    if(y0 > s[i].y) {		/* expand search upwards */
	       i0 = i; i1 = i0 + 1;
	       while(y0 >= s[i1].y) {
		  i0 = i1;
		  step += step;		/* double step size */
		  i1 = i0 + step;
		  if(i1 >= nspan) {	/* reached top of array */
		     i1 = nspan - 1;
		     break;
		  }
	       }
	    } else {			/* expand it downwards */
	       i1 = i; i0 = i - 1;
	       if(i0 >= 0) {
		  while(y0 < s[i0].y) {
		     i1 = i0;
		     step += step;	/* double step size */
		     i0 = i1 - step;
		     if(i0 < 0) {	/* off bottom of array */
			i0 = -1;
			break;
		     }
		  }
	       }
	    }
	 }
/*
 * OK, now a binary search to find the first SPAN with the desired row
 */
	 i = (i0 + i1)/2;
	 for(;;) {
	    if(y0 > s[i].y) {
	       i0 = i;
	    } else {
	       i1 = i;
	    }

	    if(i1 - i0 <= 1) {
	       if(i0 == -1 || i1 == nspan) { /* may be off end of array */
		  if(s[i].y == y0) {
		     break;
		  }
		  return(0);
	       }
	       
	       i = i0;
	       break;
	    }
	    i = (i0 + i1)/2;
	 }

	 do {				/* search all spans with y == y0 */
	    if(x0 >= s[i].x1 && x0 <= s[i].x2) {
	       return 1;			/* Ah! It is in mask */
	    }
	 } while(++i < nspan && s[i].y == y0);
      }
   }

   return 0;
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Extract all spanvecs in the CHAIN whose bounding boxes overlap the
 * given bounding box.
 */
CHAIN *
phFindSVChainBB(
		const CHAIN *chain,	/* array of OBJMASK */
		int x1,			/* boundaries of rectangle */
		int y1,
		int x2,
		int y2
		)
{
   CURSOR_T curs;			/* cursor for chain */
   OBJMASK *sv;
   CHAIN *oChain;
   char *cobjmask_type = (char *)objmask_type;

   shAssert(chain != NULL && chain->type == objmask_type);
   shAssert(x1 <= x2);
   shAssert(y1 <= y2);

   oChain = shChainNew(cobjmask_type);
   curs = shChainCursorNew(chain);
   while((sv = shChainWalk(chain,curs,NEXT)) != NULL) {
       if(x2 < sv->cmin || sv->cmax < x1 /* check bounding boxes */
	  || y2 < sv->rmin || sv->rmax < y1)
	   continue;
      (void)shChainElementAddByPos(oChain,sv,cobjmask_type,TAIL,AFTER);
   }
   shChainCursorDel(chain,curs);
   return oChain;
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Extract all spanvecs in the CHAIN whose bounding boxes overlap the
 * given circle.
 */
CHAIN *
phFindSVChainBBCircle(
		      const CHAIN *chain, /* array of OBJMASK */
		      int x,		/* centre of circle */
		      int y,
		      int rad		/* radius of circle */
		      )
{
   CURSOR_T curs;			/* cursor for chain */
   OBJMASK *sv;
   CHAIN *oChain;
   int rad2;
   char *cobjmask_type = (char *)objmask_type;

   shAssert(chain != NULL && chain->type == objmask_type);

   rad2 = rad*rad;
   oChain = shChainNew(cobjmask_type);
   curs = shChainCursorNew(chain);
   while((sv = shChainWalk(chain,curs,NEXT)) != NULL) {
       int dx, dx1, dy, dy1;
       int dist2;
       
       dx = sv->cmin - x;
       dx1 = x - sv->cmax;
       if(dx > 0) {			/* OBJMASK is to left of centre */
	  dist2 = dx*dx;
       } else if(dx1 > 0) {		/* OBJMASK is to right of centre */
	  dist2 = dx1*dx1;
       } else {				/* OBJMASK includes center */
	  dist2 = 0;
       }
       
       dy = sv->rmin - y;		/* similarily in y */
       dy1 = y - sv->rmax;
       if(dy > 0) {
	  dist2 += dy*dy;
       } else if(dy1 > 0) {
	  dist2 += dy1*dy1;
       } else {
	  ;
       }
       
       if(dist2 > rad2)
	   continue;
       
      (void)shChainElementAddByPos(oChain,sv,cobjmask_type,TAIL,AFTER);
   }
   shChainCursorDel(chain,curs);
   return oChain;
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Trim the mask so that it fits inside the rectangle.
 */
CHAIN *
phTrimMaskToRect (
		  const CHAIN *iMask,	/* array of OBJMASK */
		  int rx1,		/* boundaries of rectangle */
		  int ry1,
		  int rx2,
		  int ry2
		  )
{
   CHAIN *oMask;			/* output mask */
   OBJMASK *isv;			/* spanvecs from input mask */
   OBJMASK *osv;			/* spanvecs from output mask */
   int i, j;
   int inew;
   int nmask;				/* number of masks in iMask */
   int npix;				/* number of pixels in osv */
   char *cobjmask_type = (char *)objmask_type;

   shAssert(iMask != NULL && iMask->type == objmask_type);
   shAssert(rx1 <= rx2);
   shAssert(ry1 <= ry2);

   oMask = shChainNew(cobjmask_type);
   nmask = shChainSize(iMask);
   for(j = 0; j < nmask; j++) {
      isv = shChainElementGetByPos(iMask,j);
      if(rx2 < isv->cmin || isv->cmax < rx1 ||
					  ry2 < isv->rmin || isv->rmax < ry1) {
	 continue;			/* OBJMASK is outside rectangle */
      }

      osv = phObjmaskNew(isv->nspan);
      npix = inew = 0;
      for(i = 0; i < isv->nspan; i++) {
	 int y, x1, x2;			/* unaliased from isv->s[i] */
	 
	 y = isv->s[i].y;
	 if(y < ry1 || y > ry2) {	/* span's above or below */
	    continue;
	 }

	 x1 = isv->s[i].x1;
	 x2 = isv->s[i].x2;
	 if(x2 < rx1 || x1 > rx2) {	/* span's to left or right */
	    continue;
	 }

	 osv->s[inew].y = y;		/* trim span to fit */
	 x1 = (x1 < rx1) ? rx1 : x1;
	 x2 = (x2 > rx2) ? rx2 : x2;
	 osv->s[inew].x1 = x1;
	 osv->s[inew].x2 = x2;
	 npix += (x2 - x1 + 1);
	 inew++;
      }
      if(inew == 0) {
	 phObjmaskDel(osv);		/* There was, in fact, no overlap */
      } else {
	 osv->nspan = inew;
	 osv->cmin = (isv->cmin > rx1) ? isv->cmin : rx1;
	 osv->cmax = (isv->cmax < rx2) ? isv->cmax : rx2;
	 osv->rmin = (isv->rmin > ry1) ? isv->rmin : ry1;
	 osv->rmax = (isv->rmax < ry2) ? isv->rmax : ry2;
	 osv->npix = npix;

	 (void)shChainElementAddByPos(oMask,osv,cobjmask_type,TAIL,AFTER);
      }
   }

   return oMask;
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Return an OBJMASK that is the union of two OBJMASKs; the second mask
 * may be NULL
 *
 * See also phObjmaskMerge() and phObjmaskOrObjmask()
 */
OBJMASK *
phObjmaskUnion(const OBJMASK *om1,
	       const OBJMASK *om2)
{
   int i1, i2, j;
   OBJMASK *uni;			/* the desired union */
   int y;

   shAssert(om1 != NULL);
   if(om2 == NULL) {
      return(phObjmaskCopy(om1, 0, 0));
   }
   shAssert(om1->col0 == om2->col0 && om1->row0 == om2->row0);

   uni = phObjmaskNew(om1->nspan + om2->nspan);

   i1 = i2 = j = 0;
   while(i1 < om1->nspan) {
      y = om1->s[i1].y;
      while(i2 < om2->nspan && om2->s[i2].y <= y) {
	 uni->s[j++] = om2->s[i2++];
      }
      
      uni->s[j++] = om1->s[i1++];
   }
   while(i2 < om2->nspan) {
      uni->s[j++] = om2->s[i2++];
   }

   uni->nspan = j;
   shAssert(uni->size == j);
   phCanonizeObjmask(uni,1);

   return(uni);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Construct the mask which is the intersection of two masks.
 * This is currently an order (sizeof(mask1)*sizeof(mask2)) operation.
 * We may need to be more clever about this.
 */
CHAIN *
phMaskIntersection(
		const CHAIN *mask1, 
		const CHAIN *mask2
		)
{
   CURSOR_T curs1;			/* cursor for chain */
   CURSOR_T curs2;			/* cursor for chain */
   OBJMASK *sv1;
   OBJMASK *sv2;
   OBJMASK *svnew;
   CHAIN *imask;			/* intersection mask */
   char *cobjmask_type = (char *)objmask_type;
    
   shAssert(mask1 != NULL && mask1->type == objmask_type);
   shAssert(mask2 != NULL && mask2->type == objmask_type);

   curs1 = shChainCursorNew(mask1);
   curs2 = shChainCursorNew(mask2);
   imask = shChainNew(cobjmask_type);

   while((sv1 = shChainWalk(mask1,curs1,NEXT)) != NULL) {
      int col0 = sv1->col0, row0 = sv1->row0;
      while((sv2 = shChainWalk(mask2,curs2,NEXT)) != NULL) {
	 int irow1, irow2;
	 int inew;
	 SPAN *s1, *s2;
	 
	 shAssert(sv2->col0 == col0 && sv2->row0 == row0);
/*
 * check bounding box
 */
	 if(sv2->cmax < sv1->cmin || sv1->cmax < sv2->cmin ||
	    sv2->rmax < sv1->rmin || sv1->rmax < sv2->rmin) {
	    continue;
	 }

	 svnew = phObjmaskNew(sv1->nspan > sv2->nspan ?
			      sv1->nspan : sv2->nspan);
	 s1 = sv1->s;
	 s2 = sv2->s;
	 
	 irow1 = 0;
	 irow2 = 0;
	 inew = 0;
	 while(irow1 < sv1->nspan && irow2 < sv2->nspan) {
	    int icol1, icol2;
/*
 * find a span with a common y
 */
	    if(s1[irow1].y < s2[irow2].y) {
	       irow1++;
	       continue;
	    }
	    if(s2[irow2].y < s1[irow1].y) {
	       irow2++;
	       continue;
	    }
/*
 * loop over spans with this y
 */
	    for(icol1 = irow1, icol2 = irow2;
		icol1 < sv1->nspan && s1[icol1].y == s1[irow1].y
		&& s2[icol2].y == s2[irow2].y; icol1++) {
	       
	       while(icol2 < sv2->nspan &&
		     s2[icol2].y == s2[irow2].y &&
		     s2[icol2].x1 <= s1[icol1].x2) {
		  if(s2[icol2].x2 >= s1[icol1].x1) {
		     svnew->s[inew].y = s1[icol1].y;
		     svnew->s[inew].x1 = (s1[icol1].x1 > s2[icol2].x1) ? 
		       s1[icol1].x1 : s2[icol2].x1;
		     svnew->s[inew].x2 = (s1[icol1].x2 < s2[icol2].x2) ? 
		       s1[icol1].x2 : s2[icol2].x2;
		     inew++;
		  }
		  icol2++;
	       }
	    }
/*
 * No intersection at this y value; proceed to next spans
 */
	    irow1 = icol1; irow2 = icol2;
	 }
	 if(inew > 0) {
	    int k, npix;
	    int x1, x2;			/* unpacked from svnew */
					/* Find new bounding box */
	    svnew->nspan = inew;
	    svnew->cmin = MAXSHORT;
	    svnew->cmax = -MAXSHORT;
	    svnew->rmin = svnew->s[0].y;
	    svnew->rmax = svnew->s[inew-1].y;
	    npix = 0;
	    for(k = 0; k < inew; k++) {
	       x1 = svnew->s[inew].x1;
	       x2 = svnew->s[inew].x2;
	       if(x1 < svnew->cmin) svnew->cmin = x1;
	       if(x2 > svnew->cmax) svnew->cmax = x2;
	       npix += (x2 - x1 + 1);
	    }
	    svnew->npix = npix;
	    for(k = 0; k < inew; k++) {
	       svnew->s[k].y -= svnew->rmin;
	       svnew->s[k].x1 -= svnew->cmin;
	       svnew->s[k].x2 -= svnew->cmin;
	    }
	    
	    (void)shChainElementAddByPos(imask,svnew,cobjmask_type,TAIL,
					 AFTER);
	 } else {
	    phObjmaskDel(svnew);
	 }
      }
   }
   shChainCursorDel(mask1, curs1);
   shChainCursorDel(mask2, curs2);
   return imask;
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Replace sv1 with the intersection of sv1 and sv2, i.e. sv1 &= sv2;
 * see also phObjmaskIntersection(), phObjmaskNotIntersectionObjmask(),
 * and phObjmaskAndObjmask().
 *
 * The mask sv1 is returned.
 *
 * This is currently an order sizeof(mask2) operation.
 * We may need to be cleverer about this.
 */
OBJMASK *
phObjmaskAndObjmask(OBJMASK *sv1,	/* mask to clip */
		    const OBJMASK *sv2)	/* only keep pixels in this mask */
{
   int col0_1, row0_1;			/* == sv1->{col,row}0 */
   int cmin_1, cmax_1, rmin_1, rmax_1;	/* == sv1->{c,r}{min,max} */
   int col0_2, row0_2;
   int irow1, irow2;
   int inew;				/* index of new spans */
   int npix;				/* number of pixels in svnew */
   SPAN *s1;				/* == sv1->s */
   const SPAN *s2;			/* == sv2->s */
   OBJMASK *svnew;
   int x1, x2;				/* limits of overlapping span */
    
   shAssert(sv1 != NULL && sv2 != NULL);

   col0_1 = sv1->col0; row0_1 = sv1->row0;
   cmin_1 = sv1->cmin; cmax_1 = sv1->cmax;
   rmin_1 = sv1->rmin; rmax_1 = sv1->rmax;
   col0_2 = sv2->col0; row0_2 = sv2->row0;
/*
 * check bounding box
 */
   if(sv2->cmax + col0_2 < cmin_1 + col0_1 ||
					cmax_1 + col0_1 < sv2->cmin + col0_2 ||
      sv2->rmax + row0_2 < rmin_1 + row0_1 ||
					rmax_1 + row0_1 < sv2->rmin + row0_2) {
      return(phObjmaskRealloc(sv1, 0));
   }
   
   s1 = sv1->s;
   s2 = sv2->s;
   npix = inew = 0;
   svnew = phObjmaskNew(sv1->nspan);
   
   irow1 = 0;
   irow2 = 0;
   while(irow1 < sv1->nspan && irow2 < sv2->nspan) {
      int icol1, icol2;
/*
 * find a span with a common y
 */
      if(s1[irow1].y + row0_1 < s2[irow2].y + row0_2) {
	 irow1++;
	 continue;
      }
      if(s2[irow2].y + row0_2 < s1[irow1].y + row0_1) {
	 irow2++;
	 continue;
      }
				/* loop over spans with this y */
      for(icol1 = irow1, icol2 = irow2;
	  icol1 < sv1->nspan && s1[icol1].y == s1[irow1].y; icol1++) {
	 for(icol2 = irow2; icol2 < sv2->nspan && s2[icol2].y == s2[irow2].y &&
	     s2[icol2].x1 + col0_2 <= s1[icol1].x2 + col0_1; icol2++) {
	    if(s2[icol2].x2 + col0_2 >= s1[icol1].x1 + col0_1) {
	       /* the spans overlap */
	       if(s2[icol2].x1 + col0_2 <= s1[icol1].x1 + col0_1) {
		  x1 = s1[icol1].x1;
	       } else {
		  x1 = s2[icol2].x1 + col0_2 - col0_1;
	       }
	       if(s2[icol2].x2 + col0_2 < s1[icol1].x2 + col0_1) {
		  x2 = s2[icol2].x2 + col0_2 - col0_1;
	       } else {
		  x2 = s1[icol1].x2;
	       }
	       
	       if(inew == svnew->size) {
		  phObjmaskRealloc(svnew, inew*1.4 + 2);
	       }
	       svnew->s[inew].y = s1[icol1].y;
	       svnew->s[inew].x1 = x1;
	       svnew->s[inew].x2 = x2;
	       npix += (x2 - x1 + 1);
	       inew++;
	    }
	 }
      }
/*
 * No intersection at this y value; proceed to next spans
 */
      irow1 = icol1; irow2 = icol2;
   }

   if(inew == 0) {			/* no intersection */
      phObjmaskDel(svnew);
      return(phObjmaskRealloc(sv1, 0));
   } else {
      svnew->nspan = inew;
      svnew->npix = npix;
      svnew->col0 = sv1->col0; svnew->row0 = sv1->row0;
      phCanonizeObjmask(svnew,0);

      shFree(sv1->s);			/* copy svnew's spans to sv1 */
      sv1->s = svnew->s;
      sv1->npix = svnew->npix;
      sv1->nspan = svnew->nspan;
      sv1->size = svnew->size;
      svnew->s = NULL;
      phObjmaskDel(svnew);
   }

   return(sv1);
}

/*****************************************************************************/
/*
 * This function does the real work for phObjmaskIntersection() and
 * phObjmaskIntersectMask(). If intersect is NULL, only return a
 * 1 for intersection 0 otherwise; otherwise find the intersection
 * of the OBJMASK and CHAIN and return it as *intersect
 *
 * If returned, *intersect has row0, col0 taken from sv1 rather than mask2
 *
 * This is currently an order sizeof(mask2) operation.
 * We may need to be more clever about this.
 */
static int
find_mask_intersection(const OBJMASK *sv1, 
		       const CHAIN *mask2,
		       OBJMASK **intersect
		       )
{
   int col0_1, row0_1;			/* == sv1->{col,row}0 */
   int cmin_1, cmax_1, rmin_1, rmax_1;	/* == sv1->{c,r}{min,max} */
   int i;
   int inew;				/* index of new spans */
   int nmask2;				/* number of masks on mask2 */
   int npix;				/* number of pixels in svnew */
   OBJMASK *sv2;
   OBJMASK *svnew;
   int x1, x2;				/* limits of overlapping span */
    
   shAssert(sv1 != NULL);
   shAssert(mask2 != NULL && mask2->type == objmask_type);

   col0_1 = sv1->col0; row0_1 = sv1->row0;
   cmin_1 = sv1->cmin; cmax_1 = sv1->cmax;
   rmin_1 = sv1->rmin; rmax_1 = sv1->rmax;

   svnew = phObjmaskNew(intersect == NULL ? 0 : sv1->nspan);
   npix = inew = 0;

   nmask2 = shChainSize(mask2);
   for(i = 0; i < nmask2; i++) {
      int col0_2, row0_2;
      int irow1, irow2;
      SPAN *s1, *s2;
      sv2 = shChainElementGetByPos(mask2,i);
/*
 * check bounding box
 */
      if(sv2->cmax + sv2->col0 < cmin_1 + col0_1 ||
	 cmax_1 + col0_1 < sv2->cmin + sv2->col0 ||
	 sv2->rmax + sv2->row0 < rmin_1 + row0_1 ||
	 rmax_1 + row0_1 < sv2->rmin + sv2->row0)
	 continue;
      
      s1 = sv1->s;
      s2 = sv2->s;
      col0_2 = sv2->col0; row0_2 = sv2->row0;
      
      irow1 = 0;
      irow2 = 0;
      while(irow1 < sv1->nspan && irow2 < sv2->nspan) {
	 int icol1, icol2;
	 
	 /* find a span with a common y */
	 if(s1[irow1].y + row0_1 < s2[irow2].y + row0_2) {
	    irow1++;
	    continue;
	 }
	 if(s2[irow2].y + row0_2 < s1[irow1].y + row0_1) {
	    irow2++;
	    continue;
	 }
				/* loop over spans with this y */
	 for(icol1 = irow1, icol2 = irow2;
	     icol1 < sv1->nspan && s1[icol1].y == s1[irow1].y; icol1++) {
	    for(icol2 = irow2; icol2 < sv2->nspan && s2[icol2].y == s2[irow2].y
		&& s2[icol2].x1 + col0_2 <= s1[icol1].x2 + col0_1; icol2++) {
	       if(s2[icol2].x2 + col0_2 >= s1[icol1].x1 + col0_1) {
		  /* the spans overlap */
		  if(s2[icol2].x1 + col0_2 <= s1[icol1].x1 + col0_1) {
		     x1 = s1[icol1].x1;
		  } else {
		     x1 = s2[icol2].x1 + col0_2 - col0_1;
		  }
		  if(s2[icol2].x2 + col0_2 < s1[icol1].x2 + col0_1) {
		     x2 = s2[icol2].x2 + col0_2 - col0_1;
		  } else {
		     x2 = s1[icol1].x2;
		  }

		  if(inew == svnew->size) {
		     if(intersect == NULL) { /* we just want the return value*/
			phObjmaskDel(svnew);
			return(1);
		     }
		     phObjmaskRealloc(svnew, inew*1.4 + 2);
		  }
		  svnew->s[inew].y = s1[icol1].y;
		  svnew->s[inew].x1 = x1;
		  svnew->s[inew].x2 = x2;
		  npix += (x2 - x1 + 1);
		  inew++;
	       }
	    }
	 }
/*
 * No intersection at this y value; proceed to next spans
 */
	 irow1 = icol1; irow2 = icol2;
      }
   }

   if(intersect != NULL) {
      if(inew == 0) {			/* no intersection */
	 phObjmaskDel(svnew);
	 svnew = NULL;
      } else {
	 svnew->nspan = inew;
	 svnew->npix = npix;
	 svnew->col0 = sv1->col0; svnew->row0 = sv1->row0;
	 phCanonizeObjmask(svnew,0);
      }
      *intersect = svnew;
   } else {
      shAssert(inew == 0);
      phObjmaskDel(svnew);		/* mask2 and sv1 don't intersect */
   }
   
   return(inew == 0 ? 0 : 1);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Return 1 if OBJMASK and the second OBJMASK chain have pixels in common,
 * else 0.
 */
int
phObjmaskIntersectMask(const CHAIN *chain, /* array of OBJMASK */
		       const OBJMASK *sv) /* the OBJMASK in question */
{
   return(find_mask_intersection(sv, chain, NULL));
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Construct the OBJMASK consisting of all pixels in the first OBJMASK
 * and also in the second OBJMASK chain.
 *
 * See also phObjmaskIntersectionChain() and phObjmaskAndObjmask()
 */
OBJMASK *
phObjmaskIntersection(const OBJMASK *sv1, 
		      const CHAIN *mask2)
{
   OBJMASK *intersect = NULL;
   
   (void)find_mask_intersection(sv1, mask2, &intersect);

   return(intersect);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Return an OBJMASK chain consisting of all the OBJMASKs in mask2 that
 * intersect with om1; note that we don't copy the OBJMASKs, merely
 * return pointers (but with their reference counts incremented)
 *
 * This is currently an order sizeof(mask2) operation.
 * We may need to be more clever about this.
 *
 * See also phObjmaskIntersection()
 */
CHAIN *
phObjmaskIntersectionChain(const OBJMASK *om1, 
			   const CHAIN *mask2)
{
   int i;
   int nmask2;				/* number of masks on mask2 */
   CHAIN *result = NULL;		/* the desired result */
   OBJMASK *om2;
    
   shAssert(om1 != NULL);
   shAssert(mask2 != NULL && mask2->type == objmask_type);

   nmask2 = shChainSize(mask2);
   for(i = 0; i < nmask2; i++) {
      om2 = shChainElementGetByPos(mask2,i);

      if(phObjmaskIntersect(om1, om2,0,0)) {
	 if(result == NULL) {
	    result = shChainNew((char *)objmask_type);
	 }
	 (void)shChainElementAddByPos(result, om2, (char *)objmask_type,
				      TAIL, AFTER);
	 om2->refcntr++;
      }
   }
   
   return(result);
}

/*****************************************************************************/
/*
 * here's the workhorse for phObjmaskNotIntersectionObjmask
 *
 * Note that it does _not_ set its result's bounding box
 */
static void
i_phObjmaskNotIntersectionObjmask(OBJMASK *svnew, 
				  const OBJMASK *sv2)
{
   int irow1, irow2;
   int inew;
   int ndel;				/* number of spans deleted */
   SPAN *s1, *s2;
   int y;				/* the y that we are examining */
/*
 * check bounding box
 */
   if(sv2->cmax < svnew->cmin || svnew->cmax < sv2->cmin ||
      sv2->rmax < svnew->rmin || svnew->rmax < sv2->rmin) {
      return;
   }
/*
 * OK, there may be some (partial?) spans to remove from svnew
 */
   ndel = 0;
      
   s1 = svnew->s;
   s2 = sv2->s;
      
   irow1 = 0;
   irow2 = 0;
   inew = svnew->nspan;
   while(irow1 < svnew->nspan && irow2 < sv2->nspan) {
      int icol1, icol2;
      /* find a span with a common y */
      if(s1[irow1].y < s2[irow2].y) {
	 irow1++;
	 continue;
      }
      if(s2[irow2].y < s1[irow1].y) {
	 irow2++;
	 continue;
      }
      /* loop over spans with this y */
      y = s1[irow1].y;
      shAssert(s2[irow2].y == y);
      for(icol1 = irow1, icol2 = irow2; icol1 < svnew->nspan &&
			       s1[icol1].y == y && s2[icol2].y == y; icol1++) {
	 
	 while(icol2 < sv2->nspan &&
			    s2[icol2].y == y && s2[icol2].x1 <= s1[icol1].x2) {
	    
	    if(s2[icol2].x2 >= s1[icol1].x1) {
	       /* We've got some kind of overlap, find */
	       /* out what it is. */
	       if(s2[icol2].x1 <= s1[icol1].x1) {
		  if(s2[icol2].x2 < s1[icol1].x2) { /* overlap on left */
		     s1[icol1].x1 = s2[icol2].x2 + 1;
		  } else {		/* complete overlap: delete span */
		     s1[icol1].y = -MAXSHORT + 1;
		     ndel++;
		     break;		/* don't increment icol2 */
		  }
	       } else {
		  if(s2[icol2].x2 >= s1[icol1].x2) { /* overlap on right */
		     s1[icol1].x2 = s2[icol2].x1 - 1;
		     break;		/* don't increment icol2 */
		  } else {		/* in the middle: a new span */
		     if(inew == svnew->size) {
			phObjmaskRealloc(svnew, inew*1.4 + 2);
			s1 = svnew->s;
		     }
		     svnew->s[inew] = s1[icol1];
		     s1[icol1].x1 = s2[icol2].x2 + 1;
		     s1[inew].x2 = s2[icol2].x1 - 1;
		     inew++;
		  }
	       }
	    }
	    icol2++;
	 }
      }
/*
 * No intersection at this y value; proceed to next spans
 */
      irow1 = icol1; irow2 = icol2;
   }
   
   if(ndel > 0 || inew > svnew->nspan) {
      svnew->nspan = inew;
      phCanonizeObjmask(svnew,1);
   }
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Set all pixels in the second OBJMASK in the first
 * (see also phObjmaskUnion()). The mask om1 is returned.
 */
OBJMASK *
phObjmaskOrObjmask(OBJMASK *om1,
		   const OBJMASK *om2)
{
   int i;
   int nspan1, nspan2;
   SPAN *s1, *s2;			/* == om[12]->s */

   shAssert(om1 != NULL);
   if(om2 == NULL) {
      return(om1);
   }
   nspan1 = om1->nspan; nspan2 = om2->nspan;
   
   phObjmaskRealloc(om1, nspan1 + nspan2);
   om1->nspan += nspan2;

   s1 = om1->s; s2 = om2->s;
   for(i = 0; i < nspan2; i++) {
      s1[nspan1 + i] = s2[i];
   }

   phCanonizeObjmask(om1, 1);

   return(om1);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Remove all pixels that are set in the second OBJMASK from the first
 * (see also phObjmaskNotIntersectionObjmask() and phObjmaskAndObjmask()).
 *
 * The mask sv1 is returned.
 *
 * This is currently an order sizeof(sv2) operation.
 * We may need to be cleverer about this.
 */
OBJMASK *
phObjmaskAndNotObjmask(OBJMASK *sv1, 
		       const OBJMASK *sv2)
{
   shAssert(sv1 != NULL);
   shAssert(sv2 == NULL || (sv2->col0 == sv1->col0 && sv2->row0 == sv1->row0));
/*
 * see if there are some spans to remove
 */
   if(sv2 != NULL) {
      i_phObjmaskNotIntersectionObjmask(sv1, sv2);
      phObjmaskBBSet(sv1);
   }
   
   return(sv1);
}

/*
 * <AUTO EXTRACT>
 *
 * Construct and return the OBJMASK consisting of all pixels in the
 * first OBJMASK, but NOT in the second OBJMASK (which may be NULL).
 * (See also phObjmaskAndNotObjmask() and phObjmaskAndObjmask())
 *
 * This is currently an order sizeof(sv2) operation.
 * We may need to be more clever about this.
 */
OBJMASK *
phObjmaskNotIntersectionObjmask(const OBJMASK *sv1, 
				const OBJMASK *sv2)
{
   int i;
   OBJMASK *svnew;
    
   shAssert(sv1 != NULL);
   shAssert(sv2 == NULL || (sv2->col0 == sv1->col0 && sv2->row0 == sv1->row0));
/*
 * start by simply copying sv1
 */
   svnew = phObjmaskNew(sv1->nspan);
   for(i = 0; i < sv1->nspan; i++) {
      svnew->s[i] = sv1->s[i];
   }
   svnew->nspan = sv1->nspan;
   svnew->npix = sv1->npix;
   svnew->rmin = sv1->rmin; svnew->rmax = sv1->rmax;
   svnew->cmin = sv1->cmin; svnew->cmax = sv1->cmax;
/*
 * now see if there are some spans to remove
 */
   if(sv2 != NULL) {
      i_phObjmaskNotIntersectionObjmask(svnew, sv2);
   }
   
   phObjmaskBBSet(svnew);
   return(svnew);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Construct and return the OBJMASK consisting of all pixels in the
 * first OBJMASK, but NOT in the second OBJMASK chain.
 *
 * This is currently an order sizeof(mask2) operation.
 * We may need to be more clever about this.
 */
OBJMASK *
phObjmaskNotIntersectionObjmaskChain(const OBJMASK *sv1, 
				     const CHAIN *mask2)
{
   int col0, row0;
   CURSOR_T curs2;			/* cursor for chain */
   const OBJMASK *sv2;
   OBJMASK *svnew;
   int i;
    
   shAssert(sv1 != NULL);
   shAssert(mask2 != NULL && mask2->type == objmask_type);

   col0 = sv1->col0;
   row0 = sv1->row0;
   curs2 = shChainCursorNew(mask2);

   svnew = phObjmaskNew(sv1->nspan);
   for(i = 0; i < sv1->nspan; i++) {
      svnew->s[i] = sv1->s[i];
   }
   svnew->nspan = sv1->nspan;
   svnew->npix = sv1->npix;
   svnew->rmin = sv1->rmin; svnew->rmax = sv1->rmax;
   svnew->cmin = sv1->cmin; svnew->cmax = sv1->cmax;

   while((sv2 = shChainWalk(mask2,curs2,NEXT)) != NULL) {
      shAssert(sv2->col0 == col0 && sv2->row0 == row0);

      i_phObjmaskNotIntersectionObjmask(svnew, sv2);
   }

   shChainCursorDel(mask2, curs2);
   
   phObjmaskBBSet(svnew);
   return(svnew);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Construct the mask consisting of all pixels in the first mask, but
 * NOT in the second mask.
 * This is currently an order (sizeof(mask1)*sizeof(mask2)) operation.
 * We may need to be more clever about this.
 */
CHAIN *
phMaskNotIntersection(
		   const CHAIN *mask1, 
		   const CHAIN *mask2
		   )
{
   CURSOR_T curs1;			/* cursor for chain1 */
   OBJMASK *sv1;
   OBJMASK *svnew;
   CHAIN *imask;			/* mask of pixels in 1 but not in 2 */
   char *cobjmask_type = (char *)objmask_type;
    
   shAssert(mask1 != NULL && mask1->type == objmask_type);
   shAssert(mask2 != NULL && mask2->type == objmask_type);

   curs1 = shChainCursorNew(mask1);

   imask = shChainNew(cobjmask_type);
   while((sv1 = shChainWalk(mask1,curs1,NEXT)) != NULL) {
      svnew = phObjmaskNotIntersectionObjmaskChain(sv1, mask2);

      if(svnew->nspan > 0) {
	 (void)shChainElementAddByPos(imask,svnew,cobjmask_type,TAIL, AFTER);
      } else {
	 phObjmaskDel(svnew);
      }
   }
   shChainCursorDel(mask1, curs1);

   return imask;
}
#endif /* !defined(STAND_ALONE) */

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Set pixels in a MASK from an OBJMASK
 */
void
phMaskSetFromObjmask(const OBJMASK *om, /* the OBJMASK */
		     MASK *mask,	/* mask to set */
		     const int val)	/* value to OR into mask */
{
   char cval = val;			/* unpacked from val */
   int i;
   int nspan;				/* unpack obj->nspan for compiler */
   SPAN *spans;				/* SPANs in this OBJMASK */
   int col0, row0, ncol, nrow;		/* unpacked from mask (and corrected
					   for om->{row,col}0) */
   int y, x1, x2;			/* unpacked from a SPAN */
   unsigned char **mask_rows;		/* unpacked from obj1->mask->rows */

   shAssert(om != NULL && mask != NULL);
/*
 * Examine list of spans, setting bits in the MASK.
 */
   mask_rows = mask->rows;
   nrow = mask->nrow; ncol = mask->ncol;
   row0 = mask->row0 - om->row0; col0 = mask->col0 - om->col0;
   
   nspan = om->nspan;
   spans = om->s;
   
   for(i = 0;i < nspan;i++) {
      int j;
      unsigned char *mask_row;
      
      y = spans[i].y - row0;
      if(y < 0 || y >= nrow) {
	 continue;
      }
      x1 = spans[i].x1 - col0; x2 = spans[i].x2 - col0;
      if(x1 < 0) {
	 x1 = 0;
      }
      if(x2 >= ncol) {
	 x2 = ncol - 1;
      }
      
      mask_row = mask_rows[y];
      for(j = x1;j <= x2;j++) {
	 mask_row[j] |= cval;
      }
   }
}

#if !defined(STAND_ALONE)

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Set pixels in a MASK from a CHAIN of OBJMASK
 */
void
phMaskSetFromObjmaskChain(
			  const CHAIN *chain, /* array of OBJMASK */
			  MASK *mask,	/* mask of REGION objs were found in */
			  const int val	/* value to OR into mask */
			 )
{
   char cval = val;			/* unpacked from val */
   CURSOR_T curs;			/* cursor for chain */
   int i;
   int nspan;				/* unpack obj->nspan for compiler */
   SPAN *spans;				/* SPANs in this OBJMASK */
   int col0, row0, ncol, nrow;		/* unpacked from mask (and corrected
					   for sv->{row,col}0) */
   int y, x1, x2;			/* unpacked from a SPAN */
   unsigned char **mask_rows;		/* unpacked from obj1->mask->rows */
   OBJMASK *sv;

   shAssert(chain != NULL && chain->type == objmask_type);
   shAssert(mask != NULL);
/*
 * Examine list of spans, setting the parent REGION's MASK.
 */

   mask_rows = mask->rows;
   nrow = mask->nrow; ncol = mask->ncol;
   
   curs = shChainCursorNew(chain);
   while((sv = shChainWalk(chain,curs,NEXT)) != NULL) {
      
      row0 = mask->row0 - sv->row0; col0 = mask->col0 - sv->col0;
      nspan = sv->nspan;
      spans = sv->s;
      
      for(i = 0;i < nspan;i++) {
	 int j;
	 unsigned char *mask_row;
	 
	 y = spans[i].y - row0;
	 if(y < 0 || y >= nrow) {
	    continue;
	 }
	 x1 = spans[i].x1 - col0; x2 = spans[i].x2 - col0;
	 if(x1 < 0) {
	    x1 = 0;
	 }
	 if(x2 >= ncol) {
	    x2 = ncol - 1;
	 }
	 
	 mask_row = mask_rows[y];
	 for(j = x1;j <= x2;j++) {
	    mask_row[j] |= cval;
	 }
      }
   }
   shChainCursorDel(chain,curs);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Given a possible non-simply-connected objmask, return an array of the
 * simply connected components (followed by a NULL)
 *
 * Usually a set of new masks will be allocated, but, as a special case,
 * if alloc_one is false and the initial mask is simply connected,
 * and thus *nmask == 1, the routine'll return NULL (saving the expense
 * of allocating an array and copying om)
 */
OBJMASK **
phObjmaskSplitSimple(const OBJMASK *om,	/* initial OBJMASK */
		     int *pnmask,	/* number of resulting masks or NULL */
		     int alloc_one)	/* allocate a copy of om to return
					   if *nmask == 1? */
{
   int i, j, k;
   int *ids;				/* ID values for spans */
   int id;				/* ID of current span */
   int i0, i1;				/* range of spans with current y */
   int im0, im1;			/* range of spans with previous y */
   OBJMASK **masks = NULL;		/* return array */
   int next_id = 1;			/* next id to use */
   int nmask = 0;			/* number of masks; *pnmask == nmask */
   int nspan;				/* == om->nspan */
   int *nspans;				/* number of spans in each OBJMASK */
   SPAN *spans;				/* == om->s */
   int x1, x2;				/* unpacked from a span */
   int y;				/* the row of current interest */

   shAssert(om != NULL);
   nspan = om->nspan;
   spans = om->s;

   if(pnmask == NULL) pnmask = &nmask;
/*
 * handle trivial special case first; no spans in mask
 */
   if(nspan == 0) {
      *pnmask = nmask = 0;
      if(alloc_one) {
	 masks = shMalloc((nmask + 1)*sizeof(OBJMASK *));
	 masks[0] = NULL;
      }

      return(masks);
   }
/*
 * we have to do some work. Work through om seeing which spans touch,
 * and assigning each connected mask an id (we'll sometimes have to merge
 * apparently disconnected masks of course, the classic \Lambda shaped
 * object problem)
 */
   ids = alloca(nspan*sizeof(int));
   memset(ids, '\0', nspan*sizeof(int)); shAssert(ids[0] == 0);
   nspans = alloca((nspan + 1)*sizeof(int));
   y = spans[0].y - 2;
   i0 = 0; i1 = -1;			/* range of spans with s.y == y-1 */

   i = 0;
   while(i < nspan) {
      if(spans[i].y == y + 1) {		/* may touch previous spans */
	 im0 = i0; im1 = i1;
      } else {
	 im0 = 0; im1 = -1;
      }
      y = spans[i].y;

      for(i0 = i++; i < nspan && spans[i].y == y; i++) continue;
      i1 = i - 1;

      for(j = i0; j <= i1; j++) {
	 x1 = spans[j].x1; x2 = spans[j].x2;
	 id = ids[j];
	 for(k = im0; k <= im1; k++) {
	    if(spans[k].x1 <= x2 + 1 &&
	       spans[k].x2 >= x1 - 1) { /* spans touch */
	       if(id == 0) {		/* not already assigned */
		  id = ids[j] = ids[k];
		  nspans[id]++;
	       } else if(id == ids[k]) { /* consistent assignment */
		  ;
	       } else {			/* inconsistent assignment */
		  int l;
		  const int idk = ids[k];
		  for(l = 0; l < j; l++) {
		     if(ids[l] == idk) {
			ids[l] = id;
		     }
		  }
		  nspans[id] += nspans[idk];
		  nmask--;		/* we merged two */
		  nspans[idk] = 0;
	       }
	    }
	 }
	 if(id == 0) {			/* a new object */
	    ids[j] = next_id++;
	    nspans[ids[j]] = 1;
	    nmask++;
	 }
      }
   }
/*
 * Is this the special case of a simply connected OBJMASK?
 */
   if(nmask == 1) {
      *pnmask = nmask;
      if(alloc_one) {
	 masks = shMalloc((nmask + 1)*sizeof(OBJMASK *));
	 masks[0] = phObjmaskCopy(om, 0, 0);
	 masks[1] = NULL;
      }
      return(masks);
   }
/*
 * pack those spans into OBJMASKs of their own. Start by allocating the
 * desired OBJMASKs; note that the array masks is initially sparse if
 * we've had to merge OBJMASKs in the previous step
 */
   *pnmask = nmask;
   masks = shMalloc(next_id*sizeof(OBJMASK *));
   masks[0] = NULL;
   j = k = 0;
   for(i = 1; i < next_id; i++) {
      if(nspans[i] > 0) {
	 j++; k += nspans[i];		/* used in shAssertions */
	 masks[i] = phObjmaskNew(nspans[i]);
	 masks[i]->nspan = nspans[i];
	 masks[i]->row0 = om->row0; masks[i]->col0 = om->col0;
	 nspans[i] = 0;
      } else {
	 masks[i] = NULL;
      }
   }

   shAssert(j == nmask);
   shAssert(k == om->nspan);
/*
 * copy the spans into the proper places
 */
   for(i = 0; i < nspan; i++) {
      id = ids[i];
      masks[id]->s[nspans[id]++] = spans[i];
   }
/*
 * make the masks array compact, and set bounding boxes
 */
   for(i = j = 0; i < next_id; i++) {
      if(masks[i] != NULL) {
	 shAssert(masks[i]->nspan == nspans[i]);
	 phObjmaskBBSet(masks[i]);
	 masks[j++] = masks[i];
      }
   }
   masks[j] = NULL;

   return(masks);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Set a given pixel in a SPANMASK
 */
void
phSpanmaskSetAsPix(
		   SPANMASK *sm, /* SPANMASK to set */
		   int row,	/* pixel location */
		   int col,
		   const S_MASKTYPE val	/* type of mask */
		   )
{
    OBJMASK *om;
    char *cobjmask_type = (char *)objmask_type;
    
    shAssert(sm != NULL);
    shAssert(sm->cookie == SPAN_COOKIE);

    om = phObjmaskNew(1);
    om->npix = om->nspan = 1;
    om->s[0].y = row;
    om->s[0].x1 = col;
    om->s[0].x2 = col;
    om->cmin = col;
    om->cmax = col;
    om->rmin = row;
    om->rmax = row;
    (void)shChainElementAddByPos(sm->masks[val],om,cobjmask_type,TAIL, AFTER);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * Find how many pixels are set in an OBJMASK
 */
int phObjmaskNumPix(
		       const OBJMASK *mask
		       )
{
    int i, npix=0;
    SPAN *span;

    span = mask->s;
    for (i=0; i<mask->nspan; i++) {
	npix += span[i].x2 - span[i].x1 + 1;
    }
    return npix;
}
/*****************************************************************************/
/*
 * Here's an implementation of heapsort, specialised to sorting SPANs
 * in y and x
 *
 * Because SPAN.y and SPAN.x1 are both less than 16 bits, we can generate
 * a 32-bit sort key by putting y in the upper 16 bits, and x1 in the lower;
 * if the spans are sorted on this they are sorted first by y and then by x
 */
static void
heap_sort(SPAN *s, unsigned int n)
{
   unsigned int i, j, k, l;
   SPAN elem;
   
   if(n <= 1) return;
/*
 * Build the heap
 */
   l = n/2;
   k = n - 1;

   while(l > 0) {
      elem = s[--l];
      
      i = l;
      j = 2*l + 1;
      while(j <= k) {
	 if(j < k &&
	    SORTKEY(s[j].x1, s[j].y) < SORTKEY(s[j + 1].x1, s[j + 1].y)) {
	    j++;
	 }
	 if(SORTKEY(elem.x1, elem.y) < SORTKEY(s[j].x1, s[j].y)) {
	    s[i] = s[j];
	    i = j;
	    j = 2*j + 1;
	 } else {
	    j = k + 1;
	 }
      }
      s[i] = elem;
   }
/*
 * and destroy it again, resulting in a sorted array
 */
   for(;;) {
      elem = s[k];
      s[k] = s[0];
      if(--k == 0) {
	 s[0] = elem;
	 break;
      }
      
      i = 0;
      j = 1;
      while(j <= k) {
	 if(j < k &&
	    SORTKEY(s[j].x1, s[j].y) < SORTKEY(s[j + 1].x1, s[j + 1].y)) {
	    j++;
	 }
	 if(SORTKEY(elem.x1, elem.y) < SORTKEY(s[j].x1, s[j].y)) {
	    s[i] = s[j];
	    i = j;
	    j = 2*j + 1;
	 } else {
	    j = k + 1;
	 }
      }
      s[i] = elem;
   }
}

/*****************************************************************************/
/*
 * And here's an insertion sort. This is a Good Choice as the data to be
 * sorted is nearly sorted (and thus insertion is linear with a small
 * coefficient)
 */
static void
insertion_sort(SPAN *s, unsigned int n)
{
   int i, j;
   int nswap = 0;
   SPAN elem;
   SPAN s0;
   
   if(n <= 1) return;

   s0 = s[0];
   s[0].y = -MAXSHORT;			/* use as a sentinel */
   for(j = 1;j < n;j++) {
      elem = s[j];
      for(i = j - 1;SORTKEY(s[i].x1, s[i].y) > SORTKEY(elem.x1, elem.y);i--) {
	 s[i + 1] = s[i];
	 nswap++;
      }
      if(i == 0) {
	 if(SORTKEY(s0.x1, s0.y) > SORTKEY(elem.x1, elem.y)) {
	    s[1] = s0;
	    s0 = elem;
	    nswap++;
	    continue;
	 }
      }
      s[i + 1] = elem;
   }
   s[0] = s0;

#if 0					/* was this a good choice of sort? */
   if(nswap > n) {
      fprintf(stderr,"insertion_sort: n == %d > nswap == %d\n",n,nswap);
   }
#endif
}

/*****************************************************************************/
/*
 * See if an OBJMASK is canonical, i.e.
 *   Each span has x2 >= x1
 *   All the spans are sorted in y-then-x
 *   npix is correct
 */
int
phObjmaskIsCanonical(const OBJMASK *om)
{
   int i;
   long key;				/* key in sort; >= 32bits */
   int npix;				/* number of pixels in om */
   long okey;				/* old value of key */
   int nspan = om->nspan;		/* unalias nspan and s */
   SPAN *s = om->s;

   if(nspan <= 1) {
      return(1);
   }

   npix = s[0].x2 - s[0].x1 + 1;
   key = SORTKEY(s[0].x1, s[0].y);
   if(s[0].x2 < s[0].x1) {
      return(0);
   }
   
   for(i = 1;i < nspan;i++) {
      if(s[i].x2 < s[i].x1) {
	 return(0);
      }
      okey = key;
      key = SORTKEY(s[i].x1, s[i].y);
      if(key < okey) {
	 return(0);
      }
      npix += s[i].x2 - s[i].x1 + 1;
   }

   return(om->npix == npix ? 1 : 0);
}

/*****************************************************************************/
/*
 * See if a CHAIN of OBJMASKs is canonical, i.e.
 *   All the OBJMASKs have the same row0, col0
 *   The OBJMASKs rmin, cmin values are sorted
 *   None of the OBJMASKs overlap -- note that this is an N^2 check
 */
static int
is_canonical_objmaskChain(const CHAIN *chain)
{
   CURSOR_T crsr1, crsr2;
   long key;				/* key in sort; >= 32bits */
   long okey;				/* old value of key */
   int row0, col0;
   const OBJMASK *obj1, *obj2;
   
   shAssert(chain != NULL && chain->type == objmask_type);
   
   crsr1 = shChainCursorNew(chain);
   if((obj1 = shChainWalk(chain, crsr1, NEXT)) == NULL) {
      shChainCursorDel(chain,crsr1);
      return(1);
   }
   
   row0 = obj1->row0; col0 = obj1->col0;
   key = SORTKEY(obj1->cmin, obj1->rmin);
   
   crsr2 = shChainCursorNew(chain);
   while((obj1 = shChainWalk(chain, crsr1, NEXT)) != NULL) {
      okey = key;
      key = SORTKEY(obj1->cmin, obj1->rmin);
      if(obj1->row0 != row0 || obj1->col0 != col0) {
	 break;
      }
      if(key < okey) {
	 break;
      }

      shChainCursorSet(chain,crsr2,HEAD);
      while((obj2 = shChainWalk(chain, crsr2, NEXT)) != obj1) {
	 if(phObjmaskIntersect(obj1,obj2,0,0)) {
	    break;
	 }
      }
   }
   shChainCursorDel(chain, crsr1);
   shChainCursorDel(chain, crsr2);

#if 1
   if(obj1 != NULL) {
      shError("canonical test fails");
      obj1 = NULL;
   }
#endif

   return(obj1 == NULL ? 1 : 0);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Grow an OBJMASK to NxN super-pixel boundaries. For example, if the
 * superpixels are 4x4, the OBJMASK
 *  (y: x1, x2) = (1: 2, 3) (2: 3, 5) (3: 2, 4) (4: 1, 3)
 * would become
 *                (0: 0, 7) (1: 0, 7) (2: 0, 7) (3: 0, 7)
 *                (4: 0, 3) (5: 0, 3) (6: 0, 3) (7: 0, 3)
 * i.e. it covers all the superpixels that have at least one pixel in
 * the original OBJMASK
 *
 * This is indended for use in matching master_masks to the binned
 * version of the data region
 */
OBJMASK *
phObjmaskGrowToSuperpixel(const OBJMASK *om, /* OBJMASK to grow */
			  int n)	/* size of superpixels */
{
   int i,j,k;
   int nspan;				/* == om->nspan */
   OBJMASK *onew;			/* returned OBJMASK */
   OBJMASK *obinned;			/* om with binned coordinates */
   int x1, x2, y;

   shAssert(om != NULL);
   shAssert(phObjmaskIsCanonical(om));
   shAssert(n > 0);

/*
 * Start by making a copy of the initial OBJMASK, which we then bin; 
 * we'll later expand it again
 */
   obinned = phObjmaskCopy(om,0,0);
   
   nspan = om->nspan;
   for(i = 0;i < nspan;i++) {
      obinned->s[i].y /= n;
      obinned->s[i].x1 /= n;
      obinned->s[i].x2 /= n;
   }
   phCanonizeObjmask(obinned, 1);
/*
 * Each of those spans will expand to n spans when we unbin, so allocate
 * the objmask and expand obinned into it
 */
   onew = phObjmaskNew(n*obinned->nspan);

   nspan = obinned->nspan;
   for(i = j = 0;i < nspan;i++) {
      y = n*obinned->s[i].y;
      x1 = n*obinned->s[i].x1;
      x2 = n*obinned->s[i].x2 + n - 1;
      for(k = 0;k < n;k++) {
	 onew->s[j].y = y + k;
	 onew->s[j].x1 = x1;
	 onew->s[j++].x2 = x2;
      }
   }

   onew->nspan = n*obinned->nspan;
   onew->rmin = n*obinned->rmin;
   onew->cmin = n*obinned->cmin;
   onew->rmax = n*obinned->rmax + n - 1;
   onew->cmax = n*obinned->cmax + n - 1;

   shAssert(j == onew->nspan);
   phCanonizeObjmask(onew, 1);

   phObjmaskDel(obinned);

   return(onew);
}
#endif /* !defined(STAND_ALONE) */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "dervish.h"
#include "phConsts.h"

/*****************************************************************************/
/*
 * A routine to return this version of photo's name
 */
const char *
phPhotoVersion(void)
{
   static const char version[] = "$Name: v5_4_11 $";

   if(strlen(version) <= 9) {
      return("NOCVS");
   } else {
      return(version);
   }
}

/*****************************************************************************/
/*
 * utilities that dervish usually provides
 */
#if !defined(DERVISH_H)			/* we haven't got the _real_ dervish */

TYPE
shTypeGetFromName(const char *name)	/* NOTUSED */
{
   return(UNKNOWN);
}

/*****************************************************************************/

void
shError(char *fmt, ...)
{
   va_list args;
   char buff[1000];
   int len;

   va_start(args,fmt);
   vsprintf(buff,fmt,args);
   va_end(args);

   if(buff[len = strlen(buff)] == '\n') {
      buff[len] = '\0';
   }
   fprintf(stderr,"Error: %s\n",buff);
   fflush(stderr);
}

/*****************************************************************************/
/*
 * This is the same as shError for current purposes
 */
void
shErrStackPush(char *fmt, ...)
{
   va_list args;
   char buff[1000];
   int len;

   va_start(args,fmt);
   vsprintf(buff,fmt,args);
   va_end(args);

   if(buff[len = strlen(buff)] == '\n') {
      buff[len] = '\0';
   }
   shErrStackPush("%s",buff);
   fprintf(stderr,"Error: %s\n",buff);
   fflush(stderr);
}

/*
 * and here's the fatal handler
 */
void
shFatal(char *fmt, ...)
{
   va_list args;

   va_start(args,fmt);
   fprintf(stderr,"Fatal error: ");
   vfprintf(stderr,fmt,args);
   fprintf(stderr,"\n");
   fflush(stderr);
   va_end(args);
   abort();
}

/*****************************************************************************/
/*
 * memory management
 */
void *
shMalloc(size_t n)
{
   void *ptr = malloc(n);

   if(ptr == NULL) {
      shFatal("failed to allocate %ld bytes", (long)n);
   }

   return(ptr);
}

void *
shRealloc(void *ptr, size_t n)
{
   ptr = realloc(ptr, n);

   if(ptr == NULL) {
      shFatal("failed to reallocate %ld bytes", (long)n);
   }

   return(ptr);
}

void
shFree(void *ptr)
{
   free(ptr);
}

int
p_shMemRefCntrGet(void *ptr)		/* NOTUSED */
{
   return(0);
}

void
p_shMemRefCntrDecr(void *ptr)		/* NOTUSED */
{
   ;
}

/*****************************************************************************/
/*
 * regions
 */
REGION *
shRegNew(const char *name,		/* NOTUSED */
	 int nrow,
	 int ncol,
	 int type)
{
   int i;
   REGION *reg = shMalloc(sizeof(REGION));
   
   reg->type = type;
   reg->nrow = nrow; reg->ncol = ncol;

   reg->rows_u8 = NULL; reg->rows_s8 = NULL;
   reg->rows = reg->rows_u16 = NULL; reg->rows_s16 = NULL;
   reg->rows_u32 = NULL; reg->rows_s32 = NULL;
   reg->rows_fl32 = NULL;
   reg->mask = NULL;
   reg->row0 = reg->col0 = 0;
   
   switch (reg->type) {
    case TYPE_U16:
      reg->rows = reg->rows_u16 = shMalloc(nrow*sizeof(U16 *));
      reg->rows[0] = shMalloc(nrow*ncol*sizeof(U16));

      for(i = 1; i < nrow; i++) {
	 reg->rows[i] = reg->rows[i - 1] + ncol;
      }
      
      break;
    case TYPE_FL32:
      reg->rows_fl32 = shMalloc(nrow*sizeof(FL32 *));
      reg->rows_fl32[0] = shMalloc(nrow*ncol*sizeof(FL32));

      for(i = 1; i < nrow; i++) {
	 reg->rows_fl32[i] = reg->rows_fl32[i - 1] + ncol;
      }
      
      break;
    default:
      shFatal("Impossible type of REGION: %d", reg->type);
      break;				/* NOTREACHED */
   }
   
   return(reg);
}

void
shRegDel(REGION *reg)
{
   if(reg == NULL) {
      return;
   }
   
   switch (reg->type) {
    case TYPE_U16:
      if(reg->rows != NULL) {
	 shFree(reg->rows[0]);
	 shFree(reg->rows);
      }
      break;
    case TYPE_FL32:
      if(reg->rows_fl32 != NULL) {
	 shFree(reg->rows_fl32[0]);
	 shFree(reg->rows_fl32);
      }
      break;
    default:
      shFatal("Impossible type of REGION: %d", reg->type);
      break;				/* NOTREACHED */
   }

   shFree(reg);
}

/*****************************************************************************/
/*
 * masks
 */
MASK *
shMaskNew(const char *name,		/* NOTUSED */
	 int nrow,
	 int ncol)
{
   int i;
   MASK *mask = shMalloc(sizeof(MASK));
   
   mask->rows = shMalloc(nrow*sizeof(unsigned char *));
   mask->rows[0] = shMalloc(nrow*ncol);
   mask->nrow = nrow; mask->ncol = ncol;
   mask->row0 = mask->col0 = 0;

   for(i = 1; i < nrow; i++) {
      mask->rows[i] = mask->rows[i - 1] + ncol;
   }

   return(mask);
}

void
shMaskDel(MASK *mask)
{
   if(mask != NULL) {
      if(mask->rows != NULL) {
	 shFree(mask->rows[0]);
	 shFree(mask->rows);
      }
      shFree(mask);
   }
}

void
shMaskClear(MASK *mask)
{
   int i;
   
   shAssert(mask != NULL && mask->rows != NULL && mask->rows[0] != NULL);
   shAssert(mask->nrow >= 1);

   for(i = 0; i < mask->ncol; i++) {
      mask->rows[0][i] = '\0';
   }
   for(i = 1; i < mask->nrow; i++) {
      memcpy(mask->rows[i], mask->rows[0], mask->ncol);
   }
}

/*****************************************************************************/
/*
 * CHAINs
 */
CHAIN *
shChainNew(char *type)			/* NOTUSED */
{
   return(NULL);
}

void
shChainDel(CHAIN *ch)			/* NOTUSED */
{
   ;
}

void
shChainElementAddByPos(CHAIN *ch,	/* NOTUSED */
		       void *el,	/* NOTUSED */
		       char *type,	/* NOTUSED */
		       int w,		/* NOTUSED */
		       int how)		/* NOTUSED */
{
   ;
}

void *
shChainElementGetByPos(const CHAIN *ch,	/* NOTUSED */
		       int el)		/* NOTUSED */
{
   return(NULL);
}

void *
shChainElementRemByPos(const CHAIN *ch,	/* NOTUSED */
		       int el)		/* NOTUSED */
{
   return(NULL);
}

/*****************************************************************************/
/*
 * Misc utilities
 */

/*
 * Clear a REGION
 */
int
shRegClear(REGION *reg)			/* region to clear */
{
   int i, j;
   int ncol,nrow;			/* unpacked from reg */

   shAssert(reg != NULL);

   ncol = reg->ncol;
   nrow = reg->nrow;

   switch (reg->type) {
    case TYPE_U16:
      for(i = 0;i < nrow;i++) {
	 memset(reg->rows_u16[i],ncol*sizeof(U16), '\0');
      }
      shAssert(reg->rows_u16[0][0] == 0); /* check 0 is all-bits-0 */
      break;
    case TYPE_FL32:
      for(i = 0;i < nrow;i++) {
	 memset(reg->rows_fl32[i],ncol*sizeof(FL32), '\0');
      }
      shAssert(reg->rows_fl32[0][0] == 0.0); /* check 0.0 is all-bits-0 */
      break;
    default:
      shFatal("shRegClear doesn't handle regions of type %d\n", reg->type);
      return(SH_GENERIC_ERROR);		/* NOTREACHED */
   }

   return(SH_SUCCESS);
}

#endif

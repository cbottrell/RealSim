#if !defined(PHDERVISH_H)		/* not DERVISH_H -- this is a fake */
#define PHDERVISH_H

const char *phPhotoVersion(void);
/*
 * Try to include the real dervish.h; if we succeed it'll define DERVISH_H
 * and we'll know not to provide our own fake version.
 */
#include <dervish.h>

#if !defined(DERVISH_H)			/* we didn't find the real one */
/*
 * functions usually provided by dervish
 */
typedef int RET_CODE;
typedef int TYPE;
#define UNKNOWN 1

#define SH_SUCCESS 0
#define SH_GENERIC_ERROR -1

TYPE shTypeGetFromName(const char *type);
/*
 * error reporting
 */
void shError(char *fmt, ...);
void shErrStackPush(char *fmt, ...);
void shFatal(char *fmt, ...);
/*
 * memory
 */
void *shMalloc(size_t n);
void *shRealloc(void *ptr, size_t n);
void shFree(void *ptr);
int p_shMemRefCntrGet(void *ptr);
void p_shMemRefCntrDecr(void *ptr);
/*
 * assertions
 */
#include <assert.h>
#define shAssert assert
/*
 * REGIONs
 *
 * TYPE_... must agree with dervish's region.h
 */
#define TYPE_U8 2
#define TYPE_S8 4
#define TYPE_U16 8
#define TYPE_S16 16
#define TYPE_U32 32
#define TYPE_S32 64
#define TYPE_FL32 128

typedef unsigned char U8;
typedef char S8;
typedef unsigned short U16;
typedef short S16;
typedef unsigned int U32;
typedef int S32;
typedef float FL32;

typedef struct {
   int nrow, ncol;			/* size of mask */
   unsigned char **rows;		/* data in mask */
   int row0, col0;			/* origin of mask in larger mask */
} MASK;

typedef int PIXDATATYPE;		/* type of REGION */

typedef struct {			/* must agree with dervish up to col0*/
   char *name;				/* a unique identifier */
   int nrow;				/* number of rows in image */
   int ncol;				/* number of columns in image */
   PIXDATATYPE type;			/* pixel data type */
   U16	          **rows;		/* pointer to pointers to rows */
   U8             **rows_u8;
   S8             **rows_s8;
   U16 	          **rows_u16;
   S16 	          **rows_s16;
   U32            **rows_u32;
   S32		  **rows_s32;
   FL32	          **rows_fl32;
   MASK *mask;				/* associated bad pixel mask */
   int row0,col0;			/* location of LLH corner of child */
} REGION;

REGION *shRegNew(const char *, int nrow, int ncol, int type);
void shRegDel(REGION *reg);

MASK *shMaskNew(const char *name, int nrow, int ncol);
void shMaskDel(MASK *mask);
void shMaskClear(MASK *mask);

/*
 * Chains
 */
#define AFTER 0
#define BEFORE 1
#define TAIL 1

typedef struct chain_elem {
   struct chain_elem *pNext;
   void *pElement;
} CHAIN_ELEM;

typedef struct chain {
   int nElements;
   CHAIN_ELEM *pFirst;
   TYPE type;
} CHAIN;

CHAIN *shChainNew(char *type);
void shChainDel(CHAIN *ch);
void shChainElementAddByPos(CHAIN *ch, void *el, char *type, int w, int);
void *shChainElementGetByPos(const CHAIN *ch, int el);
void *shChainElementRemByPos(const CHAIN *ch, int el);

int shRegClear(REGION *reg);

#endif
#endif

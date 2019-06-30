#if !defined(PH_RICE_H)
#define PH_RICE_H

#include "phConsts.h"

#define NOT_COMPRESSED 0		/* compression protocol used */
#define RICE_COMPRESSION 1		/* Rice compression, pre Feb 1999 */
#define RICE_COMPRESSION_V2 2		/* Rice Compression, Feb 1999 --- */

#define RICE_J 128			/* number of pixels to treat together*/
/*
 * types for RICE_BLOCK->option
 */
#define RICE_OPTION_LOW	0
#define RICE_OPTION_RUN 0x10		/* modifiers for */
#define RICE_OPTION_EXTENDED_LOW 0x20	/*               RICE_OPTION_LOW */
#define RICE_OPTION_FS	1		/* RICE_OPTION_Kn == RICE_OPTION_FS+n*/
#define RICE_OPTION_K1	2
#define RICE_OPTION_K2	3
#define RICE_OPTION_K3	4
#define RICE_OPTION_K4	5
#define RICE_OPTION_K5	6
#define RICE_OPTION_K6	7
#define RICE_OPTION_K7	8
#define RICE_OPTION_K8	9
#define RICE_OPTION_K9	10
#define RICE_OPTION_K10	11
#define RICE_OPTION_K11	12
#define RICE_OPTION_K12	13
#define RICE_OPTION_K13	14
#define RICE_OPTION_DEF	15
#define RICE_OPTION_MASK 0xf		/* mask for option without modifiers */
/*
 * a block of RICE_J pixels of compressed data
 */
typedef struct {
   unsigned char option;		/* type of block; see above */
   U16 reference;			/* reference value */
   U16 var[RICE_J - 1];			/* variable part of block */
} RICE_BLOCK;

typedef struct {
   int nblock;				/* number of blocks */
   int npix;				/* number of pixels in segment */
   short packed;			/* are the BLOCKs packed? */
   short J;				/* number of pixels per block */
   RICE_BLOCK *blocks;			/* compressed data blocks */
   int *len;				/* used length of blocks[].var,
					   in bytes */
} RICE_SEG;

RICE_SEG *phRiceSegNew(int size);
RICE_SEG * phRiceSegRealloc(RICE_SEG *rseg, int size);
void phRiceSegDel(RICE_SEG *rseg);

void
phRiceSegCompress(const U16 *data,	/* data to compress */
		  int ndata,		/* dimen of data */
		  RICE_SEG *seg,	/* where to write the compressed data*/
		  int J);		/* number of pixels to treat together*/
void
phRiceSegUncompress(U16 *data,		/* where to write the data*/
		    const RICE_SEG *seg); /* data to uncompress */


int
phRiceSegFlatten(const RICE_SEG *rseg,	/* RICE_SEG to flatten */
		 unsigned char *buff,	/* into this buffer; can be NULL */
		 int len);		/* size of buff */
RICE_SEG *
phRiceSegInflate(RICE_SEG *rseg,
		 unsigned char *buff,	/* buffer to inflate from */
		 int *plen);		/* if non-NULL:
					   on input, length of buff or 0
					   on output, number of bytes read */

#endif

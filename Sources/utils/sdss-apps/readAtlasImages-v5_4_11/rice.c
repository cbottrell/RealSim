#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "dervish.h"
#include "phRice.h"

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
swap4(unsigned char *buff,
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
swap4(unsigned char *buff,		/* NOTUSED */
      int n)				/* NOTUSED */
{}
#endif

/*****************************************************************************/

RICE_SEG *
phRiceSegNew(int size)			/* number of blocks desired, or 0 */
{
   RICE_SEG *rseg = shMalloc(sizeof(RICE_SEG));

   rseg->npix = -1;			/* unknown */
   rseg->packed = 0;
   rseg->J = RICE_J;
   if(size <= 0) {
      rseg->nblock = 0;
      rseg->blocks = NULL;
      rseg->len = NULL;
   } else {
      rseg->nblock = size;
      rseg->blocks = shMalloc(size*sizeof(RICE_BLOCK));
      rseg->len = shMalloc(size*sizeof(int));
   }

   return(rseg);
}

RICE_SEG *
phRiceSegRealloc(RICE_SEG *rseg,	/* RICE_SEG to reallocate */
		 int size)		/* number of blocks desired, or 0 */
{
   if(rseg == NULL) {
      rseg = phRiceSegNew(size);
   } else if(rseg->nblock == size) {
      ;					/* nothing to do */
   } else {
      shAssert(size > 0);

      rseg->nblock = size;
      rseg->blocks = shRealloc(rseg->blocks, size*sizeof(RICE_BLOCK));
      rseg->len = shRealloc(rseg->len, size*sizeof(int));
   }

   return(rseg);
}

void
phRiceSegDel(RICE_SEG *rseg)
{
   if(rseg != NULL) {
      shFree(rseg->len);
      shFree(rseg->blocks);
      shFree(rseg);
   }
}

/*****************************************************************************/
/*
 * compress J pixels into block
 */
#define BITS_PER_BYTE 8
#define U16LEN (BITS_PER_BYTE*sizeof(U16))

/*****************************************************************************/
/*
 * apply the "standard mapper" to a prediction and an observed value
 */
static unsigned int
get_delta(int pred,			/* predicted value */
	  int val)			/* observed value */
{
   int delta;				/* mapped difference of two pixels */
   int Delta = val - pred;
   int theta = (pred <= MAX_U16/2) ? pred : MAX_U16 - pred;

   if(Delta < -theta) {
      delta = theta - Delta;
   } else if(Delta < 0) {
      delta = 2*(-Delta) - 1;
   } else if(Delta <= theta) {
      delta = 2*Delta;
   } else {
      delta = theta + Delta;
   }
   shAssert(delta >= 0);      

   return(delta);
}

/*****************************************************************************/
/*
 * undo the "standard mapper" given a prediction and an mapped value
 */
static int
unget_delta(int pred,			/* predicted value */
	    int delta)			/* mapped value */
{
   int val;				/* mapped difference of two pixels */
   int theta = (pred <= MAX_U16/2) ? pred : MAX_U16 - pred;

   if(delta <= 2*theta) {
      if(delta & 0x1) {
	 val = -(delta + 1)/2;
      } else {
	 val = delta/2;
      }
   } else {
      if(pred > MAX_U16/2) {
	 val = theta - delta;
      } else {
	 val = delta - theta;
      }
   }

   return(pred + val);
}

/*****************************************************************************/
/*
 * find out how long a block would be in bytes, once encoded.
 * The length does _not_ include the reference value
 */
#if 0
#  define PREDICT(DATA, I) ((3*DATA[I] - DATA[I-1])/2)
#else
#  define PREDICT(DATA, I) (DATA[I])
#endif

static int
block_length(const U16 *data,		/* data to compress */
	     int ndata,			/* number of pixels to process */
	     unsigned int nbit)		/* number of bits of noise */
{
   int i;
   int len;				/* length if FS coded */
   int pred;				/* predicted value */
/*
 * Find out how many bits the non-noise bits would take up if encoded as
 * FS codes. Note that we don't encode data[0], as it's the reference
 */
   len = 0;
   pred = data[0];
   for(i = 1;i < ndata;i++) {
      len += get_delta(pred, data[i]) >> nbit;
      pred = PREDICT(data, i);
   }
   len += (ndata - 1);			/* encoding 1 costs 2 bits */
   len += (ndata - 1)*nbit;		/* add in the noise bits */

   return(len);
}

/*****************************************************************************/
/*
 * Actually process a block. Return the number of bytes written to block->var
 */
static int
do_block(const U16 *data,		/* data to compress */
	 int ndata,			/* number of pixels to process */
	 RICE_BLOCK *block,		/* block to set */
	 int blen,			/* length of compressed data stream
					   (in bits) */
	 unsigned int nbit,		/* number of bits of noise */
	 unsigned int data_mask)	/* mask to remove nbit */
{
   unsigned int delta;			/* mapped difference of two pixels */
   int i;
   int len;				/* length of compressed data in bytes*/
   int nleft;				/* how many more bits fit in word */
   U16 noise[RICE_J];			/* buffer for noise bits */
#define DEBUG 0
#if DEBUG
   U16 darray[RICE_J];			/* array of values of delta */
#endif
   const unsigned int noise_mask = ~data_mask; /* mask to isolate nbit */
   int pred;				/* predicted value */
   int skip_words = 0;			/* how many 0 U16 words to skip over */
   U16 *vptr;				/* FS coded (trimmed) values of data */

   len = blen/BITS_PER_BYTE;
   if(blen != len*BITS_PER_BYTE) len++;
/*
 * We have to ensure that we have an integral number of U16s, or we won't
 * be able to byteswap. This wastes about 0.5bytes; a shame but not a disaster
 *
 * Furthermore, on little-endian machines we need to clear all the bytes in
 * our U16 buffers; clearing the number used is not good enough as we'll be
 * using the "wrong" end of the last short-word if len is odd
 */
   if(len%2 == 1) {
      len++;
   }

   if(blen > (ndata - 1)*U16LEN) {
      nbit = 16;
   }

   if(nbit <= 13) {
      block->option = RICE_OPTION_FS + nbit;
   } else if(nbit == 16) {
      block->option = RICE_OPTION_DEF;
   } else {
      shFatal("do_block: impossible value of nbit %d", nbit);
   }
   block->reference = data[0];

   if(block->option == RICE_OPTION_DEF) {
      len = (ndata - 1)*sizeof(U16);	/* ndata-1 => not the reference value*/
      memcpy(block->var, &data[1], len);
      swap2((void *)&block->reference, sizeof(U16));
      swap2((void *)block->var, len);
      return(len);
   }
/*
 * OK, now do the clipping and encoding of the non-noise bits.
 *
 * Note that FS(0) == 0x1, requiring us to write 1 bit,
 * so (delta + 1) is the number of bits required
 */
   memset(block->var, '\0', len + len%2); /* an even number  */

   pred = data[0];
   vptr = block->var;
   nleft = U16LEN;
   for(i = 0;i < ndata - 1;i++) {
      delta = get_delta(pred, data[i+1]);
#if DEBUG
      if(unget_delta(pred, delta) != data[i + 1]) {
	 fprintf(stderr,"Problem with delta: %d %d %d\n",
		 pred, delta, data[i + 1]);
      }
#endif
   
      pred = PREDICT(data, i + 1);
#if DEBUG
      darray[i] = delta;
#endif

      noise[i] = delta & noise_mask;
      delta = (delta & data_mask) >> nbit;

      if(delta + 1 > nleft) {		/* it won't all fit in this word */
	 delta -= nleft;
	 vptr++; nleft = U16LEN;

	 skip_words = delta/U16LEN;
	 vptr += skip_words;

	 delta -= skip_words*U16LEN;
      }

      shAssert(delta + 1 <= nleft);
      *vptr |= (1 << (nleft - (delta + 1)));
      nleft -= delta + 1;
      shAssert(nleft >= 0);
      if(nleft == 0) {
	 vptr++; nleft = U16LEN;
      }
   }
   shAssert(vptr - block->var < sizeof(block->var));
/*
 * and now the noise bits, if any
 */
   if(nbit > 0) {
      U16 n;				/* == noise[] */
      int nb;				/* number of bits still to write */
      for(i = 0;i < ndata - 1;i++) {
	 n = noise[i];
	 nb = nbit;
	 
	 if(nb > nleft) {		/* won't all fit in this word */
	    *vptr++ |= (n >> (nb - nleft));
	    nb -= nleft;
	    nleft = U16LEN;
	 }
	 *vptr |= n << (nleft - nb);
	 nleft -= nb;
      }
   }
/*
 * byte swap if required. The compressed data is just a string of bits,
 * but it's treated as a series of U16 numbers for efficiency of manipulation.
 */
   swap2((void *)&block->reference, sizeof(U16));
   swap2((void *)block->var, len);
/*
 * Done. Check that we wrote the right number of bytes if shAssert is active
 */
#if !defined(NDEBUG)
   i = (char *)vptr - (char *)block->var; /* number of complete words */
   
   if(nleft < U16LEN/2) {
      i += 2;
   } else if(nleft < U16LEN) {
      i++;
   }
   if(i%2 == 1) i++;			/* as was done for len a moment ago */
   shAssert(i == len);
   shAssert(len <= sizeof(block->var));
#endif

   return(len);
}

/*
 * Actually process a block using one of the low-entropy algorithms,
 * as specified in block->option on input
 *
 * Return the number of bytes written to block->var
 */
static int
do_low_entropy_block(const U16 *data,	/* data to compress */
		     int ndata,		/* no. of pixels to process NOTUSED */
		     RICE_BLOCK *block)	/* block to set */
{
   switch(block->option) {
    case (RICE_OPTION_LOW | RICE_OPTION_RUN):
      block->reference = data[0];
      swap2((void *)&block->reference, sizeof(U16));
      return(0);
    case (RICE_OPTION_LOW | RICE_OPTION_EXTENDED_LOW):
      shFatal("do_low_entropy_block doesn't yet support "
	      "RICE_OPTION_EXTENDED_LOW");
      break;				/* NOTREACHED */
    default:
      shError("Unknown option to do_low_entropy_block: %d", block->option);
      break;				/* NOTREACHED */
   }
   
   return(0);				/* NOTREACHED */
}

/*****************************************************************************/
/*
* c & data_mask[n] => remove last n bits
 */
static const U16 data_mask[U16LEN] = {
   0xffff, 0xfffe, 0xfffc, 0xfff8,
   0xfff0, 0xffe0, 0xffc0, 0xff80,
   0xff00, 0xfe00, 0xfc00, 0xf800,
   0xf000, 0xe000, 0xc000, 0x8000,
};
/*
 * Compress a segment of data (e.g. a row of an image)
 */
void
phRiceSegCompress(const U16 *data,	/* data to compress */
		  int ndata,		/* dimen of data */
		  RICE_SEG *seg,	/* where to write the compressed data*/
		  int J)		/* number of pixels to treat together*/
{
   RICE_BLOCK *blocks;			/* == seg->blocks */
   int i, j, k;
   int len;				/* length of compressed block */
   int len_best;			/* best value of len */
   int n;				/* number of pixels to compress */
   int nbit;				/* number of bits of noise to trim */
   int nbit_best;			/* best value of nbit */
   
   if(J == 0) {
      J = RICE_J;
   } else if(J > RICE_J) {
      fprintf(stderr,"Maximum number of pixels in a block is %d\n",RICE_J);
      J = RICE_J;
   }
   seg->J = J;

   seg->npix = ndata;
   blocks = seg->blocks;
/*
 * Process the data
 */
   for(i = j = 0; i < ndata; i += J, j++) {
      n = (ndata - i > J ? J : ndata - i);

      nbit = 0;
      len = block_length(&data[i], n, nbit);
      len_best = len; nbit_best = nbit; 

      for(nbit = 1; nbit <= 13; nbit++) {
	 len = block_length(&data[i], n, nbit);
	 if(len < len_best) {	
	    len_best = len; nbit_best = nbit; 
	 }
      }
/*
 * if we want to trim 0 bits of noise, this is a candidate for a
 * low entropy option. The Rice standard defines two such options,
 * zero-run and "extended low entropy"; we add a third, a soft-bias run.
 */
      if(nbit_best == 0) {			/* may be a run */
	 for(k = 1; k < n; k++) {
	    if(data[i + k] != data[i + k-1]) {
	       break;
	    }
	 }
	 if(k == n) {		/* yes; all the values are identical */
	    blocks[j].option = (RICE_OPTION_LOW | RICE_OPTION_RUN);
	    seg->len[j] = do_low_entropy_block(&data[i], n, &blocks[j]);
	       
	    continue;
	 }
      }
      
      seg->len[j] = do_block(&data[i], n, &blocks[j],
			     len_best, nbit_best, data_mask[nbit_best]);
   }
}

/*****************************************************************************/
/*
 * uncompress a block; return the total number of bytes read from the variable
 * part; this excludes the option and reference value
 */
static int
undo_block(const RICE_BLOCK *block,	/* block to uncompress */
	   int ndata,			/* number of pixels to process */
	   U16 *data)			/* put data here */
{
   int delta;				/* mapped difference of two pixels */
   int i, j;
   int len;				/* length of a section of data, bytes*/
   int nbit;				/* number of bits of noise */
   int nleft;				/* how many more bits fit in word */
   int pred;				/* predicted value */
   U16 vval;				/* == block->var[] */

   data[0] = block->reference;

   switch(block->option) {
    case RICE_OPTION_DEF:
      len = ndata*sizeof(U16);
      memcpy(data, block->var, len);	/* including reference, data[0] */
      swap2((void *)data, len);
      
      return(len - sizeof(U16));	/* excludes option and reference */
    case (RICE_OPTION_LOW | RICE_OPTION_RUN):
      swap2((void *)data, sizeof(U16));
      vval = data[0];
      for(i = 1; i < ndata; i++) {
	 data[i] = vval;
      }
      
      return(0);			/* excludes option and reference */
    case (RICE_OPTION_LOW | RICE_OPTION_EXTENDED_LOW):
      shFatal("undo_block doesn't yet support "
	      "RICE_OPTION_EXTENDED_LOW");
      break;				/* NOTREACHED */
    default:
      break;
   }
/*
 * We have to unpack the compressed data
 */
   shAssert(block->option <= RICE_OPTION_DEF);
   nbit = block->option - RICE_OPTION_FS;

   j = 0;
   vval = block->var[j++];
   swap2((void *)&vval, sizeof(U16));

   nleft = U16LEN;
   for(i = 1; i < ndata; i++) {
      delta = 0;

      while(vval == 0) {
	 delta += nleft;
	 vval = block->var[j++];
	 swap2((void *)&vval, sizeof(U16));
	 nleft = U16LEN;
      }

      while((vval & (1 << (U16LEN - 1))) == 0) {
	 delta++;
	 vval <<= 1;
	 nleft--;
      }
      vval <<= 1; nleft--;

      data[i] = delta;
   }
/*
 * byte swap if required.
 */
   swap2((void *)&data[0], sizeof(U16));
/*
 * and the noise bits back
 */
   if(nbit > 0) {
      U16 n;				/* == noise */
      int nb;				/* number of bits still to read */
      const unsigned int noise_mask = data_mask[U16LEN - nbit]; /* get nbit
								   left bits */

      for(i = 1;i < ndata;i++) {
	 n = 0;
	 nb = nbit;
	 
	 if(nb > nleft) {		/* it isn't all in this word */
	    n = vval & noise_mask;	/* here's part that is */
	    nb -= nleft;

	    vval = block->var[j++];
	    swap2((void *)&vval, sizeof(U16));
	    nleft = U16LEN;
	    
	    n |= (vval & noise_mask) >> (nbit - nb);
	 } else {
	    n = vval & noise_mask;
	 }
	 vval <<= nb;
	 nleft -= nb;

	 data[i] = (data[i] << nbit) | (n >> (U16LEN - nbit));
      }
   }
/*
 * and undo the preprocessor
 */
   pred = data[0];
   for(i = 1; i < ndata; i++) {
      data[i] = unget_delta(pred, data[i]);
      pred = PREDICT(data, i);
   }
/*
 * find how many bytes of padded[] we used; note that we always write an
 * multiple of sizeof(U16) bytes, so we may waste a single byte at the end
 */
   j = j*sizeof(U16);			/* number of words we processed */

   return(j);
}

/*****************************************************************************/
/*
 * uncompress a packed RICE_BLOCK
 */
static int
undo_pblock(const char *packed,		/* packed block to uncompress */
	   int ndata,			/* number of pixels to process */
	   U16 *data)			/* put data here */
{
   RICE_BLOCK block;			/* block to uncompress */
   int len;				/* length of variable part */

   block.option = packed[0]; packed++;
   memcpy(&block.reference, packed, sizeof(U16)); packed += sizeof(U16);
   memcpy(block.var, packed, sizeof(block.var));
   
   len = 1 + sizeof(U16) + undo_block(&block, ndata, data);

   return(len);
}

/*****************************************************************************/
/*
 * Uncompress a Rice-compressed segment
 */
void
phRiceSegUncompress(U16 *data,		/* where to write the data*/
		    const RICE_SEG *seg) /* data to uncompress */
{
   const char *bptr = (char *)seg->blocks; /* pointer to blocks in segment */
   int i;
   int J = seg->J;			/* number of pixels per segment */
   int n;				/* number of pixels to uncompress
					   from a block */
   int ndata = seg->npix;
   int packed = seg->packed;		/* is segment packed? */

   for(i = 0; i < seg->nblock; i++) {
      n = (ndata - i*J > J ? J : ndata - i*J);
      if(packed) {
	 bptr += undo_pblock(bptr, n, &data[i*J]);
      } else {
	 (void)undo_block((const RICE_BLOCK *)bptr, n, &data[i*J]);
	 bptr += sizeof(RICE_BLOCK);
      }
   }
}

/*****************************************************************************/
/*
 * flatten/inflate a RICE_SEG to/from a network-byte-order array of bytes
 */
int
phRiceSegFlatten(const RICE_SEG *rseg,	/* RICE_SEG to flatten */
		 unsigned char *buff,	/* into this buffer; can be NULL */
		 int len)		/* size of buff */
{
   unsigned char *const buff0 = buff;	/* initial value of buff */
   int j;
   int nbyte;				/* no. of bytes in compressed blocks */
   int size;
   
   nbyte = 0;
   for(j = 0; j < rseg->nblock; j++) {
      nbyte++;				/* option byte */
      nbyte += sizeof(U16);		/* reference value */
      nbyte += rseg->len[j];		/* the compressed data */
   }
   size = 3*sizeof(int) + nbyte;   

   if(buff == NULL) {
      return(size);
   }
   shAssert(len >= size);
/*
 * first the ints, describing e.g. the number of blocks
 */
   memcpy(buff, &rseg->nblock, sizeof(int)); buff += sizeof(int);
   memcpy(buff, &rseg->npix, sizeof(int)); buff += sizeof(int);
   memcpy(buff, &nbyte, sizeof(int)); buff += sizeof(int);
   swap4(buff - 3*sizeof(int), 3*sizeof(int));
/*
 * now the compressed blocks themselves
 */
   for(j = 0; j < rseg->nblock; j++) {
      memcpy(buff,&rseg->blocks[j].option, 1); buff++;
      memcpy(buff,&rseg->blocks[j].reference,sizeof(U16)); buff += sizeof(U16);
      memcpy(buff,rseg->blocks[j].var, rseg->len[j]); buff += rseg->len[j];
   }

   shAssert(buff <= buff0 + len);

   return(size);
}

RICE_SEG *
phRiceSegInflate(RICE_SEG *rseg,
		 unsigned char *buff,	/* buffer to inflate from */
		 int *plen)		/* if non-NULL:
					   on input, length of buff or 0
					   on output, number of bytes read */
{
   unsigned char *const buff0 = buff;	/* initial value of buff */
   int len = 0;				/* length of a block */
   int nblock = 0;			/*  */

   shAssert(buff != NULL);
/*
 * first the ints, 
 */
   swap4(buff,3*sizeof(int));
   memcpy((void *)&nblock, buff, sizeof(int)); buff += sizeof(int);

   rseg = phRiceSegRealloc(rseg, nblock);
   memcpy((void *)&rseg->npix, buff, sizeof(int)); buff += sizeof(int);
   memcpy((void *)&len, buff, sizeof(int)); buff += sizeof(int);
   
   memcpy((void *)&rseg->blocks[0], buff, len); buff += len;
   rseg->packed = 1;

   if(plen != NULL) {
      if(*plen > 0) {
	 shAssert(buff - buff0 <= *plen);
      }
      *plen = buff - buff0;
   }

   return(rseg);
}

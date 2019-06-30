/*
 * Write a FITS file; BITPIX == 8 (MASK) or 16 (U16 REGION) only
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "phFits.h"

#define BITS_PER_BYTE 8			/* number of bits in a byte */
#define FITSIZE 2880			/* size of FITS record
					   Note that we need FITSIZE+1
					   in write_card() */
static char record[FITSIZE + 1];

static void flip_high_bits(char *record, int nbyte);
#if defined(SDSS_LITTLE_ENDIAN)
static void swab2(char *buff, int n);
static void swab4(char *buff, int n);
#endif
static void write_card_d(char *card, char *keyword, int value, char *comment);
static void write_card_f(char *card, char *keyword, float value,char *comment);
static void write_card_l(char *card, char *keyword, char *value,char *comment);
static void write_card_s(char *card, char *keyword, char *value,char *comment);

void
write_fits_file(const char *name,	/* name of the file */
		const void *what,	/* what we are supposed to write */
		int bitpix)		/* BITPIX */
{
   const float bscale = 1.0, bzero = 32768; /* from the FITS header */
   unsigned char *data;			/* U16 data to write */
   int fd;
   const MASK *mask;			/* == data if bitpix == 8 */
   int nbytes;				/* number of bytes of data to write */
   int ncard;
   int nrow, ncol;			/* size of data region */
   const REGION *reg;			/* == data if bitpix == 16/-32 */

   shAssert(what != NULL);

   switch (bitpix) {
    case 8:
      mask = what;
      nrow = mask->nrow; ncol = mask->ncol;
      data = mask->rows[0];

      shAssert(mask->rows[1] == mask->rows[0] + ncol); /* data's contiguous */
      break;
    case 16:
      reg = what;
      nrow = reg->nrow; ncol = reg->ncol;
      data = (unsigned char *)reg->rows[0];

      shAssert(reg->type == TYPE_U16);
      shAssert(reg->rows[1] == reg->rows[0] + ncol); /* data's contiguous */
      break;
    case -32:
      reg = what;
      nrow = reg->nrow; ncol = reg->ncol;
      data = (unsigned char *)reg->rows_fl32[0];

      shAssert(reg->type == TYPE_FL32);
      shAssert(reg->rows_fl32[1] == \
	       reg->rows_fl32[0] + ncol); /* data's contiguous */
      break;
    default:
      data = NULL; nrow = ncol = 0;	/* supress compiler warnings */
      shFatal("Illegal value of bitpix: %d", bitpix);
   }
   
   if((fd = creat(name,0644)) < 0) {
      fprintf(stderr,"Can't create %s\n",name);
      return;
   }

   ncard = 0;
   write_card_l(&record[80*ncard++],"SIMPLE","T","");
   write_card_d(&record[80*ncard++],"BITPIX",bitpix,"");
   write_card_d(&record[80*ncard++],"NAXIS",2,"");
   write_card_d(&record[80*ncard++],"NAXIS1",ncol,"");
   write_card_d(&record[80*ncard++],"NAXIS2",nrow,"");
   if(bitpix == 16) {
      write_card_f(&record[80*ncard++],"BZERO",bzero,"");
      write_card_f(&record[80*ncard++],"BSCALE",bscale,"");
   }
   write_card_s(&record[80*ncard++],"END","","");

   while(ncard%36 != 0) {
      write_card_s(&record[80*ncard++],"","","");
   }
   write(fd,record,FITSIZE);
/*
 * Now write data
 */
   nbytes = ncol*nrow*(bitpix > 0 ? bitpix : -bitpix)/BITS_PER_BYTE;

   while(nbytes >= FITSIZE) {
      memcpy(record,data,FITSIZE);
      if(bitpix == 16) {
	 flip_high_bits(record, FITSIZE);
#if defined(SDSS_LITTLE_ENDIAN)
	 swab2(record, FITSIZE);
#endif
      } else if(bitpix == -32) {
#if defined(SDSS_LITTLE_ENDIAN)
	 swab4(record, FITSIZE);
#endif
      }
      
      if(write(fd,record,FITSIZE) != FITSIZE) {
	 perror("Writing data");
	 break;
      }
      data += FITSIZE;
      nbytes -= FITSIZE;
   }
   if(nbytes > 0) {
      memcpy(record,data,nbytes);
      if(bitpix == 16) {
	 flip_high_bits(record, nbytes);
#if defined(SDSS_LITTLE_ENDIAN)
	 swab2(record, nbytes);
#endif
      } else if(bitpix == -32) {
#if defined(SDSS_LITTLE_ENDIAN)
	 swab4(record, FITSIZE);
#endif
      }
      memset(record + nbytes,' ',FITSIZE - nbytes);
      if(write(fd,record,FITSIZE) != FITSIZE) {
	 perror("Writing final record");
      }
   }

   close(fd);
}

#if defined(SDSS_LITTLE_ENDIAN)
static void
swab2(char *buff,			/* buffer to swap ABAB --> BABA */
      int n)				/* number of _bytes_ */
{
   int i;
   char tmp;

   for(i = 0; i < n; i += 2) {
      tmp = buff[i]; buff[i] = buff[i + 1]; buff[i + 1] = tmp;
   }
}

static void
swab4(char *buff,			/* buffer to swap ABCD --> DCBA */
      int n)				/* number of _bytes_ */
{
   int i;
   char tmp;

   for(i = 0; i < n; i += 4) {
      tmp = buff[i];     buff[i] = buff[i + 3];     buff[i + 3] = tmp;
      tmp = buff[i + 1]; buff[i + 1] = buff[i + 2]; buff[i + 2] = tmp;
   }
}
#endif

/*****************************************************************************/
/*
 * Go through a record flipping the high order bit on every short. This
 * is required by the idiocy of the FITS committee in refusing to allow
 * unsigned short as a valid fits data type
 */
static void
flip_high_bits(char *record,		/* data to flip */
	       int nbyte)		/* length of record in bytes */
{
   int i;
   int n = nbyte/2;			/* number of shorts */
   unsigned short *ptr = (unsigned short *)record;

   for(i = 0; i < n; i++) {
      ptr[i] ^= 0x8000;
   }
}


/*****************************************************************************/
/*
 * Write a FITS card image. Various calls for character, int, logical, float
 *
 * Write a string value
 */
static void
write_card_s( char *card, char *keyword, char *value, char *comment )
{
   if(*keyword == '\0' || !strcmp(keyword,"COMMENT") || !strcmp(keyword,"END")
		       || !strcmp(keyword,"HISTORY")) {
      if(comment[0] != '\0') {
	 fprintf(stderr,
       	       "You can't add a comment to a COMMENT, END, or HISTORY card\n");
      }
      sprintf(card,"%-8.8s%-72s",keyword,value);
   } else {
      sprintf(card,"%-8.8s= '%-8s' %c%-*s", keyword, value,
	      (comment[0] == '\0' ? ' ' : '/'),
	      (int)(80 - 14 - strlen(value)), comment);
   }
}   

/*****************************************************************************/
/*
 * Write an integer (%d)
 */
static void
write_card_d( char *card, char *keyword, int value, char *comment )
{
   sprintf(card,"%-8.8s= %20d %c%-48s",keyword,value,
	   		(comment[0] == '\0' ? ' ' : '/'),comment);
}   

/*****************************************************************************/
/*
 * Write a logical value
 */
static void
write_card_l( char *card, char *keyword, char *value, char *comment )
{
   if(strcmp(value,"T") != 0 && strcmp(value,"F") != 0) {
      fprintf(stderr,"Invalid logical %s for keyword %s\n",value,keyword);
      sprintf(value,"?");
   }
   sprintf(card,"%-8.8s= %20s %c%-48s",keyword,value,
	   		(comment[0] == '\0' ? ' ' : '/'),comment);
}   

/*****************************************************************************/
/*
 * Write a floating value
 */
static void
write_card_f( char *card, char *keyword, float value, char *comment )
{
   sprintf(card,"%-8.8s= %20g %c%-48s",keyword,value,
	   		(comment[0] == '\0' ? ' ' : '/'),comment);
}   

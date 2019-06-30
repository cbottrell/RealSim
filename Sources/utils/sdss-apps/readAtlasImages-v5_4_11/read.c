/*
 * Read a fits binary table
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "dervish.h"
#include "phFits.h"
#include "phDataIo.h"
#include "phVariablePsf.h"

#define FITS_SIZE 2880			/* Fits record length */

#define TYPE_OBJMASK 0			/* type of object to read */
#define TYPE_ATLAS_IMAGE 1
#define TYPE_PSF_KL_COMP 2

#if defined(DARWIN)			/* the OS/X loader has trouble wanting
					   to pull these out of the libraries,
					   but failing to do so */
void *g_saoCmdHandle;
int f_BYTESWAP, f_WORDSWAP;
#endif

extern int verbose;

static char *key;
static char *value;
static char *comment;
static char ccard[81];

/*****************************************************************************/

FITS *
phFitsNew(const char *file)		/* file to open */
{
   int fd, hfd;
   FITS *fits;

   if((fd = open(file,0)) < 0) {
      shError("Can't open %s",file);
      return(NULL);
   }
   if((hfd = open(file,0)) < 0) {
      close(fd);
      shError("Can't open %s",file);
      return(NULL);
   }

   fits = shMalloc(sizeof(FITS));
   fits->fd = fd; fits->hfd = hfd;

   fits->maskrows = fits->maskcols = 0;
   
   return(fits);
}

void
phFitsDel(FITS *fits)
{
   if(fits != NULL) {
      close(fits->fd);
      close(fits->hfd);
      shFree(fits);
   }
}

/*****************************************************************************/

static void
parse_card(const char *card)
{
   char *ptr;
   
   strncpy(ccard,card,80); ccard[80] = '\0';

   key = ccard;
   for(ptr = ccard;*ptr != '\0' && !isspace(*ptr);ptr++) {
      if(*ptr == '=') {
	 break;
      }
   }
   *ptr++ = '\0';
   
   while(isspace(*ptr)) ptr++;
   if(*ptr == '=') {
      ptr++;
      while(isspace(*ptr)) ptr++;
   }

   value = ptr;
   if(*ptr == '\'') {			/* string valued */
      ptr++; value++;
      for(;*ptr != '\0';ptr++) {
	 if(*ptr == '\'') {
	    if(*(ptr + 1) == '\'') {	/* escaped ' */
	       ptr++;
	    } else {
	       break;
	    }
	 }
      }
   } else {
      while(*ptr != '\0' && !isspace(*ptr)) ptr++;
   }
   *ptr++ = '\0';

   while(isspace(*ptr)) ptr++;
   if(*ptr == '/') {
      ptr++;
      while(isspace(*ptr)) ptr++;
   }
   comment = ptr;
}

/*****************************************************************************/

static int
read_header(FITS *fits,			/* the FITS in question */
	    const char *hdrtype,	/* type of header */
	    const char *hdrval)		/* desired value of hdrtype */
{
   char cbuff[82];			/* full filename/card buffer */
   int extended;			/* is the EXTEND keyword in header? */
   int fd = fits->fd;			/* the file's file descriptor */
   char hdr_value[81];			/* value associated with hdrtype */
   int ncard;				/* number of card images read */
   int naxis3;				/* value of NAXIS3 */
   
   *hdr_value = '\0';
   fits->naxis1 = fits->naxis2 = 0;
   fits->naxis = fits->pcount = 0;
   naxis3 = 1;
   extended = 0;
   for(ncard = 0;read(fd,cbuff,80) == 80;ncard++) {
      parse_card(cbuff);
      if(verbose > 1) printf("%s | %s | %s\n",key,value,comment);
      
      if(strcmp(key,hdrtype) == 0) {
	 strcpy(hdr_value, value);
      } else if(strcmp(key,"BITPIX") == 0) {
	 fits->bitpix = atoi(value);
      } else if(strcmp(key,"NAXIS") == 0) {
	 fits->naxis = atoi(value);
      } else if(strcmp(key,"NAXIS1") == 0) {
	 fits->naxis1 = atoi(value);
      } else if(strcmp(key,"NAXIS2") == 0) {
	 fits->naxis2 = atoi(value);
      } else if(strcmp(key,"NAXIS3") == 0) {
	 naxis3 = atoi(value);
      } else if(strcmp(key,"EXTEND") == 0) {
	 extended = (strcmp(value,"T") == 0) ? 1 : 0;
      } else if(strcmp(key,"PCOUNT") == 0) {
	 fits->pcount = atoi(value);
      } else if(strcmp(key,"THEAP") == 0) {
	 fits->theap = atoi(value);
      } else if(strcmp(key,"MASKCOLS") == 0) {
	 fits->maskcols = atoi(value);
      } else if(strcmp(key,"MASKROWS") == 0) {
	 fits->maskrows = atoi(value);
      } else if(strcmp(key,"END") == 0) {
	 break;
      } else {
	 ;				/* should save the entry somewhere. */
      }
   }
   if(ncard == 0) {			/* failed to even read first card */
      shError("Failed to read card %d",ncard + 1);
      return(-1);
   }
   
   while(++ncard % 36 != 0) {		/* read rest of header */
      if(read(fd,cbuff,80) != 80) {
	 shError("Failed to read card %d",ncard + 1);
	 return(-1);
      }
   }
/*
 * if it's a primary header, SIMPLE must be present, and true
 */
   if(*hdr_value == '\0') {
      shError("failed to find keyword %s", hdrtype);
      return(-1);
   }

   if(strcmp(hdr_value, hdrval) != 0) {
      shError("expected to find value %s for keyword %s; saw %s",
						   hdrval, hdrtype, hdr_value);
      return(-1);
   }
/*
 * if there's no data, extended had better be true
 */
   if(naxis3 != 1) {
      shError("We cannot handle 3-d tables");
      return(-1);
   }
   if(strcmp(hdrtype,"SIMPLE") == 0 &&
      ((fits->naxis == 0) ||
       (fits->naxis == 1 && fits->naxis1 == 0) ||
       (fits->naxis == 2 && fits->naxis1*fits->naxis2 == 0))) {
      if(!extended) {
	 shError("Header has no associated data, and is not extended");
	 return(-1);
      }
   }

   return(0);
}

/*****************************************************************************/
/*
 * open and return a fits table. The file pointer will point to the
 * start of the data in the table
 *
 * This function is pretty general, and should be easily reused for
 * other purposes than reading atlas images
 */
FITS *
open_fits_table(const char *name,	/* name of the file */
		int hdu)		/* desired HDU */
{
   FITS *fits;				/* open a fits file */
   int i;
   int len;				/* total length of main table data */
/*
 * Start by reading the primary header
 */
   if((fits = phFitsNew(name)) == NULL) {
      return(NULL);
   }

   if(read_header(fits, "SIMPLE", "T") < 0) {
      phFitsDel(fits);
      return(NULL);
   }
   
   if(fits->naxis != 0 || fits->bitpix != 8) {
      shError("open_fits_file: %s's PDU must have no data (of type byte)",
									 name);
      phFitsDel(fits);
      return(NULL);
   }
/*
 * Now the secondary header
 */
   for(i = 1; i <= hdu; i++) {
      if(read_header(fits, "XTENSION", "BINTABLE") < 0) {
	 phFitsDel(fits);
	 return(NULL);
      }
      if(i == hdu) {			/* we're there */
	 break;
      }
/*
 * skip over main table data, allowing for padding at the end of the table
 */
      len = (fits->naxis1*fits->naxis2 + fits->pcount);
      if(len%FITS_SIZE != 0) {
	 len = (len/FITS_SIZE + 1)*FITS_SIZE;
      }
      
      if(lseek(fits->fd, len, SEEK_CUR) == -1) {
	 shError("open_fits_table: cannot skip %dth HDU", i);
	 phFitsDel(fits);
	 return(NULL);
      }
   }
/*
 * remember the start of the data, and move heap pointer to start of heap
 */
   if((fits->data_start = lseek(fits->fd,0,SEEK_CUR)) == -1) {
      shError("open_fits_table: cannot find where file pointer is");
      phFitsDel(fits);
      return(NULL);
   }
   
   fits->heap_start = fits->data_start + fits->theap;
   if(lseek(fits->hfd,fits->heap_start,SEEK_SET) == -1) {
      shError("open_fits_table: cannot seek to start of heap");
      phFitsDel(fits);
      return(NULL);
   }
   
   return(fits);
}

/*****************************************************************************/
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
swab4(unsigned char *buff,
      int n)
{
   int i;
   char tmp;

   for(i = 0; i < n; i += 4) {
      tmp = buff[i];     buff[i] =     buff[i + 3]; buff[i + 3] = tmp;
      tmp = buff[i + 1]; buff[i + 1] = buff[i + 2]; buff[i + 2] = tmp;
   }
}
#endif

/*****************************************************************************/
/*
 * read a row from the table. We'll specialize to the cases of ATLAS_IMAGEs,
 * OBJMASKs, and KL_PSF_COMPs
 */
#define NINT 36				/* we need to read NINT ints */

static void *
read_fits_row(const FITS *fits,		/* the table */
	      int row,			/* the desired row */
	      int type)			/* what are we reading? */
{
   unsigned char *buff;			/* buffer to read heap data */
   int line[NINT];			/* buffer to read a line */
   int nbyte = 0, offset = 0;		/* size and offset of AI in heap
					   (initially indices into line[]) */
   int nint = 0;			/* number of ints to read */
   void *val = NULL;			/* the object to return */

   switch (type) {
    case TYPE_ATLAS_IMAGE:
#if 0
      nint = 3;				/* v5_2 */
      nbyte = 1; offset = 2;
#elif 0
      nint = 8;				/* v5_3 */
      nbyte = 6; offset = 7;
#else
      nint = 0;				/* compatibility */
      nbyte = -2; offset = -1;
#endif
      break;
    case TYPE_OBJMASK:
      nint = 11;
      nbyte = 9; offset = 10;
      break;
    case TYPE_PSF_KL_COMP:
      nint = 36;
      nbyte = 31; offset = 32;
      break;
    default:
      shFatal("Unknown type to read: %d", type);
   }

   shAssert(fits != NULL);
   shAssert(nint <= NINT);
   shAssert(sizeof(float) == sizeof(int));

   if(nint == 0) {
      nint = fits->naxis1/sizeof(int);
   }
   
   if(nbyte < 0) {
      nbyte += nint;
   }
   if(offset < 0) {
      offset += nint;
   }
     
   if(fits->naxis1 != nint*sizeof(int)) {
      shError("read_fits_row: expected naxis1 == %d, saw %d",
	      nint*sizeof(int), fits->naxis1);
   }

   if(row <= 0 || row > fits->naxis2) {
      shError("read_fits_row: invalid row %d (must be in range 1..%d",
	      row, fits->naxis2);
      return(NULL);
   }
   
   if(lseek(fits->fd,fits->data_start + (row-1)*fits->naxis1 ,SEEK_SET) == -1){
      shError("read_fits_row: cannot seek to start of row %d", row);
      return(NULL);
   }

   if(read(fits->fd, line, nint*sizeof(int)) != nint*sizeof(int)) {
      shError("read_fits_row: failed to read row %d", row);
      return(NULL);
   }
#if defined(SDSS_LITTLE_ENDIAN)
   swab4((void *)line, nint*sizeof(int));
#endif

   nbyte = line[nbyte];
   offset = line[offset];

   switch (type) {
    case TYPE_ATLAS_IMAGE:
      buff = shMalloc(nbyte);
      break;
    case TYPE_OBJMASK:
      {
	 OBJMASK *om = phObjmaskNew(line[1]); /* line[1] == nspan */
	 om->nspan = om->size;
	 om->row0 = line[2]; om->col0 = line[3];
	 om->rmin = line[4]; om->rmax = line[5];
	 om->cmin = line[6]; om->cmax = line[7];
	 om->npix = line[8];
	 
	 val = om;
	 buff = (unsigned char *)om->s;
      }
      break;
    case TYPE_PSF_KL_COMP:
      {
	 int i = 0;
	 PSF_KL_COMP *klc = phPsfKLCompNew(NULL);

	 klc->nrow_b = line[i]; i++;
	 klc->ncol_b = line[i]; i++;
	 memcpy(klc->c, &line[i], 25*sizeof(float)); i += 25;
	 klc->lambda = *(float *)&line[i]; i++;
	 {
	    int nrow = line[i++];
	    int ncol = line[i++];
	    PIXDATATYPE pixType = (PIXDATATYPE)line[i++];
	    klc->reg = shRegNew("", nrow, ncol, pixType);
	    klc->reg = shRegNew("", nrow, ncol, pixType);

	    i += 2;			/* offset and nbyte */
	    klc->reg->row0 = line[i++];
	    klc->reg->col0 = line[i++];
	 }
	 i++;				/* counts */
	 
	 shAssert(i == nint);

	 val = klc;
	 buff = (unsigned char *)klc->reg->rows_fl32[0];
	 nbyte *= sizeof(float);
      }
      break;
    default:
      shError("You cannot get here: %d", type);
      buff = NULL; offset = nint = 0;	/* NOTREACHED */
      break;				/* NOTREACHED */
   }

   if(lseek(fits->hfd, fits->heap_start + offset, SEEK_SET) == -1) {
      shError("read_fits_row: cannot seek to heap data for row %d", row);
      if(type) {
	 free(buff);
      } else {
	 phObjmaskDel(val);
      }
      return(NULL);
   }
   
   if(read(fits->hfd, buff, nbyte) != nbyte) {
      shError("read_fits_row: failed to read heap for row %d", row);
      if(type) {
	 free(buff);
      } else {
	 phObjmaskDel(val);
      }
      return(NULL);
   }

   switch (type) {
    case TYPE_ATLAS_IMAGE:
      val = phAtlasImageInflate(NULL, buff, &nbyte);
      shFree(buff);
      {
	 ATLAS_IMAGE *ai = val;
      	 ai->id = line[0];
      	 ai->run = line[1];
	 ai->rerun = line[2];
	 ai->camCol = line[3];
	 ai->field = line[4];
	 ai->parent = line[5];
      }

      break;
    case TYPE_OBJMASK:
#if defined(SDSS_LITTLE_ENDIAN)
      swab2(buff, nbyte);
#endif
      break;
    case TYPE_PSF_KL_COMP:
#if defined(SDSS_LITTLE_ENDIAN)
      swab4(buff, nbyte);
#endif
      break;
   }

   return(val);
}

/*****************************************************************************/
/*
 * wrappers to read atlas images or objmasks
 */
ATLAS_IMAGE *
read_atlas_image(const FITS *fits,	/* the table */
		 int row)		/* the desired row */
{
   return(read_fits_row(fits, row, TYPE_ATLAS_IMAGE));
}

OBJMASK *
read_objmask(const FITS *fits,		/* the table */
		 int row)		/* the desired row */
{
   return(read_fits_row(fits, row, TYPE_OBJMASK));
}

PSF_KL_COMP *
read_KLC(const FITS *fits,		/* the table */
	 int row)			/* the desired row */
{
   return(read_fits_row(fits, row, TYPE_PSF_KL_COMP));
}

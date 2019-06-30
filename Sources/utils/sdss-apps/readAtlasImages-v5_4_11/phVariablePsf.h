#if !defined(PH_VARIABLE_PSF_H)
#define PH_VARIABLE_PSF_H
#include "dervish.h"

/*
 * a struct to contain coefficients describing the kernel required
 * to convert one PSF to another
 */
#define MAX_ORDER_B 5			/* maximum spatial order */
typedef struct {
   int shMalloc;			/* was this ACOEFF shMalloced? */
   int nrow;				/* order of polynomial in row */
   int ncol;				/*                        and column */
   float c[MAX_ORDER_B][MAX_ORDER_B];	/* Kernel coefficient */
} ACOEFF;				/* pragma SCHEMA */

typedef struct {
   REGION *reg;				/* region containing object */
   float rowc, colc;			/* center of object */
   float counts;			/* number of counts in object */
   float countsErr;                     /* counts error */
} PSF_REG;				/* pragma SCHEMA */

typedef struct {
   int refcntr;				/* number of live references to
					   this PSF_KERNEL */
   int border;				/* border to ignore around REGIONs */
   int nsigma;				/* number of Gaussian widths */
   int nrow_a;				/* orders in row and column */
   int ncol_a;				/*   for kernel at a point in image */
   int nrow_b;				/* maximum orders in row and column */
   int ncol_b;				/*   for spatial variation of kernel */
   float *sigma;			/* values of Gaussian sigma */
   ACOEFF ***a;				/* kernel coeffs at a point */
} PSF_KERNEL;				/* pragma SCHEMA */

typedef enum {
   UNKNOWN_BASIS,
   GAUSSIAN_BASIS,
   KL_BASIS
} BASIS_TYPE;				/* pragma SCHEMA */

typedef struct psf_basis {
   BASIS_TYPE type;			/* type of PSF_BASIS */
   const PSF_KERNEL *kern;		/* description of the kernel in use */
   PSF_REG ****regs;			/* basis functions */
   int soft_bias;			/* software bias in regions */
   int Nstars_basis;                    /* # of stars for basis functions */
   int Nframes_basis;                   /* # of frames for collecting
					   basis stars */
   int Nstars_coeffs;                   /* # of stars used for coefficients */
   int Nframes_coeffs;                  /* # of frames used for collecting
					   coeff. stars */
   struct psf_basis *deriv[2 + 3];	/* 1st and 2nd derivatives of basis */
} PSF_BASIS;				/* pragma SCHEMA */

typedef enum {				/* indices into PSF_BASIS->deriv */
   DBASIS_DROW = 0,
   DBASIS_DCOL,
   DBASIS_DROW2,
   DBASIS_DROWDCOL,
   DBASIS_DCOL2
} PSF_BASIS_DERIV;			/* pragma SCHEMA */

typedef struct {
   int nrow_b, ncol_b;			/* how many of c's coeffs are real */
   float c[MAX_ORDER_B][MAX_ORDER_B];	/* spatial variation of coefficients */
   float lambda;			/* eigenvalue */
   REGION *reg;				/* an eigenimage */
   float counts;			/* sum of pixels in reg */
} PSF_KL_COMP;				/* pragma SCHEMA */

/*****************************************************************************/

ACOEFF *phAcoeffNew(int nrow, int ncol);
void phAcoeffDel(ACOEFF *acoeff);

PSF_REG *phPsfRegNew(int nrow, int ncol, PIXDATATYPE type);
void phPsfRegDel(PSF_REG *preg);

PSF_KERNEL *
phPsfKernelNew(int nsigma,		/* number of Gaussian widths */
	       int nrow_a,		/* orders in row and column */
	       int ncol_a,		/*   for kernel at a point in image */
	       int nrow_b,		/* maximum orders in row and column */
	       int ncol_b);		/*   for spatial variation of kernel */
void phPsfKernelDel(PSF_KERNEL *kern);

PSF_BASIS *
phPsfBasisNew(const PSF_KERNEL *kern,	/* description of kernel */
	      BASIS_TYPE basis_type,	/* desired type of PSF_BASIS */
	      int nrow, int ncol,	/* size of desired REGIONs */
	      PIXDATATYPE type);		/* type of desired REGIONs */
void
phPsfBasisDel(PSF_BASIS *basis);

void
phPsfBasisSet(PSF_BASIS *basis);

PSF_BASIS *
phRegConvolveWithPsfKernel(const REGION *reg, /* region to convolve */
			   PSF_KERNEL *kern, /* kernel to use */
			   PSF_BASIS *basis, /* basis functions, or NULL */
			   PIXDATATYPE type); /* type of basis REGIONs or -1 */

int
phRegionDotRegion(double *dot,		/* desired value of dot product */
		  const REGION *reg1,	/* one region */
		  const REGION *reg2,	/* the other region */
		  int border);		/* width of border to ignore */

PSF_KERNEL *
phPsfKernelSet(PSF_KERNEL *kern,	/* desired kernel, or NULL */
	       const CHAIN *regs,	/* regions to be fit */
	       const PSF_BASIS *basis);	/* basis functions */
REGION *
phPsfKernelApply(const PSF_REG *preg,	/* object to be transformed */
		 const PSF_KERNEL *kern, /* kernel to apply */
		 REGION *iscr1,		/* scratch FL32 region, or NULL */
		 REGION *iscr2,		/* scratch FL32 region, or NULL */
		 REGION *iscr3);	/* scratch FL32 region, or NULL */

PSF_BASIS *
phPsfKLDecomp(const CHAIN *regs,	/* the PSF_REGs with input stars */
	      int border,		/* how many pixels to ignore
					   around regions */
	      int ncomp,		/* number of components to keep */
	      int nrow_b,		/* maximum orders in row and column */
	      int ncol_b);		/*   for spatial variation of kernel */


PSF_REG *
phPsfKLReconstruct(const PSF_BASIS *basis, /* basis to use */
		   float rowc,		/* location of */
		   float colc,		/*    desired PSF */
		   PIXDATATYPE regType); /* type of desired REGION */
		   

PSF_KL_COMP *phPsfKLCompNew(REGION *reg);
void phPsfKLCompDel(PSF_KL_COMP *klc);

PSF_KL_COMP *
phPsfKLCompSetFromBasis(PSF_KL_COMP *klc, /* PSF_KL_COMP to set, or NULL */
			const PSF_BASIS *basis,	/* basis with desired data */
			int comp);		/* desired component */
void
phPsfBasisSetFromKLComp(PSF_BASIS *basis, /* basis to be set */
			int comp,	/* component to set*/
			const PSF_KL_COMP *klc, /* PSF_KL_COMP with data */
			int copy_reg);	/* copy region? */

PSF_BASIS *
phPsfBasisDifferentiate(const PSF_BASIS *basis,	/* input basis */
			REGION *scr,	/* scratch region, or NULL */
			int nderiv_row,	/* order of derivative in row */
			int nderiv_col, /*            and column directions */
			int lbell);	/* desired half-length of cosbell */

void
phPsfBasisSetDerivatives(PSF_BASIS *basis, /* basis to differentiate */
			 REGION *scr,	/* scratch region, or NULL */
			 int lbell);	/* desired half-length of cosbell */
#endif

/*
 * Code to handle variable PSFs
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <alloca.h>
#include "dervish.h"
#if defined(STAND_ALONE)
   static int shRegIntCopy(REGION *tmp_reg, const REGION *reg);
#else
#  include "atConversions.h"		/* for M_PI */
#  include "phUtils.h"
#  include "phSkyUtils.h"
#  include "phMeschach.h"
#endif
#include "phConsts.h"
#include "phVariablePsf.h"
/*
 * Set the pointers for a 3-dimensional array BASE of objects of type TYPE,
 * with dimensions NSIGMA, NROW_A, and NCOL_A
 */
#define RC_SCALE 1e-3			/* scale factor for rowc/colc coeffs */

#define SET_POINTERS(BASE, TYPE, NSIGMA, NROW_A, NCOL_A) \
   BASE = shMalloc((NSIGMA + 1)*sizeof(TYPE **)); \
   BASE[0] = shMalloc((1 + NSIGMA*NROW_A)*sizeof(TYPE *)); \
   BASE[0][0] = shMalloc((1 + NSIGMA*NROW_A*NCOL_A)*sizeof(TYPE)); \
   \
   for(is = 1; is <= NSIGMA; is++) { \
      BASE[is] = BASE[0] + 1 + NROW_A*(is - 1); \
      \
      BASE[is][0] = BASE[0][0] + 1 + NROW_A*NCOL_A*(is - 1); \
      for(ira = 0; ira < NROW_A; ira++) { \
	 BASE[is][ira] = BASE[is][0] + NCOL_A*ira; \
      } \
   }

/*****************************************************************************/
/*
 * Create/Delete an ACOEFF; not used by PSF_KERNEL for historical/efficiency
 * regions
 */
ACOEFF *
phAcoeffNew(int nrow, int ncol)
{
   ACOEFF *acoeff = shMalloc(sizeof(ACOEFF));

   shAssert(nrow <= MAX_ORDER_B && ncol <= MAX_ORDER_B);

   acoeff->shMalloc = 1;
   acoeff->nrow = nrow;
   acoeff->ncol = ncol;

   return(acoeff);
}

void
phAcoeffDel(ACOEFF *acoeff)
{
   if(acoeff == NULL) return;

   shAssert(acoeff->shMalloc);
   shFree(acoeff);
}

/*****************************************************************************/
/*
 * Create a PSF_REG
 */
PSF_REG *
phPsfRegNew(int nrow, int ncol,		/* size of desired region */
	    PIXDATATYPE type)		/* type of desired REGION */
{
   PSF_REG *preg = shMalloc(sizeof(PSF_REG));
   preg->reg = (nrow < 0) ? NULL : shRegNew("PSF_REG", nrow, ncol, type);
   preg->rowc = preg->colc = VALUE_IS_BAD;
   preg->counts = preg->countsErr = VALUE_IS_BAD;

   return(preg);
}

void
phPsfRegDel(PSF_REG *preg)
{
   if(preg == NULL) {
      return;
   }
   
   shRegDel(preg->reg);
   shFree(preg);
}

/*****************************************************************************/
/*
 * Create a new PSF_KERNEL.  Note that the arrays start at -1, corresponding
 * to a delta function; there are no associated row^i col^j N(0,s^2) terms
 */
PSF_KERNEL *
phPsfKernelNew(int nsigma,		/* number of Gaussian widths */
	       int nrow_a,		/* orders in row and column */
	       int ncol_a,		/*   for kernel at a point in image */
	       int nrow_b,		/* maximum orders in row and column */
	       int ncol_b)		/*   for spatial variation of kernel */
{
   PSF_KERNEL *kern = shMalloc(sizeof(PSF_KERNEL));
   int is, ira, ica;

   shAssert(nsigma >= 0);
   shAssert(nrow_a > 0 && ncol_a > 0);
   shAssert(nrow_b > 0 && ncol_b > 0 && \
	    nrow_b <= MAX_ORDER_B && ncol_b <= MAX_ORDER_B);

   kern->refcntr = 1;
   kern->nsigma = nsigma;
   kern->nrow_a = nrow_a;
   kern->ncol_a = ncol_a;
   kern->nrow_b = nrow_b;
   kern->ncol_b = ncol_b;
   
   kern->border = 0;
   kern->sigma = shMalloc((nsigma + 1)*sizeof(float));
   kern->sigma[0] = 0.0;		/* a delta function */

   SET_POINTERS(kern->a, ACOEFF, nsigma, nrow_a, ncol_a);

   kern->sigma++;			/* put delta-fn at [-1] */
   kern->a++;				/* put delta-fn at [-1] */

   kern->a[-1][0][0].shMalloc = 0;	/* _not_ shMalloced */
   kern->a[-1][0][0].c[0][0] = 0;
   kern->a[-1][0][0].nrow = kern->a[-1][0][0].ncol = 1;
   for(is = 0; is < nsigma; is++) {
      kern->sigma[is] = 0.0;
      for(ira = 0; ira < nrow_a; ira++) {
	 for(ica = 0; ica < ncol_a; ica++) {
	    kern->a[is][0][0].shMalloc = 0; /* _not_ shMalloced */
	    kern->a[is][ira][ica].nrow = nrow_a;
	    kern->a[is][ira][ica].ncol = ncol_a;
	    kern->a[is][ira][ica].c[0][0] = 0;
	 }
      }
   }
   
   return(kern);
}

void
phPsfKernelDel(PSF_KERNEL *kern)
{
   if(kern == NULL) {
      return;
   }

   if(--kern->refcntr > 0) {		/* still in use somewhere else */
      return;
   }
   
   shAssert(kern->sigma != NULL);
   kern->sigma--;			/* we ++ed it */
   shFree(kern->sigma);

   shAssert(kern->a != NULL && kern->a[-1] != NULL && kern->a[-1][0] != NULL);
   kern->a--;				/* we ++ed it */
   shFree(kern->a[0][0]);
   shFree(kern->a[0]);
   shFree(kern->a);

   shFree(kern);
}

/*****************************************************************************/
/*
 * Now for a set of REGIONs corresponding to all those basis functions
 * 
 * Note that there's a REGION at regs[-1], corresponding to a delta function;
 * there are no associated row^i col^j N(0,s^2) REGIONs
 */
PSF_BASIS *
phPsfBasisNew(const PSF_KERNEL *kern,	/* description of kernel */
	      BASIS_TYPE basis_type,	/* desired type of PSF_BASIS */
	      int nrow, int ncol,	/* size of desired REGIONs */
	      PIXDATATYPE type)		/* type of desired REGIONs */
{
   PSF_BASIS *basis = shMalloc(sizeof(PSF_BASIS));
   int is, ira, ica;

   shAssert(kern != NULL);
   shAssert(nrow >= 0 || ncol < 0);

   if(basis_type == UNKNOWN_BASIS ||
      basis_type == GAUSSIAN_BASIS ||
      basis_type == KL_BASIS) {
   } else {
      shError("phPsfBasisNew: unknown BASIS_TYPE: %d", basis_type);
      basis_type = UNKNOWN_BASIS;
   }
   basis->type = basis_type;
   basis->kern = kern; ((PSF_KERNEL *)kern)->refcntr++;
   basis->Nstars_basis = -1;
   basis->Nframes_basis = -1;      
   basis->Nstars_coeffs = -1;                  
   basis->Nframes_coeffs = -1;                
   
   SET_POINTERS(basis->regs, PSF_REG *, \
		kern->nsigma, kern->nrow_a, kern->ncol_a);
   basis->regs++;			/* put delta-fn at [-1] */
   
   switch (type) {
    case TYPE_U8: basis->soft_bias = 100; break;
    case TYPE_S8: basis->soft_bias = 0; break;
    case TYPE_U16: basis->soft_bias = 1000; break;
    case TYPE_S16: basis->soft_bias = 0; break;
    case TYPE_S32: basis->soft_bias = 0; break;
    case TYPE_FL32: basis->soft_bias = 0; break;
    default:
      shError("phPsfBasisNew: unknown REGION type %d", type);
      basis->soft_bias = 0; break;
   }

   basis->regs[-1][0][0] = phPsfRegNew(nrow, ncol, type);
   for(is = 0; is < kern->nsigma; is++) {
      for(ira = 0; ira < kern->nrow_a; ira++) {
	 for(ica = 0; ica < kern->ncol_a; ica++) {
	    basis->regs[is][ira][ica] = phPsfRegNew(nrow, ncol, type);
	 }
      }
   }

   for(is = 0; is < 5; is++) {
      basis->deriv[is] = NULL;
   }
   
   return(basis);
}

void
phPsfBasisDel(PSF_BASIS *basis)
{
   int i;
   int nreg;
   
   if(basis == NULL) {
      return;
   }

   shAssert(basis->kern != NULL);   
   shAssert(basis->regs != NULL && basis->regs[-1] != NULL && \
	    basis->regs[-1][0] != NULL);

   nreg = 1 + basis->kern->nsigma*basis->kern->nrow_a*basis->kern->ncol_a;
   for(i = 0; i < nreg; i++) {
      phPsfRegDel(basis->regs[-1][0][i]);
   }
   
   basis->regs--;			/* undo the ++ in the constructor */
   shFree(basis->regs[0][0]);
   shFree(basis->regs[0]);
   shFree(basis->regs);

   phPsfKernelDel((PSF_KERNEL *)basis->kern);

   for(i = 0; i < 5; i++) {
      phPsfBasisDel(basis->deriv[i]);
   }

   shFree(basis);
}

#if !defined(STAND_ALONE)
/*****************************************************************************/
/*
 * set a REGION to a row^ir * col^ic * Gaussian
 */
static void
kernel_set(REGION *reg,			/* REGION to set */
	   int amp,			/* central amplitude */
	   int bkgd,			/* background level */
	   float sigma,			/* width parameter for Gaussian */
	   int ir,			/* order in row direction */
	   int ic)			/* order in column direction */
{
   float drow, drow2, dcol, dcol2;	/* (row-rowc), drow^2, etc. */
   PIX **rows = reg->ROWS;		/* == reg->ROWS */
   PIX *row;				/* == reg->ROWS[] */
   const int ncol = reg->ncol;
   const int nrow = reg->nrow;
   int r, c;
   float rmult, cmult;			/* coefficients for Gaussian */
   const float isigma = 1.0/sigma;
   const float rowc = nrow/2, colc = ncol/2; /* centre of REGION */

   for(r = 0; r < nrow; r++) {
      row = rows[r];
      drow = r - rowc;
      drow2 = drow*drow;
      rmult = amp*pow(isigma*drow, ir);
      for(c = 0; c < ncol; c++) {
	 dcol = c - colc;
	 dcol2 = dcol*dcol;
	 cmult = pow(isigma*dcol, ic);
	 row[c] = bkgd + rmult*cmult*exp(-0.5*isigma*isigma*(drow2 + dcol2));
      }
   }
}

/*****************************************************************************/
/*
 * calculate a set of basis functions
 */
void
phPsfBasisSet(PSF_BASIS *basis)
{
   float amp = 30000;
   int is;				/* counter for sigma */
   int ira, ica;			/* counter for row/column coeffs */
   float sigma;				/* == basis->kern->sigma[] */

   shAssert(basis != NULL && basis->kern != NULL);
   shAssert(basis->regs != NULL &&basis->regs[0][0][0]->reg->type == TYPE_PIX);

   shRegIntSetVal(basis->regs[-1][0][0]->reg, basis->soft_bias);
   {
      int nrow = basis->regs[-1][0][0]->reg->nrow/2;
      int ncol = basis->regs[-1][0][0]->reg->ncol/2;
      basis->regs[-1][0][0]->reg->rows_fl32[nrow/2][ncol/2] += 10000;
   }
   for(is = 0; is < basis->kern->nsigma; is++) {
      sigma = basis->kern->sigma[is];
      for(ira = 0; ira < basis->kern->nrow_a; ira++) {
	 for(ica = 0; ica < basis->kern->ncol_a; ica++) {
	    shRegIntSetVal(basis->regs[is][ira][ica]->reg, basis->soft_bias);
	    kernel_set(basis->regs[is][ira][ica]->reg, basis->soft_bias, amp,
		       sigma, ira, ica);
	 }
      }
   }
}

/*****************************************************************************/
/*
 * Convolve a REGION with a set of basis functions, returning a PSF_BASIS;
 * if the provided PSF_BASIS is non-NULL it'll be reused. In this case,
 * its regions must be of the same size and type as reg
 */
PSF_BASIS *
phRegConvolveWithPsfKernel(const REGION *reg, /* region to convolve */
			   PSF_KERNEL *kern, /* kernel to use */
			   PSF_BASIS *basis, /* basis functions, or NULL */
			   PIXDATATYPE type) /* type of basis REGIONs or -1 */
{
   int filtsize = 11;			/* size of convolution filter */
   int is;				/* counter for sigma */
   int ira, ica;			/* counter for row/column coeffs */
   REGION *scr;				/* scratch region for convolutions */
   float sigma;				/* == basis->kern->sigma[] */

   shAssert(reg != NULL);

   if(basis == NULL) {
      basis = phPsfBasisNew(kern, GAUSSIAN_BASIS, reg->nrow, reg->ncol, type);
   } else {
      shAssert(basis->regs != NULL);
      shAssert(type == (PIXDATATYPE)-1 ||
	       type == basis->regs[0][0][0]->reg->type);
      shAssert(basis->regs[0][0][0]->reg->nrow == reg->nrow && \
	       basis->regs[0][0][0]->reg->ncol == reg->ncol && \
	       (type == TYPE_U16 || type == TYPE_FL32));

      phPsfKernelDel((PSF_KERNEL *)basis->kern);
      basis->kern = kern; kern->refcntr++;
      type = basis->regs[0][0][0]->reg->type;
   }

   scr = shRegNew("scratch", reg->nrow, reg->ncol, type);

   shRegIntCopy(basis->regs[-1][0][0]->reg, reg);
   for(is = 0; is < basis->kern->nsigma; is++) {
      sigma = basis->kern->sigma[is];
      for(ira = 0; ira < basis->kern->nrow_a; ira++) {
	 for(ica = 0; ica < basis->kern->ncol_a; ica++) {
	    phConvolveWithGaussianPowers(basis->regs[is][ira][ica]->reg,
					 reg, scr, filtsize, sigma, ira, ica);
	 }
      }
   }

   shRegDel(scr);

   return(basis);
}

/*****************************************************************************/
/*
 * Calculate the dot product of two REGIONs of the same size; for U16
 * REGIONs, subtract the SOFT_BIAS that's assumed to be present
 *
 * Return 0, or -1 in case of problems
 */
int
phRegionDotRegion(double *dot,		/* desired value */
		  const REGION *reg1,	/* one region */
		  const REGION *reg2,	/* the other region */
		  int border)		/* width of border to ignore */
{
   int i, j;
   int nrow, ncol;			/* == reg->n{row,col} */
   int row_0 = border, col_0 = border;	/* LLC and */
   int row_1 = -border, col_1 = -border; /*        URC of region to use */
   double sum = 0.0;			/* the desired value */
   
   shAssert(reg1 != NULL && reg2 != NULL);
   shAssert(reg1->nrow == reg2->nrow && reg1->ncol == reg2->ncol);
   nrow = reg1->nrow; ncol = reg1->ncol;
   row_1 += nrow - 1;			/* specified wrt top right */
   col_1 += ncol - 1;
   shAssert(row_0 >= 0 && col_0 >= 0 && row_1 < nrow && col_1 < ncol);

   *dot = 0.0;
   if(reg1->type == TYPE_U16) {
      if(reg2->type == TYPE_U16) {
	 U16 **rows1 = reg1->rows_u16, *row1;
	 U16 **rows2 = reg2->rows_u16, *row2;
	 for(i = row_0; i < row_1; i++) {
	    row1 = rows1[i]; row2 = rows2[i];
	    for(j = col_0; j < col_1; j++) {
	       sum += (row1[j] - SOFT_BIAS)*(row2[j] - SOFT_BIAS);
	    }
	 }
      } else if(reg2->type == TYPE_FL32) {
	 U16 **rows1 = reg1->rows_u16, *row1;
	 FL32 **rows2 = reg2->rows_fl32, *row2;
	 for(i = row_0; i < row_1; i++) {
	    row1 = rows1[i]; row2 = rows2[i];
	    for(j = col_0; j < col_1; j++) {
	       sum += (row1[j] - SOFT_BIAS)*row2[j];
	    }
	 }
      } else {
	 return(-1);
      }
   } else if(reg1->type == TYPE_FL32) {
      if(reg2->type == TYPE_U16) {
	 return(phRegionDotRegion(dot, reg2, reg1, border)); /* switch order */
      } else if(reg2->type == TYPE_FL32) {
	 FL32 **rows1 = reg1->rows_fl32, *row1;
	 FL32 **rows2 = reg2->rows_fl32, *row2;
	 for(i = row_0; i < row_1; i++) {
	    row1 = rows1[i]; row2 = rows2[i];
	    for(j = col_0; j < col_1; j++) {
	       sum += row1[j]*row2[j];
	    }
	 }
      } else {
	 return(-1);
      }
   }

   *dot = sum;

   return(0);
}

/*****************************************************************************/
/*
 * Given a REGION, a chain of PSF_REGs, and a PSF_BASIS with its regs fields
 * filled out, set the a coefficients in a PSF_KERNEL
 */
PSF_KERNEL *
phPsfKernelSet(PSF_KERNEL *kern,	/* desired kernel, or NULL */
	       const CHAIN *regs,	/* the PSF_REGs with input stars */
	       const PSF_BASIS *basis)	/* basis functions */
{
   MAT *A;				/* normal equations for LSQ fit are */
   VEC *b;				/*      A*w = b */
   MAT *A0;				/* A without kernel variation */
   VEC *b0;				/* b without kernel variation */
   int border;				/* == kern->border */
   float *coeffs;			/* coefficients for spatial variation*/
   double dot;				/* inner product of two regions */
   int i, j;
   int k0, l0;				/* starting indices for i and j
					   in a sub-block of A */
   int k, l;				/* counters over spatial terms */
   int ireg;				/* index for elements of regs */
   VEC *lambda;				/* eigen values */
   float norm;				/* normalisation constant */
   int nparam;				/* total number of parameters to fit */
   int nparam0;				/* no. of parameters to fit per star */
   int nreg;				/* number of regions in regs */
   int nsigma, nrow_a, ncol_a, nrow_b, ncol_b; /* == kern->nsigma etc. */
   int nrow, ncol;			/* == reg->n{row,col} */
   MAT *Q;				/* eigen vectors */
   const PSF_REG *psf_reg;		/* an element of regs */
   const REGION *reg;			/* region containing an object */
   float rowc, colc;			/* position of a star in regs */
#if 0
   FL32 *row;				/* a row of the image */
   int row_0, col_0;			/* LLC and URC of part of regions */
   int row_1, col_1;			/*         use from PSF_BASIS->regs */
   double sum;				/* sum of weight*pixel forall pixels*/
   double sum_reg;			/* sum of all pixels in a region */
#endif
   VEC *w;				/* desired weights */

   shAssert(regs != NULL && regs->type == shTypeGetFromName("PSF_REG"));
   nreg = regs->nElements;
   if(nreg == 0) {
      return(NULL);
   }

   reg = ((PSF_REG *)shChainElementGetByPos(regs, 0))->reg;
   shAssert(reg != NULL && (reg->type == TYPE_U16 || reg->type == TYPE_FL32));
   nrow = reg->nrow; ncol = reg->ncol;
   shAssert(basis != NULL && basis->kern != NULL && basis->regs != NULL);
   shAssert(basis->regs[-1][0][0]->reg->type == TYPE_FL32);
   shAssert(basis->regs[-1][0][0]->reg->nrow == nrow && \
	    basis->regs[-1][0][0]->reg->ncol == ncol);

   if(kern == NULL) {
      kern = phPsfKernelNew(basis->kern->nsigma,
			    basis->kern->nrow_a, basis->kern->ncol_a,
			    basis->kern->nrow_b, basis->kern->ncol_b);
      memcpy(kern->sigma, basis->kern->sigma,
				      (basis->kern->nsigma + 1)*sizeof(float));
   } else {
      shAssert(kern->nsigma == basis->kern->nsigma && \
	       kern->nrow_a == basis->kern->nrow_a && \
	       kern->ncol_a == basis->kern->ncol_a && \
	       kern->nrow_b == basis->kern->nrow_b && \
	       kern->ncol_b == basis->kern->ncol_b);
   }
   nsigma = kern->nsigma;
   nrow_a = kern->nrow_a;
   ncol_a = kern->ncol_a;
   nrow_b = kern->nrow_b;
   ncol_b = kern->ncol_b;
/*
 * prepare to setup the LSQ problem
 */
   nparam0 = 1 + nsigma*nrow_a*ncol_a;
   nparam = nparam0*nrow_b*ncol_b;
   
   A = phMatNew(nparam, nparam);
   b = phVecNew(nparam);
   A0 = phMatNew(nparam0, nparam0);
   b0 = phVecNew(nparam0);
   lambda = NULL;
   Q = phMatNew(nparam, nparam);
/*
 * and fill out the matrices for the normal equations, Aw = b.
 * Note that it is required that the memory for all the basis->regs
 * be contiguous, as we check with an assertion
 */
   shAssert(basis->regs[-1][0][nparam0 - 1] == \
	    basis->regs[nsigma - 1][nrow_a - 1][ncol_a - 1]);

   border = basis->kern->border;
/*
 * Clear A and b as we'll be accumulating into them
 */
   for(i = 0; i < nparam; i++) {
      b->ve[i] = 0;
      for(j = i; j < nparam; j++) {
	 A->me[i][j] = A->me[j][i] = 0;
      }
   }
/*
 * First set A0, whose elements are the (K_i x P_0).(K_j x P_0) terms,
 * i.e. the matrix A in the absence of spatial variation of the psf 
 */
   for(i = 0; i < nparam0; i++) {
      for(j = i; j < nparam0; j++) {
	 phRegionDotRegion(&dot, basis->regs[-1][0][i]->reg,
			   basis->regs[-1][0][j]->reg, border);
	 A0->me[j][i] = A0->me[i][j] = dot;
      }
   }
/*
 * Set a convenience array coeffs with the powers of rowc and colc to be used
 * for each of the nrow_b*ncol_b spatial terms
 */
   coeffs = alloca(nrow_b*ncol_b*sizeof(float));
/*
 * We need to sum over each region, filling out A and b. The elements of
 * A only depend on the position and the basis functions being used, so
 * we can derive them from the A0 matrix; b's elements also depend upon
 * the PSF at that point, so we have to recalculate b0 for each position
 */
   for(ireg = 0; ireg < nreg; ireg++) {   
      psf_reg = shChainElementGetByPos(regs, ireg);
      reg = psf_reg->reg;
      rowc = psf_reg->rowc;
      colc = psf_reg->colc;

      for(k = 0; k < nrow_b*ncol_b; k++) {
	 coeffs[k] = pow(rowc*RC_SCALE, k%nrow_b)*pow(colc*RC_SCALE, k/nrow_b);
      }

      for(i = 0; i < nparam0; i++) {
	 phRegionDotRegion(&dot, reg, basis->regs[-1][0][i]->reg, border);
	 b0->ve[i] = dot/psf_reg->counts;
      }

      for(i = k0 = 0; i < nparam0; i++, k0 += nrow_b*ncol_b) {
	 for(k = 0; k < nrow_b*ncol_b; k++) {
	    b->ve[k0 + k] += coeffs[k]*b0->ve[i];
	 }
      }
      
      for(i = k0 = 0; i < nparam0; i++, k0 += nrow_b*ncol_b) {
	 for(j = l0 = 0; j < nparam0; j++, l0 += nrow_b*ncol_b) {
	    for(k = 0; k < nrow_b*ncol_b; k++) {
	       for(l = 0; l < nrow_b*ncol_b; l++) {
		  A->me[k0 + k][l0 + l] += coeffs[k]*coeffs[l]*A0->me[i][j];
	       }
	    }
	 }
      }
   }
/*
 * solve the system, replacing any eigenvalues that are too small with 0,
 * which phEigenBackSub() will interpret as infinity.
 */
   lambda = phEigen(A, Q, lambda);

   for(i = 0; i < nparam; i++) {
      if(fabs(lambda->ve[i]) < 1e-6) {
	 lambda->ve[i] = 0.0;
      }
   }

   w = phEigenBackSub(Q, lambda, b);
#if 0
/*
 * normalise the resulting transformation kernel to unit sum
 */
   row_0 = border; col_0 = border;
   row_1 = nrow - border - 1; col_1 = ncol - border - 1;

   sum = 0;
   for(i = 0; i < nparam0; i++) {
      sum_reg = 0;
      for(j = row_0; j <= row_1; j++) {
	 row = basis->regs[-1][0][i]->reg->rows_fl32[j];
	 for(k = col_0; k <= col_1; k++) {
	    sum_reg += row[k];
	 }
      }
      coeff_k = 0;
      for(k = 0; k < nrow_b*ncol_b; k++) {
	 coeff_k += w->ve[k*nparam0 + i];
      }
      sum += coeff_k*sum_reg;
   }
/*
 * now find the normalisation of the initial image
 */
   sum_reg = 0;
   if(reg->type == TYPE_U16) {
      U16 *row_u16;
      for(j = row_0; j <= row_1; j++) {
	 row_u16 = reg->rows_u16[j];
	 for(k = col_0; k <= col_1; k++) {
	    sum_reg += row_u16[k];
	 }
      }
      sum_reg -= SOFT_BIAS*(row_1 - row_0 + 1)*(col_1 - col_0 + 1);
   } else if(reg->type == TYPE_FL32) {
      for(j = row_0; j <= row_1; j++) {
	 row = reg->rows_fl32[j];
	 for(k = col_0; k <= col_1; k++) {
	    sum_reg += row[k];
	 }
      }
   } else {
      shFatal("You cannot get here");
   }

   shAssert(sum != 0.0);
   norm = sum_reg/sum;
#else
   norm = 1;
#endif
/*
 * pack the result into the PSF_KERNEL
 */
   kern->border = basis->kern->border;
   for(i = -1; i < nsigma; i++) {
      kern->sigma[i] = basis->kern->sigma[i];
   }

   for(i = k0 = 0; i < nparam0; i++, k0 += nrow_b*ncol_b) {
      for(k = 0; k < nrow_b*ncol_b; k++) {
	 kern->a[-1][0][i].c[k%nrow_b][k/nrow_b] = w->ve[k0 + k]*norm;
      }
   }
/*
 * clean up
 */
   phMatDel(A);
   phVecDel(b);
   phMatDel(A0);
   phVecDel(b0);
   phVecDel(lambda);
   phMatDel(Q);
   phVecDel(w);

   return(kern);
}

/*****************************************************************************/
/*
 * Given a PSF_REG <preg> and a PSF_KERNEL with its coefficients set, return
 * a REGION consisting of <in> convolved with that kernel
 */
REGION *
phPsfKernelApply(const PSF_REG *preg,	/* object to be transformed */
		 const PSF_KERNEL *kern, /* kernel to apply */
		 REGION *iscr1,		/* scratch FL32 region, or NULL */
		 REGION *iscr2,		/* scratch FL32 region, or NULL */
		 REGION *iscr3)		/* scratch FL32 region, or NULL */
{
   float acoeff;			/* the "a" coefficient for component */
   float *coeffs;			/* coefficients for spatial variation*/
   int filtsize = 11;			/* size of convolution filter */
   REGION *flout;			/* FL32 version of returned region */
   int is, ira, ica;			/* counters in sigma etc. direction */
   int i, j, k;
   int nrow, ncol;			/* == reg->n{rol,col} */
   int nsigma, nrow_a, ncol_a, nrow_b, ncol_b; /* == kern->nsigma etc. */
   REGION *out;				/* returned REGION */
   const REGION *reg;			/* the input REGION itself */
   float rowc, colc;			/* == preg->{row,col}c */
   FL32 **rows_flout, *row_flout;	/* == flout->rows, flout->rows[] */
   REGION *smreg;			/* Convolution of reg with filter */
   FL32 **rows_smreg, *row_smreg;	/* == smreg->rows, smreg->rows[] */
   REGION *scr;				/* scratch FL32 REGION */
   float sigma;				/* == kern->sigma[] */
   REGION *tmp_reg = NULL;		/* FL32 copy of reg, or NULL */

   shAssert(preg != NULL);
   reg = preg->reg; rowc = preg->rowc; colc = preg->colc;
   shAssert(reg != NULL);
   if(reg->type == TYPE_U16) {		/* make an FL32 copy */
      tmp_reg = shRegNew("", reg->nrow, reg->ncol, TYPE_FL32);
      shRegIntCopy(tmp_reg, reg);
      reg = tmp_reg;
   }
   
   shAssert(reg != NULL && reg->type == TYPE_FL32);
   nrow = reg->nrow; ncol = reg->ncol;
   shAssert(kern != NULL);
   nsigma = kern->nsigma;
   nrow_a = kern->nrow_a;
   ncol_a = kern->ncol_a;
   nrow_b = kern->nrow_b;
   ncol_b = kern->ncol_b;

   out = shRegNew("", nrow, ncol, reg->type);
   if(iscr1 == NULL) {
      flout = shRegNew("", nrow, ncol, TYPE_FL32);
   } else {
      shAssert(iscr1->nrow == nrow && iscr1->ncol == ncol &&
	       iscr1->type == TYPE_FL32);
      flout = iscr1;
   }
   rows_flout = flout->rows_fl32;

   if(iscr2 == NULL) {
      smreg = shRegNew("", nrow, ncol, TYPE_FL32);
   } else {
      shAssert(iscr2->nrow == nrow && iscr2->ncol == ncol &&
	       iscr2->type == TYPE_FL32);
      smreg = iscr2;
   }
   rows_smreg = smreg->rows_fl32;
   
   if(iscr3 == NULL) {
      scr = shRegNew("", nrow, ncol, TYPE_FL32);
   } else {
      shAssert(iscr3->nrow == nrow && iscr3->ncol == ncol &&
	       iscr3->type == TYPE_FL32);
      scr = iscr3;
   }
/*
 * Set a convenience array coeffs with the powers of rowc and colc to be used
 * for each of the nrow_b*ncol_b spatial terms
 */
   coeffs = alloca(nrow_b*ncol_b*sizeof(float));

   for(k = 0; k < nrow_b*ncol_b; k++) {
      coeffs[k] = pow(rowc*RC_SCALE, k%nrow_b)*pow(colc*RC_SCALE, k/nrow_b);
   }

   shRegClear(flout);
   for(is = -1; is < nsigma; is++) {
      sigma = kern->sigma[is];
      for(ira = 0; ira < nrow_a; ira++) {
	 for(ica = 0; ica < ncol_a; ica++) {
	    if(is < 0) {
	       shRegIntCopy(smreg, reg);
	    } else {
	       phConvolveWithGaussianPowers(smreg,
					 reg, scr, filtsize, sigma, ira, ica);
	    }
/*
 * add that into (floating) output region
 */
	    acoeff = 0;
	    for(k = 0; k < nrow_b*ncol_b; k++) {
	       acoeff += kern->a[is][ira][ica].c[k%nrow_b][k/nrow_b]*coeffs[k];
	    }
	    
	    for(i = 0; i < nrow; i++) {
	       row_flout = rows_flout[i];
	       row_smreg = rows_smreg[i];
	       for(j = 0; j < ncol; j++) {
		  row_flout[j] += acoeff*row_smreg[j];
	       }
	    }

	    if(is < 0) break;		/* only one component for is==0 */
	 }
	 
	 if(is < 0) break;		/* only one component for is==0 */
      }
   }

   shRegIntCopy(out, flout);

   shRegDel(tmp_reg);
   if(flout != iscr1) {
      shRegDel(flout);
   }
   if(smreg != iscr2) {
      shRegDel(smreg);
   }
   if(scr != iscr3) {
      shRegDel(scr);
   }

   return(out);
}

/*****************************************************************************/
/*
 * Return a PSF_BASIS consisting of the first ncomp elements of the
 * Karhunen-Lo\`eve basis of regs
 *
 * The notation is that in chapter 7 of Gyula Szokoly's thesis at JHU
 */
PSF_BASIS *
phPsfKLDecomp(const CHAIN *regs,	/* the PSF_REGs with input stars */
	      int border,		/* how many pixels to ignore
					   around regions */
	      int ncomp,		/* number of components to keep */
	      int nrow_b,		/* maximum orders in row and column */
	      int ncol_b)		/*   for spatial variation of kernel */
{
   PSF_BASIS *basis;			/* the desired basis */
   double dot;				/* inner product of two regions */
   int i, j, k;
   VEC *lambda;				/* eigenvalues */
   int nreg;				/* == regs->nElements */
   int nrow, ncol;			/* == reg_i->n{rol,col} */
   PSF_REG *preg;			/* an element of regs */
   MAT *Q;				/* eigen vectors */
   MAT *R;				/* residuals' inner products */
   REGION *reg_i, *reg_j;		/* == preg->reg */
  
   shAssert(regs != NULL);
   nreg = regs->nElements;
   shAssert(regs->type == shTypeGetFromName("PSF_REG") && nreg > 0);
   preg = shChainElementGetByPos(regs, 0);
   reg_i = preg->reg;
   shAssert(reg_i != NULL);

   nrow = reg_i->nrow; ncol = reg_i->ncol;
/*
 * Create the PSF_BASIS
 */
   {
      PSF_KERNEL *kern = phPsfKernelNew(ncomp - 1, 1, 1, nrow_b, ncol_b);
      kern->border = border;
      basis = phPsfBasisNew(kern, KL_BASIS, nrow, ncol, TYPE_FL32);
      phPsfKernelDel(kern);		/* decrement reference counter */
   }
#if 0
/*
 * Estimate the mean background level; value is not used except to check
 * sky subtraction for input regions (using gdb)
 */
      {
	 float sky = 0;			/* estimate of sky level */

	 for(i = 0; i < nreg; i++) {
	    preg = shChainElementGetByPos(regs, i);
	    reg_i = preg->reg;
	    shAssert(reg_i->type == TYPE_PIX);
	    
	    for(j = border + 1; j < nrow - border - 1; j++) {
	       sky += reg_i->ROWS[j][border] +
					     reg_i->ROWS[j][ncol - border - 1];
	    }
	    for(k = border; k < ncol - border; k++) {
	       sky += reg_i->ROWS[border][k] +
					     reg_i->ROWS[nrow - border - 1][k];
	    }
	 }
	 sky /= nreg*(2*((nrow - 2*border) + (ncol - 2*border)) - 4);
	 sky -= SOFT_BIAS;
      }
#endif
/*
 * Find the eigenvectors/values of the scalar product matrix, R' (Eq. 7.4)
 */
   lambda = NULL;
   R = phMatNew(nreg, nreg);
   Q = phMatNew(nreg, nreg);

   for(i = 0; i < nreg; i++) {
      preg = shChainElementGetByPos(regs, i);
      reg_i = preg->reg;
#if 0					/* dump input regions */
      {
	 char dump_filename[100];
	 sprintf(dump_filename, "reg%d.fts", i);
      	 shRegWriteAsFits(reg_i,
			  dump_filename, STANDARD, 2, DEF_NONE, NULL, 0);
      }
#endif
      for(j = i; j < nreg; j++) {
	 preg = shChainElementGetByPos(regs, j);
	 reg_j = preg->reg;

	 phRegionDotRegion(&dot, reg_i, reg_j, border);
	 R->me[i][j] = R->me[j][i] = dot/nreg;
      }
   }

   lambda = phEigen(R, Q, lambda);
/*
 * Contruct the first ncomp eigenimages in basis
 */
   for(i = 0; i < ncomp; i++) {
      basis->kern->sigma[i-1] = lambda->ve[i];
      shRegClear(basis->regs[i-1][0][0]->reg);
      if(i >= nreg) {
	 continue;
      }

#if 0					/* debugging only; show the stars */
      preg = shChainElementGetByPos(regs, i);
      shRegIntCopy(basis->regs[i-1][0][0]->reg, preg->reg);
      continue;      
#endif
      
      for(j = 0; j < nreg; j++) {
	 preg = shChainElementGetByPos(regs, j);
	 reg_j = preg->reg;
	 
	 shRegIntLincom(basis->regs[i-1][0][0]->reg, reg_j, 0, 1, Q->me[j][i]);
      }
#define SUBTRACT_SKY 0
#if SUBTRACT_SKY
/*
 * Estimate and subtract the mean background level
 *
 * It is not at all clear that doing this is a good idea; it'd be
 * better to get the sky level right in the first place.
 */
      if(i == 0) {			/* only zeroth KL component */
	 float sky = 0;			/* estimate of sky level */
	 REGION *sreg;			/* reg_i minus the border */

	 reg_i = basis->regs[i-1][0][0]->reg;
	 shAssert(reg_i->type == TYPE_FL32);
	 
	 for(j = border + 1; j < nrow - border - 1; j++) {
	    sky += reg_i->rows_fl32[j][border] +
					reg_i->rows_fl32[j][ncol - border - 1];
	 }
	 for(k = border; k < ncol - border; k++) {
	    sky += reg_i->rows_fl32[border][k] +
					reg_i->rows_fl32[nrow - border - 1][k];
	 }
	 sky /= 2*((nrow - 2*border) + (ncol - 2*border)) - 4;
	 
	 sreg = shSubRegNew("", basis->regs[i-1][0][0]->reg,
			    nrow - 2*border, ncol - 2*border,
			    border, border, NO_FLAGS);
	 shAssert(sreg != NULL);
	 
	 shRegIntConstAdd(sreg, -sky, 0);
	 shRegDel(sreg);
      }
#endif

      {					/* find eigenimage's normalisation */
	 FL32 *row;			/* a row of the image */
	 int row_0, col_0;		/* LLC and URC of part of regions */
	 int row_1, col_1;		/*         use from PSF_BASIS->regs */
	 double sum = 0;

	 row_0 = border; col_0 = border;
	 row_1 = nrow - border - 1; col_1 = ncol - border - 1;

         /* N.B. PSF_BASIS->counts is unweighted sum of all pixel values */
	 /* use the same loop to set pixels in the "border" regions to 0 */
	 for(j = 0; j < nrow; j++) {
	    row = basis->regs[i-1][0][0]->reg->rows_fl32[j];
	    for(k = 0; k < ncol; k++) {
	       if (j >= row_0 && j <= row_1 && k >= col_0 && k <= col_1) {
	           sum += row[k];
	       } else {
                   row[k] = 0;
               }
	    }
	 }
	 basis->regs[i-1][0][0]->counts = sum;
      }
   }
/*
 * clean up
 */
   phVecDel(lambda);
   phMatDel(Q);
   phMatDel(R);

   return(basis);
}
#endif

/*****************************************************************************/
/*
 * Given the position of an object (<rowc>, <colc>) and a KL PSF_BASIS
 * with its coefficients set, return a PSF_REG consisting of the PSF
 * reconstructed at that point
 */
PSF_REG *
phPsfKLReconstruct(const PSF_BASIS *basis, /* basis to use */
		   float rowc,		/* location of */
		   float colc,		/*    desired PSF */
		   PIXDATATYPE regType)	/* type of desired REGION */

{
   float acoeff;			/* the "a" coefficient for component */
   float *coeffs;			/* coefficients for spatial variation*/
   const REGION *comp;			/* component of KL decomposition */
   double counts;			/* psf_out->counts (from basis) */
   const PSF_KERNEL *kern;		/* == basis->kern */
   int is;				/* counter in sigma (== comp) dirn. */
   int i, j, k;
   int ncomp;				/* number of KL components in use */
   int nrow, ncol;			/* == reg->n{rol,col} */
   int nrow_b, ncol_b;			/* == kern->n{row,col}_b */
   PSF_REG *psf_out;			/* returned PSF_REG */
   REGION *out;				/* == psf_out->reg */
   FL32 **rows_comp, *row_comp;		/* == comp->rows, comp->rows[] */
   FL32 **rows_out, *row_out;		/* == out->rows, out->rows[] */

   kern = basis->kern;
   shAssert(kern != NULL);   
   shAssert(kern->nrow_a == 1 && kern->ncol_a == 1);
   ncomp = kern->nsigma + 1;
   nrow_b = kern->nrow_b;
   ncol_b = kern->ncol_b;

   shAssert(basis != NULL && basis->regs != NULL && basis->type == KL_BASIS);
   comp = basis->regs[-1][0][0]->reg;
   shAssert(comp != NULL && comp->type == TYPE_FL32);
   nrow = comp->nrow; ncol = comp->ncol;
   psf_out = phPsfRegNew(nrow, ncol, TYPE_FL32);
   out = psf_out->reg;
   rows_out = out->rows_fl32;
/*
 * Set a convenience array coeffs with the powers of rowc and colc to be used
 * for each of the nrow_b*ncol_b spatial terms
 */
   coeffs = alloca(nrow_b*ncol_b*sizeof(float));

   for(k = 0; k < nrow_b*ncol_b; k++) {
      coeffs[k] = pow(rowc*RC_SCALE, k%nrow_b)*pow(colc*RC_SCALE, k/nrow_b);
   }

   shRegClear(out);
   counts = 0;
   for(is = -1; is < ncomp - 1; is++) {
      acoeff = 0;
      for(k = 0; k < nrow_b*ncol_b; k++) {
	 acoeff += kern->a[is][0][0].c[k%nrow_b][k/nrow_b]*coeffs[k];
      }
      
      counts += acoeff*basis->regs[is][0][0]->counts;
      comp = basis->regs[is][0][0]->reg;
      rows_comp = comp->rows_fl32;
	
      for(i = 0; i < nrow; i++) {
	 row_out = rows_out[i];
	 row_comp = rows_comp[i];
	 for(j = 0; j < ncol; j++) {
	    row_out[j] += acoeff*row_comp[j];
	 }
      }
   }
/*
 * Sky subtract?
 */
#define SKY_SUBTRACT_PSF 0
#if SKY_SUBTRACT_PSF
   {
      const int border = basis->kern->border;
      int nel = 0;
      int nval = 2*(nrow + ncol - 4*border) - 4;
      float *arr = alloca(nval*sizeof(float));
      float sky = 0;

      for(i = border + 1; i < nrow - border - 1; i++) {
	 arr[nel++] = rows_out[border][i];
	 arr[nel++] = rows_out[nrow - border - 1][i];
      }
      for(i = border; i < ncol - border; i++) {
	 arr[nel++] = rows_out[i][border];
	 arr[nel++] = rows_out[i][ncol - border - 1];
      }
      shAssert(nel == nval);

      phFloatArrayStats(arr, nel, 0, &sky, NULL, NULL, NULL);

      for(i = border; i < nrow - border; i++) {
	 row_out = rows_out[i];
	 for(j = border; j < ncol - border; j++) {
	    row_out[j] -= sky;
	 }
      }
   }
#endif
/*
 * convert to desired type, if not the same as out
 */
   if(regType != TYPE_FL32) {
      REGION *iout = shRegNew("", nrow, ncol, regType);
      if(regType == TYPE_U16) {
	 U16 **rows_iout = iout->rows_u16, *row_iout;
	 int pixval;
         float scale = 0.0;
  	 if (out->rows_fl32[nrow/2][ncol/2] > 0) {
	    scale = (3e4 - SOFT_BIAS)/out->rows_fl32[nrow/2][ncol/2];
         }
	 for(i = 0; i < nrow; i++) {
	    row_out = rows_out[i];
	    row_iout = rows_iout[i];
	    for(j = 0; j < ncol; j++) {
	       pixval = scale*row_out[j] + SOFT_BIAS + 0.5;
	       row_iout[j] = (pixval < 0) ? 0 : pixval;
	    }
	 }
	 counts *= scale;
      } else {
	 shRegIntCopy(iout, out);
      }
      shRegDel(out);
      psf_out->reg = iout;
   }

   psf_out->counts = counts;
   psf_out->rowc = rowc;
   psf_out->colc = colc;

   return(psf_out);
}

/*****************************************************************************/
/*
 * Make/destroy PSF_KL_COMPs
 */
PSF_KL_COMP *
phPsfKLCompNew(REGION *reg)
{
   int i, j;
   PSF_KL_COMP *klc = shMalloc(sizeof(PSF_KL_COMP));

   klc->nrow_b = klc->ncol_b = 0;
   for(i = 0; i < MAX_ORDER_B; i++) {
      for(j = 0; j < MAX_ORDER_B; j++) {
	 klc->c[i][j] = 0;
      }
   }
   klc->lambda = 0.0;
   klc->reg = reg;

   return(klc);   
}

void
phPsfKLCompDel(PSF_KL_COMP *klc)
{
   if(klc != NULL) {
      shRegDel(klc->reg);
   }
   shFree(klc);
}

#if !defined(STAND_ALONE)
/*****************************************************************************/
/*
 * Set the fields in a PSF_KL_COMP from a PSF_BASIS
 */
PSF_KL_COMP *
phPsfKLCompSetFromBasis(PSF_KL_COMP *klc, /* PSF_KL_COMP to set, or NULL */
			const PSF_BASIS *basis,	/* basis with desired data */
			int comp)	/* desired component */
{
   int i, j;
   
   shAssert(basis != NULL && basis->kern != NULL);
   shAssert(basis->type == KL_BASIS && basis->regs != NULL);
   shAssert(comp >= 0 && comp <= basis->kern->nsigma);
   shAssert(basis->kern->a[comp - 1] != NULL);
   shAssert(basis->regs[comp - 1] != NULL);

   if(klc == NULL) {
      klc = phPsfKLCompNew(NULL);
   }
   
   klc->nrow_b = basis->kern->nrow_b;
   klc->ncol_b = basis->kern->ncol_b;
   for(i = 0; i < MAX_ORDER_B; i++) {
      for(j = 0; j < MAX_ORDER_B; j++) {
	 klc->c[i][j] = basis->kern->a[comp - 1][0][0].c[i][j];
      }
   }

   klc->lambda = basis->kern->sigma[comp - 1];
   klc->reg = basis->regs[comp - 1][0][0]->reg;
   klc->counts = basis->regs[comp - 1][0][0]->counts;

   return(klc);
}
#endif
/*
 * Set some fields in a PSF_BASIS from a PSF_KL_COMP
 */
void
phPsfBasisSetFromKLComp(PSF_BASIS *basis, /* basis to be set */
			int comp,	/* component to set*/
			const PSF_KL_COMP *klc, /* PSF_KL_COMP with data */
			int copy_reg)	/* copy region? */
{
   int i, j;
   
   shAssert(basis != NULL && basis->kern != NULL);
   shAssert(basis->type == KL_BASIS && basis->regs != NULL);
   shAssert(comp >= 0 && comp <= basis->kern->nsigma);
   shAssert(basis->kern->a[comp - 1] != NULL);
   shAssert(basis->regs[comp - 1] != NULL);
   shAssert(klc != NULL);
   
   shAssert(klc->nrow_b == basis->kern->nrow_b);
   shAssert(klc->ncol_b == basis->kern->ncol_b);
   
   for(i = 0; i < MAX_ORDER_B; i++) {
      for(j = 0; j < MAX_ORDER_B; j++) {
	 basis->kern->a[comp - 1][0][0].c[i][j] = klc->c[i][j];
      }
   }

   basis->kern->sigma[comp - 1] = klc->lambda;
   basis->regs[comp - 1][0][0]->counts = klc->counts;
   if(copy_reg) {
      REGION *reg = basis->regs[comp - 1][0][0]->reg;
      if(klc->reg == NULL) {
	 reg = basis->regs[comp - 1][0][0]->reg =
				 shRegNew("", reg->nrow, reg->ncol, reg->type);
      }
      shRegIntCopy(reg, klc->reg);
   } else {
      shAssert(basis->regs[comp - 1][0][0]->reg == NULL);
      basis->regs[comp - 1][0][0]->reg = klc->reg;
   }
}

#if !defined(STAND_ALONE)
/*****************************************************************************/
/*
 * Code to calculate derivatives of critically sampled images,
 * using derivatives of sinc functions
 */
/* 
 * Generate a cos-belled sinc filter to differentiate an image a
 * specified number of times; the filter is tapered with 
 * a cosine bell of half-length lbell (filter length 2*lbell - 1). 
 *
 * The filter is returned; either the one you pass in, or one
 * allocated for you if you pass a NULL pointer
 */
static float *
get_sinc_derivative(float *filt,	/* the filter in question,
					   dimen >= 2*lbell - 1, or NULL */
		    int lbell,		/* half-length of cosbell */
		    int nderiv)		/* order of derivative desired */
{
   float cbell;				/* value of cosbell */
   int cen = lbell - 1;			/* index of centre of filter */
   int i;
   const int len = 2*lbell - 1;		/* total length of filter */
   double sum;				/* used to normalise */
   int sgn;				/* signs alternate */
   const double flbell = M_PI/(double)lbell;

   if(filt == NULL) {
      filt = shMalloc(len*sizeof(float));
   }
/*
 * generate the coefficients
 */
    switch (nderiv) {
     case 0:				/* identity */
       filt[cen] = 1;
       for(i = 1; i < lbell; i++) {
	  filt[cen - i] = filt[cen + i] = 0;
       }

       return(filt);
       break;				/* NOTREACHED */
     case 1:				/* first derivative */
       filt[cen] = 0;
       for(sgn = -1, i = 1; i < lbell; sgn = -sgn, i++) {
	  filt[cen + i] = sgn/(float)i;
	  filt[cen - i] = -filt[cen + i];
       }

       break;
     case 2:				/* second derivative */
       filt[cen] = -M_PI*M_PI/3;	/* value for infinite series */
       for(sgn = +1, i = 1; i < lbell; sgn = -sgn, i++) {
	  filt[cen + i] = 2*sgn/(float)(i*i);
	  filt[cen - i] = filt[cen + i];
       }

       break;
     default:
       shFatal("get_sinc_derivative: I can't differentiate %d times", nderiv);
       break;				/* NOTREACHED */
    }
/*
 * Apply cosbell, and force filter to sum to zero (the derivative of
 * a constant is 0).  The case nderiv == 0 is already covered
 */
   sum = 0;
   for(i = 1; i < lbell; i++) {
      cbell = 0.5*(cos(i*flbell) + 1.0);
      filt[cen - i] *= cbell;
      filt[cen + i] *= cbell;
      sum += filt[cen - i] + filt[cen + i];
   }
   filt[cen] = -sum;
/*
 * Calculate appropriate normalisation. This is a little tricky, as the
 * coefficients have zero sum. We get around this by normalising against
 * functions with constant first/second derivatives
 */
   if(nderiv == 1) {
      sum = 0*filt[cen];
      for(i = 1; i < lbell; i++) {
	 sum += 2*i*filt[cen + i];	/* == -i*filt[cen-i] + i*filt[cen+i] */
      }
   } else {
      shAssert(nderiv == 2);

      sum = 0*filt[cen];
      for(i = 1; i < lbell; i++) {
	 sum += i*i*filt[cen + i];	/* == i^2/2*(filt[cen-i]+filt[cen+i])*/
      }
   }
/*
 * Apply that normalisation
 */
   shAssert(fabs(sum) > 1e-10);
   {
      const float isum = 1/sum;
      
      filt[cen] *= isum;
      for(i = 1; i < lbell; i++) {
	 filt[cen + i] *= isum;
	 filt[cen - i] *= isum;
      }
   }

   return(filt);
}

/*****************************************************************************/
/*
 * Calculate a PSF_BASIS by (sinc) differentiating a PSF_BASIS, returning
 * the desired basis
 */
PSF_BASIS *
phPsfBasisDifferentiate(const PSF_BASIS *basis,	/* input basis */
			REGION *scr,	/* scratch region, or NULL */
			int nderiv_row,	/* order of derivative in row */
			int nderiv_col,	/*             and column directions */
			int lbell)	/* desired half-length of cosbell */
{
   int free_scr = 0;			/* should we free REGION scr? */
   PSF_BASIS *deriv;			/* desired derivatives */
   float *filt_r;			/* filter in row direction */
   float *filt_c;			/* filter in column direction */
   int is;				/* counter for sigma */
   int ira, ica;			/* counter for row/column coeffs */
   int lbell_r, lbell_c;		/* adopted half-length in row, column*/
   int nrow, ncol;			/* == basis->regs[][][]->reg->nrow */
   PIXDATATYPE type;			/* == basis->regs[][][]->reg->type */

   shAssert(basis != NULL && basis->kern != NULL && basis->kern->nsigma > 0);
   shAssert(basis->regs != NULL);
   nrow = basis->regs[0][0][0]->reg->nrow;
   ncol = basis->regs[0][0][0]->reg->ncol;
   type = basis->regs[0][0][0]->reg->type;
   shAssert(type == TYPE_FL32);

   if(scr == NULL) {
      scr = shRegNew("", nrow, ncol, type);
      free_scr = 1;			/* we allocated it, so free it */
   } else {
      shAssert(scr->nrow == nrow && scr->ncol == ncol && scr->type == type);
   }
/*
 * allocate deriv basis, and setup the desired filters
 */
   deriv = phPsfBasisNew(basis->kern, basis->type,
			 basis->regs[0][0][0]->reg->nrow,
			 basis->regs[0][0][0]->reg->ncol,
			 basis->regs[0][0][0]->reg->type);

   lbell_r = (nderiv_row == 0) ? 1 : lbell;
   lbell_c = (nderiv_col == 0) ? 1 : lbell;
   filt_r = get_sinc_derivative(NULL, lbell_r, nderiv_row);
   filt_c = get_sinc_derivative(NULL, lbell_c, nderiv_col);
/*
 * do the work
 */
   phConvolve(deriv->regs[-1][0][0]->reg, basis->regs[-1][0][0]->reg, scr,
	      2*lbell_c - 1, 2*lbell_r - 1, filt_c, filt_r, 0, CONVOLVE_MULT);
   for(is = 0; is < basis->kern->nsigma; is++) {
      for(ira = 0; ira < basis->kern->nrow_a; ira++) {
	 for(ica = 0; ica < basis->kern->ncol_a; ica++) {
	    phConvolve(deriv->regs[is][ira][ica]->reg,
		       basis->regs[-1][0][0]->reg, scr,
		       2*lbell_c - 1, 2*lbell_r - 1, filt_c, filt_r,
		       0, CONVOLVE_MULT);
	 }
      }
   }
/*
 * clean up
 */
   shFree(filt_r); shFree(filt_c);
   if(free_scr) {
      shRegDel(scr);
   }
   return(deriv);
}

/*****************************************************************************/
void
phPsfBasisSetDerivatives(PSF_BASIS *basis, /* basis to differentiate */
			 REGION *scr,	/* scratch region, or NULL */
			 int lbell)	/* desired half-length of cosbell */
{
   int free_scr = 0;			/* should we free REGION scr? */
   int i;
   int nrow, ncol;			/* == basis->regs[][][]->reg->nrow */
   PIXDATATYPE type;			/* == basis->regs[][][]->reg->type */

   shAssert(basis != NULL && basis->kern != NULL);
   shAssert(basis != NULL && basis->kern != NULL && basis->kern->nsigma > 0);
   shAssert(basis->regs != NULL);
   nrow = basis->regs[0][0][0]->reg->nrow;
   ncol = basis->regs[0][0][0]->reg->ncol;
   type = basis->regs[0][0][0]->reg->type;
   shAssert(type == TYPE_FL32);

   if(scr == NULL) {
      scr = shRegNew("", nrow, ncol, type);
      free_scr = 1;			/* we allocated it, so free it */
   } else {
      shAssert(scr->nrow == nrow && scr->ncol == ncol && scr->type == type);
   }
/*
 * do the work
 */
   for(i = 0; i < 5; i++) {
      phPsfBasisDel(basis->deriv[i]);
   }

   basis->deriv[DBASIS_DROW] =
			      phPsfBasisDifferentiate(basis, scr, 1, 0, lbell);
   basis->deriv[DBASIS_DCOL] =
			      phPsfBasisDifferentiate(basis, scr, 0, 1, lbell);
   basis->deriv[DBASIS_DROW2] =
			      phPsfBasisDifferentiate(basis, scr, 2, 0, lbell);
   basis->deriv[DBASIS_DROWDCOL] =
			      phPsfBasisDifferentiate(basis, scr, 1, 1, lbell);
   basis->deriv[DBASIS_DCOL2] =
			      phPsfBasisDifferentiate(basis, scr, 0, 2, lbell);
/*
 * clean up
 */
   if(free_scr) {
      shRegDel(scr);
   }   
}
#endif

#if defined(STAND_ALONE)

/***************************************************************************
 * <AUTO EXTRACT>
 *
 * Copy one region to another. It is also able to convert U16 to FL32,
 * removing the SOFT_BIAS
 *
 * return: SH_SUCCESS		If region type is supported
 *         SH_GENERIC_ERROR	otherwise
 */
static int
shRegIntCopy(REGION *out,		/* output region */
	     const REGION *in)		/* input region */
{
   int i, j;
   int ncol,nrow;			/* unpacked from out */

   shAssert(out != NULL && in != NULL);
   shAssert(in->nrow == out->nrow && in->ncol == out->ncol);

   ncol = out->ncol;
   nrow = out->nrow;

   if(out->type == TYPE_U8) {
      shAssert(in->type == out->type);
      for(i = 0;i < nrow;i++) {
	 memcpy(out->rows_u8[i],in->rows_u8[i],ncol);
      }
   } else if(out->type == TYPE_S8) {
      shAssert(in->type == out->type);
      for(i = 0;i < nrow;i++) {
	 memcpy(out->rows_s8[i],in->rows_s8[i],ncol);
      }
   } else if(out->type == TYPE_U16) {
      if(in->type == TYPE_U16) {
	 for(i = 0;i < nrow;i++) {
	    memcpy(out->rows_u16[i],in->rows_u16[i],ncol*sizeof(U16));
	 }
      } else if(in->type == TYPE_FL32) {
	 FL32 *iptr;
	 U16 *optr;

	 for(i = 0;i < nrow;i++) {
	    iptr = in->rows_fl32[i];
	    optr = out->rows_u16[i];
	    for(j = 0;j < ncol; j++) {
	       optr[j] = iptr[j] + SOFT_BIAS + 0.5;
	    }
	 }
      } else {
	 shError("shRegIntCopy doesn't convert regions of type %d to U16\n",
		 in->type);
	 return(SH_GENERIC_ERROR);
      }
   } else if(out->type == TYPE_S16) {
      shAssert(in->type == out->type);
      for(i = 0;i < nrow;i++) {
	 memcpy(out->rows_s16[i],in->rows_s16[i],ncol*sizeof(S16));
      }
   } else if(out->type == TYPE_U32) {
      shAssert(in->type == out->type);
      for(i = 0;i < nrow;i++) {
	 memcpy(out->rows_u32[i],in->rows_u32[i],ncol*sizeof(U32));
      }
   } else if(out->type == TYPE_S32) {
      shAssert(in->type == out->type);
      for(i = 0;i < nrow;i++) {
	 memcpy(out->rows_s32[i],in->rows_s32[i],ncol*sizeof(S32));
      }
   } else if(out->type == TYPE_FL32) {
      if(in->type == TYPE_U16) {
	 U16 *iptr;
	 FL32 *optr;

	 for(i = 0;i < nrow;i++) {
	    iptr = in->rows_u16[i];
	    optr = out->rows_fl32[i];
	    for(j = 0;j < ncol; j++) {
	       optr[j] = (int)iptr[j] - SOFT_BIAS;
	    }
	 }
      } else if(in->type == TYPE_FL32) {
	 for(i = 0;i < nrow;i++) {
	    memcpy(out->rows_fl32[i],in->rows_fl32[i],ncol*sizeof(FL32));
	 }
      } else {
	 shError("shRegIntCopy doesn't convert regions of type %d to FL32\n",
		 in->type);
	 return(SH_GENERIC_ERROR);
      }
   } else {
      shError("shRegIntCopy doesn't handle regions of type %d\n",
	      out->type);
      return(SH_GENERIC_ERROR);
   }

   return(SH_SUCCESS);
}
#endif

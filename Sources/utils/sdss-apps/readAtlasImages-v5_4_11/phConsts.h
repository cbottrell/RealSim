#if !defined(PHCONSTS_H)
#define PHCONSTS_H
/*
 * The fundamental REGION datatype; all values must be consistent!
 *
 * N.b. Use #define not typedef for PIX so as not to confuse dervish's schema
 */
#if !defined(FLOATING_PHOTO)
#  define FLOATING_PHOTO 0
#endif

#define PHOTO_U16_COMPAT 0		/* try to make floating photo give
					   exactly the same results as int? */

#if FLOATING_PHOTO
#  error Floating point photo is not supported in this release
#  define TYPE_PIX TYPE_FL32
#  define PIX FL32			/* pixels' datatype */
#  define ROWS rows_fl32		/* Name of rows pointer in REGION */
#  define FLT2PIX(V) (V)		/* how to round a pixel value */
#  define PIX2INT(V) ((int)(V + 0.5))	/* convert a pixel to integer */
#  define PIXFMT "%g"			/* format to print a pixel */
#else
#  define TYPE_PIX TYPE_U16
#  define PIX U16			/* pixel's datatype */
#  define ROWS rows_u16			/* Name of rows pointer in REGION */
#  define FLT2PIX(V) (V + 0.5)		/* how to round a pixel value */
#  define PIX2INT(V) (V)		/* convert a pixel to integer */
#  define PIXFMT "%d"			/* format to print a pixel */
#endif

/*
 * Misc constants that crop up all over photo.  N.b. PHOTO_CONSTS() at
 * the TCL level contains ints, so don't try putting floats into it
 */
#define NCOLOR 5			/* max number of colours */

#define SOFT_BIAS 1000			/* soft bias for regions */

#define IQR_TO_SIGMA 0.741301		/* sigma = 0.741*(interquartile range)
					   for a Gaussian */

#define MIN_2N_BIAS -0.5641895835	/* The mean value of the minimum of
					   two N(0,1) variates is MIN_N2_BIAS*/

#define MIN_S8 -128			/* smallest possible S8
					   pragma typedef { */
#define MAX_S8 127			/* largest possible S8 */
#define MAX_U8 255			/* largest possible U8 */
#define MIN_S16 -32768			/* smallest possible S16 */
#define MAX_S16 32767			/* largest possible S16 */
#define MAX_U16 65535			/* largest possible U16 */
#define MIN_S32 (-MAX_S32-1)		/* smallest possible S32, as an int */
#define MAX_S32 2147483647		/* largest possible S32 */
#define MAX_U32 4294967295UL		/* largest possible U32 */
#define S32_SIGN_BIT (1U<<31)		/* sign bit for 32 bit ints */
#define ERROR_IS_BAD -1000		/* Value of error if error is unknown*/
#define VALUE_IS_BAD -9999		/* Value of error if value is unknown
					   pragma } PHOTO_CONSTS */

/*****************************************************************************/
/*
 * Softening for likelihoods
 */
#define L_SOFT 1e-10
/*
 * LPC coefficients for sigma = 1, S/N = infty
 */
#define INTERP_1_C1 0.7737
#define INTERP_1_C2 -0.2737
/*
 * LPC coefficients for sigma = 1/sqrt(2), S/N = infty. These are the coeffs
 * to use when interpolating at 45degrees to the row/column
 */
#define INTERP_IS2_C1 0.7358
#define INTERP_IS2_C2 -0.2358



#endif

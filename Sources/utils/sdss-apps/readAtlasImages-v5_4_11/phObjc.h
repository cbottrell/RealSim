#if !defined(PHOBJC_H)
#define PHOBJC_H
#include "dervish.h"
#include "phExtract.h"
#include "phSpanUtil.h"
#include "phFramestat.h"
#include "phPeaks.h"
#include "phConsts.h"

/*
 * An enum for types of object
 */
typedef enum {
   OBJ_UNK = 0, 
   OBJ_CR, 
   OBJ_DEFECT, 
   OBJ_GALAXY, 
   OBJ_GHOST, 
   OBJ_KNOWNOBJ,
   OBJ_STAR, 
   OBJ_TRAIL,
   OBJ_SKY,
   OBJ_NTYPE
} OBJ_TYPE;

/* here are the values the "flags" field can take, in bitwise fashion */
#define OBJECT1_CANONICAL_CENTER 0x1	/* used canonical, not local, centre
				   pragma typedef { */
#define OBJECT1_BRIGHT 0x2	/* detected by Bright Objects */
#define OBJECT1_EDGE 0x4	/* object is too close to edge of frame */
#define OBJECT1_BLENDED 0x8	/* object is/was blended */
#define OBJECT1_CHILD  0x10	/* object is a child */
#define OBJECT1_PEAKCENTER 0x20	/* given centre is position of peak pixel */
#define OBJECT1_NODEBLEND 0x40	/* no deblending attempted */
#define OBJECT1_NOPROFILE 0x80	/* too small to estimate a profile */
#define OBJECT1_NOPETRO 0x100	/* no Petrosian radius */
#define OBJECT1_MANYPETRO 0x200	/* more than one Petrosian radius */
#define OBJECT1_NOPETRO_BIG 0x400 /* no Petrosian radius as object is too big*/
#define OBJECT1_DEBLEND_TOO_MANY_PEAKS 0x800 /* too many peaks to deblend */
#define OBJECT1_CR 0x1000	/* contains a CR pixel */
#define OBJECT1_MANYR50 0x2000	/* more than one 50% radius */
#define OBJECT1_MANYR90 0x4000	/* more than one 90% radius */
#define OBJECT1_BAD_RADIAL 0x8000 /* some low S/N radial points */
#define OBJECT1_INCOMPLETE_PROFILE 0x10000 /* r_P includes off-frame pixels */
#define OBJECT1_INTERP 0x20000	/* object contains interpolated pixels */
#define OBJECT1_SATUR 0x40000	/* object contains saturated pixels */
#define OBJECT1_NOTCHECKED 0x80000 /* object contains NOTCHECKED pixels */
#define OBJECT1_SUBTRACTED 0x100000 /* object had wings subtracted */
#define OBJECT1_NOSTOKES 0x200000 /* object has no measured stokes params */
#define OBJECT1_BADSKY 0x400000  /* sky level is so bad that object is -ve */
#define OBJECT1_PETROFAINT 0x800000	/* >= 1 Petrosian radius too faint */
#define OBJECT1_TOO_LARGE 0x1000000	/* object is too large */
#define OBJECT1_DEBLENDED_AS_PSF 0x2000000 /* deblender treated obj as PSF */
#define OBJECT1_DEBLEND_PRUNED 0x4000000 /* deblender pruned peak list */
#define OBJECT1_ELLIPFAINT 0x8000000	/* Centre's fainter than desired
					   elliptical isophote */
#define OBJECT1_BINNED1 0x10000000 /* object was found in 1x1 binned image */
#define OBJECT1_BINNED2 0x20000000 /* object was found in 2x2 binned image */
#define OBJECT1_BINNED4 0x40000000 /* object was found in 4x4 binned image */
#define OBJECT1_MOVED 0x80000000 /* Object may have moved; largest flag value
				     pragma } OBJECT1_FLAGS */

#define OBJECT1_DETECTED (OBJECT1_BINNED1|OBJECT1_BINNED2|OBJECT1_BINNED4)
/*
 * here are values for flags2
 */
#define OBJECT2_DEBLENDED_AS_MOVING 0x1	/* deblended as a moving object
					   pragma typedef { */
#define OBJECT2_NODEBLEND_MOVING 0x2	/* no deblend of moving object */
#define OBJECT2_TOO_FEW_DETECTIONS 0x4  /* too few detections to deblend
					   as moving */
#define OBJECT2_BAD_MOVING_FIT 0x8	/* fit to moving object was too poor */
#define OBJECT2_STATIONARY 0x10		/* velocity is consistent with zero */
#define OBJECT2_PEAKS_TOO_CLOSE 0x20	/* at least some peaks were too
					   close, and thus merged */
#define OBJECT2_BINNED_CENTER 0x40	/* image was binned while centroiding*/
#define OBJECT2_LOCAL_EDGE 0x80		/* per-band centre's too near edge */
#define OBJECT2_BAD_COUNTS_ERROR 0x100	/* psf|fiberCountsErr is bad/unknown */
#define OBJECT2_BAD_MOVING_FIT_CHILD 0x200 /* moving child's fit was too poor*/
#define OBJECT2_DEBLEND_UNASSIGNED_FLUX 0x400 /* deblender failed to assign
						 enough of flux to children */
#define OBJECT2_SATUR_CENTER 0x800	/* object's centre's saturated */
#define OBJECT2_INTERP_CENTER 0x1000	/* object's centre's interpolated */
#define OBJECT2_DEBLENDED_AT_EDGE 0x2000 /* object's deblended although EDGE */
#define OBJECT2_DEBLEND_NOPEAK 0x4000	/* object had no detected peak */
#define OBJECT2_PSF_FLUX_INTERP 0x8000	/* a significant amount of
					   PSF's flux is interpolated */
#define OBJECT2_TOO_FEW_GOOD_DETECTIONS 0x10000	/* too few good detections to
						   deblend as moving */
#define OBJECT2_CENTER_OFF_AIMAGE 0x20000 /* at least one peak's centre lay off
					     the atlas image in some band */
#define OBJECT2_DEBLEND_DEGENERATE 0x40000 /* at least one potential child has
					      been pruned as being too similar
					      to some other template */
#define OBJECT2_BRIGHTEST_GALAXY_CHILD 0x80000 /* this is the brightest child
						  galaxy in a blend */
#define OBJECT2_CANONICAL_BAND 0x100000 /* This band was primary (usually r')*/
#define OBJECT2_AMOMENT_UNWEIGHTED 0x200000 /* `adaptive' moments are actually
					       unweighted */
#define OBJECT2_AMOMENT_SHIFT 0x400000	/* centre moved too far while
					   determining adaptive moments */
#define OBJECT2_AMOMENT_MAXITER 0x800000 /* Too many iterations while
					    determining adaptive moments */
#define OBJECT2_MAYBE_CR 0x1000000	/* object may be a cosmic ray */
#define OBJECT2_MAYBE_EGHOST 0x2000000	/* object may be an electronics ghost*/
#define OBJECT2_NOTCHECKED_CENTER 0x4000000 /* object's centre is NOTCHECKED */
#define OBJECT2_HAS_SATUR_DN 0x8000000	/* Counts include DN in bleed trails */
#define OBJECT2_SPARE4 0x10000000	/* unused */
#define OBJECT2_SPARE3 0x20000000	/* unused */
#define OBJECT2_SPARE2 0x40000000	/* unused */
#define OBJECT2_SPARE1 0x80000000	/* unused 
					   pragma } OBJECT2_FLAGS */
/*
 * These are book-keeping bits, and aren't written to disk
 */
#define OBJECT3_HAS_SATUR_DN 0x1	/* object has extra (saturated) DN
					   pragma typedef { */
#define OBJECT3_MEASURED 0x10000000	/* object has been measured */
#define OBJECT3_GROWN_MERGED 0x20000000	/* growing led to a merger */
#define OBJECT3_HAS_CENTER 0x40000000	/* OBJC has a canonical centre */
#define OBJECT3_MEASURE_BRIGHT 0x80000000 /* object should be measured bright
					     pragma } OBJECT3_FLAGS */
/*
 * This struct holds one object in one colour
 *
 * In the following comments, (NS) means Not Saved to disk at the end of a run
 */
typedef struct {
   const int id;		/* id number (set in constructor) */
   REGION *region;		/* REGION containing object */
   OBJMASK *mask;		/* MASK of which pixels are in this object */
   float rowc,rowcErr;		/* row position and error of center */
   float colc,colcErr;		/* column position and error of center */
   float sky, skyErr;		/* sky value and error near this object */
/*
 * PSF, aperture, and Petrosian fits/magnitudes
 */
   float psfCounts,psfCountsErr;/* Counts via PSF-fitting and error */
   float fiberCounts,fiberCountsErr;/* Counts within 3" aperture */
   float petroCounts, petroCountsErr;/* Counts defining a Petrosian magnitude*/
   float petroRad, petroRadErr;	/* Petrosian radius (pixels) */
   float petroR50, petroR50Err;	/* radius with 50% of petrosian light */
   float petroR90, petroR90Err;	/*  "  "    "  90% "    "    "   " "  */
/*
 * Shape of object
 */
   float Q, QErr, U, UErr;		/* <(col^2 - row^2)/radius^2> etc. */

   float M_e1, M_e2;			/* Adaptive E1/E2 shape measures */
   float M_e1e1Err, M_e1e2Err, M_e2e2Err; /* E1/E2 covariances */
   float M_rr_cc, M_rr_ccErr;		/* Adaptive (<r^2> + <c^2>) and error*/
   float M_cr4;				/* Adaptive fourth moment of object */
   float M_e1_psf, M_e2_psf;		/* Adaptive E1/E2 for PSF*/
   float M_rr_cc_psf;			/* Adaptive (<r^2> + <c^2>) for PSF */
   float M_cr4_psf;			/* Adaptive fourth moment of PSF */
/*
 * Properties of a certain isophote. The `Err' is the 1-sigma error,
 * the `Grad' how much the value changes if we go to an isophote
 * one magnitude fainter
 */
   float iso_rowc, iso_rowcErr, iso_rowcGrad;
   float iso_colc, iso_colcErr, iso_colcGrad;
   float iso_a, iso_aErr, iso_aGrad;
   float iso_b, iso_bErr, iso_bGrad;
   float iso_phi, iso_phiErr, iso_phiGrad;
/*
 * Model parameters for deV and exponential models
 */
   float r_deV, r_deVErr;
   float ab_deV, ab_deVErr;
   float phi_deV, phi_deVErr;
   float counts_deV, counts_deVErr;
   float r_exp, r_expErr;
   float ab_exp, ab_expErr;
   float phi_exp, phi_expErr;
   float counts_exp, counts_expErr;
/*
 * The counts for the galaxy model that fit best in the canonical band
 */
   float counts_model, counts_modelErr;
/*
 * Measures of image structure
 */
   float texture;
/*
 * Classification information
 */
   float star_L, star_lnL;	/* (log-)likelihood of being a star */
   float exp_L, exp_lnL;	/* (log-)likelihood of being an exp. disk */
   float deV_L, deV_lnL;	/* (log-)likelihood of being a deV profile */

   float fracPSF;		/* fraction of total light that's in a
				   stellar component */

   int flags;			/* flags for e.g. blends */
   int flags2;			/* More flag bits */

   OBJ_TYPE type;		/* classification: OBJ_STAR, etc. */
   float prob_psf;		/* Bayesian probability of being a PSF */
/*
 * Profile and extent of object
 */
   int nprof;			/* num of valid profile values in arrays */
   float profMean[NANN];	/* Radial profile, mean within annuli */
   float profMed[NANN];		/* Radial profile, median within annuli */
   float profErr[NANN];		/* Radial profile, stdev from mean in annuli */
/*
 * None of the quantities below this comment are to be written to disk
 */
   int flags3;			/* flags for internal photo info (NS) */
   float aratio;		/* axis ratio from central moments (NS) */
   float majaxis;		/* pos. angle of major axis (radians) (NS) */
   int comp;			/* flag for compression of REGION/MASK (NS)*/
   PEAKS *peaks;		/* all peaks in object (NS) */
   int npix;			/* number of pixels in object (NS) */
   float satur_DN;		/* DN in bleed trails associated with object */
/*
 * These are the quantities from which the likelihoods are derived
 */
   int nu_star;				/* number of d.o.f. for psf fit (NS)*/
   float chisq_star;			/* chi^2 for psf fit (NS) */
   int nu_deV;				/* number of d.o.f. for deV fit (NS) */
   float chisq_deV;			/* chi^2 for deV fit (NS) */
   int nu_exp;				/* number of d.o.f. for disk fit (NS)*/
   float chisq_exp;			/* chi^2 for disk fit (NS) */
} OBJECT1;

/*
 * Here's an atlas image
 */
typedef struct {
   int id;				/* id number for this objc */
   int run;				/* run this file was produced from */
   int rerun;				/* rerun number (-1 in photo) */
   int camCol;				/* camera column */
   int field;				/* field number */
   int parent;			        /* id of parent for deblends */
   const int ncolor;			/* number of colours */
   int shallow_copy;			/* is atlas image a shallow copy? */
   OBJMASK *master_mask;                /* master mask for OBJECT1s */

   int drow[NCOLOR], dcol[NCOLOR];	/* offsets wrt reference colour */
   SPANMASK *regmask[NCOLOR];		/* status mask from regions */
   OBJMASK *mask[NCOLOR];		/* object finder's masks */
   PIX *pix[NCOLOR];			/* image values in images for pixels in
					   master_mask, offset by drow, dcol */
/*
 * this field is not saved to disk
 */
   int npix;				/* dimension of pix (in each band) */
} ATLAS_IMAGE;				/* pragma SCHEMA */

/*
 * Here is the internal photo representation of an object in all bands
 */
typedef struct objc {
   const int id;			/* id number for this objc */
   const int ncolor;			/* number of colours */
   OBJECT1 **color;			/* obj1 information from each band */
   OBJ_TYPE type;			/* overall classification */
   float prob_psf;			/* Bayesian probability of being PSF*/
   ATLAS_IMAGE *aimage;			/* atlas image for this object */
   struct test_info *test;		/* information for testers */
   float rowc, rowcErr;			/* row position and error of centre */
   float colc, colcErr;			/* column position and error */
   float rowv, rowvErr;			/* row velocity, in pixels/frame */
   float colv, colvErr;			/* col velocity, in pixels/frame */

   int catID;			        /* catalog id number */
   int flags;				/* c.f. OBJECT1->flags */
   int flags2;				/* c.f. OBJECT1->flags2 */
   int flags3;				/* photo book-keeping; Not Saved */
   PEAKS *peaks;			/* merged peak list */
   int nchild;				/* number of children */
   struct objc *parent;			/* parent if blended */
   struct objc *sibbs;			/* list of siblings */
   struct objc *children;		/* list of children */
} OBJC;					/* pragma SCHEMA */

/*
 * And here is an OBJC reorganised for output, including an Atlas Image
 */
typedef struct {
   int id;				/* id number for this objc */
   int parent;			        /* id of parent for deblends */
   const int ncolor;			/* number of colours */
   int nchild;				/* number of children */
   OBJ_TYPE objc_type;			/* overall classification */
   float objc_prob_psf;			/* Bayesian probability of being PSF */
   int catID;			        /* catalog id number */
   int objc_flags;			/* flags from OBJC */
   int objc_flags2;			/* flags2 from OBJC */
   float objc_rowc, objc_rowcErr;	/* row position and error of centre */
   float objc_colc, objc_colcErr;	/* column position and error */
   float rowv, rowvErr;			/* row velocity, in pixels/frame (NS)*/
   float colv, colvErr;			/* col velocity, in pixels/frame (NS)*/

   ATLAS_IMAGE *aimage;			/* the atlas image */
   struct test_info *test;		/* information for testers */
/*
 * Unpacked OBJECT1s
 */
   float rowc[NCOLOR], rowcErr[NCOLOR];
   float colc[NCOLOR], colcErr[NCOLOR];
   float sky[NCOLOR], skyErr[NCOLOR];
/*
 * PSF and aperture fits/magnitudes
 */
   float psfCounts[NCOLOR], psfCountsErr[NCOLOR];
   float fiberCounts[NCOLOR], fiberCountsErr[NCOLOR];
   float petroCounts[NCOLOR], petroCountsErr[NCOLOR];
   float petroRad[NCOLOR], petroRadErr[NCOLOR];
   float petroR50[NCOLOR], petroR50Err[NCOLOR];
   float petroR90[NCOLOR], petroR90Err[NCOLOR];
/*
 * Shape of object
 */
   float Q[NCOLOR], QErr[NCOLOR], U[NCOLOR], UErr[NCOLOR];

   float M_e1[NCOLOR], M_e2[NCOLOR];
   float M_e1e1Err[NCOLOR], M_e1e2Err[NCOLOR], M_e2e2Err[NCOLOR];
   float M_rr_cc[NCOLOR], M_rr_ccErr[NCOLOR];
   float M_cr4[NCOLOR];
   float M_e1_psf[NCOLOR], M_e2_psf[NCOLOR];
   float M_rr_cc_psf[NCOLOR];
   float M_cr4_psf[NCOLOR];
/*
 * Properties of a certain isophote.
 */
   float iso_rowc[NCOLOR], iso_rowcErr[NCOLOR], iso_rowcGrad[NCOLOR];
   float iso_colc[NCOLOR], iso_colcErr[NCOLOR], iso_colcGrad[NCOLOR];
   float iso_a[NCOLOR], iso_aErr[NCOLOR], iso_aGrad[NCOLOR];
   float iso_b[NCOLOR], iso_bErr[NCOLOR], iso_bGrad[NCOLOR];
   float iso_phi[NCOLOR], iso_phiErr[NCOLOR], iso_phiGrad[NCOLOR];
/*
 * Model parameters for deV and exponential models
 */
   float r_deV[NCOLOR], r_deVErr[NCOLOR];
   float ab_deV[NCOLOR], ab_deVErr[NCOLOR];
   float phi_deV[NCOLOR], phi_deVErr[NCOLOR];
   float counts_deV[NCOLOR], counts_deVErr[NCOLOR];
   float r_exp[NCOLOR], r_expErr[NCOLOR];
   float ab_exp[NCOLOR], ab_expErr[NCOLOR];
   float phi_exp[NCOLOR], phi_expErr[NCOLOR];
   float counts_exp[NCOLOR], counts_expErr[NCOLOR];

   float counts_model[NCOLOR], counts_modelErr[NCOLOR];
/*
 * Measures of image structure
 */
   float texture[NCOLOR];
/*
 * Classification information
 */
   float star_L[NCOLOR], star_lnL[NCOLOR];
   float exp_L[NCOLOR], exp_lnL[NCOLOR];
   float deV_L[NCOLOR], deV_lnL[NCOLOR];
   float fracPSF[NCOLOR];

   int flags[NCOLOR];
   int flags2[NCOLOR];

   OBJ_TYPE type[NCOLOR];
   float prob_psf[NCOLOR];
/*
 * Profile and extent of object
 */
   int nprof[NCOLOR];
   float profMean[NCOLOR][NANN];
   float profErr[NCOLOR][NANN];
} OBJC_IO;				/* pragma SCHEMA */

/*****************************************************************************/
/*
 * Constructors/Destructors
 */
OBJECT1 *phObject1New(void);
void phObject1Del(OBJECT1 *object1);

OBJC *phObjcNew(int ncolor);
OBJC *phObjcNewFromObjc(const OBJC *objc, /* OBJC to copy */
			int deep,	/* copy object1s, ATLAS_IMAGEs etc? */
			int copy_ai);	/* copy atlas images? */
void phObjcDel(OBJC *objc, int deep);
void phObjcChainDel(CHAIN *chain, int deep);

OBJC_IO *phObjcIoNew(int ncolor);
OBJC_IO *phObjcIoNewFromObjc(const OBJC *objc);
void phObjcIoDel(OBJC_IO *objc_io, int deep);

ATLAS_IMAGE *phAtlasImageNew(int ncolor);
ATLAS_IMAGE *phAtlasImageNewFromObjc(const OBJC *objc);
void phAtlasImageDel(ATLAS_IMAGE *aimage, int deep);
void phAtlasImageDelFromObjc(OBJC *objc, int deep);
ATLAS_IMAGE *phAtlasImageCopy(const ATLAS_IMAGE *old, int deep);
int phFindAtlasImagePixel(const ATLAS_IMAGE *aimage, int row, int col);
/*
 * remove the REGIONs from an OBJC
 */
void phRegDelFromObjc(OBJC *objc, int deep);
/*
 * an iterator for OBJCs
 */
OBJC *phObjcDescendentNext(const OBJC *objc);
/*
 * ATLAS_IMAGE utilities
 */
void phAtlasImageCut(OBJC *objc, int color, const FIELDPARAMS *fparams,
		     int val, float sigma, RANDOM *rand);
void phRegionSetFromAtlasImage(const ATLAS_IMAGE *ai,
			       int c, REGION *reg, int row0, int col0,
			       float sky);
void phRegionSetValFromAtlasImage(const ATLAS_IMAGE *ai, int c, REGION *reg,
				  int val, float sigma, RANDOM *rand,
				  int row0, int col0);
void phAtlasImageSetFromRegion(ATLAS_IMAGE *ai, int c,
			       const REGION *data);
void phAtlasImageSetInObjmask(ATLAS_IMAGE *ai, int c,
			      const OBJMASK *omask, const PIX val);
void phAtlasImageSetIfNotInObjmask(ATLAS_IMAGE *ai, int c,
				   const OBJMASK *omask, const PIX val);
void phAtlasImagesPlusEquals(ATLAS_IMAGE *ai, const ATLAS_IMAGE *cai,
			     int bias, int c);
void phAtlasImagesMinusEquals(ATLAS_IMAGE *ai, const ATLAS_IMAGE *cai,
			      int bias, int c);
void phAtlasImageTrimToRect(OBJC *objc, int rmin,int cmin, int rmax,int cmax);

/*
 * Print utilities
 */
void phObject1PrintTerse(const OBJECT1 *obj, FILE *fp);
void phObject1PrintPretty(const OBJECT1 *obj, FILE *fp);
void phObject1PrintCellProf(const OBJECT1 *obj, FILE *fp);
void phObjcListPrint(CHAIN *objclist, char *fname);
void phObjcPrintPretty(OBJC *objc, char *fname);
/*
 * Comparison functions
 */
int phObject1Compare(const OBJECT1 *obj1, const OBJECT1 *obj2);
int phObjcCompare(const OBJC *objc1, const OBJC *objc2);
OBJC *phObjcClosest(CHAIN *chain, float xc, float yc, int color);

#endif


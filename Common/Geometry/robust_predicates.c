#define ANSI_DECLARATORS

#define REAL double

#define INEXACT /* Nothing */
/* #define INEXACT volatile */

/* The vertex data structure.  Each vertex is actually an array of REALs.    */
/*   The number of REALs is unknown until runtime.  An integer boundary      */
/*   marker, and sometimes a pointer to a triangle, is appended after the    */
/*   REALs.                                                                  */

typedef REAL *vertex;

/* Global constants.                                                         */

REAL splitter;       /* Used to split REAL factors for exact multiplication. */
REAL epsilon;                             /* Floating-point machine epsilon. */
REAL resulterrbound;
REAL ccwerrboundA, ccwerrboundB, ccwerrboundC;
REAL iccerrboundA, iccerrboundB, iccerrboundC;
REAL o3derrboundA, o3derrboundB, o3derrboundC;

/********* Geometric primitives begin here                           *********/
/**                                                                         **/
/**                                                                         **/

/* The adaptive exact arithmetic geometric predicates implemented herein are */
/*   described in detail in my paper, "Adaptive Precision Floating-Point     */
/*   Arithmetic and Fast Robust Geometric Predicates."  See the header for a */
/*   full citation.                                                          */

/* Which of the following two methods of finding the absolute values is      */
/*   fastest is compiler-dependent.  A few compilers can inline and optimize */
/*   the fabs() call; but most will incur the overhead of a function call,   */
/*   which is disastrously slow.  A faster way on IEEE machines might be to  */
/*   mask the appropriate bit, but that's difficult to do in C without       */
/*   forcing the value to be stored to memory (rather than be kept in the    */
/*   register to which the optimizer assigned it).                           */

#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a))
/* #define Absolute(a)  fabs(a) */

/* Many of the operations are broken up into two pieces, a main part that    */
/*   performs an approximate operation, and a "tail" that computes the       */
/*   roundoff error of that operation.                                       */
/*                                                                           */
/* The operations Fast_Two_Sum(), Fast_Two_Diff(), Two_Sum(), Two_Diff(),    */
/*   Split(), and Two_Product() are all implemented as described in the      */
/*   reference.  Each of these macros requires certain variables to be       */
/*   defined in the calling routine.  The variables `bvirt', `c', `abig',    */
/*   `_i', `_j', `_k', `_l', `_m', and `_n' are declared `INEXACT' because   */
/*   they store the result of an operation that may incur roundoff error.    */
/*   The input parameter `x' (or the highest numbered `x_' parameter) must   */
/*   also be declared `INEXACT'.                                             */

#define Fast_Two_Sum_Tail(a, b, x, y) \
  bvirt = x - a; \
  y = b - bvirt

#define Fast_Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Fast_Two_Sum_Tail(a, b, x, y)

#define Two_Sum_Tail(a, b, x, y) \
  bvirt = (REAL) (x - a); \
  avirt = x - bvirt; \
  bround = b - bvirt; \
  around = a - avirt; \
  y = around + bround

#define Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Two_Sum_Tail(a, b, x, y)

#define Two_Diff_Tail(a, b, x, y) \
  bvirt = (REAL) (a - x); \
  avirt = x + bvirt; \
  bround = bvirt - b; \
  around = a - avirt; \
  y = around + bround

#define Two_Diff(a, b, x, y) \
  x = (REAL) (a - b); \
  Two_Diff_Tail(a, b, x, y)

#define Split(a, ahi, alo) \
  c = (REAL) (splitter * a); \
  abig = (REAL) (c - a); \
  ahi = c - abig; \
  alo = a - ahi

#define Two_Product_Tail(a, b, x, y) \
  Split(a, ahi, alo); \
  Split(b, bhi, blo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

#define Two_Product(a, b, x, y) \
  x = (REAL) (a * b); \
  Two_Product_Tail(a, b, x, y)

/* Two_Product_Presplit() is Two_Product() where one of the inputs has       */
/*   already been split.  Avoids redundant splitting.                        */

#define Two_Product_Presplit(a, b, bhi, blo, x, y) \
  x = (REAL) (a * b); \
  Split(a, ahi, alo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

/* Square() can be done more quickly than Two_Product().                     */

#define Square_Tail(a, x, y) \
  Split(a, ahi, alo); \
  err1 = x - (ahi * ahi); \
  err3 = err1 - ((ahi + ahi) * alo); \
  y = (alo * alo) - err3

#define Square(a, x, y) \
  x = (REAL) (a * a); \
  Square_Tail(a, x, y)

/* Macros for summing expansions of various fixed lengths.  These are all    */
/*   unrolled versions of Expansion_Sum().                                   */

#define Two_One_Sum(a1, a0, b, x2, x1, x0) \
  Two_Sum(a0, b , _i, x0); \
  Two_Sum(a1, _i, x2, x1)

#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b , _i, x0); \
  Two_Sum( a1, _i, x2, x1)

#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b0, _j, _0, x0); \
  Two_One_Sum(_j, _0, b1, x3, x2, x1)

#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Diff(a1, a0, b0, _j, _0, x0); \
  Two_One_Diff(_j, _0, b1, x3, x2, x1)

/* Macro for multiplying a two-component expansion by a single component.    */

#define Two_One_Product(a1, a0, b, x3, x2, x1, x0) \
  Split(b, bhi, blo); \
  Two_Product_Presplit(a0, b, bhi, blo, _i, x0); \
  Two_Product_Presplit(a1, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x1); \
  Fast_Two_Sum(_j, _k, x3, x2)

/*****************************************************************************/
/*                                                                           */
/*  exactinit()   Initialize the variables used for exact arithmetic.        */
/*                                                                           */
/*  `epsilon' is the largest power of two such that 1.0 + epsilon = 1.0 in   */
/*  floating-point arithmetic.  `epsilon' bounds the relative roundoff       */
/*  error.  It is used for floating-point error analysis.                    */
/*                                                                           */
/*  `splitter' is used to split floating-point numbers into two half-        */
/*  length significands for exact multiplication.                            */
/*                                                                           */
/*  I imagine that a highly optimizing compiler might be too smart for its   */
/*  own good, and somehow cause this routine to fail, if it pretends that    */
/*  floating-point arithmetic is too much like real arithmetic.              */
/*                                                                           */
/*  Don't change this routine unless you fully understand it.                */
/*                                                                           */
/*****************************************************************************/

void exactinit()
{
	REAL half;
	REAL check, lastcheck;
	int every_other;
#ifdef LINUX
	int cword;
#endif /* LINUX */

#ifdef CPU86
#ifdef SINGLE
	_control87(_PC_24, _MCW_PC); /* Set FPU control word for single precision. */
#else /* not SINGLE */
	_control87(_PC_53, _MCW_PC); /* Set FPU control word for double precision. */
#endif /* not SINGLE */
#endif /* CPU86 */
#ifdef LINUX
#ifdef SINGLE
								 /*  cword = 4223; */
	cword = 4210;                 /* set FPU control word for single precision */
#else /* not SINGLE */
								 /*  cword = 4735; */
	cword = 4722;                 /* set FPU control word for double precision */
#endif /* not SINGLE */
	_FPU_SETCW(cword);
#endif /* LINUX */

	every_other = 1;
	half = 0.5;
	epsilon = 1.0;
	splitter = 1.0;
	check = 1.0;
	/* Repeatedly divide `epsilon' by two until it is too small to add to      */
	/*   one without causing roundoff.  (Also check if the sum is equal to     */
	/*   the previous sum, for machines that round up instead of using exact   */
	/*   rounding.  Not that these routines will work on such machines.)       */
	do {
		lastcheck = check;
		epsilon *= half;
		if (every_other) {
			splitter *= 2.0;
		}
		every_other = !every_other;
		check = 1.0 + epsilon;
	} while ((check != 1.0) && (check != lastcheck));
	splitter += 1.0;
	/* Error bounds for orientation and incircle tests. */
	resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
	ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon;
	ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon;
	ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
	iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon;
	iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon;
	iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
	o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
	o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
	o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
}

/*****************************************************************************/
/*                                                                           */
/*  fast_expansion_sum_zeroelim()   Sum two expansions, eliminating zero     */
/*                                  components from the output expansion.    */
/*                                                                           */
/*  Sets h = e + f.  See my Robust Predicates paper for details.             */
/*                                                                           */
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
/*  properties.                                                              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
int fast_expansion_sum_zeroelim(int elen, REAL *e, int flen, REAL *f, REAL *h)
#else /* not ANSI_DECLARATORS */
int fast_expansion_sum_zeroelim(elen, e, flen, f, h)  /* h cannot be e or f. */
int elen;
REAL *e;
int flen;
REAL *f;
REAL *h;
#endif /* not ANSI_DECLARATORS */

{
	REAL Q;
	INEXACT REAL Qnew;
	INEXACT REAL hh;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	int eindex, findex, hindex;
	REAL enow, fnow;

	enow = e[0];
	fnow = f[0];
	eindex = findex = 0;
	if ((fnow > enow) == (fnow > -enow)) {
		Q = enow;
		enow = e[++eindex];
	}
	else {
		Q = fnow;
		fnow = f[++findex];
	}
	hindex = 0;
	if ((eindex < elen) && (findex < flen)) {
		if ((fnow > enow) == (fnow > -enow)) {
			Fast_Two_Sum(enow, Q, Qnew, hh);
			enow = e[++eindex];
		}
		else {
			Fast_Two_Sum(fnow, Q, Qnew, hh);
			fnow = f[++findex];
		}
		Q = Qnew;
		if (hh != 0.0) {
			h[hindex++] = hh;
		}
		while ((eindex < elen) && (findex < flen)) {
			if ((fnow > enow) == (fnow > -enow)) {
				Two_Sum(Q, enow, Qnew, hh);
				enow = e[++eindex];
			}
			else {
				Two_Sum(Q, fnow, Qnew, hh);
				fnow = f[++findex];
			}
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
		}
	}
	while (eindex < elen) {
		Two_Sum(Q, enow, Qnew, hh);
		enow = e[++eindex];
		Q = Qnew;
		if (hh != 0.0) {
			h[hindex++] = hh;
		}
	}
	while (findex < flen) {
		Two_Sum(Q, fnow, Qnew, hh);
		fnow = f[++findex];
		Q = Qnew;
		if (hh != 0.0) {
			h[hindex++] = hh;
		}
	}
	if ((Q != 0.0) || (hindex == 0)) {
		h[hindex++] = Q;
	}
	return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  scale_expansion_zeroelim()   Multiply an expansion by a scalar,          */
/*                               eliminating zero components from the        */
/*                               output expansion.                           */
/*                                                                           */
/*  Sets h = be.  See my Robust Predicates paper for details.                */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
/*  properties as well.  (That is, if e has one of these properties, so      */
/*  will h.)                                                                 */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
int scale_expansion_zeroelim(int elen, REAL *e, REAL b, REAL *h)
#else /* not ANSI_DECLARATORS */
int scale_expansion_zeroelim(elen, e, b, h)   /* e and h cannot be the same. */
int elen;
REAL *e;
REAL b;
REAL *h;
#endif /* not ANSI_DECLARATORS */

{
	INEXACT REAL Q, sum;
	REAL hh;
	INEXACT REAL product1;
	REAL product0;
	int eindex, hindex;
	REAL enow;
	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;

	Split(b, bhi, blo);
	Two_Product_Presplit(e[0], b, bhi, blo, Q, hh);
	hindex = 0;
	if (hh != 0) {
		h[hindex++] = hh;
	}
	for (eindex = 1; eindex < elen; eindex++) {
		enow = e[eindex];
		Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
		Two_Sum(Q, product0, sum, hh);
		if (hh != 0) {
			h[hindex++] = hh;
		}
		Fast_Two_Sum(product1, sum, Q, hh);
		if (hh != 0) {
			h[hindex++] = hh;
		}
	}
	if ((Q != 0.0) || (hindex == 0)) {
		h[hindex++] = Q;
	}
	return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  estimate()   Produce a one-word estimate of an expansion's value.        */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
REAL estimate(int elen, REAL *e)
#else /* not ANSI_DECLARATORS */
REAL estimate(elen, e)
int elen;
REAL *e;
#endif /* not ANSI_DECLARATORS */

{
	REAL Q;
	int eindex;

	Q = e[0];
	for (eindex = 1; eindex < elen; eindex++) {
		Q += e[eindex];
	}
	return Q;
}

/*****************************************************************************/
/*                                                                           */
/*  counterclockwise()   Return a positive value if the points pa, pb, and   */
/*                       pc occur in counterclockwise order; a negative      */
/*                       value if they occur in clockwise order; and zero    */
/*                       if they are collinear.  The result is also a rough  */
/*                       approximation of twice the signed area of the       */
/*                       triangle defined by the three points.               */
/*                                                                           */
/*  Uses exact arithmetic if necessary to ensure a correct answer.  The      */
/*  result returned is the determinant of a matrix.  This determinant is     */
/*  computed adaptively, in the sense that exact arithmetic is used only to  */
/*  the degree it is needed to ensure that the returned value has the        */
/*  correct sign.  Hence, this function is usually quite fast, but will run  */
/*  more slowly when the input points are collinear or nearly so.            */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
REAL counterclockwiseadapt(vertex pa, vertex pb, vertex pc, REAL detsum)
#else /* not ANSI_DECLARATORS */
REAL counterclockwiseadapt(pa, pb, pc, detsum)
vertex pa;
vertex pb;
vertex pc;
REAL detsum;
#endif /* not ANSI_DECLARATORS */

{
	INEXACT REAL acx, acy, bcx, bcy;
	REAL acxtail, acytail, bcxtail, bcytail;
	INEXACT REAL detleft, detright;
	REAL detlefttail, detrighttail;
	REAL det, errbound;
	REAL B[4], C1[8], C2[12], D[16];
	INEXACT REAL B3;
	int C1length, C2length, Dlength;
	REAL u[4];
	INEXACT REAL u3;
	INEXACT REAL s1, t1;
	REAL s0, t0;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j;
	REAL _0;

	acx = (REAL)(pa[0] - pc[0]);
	bcx = (REAL)(pb[0] - pc[0]);
	acy = (REAL)(pa[1] - pc[1]);
	bcy = (REAL)(pb[1] - pc[1]);

	Two_Product(acx, bcy, detleft, detlefttail);
	Two_Product(acy, bcx, detright, detrighttail);

	Two_Two_Diff(detleft, detlefttail, detright, detrighttail,
		B3, B[2], B[1], B[0]);
	B[3] = B3;

	det = estimate(4, B);
	errbound = ccwerrboundB * detsum;
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	Two_Diff_Tail(pa[0], pc[0], acx, acxtail);
	Two_Diff_Tail(pb[0], pc[0], bcx, bcxtail);
	Two_Diff_Tail(pa[1], pc[1], acy, acytail);
	Two_Diff_Tail(pb[1], pc[1], bcy, bcytail);

	if ((acxtail == 0.0) && (acytail == 0.0)
		&& (bcxtail == 0.0) && (bcytail == 0.0)) {
		return det;
	}

	errbound = ccwerrboundC * detsum + resulterrbound * Absolute(det);
	det += (acx * bcytail + bcy * acxtail)
		- (acy * bcxtail + bcx * acytail);
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	Two_Product(acxtail, bcy, s1, s0);
	Two_Product(acytail, bcx, t1, t0);
	Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
	u[3] = u3;
	C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1);

	Two_Product(acx, bcytail, s1, s0);
	Two_Product(acy, bcxtail, t1, t0);
	Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
	u[3] = u3;
	C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2);

	Two_Product(acxtail, bcytail, s1, s0);
	Two_Product(acytail, bcxtail, t1, t0);
	Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
	u[3] = u3;
	Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D);

	return(D[Dlength - 1]);
}

#ifdef ANSI_DECLARATORS
REAL counterclockwise(// my_change //struct mesh *m, struct behavior *b,
	vertex pa, vertex pb, vertex pc)
#else /* not ANSI_DECLARATORS */
REAL counterclockwise(m, b, pa, pb, pc)
struct mesh *m;
struct behavior *b;
vertex pa;
vertex pb;
vertex pc;
#endif /* not ANSI_DECLARATORS */

{
	REAL detleft, detright, det;
	REAL detsum, errbound;

	// my_change //m->counterclockcount++;

	detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
	detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
	det = detleft - detright;

	// my_change //if (b->noexact) {
	// my_change //	return det;
	// my_change //}

	if (detleft > 0.0) {
		if (detright <= 0.0) {
			return det;
		}
		else {
			detsum = detleft + detright;
		}
	}
	else if (detleft < 0.0) {
		if (detright >= 0.0) {
			return det;
		}
		else {
			detsum = -detleft - detright;
		}
	}
	else {
		return det;
	}

	errbound = ccwerrboundA * detsum;
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	return counterclockwiseadapt(pa, pb, pc, detsum);
}

/*****************************************************************************/
/*                                                                           */
/*  incircle()   Return a positive value if the point pd lies inside the     */
/*               circle passing through pa, pb, and pc; a negative value if  */
/*               it lies outside; and zero if the four points are cocircular.*/
/*               The points pa, pb, and pc must be in counterclockwise       */
/*               order, or the sign of the result will be reversed.          */
/*                                                                           */
/*  Uses exact arithmetic if necessary to ensure a correct answer.  The      */
/*  result returned is the determinant of a matrix.  This determinant is     */
/*  computed adaptively, in the sense that exact arithmetic is used only to  */
/*  the degree it is needed to ensure that the returned value has the        */
/*  correct sign.  Hence, this function is usually quite fast, but will run  */
/*  more slowly when the input points are cocircular or nearly so.           */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
REAL incircleadapt(vertex pa, vertex pb, vertex pc, vertex pd, REAL permanent)
#else /* not ANSI_DECLARATORS */
REAL incircleadapt(pa, pb, pc, pd, permanent)
vertex pa;
vertex pb;
vertex pc;
vertex pd;
REAL permanent;
#endif /* not ANSI_DECLARATORS */

{
	INEXACT REAL adx, bdx, cdx, ady, bdy, cdy;
	REAL det, errbound;

	INEXACT REAL bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
	REAL bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
	REAL bc[4], ca[4], ab[4];
	INEXACT REAL bc3, ca3, ab3;
	REAL axbc[8], axxbc[16], aybc[8], ayybc[16], adet[32];
	int axbclen, axxbclen, aybclen, ayybclen, alen;
	REAL bxca[8], bxxca[16], byca[8], byyca[16], bdet[32];
	int bxcalen, bxxcalen, bycalen, byycalen, blen;
	REAL cxab[8], cxxab[16], cyab[8], cyyab[16], cdet[32];
	int cxablen, cxxablen, cyablen, cyyablen, clen;
	REAL abdet[64];
	int ablen;
	REAL fin1[1152], fin2[1152];
	REAL *finnow, *finother, *finswap;
	int finlength;

	REAL adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;
	INEXACT REAL adxadx1, adyady1, bdxbdx1, bdybdy1, cdxcdx1, cdycdy1;
	REAL adxadx0, adyady0, bdxbdx0, bdybdy0, cdxcdx0, cdycdy0;
	REAL aa[4], bb[4], cc[4];
	INEXACT REAL aa3, bb3, cc3;
	INEXACT REAL ti1, tj1;
	REAL ti0, tj0;
	REAL u[4], v[4];
	INEXACT REAL u3, v3;
	REAL temp8[8], temp16a[16], temp16b[16], temp16c[16];
	REAL temp32a[32], temp32b[32], temp48[48], temp64[64];
	int temp8len, temp16alen, temp16blen, temp16clen;
	int temp32alen, temp32blen, temp48len, temp64len;
	REAL axtbb[8], axtcc[8], aytbb[8], aytcc[8];
	int axtbblen, axtcclen, aytbblen, aytcclen;
	REAL bxtaa[8], bxtcc[8], bytaa[8], bytcc[8];
	int bxtaalen, bxtcclen, bytaalen, bytcclen;
	REAL cxtaa[8], cxtbb[8], cytaa[8], cytbb[8];
	int cxtaalen, cxtbblen, cytaalen, cytbblen;
	REAL axtbc[8], aytbc[8], bxtca[8], bytca[8], cxtab[8], cytab[8];
	int axtbclen, aytbclen, bxtcalen, bytcalen, cxtablen, cytablen;
	REAL axtbct[16], aytbct[16], bxtcat[16], bytcat[16], cxtabt[16], cytabt[16];
	int axtbctlen, aytbctlen, bxtcatlen, bytcatlen, cxtabtlen, cytabtlen;
	REAL axtbctt[8], aytbctt[8], bxtcatt[8];
	REAL bytcatt[8], cxtabtt[8], cytabtt[8];
	int axtbcttlen, aytbcttlen, bxtcattlen, bytcattlen, cxtabttlen, cytabttlen;
	REAL abt[8], bct[8], cat[8];
	int abtlen, bctlen, catlen;
	REAL abtt[4], bctt[4], catt[4];
	int abttlen, bcttlen, cattlen;
	INEXACT REAL abtt3, bctt3, catt3;
	REAL negate;

	INEXACT REAL bvirt;
	REAL avirt, bround, around;
	INEXACT REAL c;
	INEXACT REAL abig;
	REAL ahi, alo, bhi, blo;
	REAL err1, err2, err3;
	INEXACT REAL _i, _j;
	REAL _0;

	adx = (REAL)(pa[0] - pd[0]);
	bdx = (REAL)(pb[0] - pd[0]);
	cdx = (REAL)(pc[0] - pd[0]);
	ady = (REAL)(pa[1] - pd[1]);
	bdy = (REAL)(pb[1] - pd[1]);
	cdy = (REAL)(pc[1] - pd[1]);

	Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
	Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
	Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
	bc[3] = bc3;
	axbclen = scale_expansion_zeroelim(4, bc, adx, axbc);
	axxbclen = scale_expansion_zeroelim(axbclen, axbc, adx, axxbc);
	aybclen = scale_expansion_zeroelim(4, bc, ady, aybc);
	ayybclen = scale_expansion_zeroelim(aybclen, aybc, ady, ayybc);
	alen = fast_expansion_sum_zeroelim(axxbclen, axxbc, ayybclen, ayybc, adet);

	Two_Product(cdx, ady, cdxady1, cdxady0);
	Two_Product(adx, cdy, adxcdy1, adxcdy0);
	Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
	ca[3] = ca3;
	bxcalen = scale_expansion_zeroelim(4, ca, bdx, bxca);
	bxxcalen = scale_expansion_zeroelim(bxcalen, bxca, bdx, bxxca);
	bycalen = scale_expansion_zeroelim(4, ca, bdy, byca);
	byycalen = scale_expansion_zeroelim(bycalen, byca, bdy, byyca);
	blen = fast_expansion_sum_zeroelim(bxxcalen, bxxca, byycalen, byyca, bdet);

	Two_Product(adx, bdy, adxbdy1, adxbdy0);
	Two_Product(bdx, ady, bdxady1, bdxady0);
	Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
	ab[3] = ab3;
	cxablen = scale_expansion_zeroelim(4, ab, cdx, cxab);
	cxxablen = scale_expansion_zeroelim(cxablen, cxab, cdx, cxxab);
	cyablen = scale_expansion_zeroelim(4, ab, cdy, cyab);
	cyyablen = scale_expansion_zeroelim(cyablen, cyab, cdy, cyyab);
	clen = fast_expansion_sum_zeroelim(cxxablen, cxxab, cyyablen, cyyab, cdet);

	ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
	finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

	det = estimate(finlength, fin1);
	errbound = iccerrboundB * permanent;
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
	Two_Diff_Tail(pa[1], pd[1], ady, adytail);
	Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
	Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
	Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
	Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
	if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
		&& (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)) {
		return det;
	}

	errbound = iccerrboundC * permanent + resulterrbound * Absolute(det);
	det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail)
		- (bdy * cdxtail + cdx * bdytail))
		+ 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
		+ ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail)
			- (cdy * adxtail + adx * cdytail))
			+ 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
		+ ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail)
			- (ady * bdxtail + bdx * adytail))
			+ 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
	if ((det >= errbound) || (-det >= errbound)) {
		return det;
	}

	finnow = fin1;
	finother = fin2;

	if ((bdxtail != 0.0) || (bdytail != 0.0)
		|| (cdxtail != 0.0) || (cdytail != 0.0)) {
		Square(adx, adxadx1, adxadx0);
		Square(ady, adyady1, adyady0);
		Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0, aa3, aa[2], aa[1], aa[0]);
		aa[3] = aa3;
	}
	if ((cdxtail != 0.0) || (cdytail != 0.0)
		|| (adxtail != 0.0) || (adytail != 0.0)) {
		Square(bdx, bdxbdx1, bdxbdx0);
		Square(bdy, bdybdy1, bdybdy0);
		Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, bb3, bb[2], bb[1], bb[0]);
		bb[3] = bb3;
	}
	if ((adxtail != 0.0) || (adytail != 0.0)
		|| (bdxtail != 0.0) || (bdytail != 0.0)) {
		Square(cdx, cdxcdx1, cdxcdx0);
		Square(cdy, cdycdy1, cdycdy0);
		Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, cc3, cc[2], cc[1], cc[0]);
		cc[3] = cc3;
	}

	if (adxtail != 0.0) {
		axtbclen = scale_expansion_zeroelim(4, bc, adxtail, axtbc);
		temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, 2.0 * adx,
			temp16a);

		axtcclen = scale_expansion_zeroelim(4, cc, adxtail, axtcc);
		temp16blen = scale_expansion_zeroelim(axtcclen, axtcc, bdy, temp16b);

		axtbblen = scale_expansion_zeroelim(4, bb, adxtail, axtbb);
		temp16clen = scale_expansion_zeroelim(axtbblen, axtbb, -cdy, temp16c);

		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
			temp16blen, temp16b, temp32a);
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
			temp32alen, temp32a, temp48);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
			temp48, finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (adytail != 0.0) {
		aytbclen = scale_expansion_zeroelim(4, bc, adytail, aytbc);
		temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, 2.0 * ady,
			temp16a);

		aytbblen = scale_expansion_zeroelim(4, bb, adytail, aytbb);
		temp16blen = scale_expansion_zeroelim(aytbblen, aytbb, cdx, temp16b);

		aytcclen = scale_expansion_zeroelim(4, cc, adytail, aytcc);
		temp16clen = scale_expansion_zeroelim(aytcclen, aytcc, -bdx, temp16c);

		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
			temp16blen, temp16b, temp32a);
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
			temp32alen, temp32a, temp48);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
			temp48, finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (bdxtail != 0.0) {
		bxtcalen = scale_expansion_zeroelim(4, ca, bdxtail, bxtca);
		temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, 2.0 * bdx,
			temp16a);

		bxtaalen = scale_expansion_zeroelim(4, aa, bdxtail, bxtaa);
		temp16blen = scale_expansion_zeroelim(bxtaalen, bxtaa, cdy, temp16b);

		bxtcclen = scale_expansion_zeroelim(4, cc, bdxtail, bxtcc);
		temp16clen = scale_expansion_zeroelim(bxtcclen, bxtcc, -ady, temp16c);

		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
			temp16blen, temp16b, temp32a);
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
			temp32alen, temp32a, temp48);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
			temp48, finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (bdytail != 0.0) {
		bytcalen = scale_expansion_zeroelim(4, ca, bdytail, bytca);
		temp16alen = scale_expansion_zeroelim(bytcalen, bytca, 2.0 * bdy,
			temp16a);

		bytcclen = scale_expansion_zeroelim(4, cc, bdytail, bytcc);
		temp16blen = scale_expansion_zeroelim(bytcclen, bytcc, adx, temp16b);

		bytaalen = scale_expansion_zeroelim(4, aa, bdytail, bytaa);
		temp16clen = scale_expansion_zeroelim(bytaalen, bytaa, -cdx, temp16c);

		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
			temp16blen, temp16b, temp32a);
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
			temp32alen, temp32a, temp48);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
			temp48, finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (cdxtail != 0.0) {
		cxtablen = scale_expansion_zeroelim(4, ab, cdxtail, cxtab);
		temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, 2.0 * cdx,
			temp16a);

		cxtbblen = scale_expansion_zeroelim(4, bb, cdxtail, cxtbb);
		temp16blen = scale_expansion_zeroelim(cxtbblen, cxtbb, ady, temp16b);

		cxtaalen = scale_expansion_zeroelim(4, aa, cdxtail, cxtaa);
		temp16clen = scale_expansion_zeroelim(cxtaalen, cxtaa, -bdy, temp16c);

		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
			temp16blen, temp16b, temp32a);
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
			temp32alen, temp32a, temp48);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
			temp48, finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}
	if (cdytail != 0.0) {
		cytablen = scale_expansion_zeroelim(4, ab, cdytail, cytab);
		temp16alen = scale_expansion_zeroelim(cytablen, cytab, 2.0 * cdy,
			temp16a);

		cytaalen = scale_expansion_zeroelim(4, aa, cdytail, cytaa);
		temp16blen = scale_expansion_zeroelim(cytaalen, cytaa, bdx, temp16b);

		cytbblen = scale_expansion_zeroelim(4, bb, cdytail, cytbb);
		temp16clen = scale_expansion_zeroelim(cytbblen, cytbb, -adx, temp16c);

		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
			temp16blen, temp16b, temp32a);
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
			temp32alen, temp32a, temp48);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
			temp48, finother);
		finswap = finnow; finnow = finother; finother = finswap;
	}

	if ((adxtail != 0.0) || (adytail != 0.0)) {
		if ((bdxtail != 0.0) || (bdytail != 0.0)
			|| (cdxtail != 0.0) || (cdytail != 0.0)) {
			Two_Product(bdxtail, cdy, ti1, ti0);
			Two_Product(bdx, cdytail, tj1, tj0);
			Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
			u[3] = u3;
			negate = -bdy;
			Two_Product(cdxtail, negate, ti1, ti0);
			negate = -bdytail;
			Two_Product(cdx, negate, tj1, tj0);
			Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
			v[3] = v3;
			bctlen = fast_expansion_sum_zeroelim(4, u, 4, v, bct);

			Two_Product(bdxtail, cdytail, ti1, ti0);
			Two_Product(cdxtail, bdytail, tj1, tj0);
			Two_Two_Diff(ti1, ti0, tj1, tj0, bctt3, bctt[2], bctt[1], bctt[0]);
			bctt[3] = bctt3;
			bcttlen = 4;
		}
		else {
			bct[0] = 0.0;
			bctlen = 1;
			bctt[0] = 0.0;
			bcttlen = 1;
		}

		if (adxtail != 0.0) {
			temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, adxtail, temp16a);
			axtbctlen = scale_expansion_zeroelim(bctlen, bct, adxtail, axtbct);
			temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, 2.0 * adx,
				temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;
			if (bdytail != 0.0) {
				temp8len = scale_expansion_zeroelim(4, cc, adxtail, temp8);
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail,
					temp16a);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
					temp16a, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (cdytail != 0.0) {
				temp8len = scale_expansion_zeroelim(4, bb, -adxtail, temp8);
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail,
					temp16a);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
					temp16a, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}

			temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, adxtail,
				temp32a);
			axtbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adxtail, axtbctt);
			temp16alen = scale_expansion_zeroelim(axtbcttlen, axtbctt, 2.0 * adx,
				temp16a);
			temp16blen = scale_expansion_zeroelim(axtbcttlen, axtbctt, adxtail,
				temp16b);
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32b);
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
				temp32blen, temp32b, temp64);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
				temp64, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (adytail != 0.0) {
			temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, adytail, temp16a);
			aytbctlen = scale_expansion_zeroelim(bctlen, bct, adytail, aytbct);
			temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, 2.0 * ady,
				temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;


			temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, adytail,
				temp32a);
			aytbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adytail, aytbctt);
			temp16alen = scale_expansion_zeroelim(aytbcttlen, aytbctt, 2.0 * ady,
				temp16a);
			temp16blen = scale_expansion_zeroelim(aytbcttlen, aytbctt, adytail,
				temp16b);
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32b);
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
				temp32blen, temp32b, temp64);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
				temp64, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
	}
	if ((bdxtail != 0.0) || (bdytail != 0.0)) {
		if ((cdxtail != 0.0) || (cdytail != 0.0)
			|| (adxtail != 0.0) || (adytail != 0.0)) {
			Two_Product(cdxtail, ady, ti1, ti0);
			Two_Product(cdx, adytail, tj1, tj0);
			Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
			u[3] = u3;
			negate = -cdy;
			Two_Product(adxtail, negate, ti1, ti0);
			negate = -cdytail;
			Two_Product(adx, negate, tj1, tj0);
			Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
			v[3] = v3;
			catlen = fast_expansion_sum_zeroelim(4, u, 4, v, cat);

			Two_Product(cdxtail, adytail, ti1, ti0);
			Two_Product(adxtail, cdytail, tj1, tj0);
			Two_Two_Diff(ti1, ti0, tj1, tj0, catt3, catt[2], catt[1], catt[0]);
			catt[3] = catt3;
			cattlen = 4;
		}
		else {
			cat[0] = 0.0;
			catlen = 1;
			catt[0] = 0.0;
			cattlen = 1;
		}

		if (bdxtail != 0.0) {
			temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, bdxtail, temp16a);
			bxtcatlen = scale_expansion_zeroelim(catlen, cat, bdxtail, bxtcat);
			temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, 2.0 * bdx,
				temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;
			if (cdytail != 0.0) {
				temp8len = scale_expansion_zeroelim(4, aa, bdxtail, temp8);
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail,
					temp16a);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
					temp16a, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (adytail != 0.0) {
				temp8len = scale_expansion_zeroelim(4, cc, -bdxtail, temp8);
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail,
					temp16a);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
					temp16a, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}

			temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, bdxtail,
				temp32a);
			bxtcattlen = scale_expansion_zeroelim(cattlen, catt, bdxtail, bxtcatt);
			temp16alen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, 2.0 * bdx,
				temp16a);
			temp16blen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, bdxtail,
				temp16b);
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32b);
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
				temp32blen, temp32b, temp64);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
				temp64, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (bdytail != 0.0) {
			temp16alen = scale_expansion_zeroelim(bytcalen, bytca, bdytail, temp16a);
			bytcatlen = scale_expansion_zeroelim(catlen, cat, bdytail, bytcat);
			temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, 2.0 * bdy,
				temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;


			temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, bdytail,
				temp32a);
			bytcattlen = scale_expansion_zeroelim(cattlen, catt, bdytail, bytcatt);
			temp16alen = scale_expansion_zeroelim(bytcattlen, bytcatt, 2.0 * bdy,
				temp16a);
			temp16blen = scale_expansion_zeroelim(bytcattlen, bytcatt, bdytail,
				temp16b);
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32b);
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
				temp32blen, temp32b, temp64);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
				temp64, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
	}
	if ((cdxtail != 0.0) || (cdytail != 0.0)) {
		if ((adxtail != 0.0) || (adytail != 0.0)
			|| (bdxtail != 0.0) || (bdytail != 0.0)) {
			Two_Product(adxtail, bdy, ti1, ti0);
			Two_Product(adx, bdytail, tj1, tj0);
			Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
			u[3] = u3;
			negate = -ady;
			Two_Product(bdxtail, negate, ti1, ti0);
			negate = -adytail;
			Two_Product(bdx, negate, tj1, tj0);
			Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
			v[3] = v3;
			abtlen = fast_expansion_sum_zeroelim(4, u, 4, v, abt);

			Two_Product(adxtail, bdytail, ti1, ti0);
			Two_Product(bdxtail, adytail, tj1, tj0);
			Two_Two_Diff(ti1, ti0, tj1, tj0, abtt3, abtt[2], abtt[1], abtt[0]);
			abtt[3] = abtt3;
			abttlen = 4;
		}
		else {
			abt[0] = 0.0;
			abtlen = 1;
			abtt[0] = 0.0;
			abttlen = 1;
		}

		if (cdxtail != 0.0) {
			temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, cdxtail, temp16a);
			cxtabtlen = scale_expansion_zeroelim(abtlen, abt, cdxtail, cxtabt);
			temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, 2.0 * cdx,
				temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;
			if (adytail != 0.0) {
				temp8len = scale_expansion_zeroelim(4, bb, cdxtail, temp8);
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail,
					temp16a);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
					temp16a, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (bdytail != 0.0) {
				temp8len = scale_expansion_zeroelim(4, aa, -cdxtail, temp8);
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail,
					temp16a);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
					temp16a, finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}

			temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, cdxtail,
				temp32a);
			cxtabttlen = scale_expansion_zeroelim(abttlen, abtt, cdxtail, cxtabtt);
			temp16alen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, 2.0 * cdx,
				temp16a);
			temp16blen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, cdxtail,
				temp16b);
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32b);
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
				temp32blen, temp32b, temp64);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
				temp64, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (cdytail != 0.0) {
			temp16alen = scale_expansion_zeroelim(cytablen, cytab, cdytail, temp16a);
			cytabtlen = scale_expansion_zeroelim(abtlen, abt, cdytail, cytabt);
			temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, 2.0 * cdy,
				temp32a);
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp32alen, temp32a, temp48);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
				temp48, finother);
			finswap = finnow; finnow = finother; finother = finswap;


			temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, cdytail,
				temp32a);
			cytabttlen = scale_expansion_zeroelim(abttlen, abtt, cdytail, cytabtt);
			temp16alen = scale_expansion_zeroelim(cytabttlen, cytabtt, 2.0 * cdy,
				temp16a);
			temp16blen = scale_expansion_zeroelim(cytabttlen, cytabtt, cdytail,
				temp16b);
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
				temp16blen, temp16b, temp32b);
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
				temp32blen, temp32b, temp64);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
				temp64, finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
	}

	return finnow[finlength - 1];
}

#ifdef ANSI_DECLARATORS
REAL incircle(// my_change //struct mesh *m, struct behavior *b,
	vertex pa, vertex pb, vertex pc, vertex pd)
#else /* not ANSI_DECLARATORS */
REAL incircle(m, b, pa, pb, pc, pd)
struct mesh *m;
struct behavior *b;
vertex pa;
vertex pb;
vertex pc;
vertex pd;
#endif /* not ANSI_DECLARATORS */

{
	REAL adx, bdx, cdx, ady, bdy, cdy;
	REAL bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
	REAL alift, blift, clift;
	REAL det;
	REAL permanent, errbound;

	// my_change //m->incirclecount++;

	adx = pa[0] - pd[0];
	bdx = pb[0] - pd[0];
	cdx = pc[0] - pd[0];
	ady = pa[1] - pd[1];
	bdy = pb[1] - pd[1];
	cdy = pc[1] - pd[1];

	bdxcdy = bdx * cdy;
	cdxbdy = cdx * bdy;
	alift = adx * adx + ady * ady;

	cdxady = cdx * ady;
	adxcdy = adx * cdy;
	blift = bdx * bdx + bdy * bdy;

	adxbdy = adx * bdy;
	bdxady = bdx * ady;
	clift = cdx * cdx + cdy * cdy;

	det = alift * (bdxcdy - cdxbdy)
		+ blift * (cdxady - adxcdy)
		+ clift * (adxbdy - bdxady);

	// my_change //if (b->noexact) {
	// my_change //	return det;
	// my_change //}

	permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * alift
		+ (Absolute(cdxady) + Absolute(adxcdy)) * blift
		+ (Absolute(adxbdy) + Absolute(bdxady)) * clift;
	errbound = iccerrboundA * permanent;
	if ((det > errbound) || (-det > errbound)) {
		return det;
	}

	return incircleadapt(pa, pb, pc, pd, permanent);
}

/**                                                                         **/
/**                                                                         **/
/********* Geometric primitives end here                             *********/
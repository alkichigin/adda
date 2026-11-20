/* Inline complex functions, functions on length-3 real and complex vectors, and several auxiliary functions
 *
 * Copyright (C) ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#ifndef __cmplx_h
#define __cmplx_h

// project headers
#include "const.h"    // for math constants
#include "types.h"    // for doublecomplex
// system headers
#include <math.h>
#include <string.h>   // for memcpy

// Uncomment this to turn off calculation of imExp using tables
//#define NO_IMEXP_TABLE

#ifdef USE_SSE3
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#endif

// useful macro for printing complex numbers and matrices
#define REIM(a) creal(a),cimag(a)
#define REIM3V(a) REIM((a)[0]),REIM((a)[1]),REIM((a)[2])

extern void cik01_(doublecomplex *z, doublecomplex *cbi0, doublecomplex *cdi0, doublecomplex *cbi1, doublecomplex *cdi1, doublecomplex *cbk0, doublecomplex *cdk0, doublecomplex *cbk1, doublecomplex *cdk1);

#ifndef NO_IMEXP_TABLE
void imExpTableInit(void);
doublecomplex imExpTable(double arg);
#endif
void imExp_arr(doublecomplex arg,int size,doublecomplex *c);

/* We do not use 'restrict' in the following functions since they are all inline - compiler will optimize the code
 * inside the calling function and decide whether the arrays can alias or not.
 */

//======================================================================================================================
// operations on complex numbers

static inline double cAbs2(const doublecomplex a)
// square of absolute value of complex number; |a|^2
{
	return creal(a)*creal(a) + cimag(a)*cimag(a);
}

//======================================================================================================================

static inline doublecomplex cSqrtCut(const doublecomplex a)
// square root of complex number, with explicit handling of branch cut (not to depend on sign of zero of imaginary part)
/* It is designed for calculating normal component of the transmitted wavevector when passing through the plane
 * interface. However, such choice of branch cut (while physically correct) leads to all kind of weird consequences.
 *
 * For instance, the electric field above the interface for plane wave propagating from a slightly absorbing substrate
 * at large incident angle (larger than critical angle for purely real refractive index) is unexpectedly large. This
 * happens because the wave in the vacuum is inhomogeneous and the real part of wavevector is almost parallel to the
 * surface. So the field above the surface actually comes from distant points on the surface, which has much larger
 * amplitude of the incident wave from below (compared to that under the observation point). Since the distance along
 * the surface (or the corresponding slope) is inversely proportional to the imaginary part of the substrate refractive
 * index, the effect remains finite even in the limit of absorption going to zero. Therefore, in this case there exists
 * a discontinuity when switching from non-absorbing to absorbing substrate. Physically, this fact is a consequence of
 * the infinite lateral extent of the plane wave.
 *
 * Exactly the same issue exist when scattering into the absorbing medium is calculated. At large scattering angles the
 * amplitude becomes very large, which also amplifies a lot the calculated Csca.
 */
{
	if (cimag(a)==0) {
		if (creal(a)>=0) return sqrt(creal(a));
		else return I*sqrt(-creal(a));
	}
	else return csqrt(a);
}

//======================================================================================================================

static inline doublecomplex imExpReal(const double arg)
/* exponent of imaginary argument Exp(i*arg)
 * !!! should not be used in parameter parsing (table is initialized in VariablesInterconnect())
 */
{
#ifdef NO_IMEXP_TABLE
	/* We tried different standard options. (cos + I*sin) is almost twice slower than cexp, while sincos (GNU extension)
	 * is slightly faster (3.52 - 2.39 - 2.29 for matvec in test sparse runs, where about 1.23 is for non-exp part -
	 * median values over 10 runs). So we prefer to use standard cexp.
	 * When using table (below) the corresponding timing is 1.70.
	 */
	return cexp(I*arg);
#else
	return imExpTable(arg);
#endif
}

//======================================================================================================================

static inline doublecomplex imExp(const doublecomplex arg)
/* exponent of complex argument Exp(i*arg), arg = a + i*b
 * !!! should not be used in parameter parsing (table is initialized in VariablesInterconnect())
 */
{
	if(cimag(arg)==0) return imExpTable(creal(arg)); //this case is needed to make things faster for real argument
	else return imExpTable(creal(arg))*exp(-cimag(arg));
}

static inline doublecomplex imExpM1(const double arg)
/* exp(i*arg) - 1 (should be used for small argument to avoid precision loss
 * We employ special code only for small arguments, ignoring the case when arg is close to 2piN. The latter can,
 * in principle, be handled by preliminary range reduction as in imExpTable. We do not implement it here, because such
 * case is a "coincidence" - may happen for a single voxel (or a plane of voxels), while the case of small arg may
 * happen for all voxels. In the latter case loss of precision affects all computed quantities.
 * The used expression through tan(arg/2) follows from general expression for cexpm1 below.
 */
{
	if (fabs(arg)<1) {
		double t=tan(0.5*arg);
		return -2*t/(I+t); // Alternatively, (I-t)*2t/(1+t^2), but we leave the optimization to compiler
	}
	else return imExp(arg)-1;
}

//======================================================================================================================

static inline int sgn(double v){
  return (v > 0) - (v < 0);
}

/*
static inline doublecomplex cсSinh(doublecomplex arg){
	return (cexp(arg)-cexp(-arg))/2;
}

static inline doublecomplex cсCosh(doublecomplex arg){
	return (cexp(arg)+cexp(-arg))/2;
}

static inline doublecomplex cсCoth(doublecomplex arg){
	return (cexp(arg)+cexp(-arg))/(cexp(arg)-cexp(-arg));
}
*/

static inline doublecomplex ArcCoth(const doublecomplex arg){
	return .5*clog((arg+1)/(arg-1));
}

static inline doublecomplex incJintegval(const doublecomplex w, const int nun2m){
  if(nun2m == 0){
	  return w;
  }
  else{
	  return csinh(nun2m*w)/nun2m;
  }
}

static inline double fac(int n)
{
	unsigned int i;
	double result = 1;
	for(i = 2; i<=n; i++){
		result *= i;
		//printf("%u %lf\n", i, result);
	}
	return result;
}

static inline doublecomplex incJs(const int nu, const doublecomplex arg, const doublecomplex omega, const double tol)
/* incomplete Bessel function J of complex argument incJ(order, arg, omega, tolerance) through series calculation for small arguments
 */
{
	doublecomplex Sn, Sm;
	doublecomplex incn0,incn,incm;
	int n,m;
	double fc;
	Sn = 0;
	incn0 = incJintegval(omega, nu);
	incn = incn0;
	n = 0;
	//printf("---Calculating incJ0\n");
	//printf("tol %e\n",tol);
	//printf("abs(creal(incn0))*tol %e\n",fabs(creal(incn0))*tol);
	//printf("abs(cimag(incn0))*tol %e\n",fabs(cimag(incn0))*tol);
	//printf("incn0\t=\t"CFORM"\n",REIM(incn0));
	while(fabs(creal(incn)) >= fabs(creal(incn0))*tol || fabs(cimag(incn)) >= fabs(cimag(incn0))*tol){
		Sm = 0;

		for(m = 0; m <= n; m++){
			//fc = 1/(fac(m)*fac(n - m));
			incm = incJintegval(omega, (nu + n - 2*m))/(fac(m)*fac(n - m));
			//printf("%i %i integval "CFORM" fc=%e incm "CFORM"\n",n,m,REIM(incJintegval(omega, (nu + n - 2*m))),fc,REIM(incm));
			Sm += incm;
		}
		incn = cpow((-arg/2), n)*Sm;
		Sn += incn;
		//printf("n %i, incn = "CFORM", Sm = "CFORM"\n",n,REIM(incn),REIM(Sm));

		n++;
		//if(n>=15) return Sn;
	}
	return Sn;
}

static inline doublecomplex incKpos(const int nu, const doublecomplex arg, const doublecomplex Omega)
{
	doublecomplex b0 = ccosh(nu*Omega)/csinh(Omega);
	doublecomplex b1 = (nu*csinh(nu*Omega) - ccosh(nu*Omega)*(ccosh(Omega)/csinh(Omega)))/(csinh(Omega)*csinh(Omega));
	doublecomplex b2 = (nu*nu*ccosh(nu*Omega) - 3*nu*csinh(nu*Omega)*(ccosh(Omega)/csinh(Omega)) + (2 + 3/(csinh(Omega)*csinh(Omega)))*ccosh(nu*Omega))/(2*csinh(Omega)*csinh(Omega)*csinh(Omega));

	return cexp(-arg*ccosh(Omega))*((1*b0/arg) + (1*b1/(arg*arg)) + (2*b2/(arg*arg*arg)));
}

static inline doublecomplex incJl(const int nu, const doublecomplex arg, const doublecomplex Omega)
{
	doublecomplex t1,t2,t3,t4,t5,t6,K0,K1;
	cik01_(&arg, &t1, &t2, &t3, &t4, &K0, &t5, &K1, &t6);
	if(nu == 0){
		return sgn(creal(Omega))*K0 - incKpos(nu,arg,Omega);
	}
	else if(nu == 1){
		return sgn(creal(Omega))*K1 - incKpos(nu,arg,Omega);
	}
	return 0;
}

static inline doublecomplex incJ(const int nu, const doublecomplex arg, const doublecomplex Omega, const double tol)
/* incomplete Bessel function J of complex argument incJ(order, arg, omega, tolerance)
 */
{
	//return incJs(nu,arg,Omega,tol);
	if(cabs(arg) < 20){
		//print("KAWABUNGA!!");
		return incJs(nu,arg,Omega,tol);
	}
	else
	{
		return incJl(nu,arg,Omega);
	}
}

//======================================================================================================================
/* Generalized exponential integral and incomplete modified Bessel functions K_nu(zeta,eta) of the substrate-aware
 * electron field. Equation numbers refer to the Supporting Information of the substrate EELS/CL paper (S-equations).
 */
#define EULER_GAMMA 0.57721566490153286061

static inline doublecomplex cExpIntE(const int n, const doublecomplex x)
/* generalized exponential integral E_n(x) = int_1^inf exp(-x t)/t^n dt for integer n>=0 and Re(x)>=0, x!=0;
 * power series for |x|<=1 and modified Lentz continued fraction otherwise (the standard algorithm, cf. DLMF 8.19.8,
 * 8.19.17). For n<=0 the closed form is used (eq. S38 for n=0).
 */
{
	const double eps=1e-16;
	const int maxit=5000;
	int i,ii;
	if (n==0) return cexp(-x)/x;
	if (cabs(x)<=1) {
		const int nm1=n-1;
		doublecomplex ans,fact=1,del,psi;
		ans = (nm1!=0) ? 1.0/nm1 : -clog(x)-EULER_GAMMA;
		for (i=1;i<=maxit;i++) {
			fact *= -x/i;
			if (i!=nm1) del = -fact/(i-nm1);
			else {
				psi = -EULER_GAMMA;
				for (ii=1;ii<=nm1;ii++) psi += 1.0/ii;
				del = fact*(-clog(x)+psi);
			}
			ans += del;
			if (cabs(del)<cabs(ans)*eps) return ans;
		}
		return ans;
	}
	else {
		const double tiny=1e-300;
		doublecomplex b=x+n,c=1/tiny,d=1/b,h=d,an,del;
		for (i=1;i<=maxit;i++) {
			an = -i*(double)(n-1+i);
			b += 2;
			d = an*d+b;
			if (cabs(d)<tiny) d=tiny;
			c = b+an/c;
			if (cabs(c)<tiny) c=tiny;
			d = 1/d;
			del = c*d;
			h *= del;
			if (cabs(del-1)<eps) break;
		}
		return h*cexp(-x);
	}
}

static inline doublecomplex incKS33(const int nu, const doublecomplex zeta, const doublecomplex eta, const double tol)
/* K_nu(zeta,eta) by the series in generalized exponential integrals, eq. (S33); efficient for small |zeta exp(-eta)|,
 * i.e. for large Re(eta) - exactly the case of grazing beams, where eq. (S22) overflows
 */
{
	const int maxit=500;
	const doublecomplex x=zeta*cexp(eta)/2, y=zeta*cexp(-eta)/2, ep=cexp(nu*eta), em=cexp(-nu*eta);
	doublecomplex sum=0,fn=1,term; // fn = (-y)^n/n!
	int n;
	for (n=0;n<maxit;n++) {
		if (n>0) fn *= -y/n;
		term = fn/2*(ep*cExpIntE(n+1-nu,x) + em*cExpIntE(n+1+nu,x));
		sum += term;
		if (n>=1 && cabs(term)<=tol*cabs(sum)) break;
	}
	return sum;
}

static inline doublecomplex incK(const int nu, const doublecomplex zeta, const doublecomplex eta, const double tol)
/* incomplete modified Bessel function K_nu(zeta,eta) = K_nu(zeta) - J(zeta,nu,eta), eq. (S17), for nu=0,1:
 * Re(eta)<0 is reduced to Re(eta)>0 by eq. (S18); |zeta|>=20 - large-argument expansion (S25), with n! b_n;
 * small |zeta exp(eta)| - series (S22); otherwise - series (S33).
 */
{
	doublecomplex t1,t2,t3,t4,t5,t6,K0,K1,Kfull;
	if (creal(eta)<0) {
		cik01_((doublecomplex *)&zeta,&t1,&t2,&t3,&t4,&K0,&t5,&K1,&t6);
		Kfull = (nu==0) ? K0 : K1;
		return 2*Kfull - incK(nu,zeta,-eta,tol);
	}
	if (cabs(zeta)>=20) return incKpos(nu,zeta,eta);
	if (cabs(zeta*cexp(eta))<=2) {
		cik01_((doublecomplex *)&zeta,&t1,&t2,&t3,&t4,&K0,&t5,&K1,&t6);
		Kfull = (nu==0) ? K0 : K1;
		return Kfull - incJs(nu,zeta,eta,tol);
	}
	return incKS33(nu,zeta,eta,tol);
}

//======================================================================================================================

static inline doublecomplex cExpM1(const doublecomplex a)
/* Complex analogue of expm1 function ( exp(a) - 1 ), should be used for small arguments to avoid precision loss
 * The algorithm is a simplified version of the one published in Section 17.7 of Beebe N.H.F., The Mathematical-Function
 * Computation Handbook: Programming Using the MathCW Portable Software Library. Springer; 2017.
 *  */
{
	if (fabs(creal(a))+fabs(cimag(a))<1) { // uses faster L1-norm instead of L2-norm
		doublecomplex t=ctanh(0.5*a);
		return 2*t/(1-t);
	}
	else return cexp(a)-1;
}

//======================================================================================================================
// operations on complex vectors

static inline void cvInit(doublecomplex a[static 3])
// set complex vector[3] to zero; a=0
{
	a[0] = 0;
	a[1] = 0;
	a[2] = 0;
}

//======================================================================================================================

static inline void vConj(const doublecomplex a[static 3],doublecomplex b[static 3])
// complex conjugate of vector[3]; b=a*
{
	b[0] = conj(a[0]);
	b[1] = conj(a[1]);
	b[2] = conj(a[2]);
}

//======================================================================================================================

static inline void vReal(const doublecomplex a[static 3],double b[static 3])
// takes real part of the complex vector; b=Re(a)
{
	b[0]=creal(a[0]);
	b[1]=creal(a[1]);
	b[2]=creal(a[2]);
}

//======================================================================================================================

static inline void cvBuildRe(const double a[static 3],doublecomplex b[static 3])
// builds complex vector from real part; b=a + i*0
{
	b[0]=a[0];
	b[1]=a[1];
	b[2]=a[2];
}

//======================================================================================================================

static inline void vInvRefl_cr(const double a[static 3],doublecomplex b[static 3])
/* reflects real vector with respect to the xy-plane and then inverts it, equivalent to reflection about the z-axis
 * result is stored into the complex vector
 */
{
	b[0]=-a[0];
	b[1]=-a[1];
	b[2]=a[2];
}

//======================================================================================================================

static inline void crCrossProd(const double a[static 3],const doublecomplex b[static 3],doublecomplex c[static 3])
// cross product of real and complex vector; c = a x b; !!! vectors must not alias
{
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}

//======================================================================================================================

static inline void cvMultScal(const double a,const doublecomplex b[static 3],doublecomplex c[static 3])
// multiplication of vector[3] by real scalar; c=ab;
{
	c[0] = a*b[0];
	c[1] = a*b[1];
	c[2] = a*b[2];
}

//======================================================================================================================

static inline void cvMultScal_RVec(const doublecomplex a,const double b[static 3],doublecomplex c[static 3])
// complex scalar - real vector[3] multiplication; c=ab
{
	c[0] = a*b[0];
	c[1] = a*b[1];
	c[2] = a*b[2];
}

//======================================================================================================================

static inline void cvMultScal_cmplx(const doublecomplex a,const doublecomplex b[static 3],doublecomplex c[static 3])
// multiplication of vector[3] by complex scalar; c=ab
{
	c[0] = a*b[0];
	c[1] = a*b[1];
	c[2] = a*b[2];
}

//======================================================================================================================

static inline void cvMult(const doublecomplex a[static 3],const doublecomplex b[static 3],doublecomplex c[static 3])
// multiplication of real vector by complex vector (by elements); c[i]=a[i]*b[i];
{
	c[0]=a[0]*b[0];
	c[1]=a[1]*b[1];
	c[2]=a[2]*b[2];
}

//======================================================================================================================

static inline double cvNorm2(const doublecomplex a[static 3])
// square of the norm of a complex vector[3]
{
	return cAbs2(a[0]) + cAbs2(a[1]) + cAbs2(a[2]);
}


//======================================================================================================================

static inline doublecomplex cDotProd(const doublecomplex a[static 3],const doublecomplex b[static 3])
// conjugate dot product of two complex vector[3]; c=a.b = a[0]*b*[0]+...+a[2]*b*[2]; for one vector use cvNorm2
{
	return a[0]*conj(b[0]) + a[1]*conj(b[1]) + a[2]*conj(b[2]);
}

//======================================================================================================================

static inline double cDotProd_Im(const doublecomplex a[static 3],const doublecomplex b[static 3])
/* imaginary part of dot product of two complex vector[3]; c=Im(a.b)
 * It is not clear whether this way is faster than cimag(cDotProd)
 */
{
	return ( cimag(a[0])*creal(b[0]) - creal(a[0])*cimag(b[0])
	       + cimag(a[1])*creal(b[1]) - creal(a[1])*cimag(b[1])
	       + cimag(a[2])*creal(b[2]) - creal(a[2])*cimag(b[2]) );
}

//======================================================================================================================

static inline doublecomplex cDotProd_conj(const doublecomplex a[static 3],const doublecomplex b[static 3])
// dot product of two complex vector[3] without conjugation; a.(b*) = a[0]*b[0]+...+a[2]*b[2]
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

//======================================================================================================================

static inline void cvAdd(const doublecomplex a[static 3],const doublecomplex b[static 3],doublecomplex c[static 3])
// add two complex vector[3]; c=a+b;
{
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
}

//======================================================================================================================

static inline void cvSubtr(const doublecomplex a[static 3],const doublecomplex b[static 3],doublecomplex c[static 3])
// add two complex vector[3]; c=a-b;
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}

//======================================================================================================================

static inline void cvInvSign(doublecomplex a[static 3])
// inverts the sign in the complex double vector[3]
{
	a[0]=-a[0];
	a[1]=-a[1];
	a[2]=-a[2];
}

//======================================================================================================================

static inline void cvAdd2Self(doublecomplex a[static 3],const doublecomplex b[static 3],const doublecomplex c[static 3])
// increment one complex vector[3] by sum of other two; a+=b+c
{
	a[0] += b[0] + c[0];
	a[1] += b[1] + c[1];
	a[2] += b[2] + c[2];
}

//======================================================================================================================

static inline doublecomplex crDotProd(const doublecomplex a[static 3],const double b[static 3])
// dot product of complex and real vectors[3]; a.b
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

//======================================================================================================================

static inline double crDotProd_Re(const doublecomplex a[static 3],const double b[static 3])
// real part of dot product of complex and real vectors[3]; Re(a.b)
{
	return creal(a[0])*b[0] + creal(a[1])*b[1] + creal(a[2])*b[2];
}

//======================================================================================================================

static inline void cvLinComb1(const doublecomplex a[static 3],const doublecomplex b[static 3],const double c1,
	doublecomplex c[static 3])
// linear combination of complex vectors[3] with real coefficient; second coefficient is unity; c=c1*a+b
{
	c[0] = c1*a[0] + b[0];
	c[1] = c1*a[1] + b[1];
	c[2] = c1*a[2] + b[2];
}

//======================================================================================================================

static inline void cvLinComb1_cmplx(doublecomplex a[static 3],doublecomplex b[static 3],const doublecomplex c1,
	doublecomplex c[static 3])
// linear combination of complex vectors[3] with complex coefficients; second coefficient is unity; c=c1*a+b
{
	c[0] = c1*a[0] + b[0];
	c[1] = c1*a[1] + b[1];
	c[2] = c1*a[2] + b[2];
}

//======================================================================================================================

static inline void cSymMatrVec(const doublecomplex matr[static 6],const doublecomplex vec[static 3],
	doublecomplex res[static 3])
// multiplication of complex symmetric matrix[6] by complex vec[3]; res=matr.vec; !!! vec and res must not alias
{
	res[0] = matr[0]*vec[0] + matr[1]*vec[1] + matr[2]*vec[2];
	res[1] = matr[1]*vec[0] + matr[3]*vec[1] + matr[4]*vec[2];
	res[2] = matr[2]*vec[0] + matr[4]*vec[1] + matr[5]*vec[2];
}

//======================================================================================================================

static inline void cSymMatrVecReal(const doublecomplex matr[static 6],const double vec[static 3],
	doublecomplex res[static 3])
// multiplication of complex symmetric matrix[6] by  vec[3]; res=matr.vec; !!! vec and res must not alias
{
	res[0] = matr[0]*vec[0] + matr[1]*vec[1] + matr[2]*vec[2];
	res[1] = matr[1]*vec[0] + matr[3]*vec[1] + matr[4]*vec[2];
	res[2] = matr[2]*vec[0] + matr[4]*vec[1] + matr[5]*vec[2];
}

//======================================================================================================================

static inline void cReflMatrVec(const doublecomplex matr[static 6],const doublecomplex vec[static 3],
	doublecomplex res[static 3])
/* multiplication of matrix[6] by complex vec[3]; res=matr.vec; passed components are the same as for symmetric matrix:
 * 11,12,13,22,23,33, but the matrix has the following symmetry - M21=M12, M31=-M13, M32=-M23
 * !!! vec and res must not alias
 */
{
	res[0] = matr[0]*vec[0] + matr[1]*vec[1] + matr[2]*vec[2];
	res[1] = matr[1]*vec[0] + matr[3]*vec[1] + matr[4]*vec[2];
	res[2] = - matr[2]*vec[0] - matr[4]*vec[1] + matr[5]*vec[2];
}

//======================================================================================================================

static inline void cReflMatrVecReal(const doublecomplex matr[static 6],const double vec[static 3],
	doublecomplex res[static 3])
/* multiplication of matrix[6] by complex vec[3]; res=matr.vec; passed components are the same as for symmetric matrix:
 * 11,12,13,22,23,33, but the matrix has the following symmetry - M21=M12, M31=-M13, M32=-M23
 * !!! vec and res must not alias
 */
{
	res[0] = matr[0]*vec[0] + matr[1]*vec[1] + matr[2]*vec[2];
	res[1] = matr[1]*vec[0] + matr[3]*vec[1] + matr[4]*vec[2];
	res[2] = - matr[2]*vec[0] - matr[4]*vec[1] + matr[5]*vec[2];
}

//======================================================================================================================
// operations on real vectors

static inline void vInit(double a[static 3])
// set real vector[3] to zero; a=0
{
	a[0] = 0;
	a[1] = 0;
	a[2] = 0;
}

//======================================================================================================================

static inline void vCopy(const double a[static 3],double b[static 3])
// copies one vector into another; b=a
{
	// can be rewritten through memcpy, but compiler should be able to do it itself if needed
	b[0]=a[0];
	b[1]=a[1];
	b[2]=a[2];
}

//======================================================================================================================

static inline void vAdd(const double a[static 3],const double b[static 3],double c[static 3])
// adds two real vectors; c=a+b
{
	c[0]=a[0]+b[0];
	c[1]=a[1]+b[1];
	c[2]=a[2]+b[2];
}

//======================================================================================================================

static inline void vSubtr(const double a[static 3],const double b[static 3],double c[static 3])
// subtracts two real vectors; c=a-b
{
	c[0]=a[0]-b[0];
	c[1]=a[1]-b[1];
	c[2]=a[2]-b[2];
}

//======================================================================================================================

static inline void vInvSign(double a[static 3])
// inverts the sign in the double vector[3]
{
	a[0]=-a[0];
	a[1]=-a[1];
	a[2]=-a[2];
}

//======================================================================================================================

static inline void vRefl(const double inc[static 3],double ref[static 3])
// reflects the incident vector 'inc' with respect to the xy-plane (inverts z-component)
{
	ref[0]=inc[0];
	ref[1]=inc[1];
	ref[2]=-inc[2];
}

//======================================================================================================================

static inline void vMultScal(const double a,const double b[static 3],double c[static 3])
// multiplication of real vector by scalar; c=a*b;
{
	c[0]=a*b[0];
	c[1]=a*b[1];
	c[2]=a*b[2];
}

//======================================================================================================================

static inline void vMult(const double a[static 3],const double b[static 3],double c[static 3])
// multiplication of two vectors (by elements); c[i]=a[i]*b[i];
{
	c[0]=a[0]*b[0];
	c[1]=a[1]*b[1];
	c[2]=a[2]*b[2];
}

//======================================================================================================================

static inline bool vAlongZ(const double a[static 3])
// a robust (with respect to round-off errors) way to test that vector is along the z-axis (+ or -)
{
	return fabs(a[0])<ROUND_ERR && fabs(a[1])<ROUND_ERR;
}

//======================================================================================================================

static inline double DotProd(const double a[static 3],const double b[static 3])
// dot product of two real vectors[3]; use DotProd(x,x) to get squared norm
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

//======================================================================================================================

static inline double DotProdSquare(const double a[static 3],const double b[static 3])
// dot product of element-wise squares of two real vectors[3]
{
	return a[0]*a[0]*b[0]*b[0] + a[1]*a[1]*b[1]*b[1] + a[2]*a[2]*b[2]*b[2];
}

//======================================================================================================================

static inline double vNorm(const double a[static 3])
// norm of a real vector[3]
{
	return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

//======================================================================================================================

static inline void CrossProd(const double a[static 3],const double b[static 3],double c[static 3])
// cross product of two real vectors; c = a x b; !!! vectors must not alias
{
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}

//======================================================================================================================

static inline void vNormalize(double a[static 3])
// normalize real vector to have unit norm
{
	double c;
	c=1./sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	a[0]*=c;
	a[1]*=c;
	a[2]*=c;
}

//======================================================================================================================

static inline void LinComb(const double a[static 3],const double b[static 3],const double c1,const double c2,
	double c[static 3])
// linear combination of real vectors[3]; c=c1*a+c2*b
{
	c[0]=c1*a[0]+c2*b[0];
	c[1]=c1*a[1]+c2*b[1];
	c[2]=c1*a[2]+c2*b[2];
}

//======================================================================================================================

static inline void OuterSym(const double a[static 3],double matr[static 6])
// outer product of real vector a with itself, stored in symmetric matrix matr
{
	matr[0] = a[0]*a[0];
	matr[1] = a[0]*a[1];
	matr[2] = a[0]*a[2];
	matr[3] = a[1]*a[1];
	matr[4] = a[1]*a[2];
	matr[5] = a[2]*a[2];
}

//======================================================================================================================

static inline double TrSym(const double a[static 6])
// trace of a symmetric matrix stored as a vector of size 6
{
	return (a[0]+a[2]+a[5]);
}

//======================================================================================================================

static inline double QuadForm(const double matr[static 6],const double vec[static 3])
// value of a quadratic form matr (symmetric matrix stored as a vector of size 6) over a vector vec;
{
	return ( vec[0]*vec[0]*matr[0] + vec[1]*vec[1]*matr[2] + vec[2]*vec[2]*matr[5]
	       + 2*(vec[0]*vec[1]*matr[1] + vec[0]*vec[2]*matr[3] + vec[1]*vec[2]*matr[4]) );
}

//======================================================================================================================

static inline void MatrVec(double matr[static 3][3],const double vec[static 3],double res[static 3])
// multiplication of matrix[3][3] by vec[3] (all real); res=matr.vec;
{
	res[0]=matr[0][0]*vec[0]+matr[0][1]*vec[1]+matr[0][2]*vec[2];
	res[1]=matr[1][0]*vec[0]+matr[1][1]*vec[1]+matr[1][2]*vec[2];
	res[2]=matr[2][0]*vec[0]+matr[2][1]*vec[1]+matr[2][2]*vec[2];
}

//======================================================================================================================

static inline void MatrColumn(double matr[static 3][3],const int ind,double vec[static 3])
// get ind's column of matrix[3][3] and store it into vec[3] (all real, ind starts from zero); vec=matr[.][ind];
{
	vec[0]=matr[0][ind];
	vec[1]=matr[1][ind];
	vec[2]=matr[2][ind];
}

//======================================================================================================================

static inline void Permutate(double vec[static 3],const int ord[static 3])
// permutate double vector vec using permutation ord
{
	double buf[3];

	memcpy(buf,vec,3*sizeof(double));
	vec[0]=buf[ord[0]];
	vec[1]=buf[ord[1]];
	vec[2]=buf[ord[2]];
}

//======================================================================================================================

static inline void Permutate_i(int vec[static 3],const int ord[static 3])
// permutate int vector vec using permutation ord
{
	int buf[3];

	memcpy(buf,vec,3*sizeof(int));
	vec[0]=buf[ord[0]];
	vec[1]=buf[ord[1]];
	vec[2]=buf[ord[2]];
}

//======================================================================================================================
// Auxiliary functions

static inline double Deg2Rad(const double deg)
// transforms angle in degrees to radians
{
	return (PI_OVER_180*deg);
}

//======================================================================================================================

static inline double Rad2Deg(const double rad)
// transforms angle in radians to degrees
{
	return (INV_PI_180*rad);
}

//======================================================================================================================

static inline bool TestBelowDeg(const double deg)
/* tests if the direction is below the substrate using the degree theta in degrees;
 * if unsure (within rounded error) returns false (above)
 */
{
	return fabs(fmod(fabs(deg),360)-180) < 90*(1-ROUND_ERR);
}

//======================================================================================================================
// functions used for substrate

//======================================================================================================================

static inline doublecomplex FresnelRS(const doublecomplex ki,const doublecomplex kt)
/* reflection coefficient for s-polarized wave (E perpendicular to the main plane),
 * ki,kt are normal (positive) components of incident and transmitted wavevector (with arbitrary mutual scaling)
 */
{
	return (ki-kt)/(ki+kt);
}

//======================================================================================================================

static inline doublecomplex FresnelTS(const doublecomplex ki,const doublecomplex kt)
/* transmission coefficient for s-polarized wave (E perpendicular to the main plane),
 * ki,kt are normal (positive) components of incident and transmitted wavevector (with arbitrary mutual scaling)
 */
{
	return 2*ki/(ki+kt);
}

//======================================================================================================================

static inline doublecomplex FresnelRP(const doublecomplex ki,const doublecomplex kt,const doublecomplex mr)
/* reflection coefficient for p-polarized wave (E parallel to the main plane),
 * ki,kt are normal (positive) components of incident and transmitted wavevector (with arbitrary mutual scaling)
 * mr is the ratio of refractive indices (mt/mi)
 */
{
	return (mr*mr*ki-kt)/(mr*mr*ki+kt);
}

//======================================================================================================================

static inline doublecomplex FresnelTP(const doublecomplex ki,const doublecomplex kt,const doublecomplex mr)
/* transmission coefficient for p-polarized wave (E parallel to the main plane),
 * ki,kt are normal (positive) components of incident and transmitted wavevector (with arbitrary mutual scaling)
 * mr is the ratio of refractive indices (mt/mi)
 */
{
	return 2*mr*ki/(mr*mr*ki+kt);
}

#ifdef USE_SSE3

//======================================================================================================================

static inline __m128d cmul(__m128d a,__m128d b)
// complex multiplication
{
	__m128d t;
	t = _mm_movedup_pd(a);
	a = _mm_shuffle_pd(a,a,3);
	t = _mm_mul_pd(b,t);
	b = _mm_shuffle_pd(b,b,1);
	b = _mm_mul_pd(a,b);
	return _mm_addsub_pd(t,b);
}

//======================================================================================================================

static inline __m128d cadd(__m128d a,__m128d b)
// complex addition
{
	return _mm_add_pd(a,b);
}

#endif // USE_SSE3

#endif // __cmplx_h

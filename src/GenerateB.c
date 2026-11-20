/* Generates the incident beam
 *
 * Lminus beam is based on: G. Gouesbet, B. Maheu, G. Grehan, "Light scattering from a sphere arbitrary located
 * in a Gaussian beam, using a Bromwhich formulation", J.Opt.Soc.Am.A 5,1427-1443 (1988).
 * Eq.(22) - complex conjugate.
 *
 * Davis beam is based on: L. W. Davis, "Theory of electromagnetic beams," Phys.Rev.A 19, 1177-1179 (1979).
 * Eqs.(15a),(15b) - complex conjugate; in (15a) "Q" changed to "Q^2" (typo).
 *
 * Barton beam is based on: J. P. Barton and D. R. Alexander, "Fifth-order corrected electromagnetic-field components
 * for a fundamental Gaussian-beam," J.Appl.Phys. 66,2800-2802 (1989).
 * Eqs.(25)-(28) - complex conjugate.
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
#include "const.h" // keep this first
// project headers
#include "cmplx.h"
#include "comm.h"
#include "io.h"
#include "interaction.h"
#include "param.h"
#include "vars.h"
// system headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

// SEMI-GLOBAL VARIABLES

// defined and initialized in param.c
extern const int beam_Npars;
extern const double beam_pars[];
extern const char *beam_fnameY;
extern const char *beam_fnameX;
extern const opt_index opt_beam;

extern void cik01_(doublecomplex *z, doublecomplex *cbi0, doublecomplex *cdi0, doublecomplex *cbi1, doublecomplex *cdi1, doublecomplex *cbk0, doublecomplex *cdk0, doublecomplex *cbk1, doublecomplex *cdk1);

// used in CalculateE.c
double C0dipole,C0dipole_refl; // inherent cross sections of exciting dipole (in free space and addition due to surface)
int vorticity;                 // Vorticity of vortex beams (besN for Bessel beams)
// used in param.c
const char *beam_descr; // string for log file with beam parameters
/* Propagation (phase) directions of secondary incident beams above (refl) and below (tran) the surface (unit vectors)
 * When msub is complex, one of this doesn't tell the complete story, since the corresponding wave is inhomogeneous,
 * given by the complex wavenumber ktVec
 */
double prIncRefl[3],prIncTran[3];

// LOCAL VARIABLES
static double s,s2;            // beam confinement factor and its square
static double scale_x,scale_z; // multipliers for scaling coordinates
static doublecomplex ki,kt;    // abs of normal components of k_inc/k0, and ktran/k0
static doublecomplex ktVec[3]; // k_tran/k0
static double p0;              // amplitude of the incident dipole moment
static doublecomplex beta1, gamma1_inv;// v*m/c, 1/gamma
static doublecomplex e_pref1; // prefactor of the field of the electron
static double e_w_v;   // prefactor in an argument of a phase exponent in the incident field of the electron
static doublecomplex e_w_gv1, e_w_gv2;  // prefactor in an argument of the BesselK, incBesselJ in the incident field of the electron
static doublecomplex eps1, eps2; //dielectric permittivity of the media above and below the plane surface
static doublecomplex vgw, vgw2; //Some constants for the field of an electron
static const bool tilted_regularization_enabled = false; // toggle projected-grid regularization for tilted electron beams
#ifndef NO_FORTRAN
static int besN;                  // Bessel beam order
static double besAlpha;           // half-cone angle (in radians)
static doublecomplex besKt,besKz; // wave-vector components (transverse and longitudinal)
static doublecomplex besM[4];     // Matrix M defining the generalized Bessel beam
#endif
/* TO ADD NEW BEAM
 * Add here all internal variables (beam parameters), which you initialize in InitBeam() and use in GenerateB()
 * afterwards. If you need local, intermediate variables, put them into the beginning of the corresponding function.
 * Add descriptive comments, use 'static'.
 */

// EXTERNAL FUNCTIONS

#ifndef NO_FORTRAN
void bjndd_(const int *n,const double *x,double *bj,double *dj,double *fj);
#endif



static void ElectronFullLineField(const double obs[static 3],const double center[static 3],const double dir[static 3],
	const bool auxiliary,doublecomplex out[static 3])
// corrected field of the full infinite electron trajectory in a homogeneous medium;
// this helper is consistent with the perpendicular thesis medium field in the limit dir->-ez.
{
	double r1[3],r1par[3],r1per[3];
	double ro,temp,s;
	doublecomplex arg,t1,t2,t3,t4,t5,t6,K0,K1,phase;
	doublecomplex v1[3],v2[3],v3[3];

	LinComb(obs,center,1,-1,r1);
	temp = DotProd(r1,dir);
	vMultScal(temp,dir,r1par);
	vSubtr(r1,r1par,r1per);
	ro = vNorm(r1per);
	if (ro != 0) vNormalize(r1per);
	if (ro < dsX/2) ro = dsX/2;
	s = -temp;
	arg = e_w_gv1*ro;
	cik01_(&arg, &t1, &t2, &t3, &t4, &K0, &t5, &K1, &t6);
	phase = imExp(-e_w_v*s);
	if (auxiliary) {
		phase = conj(phase);
		cvMultScal_RVec(phase*K1,r1per,v1);
		cvMultScal_RVec((I)*gamma1_inv*phase*K0,dir,v2);
	}
	else {
		cvMultScal_RVec(-phase*K1,r1per,v1);
		cvMultScal_RVec((I)*gamma1_inv*phase*K0,dir,v2);
	}
	cvAdd(v1,v2,v3);
	cvMultScal_cmplx(2*e_pref1,v3,out);
}

static void ElectronFullLineFieldNoFloor(const double obs[static 3],const double center[static 3],const double dir[static 3],
	const bool auxiliary,doublecomplex out[static 3])
// corrected exact full infinite trajectory field with only a tiny numerical safeguard; no dX/2 regularization
{
	double r1[3],r1par[3],r1per[3];
	double ro,temp,s;
	doublecomplex arg,t1,t2,t3,t4,t5,t6,K0,K1,phase;
	doublecomplex v1[3],v2[3],v3[3];
	const double tiny = 1e-12*dsX;

	LinComb(obs,center,1,-1,r1);
	temp = DotProd(r1,dir);
	vMultScal(temp,dir,r1par);
	vSubtr(r1,r1par,r1per);
	ro = vNorm(r1per);
	if (ro != 0) vNormalize(r1per);
	if (ro < tiny) ro = tiny;
	s = -temp;
	arg = e_w_gv1*ro;
	cik01_(&arg, &t1, &t2, &t3, &t4, &K0, &t5, &K1, &t6);
	phase = imExp(-e_w_v*s);
	if (auxiliary) {
		phase = conj(phase);
		cvMultScal_RVec(phase*K1,r1per,v1);
		cvMultScal_RVec((I)*gamma1_inv*phase*K0,dir,v2);
	}
	else {
		cvMultScal_RVec(-phase*K1,r1per,v1);
		cvMultScal_RVec((I)*gamma1_inv*phase*K0,dir,v2);
	}
	cvAdd(v1,v2,v3);
	cvMultScal_cmplx(2*e_pref1,v3,out);
}

static void ElectronBranchPlusRealNoFloor(const double obs[static 3],const double cross[static 3],const double dir[static 3],
	const double s0,const bool auxiliary,doublecomplex out[static 3])
// corrected lower real branch: E_+^(r) for incident (auxiliary=false) and auxiliary (auxiliary=true)
{
	double r1[3],r1par[3],r1per[3],rhat[3]={0,0,0};
	double rho,s,temp,R;
	doublecomplex arg,t1,t2,t3,t4,t5,t6,K0,K1,Omega,incJ0,incJ1,F,K1branch;
	doublecomplex phase0,phase_line,expikr,A,B;
	const double k = WaveNum0;
	const double alpha = e_w_v;
	const double tol = 1e-8;
	const double tiny = 1e-12*dsX;
	doublecomplex core_pref,ep_pref;

	LinComb(obs,cross,1,-1,r1);
	temp = DotProd(r1,dir);
	vMultScal(temp,dir,r1par);
	vSubtr(r1,r1par,r1per);
	rho = vNorm(r1per);
	if (rho != 0) { vNormalize(r1per); vCopy(r1per,rhat); }
	if (rho < tiny) rho = tiny;
	s = -temp;
	R = sqrt(rho*rho+s*s);
	arg = e_w_gv1*rho;
	cik01_(&arg, &t1, &t2, &t3, &t4, &K0, &t5, &K1, &t6);
	expikr = imExp(k*R);

	if (auxiliary) {
		phase0 = imExp(-alpha*s0);
		phase_line = imExp(alpha*s);
		Omega = asinh(s/rho) - ArcCoth(beta1);
		incJ0 = incJ(0,arg,Omega,tol);
		incJ1 = incJ(1,arg,Omega,tol);
		F = phase_line*(K0-incJ0);
		K1branch = phase_line*(K1-incJ1);
		A = s/(R*R*R)*(1-I*k*R) - I*alpha/R;
		B = rho/(R*R*R)*(1-I*k*R) + I*alpha*s/(R*rho);
		core_pref = -e_pref1*phase0;
		ep_pref = (I)*e_pref1*(beta1*beta1)/gamma1_inv * phase0;
		cvMultScal_RVec(core_pref*K1branch,rhat,out);
		doublecomplex tmpv[3],tmp2[3];
		cvMultScal_RVec(core_pref*(I)*gamma1_inv*F,dir,tmpv);
		cvAdd(out,tmpv,out);
		cvMultScal_RVec(ep_pref*expikr*A/(k*k),dir,tmpv);
		cvMultScal_RVec(-ep_pref*expikr*B/(k*k),rhat,tmp2);
		cvAdd(out,tmpv,out);
		cvAdd(out,tmp2,out);
	}
	else {
		phase0 = imExp(alpha*s0);
		phase_line = imExp(-alpha*s);
		Omega = asinh(s/rho) + ArcCoth(beta1);
		incJ0 = incJ(0,arg,Omega,tol);
		incJ1 = incJ(1,arg,Omega,tol);
		F = phase_line*(K0-incJ0);
		K1branch = phase_line*(K1-incJ1);
		A = s/(R*R*R)*(1-I*k*R) + I*alpha/R;
		B = rho/(R*R*R)*(1-I*k*R) - I*alpha*s/(R*rho);
		core_pref = e_pref1*phase0;
		ep_pref = (I)*e_pref1*(beta1*beta1)/gamma1_inv * phase0;
		cvMultScal_RVec(core_pref*K1branch,rhat,out);
		doublecomplex tmpv[3],tmp2[3];
		cvMultScal_RVec(core_pref*(-I)*gamma1_inv*F,dir,tmpv);
		cvAdd(out,tmpv,out);
		cvMultScal_RVec(ep_pref*expikr*A/(k*k),dir,tmpv);
		cvMultScal_RVec(-ep_pref*expikr*B/(k*k),rhat,tmp2);
		cvAdd(out,tmpv,out);
		cvAdd(out,tmp2,out);
	}
}

static void ElectronBranchMinusImageNoFloor(const double obs[static 3],const double cross[static 3],const double dir[static 3],
	const double s0,const bool auxiliary,doublecomplex out[static 3])
// corrected lower image branch: E_-^(i) for incident (auxiliary=false) and auxiliary (auxiliary=true)
{
	double r1[3],r1par[3],r1per[3],rhat[3]={0,0,0};
	double rho,s,temp,R;
	doublecomplex arg,t1,t2,t3,t4,t5,t6,K0,K1,Omega,incJ0,incJ1,F,K1branch;
	doublecomplex phase0,phase_line,expikr,A,B;
	const double k = WaveNum0;
	const double alpha = e_w_v;
	const double tol = 1e-8;
	const double tiny = 1e-12*dsX;
	doublecomplex core_pref,ep_pref;

	LinComb(obs,cross,1,-1,r1);
	temp = DotProd(r1,dir);
	vMultScal(temp,dir,r1par);
	vSubtr(r1,r1par,r1per);
	rho = vNorm(r1per);
	if (rho != 0) { vNormalize(r1per); vCopy(r1per,rhat); }
	if (rho < tiny) rho = tiny;
	s = -temp;
	R = sqrt(rho*rho+s*s);
	arg = e_w_gv1*rho;
	cik01_(&arg, &t1, &t2, &t3, &t4, &K0, &t5, &K1, &t6);
	expikr = imExp(k*R);

	if (auxiliary) {
		phase0 = imExp(-alpha*s0);
		phase_line = imExp(alpha*s);
		// unstable direct branch would use Omega = -asinh(s/rho) + ArcCoth(beta1)
		// Evaluate it by the complement identity K_nu(-eta)=2K_nu-K_nu(eta),
		// with the stable branch eta = asinh(s/rho) - ArcCoth(beta1).
		Omega = asinh(s/rho) - ArcCoth(beta1);
		incJ0 = incJ(0,arg,Omega,tol);
		incJ1 = incJ(1,arg,Omega,tol);
		F = phase_line*(K0+incJ0);
		K1branch = phase_line*(K1+incJ1);
		A = s/(R*R*R)*(1-I*k*R) - I*alpha/R;
		B = rho/(R*R*R)*(1-I*k*R) + I*alpha*s/(R*rho);
		core_pref = -e_pref1*phase0;
		ep_pref = -(I)*e_pref1*(beta1*beta1)/gamma1_inv * phase0;
		cvMultScal_RVec(core_pref*K1branch,rhat,out);
		doublecomplex tmpv[3],tmp2[3];
		cvMultScal_RVec(core_pref*(I)*gamma1_inv*F,dir,tmpv);
		cvAdd(out,tmpv,out);
		cvMultScal_RVec(ep_pref*expikr*A/(k*k),dir,tmpv);
		cvMultScal_RVec(-ep_pref*expikr*B/(k*k),rhat,tmp2);
		cvAdd(out,tmpv,out);
		cvAdd(out,tmp2,out);
	}
	else {
		phase0 = imExp(alpha*s0);
		phase_line = imExp(-alpha*s);
		// unstable direct branch would use Omega = -asinh(s/rho) - ArcCoth(beta1)
		// Evaluate it by the complement identity K_nu(-eta)=2K_nu-K_nu(eta),
		// with the stable branch eta = asinh(s/rho) + ArcCoth(beta1).
		Omega = asinh(s/rho) + ArcCoth(beta1);
		incJ0 = incJ(0,arg,Omega,tol);
		incJ1 = incJ(1,arg,Omega,tol);
		F = phase_line*(K0+incJ0);
		K1branch = phase_line*(K1+incJ1);
		A = s/(R*R*R)*(1-I*k*R) + I*alpha/R;
		B = rho/(R*R*R)*(1-I*k*R) - I*alpha*s/(R*rho);
		core_pref = e_pref1*phase0;
		ep_pref = -(I)*e_pref1*(beta1*beta1)/gamma1_inv * phase0;
		cvMultScal_RVec(core_pref*K1branch,rhat,out);
		doublecomplex tmpv[3],tmp2[3];
		cvMultScal_RVec(core_pref*(-I)*gamma1_inv*F,dir,tmpv);
		cvAdd(out,tmpv,out);
		cvMultScal_RVec(ep_pref*expikr*A/(k*k),dir,tmpv);
		cvMultScal_RVec(-ep_pref*expikr*B/(k*k),rhat,tmp2);
		cvAdd(out,tmpv,out);
		cvAdd(out,tmp2,out);
	}
}

static bool LineIntersectsDipoleCube(const double line_point[static 3],const double dir[static 3],const double cube_center[static 3],const double half)
// infinite-line / axis-aligned cube intersection test by slab clipping
{
	double tmin = -DBL_MAX, tmax = DBL_MAX;
	for (int a=0; a<3; ++a) {
		double p = line_point[a] - cube_center[a];
		double d = dir[a];
		if (fabs(d) < 1e-15) {
			if (fabs(p) > half) return false;
		}
		else {
			double t1 = (-half - p)/d;
			double t2 = ( half - p)/d;
			if (t1 > t2) { double tmp=t1; t1=t2; t2=tmp; }
			if (t1 > tmin) tmin = t1;
			if (t2 < tmax) tmax = t2;
			if (tmin > tmax) return false;
		}
	}
	return true;
}





#define EULER_GAMMA_CONST 0.5772156649015328606

static inline doublecomplex K0_asym(const doublecomplex x)
// K0(x) = L + x^2/4 (1-L) + x^4/128 (3-2L) + O(x^6 L),  L=-ln(x/2)-gamma_E
{
	doublecomplex L = -clog(x/2.0)-EULER_GAMMA_CONST;
	doublecomplex x2 = x*x;
	return L + 0.25*x2*(1.0-L) + (x2*x2)*(3.0-2.0*L)/128.0;
}

static inline doublecomplex K1_asym(const doublecomplex x)
// K1(x) = 1/x - x/2 L - x/4 - x^3/16 L - 5 x^3/64 + O(x^5 L)
{
	doublecomplex L = -clog(x/2.0)-EULER_GAMMA_CONST;
	doublecomplex x2 = x*x;
	return 1.0/x - 0.5*x*L - 0.25*x - (x*x2)*L/16.0 - 5.0*(x*x2)/64.0;
}

static inline doublecomplex K0inc_asym(const doublecomplex x,const doublecomplex Omega)
// K_0(x,Omega)=K0(x)-J0, with J0 expanded through x^3
{
	doublecomplex L = -clog(x/2.0)-EULER_GAMMA_CONST;
	doublecomplex x2 = x*x;
	doublecomplex s1 = csinh(Omega);
	doublecomplex s2 = csinh(2.0*Omega);
	doublecomplex s3 = csinh(3.0*Omega);
	return (L - Omega)
		+ x*s1
		+ x2*( (1.0-L)/4.0 - (2.0*Omega + s2)/8.0 )
		+ (x*x2)*(s3 + 9.0*s1)/72.0;
}

static inline doublecomplex K1inc_asym(const doublecomplex x,const doublecomplex Omega)
// K_1(x,Omega)=K1(x)-J1, with J1 expanded through x^3
{
	doublecomplex L = -clog(x/2.0)-EULER_GAMMA_CONST;
	doublecomplex x2 = x*x;
	doublecomplex s1 = csinh(Omega);
	doublecomplex s2 = csinh(2.0*Omega);
	doublecomplex s4 = csinh(4.0*Omega);
	doublecomplex c2 = ccosh(Omega)*ccosh(Omega);
	return 1.0/x
		- s1
		+ x*( (2.0*Omega + s2)/4.0 - 0.5*L - 0.25 )
		- x2*( (c2 + 2.0)*s1 )/6.0
		+ (x*x2)*( (12.0*Omega + 8.0*s2 + s4)/192.0 - L/16.0 - 5.0/64.0 );
}

static void ElectronTiltedFieldPointNoFloor(const double obs[static 3],const double beam_center[static 3],
	const double cross[static 3],const double dir[static 3],const double dir_image[static 3],const double s0,const double s0_image,
	const doublecomplex refl,const bool auxiliary,doublecomplex out[static 3])
// re-implemented exact tilted field at one point with only tiny safeguards in constituent kernels; no dX/2 regularization.
// Sign convention is fixed so that the theta->0 limit reproduces the validated perpendicular implementation:
// E_tilt = E_full^(r) - r * ( E_+^(r) + E_-^(i) ).
{
	doublecomplex vfull_real[3],vplus_real[3],vminus_img[3],tmp[3];
	ElectronFullLineFieldNoFloor(obs,beam_center,dir,auxiliary,vfull_real);
	ElectronBranchPlusRealNoFloor(obs,cross,dir,s0,auxiliary,vplus_real);
	ElectronBranchMinusImageNoFloor(obs,cross,dir_image,s0_image,auxiliary,vminus_img);
	cvAdd(vplus_real,vminus_img,tmp);
	cvMultScal_cmplx(-refl,tmp,tmp);
	cvAdd(vfull_real,tmp,out);
}



static double SnapToNearestHalfInteger(const double c)
// nearest value in Z+1/2
{
	return floor(c + 0.5) - 0.5;
}


static inline void Cross3(const double a[static 3],const double b[static 3],double out[static 3])
{
	out[0] = a[1]*b[2] - a[2]*b[1];
	out[1] = a[2]*b[0] - a[0]*b[2];
	out[2] = a[0]*b[1] - a[1]*b[0];
}

static inline double NormSq3(const double a[static 3])
{
	return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
}

static bool BuildTiltedGridRegularizedPoint(const double obs[static 3],const double cross[static 3],const double dir[static 3],
	double obs_reg[static 3],double *rho_real_out)
// projected-grid regularization for arbitrary beam direction.
// For points inside the beam-hit cylinder, project the transverse offset to the nearest admissible
// face of the projected Wigner-Seitz cell of the cubic lattice. This guarantees the correct dX/2
// limit in the perpendicular/parallel cases and avoids corner-snapping to dX/sqrt(2).
{
	double r1[3],r1par[3],r1per[3],a[3][3],g[6][3],ghat[3],qreg[3],cand[3];
	double rho_real,temp,gnorm,u,dface,best_dist=DBL_MAX,best_z=-DBL_MAX;
	int j,k,m; bool found=false;
	LinComb(obs,cross,1,-1,r1);
	temp = DotProd(r1,dir);
	vMultScal(temp,dir,r1par);
	vSubtr(r1,r1par,r1per);
	rho_real = vNorm(r1per);
	if (rho_real_out) *rho_real_out = rho_real;
	if (rho_real >= dsX/2) return false;

	for (k=0;k<3;k++) {
		for (j=0;j<3;j++) a[k][j] = -dsX*dir[k]*dir[j];
		a[k][k] += dsX;
	}
	vCopy(a[0],g[0]);
	vCopy(a[1],g[1]);
	vCopy(a[2],g[2]);
	vSubtr(a[0],a[1],g[3]);
	vSubtr(a[0],a[2],g[4]);
	vSubtr(a[1],a[2],g[5]);

	for (m=0; m<6; ++m) {
		gnorm = vNorm(g[m]);
		if (gnorm < 1e-14*dsX) continue;
		vMultScal(1.0/gnorm,g[m],ghat);
		u = DotProd(r1per,ghat);
		for (int sgn=-1; sgn<=1; sgn+=2) {
			dface = sgn*0.5*gnorm - u;
			LinComb(r1per,ghat,1,dface,qreg);
			LinComb(cross,r1par,1,1,cand);
			LinComb(cand,qreg,1,1,cand);
			if (cand[2] + 1e-12 < dsX/2) continue;
			/* candidate must lie on or inside all projected Wigner-Seitz half-space boundaries */
			bool inside=true;
			for (int n=0; n<6; ++n) {
				double gn2 = NormSq3(g[n]);
				if (gn2 < 1e-28) continue;
				if (fabs(DotProd(qreg,g[n])) > 0.5*gn2 + 1e-10*dsX*dsX) { inside=false; break; }
			}
			if (!inside) continue;
			if (!found || fabs(dface) < best_dist - 1e-14 || (fabs(fabs(dface)-best_dist) < 1e-14 && cand[2] > best_z)) {
				best_dist = fabs(dface);
				best_z = cand[2];
				obs_reg[0] = cand[0];
				obs_reg[1] = cand[1];
				obs_reg[2] = cand[2];
				found = true;
			}
		}
	}
	if (found) return true;

	/* conservative fallback: move along the shortest projected lattice direction by half a step */
	{
		int best_m=-1; double best_len=DBL_MAX;
		for (m=0; m<6; ++m) {
			gnorm = vNorm(g[m]);
			if (gnorm > 1e-14*dsX && gnorm < best_len) { best_len=gnorm; best_m=m; }
		}
		if (best_m < 0) return false;
		vMultScal(1.0/best_len,g[best_m],ghat);
		for (int sgn=-1; sgn<=1; sgn+=2) {
			LinComb(r1per,ghat,1,sgn*0.5*best_len,qreg);
			LinComb(cross,r1par,1,1,cand);
			LinComb(cand,qreg,1,1,cand);
			if (cand[2] >= dsX/2 - 1e-12) { vCopy(cand,obs_reg); return true; }
		}
		vCopy(obs,obs_reg);
		if (obs_reg[2] < dsX/2) obs_reg[2] = dsX/2;
		return true;
	}
}


static void ElectronPerpendicularSurfaceField(const double obs[static 3],const double beam_center[static 3],
	const double hsub,const bool auxiliary,doublecomplex out[static 3])
// validated perpendicular -surf field used as regularized fallback inside the near-hit cylinder of the tilted branch
{
	double r1[3],ro,z,z0,R1;
	doublecomplex e_wb_gv,t1,t2,t3,t4,t5,t6,K0,K1,omp,omm,incJ0p,incJ0m,incJ1p,incJ1m;
	doublecomplex phase_z0,phase_med,expikr,med_xy,med_z,sub_pref,sub_xy,sub_z,phase_a;

	LinComb(obs,beam_center,1,-1,r1);
	ro = sqrt(r1[0]*r1[0]+r1[1]*r1[1]);
	z = obs[2] + hsub;
	if (ro < dsX/2) ro = dsX/2;
	z0 = beam_center[2] + hsub;
	e_wb_gv = e_w_gv1*ro;
	cik01_(&e_wb_gv, &t1, &t2, &t3, &t4, &K0, &t5, &K1, &t6);
	t4 = imExp(e_w_v*z);
	phase_z0 = imExp(e_w_v*z0);
	phase_med = conj(t4)*phase_z0;
	omp = asinh(z/ro)+ArcCoth(beta1);
	omm = asinh(z/ro)-ArcCoth(beta1);
	incJ0p = incJ(0,e_wb_gv,omp,1e-8);
	incJ0m = incJ(0,e_wb_gv,omm,1e-8);
	incJ1p = incJ(1,e_wb_gv,omp,1e-8);
	incJ1m = incJ(1,e_wb_gv,omm,1e-8);
	R1 = 1/sqrt(ro*ro+z*z);
	doublecomplex K1p = conj(t4)*(K1-incJ1p), K1m = t4*(K1-incJ1m), K0p = conj(t4)*(K0-incJ0p), K0m = t4*(K0-incJ0m);
	expikr = imExp(WaveNum0*sqrt(ro*ro+z*z));
	med_xy = -2*e_pref1*phase_med*K1;
	med_z = (-I)*gamma1_inv*2*e_pref1*phase_med*K0;
	sub_pref = -e_pref1*((eps1-eps2)/(eps1+eps2))*phase_z0;
	sub_xy = sub_pref*(K1p+K1m-2*z*R1/e_wb_gv*expikr);
	sub_z = sub_pref*(I*gamma1_inv*(K0p-K0m)+2*R1/e_w_gv1*expikr);
	if (!auxiliary) {
		out[0] = (r1[0]/ro)*(med_xy+sub_xy);
		out[1] = (r1[1]/ro)*(med_xy+sub_xy);
		out[2] = med_z+sub_z;
	}
	else {
		phase_a = imExp(2*e_w_v*(z-z0));
		doublecomplex Einc_sub[3]={(r1[0]/ro)*sub_xy,(r1[1]/ro)*sub_xy,sub_z};
		out[0] = -phase_a*(r1[0]/ro)*med_xy;
		out[1] = -phase_a*(r1[1]/ro)*med_xy;
		out[2] = phase_a*med_z;
		doublecomplex suba[3];
		cvMultScal_cmplx(-imExp(-2*e_w_v*z0),Einc_sub,suba);
		cvAdd(out,suba,out);
	}
}

//======================================================================================================================

void InitBeam(void)
// initialize beam; produce description string
{
	double w0; // beam width
	//CASE: B_ELECTRON
	static double e_energy;        // kinetic energy of the electron
	static double e_v;      // speed of the electron
	double q_electron = -4.803204673e-10; //electric charge of an electron, cm^(3/2)*g^(1/2)*s^(-1)
	double c_light = 29979245800; //speed of light in vacuum, cm/s
	const double e_energy_rest = 510.99895; //Electron rest mass, keV
	const char *tmp_str; // temporary string
	/* TO ADD NEW BEAM
	 * Add here all intermediate variables, which are used only inside this function.
	 */

	// initialization of global option index for error messages
	opt=opt_beam;
	// beam initialization
	beam_asym=(beam_center_0[0]!=0 || beam_center_0[1]!=0 || beam_center_0[2]!=0);
	if (!beam_asym) vInit(beam_center); // not to calculate it for each orientation
	vorticity = 0;
	switch (beamtype) {
		case B_PLANE:
			if (IFROOT) beam_descr="plane wave";
			if (surface) {
				/* here we assume that prop_0 will not change further (e.g., by rotation of particle),
				 * i.e. prop=prop_0 in GenerateBeam() below
				 */
				if (prop_0[2]==0) PrintError("Ambiguous setting of beam propagating along the surface. Please specify "
					"the incident direction to have (arbitrary) small positive or negative z-component");
				if (msubInf && prop_0[2]>0) PrintError("Perfectly reflecting surface ('-surf ... inf') is incompatible "
					"with incident direction from below (including the default one)");
				// Here we set ki,kt,ktVec and propagation directions prIncRefl,prIncTran
				if (prop_0[2]>0) { // beam comes from the substrate (below)
					// here msub should always be defined
					inc_scale=1/creal(msub);
					ki=msub*prop_0[2];
					/* Special case for msub near 1 to remove discontinuities for near-grazing incidence. The details
					 * are discussed in CalcFieldSurf() in crosssec.c.
					 */
					if (cabs(msub-1)<ROUND_ERR && cabs(ki)<SQRT_RND_ERR) kt=ki;
					else kt=cSqrtCut(1 - msub*msub*(prop_0[0]*prop_0[0]+prop_0[1]*prop_0[1]));
					// determine propagation direction and full wavevector of wave transmitted into substrate
					ktVec[0]=msub*prop_0[0];
					ktVec[1]=msub*prop_0[1];
					ktVec[2]=kt;
				}
				else if (prop_0[2]<0) { // beam comes from above the substrate
					inc_scale=1;
					ki=-prop_0[2]; // always real
					if (!msubInf) {
						// same special case as above
						if (cabs(msub-1)<ROUND_ERR && cabs(ki)<SQRT_RND_ERR) kt=ki;
						else kt=cSqrtCut(msub*msub - (prop_0[0]*prop_0[0]+prop_0[1]*prop_0[1]));
						// determine propagation direction of wave transmitted into substrate
						ktVec[0]=prop_0[0];
						ktVec[1]=prop_0[1];
						ktVec[2]=-kt;
					}
				}
				else LogError(ONE_POS,"Ambiguous setting of beam propagating along the surface. Please specify the"
					"incident direction to have (arbitrary) small positive or negative z-component");
				vRefl(prop_0,prIncRefl);
				if (!msubInf) {
					vReal(ktVec,prIncTran);
					vNormalize(prIncTran);
				}
			}
			return;
		case B_DIPOLE:
			if (surface) {
				if (beam_center_0[2]<=-hsub)
					PrintErrorHelp("External dipole should be placed strictly above the surface");
				inc_scale=1; // but scaling of Mueller matrix is weird anyway
			}
			/* definition of p0 is important for scaling of many scattering quantities (that are normalized to incident
			 * irradiance). Alternative definition is p0=1, but then the results will scale with unit of length
			 * (breaking scale invariance)
			 */
			p0=1/(WaveNum*WaveNum*WaveNum); //Is it valid if WaveNum is complex?
			if (IFROOT) beam_descr="point dipole";
			return;
		case B_LMINUS:
		case B_DAVIS3:
		case B_BARTON5:
			if (surface) PrintError("Currently, Gaussian incident beam is not supported for '-surf'");
			// initialize parameters
			w0=beam_pars[0];
			TestPositive(w0,"beam width");
			s=1/(WaveNum*w0); //Is it valid if WaveNum is complex?
			s2=s*s;
			scale_x=1/w0;
			scale_z=s*scale_x; // 1/(k*w0^2)
			// beam info
			if (IFROOT) {
				switch (beamtype) {
					case B_LMINUS:
						tmp_str="L- approximation";
						break;
					case B_DAVIS3:
						tmp_str="3rd order approximation, by Davis";
						break;
					case B_BARTON5:
						tmp_str="5th order approximation, by Barton";
						break;
					default: LogError(ONE_POS,"Incompatibility error in GenerateB");
				}
				beam_descr=dyn_sprintf("Gaussian beam (%s)\n"
				                       "\tWidth="GFORMDEF" (confinement factor s="GFORMDEF")",tmp_str,w0,s);
			}
			return;
		case B_ELECTRON:
			//Electron field is in CGS. Electron field in SI would look the same, except multiplied by 1/(4*pi*eps0)
			//if (surface) PrintError("Currently, electron incident beam is not supported for '-surf'");
			// initialize parameters
			e_energy=beam_pars[0];
			TestPositive(e_energy,"kinetic energy of the electron");
			beam_asym=(beam_center_0[0]!=0 || beam_center_0[1]!=0 || beam_center_0[2]!=0);
			//symX=symY=symZ=symR=false;
			if (!beam_asym) vInit(beam_center);
			e_v = c_light*sqrt(1-pow((e_energy_rest/(e_energy+e_energy_rest)),2));
			beta1 = e_v*mhost/c_light;
			//gamma_eps_inv = cSqrtCut(1-beta_eps*beta_eps);
			gamma1_inv=1/cSqrtCut(1/(1-beta1*beta1));
			//printf("gamma_eps_inv\t=\t"CFORM"\n",REIM(gamma_eps_inv));
			//printf("w = "EFORM"\n",WaveNum0*c_light);
			//printf("v = "EFORM"\n",e_v);
			e_w_v = WaveNum0*c_light/e_v;
			//printf("e_w_v\t=\t"CFORM"\n",REIM(e_w_v));
			e_w_gv1 = e_w_v*gamma1_inv;
			//printf("e_w_gv\t=\t"CFORM"\n",REIM(e_w_gv));
			e_pref1 = q_electron*(WaveNum0*1e7)*gamma1_inv/(beta1*beta1*c_light); // (q*k0/c) must be in statV/cm => multiplying k0*1e7
			//printf("e_pref = "CFORM"\n",REIM(e_pref));
			if(surface){
				if(mhost != 1) PrintError("Currently, '-mhost' is not supported for '-surf'");
				inc_scale=1;
				vgw = e_v/(gamma1_inv*WaveNum0*c_light);
				vgw2 = vgw*vgw;
				//printf("vgw\t=\t"CFORM"\n",REIM(vgw));
				//printf("vgw2\t=\t"CFORM"\n",REIM(vgw2));
				eps1 = mhost*mhost;
				eps2 = msub*msub;
				//printf("msub\t=\t"CFORM"\n",REIM(msub));
			}
			if (IFROOT) beam_descr=dyn_sprintf("The electron with the %g keV energy moving through "GFORM3V,e_energy,COMP3V(beam_center_0));


			return;
#ifndef NO_FORTRAN
		case B_BES_CS:
		case B_BES_CSp:
		case B_BES_M:
		case B_BES_LE:
		case B_BES_LM:
		case B_BES_TEL:
		case B_BES_TML:
			if (surface) PrintError("Currently, Bessel incident beam is not supported for '-surf'");
			// initialize parameters
			ConvertToInteger(beam_pars[0],"beam order",&besN);
			TestRangeII(besN,"beam order (might cause the incorrect calculation of Bessel function)",-50,50);
			vorticity=besN;
			besAlpha=Deg2Rad(beam_pars[1]);
			besKt=WaveNum*sin(besAlpha);
			besKz=WaveNum*cos(besAlpha);
			switch (beamtype) { // definition of elements of matrix M ((Mex,Mey),(Mmx,Mmy))
				case B_BES_CS:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					besM[0]=0.5;  besM[1]=0;
					besM[2]=0;    besM[3]=0.5;
					if (IFROOT) tmp_str="circularly symmetric energy density";
					break;
				case B_BES_CSp:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					besM[0]=0.5;  besM[1]=0;
					besM[2]=0;    besM[3]=-0.5;
					if (IFROOT) tmp_str="circularly symmetric energy density, alternative type";
					break;
				case B_BES_M:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					if (beam_Npars==6) {
						besM[0]=beam_pars[2];
						besM[1]=beam_pars[3];
						besM[2]=beam_pars[4];
						besM[3]=beam_pars[5];
					}
					else {
						besM[0]=beam_pars[2]+I*beam_pars[6];
						besM[1]=beam_pars[3]+I*beam_pars[7];
						besM[2]=beam_pars[4]+I*beam_pars[8];
						besM[3]=beam_pars[5]+I*beam_pars[9];
					}
					if (IFROOT) tmp_str="generalized";
					break;
				case B_BES_LE:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					besM[0]=0;  besM[1]=0;
					besM[2]=0;  besM[3]=1;
					if (IFROOT) tmp_str="linear electric field";
					break;
				case B_BES_LM:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					besM[0]=0;  besM[1]=1;
					besM[2]=0;  besM[3]=0;
					if (IFROOT) tmp_str="linear magnetic field";
					break;
				/* TODO: for the following two types, both 0 and 90 degrees should be fine, but may require some
				 * rearrangement of formulae. Then the tests for range of alpha should be moved outside of the case
				 */
				case B_BES_TEL:
					TestRangeNN(beam_pars[1],"half-cone angle for TEL type",0,90);
					besM[0]=-WaveNum/besKt;  besM[1]=0;
					besM[2]=0;               besM[3]=besKz/besKt;
					if (IFROOT) tmp_str="linear component of the TE";
					break;
				case B_BES_TML:
					TestRangeNN(beam_pars[1],"half-cone angle for TML type",0,90);
					besM[0]=0;              besM[1]=besKz/besKt;
					besM[2]=WaveNum/besKt;  besM[3]=0;
					if (IFROOT) tmp_str="linear component of the TM";
					break;
				default: LogError(ONE_POS,"Incompatibility error in GenerateB");
			}
			// TODO: some symmetries can be retained in some special cases
			symR=symX=symY=symZ=false;
			// beam info
			if (IFROOT) beam_descr=dyn_sprintf("Bessel beam (%s)\n"
				                               "\tOrder: %d, half-cone angle: "GFORMDEF" deg",
				                               tmp_str,besN,beam_pars[1]);
			return;
#endif // !NO_FORTRAN
		case B_READ:
			// the safest is to assume cancellation of all symmetries
			symX=symY=symZ=symR=false;
			if (surface) inc_scale=1; // since we can't know it, we assume the default case
			if (IFROOT) {
				if (beam_Npars==1) beam_descr=dyn_sprintf("specified by file '%s'",beam_fnameY);
				else beam_descr=dyn_sprintf("specified by files '%s' and '%s'",beam_fnameY,beam_fnameX);
			}
			// this type is unaffected by beam_center
			return;
	}
	LogError(ONE_POS,"Unknown type of incident beam (%d)",(int)beamtype);
	/* TO ADD NEW BEAM
	 * add a case above. Identifier ('B_...') should be defined inside 'enum beam' in const.h. The case should
	 * 1) save all the input parameters from array 'beam_pars' to local variables (defined in the beginning of this
	 *    source files)
	 * 2) test all input parameters (for that you're encouraged to use functions from param.h since they would
	 *    automatically produce informative output in case of error).
	 * 3) the symmetry breaking due to prop or beam_center is taken care of in VariablesInterconnect() in param.c.
	 *    But if there are other reasons why beam would break any symmetry, corresponding variable should be set to
	 *    false here. Do not set any of them to true, as they can be set to false by other factors.
	 *    symX, symY, symZ - symmetries of reflection over planes YZ, XZ, XY respectively.
	 *    symR - symmetry of rotation for 90 degrees over the Z axis
	 * 4) initialize the following:
	 *    beam_descr - descriptive string, which will appear in log file  (should NOT end with \n).
	 *    vorticity - (only for vortex beam) integer value, how many turns the phase experience, when one makes a full
	 *                turn around the beam axis.
	 * 5) Consider the case of surface (substrate near the particle). If the new beam type is incompatible with it, add
	 *    an explicit exception, like "if (surface) PrintErrorHelp(...);". Otherwise, you also need to define inc_scale.
	 * All other auxiliary variables, which are used in beam generation (GenerateB(), see below), should be defined in
	 * the beginning of this file. If you need temporary local variables (which are used only in this part of the code),
	 * define them in the beginning of this function.
	 */
}

//======================================================================================================================

void GenerateB (const enum incpol which,   // x - or y polarized incident light
                doublecomplex *restrict b) // the b vector for the incident field
// generates incident beam at every dipole
{
	size_t i,j;
	doublecomplex psi0,Q,Q2;
	doublecomplex v1[3],v2[3],v3[3],gt[6],Einc_charge[3],Einc_image[3],Ea_charge[3],Ea_image[3];
	double ro,ro2,ro4;
	double x,y,z,z0,s0,R1,scale,x2_s,xy_s,temp;
	double tol; //incJ calculation tolerance
	doublecomplex t1,t2,t3,t4,t5,t6,t7,t8,ctemp,e_wb_gv;
	doublecomplex epp,epm,emp,emm,K0,K1,incJ0p,incJ0m,incJ1p,incJ1m,omp,omm,hzp,hzm,hzzp,hzzm,hzxp,hzxm;
	doublecomplex K0p,K0m,K1p,K1m,expikr,phase_med,phase_a,phase_z0,sub_xy,sub_z,med_xy,med_z,sub_pref;
	const double *ex; // coordinate axis of the beam reference frame
	double ey[3];
	double r1[3],r1par[3],r1per[3],beam_center_image[3],beam_center_proj[3],prop_image[3];
	/* complex wave amplitudes of transmitted wave (with phase relative to beam center);
	 * The transmitted wave can be inhomogeneous wave (when msub is complex), then eIncTran (e) is normalized
	 * counter-intuitively. Before multiplying by tc, it satisfies (e,e)=1!=||e||^2. This normalization is consistent
	 * with used formulae for transmission coefficients. So this transmission coefficient is not (generally) equal to
	 * the ratio of amplitudes of the electric fields. In particular, when E=E0*e, ||E||!=|E0|*||e||, where
	 * ||e||^2=(e,e*)=|e_x|^2+|e_y|^2+|e_z|^2=1
	 */
	doublecomplex eIncTran[3];
#ifndef NO_FORTRAN
	// for Bessel beams
	int n1,q;
	doublecomplex vort;  // vortex phase of Bessel beam rotated by 90 deg
	doublecomplex fn[5]; // for general functions f(n,ro,phi) of Bessel beams (fn-2, fn-1, fn, fn+1, fn+2, respectively)
	double phi,arg,td1[abs(besN)+3],td2[abs(besN)+3],jn1[abs(besN)+3]; // for Bessel beams
#endif
	const char *fname;

	/* TO ADD NEW BEAM
	 * Add here all intermediate variables, which are used only inside this function. You may as well use 't1'-'t8'
	 * variables defined above.
	 */

	// set reference frame of the beam; ez=prop, ex - incident polarization
	if (which==INCPOL_Y) {
		ex=incPolY;
		vMultScal(-1,incPolX,ey);
	}
	else { // which==INCPOL_X
		ex=incPolX;
		vCopy(incPolY,ey);
	}

	//printf("beam_center_0 = "EFORM3V"\n",COMP3V(beam_center_0));
	//printf("beam_center = "EFORM3V"\n",COMP3V(beam_center));

	switch (beamtype) {
		case B_PLANE: // plane is separate to be fast (for non-surface)
			if (surface) {
				/* With respect to normalization we use here the same assumption as in the free space - the origin is in
				 * the particle center, and amplitude of incoming plane wave is equal to 1. Then irradiance of the beam
				 * coming from below is c*Re(msub)/(8pi), different from that coming from above.
				 * Original incident (incoming) beam propagating from the vacuum (above) is Exp(i*k*r.a), while - from
				 * the substrate (below) is Exp(i*k*msub*r.a). We assume that the incoming beam is homogeneous in its
				 * original medium.
				 */
				double hbeam=hsub+beam_center[2]; // height of beam center above the surface
				if (prop[2]>0) { // beam comes from the substrate (below)
					doublecomplex tc; // transmission coefficients
					//  determine amplitude of the transmitted wave; here msub is always defined
					if (which==INCPOL_Y) { // s-polarized
						cvBuildRe(ex,eIncTran);
						tc=FresnelTS(ki,kt);
					}
					else { // p-polarized
						crCrossProd(ey,ktVec,eIncTran);
						tc=FresnelTP(ki,kt,1/msub);
					}
					// phase shift due to the beam center relative to the origin and surface
					cvMultScal_cmplx(tc*cexp(I*WaveNum*((kt-ki)*hbeam-crDotProd(ktVec,beam_center))),eIncTran,eIncTran);
					// main part
					for (i=0;i<local_nvoid_Ndip;i++) {
						j=3*i;
						// b[i] = eIncTran*exp(ik*kt.r); subtraction of beam_center is avoided by the factor above
						cvMultScal_cmplx(cexp(I*WaveNum*crDotProd(ktVec,DipoleCoord+j)),eIncTran,b+j);
					}
				}
				else if (prop[2]<0) { // beam comes from above the substrate
					/* The following code takes extra care to be stable (not to lose precision) for grazing incidence.
					 * While it seems more complicated for general incidence, it requires only a few more complex
					 * arithmetic operations that are negligible compared to two complex exponents. For grazing
					 * incidence, it is not only stable but may also be faster than the straightforward computation,
					 * since one of the complex exponents is computed from a very small argument
					 */
					// determine reflection coefficient + 1
					doublecomplex rcp1;
					if (which==INCPOL_Y) { // s-polarized
						if (msubInf) rcp1=0;
						else rcp1=FresnelTS(ki,kt);
					}
					else { // p-polarized
						if (msubInf) rcp1=2;
						else rcp1=msub*FresnelTP(ki,kt,msub);
					}
					// main part
					for (i=0;i<local_nvoid_Ndip;i++) {
						j=3*i;
						vSubtr(DipoleCoord+j,beam_center,r1); // beam_center affects only the common phase factor
						/* b[i] = ex*exp(ik*r.a) + eIncRefl*exp[ik(prIncRefl.r-2az*hbeam)],
						 * where prIncRefl is prop with reflected z-component, while eIncRefl=rc*ex for s-polarization
						 * and (additionally) with reflected transverse (x,y) components for p-polarization
						 */
						ctemp=imExp(WaveNum*DotProd(r1,prop)); // exp(ik*r.a)
						t1=imExpM1(-2*WaveNum*prop[2]*(r1[2]+hbeam));
						/* t2=exp(ik*r.a){1+rc*exp[-2i*kz(z+hbeam)]}, but the terms in parentheses (rc and exp(...)) are
						 * replaced by their differences with -1 and 1, respectively (eliminating 1 inside)
						 */
						t2=ctemp*(rcp1*(t1+1)-t1);
						cvMultScal_RVec(t2,ex,b+j);
						if (which==INCPOL_X) {
							/* here we add (eIncRefl-ex*rc)*exp[ik(prIncRefl.r-2az*hbeam)] (non-zero only for
							 * p-polarization), expressing rc*exp(...) as t2-ctemp
							 */
							t3=2*(t2-ctemp);
							b[j]-=t3*ex[0];
							b[j+1]-=t3*ex[1];
						}
					}
				}
			}
			else for (i=0;i<local_nvoid_Ndip;i++) { // standard (non-surface) plane wave
				j=3*i;
				// can be replaced by complex multiplication (by precomputed phase factor), does not seem beneficial
				vSubtr(DipoleCoord+j,beam_center,r1);
				ctemp=imExp(WaveNum*DotProd(r1,prop)); // ctemp=exp(ik*r.a)
				cvMultScal_RVec(ctemp,ex,b+j); // b[i]=ctemp*ex
			}
			return;
		case B_DIPOLE: {
			double dip_p[3]; // dipole moment, = p0*prop
			vMultScal(p0,prop,dip_p);
			for (i=0;i<local_nvoid_Ndip;i++) { // here we explicitly use that dip_p is real
				j=3*i;
				vSubtr(DipoleCoord+j,beam_center,r1);
				if(r1[0] < dsX) r1[0] = dsX;
				if(r1[1] < dsX) r1[1] = dsX;
				if(r1[2] < dsX) r1[2] = dsX;
				(*InterTerm_real)(r1,gt);
				cSymMatrVecReal(gt,dip_p,b+j);
				if (surface) { // add reflected field
					r1[2]=DipoleCoord[j+2]+beam_center[2]+2*hsub;
					(*ReflTerm_real)(r1,gt);
					cReflMatrVecReal(gt,dip_p,v1);
					cvAdd(v1,b+j,b+j);
				}
			}
			/* calculate dipole inherent cross sections (for decay rate enhancements)
			 * General formula is C0=4pi*k*Im(p0*.G(r0,r0).p0) and it is used for reflection; but for direct
			 * interaction it is expected to be singular, so we use an explicit formula for point dipole:
			 * C0=(8pi/3)*k^4*|p0|^2. Thus we discard the choice of "-int ...", but still the used value should be
			 * correct for most formulations, e.g. poi, fcd, fcd_st, igt_so. Moreover, it is also logical, since the
			 * exciting dipole is really point one, in contrast to the dipoles composing the particle.
			 */
			double temp=p0*WaveNum*WaveNum; // in principle, t1=1/k, but we keep a general formula
			C0dipole=2*FOUR_PI_OVER_THREE*temp*temp;
			if (surface) {
				r1[0]=r1[1]=0;
				r1[2]=2*(beam_center[2]+hsub);
				(*ReflTerm_real)(r1,gt);
				double tmp;
				/* the following expression uses that dip_p is real and a specific (anti-)symmetry of the gt
				 * a general expression is commented out below
				 */
				tmp=dip_p[0]*dip_p[0]*cimag(gt[0])+2*dip_p[0]*dip_p[1]*cimag(gt[1])+dip_p[1]*dip_p[1]*cimag(gt[3])
				   +dip_p[2]*dip_p[2]*cimag(gt[5]);
//				v1[0]=dip_p[0]; v1[1]=dip_p[1]; v1[2]=dip_p[2];
//				cReflMatrVec(gt,v1,v2);
//				tmp=cDotProd_Im(v2,v1);
				C0dipole_refl=FOUR_PI*WaveNum*tmp;
			}
			return;
		}
		case B_LMINUS:
		case B_DAVIS3:
		case B_BARTON5:
			for (i=0;i<local_nvoid_Ndip;i++) {
				j=3*i;
				// set relative coordinates (in beam's coordinate system)
				vSubtr(DipoleCoord+j,beam_center,r1);
				x=DotProd(r1,ex)*scale_x;
				y=DotProd(r1,ey)*scale_x;
				z=DotProd(r1,prop)*scale_z;
				ro2=x*x+y*y;
				Q=1/(2*z-I);
				psi0=-I*Q*imExp(Q*ro2);
				// ctemp=exp(ik*z0)*psi0, z0 - non-scaled coordinate (z/scale_z)
				ctemp=imExpReal(WaveNum*z/scale_z)*psi0; //or should be imExp?
				// the following logic (if-else-if...) is hard to replace by a simple switch
				if (beamtype==B_LMINUS) cvMultScal_RVec(ctemp,ex,b+j); // b[i]=ctemp*ex
				else {
					/* It is possible to rewrite the formulae below to avoid division by ro2, but we prefer
					 * dimensionless variables. The value for ro2=0 doesn't really matter (cancels afterwards).
					 * The current code should work OK even for very small ro2
					 */
					if (ro2==0) x2_s=0;
					else x2_s=x*x/ro2;
					Q2=Q*Q;
					ro4=ro2*ro2;
					// some combinations that are used more than once
					t4=s2*ro2*Q2; // t4=(s*ro*Q)^2
					t5=I*ro2*Q;   // t5=i*Q*ro^2
					t6=ro4*Q2;    // t6=ro^4*Q^2
					t7=x*s*Q;     // t7=x*s*Q
					if (beamtype==B_DAVIS3) {
						// t1=1+s^2(-4Q^2*x^2-iQ^3*ro^4)=1-t4(4x2_s+t5)
						t1 = 1 - t4*(4*x2_s+t5);
						// t2=0
						t2=0;
						// t3=-s(2Qx)+s^3(8Q^3*ro^2*x+2iQ^4*ro^4*x-4iQ^2x)=2t7[-1+iQ*s2*(-4t5+t6-2)]
						t3 = 2*t7*(-1 + I*Q*s2*(-4*t5+t6-2));
					}
					else if (beamtype==B_BARTON5) {
						if (ro2==0) xy_s=0; // see comment for x2_s above
						else xy_s=x*y/ro2;
						t8=8+2*t5; // t8=8+2i*Q*ro^2
						/* t1 = 1 + s^2(-ro^2*Q^2-i*ro^4*Q^3-2Q^2*x^2)
						 *    + s^4[2ro^4*Q^4+3iro^6*Q^5-0.5ro^8*Q^6+x^2(8ro^2*Q^4+2iro^4*Q^5)]
						 *    = 1 + t4*{-1-2x2_s-t5+t4*[2+3t5-0.5t6+x2_s*t8]}
						 */
						t1 = 1 + t4*(-1 - 2*x2_s - t5 + t4*(2+3*t5-0.5*t6+x2_s*t8));
						// t2=s^2(-2Q^2*xy)+s^4[xy(8ro^2*Q^4+2iro^4*Q^5)]=xy_s*t4(-2+t4*t8)
						t2=xy_s*t4*(-2+t4*t8);
						/* t3 = s(-2Qx) + s^3[(6ro^2*Q^3+2iro^4*Q^4)x] + s^5[(-20ro^4*Q^5-10iro^6*Q^6+ro^8*Q^7)x]
						 *    = t7{-2+t4[6+2t5+t4(-20-10t5+t6)]}
						 */
						t3 = t7*(-2 + t4*(6 + 2*t5 + t4*(-20-10*t5+t6)));
					}
					else LogError(ONE_POS,"Inconsistency in beam definition"); // to remove warnings
					// b[i]=ctemp(ex*t1+ey*t2+ez*t3)
					cvMultScal_RVec(t1,ex,v1);
					cvMultScal_RVec(t2,ey,v2);
					cvMultScal_RVec(t3,prop,v3);
					cvAdd2Self(v1,v2,v3);
					cvMultScal_cmplx(ctemp,v1,b+j);
				}
			}
			return;

		//new code
		case B_ELECTRON:
			if(!surface) hsub = 0;
			cvMultScal_cmplx(1,beam_center,beam_center_proj);
			if(fabs(prop[2]) > DBL_EPSILON){ //For tilted beam: projection of beam center on z=-hsub
				if(fabs(prop[0]) > DBL_EPSILON) beam_center_proj[0] = -(prop[0]/prop[2])*(beam_center[2]+hsub) + beam_center[0];
				if(fabs(prop[1]) > DBL_EPSILON) beam_center_proj[1] = -(prop[1]/prop[2])*(beam_center[2]+hsub) + beam_center[1];
				beam_center_proj[2] = -hsub;
				printf("Beam center = (%g %g %g)\n",beam_center[0],beam_center[1],beam_center[2]);
				printf("Beam center proj = (%g %g %g)\n",beam_center_proj[0],beam_center_proj[1],beam_center_proj[2]);
				//LogError(ONE_POS,"STOP");
			}

			for (i=0;i<local_nvoid_Ndip;i++) {
				j=3*i;

				if(!surface) {
					LinComb(DipoleCoord+j,beam_center,1,-1,r1); // set relative coordinates (in beam's coordinate system)
					temp = DotProd(r1,prop);
					vMultScal(temp,prop,r1par);
					vSubtr(r1,r1par,r1per);
					ro = vNorm(r1per);
					z = temp;
					if(ro != 0) vNormalize(r1per);
					if(ro < dsX) ro = dsX;
					e_wb_gv = e_w_gv1*ro;
					cik01_(&e_wb_gv, &t1, &t2, &t3, &t4, &t7, &t5, &t8, &t6);

					t4 = imExp(e_w_v*z);
					cvMultScal_RVec(t4*t8,r1per,v1);
					cvMultScal_RVec((-I)*gamma1_inv*t4*t7,prop,v2);
					cvAdd(v1,v2,v3);
					cvMultScal_cmplx(2*e_pref1,v3,Einc_charge); //Einc from the charge

					t4 = conj(t4);
					cvMultScal_RVec(-t4*t8,r1per,v1);
					cvMultScal_RVec((-I)*gamma1_inv*t4*t7,prop,v2);
					cvAdd(v1,v2,v3);
					cvMultScal_cmplx(2*e_pref1,v3,Ea_charge); //Ea from the charge

					cvMultScal_cmplx(1,Einc_charge,b+j); //E_inc
					cvMultScal_cmplx(1,Ea_charge,E1+j); //E_a
				}
				else { //substrate
					if(fabs(prop[2]) < 2e-3) { //beam parallel or numerically near-parallel to the surface
						if(beam_center[2] >= -hsub) { //beam above the surface
							//charge
							LinComb(DipoleCoord+j,beam_center,1,-1,r1); // set relative coordinates (in beam's coordinate system)
							temp = DotProd(r1,prop);
							vMultScal(temp,prop,r1par);
							vSubtr(r1,r1par,r1per);
							ro = vNorm(r1per);
							z = temp;
							if(ro != 0) vNormalize(r1per);
							if(ro < dsX) ro = dsX;
							e_wb_gv = e_w_gv1*ro;
							cik01_(&e_wb_gv, &t1, &t2, &t3, &t4, &t7, &t5, &t8, &t6);

							t4 = imExp(e_w_v*z);
							cvMultScal_RVec(t4*t8,r1per,v1);
							cvMultScal_RVec((-I)*gamma1_inv*t4*t7,prop,v2);
							cvAdd(v1,v2,v3);
							cvMultScal_cmplx(2*e_pref1,v3,Einc_charge); //Einc from the charge

							t4 = conj(t4);
							cvMultScal_RVec(-t4*t8,r1per,v1);
							cvMultScal_RVec((-I)*gamma1_inv*t4*t7,prop,v2);
							cvAdd(v1,v2,v3);
							cvMultScal_cmplx(2*e_pref1,v3,Ea_charge); //Ea from the charge


							//image charge
							cvMultScal_cmplx(1,beam_center,beam_center_image);
							beam_center_image[2] = -beam_center[2]-2*hsub;

							LinComb(DipoleCoord+j,beam_center_image,1,-1,r1); // set relative coordinates (in beam's coordinate system)
							temp = DotProd(r1,prop);
							vMultScal(temp,prop,r1par);
							vSubtr(r1,r1par,r1per);
							ro = vNorm(r1per);
							z = temp;
							if(ro != 0) vNormalize(r1per);
							if(ro < dsX) ro = dsX;
							e_wb_gv = e_w_gv1*ro;
							cik01_(&e_wb_gv, &t1, &t2, &t3, &t4, &t7, &t5, &t8, &t6);

							t4 = imExp(e_w_v*z);
							cvMultScal_RVec(t4*t8,r1per,v1);
							cvMultScal_RVec((-I)*gamma1_inv*t4*t7,prop,v2);
							cvAdd(v1,v2,v3);
							cvMultScal_cmplx(((eps1-eps2)/(eps1+eps2))*2*e_pref1,v3,Einc_image); //Einc from the image

							t4 = conj(t4);
							cvMultScal_RVec(-t4*t8,r1per,v1);
							cvMultScal_RVec((-I)*gamma1_inv*t4*t7,prop,v2);
							cvAdd(v1,v2,v3);
							cvMultScal_cmplx(((eps1-eps2)/(eps1+eps2))*2*e_pref1,v3,Ea_image); //Ea from the image

							cvAdd(Einc_charge,Einc_image,b+j);
							cvAdd(Ea_charge,Ea_image,E1+j);
						}
						else { //beam below the surface
							LinComb(DipoleCoord+j,beam_center,1,-1,r1); // set relative coordinates (in beam's coordinate system)
							temp = DotProd(r1,prop);
							vMultScal(temp,prop,r1par);
							vSubtr(r1,r1par,r1per);
							ro = vNorm(r1per);
							z = temp;
							if(ro != 0) vNormalize(r1per);
							if(ro < dsX) ro = dsX;
							e_wb_gv = e_w_gv1*ro;
							cik01_(&e_wb_gv, &t1, &t2, &t3, &t4, &t7, &t5, &t8, &t6);

							t4 = imExp(e_w_v*z);
							cvMultScal_RVec(t4*t8,r1per,v1);
							cvMultScal_RVec((-I)*gamma1_inv*t4*t7,prop,v2);
							cvAdd(v1,v2,v3);
							cvMultScal_cmplx(2*e_pref1,v3,Einc_charge); //Einc from the charge

							t4 = conj(t4);
							cvMultScal_RVec(-t4*t8,r1per,v1);
							cvMultScal_RVec((-I)*gamma1_inv*t4*t7,prop,v2);
							cvAdd(v1,v2,v3);
							cvMultScal_cmplx(2*e_pref1,v3,Ea_charge); //Ea from the charge

							cvMultScal_cmplx(2*eps1/(eps1+eps2),Einc_charge,b+j); //E_inc
							cvMultScal_cmplx(2*eps1/(eps1+eps2),Ea_charge,E1+j); //E_a
						}
					}


					//else if(prop[0] < DBL_EPSILON && prop[1] < DBL_EPSILON) { //perpendicular to the substrate
					else if(fabs(prop[0]) < DBL_EPSILON && fabs(prop[1]) < DBL_EPSILON && prop[2] == -1) { //perpendicular to the substrate
						//LogError(ONE_POS,"old");
						//printf("--- j = %i\n",j);
						//printf("prop\t=\t%g %g %g\n",prop[0],prop[1],prop[2]);

						LinComb(DipoleCoord+j,beam_center,1,-1,r1);

						ro = sqrt(r1[0]*r1[0]+r1[1]*r1[1]);
						z = (DipoleCoord+j)[2] + hsub;
						if(ro < dsX/2) ro = dsX/2;
						//printf("b = %lf\n",ro);
						//printf("z = %lf\n",z);
						if(ro == 0) LogError(ONE_POS,"electron hit a dipole, this is currently not supported, ro = "EFORM, ro);
						z0 = beam_center[2] + hsub;
						e_wb_gv = e_w_gv1*ro;
						cik01_(&e_wb_gv, &t1, &t2, &t3, &t4, &K0, &t5, &K1, &t6);
						//printf("e_wb_gv\t=\t"CFORM"\n",REIM(e_wb_gv));
						//printf("K0\t=\t"CFORM"\n",REIM(K0));
						//printf("K1\t=\t"CFORM"\n",REIM(K1));

						/* Simplified expressions from Eqs. (53), (54), and (56) of the thesis.
						 * They are algebraically equivalent to the older implementation above, but avoid
						 * explicit h_z/h_zz/h_zx combinations. */
						t4 = imExp(e_w_v*z);    // exp(i*omega*z/v)
						phase_z0 = imExp(e_w_v*z0); // exp(i*omega*z0/v)
						phase_med = conj(t4)*phase_z0; // exp(-i*omega*(z-z0)/v)
						omp = asinh(z/ro)+ArcCoth(beta1);
						omm = asinh(z/ro)-ArcCoth(beta1);

						tol = 1e-4; // precision for incomplete Bessel functions
						incJ0p = incJ(0,e_wb_gv,omp,tol);
						incJ0m = incJ(0,e_wb_gv,omm,tol);
						incJ1p = incJ(1,e_wb_gv,omp,tol);
						incJ1m = incJ(1,e_wb_gv,omm,tol);

						R1 = 1/sqrt(ro*ro+z*z); // 1/r_tilde
						K1p = conj(t4)*(K1-incJ1p); // K_1^+
						K1m = t4*(K1-incJ1m);       // K_1^-
						K0p = conj(t4)*(K0-incJ0p); // K_0^+
						K0m = t4*(K0-incJ0m);       // K_0^-
						expikr = imExp(WaveNum0*sqrt(ro*ro+z*z)); // exp(i*k*r_tilde); direct form avoids ArcCoth/ccosh branch issues

						med_xy = -2*e_pref1*phase_med*K1;
						med_z = (-I)*gamma1_inv*2*e_pref1*phase_med*K0;
						sub_pref = -e_pref1*((eps1-eps2)/(eps1+eps2))*phase_z0;
						sub_xy = sub_pref*(K1p+K1m-2*z*R1/e_wb_gv*expikr);
						sub_z = sub_pref*(I*gamma1_inv*(K0p-K0m)+2*R1/e_w_gv1*expikr);

						Einc_charge[0] = (r1[0]/ro)*med_xy;
						Einc_charge[1] = (r1[1]/ro)*med_xy;
						Einc_charge[2] = med_z;
						Einc_image[0] = (r1[0]/ro)*sub_xy;
						Einc_image[1] = (r1[1]/ro)*sub_xy;
						Einc_image[2] = sub_z;
						cvAdd(Einc_charge,Einc_image,b+j);

						phase_a = imExp(2*e_w_v*(z-z0));
						Ea_charge[0] = -phase_a*Einc_charge[0];
						Ea_charge[1] = -phase_a*Einc_charge[1];
						Ea_charge[2] = phase_a*Einc_charge[2];
						cvMultScal_cmplx(-imExp(-2*e_w_v*z0),Einc_image,Ea_image);
						cvAdd(Ea_charge,Ea_image,E1+j);
					}

					else { //tilted
						doublecomplex refl;
						double rho_real,temp_real,scale_core;
						double core_obs[3];
						if(prop[2] > 0) LogError(ONE_POS,"Beam propagating from below the substrate is currently not supported");

						// reflected trajectory used for R_img term
						cvMultScal_cmplx(1,beam_center,beam_center_image);
						beam_center_image[2] = -beam_center[2]-2*hsub;
						vCopy(prop,prop_image);
						prop_image[2] = -prop[2];
						refl = (eps1-eps2)/(eps1+eps2);

						// distance from beam center to the interface crossing along uhat=-prop
						LinComb(beam_center,beam_center_proj,1,-1,r1);
						s0 = -DotProd(r1,prop);
						LinComb(beam_center_image,beam_center_proj,1,-1,r1);
						z0 = -DotProd(r1,prop_image); // reuse z0 as s0 for image frame

						// real-beam distance and line/cube intersection decide whether we regularize the singular real core
						LinComb(DipoleCoord+j,beam_center_proj,1,-1,r1);
						temp_real = DotProd(r1,prop);
						vMultScal(temp_real,prop,r1par);
						vSubtr(r1,r1par,r1per);
						rho_real = vNorm(r1per);

						if (tilted_regularization_enabled && BuildTiltedGridRegularizedPoint(DipoleCoord+j,beam_center_proj,prop,core_obs,&rho_real)) {
							// projected-grid regularization of the whole exact tilted field inside the beam-hit cylinder
							ElectronTiltedFieldPointNoFloor(core_obs,beam_center,beam_center_proj,prop,prop_image,s0,z0,refl,false,b+j);
							ElectronTiltedFieldPointNoFloor(core_obs,beam_center,beam_center_proj,prop,prop_image,s0,z0,refl,true,E1+j);
						}
						else {
							// exact tilted field outside the regularized beam-hit cylinder (or for unsupported out-of-plane tilt)
							ElectronTiltedFieldPointNoFloor(DipoleCoord+j,beam_center,beam_center_proj,prop,prop_image,s0,z0,refl,false,b+j);
							ElectronTiltedFieldPointNoFloor(DipoleCoord+j,beam_center,beam_center_proj,prop,prop_image,s0,z0,refl,true,E1+j);
						}
					}
				}
			}
			return;

		// old code
		/*case B_ELECTRON:
					for (i=0;i<local_nvoid_Ndip;i++) {
						j=3*i;
						// set relative coordinates (in beam's coordinate system)
						LinComb(DipoleCoord+j,beam_center,1,-1,r1);


						if(surface){
							//printf("--- j = %i\n",j);
							//printf("prop\t=\t%g %g %g\n",prop[0],prop[1],prop[2]);
							if(prop[0] != 0 || prop[1] != 0 || prop[2] != -1) PrintError("Currently, only '-prop 0 0 -1' is supported for '-surf'");
							ro = sqrt(r1[0]*r1[0]+r1[1]*r1[1]);
							z = (DipoleCoord+j)[2] + hsub;
							//printf("b = %lf\n",ro);
							//printf("z = %lf\n",z);
							if(ro == 0) LogError(ONE_POS,"electron hit a dipole, this is currently not supported, ro = "EFORM, ro);
							z0 = beam_center[2] + hsub;
							e_wb_gv = e_w_gv1*ro;
							cik01_(&e_wb_gv, &t1, &t2, &t3, &t4, &K0, &t5, &K1, &t6);
							//printf("e_wb_gv\t=\t"CFORM"\n",REIM(e_wb_gv));
							//printf("K0\t=\t"CFORM"\n",REIM(K0));
							//printf("K1\t=\t"CFORM"\n",REIM(K1));

							epp = imExp(e_w_v*(z+z0));
							epm = imExp(e_w_v*(z-z0));
							emp = imExp(-e_w_v*(z+z0));
							emm = imExp(-e_w_v*(z-z0));
							//printf("epp\t=\t"CFORM"\n",REIM(epp));
							//printf("epm\t=\t"CFORM"\n",REIM(epm));
							//printf("emp\t=\t"CFORM"\n",REIM(emp));
							//printf("emm\t=\t"CFORM"\n",REIM(emm));
							omp = asinh(z/ro)+ArcCoth(beta1);
							omm = asinh(z/ro)-ArcCoth(beta1);
							//printf("Omega+\t=\t"CFORM"\n",REIM(omp));
							//printf("Omega-\t=\t"CFORM"\n",REIM(omm));

							tol = 1e-8; //how many digits of precision we will have in incJ calculation
							incJ0p = incJ(0,e_wb_gv,omp,tol);
							//printf("incJ0+\t=\t"CFORM"\n",REIM(incJ0p));
							incJ0m = incJ(0,e_wb_gv,omm,tol);
							//printf("incJ0-\t=\t"CFORM"\n",REIM(incJ0m));
							incJ1p = incJ(1,e_wb_gv,omp,tol);
							//printf("incJ1+\t=\t"CFORM"\n",REIM(incJ1p));
							incJ1m = incJ(1,e_wb_gv,omm,tol);
							//printf("incJ1-\t=\t"CFORM"\n",REIM(incJ1m));

							R1 = (1/sqrt(ro*ro+z*z));
							hzp = R1*cexp(-e_wb_gv*ccosh(omp));
							hzm = R1*cexp(-e_wb_gv*ccosh(omm));
							//printf("hz+\t=\t"CFORM"\n",REIM(hzp));
							//printf("hz-\t=\t"CFORM"\n",REIM(hzm));
							hzzp = -R1*hzp*(z*R1+e_wb_gv*csinh(omp));
							hzzm = -R1*hzm*(z*R1+e_wb_gv*csinh(omm));
							//printf("hzz+\t=\t"CFORM"\n",REIM(hzzp));
							//printf("hzz-\t=\t"CFORM"\n",REIM(hzzm));
							hzxp = hzp*((1./vgw)*(z*R1*csinh(omp)-ccosh(omp)) - ro*R1*R1);
							hzxm = hzm*((1./vgw)*(z*R1*csinh(omm)-ccosh(omm)) - ro*R1*R1);
							//printf("hzx+\t=\t"CFORM"\n",REIM(hzxp));
							//printf("hzx-\t=\t"CFORM"\n",REIM(hzxm));



							//Einc: x,y,z
							(b+j)[0] = -e_pref1*(r1[0]/ro)*(1/(eps1+eps2))*(
									K1*(emm*(3*eps1+eps2)+epp*(eps1-eps2))
									- (eps1-eps2)*(
											emm*(incJ1p + vgw*(z/ro)*hzp - I*(vgw/e_w_v)*hzxp) +
											epp*(incJ1m + vgw*(z/ro)*hzm + I*(vgw/e_w_v)*hzxm)
											));
							//printf("EincX\t=\t"CFORM"\n",REIM((b+j)[0]));
							(b+j)[1] = -e_pref1*(r1[1]/ro)*(1/(eps1+eps2))*(
									K1*(emm*(3*eps1+eps2)+epp*(eps1-eps2))
									- (eps1-eps2)*(
											emm*(incJ1p + vgw*(z/ro)*hzp - I*(vgw/e_w_v)*hzxp) +
											epp*(incJ1m + vgw*(z/ro)*hzm + I*(vgw/e_w_v)*hzxm)
											));
							//printf("EincY\t=\t"CFORM"\n",REIM((b+j)[1]));
							(b+j)[2] = (-I*gamma1_inv)*e_pref1*(1/(eps1+eps2))*(
									K0*(emm*(3*eps1+eps2)-epp*(eps1-eps2))
									- (eps1-eps2)*(
											emm*(incJ0p + vgw2*(2*I*e_w_v*hzp - hzzp)) -
											epp*(incJ0m - vgw2*(2*I*e_w_v*hzm + hzzm))
											));
							//printf("EincZ\t=\t"CFORM"\n",REIM((b+j)[2]));

							//E1: x,y,z
							(E1+j)[0] = e_pref1*(r1[0]/ro)*(1/(eps1+eps2))*(
									K1*(epm*(3*eps1+eps2) + emp*(eps1-eps2))
									- (eps1-eps2)*(
											epm*(incJ1m + vgw*(z/ro)*hzm + I*(vgw/e_w_v)*hzxm) +
											emp*(incJ1p + vgw*(z/ro)*hzp - I*(vgw/e_w_v)*hzxp)
											));
							//printf("E1X\t=\t"CFORM"\n",REIM((E1+j)[0]));
							(E1+j)[1] = e_pref1*(r1[1]/ro)*(1/(eps1+eps2))*(
									K1*(epm*(3*eps1+eps2) + emp*(eps1-eps2))
									- (eps1-eps2)*(
											epm*(incJ1m + vgw*(z/ro)*hzm + I*(vgw/e_w_v)*hzxm) +
											emp*(incJ1p + vgw*(z/ro)*hzp - I*(vgw/e_w_v)*hzxp)
											));
							//printf("E1Y\t=\t"CFORM"\n",REIM((E1+j)[1]));
							(E1+j)[2] = (-I*gamma1_inv)*e_pref1*(1/(eps1+eps2))*(
									K0*(epm*(3*eps1+eps2) - emp*(eps1-eps2))
									- (eps1-eps2)*(
											epm*(incJ0m - vgw2*(2*I*e_w_v*hzm + hzzm)) -
											emp*(incJ0p + vgw2*(2*I*e_w_v*hzp - hzzp))
											));
							//printf("E1Z\t=\t"CFORM"\n",REIM((E1+j)[2]));

						}
						else{
							temp = DotProd(r1,prop);
							vMultScal(temp,prop,r1par);
							vSubtr(r1,r1par,r1per);
							ro = vNorm(r1per);
							z = temp;
							if(ro != 0) vNormalize(r1per);
							//else LogError(ONE_POS,"electron hit a dipole, this is currently not supported, ro = "EFORM, ro);
							if(ro < dsX/2) ro = dsX/2;

							//if(prop[2] < DBL_EPSILON && surface) LogError(ONE_POS,"ZAEBIS");

							//printf("z = "EFORM"\n",z);
							//printf("b = "EFORM"\n",ro);
							//printf("gamma = "EFORM"\n",creal(1./gamma_eps_inv));

							e_wb_gv = e_w_gv1*ro;
							cik01_(&e_wb_gv, &t1, &t2, &t3, &t4, &t7, &t5, &t8, &t6);
							//printf("e_wb_gv\t=\t"CFORM"\n",REIM(e_wb_gv));
							//printf("besselK0re\t=\t"EFORM"\n",creal(t7));
							//printf("besselK0im\t=\t"EFORM"\n",cimag(t7));
							//printf("besselK1re\t=\t"EFORM"\n",creal(t8));
							//printf("besselK1im\t=\t"EFORM"\n",cimag(t8));

							t4 = imExp(e_w_v*z);
							//printf("imExp = "CFORM"\n",t4);
							cvMultScal_RVec(t4*t8,r1per,v1);
							cvMultScal_RVec((-I)*gamma1_inv*t4*t7,prop,v2);
							cvAdd(v1,v2,v3);
							cvMultScal_cmplx(2*e_pref1,v3,b+j); //E_inc
							//printf("Einc\t=\t"CFORM3V"\n",REIM3V(b+j));

							t4 = conj(t4);
							cvMultScal_RVec(-t4*t8,r1per,v1);
							cvMultScal_RVec((-I)*gamma1_inv*t4*t7,prop,v2);
							cvAdd(v1,v2,v3);
							cvMultScal_cmplx(2*e_pref1,v3,E1+j); //E_1
							//printf("E1\t=\t"CFORM3V"\n",REIM3V(E1+j));
						}

					}
					return;*/

#ifndef NO_FORTRAN
		case B_BES_CS:
		case B_BES_CSp:
		case B_BES_M:
		case B_BES_LE:
		case B_BES_LM:
		case B_BES_TEL:
		case B_BES_TML:
			/* we assume that matrix M determines x-polarization (in the beam reference frame), while y-one is obtained
			 * by rotation with additional phase shift
			 */
			vort=(which==INCPOL_Y) ? cpow(I,besN) : 1;
			for (i=0;i<local_nvoid_Ndip;i++) {
				j=3*i;
				vSubtr(DipoleCoord+j,beam_center,r1);
				x=DotProd(r1,ex);
				y=DotProd(r1,ey);
				z=DotProd(r1,prop);
				phi=atan2(y,x); // angular coordinate in a cylindrical coordinate system
				ctemp=imExp(besN*phi)*cexp(I*besKz*z)*vort/(WaveNum*WaveNum); // common factor
				arg=besKt*sqrt(x*x+y*y); // argument of Bessel functions
				if (arg<ROUND_ERR) {
					// TODO: the following seems incorrect for n=+-1 or +-2 (other fn will be non-zero)
					if (besN==0) fn[2]=1.;
					else fn[2]=0;
					fn[0]=fn[1]=fn[3]=fn[4]=0;
				}
				else {
					// TODO: the following looks very complicated, try to simplify
					n1=abs(besN)+2;
					if (besN<=-3) {
						bjndd_(&n1,&arg,jn1,td1,td2);
						if (n1%2==0) q=1;
						else q=-1;
						fn[0] =  q*jn1[-besN+2]*imExp(-2*phi);
						fn[1] = -q*jn1[-besN+1]*imExp(-phi);
						fn[2] =  q*jn1[-besN];
						fn[3] = -q*jn1[-besN-1]*imExp(phi);
						fn[4] =  q*jn1[-besN-2]*imExp(2*phi);
					}
					if (besN >= 2) {
						bjndd_(&n1,&arg,jn1,td1,td2);
						fn[0] = jn1[besN-2]*imExp(-2*phi);
						fn[1] = jn1[besN-1]*imExp(-phi);
						fn[2] = jn1[besN];
						fn[3] = jn1[besN+1]*imExp(phi);
						fn[4] = jn1[besN+2]*imExp(2*phi);
					}
					if (besN == -2) {
						bjndd_(&n1,&arg,jn1,td1,td2);
						fn[0] =  jn1[4]*imExp(-2*phi);
						fn[1] = -jn1[3]*imExp(-phi);
						fn[2] =  jn1[2];
						fn[3] = -jn1[1]*imExp(phi);
						fn[4] =  jn1[0]*imExp(2*phi);
					}
					if (besN == -1) {
						bjndd_(&n1,&arg,jn1,td1,td2);
						fn[0] = -jn1[3]*imExp(-2*phi);
						fn[1] =  jn1[2]*imExp(-phi);
						fn[2] = -jn1[1];
						fn[3] =  jn1[0]*imExp(phi);
						fn[4] =  jn1[1]*imExp(2*phi);
					}
					if (besN == 0) {
						bjndd_(&n1,&arg,jn1,td1,td2);
						fn[0] =  jn1[2]*imExp(-2*phi);
						fn[1] = -jn1[1]*imExp(-phi);
						fn[2] =  jn1[0];
						fn[3] =  jn1[1]*imExp(phi);
						fn[4] =  jn1[2]*imExp(2*phi);
					}
					if (besN == 1) {
						bjndd_(&n1,&arg,jn1,td1,td2);
						fn[0] = -jn1[1]*imExp(-2*phi);
						fn[1] =  jn1[0]*imExp(-phi);
						fn[2] =  jn1[1];
						fn[3] =  jn1[2]*imExp(phi);
						fn[4] =  jn1[3]*imExp(2*phi);
					}
				}
				/* TODO: all factors before f[n] should be calculated in InitBeam beforehand
				 * this will improve speed and allow robust calculation in the limit, e.g., of kt->0
				 */
				t1 = (((WaveNum*WaveNum+besKz*besKz)/2.*besM[0] +  WaveNum*besKz*besM[3])*fn[2] +
					  besKt*besKt/4.*((besM[0]+I*besM[1])*fn[0] + (besM[0]-I*besM[1])*fn[4])); // Ex
				t2 = (((WaveNum*WaveNum+besKz*besKz)/2.*besM[1] -  WaveNum*besKz*besM[2])*fn[2] +
					  I*besKt*besKt/4.*((besM[0]+I*besM[1])*fn[0] - (besM[0]-I*besM[1])*fn[4])); // Ey
				t3 = ((I*besKz*(besM[0]+I*besM[1]) + WaveNum*(besM[2]+I*besM[3]))*fn[1] -
					  (I*besKz*(besM[0]-I*besM[1]) - WaveNum*(besM[2]-I*besM[3]))*fn[3])*besKt/2.; // Ez

				cvMultScal_RVec(t1,ex,v1);
				cvMultScal_RVec(t2,ey,v2);
				cvMultScal_RVec(t3,prop,v3);
				cvAdd2Self(v1,v2,v3);
				cvMultScal_cmplx(ctemp,v1,b+j);
			}
			return;

#endif // !NO_FORTRAN		

		case B_READ:
			if (which==INCPOL_Y) fname=beam_fnameY;
			else fname=beam_fnameX; // which==INCPOL_X
			ReadField(fname,b);
			return;
	}
	LogError(ONE_POS,"Unknown type of incident beam (%d)",(int)beamtype);
	/* TO ADD NEW BEAM
	 * add a case above. Identifier ('B_...') should be defined inside 'enum beam' in const.h. This case should set
	 * complex vector 'b', describing the incident field in the particle reference frame. It is set inside the cycle for
	 * each dipole of the particle and is calculated using
	 * 1) 'DipoleCoord' - array of dipole coordinates;
	 * 2) 'prop' - propagation direction of the incident field;
	 * 3) 'ex' - direction of incident polarization;
	 * 4) 'ey' - complementary unity vector of polarization (orthogonal to both 'prop' and 'ex');
	 * 5) 'beam_center' - beam center in the particle reference frame (automatically calculated from 'beam_center_0'
	 *                    defined by '-beam_center' command line option).
	 * If the new beam type is compatible with '-surf', include here the corresponding code. For that you will need
	 * the variables, related to surface - see vars.c after "// related to a nearby surface".
	 * If you need temporary local variables (which are used only in this part of the code), either use 't1'-'t8' or
	 * define your own (with more informative names) in the beginning of this function.
	 */
}

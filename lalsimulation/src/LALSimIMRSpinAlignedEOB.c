/*
*  Copyright (C) 2011 Craig Robinson, Enrico Barausse, Yi Pan, Prayush Kumar
*  (minor changes)
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**
 * \author Craig Robinson, Yi Pan, Prayush Kumar
 *
 * \file
 *
 * \brief Functions for producing SEOBNRv1 waveforms for
 * spinning binaries, as described in
 * Taracchini et al. ( PRD 86, 024011 (2012), arXiv 1202.0790 ).
 * All equation numbers in this file refer to equations of this paper,
 * unless otherwise specified.
 *
 *
 * \brief Functions for producing SEOBNRv2 waveforms for
 * spinning binaries, as described in
 * Taracchini et al. ( arXiv 1311.2544 ).
 * 
 * 
 * \brief Functions for producing SEOBNRv3 waveforms for 
 * precessing binaries of spinning compact objects, as described
 * in Pan et al. ( arXiv 1307.6232 ).
 */

#include <math.h>
#include <complex.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKutta4.h>
#include <lal/SphericalHarmonics.h>
#include <gsl/gsl_sf_gamma.h>

#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"
#include "LALSimInspiralPrecess.h"

/* Include all the static function files we need */
#include "LALSimIMREOBHybridRingdown.c"
#include "LALSimIMREOBFactorizedWaveform.c"
#include "LALSimIMREOBNewtonianMultipole.c"
#include "LALSimIMREOBNQCCorrection.c"
#include "LALSimIMRSpinEOBInitialConditions.c"
#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinAlignedEOBHcapDerivative.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBFactorizedWaveformCoefficients.c"
#include "LALSimIMRSpinEOBFactorizedWaveform.c"
#include "LALSimIMRSpinEOBFactorizedFlux.c"


#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static int
XLALEOBSpinStopCondition(double UNUSED t,
                           const double values[],
                           double dvalues[],
                           void *funcParams
                          )
{

  SpinEOBParams *params = (SpinEOBParams *)funcParams;
  double omega_x, omega_y, omega_z, omega;
  double r2;

  omega_x = values[1]*dvalues[2] - values[2]*dvalues[1];
  omega_y = values[2]*dvalues[0] - values[0]*dvalues[2];
  omega_z = values[0]*dvalues[1] - values[1]*dvalues[0];

  r2 = values[0]*values[0] + values[1]*values[1] + values[2]*values[2];
  omega = sqrt( omega_x*omega_x + omega_y*omega_y + omega_z*omega_z )/r2;

  /* Terminate when omega reaches peak, and separation is < 6M */
  //if ( omega < params->eobParams->omega )
  if ( r2 < 16. && omega < params->eobParams->omega )
  {
    return 1;
  }

  params->eobParams->omega = omega;
  return GSL_SUCCESS;
}

UNUSED static int
XLALEOBSpinStopConditionBasedOnPR(double UNUSED t,
                           const double values[],
                           double dvalues[],
                           void *funcParams
                          )
{
  SpinEOBParams UNUSED *params = (SpinEOBParams *)funcParams;
  
  REAL8 r2, pDotr = 0;
  REAL8 p[3], r[3];

  memcpy( r, values, 3*sizeof(REAL8));
  memcpy( p, values+3, 3*sizeof(REAL8));
  
  r2 = inner_product(r,r);
  pDotr = inner_product( p, r ) / sqrt(r2);
  
  /* Terminate when omega reaches peak, and separation is < 6M */
  //if ( omega < params->eobParams->omega )
  if ( ( r2 < 36 && pDotr >= 0 ) || isnan( dvalues[3] ) || isnan (dvalues[2]) || isnan (dvalues[1]) || isnan (dvalues[0]) )
  {
    return 1;
  }

  return GSL_SUCCESS;
}


/**
 * Stopping condition for the regular resolution EOB orbital evolution
 * -- stop when reaching max orbital frequency in strong field.
 * At each test,
 * if omega starts to decrease, return 1 to stop evolution;
 * if not, update omega with current value and return GSL_SUCCESS to continue evolution.
 */
static int
XLALEOBSpinAlignedStopCondition(double UNUSED t,  /**< UNUSED */
                           const double values[], /**< dynamical variable values */
                           double dvalues[],      /**< dynamical variable time derivative values */
                           void *funcParams       /**< physical parameters */
                          )
{

  REAL8 omega, r;
  SpinEOBParams *params = (SpinEOBParams *)funcParams;

  r     = values[0];
  omega = dvalues[1];

  //if ( omega < params->eobParams->omega )
  if ( r < 6. && omega < params->eobParams->omega )
  {
    return 1;
  }

  params->eobParams->omega = omega;
  return GSL_SUCCESS;
}

/**
 * Stopping condition for the high resolution EOB orbital evolution
 * -- stop when reaching a minimum radius 0.3M out of the EOB horizon (Eqs. 9b, 37)
 * or when getting nan in any of the four ODE equations
 * At each test,
 * if conditions met, return 1 to stop evolution;
 * if not, return GSL_SUCCESS to continue evolution.
 */
static int
XLALSpinAlignedHiSRStopCondition(double UNUSED t,  /**< UNUSED */
                           const double values[], /**< dynamical variable values */
                           double dvalues[],      /**< dynamical variable time derivative values */
                           void *funcParams       /**< physical parameters */
                          )
{
  SpinEOBParams *params = (SpinEOBParams *)funcParams;
  REAL8 K, eta, chiK;
  REAL8 rshift = 0.6;
  eta  = params->eobParams->eta;
  chiK = params->sigmaKerr->data[2] / (1.-2.*eta);
  K = 1.4467 -  1.7152360250654402 * eta - 3.246255899738242 * eta * eta;

  if ( chiK < 0.8 ) rshift = 0.5;
  if ( chiK < 0.78) rshift = 0.475;
  if ( chiK < 0.72) rshift = 0.45;
  if ( chiK > -0.8 && chiK < 0.67 ) rshift = 0.35;

  if ( values[0] <= (1.+sqrt(1-params->a * params->a))*(1.-K*eta) + rshift+0.02 || isnan( dvalues[3] ) || isnan (dvalues[2]) || isnan (dvalues[1]) || isnan (dvalues[0]) )
  {
    return 1;
  }
  return GSL_SUCCESS;
}

/**
 * This function generates spin-aligned SEOBNRv1 waveforms h+ and hx.
 * Currently, only the h22 harmonic is available.
 * STEP 0) Prepare parameters, including pre-computed coefficients
 * for EOB Hamiltonian, flux and waveform
 * STEP 1) Solve for initial conditions
 * STEP 2) Evolve EOB trajectory until reaching the peak of orbital frequency
 * STEP 3) Step back in time by tStepBack and volve EOB trajectory again
 * using high sampling rate, stop at 0.3M out of the "EOB horizon".
 * STEP 4) Locate the peak of orbital frequency for NQC and QNM calculations
 * STEP 5) Calculate NQC correction using hi-sampling data
 * STEP 6) Calculate QNM excitation coefficients using hi-sampling data
 * STEP 7) Generate full inspiral waveform using desired sampling frequency
 * STEP 8) Generate full IMR modes -- attaching ringdown to inspiral
 * STEP 9) Generate full IMR hp and hx waveforms
 */
int XLALSimIMRSpinAlignedEOBWaveform(
        REAL8TimeSeries **hplus,     /**<< OUTPUT, +-polarization waveform */
        REAL8TimeSeries **hcross,    /**<< OUTPUT, x-polarization waveform */
        const REAL8     phiC,        /**<< coalescence orbital phase (rad) */ 
        REAL8           deltaT,      /**<< sampling time step */
        const REAL8     m1SI,        /**<< mass-1 in SI unit */ 
        const REAL8     m2SI,        /**<< mass-2 in SI unit */
        const REAL8     fMin,        /**<< starting frequency (Hz) */
        const REAL8     r,           /**<< distance in SI unit */
        const REAL8     inc,         /**<< inclination angle */
        const REAL8     spin1z,      /**<< z-component of spin-1, dimensionless */
        const REAL8     spin2z,       /**<< z-component of spin-2, dimensionless */
        UINT4           SpinAlignedEOBversion /**<< 1 for SEOBNRv1, 2 for SEOBNRv2 */
     )
{
  /* If the EOB version flag is neither 1 nor 2, exit */
  if (SpinAlignedEOBversion != 1 && SpinAlignedEOBversion != 2)
  {
    XLALPrintError("XLAL Error - %s: SEOBNR version flag incorrectly set to %u\n",
        __func__, SpinAlignedEOBversion);
    XLAL_ERROR( XLAL_EERR );
  }
  
  Approximant SpinAlignedEOBapproximant = (SpinAlignedEOBversion == 1) ? SEOBNRv1 : SEOBNRv2;

  /* If either spin > 0.6, model not available, exit */
  if ( SpinAlignedEOBversion == 1 && ( spin1z > 0.6 || spin2z > 0.6 ) )
  {
    XLALPrintError( "XLAL Error - %s: Component spin larger than 0.6!\nSEOBNRv1 is only available for spins in the range -1 < a/M < 0.6.\n", __func__);
    XLAL_ERROR( XLAL_EINVAL );
  }

  INT4 i;

  REAL8Vector *values = NULL;

  /* EOB spin vectors used in the Hamiltonian */
  REAL8Vector *sigmaStar = NULL;
  REAL8Vector *sigmaKerr = NULL;
  REAL8       a, tplspin;
  REAL8       chiS, chiA;

  /* Wrapper spin vectors used to calculate sigmas */
  REAL8Vector s1Vec, s1VecOverMtMt;
  REAL8Vector s2Vec, s2VecOverMtMt;
  REAL8       spin1[3] = {0, 0, spin1z};
  REAL8       spin2[3] = {0, 0, spin2z};
  REAL8       s1Data[3], s2Data[3], s1DataNorm[3], s2DataNorm[3];

  /* Parameters of the system */
  REAL8 m1, m2, mTotal, eta, mTScaled;
  REAL8 amp0;
  REAL8 sSub = 0.0;
  LIGOTimeGPS tc = LIGOTIMEGPSZERO;

  /* Dynamics of the system */
  REAL8Vector rVec, phiVec, prVec, pPhiVec;
  REAL8       omega, v, ham;

  /* Cartesian vectors needed to calculate Hamiltonian */
  REAL8Vector cartPosVec, cartMomVec;
  REAL8       cartPosData[3], cartMomData[3];

  /* Signal mode */
  COMPLEX16   hLM;
  REAL8Vector *sigReVec = NULL, *sigImVec = NULL;

  /* Non-quasicircular correction */
  EOBNonQCCoeffs nqcCoeffs;
  COMPLEX16      hNQC;
  REAL8Vector    *ampNQC = NULL, *phaseNQC = NULL;

  /* Ringdown freq used to check the sample rate */
  COMPLEX16Vector modefreqVec;
  COMPLEX16      modeFreq;

  /* Spin-weighted spherical harmonics */
  COMPLEX16  MultSphHarmP;
  COMPLEX16  MultSphHarmM;

  /* We will have to switch to a high sample rate for ringdown attachment */
  REAL8 deltaTHigh;
  UINT4 resampFac;
  UINT4 resampPwr;
  REAL8 resampEstimate;

  /* How far will we have to step back to attach the ringdown? */
  REAL8 tStepBack;
  INT4  nStepBack;

  /* Dynamics and details of the high sample rate part used to attach the ringdown */
  UINT4 hiSRndx;
  REAL8Vector timeHi, rHi, phiHi, prHi, pPhiHi;
  REAL8Vector *sigReHi = NULL, *sigImHi = NULL;
  REAL8Vector *omegaHi = NULL;

  /* Indices of peak frequency and final point */
  /* Needed to attach ringdown at the appropriate point */
  UINT4 peakIdx = 0, finalIdx = 0;

  /* (2,2) and (2,-2) spherical harmonics needed in (h+,hx) */
  REAL8 y_1, y_2, z1, z2;

  /* Variables for the integrator */
  ark4GSLIntegrator       *integrator = NULL;
  REAL8Array              *dynamics   = NULL;
  REAL8Array              *dynamicsHi = NULL;
  INT4                    retLen;
  REAL8  UNUSED           tMax;

  /* Accuracies of adaptive Runge-Kutta integrator */
  const REAL8 EPS_ABS = 1.0e-10;
  const REAL8 EPS_REL = 1.0e-9;

  /*
   * STEP 0) Prepare parameters, including pre-computed coefficients 
   *         for EOB Hamiltonian, flux and waveform  
   */

  /* Parameter structures containing important parameters for the model */
  SpinEOBParams           seobParams;
  SpinEOBHCoeffs          seobCoeffs;
  EOBParams               eobParams;
  FacWaveformCoeffs       hCoeffs;
  NewtonMultipolePrefixes prefixes;

  /* Initialize parameters */
  m1 = m1SI / LAL_MSUN_SI;
  m2 = m2SI / LAL_MSUN_SI;
  mTotal = m1 + m2;
  mTScaled = mTotal * LAL_MTSUN_SI;
  eta    = m1 * m2 / (mTotal*mTotal);

  amp0 = mTotal * LAL_MRSUN_SI / r;

  /* TODO: Insert potentially necessary checks on the arguments */

  /* Calculate the time we will need to step back for ringdown */
  tStepBack = 100. * mTScaled;
  nStepBack = ceil( tStepBack / deltaT );

  /* Calculate the resample factor for attaching the ringdown */
  /* We want it to be a power of 2 */
  /* If deltaT > Mtot/50, reduce deltaT by the smallest power of two for which deltaT < Mtot/50 */
  resampEstimate = 50. * deltaT / mTScaled;
  resampFac = 1;
  //resampFac = 1 << (UINT4)ceil(log2(resampEstimate));
  
  if ( resampEstimate > 1. )
  {
    resampPwr = (UINT4)ceil( log2( resampEstimate ) );
    while ( resampPwr-- )
    {
      resampFac *= 2u;
    }
  }
  

  /* Allocate the values vector to contain the initial conditions */
  /* Since we have aligned spins, we can use the 4-d vector as in the non-spin case */
  if ( !(values = XLALCreateREAL8Vector( 4 )) )
  {
    XLAL_ERROR( XLAL_ENOMEM );
  }
  memset ( values->data, 0, values->length * sizeof( REAL8 ));

  /* Set up structures and calculate necessary PN parameters */
  /* Unlike the general case, we only need to calculate these once */
  memset( &seobParams, 0, sizeof(seobParams) );
  memset( &seobCoeffs, 0, sizeof(seobCoeffs) );
  memset( &eobParams, 0, sizeof(eobParams) );
  memset( &hCoeffs, 0, sizeof( hCoeffs ) );
  memset( &prefixes, 0, sizeof( prefixes ) );

  /* Before calculating everything else, check sample freq is high enough */
  modefreqVec.length = 1;
  modefreqVec.data   = &modeFreq;

  if ( XLALSimIMREOBGenerateQNMFreqV2( &modefreqVec, m1, m2, spin1, spin2, 2, 2, 1, SpinAlignedEOBapproximant ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* If Nyquist freq < 220 QNM freq, exit */
  if ( deltaT > LAL_PI / creal(modeFreq) )
  {
    XLALPrintError( "XLAL Error - %s: Ringdown frequency > Nyquist frequency!\nAt present this situation is not supported.\n", __func__);
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EINVAL );
  }

  if ( !(sigmaStar = XLALCreateREAL8Vector( 3 )) )
  {
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  if ( !(sigmaKerr = XLALCreateREAL8Vector( 3 )) )
  {
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  seobParams.alignedSpins = 1;
  seobParams.tortoise     = 1;
  seobParams.sigmaStar    = sigmaStar;
  seobParams.sigmaKerr    = sigmaKerr;
  seobParams.seobCoeffs   = &seobCoeffs;
  seobParams.eobParams    = &eobParams;
  eobParams.hCoeffs       = &hCoeffs;
  eobParams.prefixes      = &prefixes;

  eobParams.m1  = m1;
  eobParams.m2  = m2;
  eobParams.eta = eta;

  s1Vec.length = s2Vec.length = 3;
  s1VecOverMtMt.length = s2VecOverMtMt.length = 3;
  s1Vec.data   = s1Data;
  s2Vec.data   = s2Data;
  s1VecOverMtMt.data   = s1DataNorm;
  s2VecOverMtMt.data   = s2DataNorm;

  /* copy the spins into the appropriate vectors, and scale them by the mass */
  memcpy( s1Data, spin1, sizeof( s1Data ) );
  memcpy( s2Data, spin2, sizeof( s2Data ) );
  memcpy( s1DataNorm, spin1, sizeof( s1DataNorm ) );
  memcpy( s2DataNorm, spin2, sizeof( s2DataNorm ) );

  /* Calculate chiS and chiA */

  chiS = 0.5 * (spin1[2] + spin2[2]);
  chiA = 0.5 * (spin1[2] - spin2[2]);

  for( i = 0; i < 3; i++ )
  {
    s1Data[i] *= m1*m1;
    s2Data[i] *= m2*m2;
  }
 for ( i = 0; i < 3; i++ )
  {
    s1DataNorm[i] = s1Data[i]/mTotal/mTotal;
    s2DataNorm[i] = s2Data[i]/mTotal/mTotal;
  } 
  seobParams.s1Vec    = &s1VecOverMtMt;
  seobParams.s2Vec    = &s2VecOverMtMt;
 
  cartPosVec.length = cartMomVec.length = 3;
  cartPosVec.data = cartPosData;
  cartMomVec.data = cartMomData;
  memset( cartPosData, 0, sizeof( cartPosData ) );
  memset( cartMomData, 0, sizeof( cartMomData ) );

  /* Populate the initial structures */
  if ( XLALSimIMRSpinEOBCalculateSigmaStar( sigmaStar, m1, m2, &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  if ( XLALSimIMRSpinEOBCalculateSigmaKerr( sigmaKerr, m1, m2, &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Calculate the value of a */
  /* XXX I am assuming that, since spins are aligned, it is okay to just use the z component XXX */
  /* TODO: Check this is actually the way it works in LAL */
  switch ( SpinAlignedEOBversion )
  {
     case 1:
       tplspin = 0.0;
       break;
     case 2:
       tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }
  /*for ( i = 0; i < 3; i++ )
  {
    a += sigmaKerr->data[i]*sigmaKerr->data[i];
  }
  a = sqrt( a );*/
  seobParams.a = a = sigmaKerr->data[2];
  seobParams.chi1 = spin1[2];
  seobParams.chi2 = spin2[2];

  /* Now compute the spinning H coefficients and store them in seobCoeffs */
  if ( XLALSimIMRCalculateSpinEOBHCoeffs( &seobCoeffs, eta, a, SpinAlignedEOBversion ) == XLAL_FAILURE )
  {    
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  if ( XLALSimIMREOBCalcSpinFacWaveformCoefficients( &hCoeffs, m1, m2, eta, tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  if ( XLALSimIMREOBComputeNewtonMultipolePrefixes( &prefixes, eobParams.m1, eobParams.m2 )
         == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /*
   * STEP 1) Solve for initial conditions
   */

  /* Set the initial conditions. For now we use the generic case */
  /* Can be simplified if spin-aligned initial conditions solver available. The cost of generic code is negligible though. */
  REAL8Vector *tmpValues = XLALCreateREAL8Vector( 14 );
  if ( !tmpValues )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  memset( tmpValues->data, 0, tmpValues->length * sizeof( REAL8 ) );

  /* We set inc zero here to make it easier to go from Cartesian to spherical coords */
  /* No problem setting inc to zero in solving spin-aligned initial conditions. */
  /* inc is not zero in generating the final h+ and hx */
  if ( XLALSimIMRSpinEOBInitialConditions( tmpValues, m1, m2, fMin, 0, s1Data, s2Data, &seobParams ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( tmpValues );
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /*fprintf( stderr, "ICs = %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", tmpValues->data[0], tmpValues->data[1], tmpValues->data[2],
      tmpValues->data[3], tmpValues->data[4], tmpValues->data[5], tmpValues->data[6], tmpValues->data[7], tmpValues->data[8],
      tmpValues->data[9], tmpValues->data[10], tmpValues->data[11] );*/

  /* Taken from Andrea's code */
/*  memset( tmpValues->data, 0, tmpValues->length*sizeof(tmpValues->data[0]));*/
#if 0
  tmpValues->data[0] = 19.9947984026;
  tmpValues->data[3] = -0.000433854158413;
  tmpValues->data[4] = 4.84217964546/tmpValues->data[0]; // q=1
#endif
#if 0
  tmpValues->data[0] = 19.9982539582;
  tmpValues->data[3] = -0.000390702473305;
  tmpValues->data[4] = 4.71107185264/tmpValues->data[0]; // q=1, chi1=chi2=0.98
#endif
#if 0
  tmpValues->data[0] = 19.996332305;
  tmpValues->data[3] = -0.000176807206312;
  tmpValues->data[4] = 4.84719922687/tmpValues->data[0]; // q=8
#endif
#if 0
  tmpValues->data[0] = 6.22645094958;
  tmpValues->data[3] = -0.00851784427559;
  tmpValues->data[4] = 3.09156589713/tmpValues->data[0]; // q=8 chi1=0.5 TEST DYNAMICS
#endif
#if 0
  tmpValues->data[0] = 19.9996712714;
  tmpValues->data[3] = -0.00016532905477;
  tmpValues->data[4] = 4.77661989696/tmpValues->data[0]; // q=8 chi1=0.5
#endif

  /* Now convert to Spherical */
  /* The initial conditions code returns Cartesian components of four vectors x, p, S1 and S2,
   * in the special case that the binary starts on the x-axis and the two spins are aligned
   * with the orbital angular momentum along the z-axis.
   * Therefore, in spherical coordinates the initial conditions are
   * r = x; phi = 0.; pr = px; pphi = r * py.
   */
  values->data[0] = tmpValues->data[0];
  values->data[1] = 0.;
  values->data[2] = tmpValues->data[3];
  values->data[3] = tmpValues->data[0] * tmpValues->data[4];

  //fprintf( stderr, "Spherical initial conditions: %e %e %e %e\n", values->data[0], values->data[1], values->data[2], values->data[3] );

  /*
   * STEP 2) Evolve EOB trajectory until reaching the peak of orbital frequency
   */

  /* Now we have the initial conditions, we can initialize the adaptive integrator */
  if (!(integrator = XLALAdaptiveRungeKutta4Init(4, XLALSpinAlignedHcapDerivative, XLALEOBSpinAlignedStopCondition, EPS_ABS, EPS_REL)))
  {
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  integrator->stopontestonly = 1;
  integrator->retries = 1;

  retLen = XLALAdaptiveRungeKutta4( integrator, &seobParams, values->data, 0., 20./mTScaled, deltaT/mTScaled, &dynamics );
  if ( retLen == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Set up pointers to the dynamics */
  rVec.length = phiVec.length = prVec.length = pPhiVec.length = retLen;
  rVec.data    = dynamics->data+retLen;
  phiVec.data  = dynamics->data+2*retLen;
  prVec.data   = dynamics->data+3*retLen;
  pPhiVec.data = dynamics->data+4*retLen;

  //printf( "We think we hit the peak at time %e\n", dynamics->data[retLen-1] );

  /* TODO : Insert high sampling rate / ringdown here */
  FILE *out = fopen( "saDynamics.dat", "w" );
  for ( i = 0; i < retLen; i++ )
  {
    fprintf( out, "%.16e %.16e %.16e %.16e %.16e\n", dynamics->data[i], rVec.data[i], phiVec.data[i], prVec.data[i], pPhiVec.data[i] );
  }
  fclose( out );

  /*
   * STEP 3) Step back in time by tStepBack and volve EOB trajectory again 
   *         using high sampling rate, stop at 0.3M out of the "EOB horizon".
   */

  /* Set up the high sample rate integration */
  hiSRndx = retLen - nStepBack;
  deltaTHigh = deltaT / (REAL8)resampFac;

  fprintf( stderr, "Stepping back %d points - we expect %d points at high SR\n", nStepBack, nStepBack*resampFac );
  fprintf( stderr, "Commencing high SR integration... from %.16e %.16e %.16e %.16e %.16e\n",
     (dynamics->data)[hiSRndx],rVec.data[hiSRndx], phiVec.data[hiSRndx], prVec.data[hiSRndx], pPhiVec.data[hiSRndx] );

  values->data[0] = rVec.data[hiSRndx];
  values->data[1] = phiVec.data[hiSRndx];
  values->data[2] = prVec.data[hiSRndx];
  values->data[3] = pPhiVec.data[hiSRndx];
  /* For HiSR evolution, we stop at a radius 0.3M from the deformed Kerr singularity, 
   * or when any derivative of Hamiltonian becomes nan */
  integrator->stop = XLALSpinAlignedHiSRStopCondition;

  retLen = XLALAdaptiveRungeKutta4( integrator, &seobParams, values->data, 0., 20./mTScaled, deltaTHigh/mTScaled, &dynamicsHi );
  if ( retLen == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  //fprintf( stderr, "We got %d points at high SR\n", retLen );

  /* Set up pointers to the dynamics */
  rHi.length = phiHi.length = prHi.length = pPhiHi.length = timeHi.length = retLen;
  timeHi.data = dynamicsHi->data;
  rHi.data    = dynamicsHi->data+retLen;
  phiHi.data  = dynamicsHi->data+2*retLen;
  prHi.data   = dynamicsHi->data+3*retLen;
  pPhiHi.data = dynamicsHi->data+4*retLen;

  out = fopen( "saDynamicsHi.dat", "w" );
  for ( i = 0; i < retLen; i++ )
  {
    fprintf( out, "%.16e %.16e %.16e %.16e %.16e\n", timeHi.data[i], rHi.data[i], phiHi.data[i], prHi.data[i], pPhiHi.data[i] );
  }
  fclose( out );

  /* Allocate the high sample rate vectors */
  sigReHi  = XLALCreateREAL8Vector( retLen + (UINT4)ceil( 20 / ( cimag(modeFreq) * deltaTHigh )) );
  sigImHi  = XLALCreateREAL8Vector( retLen + (UINT4)ceil( 20 / ( cimag(modeFreq) * deltaTHigh )) );
  omegaHi  = XLALCreateREAL8Vector( retLen + (UINT4)ceil( 20 / ( cimag(modeFreq) * deltaTHigh )) );
  ampNQC   = XLALCreateREAL8Vector( retLen );
  phaseNQC = XLALCreateREAL8Vector( retLen );

  if ( !sigReHi || !sigImHi || !omegaHi || !ampNQC || !phaseNQC )
  {
    XLAL_ERROR( XLAL_ENOMEM );
  }

  memset( sigReHi->data, 0, sigReHi->length * sizeof( sigReHi->data[0] ));
  memset( sigImHi->data, 0, sigImHi->length * sizeof( sigImHi->data[0] ));

  /* Populate the high SR waveform */
  REAL8 omegaOld = 0.0;
  INT4  phaseCounter = 0;
  for ( i = 0; i < retLen; i++ )
  {
    values->data[0] = rHi.data[i];
    values->data[1] = phiHi.data[i];
    values->data[2] = prHi.data[i];
    values->data[3] = pPhiHi.data[i];

    omegaHi->data[i] = omega = XLALSimIMRSpinAlignedEOBCalcOmega( values->data, &seobParams );
    v = cbrt( omega );

    /* Calculate the value of the Hamiltonian */
    cartPosVec.data[0] = values->data[0];
    cartMomVec.data[0] = values->data[2];
    cartMomVec.data[1] = values->data[3] / values->data[0];

    ham = XLALSimIMRSpinEOBHamiltonian( eta, &cartPosVec, &cartMomVec, &s1VecOverMtMt, &s2VecOverMtMt, sigmaKerr, sigmaStar, seobParams.tortoise, &seobCoeffs );

    if ( XLALSimIMRSpinEOBGetSpinFactorizedWaveform( &hLM, values, v, ham, 2, 2, &seobParams )
           == XLAL_FAILURE )
    {
      /* TODO: Clean-up */
      XLAL_ERROR( XLAL_EFUNC );
    }

    ampNQC->data[i]  = cabs( hLM );
    sigReHi->data[i] = (REAL4)(amp0 * creal(hLM));
    sigImHi->data[i] = (REAL4)(amp0 * cimag(hLM));
    phaseNQC->data[i]= carg( hLM ) + phaseCounter * LAL_TWOPI;

    if ( i && phaseNQC->data[i] > phaseNQC->data[i-1] )
    {
      phaseCounter--;
      phaseNQC->data[i] -= LAL_TWOPI;
    }

    if ( omega <= omegaOld && !peakIdx )
    {
      //printf( "Have we got the peak? omegaOld = %.16e, omega = %.16e\n", omegaOld, omega );
      peakIdx = i;
    }
    omegaOld = omega;
  }
  //printf( "We now think the peak is at %d\n", peakIdx );
  finalIdx = retLen - 1;

  /*
   * STEP 4) Locate the peak of orbital frequency for NQC and QNM calculations
   */

  /* Stuff to find the actual peak time */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;
  REAL8 omegaDeriv1; //, omegaDeriv2;
  REAL8 time1, time2;
  REAL8 timePeak, timewavePeak = 0., omegaDerivMid;
  REAL8 sigAmpSqHi = 0., oldsigAmpSqHi = 0.;
  INT4  peakCount = 0;

  spline = gsl_spline_alloc( gsl_interp_cspline, retLen );
  acc    = gsl_interp_accel_alloc();

  time1 = dynamicsHi->data[peakIdx];

  gsl_spline_init( spline, dynamicsHi->data, omegaHi->data, retLen );
  omegaDeriv1 = gsl_spline_eval_deriv( spline, time1, acc );
  if ( omegaDeriv1 > 0. )
  {
    time2 = dynamicsHi->data[peakIdx+1];
    //omegaDeriv2 = gsl_spline_eval_deriv( spline, time2, acc );
  }
  else
  {
    //omegaDeriv2 = omegaDeriv1;
    time2 = time1;
    time1 = dynamicsHi->data[peakIdx-1];
    peakIdx--;
    omegaDeriv1 = gsl_spline_eval_deriv( spline, time1, acc );
  }

  do
  {
    timePeak = ( time1 + time2 ) / 2.;
    omegaDerivMid = gsl_spline_eval_deriv( spline, timePeak, acc );

    if ( omegaDerivMid * omegaDeriv1 < 0.0 )
    {
      //omegaDeriv2 = omegaDerivMid;
      time2 = timePeak;
    }
    else
    {
      omegaDeriv1 = omegaDerivMid;
      time1 = timePeak;
    }
  }
  while ( time2 - time1 > 1.0e-5 );

  /*gsl_spline_free( spline );
  gsl_interp_accel_free( acc );
  */

  XLALPrintInfo( "Estimation of the peak is now at time %.16e\n", timePeak );

  /* Having located the peak of orbital frequency, we set time and phase of coalescence */
  XLALGPSAdd( &tc, -mTScaled * (dynamics->data[hiSRndx] + timePeak));
  gsl_spline_init( spline, dynamicsHi->data, phiHi.data, retLen );
  sSub = gsl_spline_eval( spline, timePeak, acc ) - phiC;
  gsl_spline_free( spline );
  gsl_interp_accel_free( acc );
  /* Apply phiC to hi-sampling waveforms */
  REAL8 thisReHi, thisImHi;
  REAL8 csSub2 = cos(2.0 * sSub);
  REAL8 ssSub2 = sin(2.0 * sSub);
  for ( i = 0; i < retLen; i++)
  {
    thisReHi = sigReHi->data[i];
    thisImHi = sigImHi->data[i];
    sigReHi->data[i] =   thisReHi * csSub2 - thisImHi * ssSub2;
    sigImHi->data[i] =   thisReHi * ssSub2 + thisImHi * csSub2; 
  }

  /*
   * STEP 5) Calculate NQC correction using hi-sampling data
   */

  /* Calculate nonspin and amplitude NQC coefficients from fits and interpolation table */
  switch ( SpinAlignedEOBversion )
  {
     case 1:
       if ( XLALSimIMRGetEOBCalibratedSpinNQC( &nqcCoeffs, 2, 2, eta, a ) == XLAL_FAILURE )
       {
         XLAL_ERROR( XLAL_EFUNC );
       }
       break;
     case 2:
       if ( XLALSimIMRGetEOBCalibratedSpinNQC3D( &nqcCoeffs, 2, 2, m1, m2, a, chiA ) == XLAL_FAILURE )
       {
         XLAL_ERROR( XLAL_EFUNC );
       }
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }

  /* Calculate phase NQC coefficients */
  if ( XLALSimIMRSpinEOBCalculateNQCCoefficients( ampNQC, phaseNQC, &rHi, &prHi, omegaHi,
          2, 2, timePeak, deltaTHigh/mTScaled, m1, m2, a, chiA, chiS, &nqcCoeffs, SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Calculate the time of amplitude peak. Despite the name, this is in fact the shift in peak time from peak orb freq time */
  switch ( SpinAlignedEOBversion )
  {
     case 1:
     timewavePeak = XLALSimIMREOBGetNRSpinPeakDeltaT(2, 2, eta,  a);
       break;
     case 2:
     timewavePeak = 0.0;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }

  /* Apply to the high sampled part */
  //out = fopen( "saWavesHi.dat", "w" );
  for ( i = 0; i < retLen; i++ )
  {
    values->data[0] = rHi.data[i];
    values->data[1] = phiHi.data[i] - sSub;
    values->data[2] = prHi.data[i];
    values->data[3] = pPhiHi.data[i];

    if ( XLALSimIMREOBNonQCCorrection( &hNQC, values, omegaHi->data[i], &nqcCoeffs ) == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }

    hLM = sigReHi->data[i];
    hLM += I * sigImHi->data[i];
    //fprintf( out, "%.16e %.16e %.16e %.16e %.16e\n", timeHi.data[i], creal(hLM), cimag(hLM), creal(hNQC), cimag(hNQC) );

    hLM *= hNQC;
    sigReHi->data[i] = (REAL4) creal(hLM);
    sigImHi->data[i] = (REAL4) cimag(hLM);
    sigAmpSqHi = creal(hLM)*creal(hLM)+cimag(hLM)*cimag(hLM);
    if (sigAmpSqHi < oldsigAmpSqHi && peakCount == 0 && (i-1)*deltaTHigh/mTScaled < timePeak - timewavePeak) 
    {
      timewavePeak = (i-1)*deltaTHigh/mTScaled;
      peakCount += 1;
    }
    oldsigAmpSqHi = sigAmpSqHi;
  }
  //fclose(out);
  /*printf("NQCs entering hNQC: %.16e, %.16e, %.16e, %.16e, %.16e, %.16e\n", nqcCoeffs.a1, nqcCoeffs.a2,nqcCoeffs.a3, nqcCoeffs.a3S, nqcCoeffs.a4, nqcCoeffs.a5 );
  printf("NQCs entering hNQC: %.16e, %.16e, %.16e, %.16e\n", nqcCoeffs.b1, nqcCoeffs.b2,nqcCoeffs.b3, nqcCoeffs.b4 );*/
  if (timewavePeak < 1.0e-16 || peakCount == 0)
  {
    //printf("YP::warning: could not locate mode peak, use calibrated time shift of amplitude peak instead.\n");
    /* NOTE: instead of looking for the actual peak, use the calibrated value,    */
    /*       ignoring the error in using interpolated NQC instead of iterated NQC */
    timewavePeak = timePeak - timewavePeak;
  }
  
  /*
   * STEP 6) Calculate QNM excitation coefficients using hi-sampling data
   */

  /*out = fopen( "saInspWaveHi.dat", "w" );
  for ( i = 0; i < retLen; i++ )
  {
    fprintf( out, "%.16e %.16e %.16e\n", timeHi.data[i], sigReHi->data[i], sigImHi->data[i] );
  }
  fclose( out );*/
  
  /* Attach the ringdown at the time of amplitude peak */
  REAL8 combSize = 7.5; /* Eq. 34 */
  REAL8 chi = (spin1[2] + spin2[2]) / 2. + (spin1[2] - spin2[2]) / 2. * sqrt(1. - 4. * eta) / (1. - 2. * eta);
  
  /* Modify the combsize for SEOBNRv2 */
  /* If chi1=chi2=0, comb = 11. if chi < 0.8, comb = 12. if chi >= 0.8, comb =
   * 13.5 */
  if( SpinAlignedEOBversion == 2 )
  {
    combSize = (spin1[2] == 0. && spin2[2] == 0.) ? 11. : (( chi >= 0.8 ) ? 13.5 : 12.);
  }

  REAL8 timeshiftPeak;
  timeshiftPeak = timePeak - timewavePeak;
  if ( SpinAlignedEOBversion == 2)
  {
    timeshiftPeak = (timePeak - timewavePeak) > 0. ? (timePeak - timewavePeak) : 0.;
  }

  /*printf("YP::timePeak and timewavePeak: %.16e and %.16e\n",timePeak,timewavePeak);
  printf("YP::timeshiftPeak and combSize: %.16e and %.16e\n",timeshiftPeak,combSize);
  printf("PK::chi and SpinAlignedEOBversion: %.16e and %u\n\n", chi,SpinAlignedEOBversion);*/

  REAL8Vector *rdMatchPoint = XLALCreateREAL8Vector( 3 );
  if ( !rdMatchPoint )
  {
    XLAL_ERROR( XLAL_ENOMEM );
  }
  
  if ( combSize > timePeak - timeshiftPeak )
  {
    XLALPrintError( "The comb size looks to be too big!!!\n" );
  }

  rdMatchPoint->data[0] = combSize < timePeak - timeshiftPeak ? timePeak - timeshiftPeak - combSize : 0;
  rdMatchPoint->data[1] = timePeak - timeshiftPeak;
  rdMatchPoint->data[2] = dynamicsHi->data[finalIdx];
  printf("YP::comb range: %f, %f\n",rdMatchPoint->data[0],rdMatchPoint->data[1]);
  if ( XLALSimIMREOBHybridAttachRingdown( sigReHi, sigImHi, 2, 2,
              deltaTHigh, m1, m2, spin1[0], spin1[1], spin1[2], spin2[0], spin2[1], spin2[2],
              &timeHi, rdMatchPoint, SpinAlignedEOBapproximant )
          == XLAL_FAILURE ) 
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /*
   * STEP 7) Generate full inspiral waveform using desired sampling frequency
   */

  /* Now create vectors at the correct sample rate, and compile the complete waveform */
  sigReVec = XLALCreateREAL8Vector( rVec.length + ceil( sigReHi->length / resampFac ) );
  sigImVec = XLALCreateREAL8Vector( sigReVec->length );

  memset( sigReVec->data, 0, sigReVec->length * sizeof( REAL8 ) );
  memset( sigImVec->data, 0, sigImVec->length * sizeof( REAL8 ) );
 
  /* Generate full inspiral waveform using desired sampling frequency */
  /* TODO - Check vectors were allocated */
  for ( i = 0; i < (INT4)rVec.length; i++ )
  {
    values->data[0] = rVec.data[i];
    values->data[1] = phiVec.data[i] - sSub;
    values->data[2] = prVec.data[i];
    values->data[3] = pPhiVec.data[i];

    omega = XLALSimIMRSpinAlignedEOBCalcOmega( values->data, &seobParams );
    v = cbrt( omega );

    /* Calculate the value of the Hamiltonian */
    cartPosVec.data[0] = values->data[0];
    cartMomVec.data[0] = values->data[2];
    cartMomVec.data[1] = values->data[3] / values->data[0];

    ham = XLALSimIMRSpinEOBHamiltonian( eta, &cartPosVec, &cartMomVec, &s1VecOverMtMt, &s2VecOverMtMt, sigmaKerr, sigmaStar, seobParams.tortoise, &seobCoeffs );

    if ( XLALSimIMRSpinEOBGetSpinFactorizedWaveform( &hLM, values, v, ham, 2, 2, &seobParams )
           == XLAL_FAILURE )
    {
      /* TODO: Clean-up */
      XLAL_ERROR( XLAL_EFUNC );
    }

    if ( XLALSimIMREOBNonQCCorrection( &hNQC, values, omega, &nqcCoeffs ) == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }

    hLM *= hNQC;

    sigReVec->data[i] = amp0 * creal(hLM);
    sigImVec->data[i] = amp0 * cimag(hLM);
  }

  /*
   * STEP 8) Generate full IMR modes -- attaching ringdown to inspiral
   */

  /* Attach the ringdown part to the inspiral */
  for ( i = 0; i < (INT4)(sigReHi->length / resampFac); i++ )
  {
    sigReVec->data[i+hiSRndx] = sigReHi->data[i*resampFac];
    sigImVec->data[i+hiSRndx] = sigImHi->data[i*resampFac];
  }

  /*
   * STEP 9) Generate full IMR hp and hx waveforms
   */
  
  /* For now, let us just try to create a waveform */
  REAL8TimeSeries *hPlusTS  = XLALCreateREAL8TimeSeries( "H_PLUS", &tc, 0.0, deltaT, &lalStrainUnit, sigReVec->length );
  REAL8TimeSeries *hCrossTS = XLALCreateREAL8TimeSeries( "H_CROSS", &tc, 0.0, deltaT, &lalStrainUnit, sigImVec->length );

  /* TODO change to using XLALSimAddMode function to combine modes */
  /* For now, calculate -2Y22 * h22 + -2Y2-2 * h2-2 directly (all terms complex) */
  /* Compute spin-weighted spherical harmonics and generate waveform */
  REAL8 coa_phase = 0.0;

  MultSphHarmP = XLALSpinWeightedSphericalHarmonic( inc, coa_phase, -2, 2, 2 );
  MultSphHarmM = XLALSpinWeightedSphericalHarmonic( inc, coa_phase, -2, 2, -2 );

  y_1 =   creal(MultSphHarmP) + creal(MultSphHarmM);
  y_2 =   cimag(MultSphHarmM) - cimag(MultSphHarmP);
  z1 = - cimag(MultSphHarmM) - cimag(MultSphHarmP);
  z2 =   creal(MultSphHarmM) - creal(MultSphHarmP);

  for ( i = 0; i < (INT4)sigReVec->length; i++ )
  {
    REAL8 x1 = sigReVec->data[i];
    REAL8 x2 = sigImVec->data[i];

    hPlusTS->data->data[i]  = (x1 * y_1) + (x2 * y_2);
    hCrossTS->data->data[i] = (x1 * z1) + (x2 * z2);
  }

  /* Point the output pointers to the relevant time series and return */
  (*hplus)  = hPlusTS;
  (*hcross) = hCrossTS;

  /* Free memory */
  XLALDestroyREAL8Vector( tmpValues );
  XLALDestroyREAL8Vector( sigmaKerr );
  XLALDestroyREAL8Vector( sigmaStar );
  XLALDestroyREAL8Vector( values );
  XLALDestroyREAL8Vector( rdMatchPoint );
  XLALDestroyREAL8Vector( ampNQC );
  XLALDestroyREAL8Vector( phaseNQC );
  XLALDestroyREAL8Vector( sigReVec );
  XLALDestroyREAL8Vector( sigImVec );
  XLALAdaptiveRungeKutta4Free( integrator );
  XLALDestroyREAL8Array( dynamics );
  XLALDestroyREAL8Array( dynamicsHi );
  XLALDestroyREAL8Vector( sigReHi );
  XLALDestroyREAL8Vector( sigImHi );
  XLALDestroyREAL8Vector( omegaHi );

  return XLAL_SUCCESS;
}

/** ********************************************************************
 *  THE FOLLOWING HAS FUNCTIONS FOR THE PRECESSING EOB MODEL
 *  ********************************************************************
 * */
/**
 * This function generates precessing spinning SEOBNRv3 waveforms h+ and hx.
 * Currently, only the h22 harmonic is available.
 * STEP 0) Prepare parameters, including pre-computed coefficients
 */

int XLALSimIMRSpinEOBWaveform(
        REAL8TimeSeries **hplus,
        REAL8TimeSeries **hcross,
        //LIGOTimeGPS     *tc,
        const REAL8     UNUSED phiC,
        const REAL8     deltaT,
        const REAL8     m1SI,
        const REAL8     m2SI,
        const REAL8     fMin,
        const REAL8     r,
        const REAL8     inc,
        const REAL8     INspin1[],
        const REAL8     INspin2[]
     )
{
  int importDynamicsAndGetDerivatives = 0;
  int debugPK = 0;

  INT4 i, k;
  UINT4 j;
  INT4 UNUSED status;
  LIGOTimeGPS tc = LIGOTIMEGPSZERO;
  
  /* SEOBNRv3 model */
  Approximant spinEOBApproximant = SEOBNRv3;
  /* Fix the underlying aligned spin EOB model */
  INT4 SpinAlignedEOBversion = 1;

  /* Vector to store the initial parameters */
  REAL8 spin1[3], spin2[3];
  memcpy( spin1, INspin1, 3*sizeof(REAL8));
  memcpy( spin2, INspin2, 3*sizeof(REAL8));
  
  REAL8Vector *values = NULL;
  /* Allocate the values vector to contain the ICs */
  /* For this model, it contains 12 dynamical variables: */
  /* values[0-2]  - x (Cartesian separation vector) */
  /* values[3-5]  - p (Cartesian momentum) */
  /* values[6-8]  - spin of body 1 */
  /* values[9-11] - spin of body 2 */
  if ( !(values = XLALCreateREAL8Vector( 14 )) )
  {
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  memset( values->data, 0, values->length * sizeof( REAL8 ));

  /* Vector to contain derivatives of ICs */
  REAL8Vector *dvalues = NULL; 
  if ( !(dvalues = XLALCreateREAL8Vector( 14 )) )
  {
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  REAL8       rdotvec[3];
  
  /* EOB spin vectors used in the Hamiltonian */
  REAL8Vector *sigmaStar = NULL;
  REAL8Vector *sigmaKerr = NULL;
  REAL8       a, tplspin;
  REAL8       chiS, chiA;
  REAL8       spinNQC;

  /* Spins not scaled by the mass */
  REAL8 mSpin1[3], mSpin2[3];
  
  /* Wrapper spin vectors used to calculate sigmas */
  REAL8Vector s1Vec, s1VecOverMtMt;
  REAL8Vector s2Vec, s2VecOverMtMt;
  REAL8       s1Data[3], s2Data[3], s1DataNorm[3], s2DataNorm[3];

  /* Parameters of the system */
  REAL8 m1, m2, mTotal, eta, mTScaled;
  REAL8 UNUSED amp0, amp;
  REAL8 UNUSED sSub = 0.0;

  /* Dynamics of the system */
  REAL8Vector UNUSED rVec, phiVec, prVec, pPhiVec;
  REAL8       omega, v, ham;

  /* Cartesian vectors needed to calculate Hamiltonian */
  REAL8Vector cartPosVec, cartMomVec;
  REAL8       cartPosData[3], cartMomData[3];
  REAL8 rcrossrdotNorm, rvec[3], rcrossrdot[3], s1dotLN, s2dotLN;; 
  REAL8 UNUSED rcrossp[3], rcrosspMag, s1dotL, s2dotL;

  /* Signal mode */
  COMPLEX16   hLM;
  REAL8Vector UNUSED *sigReVec = NULL, *sigImVec = NULL;

  /* Non-quasicircular correction */
  EOBNonQCCoeffs nqcCoeffs;
  COMPLEX16      UNUSED   hNQC;
  REAL8Vector    UNUSED   *ampNQC = NULL, *phaseNQC = NULL;

  /* Ringdown freq used to check the sample rate */
  COMPLEX16Vector UNUSED modefreqVec;
  COMPLEX16       UNUSED modeFreq;

  /* Spin-weighted spherical harmonics */
  COMPLEX16 UNUSED MultSphHarmP;
  COMPLEX16 UNUSED MultSphHarmM;

  /* We will have to switch to a high sample rate for ringdown attachment */
  //REAL8 deltaTHigh;
  //UINT4 resampFac;
  //UINT4 resampPwr;
  //REAL8 resampEstimate;

  /* How far will we have to step back to attach the ringdown? */
  //REAL8 tStepBack;
  //INT4  nStepBack;

  /* Dynamics and details of the high sample rate part used to attach the ringdown */
  //UINT4 hiSRndx;
  //REAL8Vector timeHi, rHi, phiHi, prHi, pPhiHi;
  //REAL8Vector *sigReHi = NULL, *sigImHi = NULL;
  //REAL8Vector *omegaHi = NULL;

  /* Indices of peak frequency and final point */
  /* Needed to attach ringdown at the appropriate point */
  UINT4 UNUSED peakIdx = 0, finalIdx = 0;

  /* (2,2) and (2,-2) spherical harmonics needed in (h+,hx) */
  REAL8 UNUSED y_1, y_2, z1, z2;

  /* Parameter structures containing important parameters for the model */
  SpinEOBParams           seobParams;
  SpinEOBHCoeffs          seobCoeffs;
  EOBParams               eobParams;
  FacWaveformCoeffs       hCoeffs;
  NewtonMultipolePrefixes prefixes;

  /* Set up structures and calculate necessary PN parameters */
  /* Due to precession, these need to get calculated in every step */
  /* TODO: Only calculate non-spinning parts once */
  memset( &seobParams, 0, sizeof(seobParams) );
  memset( &seobCoeffs, 0, sizeof(seobCoeffs) );
  memset( &eobParams, 0, sizeof(eobParams) );
  memset( &hCoeffs, 0, sizeof( hCoeffs ) );
  memset( &prefixes, 0, sizeof( prefixes ) );

  /* Variables for the integrator */
  ark4GSLIntegrator       *integrator = NULL;
  REAL8Array              *dynamics   = NULL;
  //REAL8Array              *dynamicsHi = NULL;
  INT4                    retLen;
  INT4                    retLenRDPatch;
  REAL8  UNUSED           tMax;

  /* Accuracies of adaptive Runge-Kutta integrator */
  const REAL8 EPS_ABS = 1.0e-9;
  const REAL8 EPS_REL = 1.0e-8;

  if ( !(sigmaStar = XLALCreateREAL8Vector( 3 )) )
  {
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  if ( !(sigmaKerr = XLALCreateREAL8Vector( 3 )) )
  {
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  /* Initialize parameters */
  //m1 = 25.0;//m1SI / LAL_MSUN_SI;
  //m2 = 5.0;//m2SI / LAL_MSUN_SI;
  m1 = m1SI / LAL_MSUN_SI;
  m2 = m2SI / LAL_MSUN_SI;
  mTotal = m1 + m2;
  mTScaled = mTotal * LAL_MTSUN_SI;
  eta    = m1 * m2 / (mTotal*mTotal);
  printf("YP: m1 = %f, m2 = %f\n",m1,m2);
  amp0 = mTotal * LAL_MRSUN_SI / r;
  //amp0 = 4. * mTotal * LAL_MRSUN_SI * eta / r;
  
#if 1  
  /*values->data[0] = 19.9;
  values->data[1] = -1.04337949716e-22;
  values->data[2] = -5.04250029413e-19;
  values->data[3] = -0.000201238161271;
  values->data[4] = 0.24386594182;
  values->data[5] = -0.000100920067298;
  values->data[6] = 0.26796875;
  values->data[7] = 0.30625;
  values->data[8] = -0.0765625000004;
  values->data[9] = -0.009375;
  values->data[10] = -0.00156249999995;
  values->data[11] = -0.00156250000004;*/

  /*values->data[0] = 15.87;
  values->data[1] = 0.;
  values->data[2] = 0.;
  values->data[3] = -0.0005229624453230818;
  values->data[4] = 0.2779908426296164 ;
  values->data[5] = -6.328418862992827e-05;
  values->data[6] = -0.2704529501882914 * (mTotal/m1) * (mTotal/m1);
  values->data[7] = -0.2168021314138335 * (mTotal/m1) * (mTotal/m1);
  values->data[8] = 0.001330438577632606 * (mTotal/m1) * (mTotal/m1);
  values->data[9] = 0.;
  values->data[10] = 0.;
  values->data[11] = 0.;*/
  
  values->data[0] = 25;
  values->data[1] = -2.738042987012824e-23;
  values->data[2] = -2.895313259628355e-19;
  values->data[3] = -0.0001185485542182621;
  values->data[4] = 0.2118058597394067;
  values->data[5] = -4.006015947416696e-05;
  values->data[6] = -0.4122563355465681 * (mTotal/m1) * (mTotal/m1);
  values->data[7] = -0.3304754197472309 * (mTotal/m1) * (mTotal/m1);
  values->data[8] = 0.1716761079860819 * (mTotal/m1) * (mTotal/m1);
  values->data[9] = 0.01348361657291579 * (mTotal/m2) * (mTotal/m2);
  values->data[10] = 2.117582368135751e-22 * (mTotal/m2) * (mTotal/m2);
  values->data[11] = 0.009796420871541223 * (mTotal/m2) * (mTotal/m2);

  /* These ones did not work with the cpp code -- underflow error
  values->data[0] = 29.9;
  values->data[1] = 0.;
  values->data[2] = 0.;
  values->data[3] = -6.95123137564155e-05;
  values->data[4] = 0.1925025010802479;
  values->data[5] = -3.168712384234343e-05;
  values->data[6] = -0.433471963530112 * (mTotal/m1) * (mTotal/m1);
  values->data[7] = -0.3474824199034985 * (mTotal/m1) * (mTotal/m1);
  values->data[8] = 3.40168431617327e-17 * (mTotal/m1) * (mTotal/m1);
  values->data[9] = 0.01348361657291579 * (mTotal/m2) * (mTotal/m2);
  values->data[10] = 0.;
  values->data[11] = 0.009796420871541223 * (mTotal/m2) * (mTotal/m2);*/
  
  /*values->data[0] = 15.87;
  values->data[1] = 0.;
  values->data[2] = 0.;
  values->data[3] = -0.000521675194648;
  values->data[4] = 0.278174373488;
  values->data[5] = -0.00012666165246;
  values->data[6] = -0.270452950188 * (mTotal/m1) * (mTotal/m1);
  values->data[7] = -0.216802131414 * (mTotal/m1) * (mTotal/m1);
  values->data[8] = 0.00133043857763 * (mTotal/m1) * (mTotal/m1);
  values->data[9] = 0.;
  values->data[10] = 0.;
  values->data[11] = 0.;*/
    
  /*values->data[0] = 7.;
  values->data[1] = 0.;
  values->data[2] = 0.;
  values->data[3] = -0.01181688738719029;
  values->data[4] = 0.4843132461494007;
  values->data[5] = -0.003144147080494646;
  values->data[6] = -0.2704529501882914 * (mTotal/m1) * (mTotal/m1);
  values->data[7] = -0.2168021314138335 * (mTotal/m1) * (mTotal/m1);
  values->data[8] = 0.001330438577632606 * (mTotal/m1) * (mTotal/m1);
  values->data[9] = 0.;
  values->data[10] = 0.;
  values->data[11] = 0.01666666666666667 * (mTotal/m2) * (mTotal/m2);*/

#if 0
  values->data[0] = 3.5;
  values->data[1] = 4.949747468305832;
  values->data[2] = -3.5;
  values->data[3] = -0.2502883144899609;
  values->data[4] = 0.3341053793667105;
  values->data[5] = 0.2458418190466291;
  values->data[6] = -0.03727389203571741;// * (mTotal/m1) * (mTotal/m1);
  values->data[7] = -0.49613957621055016;// * (mTotal/m1) * (mTotal/m1);
  values->data[8] = 0.03998328699948254;// * (mTotal/m1) * (mTotal/m1);
  values->data[9] = 0.42426406871192845; //* (mTotal/m2) * (mTotal/m2);
  values->data[10] = 0.01 * (mTotal/m2) * (mTotal/m2);
  values->data[11] = 0.42426406871192845;// * (mTotal/m2) * (mTotal/m2);
#endif
    

  /*values->data[0] = 0.130208309399131;
  values->data[1] = -7.60959058900954;
  values->data[2] = -2.45855499735253;
  values->data[3] = 0.430218830242973;
  values->data[4] = 0.04291632558714;
  values->data[5] = -0.0891018381333908;
  values->data[6] = -0.337735482229554;
  values->data[7] = -0.0292783444192628;
  values->data[8] = 0.0722997424189602;
  values->data[9] = 0.;
  values->data[10] = 0.;
  values->data[11] = 0.;*/
#endif
  
  for( i = 0; i < 3; i++ )
  {
    /* Store the dimensionless chi vector, i.e. \vec{S}_i/m_i^2 */
    spin1[i] = values->data[i+6];
    spin2[i] = values->data[i+9];
    
    /* Spins being evolved, i.e. \vec{S}_i/M^2 */
    values->data[i+6] *= m1*m1/(mTotal*mTotal);
    values->data[i+9] *= m2*m2/(mTotal*mTotal);
  }

  /* TODO: Insert potentially necessary checks on the arguments */

  /* Calculate the time we will need to step back for ringdown */
  //tStepBack = 100. * mTScaled;
  //nStepBack = ceil( tStepBack / deltaT );

  /* Calculate the resample factor for attaching the ringdown */
  /* We want it to be a power of 2 */
  /* If deltaT > Mtot/50, reduce deltaT by the smallest power of two for which deltaT < Mtot/50 */
  //resampEstimate = 50. * deltaT / mTScaled;
  //resampFac = 1;
  //resampFac = 1 << (UINT4)ceil(log2(resampEstimate));
  /*
  if ( resampEstimate > 1. )
  {
    resampPwr = (UINT4)ceil( log2( resampEstimate ) );
    while ( resampPwr-- )
    {
      resampFac *= 2u;
    }
  }*/
    

  /* Wrapper spin vectors used to calculate sigmas */
  s1VecOverMtMt.length = s2VecOverMtMt.length = 3;
  s1VecOverMtMt.data   = s1DataNorm;
  s2VecOverMtMt.data   = s2DataNorm;

  s1Vec.length = s2Vec.length = 3;
  s1Vec.data   = s1Data;
  s2Vec.data   = s2Data;

  memcpy( s1Data, spin1, sizeof(s1Data) );
  memcpy( s2Data, spin2, sizeof(s2Data) );
  memcpy( s1DataNorm, spin1, sizeof( s1DataNorm ) );
  memcpy( s2DataNorm, spin2, sizeof( s2DataNorm ) );

  for( i = 0; i < 3; i++ )
  {
    s1Data[i] *= m1*m1;
    s2Data[i] *= m2*m2;
  }

  for ( i = 0; i < 3; i++ )
  {
    s1DataNorm[i] = s1Data[i]/mTotal/mTotal;
    s2DataNorm[i] = s2Data[i]/mTotal/mTotal;
  }
 
  if( debugPK )
  {
    printf("Calculating SigmaStar\n");
    fflush(NULL);
  }
  /* Populate the initial structures */
  if ( XLALSimIMRSpinEOBCalculateSigmaStar( sigmaStar, m1, m2, 
                              &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  if( debugPK )
  {
    printf("Calculating SigmaKerr\n");
    fflush(NULL);
  }
  if ( XLALSimIMRSpinEOBCalculateSigmaKerr( sigmaKerr, m1, m2, 
                              &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  memcpy( mSpin1, spin1, sizeof( mSpin1 ) );
  memcpy( mSpin2, spin2, sizeof( mSpin2 ) );

  for ( i = 0; i < 3; i++ )
  {
    mSpin1[i] *= m1*m1;
    mSpin2[i] *= m2*m2;
  }

  /* Calculate the value of a */
  seobParams.a = a = sqrt( inner_product(sigmaKerr->data,sigmaKerr->data) );

  seobParams.s1Vec = &s1VecOverMtMt;
  seobParams.s2Vec = &s2VecOverMtMt;
  
  /* Cartesian vectors needed to calculate Hamiltonian */
  cartPosVec.length = cartMomVec.length = 3;
  cartPosVec.data = cartPosData;
  cartMomVec.data = cartMomData;
  memset( cartPosData, 0, sizeof( cartPosData ) );
  memset( cartMomData, 0, sizeof( cartMomData ) );

  /* TODO: Insert potentially necessary checks on the arguments */

  seobParams.alignedSpins = 0;
  seobParams.tortoise     = 1;
  seobParams.sigmaStar    = sigmaStar;
  seobParams.sigmaKerr    = sigmaKerr;
  seobParams.seobCoeffs   = &seobCoeffs;
  seobParams.eobParams    = &eobParams;
  seobParams.nqcCoeffs    = &nqcCoeffs;
  eobParams.hCoeffs       = &hCoeffs;
  eobParams.prefixes      = &prefixes;
  seobCoeffs.SpinAlignedEOBversion = SpinAlignedEOBversion;
  eobParams.m1  = m1;
  eobParams.m2  = m2;
  eobParams.eta = eta;

  /* ************************************************* */
  /* Populate the initial structures                   */
  /* ************************************************* */
  /* Pre-compute the Hamiltonian coefficients */
  if( debugPK )
  {
    printf("Populating the Hamiltonian Coefficients. a = %.12e\n",a);
    fflush(NULL);
  }
  if ( XLALSimIMRCalculateSpinEOBHCoeffs( &seobCoeffs, eta, a, 
                          SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Pre-compute the coefficients for the Newtonian factor of hLM */
  if( debugPK )
  {
    printf("Populating the Newtonian Multipole Coefficients\n");
    fflush(NULL);
  }
  if ( XLALSimIMREOBComputeNewtonMultipolePrefixes( &prefixes, eobParams.m1,
			eobParams.m2 ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* ************************************************************ */
  /* ***** Calculate the derivatives loading dynamics from a file */
  /* ************************************************************ */
if(importDynamicsAndGetDerivatives)
{  
  nqcCoeffs.a1 = nqcCoeffs.a2 = nqcCoeffs.a3 = nqcCoeffs.a3S = nqcCoeffs.a4 = nqcCoeffs.a5 = nqcCoeffs.b1 = nqcCoeffs.b2 = nqcCoeffs.b3 = nqcCoeffs.b4 = 0;
  if(debugPK)printf("\tl = %d, m = %d, NQC: a1 = %.16e, a2 = %.16e, a3 = %.16e, a3S = %.16e, a4 = %.16e, a5 = %.16e\n\tb1 = %.16e, b2 = %.16e, b3 = %.16e, b4 = %.16e\n", 
                           2, 2, nqcCoeffs.a1, nqcCoeffs.a2, nqcCoeffs.a3, nqcCoeffs.a3S, nqcCoeffs.a4, nqcCoeffs.a5, 
                           nqcCoeffs.b1, nqcCoeffs.b2, nqcCoeffs.b3, nqcCoeffs.b4 );
  //FILE* dynfin = fopen( "Data_q_5_chi1_0.247214_chi2_0.352671_dSO_-69.5_dSS_2.75P.dat", "r" );
  FILE* dynfin = fopen( "seobDynamics.dat", "r" );
  FILE* dynfout= fopen( "DynamicsDerivatives.dat", "w");
  
  double tin, xin, yin, zin, pxin, pyin, pzin, s1xin, s1yin, s1zin, s2xin, s2yin, s2zin;
  double UNUSED tmp1, tmp2;
  REAL8 valuesin[14], dvaluesout[14]; 
  double H, flux;
   
  REAL8Vector polarDynamics, cartDynamics;
  REAL8 polData[4];
  cartDynamics.data = valuesin;
  polarDynamics.data= polData;
  
  int tmpj = 0;
  if(debugPK){
	  printf("Reading dynamics file to compute derivatives\n"); }
  while(fscanf(dynfin, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le", 
					&tin, &xin, &yin, &zin, &pxin, &pyin, &pzin, 
					&s1xin, &s1yin, &s1zin, &s2xin, &s2yin, &s2zin 
					,&tmp1, &tmp2
					 ) == 15)
  {
	if(debugPK)printf("%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\n", 
					tin, xin, yin, zin, pxin, pyin, pzin, 
					s1xin, s1yin, s1zin, s2xin, s2yin, s2zin );
					
    valuesin[0] = xin;
    valuesin[1] = yin;
    valuesin[2] = zin;
    valuesin[3] = pxin;
    valuesin[4] = pyin;
    valuesin[5] = pzin;
    valuesin[6] = s1xin;
    valuesin[7] = s1yin;
    valuesin[8] = s1zin;
    valuesin[9] = s2xin;
    valuesin[10]= s2yin;
    valuesin[11]= s2zin;

    cartPosVec.data = cartPosData;
    cartMomVec.data = cartMomData;
    memcpy( cartPosData, valuesin, 3*sizeof( REAL8 ) );
    memcpy( cartMomData, valuesin+3, 3*sizeof( REAL8 ) );
     
    /* Populate the SPIN structures 
     * These set s{1,2}Vec, s{1,2}VecOverMtMt, which in turn set
     * the spins stored in seobParams */
    for( i=0; i<3; i++ )
    {
	  s1DataNorm[i] = valuesin[i+6];
	  s2DataNorm[i] = valuesin[i+9];
	  s1Data[i] = s1DataNorm[i] * mTotal * mTotal;
	  s2Data[i] = s2DataNorm[i] * mTotal * mTotal;
    }
  
    s1Vec.data   = s1Data;
    s2Vec.data   = s2Data;
   
    s1VecOverMtMt.data   = s1DataNorm;
    s2VecOverMtMt.data   = s2DataNorm;

    seobParams.s1Vec = &s1VecOverMtMt;
    seobParams.s2Vec = &s2VecOverMtMt;

    if ( XLALSimIMRSpinEOBCalculateSigmaStar( sigmaStar, m1, m2, 
                              &s1Vec, &s2Vec ) == XLAL_FAILURE )
		printf("\nSomething went wrong evaluating SigmaStar in line %d of the input file\n", 
			tmpj );
	
	if ( XLALSimIMRSpinEOBCalculateSigmaKerr( sigmaKerr, m1, m2, 
                              &s1Vec, &s2Vec ) == XLAL_FAILURE )
      printf("\nSomething went wrong evaluating SigmaKerr in line %d of the input file\n", 
			tmpj );
	
      /* Calculate the value of a */
    seobParams.a = a = sqrt( sigmaKerr->data[0]*sigmaKerr->data[0] 
		+ sigmaKerr->data[1]*sigmaKerr->data[1] 
		+ sigmaKerr->data[2]*sigmaKerr->data[2] );
    
	if ( debugPK )
    {
	  printf("\nValue of a calculated = %.16le\n", a);
	  
	  /* Print out all spin parameters */
	  printf("sigmaStar = {%.12e,%.12e,%.12e}\n, sigmaKerr = {%.12e,%.12e,%.12e}\n",
			(double) sigmaStar->data[0], (double) sigmaStar->data[1], 
			(double) sigmaStar->data[2], (double) sigmaKerr->data[0],
			(double) sigmaKerr->data[1], (double) sigmaKerr->data[2]);
	  printf("a = %.12e, tplspin = %.12e, chiS = %.12e, chiA = %.12e\n", 
			(double) a, (double) tplspin, (double) chiS, (double) chiA);
	  printf("a is used to compute Hamiltonian coefficients,\n tplspin and chiS and chiA for the multipole coefficients\n");

	  printf("s1Vec = {%.12e,%.12e,%.12e}\n", (double) s1VecOverMtMt.data[0], 
	  (double) s1VecOverMtMt.data[1], (double) s1VecOverMtMt.data[2]);
	  printf("s2Vec = {%.12e,%.12e,%.12e}\n", (double) s2VecOverMtMt.data[0],
	  (double) s2VecOverMtMt.data[1], (double) s2VecOverMtMt.data[2]);

      fflush(NULL);
    } 


    memset( dvaluesout, 0, 14*sizeof(REAL8) );
    status = XLALSpinHcapRvecDerivative( 0, valuesin, dvaluesout, (void*) &seobParams);
    memcpy( rdotvec, dvaluesout, 3*sizeof(REAL8) );
    
    /* Calculate r cross rDot */
    cross_product( cartPosData, rdotvec, rcrossrdot );
    rcrossrdotNorm = sqrt(inner_product( rcrossrdot, rcrossrdot ));
    for( i = 0; i < 3; i++ )
      rcrossrdot[i] /= rcrossrdotNorm;
   
    // Since s?Data have non-normalized spins, we need to divide by m?^2
    s1dotLN = inner_product( s1Data, rcrossrdot ) / (m1*m1); 
    s2dotLN = inner_product( s2Data, rcrossrdot ) / (m1*m1);  
    chiS = 0.5 * (s1dotLN + s2dotLN);
    chiA = 0.5 * (s1dotLN - s2dotLN);

    switch ( SpinAlignedEOBversion )
    {
     case 1:
       tplspin = 0.0;
       break;
     case 2:
       tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
    }
    if(debugPK)printf("\nCalculating Hamiltonian coefficients\n");
    /* Populate the different structures */
    if ( XLALSimIMRCalculateSpinEOBHCoeffs( &seobCoeffs, eta, a, 
                          SpinAlignedEOBversion ) == XLAL_FAILURE )
      printf("\nSomething went wrong evaluating XLALSimIMRCalculateSpinEOBHCoeffs in line %d of the input file\n", 
			tmpj );

    if(debugPK)printf("\nCalculating Waveform coefficients\n");
    if ( XLALSimIMREOBCalcSpinFacWaveformCoefficients( &hCoeffs, m1, m2, eta, 
        tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE )
      printf("\nSomething went wrong evaluating XLALSimIMRCalculateSpinEOBHCoeffs in line %d of the input file\n", 
			tmpj );

    /* CALCULATE THE DERIVATIVES */
    if(debugPK)printf("\nCalculating Derivatives\n");
    memset( dvaluesout, 0, 14*sizeof(REAL8) );
    if( XLALSpinHcapNumericalDerivative( (REAL8) tin, valuesin, dvaluesout, 
								(void*) &seobParams ) != XLAL_SUCCESS )
      printf("\nSomething went wrong evaluating derivatives in line %d of the input file\n", 
			tmpj );
    
    /* CALCULATE THE HAMILTONIAN AND FLUX */
    if(debugPK)printf("\nCalculating Hamiltonian\n");
    H = XLALSimIMRSpinEOBHamiltonian( eta, &cartPosVec, &cartMomVec, 
			&s1VecOverMtMt, &s2VecOverMtMt, sigmaKerr, sigmaStar, 
			seobParams.tortoise, seobParams.seobCoeffs ); 
    H *= mTotal;

    /* Calculate POLAR dynamics for Flux */
    polData[0] = sqrt( inner_product(cartPosData, cartPosData) );
    polData[1] = 0;
    polData[2] = inner_product(cartPosData, cartMomData) / polData[0];
    cross_product( cartPosData, cartMomData, rcrossp );
    polData[3] = sqrt(inner_product(rcrossp, rcrossp));
    omega = rcrossrdotNorm / (polData[0] * polData[0]);
    
    if(debugPK)printf("\nCalculating Flux\n");
    flux  = XLALInspiralPrecSpinFactorizedFlux( &polarDynamics, &cartDynamics, 
              &nqcCoeffs, omega, &seobParams, H/mTotal, 8, SpinAlignedEOBversion );
    flux = flux / eta;
    
    /* OUTPUT THE DERIVATIVES */
    if(debugPK)printf("\n\nDynamics imported\n");
    fprintf(dynfout, "%.16e\t", tin);
    for( i=0; i<12; i++ )
	{
		fprintf(dynfout, "%.16e\t", valuesin[i]);
		if(debugPK)printf("%e ", valuesin[i]);
	}
    
    if(debugPK)printf("\n\nDerivatives calculated:\n");
    for( i=0; i<14; i++ )
    {
	  fprintf(dynfout, "%.16e\t", dvaluesout[i]);
	  if(debugPK)printf("%e ", dvaluesout[i]);
    }
    	
	fprintf(dynfout, "%.16e\t%.16e\n", H, flux);
	if(debugPK)printf("%e\t%e\n", H, flux);
	
	tmpj += 1;
  }
  
  fclose(dynfin);
  fclose(dynfout);
  
  if(debugPK){
	  printf("Derivatives written for the dynamics file!\n"); }  
}
  /* ************************************************* */
  /* ***** Set up the INITIAL CONDITIONS               */
  /* ************************************************* */
  /*
   * STEP 1) Solve for initial conditions
   */

  REAL8 temp32;
  temp32 = fMin * inc; temp32 *= 2;

  /*
  if( debugPK )
  {
    printf("Calling the XLALSimIMRSpinEOBInitialConditions function!\n");
    fflush(NULL);
  }
  
  if ( XLALSimIMRSpinEOBInitialConditions( tmpValues2, m1, m2, fMin, inc,
                	mSpin1, mSpin2, &seobParams ) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }*/
  
  //exit(0);
  //YP::{x,y,z,px,py,pz,s1x,s1y,s1z,s2x,s2y,s2z} =
  //{15.87, 0, 0, -0.000521675194648, 0.278174373488, -0.00012666165246,
  // -0.270452950188, -0.216802131414, 0.00133043857763, 0, 0, 0};
#if 0
  values->data[0] = 0.;
  values->data[1] = 12.845228155660482;
  values->data[2] = -4.553894189373296;
  values->data[3] = -0.006165987975074341 /*/ eta*/;
  values->data[4] = 0.10049046440176972 /*/ eta*/;
  values->data[5] = 0.28877341851636174 /*/ eta*/;
  values->data[6] = 4.*m1*m1 * 0.1125;
  values->data[7] = 4.*m1*m1 * -0.09742785792574934;
  values->data[8] = 4.*m1*m1 * -0.16875;
  values->data[9] = 4.*m2*m2 *-0.1060660171779821;
  values->data[10] =4.*m2*m2 * 6.938893903907228e-18;
  values->data[11] =4.*m2*m2 * -0.10606601717798211;
#endif
#if 0
  values->data[0] = 15.87;
  values->data[1] = 0.;
  values->data[2] = 0.;
  values->data[3] = -0.000521675194648;
  values->data[4] = 0.278174373488;
  values->data[5] = -0.00012666165246;
  values->data[6] = -0.270452950188;
  values->data[7] = -0.216802131414;
  values->data[8] = 0.00133043857763;
  values->data[9] = 0.;
  values->data[10] = 0.;
  values->data[11] = 0.;
#endif

  /* Assume that initial conditions are available at this point, to 
   * compute the chiS and chiA parameters. 
   * Calculate the values of chiS and chiA, as given in Eq.16 of 
   * Precessing EOB paper. Assuming \vec{L} to be pointing in the 
   * direction of \vec{r}\times\vec{p} */
  if(debugPK)printf("\nReached the point where LN is to be calculated\n");
  /* Calculate rDot = \f$\partial Hreal / \partial p_r\f$ */
  memset( dvalues->data, 0, 14 * sizeof(REAL8) );
  status = XLALSpinHcapRvecDerivative( 0, values->data, dvalues->data, (void*) &seobParams);
  if(debugPK)printf("\nCalculated Rdot\n");
  memcpy( rdotvec, dvalues->data, 3*sizeof(REAL8) );

  /* Calculate r cross rDot */ 
  memcpy( rvec, values->data, 3*sizeof(REAL8) );
  cross_product( rvec, rdotvec, rcrossrdot );
  rcrossrdotNorm = sqrt(inner_product( rcrossrdot, rcrossrdot ));
  for( i = 0; i < 3; i++ )
    rcrossrdot[i] /= rcrossrdotNorm;
   
  s1dotLN = inner_product( spin1, rcrossrdot );
  s2dotLN = inner_product( spin2, rcrossrdot );
    
  rcrossp[0] = values->data[1]*values->data[5] - values->data[2]*values->data[4];
  rcrossp[1] = values->data[2]*values->data[3] - values->data[0]*values->data[5];
  rcrossp[2] = values->data[0]*values->data[4] - values->data[1]*values->data[3];
  rcrosspMag = sqrt(rcrossp[0]*rcrossp[0] + rcrossp[1]*rcrossp[1] + 
        rcrossp[2]*rcrossp[2]);
  
  rcrossp[0] /= rcrosspMag;
  rcrossp[1] /= rcrosspMag;
  rcrossp[2] /= rcrosspMag;

  s1dotL = spin1[0]*rcrossp[0] + spin1[1]*rcrossp[1] + spin1[2]*rcrossp[2];
  s2dotL = spin2[0]*rcrossp[0] + spin2[1]*rcrossp[1] + spin2[2]*rcrossp[2];

  chiS = 0.5 * (s1dotLN + s2dotLN);
  chiA = 0.5 * (s1dotLN - s2dotLN);

  /* Compute the test-particle limit spin of the deformed-Kerr background */
  /* TODO: Check this is actually the way it works in LAL */
  switch ( SpinAlignedEOBversion )
  {
     case 1:
       tplspin = 0.0;
       break;
     case 2:
       tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }

  /* ************************************************* */
  /* Populate the Waveform initial structures          */
  /* ************************************************* */
#if 0
  /* Pre-compute the non-spinning coefficients for hLM */
  if ( XLALSimIMREOBCalcPrecNoSpinFacWaveformCoefficients( &hCoeffs, m1, m2, eta, 
        tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Pre-compute the spinning coefficients for hLM */
  if ( XLALSimIMREOBCalcPrecSpinFacWaveformCoefficients( &hCoeffs, m1, m2, eta, 
        tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }
#endif
  
  /* Pre-compute the non-spinning and spinning coefficients for hLM factors */
  if( debugPK )
  {
    printf("Calling XLALSimIMREOBCalcSpinFacWaveformCoefficients for hlm Coefficients!\n");
    printf("tplspin = %.12e, chiS = %.12e, chiA = %.12e\n", tplspin, chiS, chiA);
    fflush(NULL);
  }
  if ( XLALSimIMREOBCalcSpinFacWaveformCoefficients( &hCoeffs, m1, m2, eta, 
        tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* ************************************************* */
  /* Compute the coefficients for the NQC Corrections  */
  /* ************************************************* */
  spinNQC = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
  switch ( SpinAlignedEOBversion )
  {
	  case 1:         
	    if(debugPK)printf("\t NQC: spin used = %.12e\n", spinNQC);
	    XLALSimIMRGetEOBCalibratedSpinNQC( &nqcCoeffs, 2, 2, eta, spinNQC );
	    break;
	  case 2:
	    if(debugPK)printf("\t NQC: spins used = %.12e, %.12e\n", spinNQC, chiA);
	    // XLALSimIMRGetEOBCalibratedSpinNQCv2( &nqcCoeffs, 2, 2, eta, ak->a );
	    XLALSimIMRGetEOBCalibratedSpinNQC3D( &nqcCoeffs, 2, 2, m1, m2, spinNQC, chiA );
	    break;
	  default:
	    XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
	    XLAL_ERROR( XLAL_EINVAL );
	    break;
  }
  nqcCoeffs.a1 = nqcCoeffs.a2 = nqcCoeffs.a3 = nqcCoeffs.a3S = nqcCoeffs.a4 = 
  nqcCoeffs.a5 = nqcCoeffs.b1 = nqcCoeffs.b2 = nqcCoeffs.b3 = nqcCoeffs.b4 = 0;
  if(debugPK)printf("\tl = %d, m = %d, NQC: a1 = %.16e, a2 = %.16e, a3 = %.16e, a3S = %.16e, a4 = %.16e, a5 = %.16e\n\tb1 = %.16e, b2 = %.16e, b3 = %.16e, b4 = %.16e\n", 
                       2, 2, nqcCoeffs.a1, nqcCoeffs.a2, nqcCoeffs.a3, 
                       nqcCoeffs.a3S, nqcCoeffs.a4, nqcCoeffs.a5, 
                       nqcCoeffs.b1, nqcCoeffs.b2, nqcCoeffs.b3, nqcCoeffs.b4 );
 

  if ( debugPK )
  {
	  /* Print out all mass parameters */
	  printf("m1SI = %.12e, m2SI = %.12e, m1 = %.12e, m2 = %.12e\n",
			(double) m1SI, (double) m2SI, (double) m1, (double) m2 );
	  printf("mTotal = %.12e, mTScaled = %.12e, eta = %.12e\n", 
			(double) mTotal, (double) mTScaled, (double) eta );
	  /* Print out all spin parameters */
	  printf("spin1 = {%.12e,%.12e,%.12e}, spin2 = {%.12e,%.12e,%.12e}\n",
			(double) spin1[0], (double) spin1[1], (double) spin1[2],
			(double) spin2[0], (double) spin2[1], (double) spin2[2]);
	  printf("mSpin1 = {%.12e,%.12e,%.12e}, mSpin2 = {%.12e,%.12e,%.12e}\n",
			(double) mSpin1[0], (double) mSpin1[1], (double) mSpin1[2],
			(double) mSpin2[0], (double) mSpin2[1], (double) mSpin2[2]);
	
	  printf("sigmaStar = {%.12e,%.12e,%.12e},\n sigmaKerr = {%.12e,%.12e,%.12e}\n",
			(double) sigmaStar->data[0], (double) sigmaStar->data[1], 
			(double) sigmaStar->data[2], (double) sigmaKerr->data[0],
			(double) sigmaKerr->data[1], (double) sigmaKerr->data[2]);
	  printf("a = %.12e, tplspin = %.12e, chiS = %.12e, chiA = %.12e\n", 
			(double) a, (double) tplspin, (double) chiS, (double) chiA);
	  printf("a is used to compute Hamiltonian coefficients,\n tplspin and chiS and chiA for the multipole coefficients\n");

	  printf("s1Vec = {%.12e,%.12e,%.12e}\n", (double) s1VecOverMtMt.data[0], 
	  (double) s1VecOverMtMt.data[1], (double) s1VecOverMtMt.data[2]);
	  printf("s2Vec = {%.12e,%.12e,%.12e}\n", (double) s2VecOverMtMt.data[0],
	  (double) s2VecOverMtMt.data[1], (double) s2VecOverMtMt.data[2]);

      fflush(NULL);
  } 
  
  FILE *out = NULL;

  /* Initialize the GSL integrator */
  if (!(integrator = XLALAdaptiveRungeKutta4Init(14, XLALSpinHcapNumericalDerivative,
				XLALEOBSpinStopCondition, EPS_ABS, EPS_REL)))
  {
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  integrator->stopontestonly = 1;
  for( i = 0; i < 3; i++ )
  {
	if( debugPK )
		printf("\n r = {%f,%f,%f}\n", 
		values->data[0], values->data[1], values->data[2]);
	fflush(NULL);
}
  
  if(debugPK)printf("\n\n BEGINNING THE EVOLUTION\n\n");
  retLen = XLALAdaptiveRungeKutta4( integrator, &seobParams, values->data, 
				0., 20./mTScaled, deltaT/mTScaled, &dynamics );
   if(debugPK)printf("\n\n FINISHED THE EVOLUTION\n\n");
 if ( retLen == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /*FILE *indat = fopen( "inputDynData.dat", "r" );
  i = 0;
  while ( !feof(indat) ) {
    fscanf(indat,"%e",dynamics->data+i);
    i++
  }
  retLen = (i+1)/14;*/ 

  printf("To be the man, you've got to beat the man! Woooooooo!!!!\n" );

  REAL8Vector tVec;
  tVec.length = retLen;
  tVec.data   = dynamics->data;

  REAL8 *posVecx = dynamics->data+retLen;
  REAL8 *posVecy = dynamics->data+2*retLen;
  REAL8 *posVecz = dynamics->data+3*retLen;
  REAL8 *momVecx = dynamics->data+4*retLen;
  REAL8 *momVecy = dynamics->data+5*retLen;
  REAL8 *momVecz = dynamics->data+6*retLen;
  REAL8 *s1Vecx = dynamics->data+7*retLen;
  REAL8 *s1Vecy = dynamics->data+8*retLen;
  REAL8 *s1Vecz = dynamics->data+9*retLen;
  REAL8 *s2Vecx = dynamics->data+10*retLen;
  REAL8 *s2Vecy = dynamics->data+11*retLen;
  REAL8 *s2Vecz = dynamics->data+12*retLen;
  //YP: change variable names to match those used in the waveform section
  //REAL8 *vphi   = dynamics->data+13*retLen;
  //REAL8 *phpart2= dynamics->data+14*retLen;
  REAL8 *phiDMod= dynamics->data+13*retLen;
  REAL8 *phiMod = dynamics->data+14*retLen;

  out = fopen( "seobDynamics.dat", "w" );
  for ( i = 0; i < retLen; i++ )
  {
    //YP: output orbital phase and phase modulation separately, instead of their sum
    fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", 
    i*deltaT/mTScaled, 
    posVecx[i], posVecy[i], posVecz[i], momVecx[i], momVecy[i], momVecz[i],
    s1Vecx[i], s1Vecy[i], s1Vecz[i], s2Vecx[i], s2Vecy[i], s2Vecz[i],
    phiDMod[i], phiMod[i] );
  }
  fclose( out );

  printf("DYNAMICS WRITTEN!!!!\n" );

  /* ==============================
   *   Waveform Generation
   * ==============================
   */

  REAL8 tPeakOmega, tAttach, combSize, /*longCombSize,*/ deltaNQC;
  REAL8 vX, vY, vZ, rCrossV_x, rCrossV_y, rCrossV_z, vOmega, omegasav, omegasav2;
  REAL8 magR, Lx, Ly, Lz, magL, LNhx, LNhy, LNhz, magLN, Jx, Jy, Jz, magJ;
  REAL8 aI2P, bI2P, gI2P, aP2J, bP2J, gP2J;
  REAL8 chi1J, chi2J, chiJ, kappaJL;
  REAL8 JframeEx[3], JframeEy[3], JframeEz[3];
  REAL8 LframeEx[3], LframeEy[3], LframeEz[3];
  combSize = 7.5;
  /* WaveStep 1
   * Locate merger point (max omega), calculate J, chi and kappa at merger, and construct final J frame
   */
  /* WaveStep 1.1: locate merger point */
  omegasav2 = -1.0;
  omegasav  = -0.5;
  omega     =  0.0;
  for ( i = 0; i < retLen; i++ )
  {
    for ( j = 0; j < values->length; j++ )
    {
      values->data[j] = *(dynamics->data+(j+1)*retLen+i);
    }
    vX = XLALSpinHcapNumDerivWRTParam( 3, values->data, &seobParams );
    vY = XLALSpinHcapNumDerivWRTParam( 4, values->data, &seobParams );
    vZ = XLALSpinHcapNumDerivWRTParam( 5, values->data, &seobParams );
    rCrossV_x = posVecy[i] * vZ - posVecz[i] * vY;
    rCrossV_y = posVecz[i] * vX - posVecx[i] * vZ;
    rCrossV_z = posVecx[i] * vY - posVecy[i] * vX;

    magR = sqrt(posVecx[i]*posVecx[i] + posVecy[i]*posVecy[i] + posVecz[i]*posVecz[i] );
    omega = sqrt(rCrossV_x*rCrossV_x + rCrossV_y*rCrossV_y + rCrossV_z*rCrossV_z ) / (magR*magR);
    if ( omega < omegasav )
      break;
    else
    {
      omegasav2 = omegasav;
      omegasav  = omega;
    }
  }
  if ( i == retLen - 1 )
  {
    printf("YP: Error! Failed to find peak of omega!\n");
    abort();
  }
  else
  {
    tPeakOmega = (i-(4.*omegasav-3.*omega-omegasav2)/(2.*omegasav-omega-omegasav2)/2.)*deltaT/mTScaled;
  }
  /* WaveStep 1.2: calculate J at merger */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;
  spline = gsl_spline_alloc( gsl_interp_cspline, retLen );
  acc    = gsl_interp_accel_alloc();

  for ( j = 0; j < values->length; j++ )
  {
    gsl_spline_init( spline, dynamics->data, dynamics->data+(j+1)*retLen, retLen );
    values->data[j] = gsl_spline_eval( spline, tPeakOmega, acc );
  }
  vX = XLALSpinHcapNumDerivWRTParam( 3, values->data, &seobParams );
  vY = XLALSpinHcapNumDerivWRTParam( 4, values->data, &seobParams );
  vZ = XLALSpinHcapNumDerivWRTParam( 5, values->data, &seobParams );
  rCrossV_x = values->data[1] * vZ - values->data[2] * vY;
  rCrossV_y = values->data[2] * vX - values->data[0] * vZ;
  rCrossV_z = values->data[0] * vY - values->data[1] * vX;
  Lx = values->data[1] * values->data[5] - values->data[2] * values->data[4];
  Ly = values->data[2] * values->data[3] - values->data[0] * values->data[5];
  Lz = values->data[0] * values->data[4] - values->data[1] * values->data[3];
  magL = sqrt( Lx*Lx + Ly*Ly + Lz*Lz );
  Jx = eta*Lx + values->data[6] + values->data[9];
  Jy = eta*Ly + values->data[7] + values->data[10];
  Jz = eta*Lz + values->data[8] + values->data[11];
  magJ = sqrt( Jx*Jx + Jy*Jy + Jz*Jz );
  /* WaveStep 1.3: calculate chi and kappa at merger */
  chi1J = values->data[6]*Jx + values->data[7] *Jy + values->data[8] *Jz;
  chi2J = values->data[9]*Jx + values->data[10]*Jy + values->data[11]*Jz;
  chi1J/= magJ*m1*m1/mTotal/mTotal;
  chi2J/= magJ*m2*m2/mTotal/mTotal;
  chiJ = (chi1J+chi2J)/2. + (chi1J-chi2J)/2.*sqrt(1. - 4.*eta)/(1. - 2.*eta);
  if (chiJ <= 0.) {
    deltaNQC = 2.5;
  }
  else {
    deltaNQC = 2.5 + 1.77*(chiJ/0.43655)*(chiJ/0.43655)*(chiJ/0.43655)*(chiJ/0.43655);
  }
  kappaJL      = (Lx*Jx + Ly*Jy + Lz*Jz) / magL / magJ;
  combSize    *= 1.0 + 9.0 * (1.0 - fabs(kappaJL));
  //longCombSize = combSize;
  deltaNQC    += 10.0 * (1.0 - fabs(kappaJL));
  tAttach = tPeakOmega - deltaNQC;
  /* WaveStep 1.4: calculate combsize and deltaNQC */
  switch ( SpinAlignedEOBversion )
  {
     case 1:
       combSize = 7.5;
       if ( chiJ <= 0. ) {
         deltaNQC = 2.5;
       }
       else {
         deltaNQC = 2.5 + 1.77*(chiJ/0.43655)*(chiJ/0.43655)*(chiJ/0.43655)*(chiJ/0.43655);
       }
       break;
     case 2:
       combSize = 12.;
       if ( chiJ > 0.8 ) combSize = 13.5;
       if ( chi1J == 0. && chi2J == 0. ) combSize = 11.;
       if ( chiJ <= 0. ) {
         deltaNQC = 2.5 + (1. + chiJ)*(-2.5 + 2.5*sqrt(1.-4.*eta));
       }
       else if ( chiJ <= 0.8 ) {
         deltaNQC = (0.75*eta*chiJ + sqrt(1. - 4.*eta)) * (2.5 + 10.*chiJ*chiJ + 24.*chiJ*chiJ*chiJ*chiJ);
       }
       else {
         deltaNQC = (0.75*eta*chiJ + sqrt(1. - 4.*eta)) * (57.1755 - 48.0564*chiJ);
       }
       break;
     default:
       XLALPrintError( "XLAL Error - %s: wrong SpinAlignedEOBversion value, must be 1 or 2!\n", __func__ );
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }

  /* WaveStep 1.5: get  */
  modefreqVec.length = 1;
  modefreqVec.data   = &modeFreq;

  if ( XLALSimIMREOBGenerateQNMFreqV2( &modefreqVec, m1, m2, spin1, spin2, 2, 2, 1, spinEOBApproximant ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }
  retLenRDPatch = (UINT4)ceil( 20 / ( cimag(modeFreq) * deltaT ));
  /* WaveStep 1.6: construct J-frame */
  JframeEz[0] = Jx / magJ;
  JframeEz[1] = Jy / magJ;
  JframeEz[2] = Jz / magJ;
  if ( 1.-JframeEz[2] < 1.0e-16 ) {
    JframeEx[0] = 1.;
    JframeEx[1] = 0.;
    JframeEx[2] = 0.;
  }
  else {
    JframeEx[0] = JframeEz[1];
    JframeEx[1] = -JframeEz[0];
    JframeEx[2] = 0.;
  }
  JframeEx[0] /= sqrt( JframeEz[0]*JframeEz[0] + JframeEz[1]*JframeEz[1] );
  JframeEx[1] /= sqrt( JframeEz[0]*JframeEz[0] + JframeEz[1]*JframeEz[1] );
  JframeEy[0] = JframeEz[1]*JframeEx[2] - JframeEz[2]*JframeEx[1];
  JframeEy[1] = JframeEz[2]*JframeEx[0] - JframeEz[0]*JframeEx[2];
  JframeEy[2] = JframeEz[0]*JframeEx[1] - JframeEz[1]*JframeEx[0];

  /* WaveStep 2
   * Calculate quasi-nonprecessing waveforms
   */
  /* WaveStep 2.1: create time-series containers for euler angles and hlm harmonics */
  SphHarmTimeSeries *hlmPTS = NULL;
  SphHarmTimeSeries *hIMRlmJTS = NULL;
  REAL8TimeSeries UNUSED *hPlusTS  = XLALCreateREAL8TimeSeries( "H_PLUS", &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries UNUSED *hCrossTS = XLALCreateREAL8TimeSeries( "H_CROSS", &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries *alphaI2PTS = XLALCreateREAL8TimeSeries( "alphaI2P", &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries  *betaI2PTS = XLALCreateREAL8TimeSeries(  "betaI2P", &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries *gammaI2PTS = XLALCreateREAL8TimeSeries( "gammaI2P", &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries *alphaP2JTS = XLALCreateREAL8TimeSeries( "alphaP2J", &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries  *betaP2JTS = XLALCreateREAL8TimeSeries(  "betaP2J", &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries *gammaP2JTS = XLALCreateREAL8TimeSeries( "gammaP2J", &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8Sequence *tlist        = NULL;
  REAL8Sequence *tlistRDPatch = NULL;
  COMPLEX16TimeSeries *h22TS   = XLALCreateCOMPLEX16TimeSeries( "H_22",   &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h21TS   = XLALCreateCOMPLEX16TimeSeries( "H_21",   &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h20TS   = XLALCreateCOMPLEX16TimeSeries( "H_20",   &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h2m1TS  = XLALCreateCOMPLEX16TimeSeries( "H_2m1",  &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h2m2TS  = XLALCreateCOMPLEX16TimeSeries( "H_2m2",  &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h22PTS  = XLALCreateCOMPLEX16TimeSeries( "HP_22",  &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h21PTS  = XLALCreateCOMPLEX16TimeSeries( "HP_21",  &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h20PTS  = XLALCreateCOMPLEX16TimeSeries( "HP_20",  &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h2m1PTS = XLALCreateCOMPLEX16TimeSeries( "HP_2m1", &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h2m2PTS = XLALCreateCOMPLEX16TimeSeries( "HP_2m2", &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h22JTS  = XLALCreateCOMPLEX16TimeSeries( "HJ_22",  &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h21JTS  = XLALCreateCOMPLEX16TimeSeries( "HJ_21",  &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h20JTS  = XLALCreateCOMPLEX16TimeSeries( "HJ_20",  &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h2m1JTS = XLALCreateCOMPLEX16TimeSeries( "HJ_2m1", &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h2m2JTS = XLALCreateCOMPLEX16TimeSeries( "HJ_2m2", &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *hJTS    = XLALCreateCOMPLEX16TimeSeries( "HJ",     &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *hIMR22JTS  = XLALCreateCOMPLEX16TimeSeries( "HIMRJ_22",  &tc, 0.0, deltaT, &lalStrainUnit,
                                                                   retLen + retLenRDPatch );
  COMPLEX16TimeSeries *hIMR21JTS  = XLALCreateCOMPLEX16TimeSeries( "HIMRJ_21",  &tc, 0.0, deltaT, &lalStrainUnit,
                                                                   retLen + retLenRDPatch );
  COMPLEX16TimeSeries *hIMR20JTS  = XLALCreateCOMPLEX16TimeSeries( "HIMRJ_20",  &tc, 0.0, deltaT, &lalStrainUnit,
                                                                   retLen + retLenRDPatch );
  COMPLEX16TimeSeries *hIMR2m1JTS = XLALCreateCOMPLEX16TimeSeries( "HIMRJ_2m1", &tc, 0.0, deltaT, &lalStrainUnit,
                                                                   retLen + retLenRDPatch );
  COMPLEX16TimeSeries *hIMR2m2JTS = XLALCreateCOMPLEX16TimeSeries( "HIMRJ_2m2", &tc, 0.0, deltaT, &lalStrainUnit,
                                                                   retLen + retLenRDPatch );
  COMPLEX16TimeSeries *hIMRJTS    = XLALCreateCOMPLEX16TimeSeries( "HIMRJ",     &tc, 0.0, deltaT, &lalStrainUnit,
                                                                   retLen + retLenRDPatch );
  if ( !(tlist = XLALCreateREAL8Vector( retLen )) || !(tlistRDPatch = XLALCreateREAL8Vector( retLen + retLenRDPatch )) )
  {
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  memset( tlist->data, 0, tlist->length * sizeof( REAL8 ));
  memset( tlistRDPatch->data, 0, tlistRDPatch->length * sizeof( REAL8 ));
  for ( i = 0; i < retLen; i++ )
  {
    tlist->data[i] = i * deltaT/mTScaled;
    tlist->data[i] = i * deltaT/mTScaled;
  }
  for ( i = retLen; i < retLen + retLenRDPatch; i++ )
  {
    tlist->data[i] = i * deltaT/mTScaled;
  }

  /* WaveStep 2.2: get Calibrated NQC coeffcients, based on underlying nonprecessing model */
  switch ( SpinAlignedEOBversion )
  {
     case 1:
       if ( XLALSimIMRGetEOBCalibratedSpinNQC( &nqcCoeffs, 2, 2, eta, a ) == XLAL_FAILURE )
       {
         XLAL_ERROR( XLAL_EFUNC );
       }
       break;
     case 2:
       if ( XLALSimIMRGetEOBCalibratedSpinNQC3D( &nqcCoeffs, 2, 2, m1, m2, a, chiA ) == XLAL_FAILURE )
       {
         XLAL_ERROR( XLAL_EFUNC );
       }
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }

  /* WaveStep 2.3: main loop for quasi-nonprecessing waveform generation */
  out = fopen( "rotDynamics.dat", "w" );
  for ( i = 0; i < retLen; i++ )
  {
    for ( j = 0; j < values->length; j++ )
    {
      values->data[j] = dynamics->data[(j+1)*retLen + i];
      if ( i == 0 ) printf("value %d: %.16e\n",j,values->data[j]);
    }

    vX = XLALSpinHcapNumDerivWRTParam( 3, values->data, &seobParams );
    vY = XLALSpinHcapNumDerivWRTParam( 4, values->data, &seobParams );
    vZ = XLALSpinHcapNumDerivWRTParam( 5, values->data, &seobParams );

    /* Cartesian vectors needed to calculate Hamiltonian */
    cartPosVec.length = cartMomVec.length = 3;
    cartPosVec.data = cartPosData;
    cartMomVec.data = cartMomData;
    memset( cartPosData, 0, sizeof( cartPosData ) );
    memset( cartMomData, 0, sizeof( cartMomData ) );

    rCrossV_x = posVecy[i] * vZ - posVecz[i] * vY;
    rCrossV_y = posVecz[i] * vX - posVecx[i] * vZ;
    rCrossV_z = posVecx[i] * vY - posVecy[i] * vX;

    magR = sqrt(posVecx[i]*posVecx[i] + posVecy[i]*posVecy[i] + posVecz[i]*posVecz[i] );
    omega = sqrt(rCrossV_x*rCrossV_x + rCrossV_y*rCrossV_y + rCrossV_z*rCrossV_z ) / (magR*magR);
    vOmega = v = cbrt( omega );
    amp = amp0 * vOmega * vOmega;

    //LNhx = posVecy[i] * momVecz[i] - posVecz[i] * momVecy[i];
    //LNhy = posVecz[i] * momVecx[i] - posVecx[i] * momVecz[i];
    //LNhz = posVecx[i] * momVecy[i] - posVecy[i] * momVecx[i];
    LNhx  = rCrossV_x;
    LNhy  = rCrossV_y;
    LNhz  = rCrossV_z;
    magLN = sqrt(LNhx*LNhx + LNhy*LNhy + LNhz*LNhz);
    LNhx  = LNhx / magLN;
    LNhy  = LNhy / magLN;
    LNhz  = LNhz / magLN;

    aI2P = atan2( LNhy, LNhx );
    bI2P = acos( LNhz );
    gI2P = -phiMod[i];
    LframeEx[0] =  cos(aI2P)*cos(bI2P)*cos(gI2P) - sin(aI2P)*sin(gI2P);
    LframeEx[1] =  sin(aI2P)*cos(bI2P)*cos(gI2P) + cos(aI2P)*sin(gI2P);
    LframeEx[2] = -sin(bI2P)*cos(gI2P);
    LframeEy[0] = -cos(aI2P)*cos(bI2P)*sin(gI2P) - sin(aI2P)*cos(gI2P);
    LframeEy[1] = -sin(aI2P)*cos(bI2P)*sin(gI2P) + cos(aI2P)*cos(gI2P);
    LframeEy[2] =  sin(bI2P)*sin(gI2P);
    LframeEz[0] =  LNhx;
    LframeEz[1] =  LNhy;
    LframeEz[2] =  LNhz;
    aP2J = atan2(JframeEz[0]*LframeEy[0]+JframeEz[1]*LframeEy[1]+JframeEz[2]*LframeEy[2],
                 JframeEz[0]*LframeEx[0]+JframeEz[1]*LframeEx[1]+JframeEz[2]*LframeEx[2]);
    bP2J = acos( JframeEz[0]*LframeEz[0]+JframeEz[1]*LframeEz[1]+JframeEz[2]*LframeEz[2]);
    gP2J = atan2(  JframeEy[0]*LframeEz[0]+JframeEy[1]*LframeEz[1]+JframeEy[2]*LframeEz[2],
                 -(JframeEx[0]*LframeEz[0]+JframeEx[1]*LframeEz[1]+JframeEx[2]*LframeEz[2]));

if (i==0||i==1900) printf("{{%f,%f,%f},{%f,%f,%f},{%f,%f,%f}}\n",JframeEx[0],JframeEx[1],JframeEx[2],JframeEy[0],JframeEy[1],JframeEy[2],JframeEz[0],JframeEz[1],JframeEz[2]);
if (i==0||i==1900) printf("{{%f,%f,%f},{%f,%f,%f},{%f,%f,%f}}\n",LframeEx[0],LframeEx[1],LframeEx[2],LframeEy[0],LframeEy[1],LframeEy[2],LframeEz[0],LframeEz[1],LframeEz[2]);
if (i==0||i==1900) printf("YP: study time = %f\n",i*deltaT/mTScaled);
if (i==1900) printf("YP: gamma: %f, %f, %f, %f\n", JframeEy[0]*LframeEz[0]+JframeEy[1]*LframeEz[1]+JframeEy[2]*LframeEz[2], JframeEx[0]*LframeEz[0]+JframeEx[1]*LframeEz[1]+JframeEx[2]*LframeEz[2], gP2J, atan2(-0.365446,-0.378524));
    /* I2P Euler angles are stored only for debugging purposes */
    alphaI2PTS->data->data[i] = aI2P;
     betaI2PTS->data->data[i] = bI2P;
    gammaI2PTS->data->data[i] = gI2P;
    alphaP2JTS->data->data[i] = -gP2J;
     betaP2JTS->data->data[i] = 0.*LAL_PI/2.-bP2J;
    gammaP2JTS->data->data[i] = -aP2J;

    /* Calculate the value of the Hamiltonian */
    cartPosVec.data[0] = values->data[0];
    cartPosVec.data[1] = values->data[1];
    cartPosVec.data[2] = values->data[2];
    cartMomVec.data[0] = values->data[3];
    cartMomVec.data[1] = values->data[4];
    cartMomVec.data[2] = values->data[5];

    /*if (i == 287)
    {
      printf("%f, %f %f %f, %f %f %f\n",eta, cartPosVec.data[0], cartPosVec.data[1], cartPosVec.data[2], cartMomVec.data[0], cartMomVec.data[1], cartMomVec.data[2]);
    }*/
    ham = XLALSimIMRSpinEOBHamiltonian( eta, &cartPosVec, &cartMomVec,
                  &s1VecOverMtMt, &s2VecOverMtMt,
                  sigmaKerr, sigmaStar, seobParams.tortoise, &seobCoeffs );

    if ( XLALSimIMRSpinEOBGetSpinFactorizedWaveform( &hLM, values, v, ham, 2, 2, &seobParams )
           == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }
    if ( XLALSimIMRSpinEOBNonQCCorrection( &hNQC, values, omega, &nqcCoeffs ) == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }
    hLM *= hNQC;
    h22TS->data->data[i]  = hLM;
    h2m2TS->data->data[i] = conjl(hLM);
    if ( XLALSimIMRSpinEOBGetSpinFactorizedWaveform( &hLM, values, v, ham, 2, 1, &seobParams )
           == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }
    h21TS->data->data[i]  = hLM;
    h2m1TS->data->data[i] = conjl(hLM);
    h20TS->data->data[i]  = 0.0;

    if (i == 95 )
    {
      printf("%.16e %.16e %.16e\n",ham, omega, v);
      printf("%.16e %.16e %.16e\n",values->data[0],values->data[1],values->data[2]);
      printf("%.16e %.16e %.16e\n",values->data[3],values->data[4],values->data[5]);
      printf("%.16e %.16e %.16e\n",values->data[6],values->data[7],values->data[8]);
      printf("%.16e %.16e %.16e\n",values->data[9],values->data[10],values->data[11]);
      printf("%.16e %.16e %.16e %.16e %.16e %.16e\n",nqcCoeffs.a1,nqcCoeffs.a2,nqcCoeffs.a3,nqcCoeffs.a3S,nqcCoeffs.a4,nqcCoeffs.a5);
      printf("%.16e %.16e %.16e %.16e\n",nqcCoeffs.b1,nqcCoeffs.b2,nqcCoeffs.b3,nqcCoeffs.b4);
      printf("%.16e %.16e, %.16e %.16e\n",creal(hLM),cimag(hLM),creal(hNQC),cimag(hNQC));
    }

    /*hPlusTS->data->data[i]  = - 0.5 * amp * cos( 2.*vphi[i]) * cos(2.*alpha) * (1. + LNhz*LNhz)
                            + amp * sin(2.*vphi[i]) * sin(2.*alpha)*LNhz;

    hCrossTS->data->data[i] = - 0.5 * amp * cos( 2.*vphi[i]) * sin(2.*alpha) * (1. + LNhz*LNhz)
                            - amp * sin(2.*vphi[i]) * cos(2.*alpha) * LNhz;*/

    fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e ",
             i*deltaT/mTScaled, aI2P, bI2P, gI2P, aP2J, bP2J, gP2J );
    fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
             vX, vY, vZ, LNhx, LNhy, LNhz, creal(hLM), cimag(hLM) );
  }
  fclose( out );
  printf("YP: quasi-nonprecessing modes generated.\n");

  /* WaveStep 2.4: add quasi-nonprecessing spherical harmonic modes to the SphHarmTimeSeries structure */
  hlmPTS = XLALSphHarmTimeSeriesAddMode( hlmPTS, h22TS, 2, 2 );
  hlmPTS = XLALSphHarmTimeSeriesAddMode( hlmPTS, h21TS, 2, 1 );
  hlmPTS = XLALSphHarmTimeSeriesAddMode( hlmPTS, h20TS, 2, 0 );
  hlmPTS = XLALSphHarmTimeSeriesAddMode( hlmPTS, h2m1TS, 2, -1 );
  hlmPTS = XLALSphHarmTimeSeriesAddMode( hlmPTS, h2m2TS, 2, -2 );
  XLALSphHarmTimeSeriesSetTData( hlmPTS, tlist );

  h22PTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 2 );
  h21PTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 1 );
  h20PTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 0 );
  h2m1PTS = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, -1);
  h2m2PTS = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, -2);
  printf("YP: SphHarmTS structures populated.\n");

  out = fopen( "PWaves.dat", "w" );
  for ( i = 0; i < retLen; i++ )
  {
    fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
             i*deltaT/mTScaled, creal(h22PTS->data->data[i]), cimag(h22PTS->data->data[i]),
                                creal(h21PTS->data->data[i]), cimag(h21PTS->data->data[i]),
                                creal(h20PTS->data->data[i]), cimag(h20PTS->data->data[i]),
                                creal(h2m1PTS->data->data[i]), cimag(h2m1PTS->data->data[i]),
                                creal(h2m2PTS->data->data[i]), cimag(h2m2PTS->data->data[i]) );
  }
  fclose( out );
  printf("YP: P-frame waveforms written to file.\n");


  /* WaveStep 3
   * Generate IMR waveforms in the merger J (~ final spin) -frame
   */
  if ( XLALSimInspiralPrecessionRotateModes( hlmPTS, alphaP2JTS, betaP2JTS, gammaP2JTS ) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  h22JTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 2 );
  h21JTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 1 );
  h20JTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 0 );
  h2m1JTS = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, -1);
  h2m2JTS = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, -2);
  printf("YP: PtoJ rotation done.\n");

  out = fopen( "JWaves.dat", "w" );
  for ( i = 0; i < retLen; i++ )
  {
    fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
             i*deltaT/mTScaled, creal(h22JTS->data->data[i]), cimag(h22JTS->data->data[i]),
                                creal(h21JTS->data->data[i]), cimag(h21JTS->data->data[i]),
                                creal(h20JTS->data->data[i]), cimag(h20JTS->data->data[i]),
                                creal(h2m1JTS->data->data[i]), cimag(h2m1JTS->data->data[i]),
                                creal(h2m2JTS->data->data[i]), cimag(h2m2JTS->data->data[i]) );
  }
  fclose( out );
  printf("YP: P-frame waveforms written to file.\n");

  REAL8Vector *rdMatchPoint = XLALCreateREAL8Vector( 3 );
  REAL8Vector *sigReHi = NULL, *sigImHi = NULL;
  sigReHi  = XLALCreateREAL8Vector( retLen + retLenRDPatch );
  sigImHi  = XLALCreateREAL8Vector( retLen + retLenRDPatch );
  if ( !sigReHi || !sigImHi || !rdMatchPoint )
  {
    XLAL_ERROR( XLAL_ENOMEM );
  }

  if ( combSize > tAttach )
  {
    XLALPrintError( "The comb size looks to be too big!!!\n" );
  }
  rdMatchPoint->data[0] = combSize < tAttach ? tAttach - combSize : 0;
  rdMatchPoint->data[1] = tAttach;
  rdMatchPoint->data[2] = (retLen-1)*deltaT/mTScaled;
  #if debugOutput
  printf("YP::comb range: %f, %f\n",rdMatchPoint->data[0],rdMatchPoint->data[1]);
  #endif
  rdMatchPoint->data[0] -= fmod( rdMatchPoint->data[0], deltaT/mTScaled );
  rdMatchPoint->data[1] -= fmod( rdMatchPoint->data[1], deltaT/mTScaled );


  for ( k = 2; k > -3; k-- )
  {
    hJTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, k );
    memset( sigReHi->data, 0, sigReHi->length * sizeof( sigReHi->data[0] ));
    memset( sigImHi->data, 0, sigImHi->length * sizeof( sigImHi->data[0] ));
    for ( i = 0; i < retLen; i++ )
    {
      sigReHi->data[i] = creal(hJTS->data->data[i]);
      sigImHi->data[i] = cimag(hJTS->data->data[i]);
    }
    if ( XLALSimIMREOBHybridAttachRingdown( sigReHi, sigImHi, 2, k,
                deltaT, m1, m2, spin1[0], spin1[1], spin1[2], spin2[0], spin2[1], spin2[2],
                &tVec, rdMatchPoint, spinEOBApproximant )
            == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }
    for ( i = 0; i < retLen + retLenRDPatch; i++ )
    {
      hIMRJTS->data->data[i] = (sigReHi->data[i]) + I * (sigImHi->data[i]);
    }
    hIMRlmJTS = XLALSphHarmTimeSeriesAddMode( hIMRlmJTS, hIMRJTS, 2, k );
  }
  XLALSphHarmTimeSeriesSetTData( hIMRlmJTS, tlistRDPatch );
  printf("YP: J wave RD attachment done.\n");

  hIMR22JTS  = XLALSphHarmTimeSeriesGetMode( hIMRlmJTS, 2, 2 );
  hIMR21JTS  = XLALSphHarmTimeSeriesGetMode( hIMRlmJTS, 2, 1 );
  hIMR20JTS  = XLALSphHarmTimeSeriesGetMode( hIMRlmJTS, 2, 0 );
  hIMR2m1JTS = XLALSphHarmTimeSeriesGetMode( hIMRlmJTS, 2, -1);
  hIMR2m2JTS = XLALSphHarmTimeSeriesGetMode( hIMRlmJTS, 2, -2);
  out = fopen( "JIMRWaves.dat", "w" );
  for ( i = 0; i < retLen + retLenRDPatch; i++ )
  {
    fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
             i*deltaT/mTScaled, creal(hIMR22JTS->data->data[i]), cimag(hIMR22JTS->data->data[i]),
                                creal(hIMR21JTS->data->data[i]), cimag(hIMR21JTS->data->data[i]),
                                creal(hIMR20JTS->data->data[i]), cimag(hIMR20JTS->data->data[i]),
                                creal(hIMR2m1JTS->data->data[i]), cimag(hIMR2m1JTS->data->data[i]),
                                creal(hIMR2m2JTS->data->data[i]), cimag(hIMR2m2JTS->data->data[i]) );
  }
  fclose( out );
  printf("YP: IMR J wave written to file.\n");


abort();
  /* Point the output pointers to the relevant time series and return */
  (*hplus)  = hPlusTS;
  (*hcross) = hCrossTS;


  return XLAL_SUCCESS;
}


#if defined(SEQUENCE_CLASS)
    #include "MrServers/MrImaging/seq/SeqDebug.h" // TRACE_PUT macros
    #include "MrServers/MrMeasSrv/MeasUtils/MeasMath.h" // M_PI
    #include "MrServers/MrImaging/seq/a_BEAT_physics/spiral_waveform/src/spiral_waveform.h"
#elif defined(GADGETRON)
    #include <cmath>
    #include "gadgetron/log.h"
    #include "spiral_waveform.h"
#else
    #include <cmath>
    #include <cstdio>
    #include "spiral_waveform.h"
#endif

#if defined(SEQUENCE_CLASS)
    #define INFO_P1(fmt, p1) \
    { \
        TRACE_PUT1(TC_INFO, TF_SEQ, fmt, p1); \
    }

    #define INFO_P2(fmt, p1, p2) \
    { \
        TRACE_PUT2(TC_INFO, TF_SEQ, fmt, p1, p2); \
    }

    #define INFO_P3(fmt, p1, p2, p3) \
    { \
        TRACE_PUT3(TC_INFO, TF_SEQ, fmt, p1, p2, p3); \
    }
#elif defined(GADGETRON)
    #define INFO_P1(fmt, p1) \
    { \
        GINFO(fmt, p1); \
    }

    #define INFO_P2(fmt, p1, p2) \
    { \
        GINFO(fmt, p1, p2); \
    }

    #define INFO_P3(fmt, p1, p2, p3) \
    { \
        GINFO(fmt, p1, p2, p3); \
    }
#else
    #define INFO_P1(fmt, p1) \
    { \
        printf(fmt, p1); \
    }

    #define INFO_P2(fmt, p1, p2) \
    { \
        printf(fmt, p1, p2); \
    }

    #define INFO_P3(fmt, p1, p2, p3) \
    { \
        printf(fmt, p1, p2, p3); \
    }
#endif


inline double min(double a, double b)
{
    return (a > b) ? b : a;
}

inline double max(double a, double b)
{
    return (a < b) ? b : a;
}

inline double lerp(double a, double b, double t) 
{
    return (1.0 - t) * a + t * b;
}

inline long ceildiv(long a, long b)
{
    return (a + b - 1) / b;
}


SpiralWaveform::SpiralWaveform()
  : m_lBaseResolution   (-1)
  , m_lSpiralArms       (-1)
  , m_lImagesPerSlab    (-1)
  , m_dFieldOfView      (-1.0)
  , m_dSlabThickness    (-1.0)
  , m_dMaxGradAmpl      (-1.0)
  , m_dMinRiseTime      (-1.0)
  , m_dDwellTime        (-1.0)
  , m_dReadoutOS        (2.0)
  , m_dGradDelay        (0.0)
  , m_dLarmorConst      (SWF_GAMMA_1H)
  , m_eSpiralType       (Arch)
  , m_eVDType           (Linear)
  , m_dVDInnerCutoff    (1.0)
  , m_dVDOuterCutoff    (1.0)
  , m_dVDOuterDensity   (1.0)
  , m_bSloppy           (false)
  , m_dSloppyPeriod     (-1.0)
  , m_dCRT              (0.0)
  , m_lGradSize         (-1)
  , m_dGradAmpl         (-1.0)
  , m_pdGradRO          (NULL)
  , m_pdGradPE          (NULL)
  , m_pdGradSS          (NULL)
  , m_lTrajSize         (-1)
  , m_lSampToSkip       (-1)
  , m_pdTrajRO          (NULL)
  , m_pdTrajPE          (NULL)
  , m_pdTrajSS          (NULL)
  , m_lRampUpTime       (-1)
  , m_lRampDownTime     (-1)
  , m_lReadOutTime      (-1)
  , m_lTotalTime        (-1)
{
    m_dCRT = static_cast<double>(SWF_GRAD_RASTER_TIME) / static_cast<double>(SWF_RASTER_OVERSAMPLING);

    // allocate memory for gradient waveforms
    m_pdGradRO = new double[SWF_MAX_WAVEFORM_SIZE];
    m_pdGradPE = new double[SWF_MAX_WAVEFORM_SIZE];
    m_pdGradSS = new double[SWF_MAX_WAVEFORM_SIZE];

    m_pdTrajRO = new double[SWF_MAX_WAVEFORM_SIZE];
    m_pdTrajPE = new double[SWF_MAX_WAVEFORM_SIZE];
    m_pdTrajSS = new double[SWF_MAX_WAVEFORM_SIZE];
}


SpiralWaveform::~SpiralWaveform()
{
    // release memory
    if (m_pdGradRO)
    {
        delete[] m_pdGradRO;
        m_pdGradRO = NULL;
    }

    if (m_pdGradPE)
    {
        delete[] m_pdGradPE;
        m_pdGradPE = NULL;
    }

    if (m_pdGradSS)
    {
        delete[] m_pdGradSS;
        m_pdGradSS = NULL;
    }

    if (m_pdTrajRO)
    {
        delete[] m_pdTrajRO;
        m_pdTrajRO = NULL;
    }

    if (m_pdTrajPE)
    {
        delete[] m_pdTrajPE;
        m_pdTrajPE = NULL;
    }

    if (m_pdTrajSS)
    {
        delete[] m_pdTrajSS;
        m_pdTrajSS = NULL;
    }
}


bool SpiralWaveform::calculate(bool bCalcTraj)
{
    static const char *ptModule = {"SpiralWaveform::calculate"};

    // * ---------------------------------------------------------------------- *
    // * Check inputs                                                           *
    // * ---------------------------------------------------------------------- *

    if (m_lBaseResolution <= 0)
    {
        INFO_P1("%s base resolution not set\n", ptModule);
        return false;
    }

    if (m_lSpiralArms <= 0)
    {
        INFO_P1("%s spiral arms not set\n", ptModule);
        return false;
    }

    if (m_dFieldOfView <= 0.0)
    {
        INFO_P1("%s field of view not set\n", ptModule);
        return false;
    }

    if (m_dMaxGradAmpl <= 0.0)
    {
        INFO_P1("%s max grad ampl not set\n", ptModule);
        return false;
    }

    if (m_dMinRiseTime <= 0.0)
    {
        INFO_P1("%s min rise time not set\n", ptModule);
        return false;
    }

    if (m_dDwellTime <= 0.0)
    {
        INFO_P1("%s dwell time not set\n", ptModule);
        return false;
    }

    if (m_bSloppy && m_dSloppyPeriod <= 0.0)
    {
        INFO_P1("%s sloppy period not set\n", ptModule);
        return false;
    }

    // member variables are stored with the most commonly used units
    // here we do some conversions to use the same units in the computation
    // we also calculate some additional derived values
    double dCRT             = 1.e-3 * m_dCRT;         // [ms]
    double dFieldOfView     = 1.e-3 * m_dFieldOfView; // [m]
    double dPixelSize       = dFieldOfView / static_cast<double>(m_lBaseResolution);
    double dSlabThickness   = 1.e-3 * m_dSlabThickness;
    double dSliceThickness  = dSlabThickness / static_cast<double>(m_lImagesPerSlab);

    // compute gradient amplitude according to the dwell time
    double dEffDwellTime    = m_dDwellTime * m_dReadoutOS;  // [us]
    double dBandwidth       = 1.e6 / dEffDwellTime;         // [Hz]
    m_dGradAmpl             = dBandwidth / (m_dLarmorConst * m_dFieldOfView);
    if (m_dGradAmpl > m_dMaxGradAmpl)
    {
        INFO_P1("%s spiral waveform would exceed gradient amplitude limits", ptModule);
        return false;  
    }

    // the next 4 variables are for variable density spirals
    // they create a transition in the radial spacing as the k-space radius goes from 0 to 1, i.e.
    //    0 < kr < m_dVDInnerCutoff : spacing = Nyquist distance
    // m_dVDInnerCutoff < kr < m_dVDOuterCutoff : spacing increases to m_dVDOuterDensity (affected by m_eVDType)
    // m_dVDOuterCutoff < kr < 1    : spacing = m_dVDOuterDensity

    // radial distance per arm to meet the Nyquist limit
    double dNyquist         = static_cast<double>(m_lSpiralArms) / dFieldOfView;
    // gamma times raster time (dGammaCRT * g = dk)
    double dGammaCRT        = m_dLarmorConst * dCRT;
    // the most the gradients can change in 1 raster period
    double dMaxGradDeltaCRT = m_dMaxSlewRate * dCRT;
    // GRT
    double dGammaGRT        = SWF_RASTER_OVERSAMPLING * dGammaCRT;

    double dRadialSpacing = 1.0;
    double dAlpha, dPhi, dTheta;
    double kr, kmx, kmy, kmz, kmr, rnorm;
    double dDirRO=0, dDirPE=0, dDirSS=0, umag;
    double gx=0, gy=0, gz=0;
    double us_i;
    double dGradMag=0,term;
    double krmax, kzmax, krmax2, kzmax2;
    double krlim;
    long i, i0, i1, i_end;
    long lJ;

    // allocate memory for t
    double* kx    = new double[SWF_RASTER_OVERSAMPLING * SWF_MAX_WAVEFORM_SIZE];
    double* ky    = new double[SWF_RASTER_OVERSAMPLING * SWF_MAX_WAVEFORM_SIZE];
    double* kz    = new double[SWF_RASTER_OVERSAMPLING * SWF_MAX_WAVEFORM_SIZE];
    double* gsign = new double[SWF_RASTER_OVERSAMPLING * SWF_MAX_WAVEFORM_SIZE];

    // initialize values
    for (i = 0; i < SWF_RASTER_OVERSAMPLING * SWF_MAX_WAVEFORM_SIZE; i++) gsign[i] = 1.;
    for (i = 0; i < SWF_RASTER_OVERSAMPLING * SWF_MAX_WAVEFORM_SIZE; i++) kx[i] = 0.;
    for (i = 0; i < SWF_RASTER_OVERSAMPLING * SWF_MAX_WAVEFORM_SIZE; i++) ky[i] = 0.;
    for (i = 0; i < SWF_RASTER_OVERSAMPLING * SWF_MAX_WAVEFORM_SIZE; i++) kz[i] = 0.;
    for (i = 0; i < SWF_MAX_WAVEFORM_SIZE; i++) m_pdGradRO[i] = 0.;
    for (i = 0; i < SWF_MAX_WAVEFORM_SIZE; i++) m_pdGradPE[i] = 0.;
    for (i = 0; i < SWF_MAX_WAVEFORM_SIZE; i++) m_pdGradSS[i] = 0.;

    krmax = 0.5/dPixelSize;
    kzmax = 0.5/dSliceThickness;
    krmax2 = krmax*krmax;
    kzmax2 = kzmax*kzmax;
    krlim = krmax*(1.-(dPixelSize/dFieldOfView));

    // start out spiral going radially at max slew-rate for 2 time-points
    kx[0] = 0.;
    ky[0] = 0.;
    kx[1] = dGammaCRT * dMaxGradDeltaCRT;
    ky[1] = 0.;
    kx[2] = 3. * dGammaCRT * dMaxGradDeltaCRT;
    ky[2] = 0.;

    // if spherical DST
    if (m_eSpiralType == SphDST)
    {
        kz[0] = kzmax;
        // stay on surface of ellipsoid
        kz[1] = sqrt(kzmax2 * (1. - ((kx[1] * kx[1] + ky[1] * ky[1]) / krmax2)));
        kz[2] = sqrt(kzmax2 * (1. - ((kx[2] * kx[2] + ky[2] * ky[2]) / krmax2)));
    }

    i = 2;
    kr = kx[2];

    /******************************/
    /* LOOP UNTIL YOU HIT MAX RES */
    /******************************/
    while (kr <= krlim)
    {
        if (i >= SWF_RASTER_OVERSAMPLING * SWF_MAX_WAVEFORM_SIZE - 1)
        {
            INFO_P1("%s spiral waveform too long\n", ptModule);
            return false;
        }

        /**************************/
        /*** STEP 1:  Determine the direction (dDirRO,dDirPE) of the gradient at ~(i+0.5) */
        /**************************/
        /* calculate dk/dCRT = m_dLarmorConst G*/
        kmx = 1.5*kx[i] - 0.5*kx[i-1];
        kmy = 1.5*ky[i] - 0.5*ky[i-1];
        kmr = sqrt(kmx*kmx + kmy*kmy);

        /////////////////////////////
        // Start dRadialSpacing logic //
        /////////////////////////////
        rnorm = 2.*dPixelSize*kmr; /* the k-space radius, normalized to go from 0 to 1 */
        
        // determine the undersample factor
        if (rnorm <= m_dVDInnerCutoff)
        {
            dRadialSpacing = 1;
        }
        else if (rnorm < m_dVDOuterCutoff)
        {
            us_i = (rnorm - m_dVDInnerCutoff) / (m_dVDOuterCutoff - m_dVDInnerCutoff); /* goes from 0 to 1 as rnorm goes from m_dVDInnerCutoff to m_dVDOuterCutoff*/
            if (m_eVDType == Linear)
            {
                dRadialSpacing = 1. + (m_dVDOuterDensity - 1.) * us_i;
            }
            else if (m_eVDType == Quadratic)
            {
                dRadialSpacing = 1. + (m_dVDOuterDensity - 1.) * us_i * us_i;
            }
            else if (m_eVDType == Hanning)
            {
                dRadialSpacing = 1. + (m_dVDOuterDensity - 1.) * 0.5 * (1. - cos(us_i * M_PI));
            }
        }
        else // rnorm > m_dVDOuterCutoff
        {
            dRadialSpacing = m_dVDOuterDensity;
        } 

        // Undersample spiral for Spherical-Distributed Spiral
        if (m_eSpiralType == SphDST)
        {
            if (rnorm < 1.0)
                dRadialSpacing = min(dSlabThickness / dSliceThickness, dRadialSpacing / sqrt(1.0 - (rnorm * rnorm)));
            else
                dRadialSpacing = dSlabThickness / dSliceThickness;
        }

        // MAKE FERMAT SPIRAL FOR FLORET
        if (m_eSpiralType == Fermat && rnorm > 0.)
            dRadialSpacing *= 1.0 / rnorm;

    /* Sloppy Spirals - add variability to dRadialSpacing for reduced aliasing coherence */ 
    // A couple different options here are commented out
    
        if (m_bSloppy)
        {
            // Lots of ways to be sloppy
    //      dRadialSpacing = max(1., (dRadialSpacing + ((dRadialSpacing-1.)*sin(2.*M_PI*(double)(i)/m_dSloppyPeriod))));
    //      dRadialSpacing += (dRadialSpacing-1.)*sin(2.*M_PI*m_dSloppyPeriod*atan2(ky[i],kx[i]));
            dRadialSpacing += (dRadialSpacing-1.)*sin(2.*M_PI*m_dSloppyPeriod*rnorm);
        }

    ///////////////////////////
    // End dRadialSpacing logic //
    ///////////////////////////

    /* See the Key Equation 4 at the beginning of the code */
        dAlpha = atan(2.*M_PI*kmr/(dRadialSpacing*dNyquist));
        dPhi = atan2(kmy,kmx);
        dTheta = dPhi + dAlpha;

        dDirRO = cos(dTheta);
        dDirPE = sin(dTheta);

    // IF SPHERICAL DST
    // u dot km is zero if moving on a sphere (km is radial, u is tangential,
    // thus km stays on the sphere)
    // We are on an ellipsoid, but can normalize u and km by krmax and kzmax to make this work
    // The final gradient vector (dDirRO dDirPE dDirSS) will be tangential to the sphere
        if (m_eSpiralType == SphDST)
        {
            kmz = 1.5*kz[i] - 0.5*kz[i-1];
            dDirSS = -((dDirRO*kmx + dDirPE*kmy)/krmax2)*(kzmax2/kmz);
            umag = sqrt(dDirRO*dDirRO + dDirPE*dDirPE + dDirSS*dDirSS);
            dDirRO = dDirRO/umag;
            dDirPE = dDirPE/umag;
            dDirSS = dDirSS/umag;
            gz = (kz[i] - kz[i-1])/dGammaCRT;
        }

    /**************************/
    /*** STEP 2: Find largest gradient magnitude with available slew */
    /**************************/

        // Current gradient
        gx = (kx[i] - kx[i-1])/dGammaCRT;
        gy = (ky[i] - ky[i-1])/dGammaCRT;

        // solve for dGradMag using the quadratic equation |dGradMag u - g| = dMaxGradDeltaCRT
        // which is
        //   (dGradMag u - g)(dGradMag u* - g*) = dMaxGradDeltaCRT^2
        // which gives
        //   dGradMag^2 (u u*) - dGradMag (g u* + u g*) + g g* - dMaxGradDeltaCRT^2 = 0

        // Replacing u u* with 1 (i.e. u is a unit vector) and
        // replacing (g u* + u g*) with 2 Real[g u*]
        // this is
        //   dGradMag^2 + dGradMag (2 b) + c = 0
        // giving
        //   dGradMag = -b +/- Sqrt(b^2 - c)
        // The variable "term" = (b^2 - c) will be positive if we can meet the desired new gradient 
        term = dMaxGradDeltaCRT * dMaxGradDeltaCRT - (gx * gx + gy * gy + gz * gz) + (dDirRO * gx + dDirPE * gy + dDirSS * gz) * (dDirRO * gx + dDirPE * gy + dDirSS * gz);

        if (term >= 0)
        {
        // Slew constraint is met! Now assign next gradient and then next k value
        // NOTE gsign is +1 or -1
        //   if gsign is positive, we are using slew to speed up (increase dGradMag) as much as possible
        //   if gsign is negative, we are using slew to slow down (decrease dGradMag) as much as possible
            dGradMag  = min((dDirRO * gx + dDirPE * gy + dDirSS * gz) + gsign[i] * sqrt(term), m_dGradAmpl);
            gx = dGradMag*dDirRO;
            gy = dGradMag*dDirPE;

            kx[i+1] = kx[i] + gx*dGammaCRT;
            ky[i+1] = ky[i] + gy*dGammaCRT;

        // If SPHERE
            if (m_eSpiralType == SphDST)
                kz[i+1] = sqrt(kzmax2*(1.-((kx[i+1]*kx[i+1]+ky[i+1]*ky[i+1])/krmax2))); // stay on surface of ellipsoid

            i++;
        } // term >= 0
        else
        {
            // We can't go further without violating the slew rate
            // This means that we've sped up too fast to turn here at the desired curvature
            // We are going to iteratively go back in time and slow down, rather than speed up, at max slew
            // Here we'll keep looking back until gsign is positive, then add another negative gsign, just far enough to make the current corner
            while ((i > 3) && (gsign[i - 1] == -1)) i--;
            gsign[i - 1] = -1;
            i = i - 2;
        } // term < 0

        kr = sqrt(kx[i]*kx[i] + ky[i]*ky[i]);

    } // MAIN kr loop

    i_end = i;

    //********************************************
    // DONE LOOPING FOR SAMPLING PORTION
    // recast k to g while subsampling by SWF_RASTER_OVERSAMPLING  
    //********************************************
    m_pdGradRO[0] = 0.;
    m_pdGradPE[0] = 0.; 
    m_pdGradSS[0] = 0.; 
    double dSumGradRO = 0.;
    double dSumGradPE = 0.;
    double dSumGradSS = 0.;
    for (lJ = 1; lJ <= (i_end / SWF_RASTER_OVERSAMPLING); lJ++)
    {
        i1 = lJ * SWF_RASTER_OVERSAMPLING;
        i0 = (lJ - 1) * SWF_RASTER_OVERSAMPLING;
        m_pdGradRO[lJ] = (kx[i1] - kx[i0]) / dGammaGRT;
        m_pdGradPE[lJ] = (ky[i1] - ky[i0]) / dGammaGRT;
        m_pdGradSS[lJ] = (kz[i1] - kz[i0]) / dGammaGRT;
        dSumGradRO += m_pdGradRO[lJ];
        dSumGradPE += m_pdGradPE[lJ];
        dSumGradSS += m_pdGradSS[lJ];
    }
    // readout time
    m_lGradSize = lJ;
    m_lReadOutTime = lJ * SWF_GRAD_RASTER_TIME;

    // recalculate these ending gradient points
    dGradMag = sqrt(m_pdGradRO[m_lGradSize - 1] * m_pdGradRO[m_lGradSize - 1] +
              m_pdGradPE[m_lGradSize - 1] * m_pdGradPE[m_lGradSize - 1] +
              m_pdGradSS[m_lGradSize - 1] * m_pdGradSS[m_lGradSize - 1]);
    dDirRO = m_pdGradRO[m_lGradSize - 1] / dGradMag;
    dDirPE = m_pdGradPE[m_lGradSize - 1] / dGradMag;
    dDirSS = m_pdGradSS[m_lGradSize - 1] / dGradMag;

    
    // * ---------------------------------------------------------------------- *
    // * Prepare ramp down                                                      *
    // * ---------------------------------------------------------------------- *

    // Ramp gradients to zero
    long lRampDownSize = 0;
    double dGradDelta = 0.0;
    double dMaxGradDelta = m_dMaxSlewRate * 1.0e-3 * SWF_GRAD_RASTER_TIME;

    if (m_lRampDownTime >= 0)
    {
        // A specific ramp time has been requested. Compute the necessary
        // gradient increment to realize this ramp time.

        if (m_lRampDownTime == 0)
        {
            INFO_P2("%s cannot realize requested ramp down time: %ld\n", ptModule, m_lRampDownTime );
            return false;
        }
        
        lRampDownSize = m_lRampDownTime / SWF_GRAD_RASTER_TIME;
        dGradDelta = dGradMag / lRampDownSize; // [(mT/m) / sample]

        if (dGradDelta > dMaxGradDelta)
        {
            INFO_P3("%s cannot realize requested ramp down time = %ld, because slew rate limit = %f (mT/m)/ms would be exceeded\n", ptModule, m_lRampDownTime, m_dMaxSlewRate );
            return false;
        }
    }
    else
    {
        // No specific ramp down time requested. Ramp down at the maximum slew
        // rate.
        dGradDelta = dMaxGradDelta;

        // Compute time required for the ramp.
        lRampDownSize = static_cast<long>(ceil(dGradMag / dGradDelta) + 0.5);
        m_lRampDownTime = lRampDownSize * SWF_GRAD_RASTER_TIME;

        // Update delta to match this duration.
        dGradDelta = dGradMag / lRampDownSize;
    }

    for (long lR = lRampDownSize - 1; lR >= 0; lJ++, lR--)
    {
        if (lJ >= SWF_MAX_WAVEFORM_SIZE)
        {
            INFO_P1("%s spiral waveform too long\n", ptModule);
            return false;
        }

        m_pdGradRO[lJ] = dGradDelta * dDirRO * lR;
        m_pdGradPE[lJ] = dGradDelta * dDirPE * lR;
        m_pdGradSS[lJ] = dGradDelta * dDirSS * lR;
        dSumGradRO += m_pdGradRO[lJ];
        dSumGradPE += m_pdGradPE[lJ];
        dSumGradSS += m_pdGradSS[lJ];
    }

    m_lGradSize += lRampDownSize;

    
    // * ---------------------------------------------------------------------- *
    // * Calculate arbitrary gradient waveforms, if necessary                   *
    // * ---------------------------------------------------------------------- *

    // normalize waveforms
    for (lJ = 0; lJ < m_lGradSize; lJ++)
    {
        m_pdGradRO[lJ] /= m_dGradAmpl;
        m_pdGradPE[lJ] /= m_dGradAmpl;
        m_pdGradSS[lJ] /= m_dGradAmpl;
    }

    // ramp-up time is always 0 for now, as we only consider spiral-out
    // trajectories
    m_lRampUpTime = 0;

    // total time
    m_lTotalTime = m_lRampUpTime + m_lReadOutTime + m_lRampDownTime;


    // * ---------------------------------------------------------------------- *
    // * Release memory                                                         *
    // * ---------------------------------------------------------------------- *

    delete[] kx;
    delete[] ky;
    delete[] kz;
    delete[] gsign;


    // * ---------------------------------------------------------------------- *
    // * Calculate k-space trajectory                                           *
    // * ---------------------------------------------------------------------- *

    long lEffDwellTime = static_cast<long>(m_dDwellTime * m_dReadoutOS * 1000.0 + 0.5); // [ns]
    m_lTrajSize = ceildiv(m_lReadOutTime * 1000, lEffDwellTime) * static_cast<long>(m_dReadoutOS + 0.5);
    
    if (!bCalcTraj)
    {
        // trajectory calculation was not requested, so we are done
        return true;
    }

    // compute trajectory at gradient raster grid, integrating with trapezoid
    // rule
    double* pdTrajRO = new double[m_lGradSize];
    double* pdTrajPE = new double[m_lGradSize];
    double* pdTrajSS = new double[m_lGradSize];

    pdTrajRO[0] = 0.0;
    pdTrajPE[0] = 0.0;
    pdTrajSS[0] = 0.0;

    double dFactor = m_dGradAmpl * m_dLarmorConst * SWF_GRAD_RASTER_TIME * 1e-6;    
    for (lJ = 1; lJ < m_lGradSize; lJ++)
    {
        pdTrajRO[lJ] = pdTrajRO[lJ - 1] + dFactor * (m_pdGradRO[lJ] + m_pdGradRO[lJ - 1]) / 2.0;
        pdTrajPE[lJ] = pdTrajPE[lJ - 1] + dFactor * (m_pdGradPE[lJ] + m_pdGradPE[lJ - 1]) / 2.0;
        pdTrajSS[lJ] = pdTrajSS[lJ - 1] + dFactor * (m_pdGradSS[lJ] + m_pdGradSS[lJ - 1]) / 2.0;
    }

    // convert trajectory from 1/mm (range +- 1 / (2 * pixel size)) to
    // radians/pixel (range +- pi)
    double dPixelSizeRO = m_dFieldOfView / m_lBaseResolution;
    double dPixelSizePE = m_dFieldOfView / m_lBaseResolution;
    double dPixelSizeSS = m_dSlabThickness / m_lImagesPerSlab;
    for (lJ = 0; lJ < m_lGradSize; lJ++)
    {
        pdTrajRO[lJ] *= dPixelSizeRO * 2.0 * M_PI;
        pdTrajPE[lJ] *= dPixelSizePE * 2.0 * M_PI;
        pdTrajSS[lJ] *= dPixelSizeSS * 2.0 * M_PI;
    }

    // now interpolate to ADC raster
    // start and spacing times for gradient and ADC
    long lT = 0, lG = 0;
    double dTime;
    double dDeltaTimeGrad = SWF_GRAD_RASTER_TIME;
    double dDeltaTimeADC = m_dDwellTime;
    double dCurrentTimeGrad = 0.0;
    double dCurrentTimeADC = -m_dGradDelay;
    m_lSampToSkip = 0;
    
    for (lT = 0; lT < m_lTrajSize; lT++, dCurrentTimeADC += dDeltaTimeADC)
    {
        if (dCurrentTimeADC < 0.0)
        {
            m_lSampToSkip++;
            m_pdTrajRO[lT] = 0.0;
            m_pdTrajPE[lT] = 0.0;
            m_pdTrajSS[lT] = 0.0;
            continue;
        }

        // find previous gradient point
        while (dCurrentTimeGrad + dDeltaTimeGrad < dCurrentTimeADC)
        {
            lG++;
            dCurrentTimeGrad += dDeltaTimeGrad;
        }

        // how far are we after previous point (0.0 is exactly at previous
        // point, 1.0 is exactly at next point)
        dTime = (dCurrentTimeADC - dCurrentTimeGrad) / dDeltaTimeGrad;
        
        m_pdTrajRO[lT] = lerp(pdTrajRO[lG], pdTrajRO[lG + 1], dTime);
        m_pdTrajPE[lT] = lerp(pdTrajPE[lG], pdTrajPE[lG + 1], dTime);
        m_pdTrajSS[lT] = lerp(pdTrajSS[lG], pdTrajSS[lG + 1], dTime);
    }

    // delete temp data
    delete[] pdTrajRO;
    delete[] pdTrajPE;
    delete[] pdTrajSS;

    return true;
}

long SpiralWaveform::getGradientWaveformSize() const
{
    return m_lGradSize;
}

double SpiralWaveform::getGradientAmplitude() const
{
    return m_dGradAmpl;
}

void SpiralWaveform::getGradientWaveformRO(double* pdWaveform) const
{
    for (long lI = 0; lI < m_lGradSize; lI++)
    {
        pdWaveform[lI] = m_pdGradRO[lI];
    }
}

void SpiralWaveform::getGradientWaveformPE(double* pdWaveform) const
{
    for (long lI = 0; lI < m_lGradSize; lI++)
    {
        pdWaveform[lI] = m_pdGradPE[lI];
    }
}

void SpiralWaveform::getGradientWaveformSS(double* pdWaveform) const
{
    for (long lI = 0; lI < m_lGradSize; lI++)
    {
        pdWaveform[lI] = m_pdGradSS[lI];
    }
}

void SpiralWaveform::getGradientWaveformRO(float* pfWaveform) const
{
    for (long lI = 0; lI < m_lGradSize; lI++)
    {
        pfWaveform[lI] = static_cast<float>(m_pdGradRO[lI]);
    }
}

void SpiralWaveform::getGradientWaveformPE(float* pfWaveform) const
{
    for (long lI = 0; lI < m_lGradSize; lI++)
    {
        pfWaveform[lI] = static_cast<float>(m_pdGradPE[lI]);
    }
}

void SpiralWaveform::getGradientWaveformSS(float* pfWaveform) const
{
    for (long lI = 0; lI < m_lGradSize; lI++)
    {
        pfWaveform[lI] = static_cast<float>(m_pdGradSS[lI]);
    }
}

long SpiralWaveform::getTrajectoryWaveformSize() const
{
    return m_lTrajSize;
}

long SpiralWaveform::getSamplesToSkip() const
{
    return m_lSampToSkip;
}

void SpiralWaveform::getTrajectoryWaveformRO(double* pdWaveform) const
{
    for (long lI = 0; lI < m_lTrajSize; lI++)
    {
        pdWaveform[lI] = m_pdTrajRO[lI];
    }
}

void SpiralWaveform::getTrajectoryWaveformPE(double* pdWaveform) const
{
    for (long lI = 0; lI < m_lTrajSize; lI++)
    {
        pdWaveform[lI] = m_pdTrajPE[lI];
    }
}

void SpiralWaveform::getTrajectoryWaveformSS(double* pdWaveform) const
{
    for (long lI = 0; lI < m_lTrajSize; lI++)
    {
        pdWaveform[lI] = m_pdTrajSS[lI];
    }
}

void SpiralWaveform::getTrajectoryWaveformRO(float* pfWaveform) const
{
    for (long lI = 0; lI < m_lTrajSize; lI++)
    {
        pfWaveform[lI] = static_cast<float>(m_pdTrajRO[lI]);
    }
}

void SpiralWaveform::getTrajectoryWaveformPE(float* pfWaveform) const
{
    for (long lI = 0; lI < m_lTrajSize; lI++)
    {
        pfWaveform[lI] = static_cast<float>(m_pdTrajPE[lI]);
    }
}

void SpiralWaveform::getTrajectoryWaveformSS(float* pfWaveform) const
{
    for (long lI = 0; lI < m_lTrajSize; lI++)
    {
        pfWaveform[lI] = static_cast<float>(m_pdTrajSS[lI]);
    }
}


extern "C"
{
    int calculate_spiral_trajectory(
        float* pfTraj,
        long* plTrajSize,
        long lBaseResolution, 
        long lSpiralArms,
        double dFieldOfView,
        double dMaxGradAmpl,
        double dMinRiseTime,
        double dDwellTime,
        double dReadoutOS,
        double dGradientDelay,
        double dLarmorConst)
    {
        SpiralWaveform wf;
        wf.setBaseResolution(lBaseResolution);
        wf.setSpiralArms(lSpiralArms);
        wf.setFieldOfView(dFieldOfView);
        wf.setMaxGradAmpl(dMaxGradAmpl);
        wf.setMinRiseTime(dMinRiseTime);
        wf.setDwellTime(dDwellTime);
        wf.setReadoutOS(dReadoutOS);
        wf.setGradientDelay(dGradientDelay);
        wf.setLarmorConst(dLarmorConst);

        if (!wf.calculate(true))
        {
            return 1;
        }

        *plTrajSize = wf.getTrajectoryWaveformSize();

        float* pfTrajRO = new float[*plTrajSize];
        float* pfTrajPE = new float[*plTrajSize];

        wf.getTrajectoryWaveformRO(pfTrajRO);
        wf.getTrajectoryWaveformPE(pfTrajPE);

        for (long lI = 0; lI < *plTrajSize; lI++)
        {
            pfTraj[2 * lI + 0] = pfTrajRO[lI];
            pfTraj[2 * lI + 1] = pfTrajPE[lI];
        }

        delete[] pfTrajRO;
        delete[] pfTrajPE;
        
        // float* pfTrajRO = new float[*plTrajSize];
        // float* pfTrajPE = new float[*plTrajSize];
        // float* pfTrajSS = new float[*plTrajSize];

        // wf.getTrajectoryWaveformRO(pfTrajRO);
        // wf.getTrajectoryWaveformPE(pfTrajPE);
        // wf.getTrajectoryWaveformSS(pfTrajSS);

        // for (long lI = 0; lI < *plTrajSize; lI++)
        // {
        //     pfTraj[3 * lI + 0] = pfTrajRO[lI];
        //     pfTraj[3 * lI + 1] = pfTrajPE[lI];
        //     pfTraj[3 * lI + 2] = pfTrajSS[lI];
        // }

        // delete[] pfTrajRO;
        // delete[] pfTrajPE;
        // delete[] pfTrajSS;

        return 0;
    }
}

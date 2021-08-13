#ifndef SPIRAL_WAVEFORM_H_
#define SPIRAL_WAVEFORM_H_

#define SWF_GRAD_RASTER_TIME        10
#define SWF_RASTER_OVERSAMPLING     10
#define SWF_GAMMA_1H                42.577478518
#define SWF_MAX_WAVEFORM_SIZE       100000

class SpiralWaveform
{
public:
    
    // main API
    SpiralWaveform();
    ~SpiralWaveform();

    bool calculate(bool bCalcTraj = false);

    long getGradientWaveformSize() const;
    double getGradientAmplitude() const;

    void getGradientWaveformRO(double* pdWaveform) const;
    void getGradientWaveformPE(double* pdWaveform) const;
    void getGradientWaveformSS(double* pdWaveform) const;

    void getGradientWaveformRO(float* pfWaveform) const;
    void getGradientWaveformPE(float* pfWaveform) const;
    void getGradientWaveformSS(float* pfWaveform) const;

    long getTrajectoryWaveformSize() const;
    long getSamplesToSkip() const;
    // double getTrajectoryAmplitude() const;

    void getTrajectoryWaveformRO(double* pdWaveform) const;
    void getTrajectoryWaveformPE(double* pdWaveform) const;
    void getTrajectoryWaveformSS(double* pdWaveform) const;

    void getTrajectoryWaveformRO(float* pfWaveform) const;
    void getTrajectoryWaveformPE(float* pfWaveform) const;
    void getTrajectoryWaveformSS(float* pfWaveform) const;

    long getRampUpTime() const;
    void setRampUpTime(long lRampUpTime);

    long getRampDownTime() const;
    void setRampDownTime(long lRampDownTime);

    long getReadOutTime() const;
    long getTotalTime() const;

    // types
    enum SpiralType     { Arch, CylDST, SphDST, Fermat };
    enum VDType         { Linear, Quadratic, Hanning };

    // getters/setters
    long getBaseResolution() const;
    void setBaseResolution(long lBaseResolution);

    long getImagesPerSlab() const;
    void setImagesPerSlab(long lImagesPerSlab);

    long getSpiralArms() const;
    void setSpiralArms(long lSpiralArms);

    double getFieldOfView() const;
    void setFieldOfView(double dFieldOfView);

    double getSlabThickness() const;
    void setSlabThickness(double dSlabThickness);

    double getMaxGradAmpl() const;
    void setMaxGradAmpl(double dMaxGradAmpl);

    double getMinRiseTime() const;
    void setMinRiseTime(double dMinRiseTime);

    double getDwellTime() const;
    void setDwellTime(double dDwellTime);

    double getReadoutOS() const;
    void setReadoutOS(double dReadoutOS);

    double getGradientDelay() const;
    void setGradientDelay(double dGradientDelay);

    double getLarmorConst() const;
    void setLarmorConst(double dLarmorConst);

    SpiralType getSpiralType() const;
    void setSpiralType(SpiralType eSpiralType);

    VDType getVDType() const;
    void setVDType(VDType eVDType);

    double getVDInnerCutoff() const;
    void setVDInnerCutoff(double dVDInnerCutoff);

    double getVDOuterCutoff() const;
    void setVDOuterCutoff(double dVDOuterCutoff);

    double getVDOuterDensity() const;
    void setVDOuterDensity(double dVDOuterDensity);

    bool getSloppy() const;
    void setSloppy(bool bSloppy);

    double getSloppyPeriod() const;
    void setSloppyPeriod(double dSloppyPeriod);
    
protected:

    long m_lBaseResolution;
    long m_lSpiralArms;
    long m_lImagesPerSlab;

    double m_dFieldOfView;      // [mm]
    double m_dSlabThickness;    // [mm]
    
    double m_dMaxGradAmpl;      // [mT/m]
    double m_dMinRiseTime;      // [us/(mT/m)]
    double m_dMaxSlewRate;      // [(mT/m)/ms]
    
    double m_dDwellTime;        // [us]
    double m_dReadoutOS;
    double m_dGradDelay;        // [us]
    double m_dLarmorConst;      // [MHz / T]

    SpiralType m_eSpiralType;
    
    // variable density configuration
    VDType m_eVDType;
    double m_dVDInnerCutoff;
    double m_dVDOuterCutoff;
    double m_dVDOuterDensity;

    // sloppy spirals
    bool m_bSloppy;
    double m_dSloppyPeriod;
    
    // calculation raster time
    double m_dCRT;

    // gradient waveforms
    long m_lGradSize;
    double m_dGradAmpl;
    double* m_pdGradRO;
    double* m_pdGradPE;
    double* m_pdGradSS;

    // k-space trajectory
    long m_lTrajSize;
    long m_lSampToSkip;
    double* m_pdTrajRO;
    double* m_pdTrajPE;
    double* m_pdTrajSS;

    // timing
    long m_lRampUpTime;     // [us]
    long m_lRampDownTime;   // [us]
    long m_lReadOutTime;    // [us]
    long m_lTotalTime;      // [us]
};

inline long SpiralWaveform::getBaseResolution() const
{
    return m_lBaseResolution;
}

inline void SpiralWaveform::setBaseResolution(long lBaseResolution)
{
    m_lBaseResolution = lBaseResolution;
}

inline long SpiralWaveform::getSpiralArms() const
{
    return m_lSpiralArms;
}

inline void SpiralWaveform::setSpiralArms(long lSpiralArms)
{
    m_lSpiralArms = lSpiralArms;
}

inline double SpiralWaveform::getFieldOfView() const
{
    return m_dFieldOfView;
}

inline void SpiralWaveform::setFieldOfView(double dFieldOfView)
{
    m_dFieldOfView = dFieldOfView;
}

inline double SpiralWaveform::getSlabThickness() const
{
    return m_dSlabThickness;
}

inline void SpiralWaveform::setSlabThickness(double dSlabThickness)
{
    m_dSlabThickness = dSlabThickness;
}

inline double SpiralWaveform::getMaxGradAmpl() const
{
    return m_dMaxGradAmpl;
}

inline void SpiralWaveform::setMaxGradAmpl(double dMaxGradAmpl)
{
    m_dMaxGradAmpl = dMaxGradAmpl;
}

inline double SpiralWaveform::getMinRiseTime() const
{
    return m_dMinRiseTime;
}

inline void SpiralWaveform::setMinRiseTime(double dMinRiseTime)
{
    m_dMinRiseTime = dMinRiseTime;
    m_dMaxSlewRate = 1.e3 / dMinRiseTime;
}

inline double SpiralWaveform::getDwellTime() const
{
    return m_dDwellTime;
}

inline void SpiralWaveform::setDwellTime(double dDwellTime)
{
    m_dDwellTime = dDwellTime;
}

inline double SpiralWaveform::getReadoutOS() const
{
    return m_dReadoutOS;
}

inline void SpiralWaveform::setReadoutOS(double dReadoutOS)
{
    m_dReadoutOS = dReadoutOS;
}

inline double SpiralWaveform::getGradientDelay() const
{
    return m_dGradDelay;
}

inline void SpiralWaveform::setGradientDelay(double dGradientDelay)
{
    m_dGradDelay = dGradientDelay;
}

inline double SpiralWaveform::getLarmorConst() const
{
    return m_dLarmorConst;
}

inline void SpiralWaveform::setLarmorConst(double dLarmorConst)
{
    m_dLarmorConst = dLarmorConst;
}

inline SpiralWaveform::SpiralType SpiralWaveform::getSpiralType() const
{
    return m_eSpiralType;
}

inline void SpiralWaveform::setSpiralType(SpiralType eSpiralType)
{
    m_eSpiralType = eSpiralType;
}

inline SpiralWaveform::VDType SpiralWaveform::getVDType() const
{
    return m_eVDType;
}

inline void SpiralWaveform::setVDType(VDType eVDType)
{
    m_eVDType = eVDType;
}

inline double SpiralWaveform::getVDInnerCutoff() const
{
    return m_dVDInnerCutoff;
}

inline void SpiralWaveform::setVDInnerCutoff(double dVDInnerCutoff)
{
    m_dVDInnerCutoff = dVDInnerCutoff;
}

inline double SpiralWaveform::getVDOuterCutoff() const
{
    return m_dVDOuterCutoff;
}

inline void SpiralWaveform::setVDOuterCutoff(double dVDOuterCutoff)
{
    m_dVDOuterCutoff = dVDOuterCutoff;
}

inline double SpiralWaveform::getVDOuterDensity() const
{
    return m_dVDOuterDensity;
}

inline void SpiralWaveform::setVDOuterDensity(double dVDOuterDensity)
{
    m_dVDOuterDensity = dVDOuterDensity;
}

inline bool SpiralWaveform::getSloppy() const
{
    return m_bSloppy;
}

inline void SpiralWaveform::setSloppy(bool bSloppy)
{
    m_bSloppy = bSloppy;
}

inline double SpiralWaveform::getSloppyPeriod() const
{
    return m_dSloppyPeriod;
}

inline void SpiralWaveform::setSloppyPeriod(double dSloppyPeriod)
{
    m_dSloppyPeriod = dSloppyPeriod;
}

inline long SpiralWaveform::getRampUpTime() const
{
    return m_lRampUpTime;
}

inline void SpiralWaveform::setRampUpTime(long lRampUpTime)
{
    m_lRampUpTime = lRampUpTime;
}

inline long SpiralWaveform::getRampDownTime() const
{
    return m_lRampDownTime;
}

inline void SpiralWaveform::setRampDownTime(long lRampDownTime)
{
    m_lRampDownTime = lRampDownTime;
}

inline long SpiralWaveform::getReadOutTime() const
{
    return m_lReadOutTime;
}

inline long SpiralWaveform::getTotalTime() const
{
    return m_lTotalTime;
}

extern "C" int calculate_spiral_trajectory(
    float* pfTraj,
    long* plTrajSize,
    long lBaseResolution,
    long lSpiralArms,
    double dFieldOfView,
    double dMaxGradAmpl,
    double dMinRiseTime,
    double dDwellTime,
    double dReadoutOS = 2.0,
    double dGradientDelay = 0.0,
    double dLarmorConst = SWF_GAMMA_1H,
    double dVDInnerCutoff = 1.0,
    double dVDOuterCutoff = 1.0,
    double dVDOuterDensity = 1.0,
    int iVDType = 0);

#endif  // SPIRAL_WAVEFORM_H_

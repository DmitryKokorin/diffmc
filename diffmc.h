#pragma once

#include <string>
#include <vector>

#include "databuff.h"


class DiffMCApp
{
public:

    DiffMCApp();

    int  run();
    bool getOpts(int argc, char ** argv);
    void printHelp();

    static const int   kSeedIncrement = 1000; //random generator seeds in threads differ by this number


private:

    enum FileOptions { Create, Load, Save };

    struct Options
    {
        Options()
          : execName()
          , workDir()
          , eeProfileOptions(Create)
          , oeProfileOptions(Create)
          , eoProfileOptions(Create)
          , oFreePathOptions(Create)
          , eFreePathOptions(Create)
          , eeProbabilityOptions(Create)
          , eeProfileName()
          , oeProfileName()
          , eoProfileName()
          , oFreePathName()
          , eFreePathName()
          , eeProbabilityName()
          , seed(1000)
          , maxPhotons(1000)
          , maxScatterings(1000)
          , maxTime(100.)
          , points(500)
        {}

        std::string execName;
        std::string workDir;

        FileOptions eeProfileOptions;
        FileOptions oeProfileOptions;
        FileOptions eoProfileOptions;
        FileOptions oFreePathOptions;
        FileOptions eFreePathOptions;
        FileOptions eeProbabilityOptions;
        std::string eeProfileName;
        std::string oeProfileName;
        std::string eoProfileName;
        std::string oFreePathName;
        std::string eFreePathName;
        std::string eeProbabilityName;

        int     seed;
        int     maxPhotons;
        int     maxScatterings;
        Float   maxTime;
        int     points;
    };


    Options options_;



    void flushBuffers(const int scatteredCount, const DataBuff& buff);
    void output();


    bool  prepareOFreePath(LinearInterpolation& length);
    bool  prepareEFreePath(LinearInterpolation& length);

    bool  prepareEChannelProb(LinearInterpolation& prob);

    bool  prepareOEPartition(Partition& partition);
    bool  prepareEOPartition(Partition& partition);
    bool  prepareEEPartition(Partition& partition);

    size_t processScattering(const Photon& ph, DataBuff& buff, size_t idx);



    LinearInterpolation m_eLength;
    LinearInterpolation m_oLength;

    LinearInterpolation m_eChannelProb;

    int m_photonCnt;
    int m_saveRate;

    //right region border and number of iterations for some partition chunk
    typedef std::pair<Float, size_t>    ChunkParam;
    typedef std::vector<ChunkParam>     ChunkParamsList;

    ChunkParamsList    m_chunkParams;

    DataBuff m_dataBuff;

private:

    DiffMCApp& operator=(const DiffMCApp& other);
    DiffMCApp(const DiffMCApp& other);
};

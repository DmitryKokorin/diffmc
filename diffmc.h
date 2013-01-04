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
    size_t processScattering(const Photon& ph, DataBuff& buff, size_t idx);


    template <typename T>
    bool  prepareFreePath(LinearInterpolation &l, const char *name, const bool options, const std::string &fileName);

    bool  prepareEChannelProb(LinearInterpolation& prob);

    template <typename T>
    bool  preparePartition(Partition &p, const char *name, const bool options, const std::string &fileName);


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

template <typename T>
bool DiffMCApp::prepareFreePath(LinearInterpolation &l, const char *name, const bool options, const std::string &fileName)
{
    using namespace std;

    if (options == Load) {

        cerr << "loading " << name << "-beam free path file..." << endl;

        if (!l.load(fileName)) {

            cerr << "can't load " << name << "-beam free path data" << endl;
            return false;
        }
    }

    if (options == Create || options == Save) {

        cerr << "calculating " << name << "-beam free path data..." << endl;
        createFreePath<T>(l);
    }

    if (options == Save) {

        cerr << "saving " << name << "-beam free path data to file..." << endl;

        if (!l.save(fileName)) {

            cerr << "can't save " << name << "-beam free path data" << endl;
            return false;
        }
    }

    cerr << "done" << endl;

    return true;
}

template <typename T>
bool DiffMCApp::preparePartition(Partition &p, const char *name, const bool options, const std::string &fileName)
{
    using namespace std;

    if (options == Load) {

        cerr << "loading " << name << " profile file..." << endl;

        if (!p.load(fileName)) {

            cerr << "can't load " << name << " profile data" << endl;
            return false;
        }
    }

    if (options == Create || options == Save) {

        cerr << "calculating " << name << " profile data..." << endl;

        const int chunksNum = m_chunkParams.size();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < chunksNum; ++i) {

            p.addChunk<T>(i == 0 ? 0. : m_chunkParams[i-1].first,
                    m_chunkParams[i].first, m_chunkParams[i].second);
        }
    }

    if (options == Save) {

        cerr << "saving " << name << " profile data to file..." << endl;

        if (!p.save(fileName)) {

            cerr << "can't save " << name << " profile data" << endl;
            return false;
        }
    }

    cerr << "done" << endl;

    return true;
}

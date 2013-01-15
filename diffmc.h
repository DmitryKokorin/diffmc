#pragma once

#include <string>
#include <vector>

#include "databuff.h"
#include "interval.h"
#include "distance.h"


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

    void flushBuffers(int &scatteredCount, DataBuff &buff);
    void output();
    size_t processScattering(const Photon& ph, DataBuff& buff, size_t idx);


    template <typename T>
    bool  prepareFreePath(LinearInterpolation &l, const char *name, const int options, const std::string &fileName);

    bool  prepareEChannelProb(LinearInterpolation& prob);

    template <typename T>
    bool  preparePartition(Partition &p, const char *name, const int options, const std::string &fileName);


    LinearInterpolation m_eLength;
    LinearInterpolation m_oLength;

    LinearInterpolation m_eChannelProb;

    int m_photonCnt;
    int m_saveRate;

    //right region border and number of iterations for some partition chunk
    //typedef std::pair<Float, size_t>    ChunkParam;
    //typedef std::vector<ChunkParam>     ChunkParamsList;
    typedef std::vector<Interval>  IntervalsContainer;
    IntervalsContainer m_intervals;

    //ChunkParamsList    m_chunkParams;

    DataBuff m_dataBuff;

private:

    DiffMCApp& operator=(const DiffMCApp& other);
    DiffMCApp(const DiffMCApp& other);
};

template <typename T>
bool DiffMCApp::prepareFreePath(LinearInterpolation &l, const char *name, const int options, const std::string &fileName)
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
bool DiffMCApp::preparePartition(Partition &p, const char *name, const int options, const std::string &fileName)
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

        cerr << "calculating profile intervals..." << endl;

        Float left = 0.001;
        Float right = 0.5* M_PI;

        const Float tolerance = 0.01;

        Float dist;
        int i = 0;

        m_intervals.push_back(Interval(0., left));
        int num_threads = omp_get_max_threads();

        #pragma omp parallel for
        for (int j = 0; j < num_threads; ++j) {

            Float local_left = left + j/(Float)num_threads*(right-left);
            Float local_right = left + (j+1)/(Float)num_threads*(right-left);
            Float local_end = local_right;

            do {

                local_right = local_end;

                T ind_left = createIndicatrix<T>(local_left);

                while (1) {

                    T ind_right  = createIndicatrix<T>(local_right);

                    dist = ::distance(ind_left, ind_right);

                    if (dist > tolerance)
                        local_right = 0.5*(local_right + local_left);
                    else
                        break;
                }

                #pragma omp critical
                {
                    m_intervals.push_back(Interval(local_left, local_right));
                    std::cout << i << " " << local_left << " " << local_right << std::endl;

                    i++;
                }

                local_left = local_right;

            }
            while (local_right != local_end);
        }

        std::sort(m_intervals.begin(), m_intervals.end());

        cerr << "calculating " << name << " profile data..." << endl;

        const int chunksNum = m_intervals.size();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < chunksNum; ++i) {

              p.addChunk<T>(m_intervals[i].left, m_intervals[i].right);
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

#pragma once

#include <string>
#include <vector>

#include "vector3.h"


class DataBuff
{
public:

    DataBuff() :
        points()
    {}

    struct Point
    {
        Point()
          : r(), r2(), r3(), r4(), r5(), r6(), time(), scatterings(), measurements()
        {
        }
    
    	Vector3 r, r2, r3, r4, r5, r6;

        Float time;

        Float   scatterings;
        int     measurements;

        void clear()
        {
            r.clear();
            r2.clear();
            r3.clear();
            r4.clear();
            r5.clear();
            r6.clear();

            scatterings  = 0.;
            measurements = 0;
        }

        Point& operator+=(const Point& rhv)
        {
            r  += rhv.r;
            r2 += rhv.r2;
            r3 += rhv.r3;
            r4 += rhv.r4;
            r5 += rhv.r5;
            r6 += rhv.r6;

            scatterings  += rhv.scatterings;
            measurements += rhv.measurements;

            return *this;
        }
    };

    typedef std::vector<Point> Points;

    void clear()
    {
        Points::iterator i = points.begin();
        for (; i != points.end(); ++i) {

            i->clear();
        }
    }

    DataBuff& operator+=(const DataBuff& rhv)
    {
        Points::iterator i = points.begin();
        Points::const_iterator j = rhv.points.begin();

        for (; i != points.end(); ++i, ++j) {

            *i += *j;
        }

        return *this;
    }

    Points points;
};


class DiffMCApp
{
public:

    DiffMCApp();

    int  run();
    bool getOpts(int argc, char ** argv);
    void printHelp();

    static const int   kSeedIncrement = 1000; //random generator seeds in threads differ by this number


private:

    int  getSeed() const { return m_seed; }

    bool isLoadOFreePath() const {return m_loadOFreePath;}
    bool isSaveOFreePath() const {return m_saveOFreePath;}
    bool isLoadEFreePath() const {return m_loadEFreePath;}
    bool isSaveEFreePath() const {return m_saveEFreePath;}

    bool isLoadEChannelProb() const {return m_loadEChannelProb;}
    bool isSaveEChannelProb() const {return m_saveEChannelProb;}


    bool isLoadEEPartition() const {return m_loadEEPartition;}
    bool isSaveEEPartition() const {return m_saveEEPartition;}
    bool isLoadOEPartition() const {return m_loadOEPartition;}
    bool isSaveOEPartition() const {return m_saveOEPartition;}
    bool isLoadEOPartition() const {return m_loadEOPartition;}
    bool isSaveEOPartition() const {return m_saveEOPartition;}

    const std::string& getOFreePathFileName() const {return m_oFreePathFileName;}
    const std::string& getEFreePathFileName() const {return m_eFreePathFileName;}

    const std::string& getEChannelProbFileName() const {return m_eChannelProbFileName;}

    const std::string& getOEPartitionFileName() const {return m_oePartitionFileName;}
    const std::string& getEOPartitionFileName() const {return m_eoPartitionFileName;}
    const std::string& getEEPartitionFileName() const {return m_eePartitionFileName;}

    const std::string& getWorkDir() const {return m_workDir;}

    void flushBuffers(const int scatteredCount, const DataBuff& buff);


    bool  prepareOFreePath(LinearInterpolation& length);
    bool  prepareEFreePath(LinearInterpolation& length);

    bool  prepareEChannelProb(LinearInterpolation& prob);

    bool  prepareOEPartition(Partition& partition);
    bool  prepareEOPartition(Partition& partition);
    bool  prepareEEPartition(Partition& partition);

    void  prepareDataBuff(DataBuff& buff);
    size_t processScattering(const Photon& ph, DataBuff& buff, size_t idx);
    void  processLastScattering(const Photon& ph);

    void output();

    void outputPositions();

    std::string m_workDir;
    std::string m_execFileName;

    std::string m_oFreePathFileName;
    std::string m_eFreePathFileName;

    std::string m_eChannelProbFileName;

    std::string m_oePartitionFileName;
    std::string m_eoPartitionFileName;
    std::string m_eePartitionFileName;


    bool m_loadOFreePath;
    bool m_saveOFreePath;
    bool m_loadEFreePath;
    bool m_saveEFreePath;

    bool m_loadEChannelProb;
    bool m_saveEChannelProb;

    bool m_loadOEPartition;
    bool m_saveOEPartition;
    bool m_loadEOPartition;
    bool m_saveEOPartition;
    bool m_loadEEPartition;
    bool m_saveEEPartition;


    LinearInterpolation m_eLength;
    LinearInterpolation m_oLength;

    LinearInterpolation m_eChannelProb;

    int m_seed;
    int m_maxPhotons;
    int m_maxScatterings;

    Float   m_maxTime;
    int     m_points;

    int m_photonCnt;
    int m_saveRate;

    //right region border and number of iterations for some partition chunk
    typedef std::pair<Float, size_t>    ChunkParam;
    typedef std::vector<ChunkParam>     ChunkParamsList;

    ChunkParamsList    m_chunkParams;

    DataBuff m_dataBuff;

    std::vector<Vector3>    m_positions;

//    Float m_meanCos;

private:

    DiffMCApp& operator=(const DiffMCApp& other);
    DiffMCApp(const DiffMCApp& other);
};

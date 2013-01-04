#include "mathcompat.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <omp.h>


#include "common.h"

#include "channel.h"
#include "partition.h"
#include "photon.h"
#include "indicatrix.h"
#include "freepath.h"
#include "optics.h"

#include "diffmc.h"



int main(int argc, char ** argv)
{
    DiffMCApp app;

    int res = 0;

    if (app.getOpts(argc, argv)) {

        res = app.run();
    }
    else {

        app.printHelp();
        res = -1;
    }

    return res;
}


/////////////////////////////////////////////


DiffMCApp::DiffMCApp()
  : options_()
  , m_eLength()
  , m_oLength()
  , m_eChannelProb()
  , m_photonCnt(0)
  , m_saveRate(0)
  , m_chunkParams()
  , m_dataBuff()
{
}

bool DiffMCApp::getOpts(int argc, char ** argv)
{
    options_.execName = argv[0];

    for (int i = 1; i < argc; ++i) {

        std::string arg = argv[i];

        if (arg == "--workdir") {

            if (++i == argc)
                return false;

            options_.workDir = argv[i];

            if (!options_.workDir.empty())
                options_.workDir += '/';
        }
        else if (arg == "--seed") {

            if (++i == argc)
                return false;

            std::stringstream stream(argv[i]);
            stream >> options_.seed;
        }
        else if (arg == "--photons") {

            if (++i == argc)
                return false;

            std::stringstream stream(argv[i]);
            stream >> options_.maxPhotons;
        }
        else if (arg == "--scatterings") {

            if (++i == argc)
                return false;

            std::stringstream stream(argv[i]);
            stream >> options_.maxScatterings;
        }
        else if (arg == "--maxtime") {

            if (++i == argc)
                return false;

            std::stringstream stream(argv[i]);
            stream >> options_.maxTime;
        }
        else if (arg == "--points") {

            if (++i == argc)
                return false;

            std::stringstream stream(argv[i]);
            stream >> options_.seed;
        }
        else if (arg == "--loadoeprofile") {

            if (++i == argc)
                return false;

            options_.oeProfileOptions = Load;
            options_.oeProfileName    = argv[i];
        }
        else if (arg == "--loadeoprofile") {

            if (++i == argc)
                return false;

            options_.eoProfileOptions = Load;
            options_.eoProfileName    = argv[i];
        }
        else if (arg == "--loadeeprofile") {

            if (++i == argc)
                return false;

            options_.eeProfileOptions = Load;
            options_.eeProfileName    = argv[i];
        }
        else if (arg == "--loadofreepath") {

            if (++i == argc)
                return false;

            options_.oFreePathOptions = Load;
            options_.oFreePathName    = argv[i];
        }
        else if (arg == "--loadefreepath") {

            if (++i == argc)
                return false;

            options_.eFreePathOptions = Load;
            options_.eFreePathName    = argv[i];
        }
        else if (arg == "--loadeeprobability") {

            if (++i == argc)
                return false;

            options_.eeProbabilityOptions = Load;
            options_.eeProbabilityName    = argv[i];
        }
        else if (arg == "--saveoeprofile") {

            if (++i == argc)
                return false;

            options_.oeProfileOptions = Save;
            options_.oeProfileName    = argv[i];
        }
        else if (arg == "--saveeoprofile") {

            if (++i == argc)
                return false;

            options_.eoProfileOptions = Save;
            options_.eoProfileName    = argv[i];
        }
        else if (arg == "--saveeeprofile") {

            if (++i == argc)
                return false;

            options_.eeProfileOptions = Save;
            options_.eeProfileName    = argv[i];
        }
        else if (arg == "--saveofreepath") {

            if (++i == argc)
                return false;

            options_.oFreePathOptions = Save;
            options_.oFreePathName    = argv[i];
        }
        else if (arg == "--saveefreepath") {

            if (++i == argc)
                return false;

            options_.eFreePathOptions = Save;
            options_.eFreePathName    = argv[i];
        }
        else if (arg == "--saveeeprobability") {

            if (++i == argc)
                return false;

            options_.eeProbabilityOptions = Save;
            options_.eeProbabilityName    = argv[i];
        }
        else
            return false;
    }

    return true;
}




int DiffMCApp::run()
{
    using namespace std;

    cout.precision(17);
    cout << scientific;

    cerr.precision(17);
    cerr << scientific;

    cerr << "# seed = "    << options_.seed    << endl;
    cerr << "# maxtime = " << options_.maxTime << endl;
    cerr << "# H = "       << Optics::H        << endl;
    cerr << "# lambda = "  << Optics::lambda   << endl;

    options_.maxTime *= Optics::c;

    m_dataBuff.prepare(options_.points, options_.maxTime);

    //free path
    if (!prepareOFreePath(m_oLength))
        return -1;

    if (!prepareEFreePath(m_eLength))
        return -1;


    //e channel probability
    if (!prepareEChannelProb(m_eChannelProb))
        return -1;


    //partition
    int numChunks = 100;
    Float chunkStep = 0.5*M_PI / numChunks;

    for (int i = 1; i <= numChunks; ++i)
        m_chunkParams.push_back(ChunkParam(i*chunkStep, 20));

    Partition pOE, pEO, pEE;

    if (!prepareOEPartition(pOE))
        return -1;

    if (!prepareEOPartition(pEO))
        return -1;

    if (!prepareEEPartition(pEE))
        return -1;


    cerr << "scattering..." << endl;
    Photon::init(&m_oLength, &m_eLength, &pOE, &pEO, &pEE, &m_eChannelProb);

    const int flushRate = 1;
    m_saveRate  = omp_get_max_threads()*flushRate;

    const Float t = 0.5*M_PI; //angle with director
    const Vector3 initVector = Vector3(cos(t), 0, sin(t)).normalize();
    cerr << "initial angle: " << t << endl;


    //main loop

    #pragma omp parallel
    {
        RngEngine rng_engine;
        rng_engine.seed(options_.seed + kSeedIncrement*omp_get_thread_num());

        bool flush = false;
        int  scatteredCount = 0;
        DataBuff buff(options_.points, options_.maxTime);

        #pragma omp for schedule (dynamic)
        for (int i = 0; i < options_.maxPhotons; ++i) {

            Photon ph(rng_engine, initVector, Optics::ECHANNEL);
            size_t timeIdx = 0;

            if (flush) {

                flush = false;
                scatteredCount = 0;
                buff.clear();
            }


            while ((ph.scatterings < options_.maxScatterings) && (ph.time <= options_.maxTime)) {

                ph.move();

                timeIdx = processScattering(ph, buff, timeIdx);

                ph.scatter();

            }

            if (++scatteredCount == flushRate) {

                flush = true;
                flushBuffers(scatteredCount, buff);
            }
        }

        if (options_.maxPhotons % flushRate)
            flushBuffers(scatteredCount, buff);
    }

    m_dataBuff.average();

    output();
    return 0;
}


size_t DiffMCApp::processScattering(const Photon& ph, DataBuff& buff, size_t timeIdx)
{
    if (timeIdx >= buff.points.size())
        return timeIdx;

    if (ph.time < buff.points[timeIdx].time)
        return timeIdx;

    const Angle a = Angle(ph.s_i, Optics::director);
    const Float nn = (Optics::OCHANNEL == ph.channel) ? Optics::OBeam::n(a) : Optics::EBeam::n(a);

    DataBuff::Points::iterator i = buff.points.begin() + timeIdx;

    while ( i != buff.points.end() && (i->time < ph.time)) {

        Vector3 r = ph.pos - ph.s_i * (ph.time - i->time)/nn;

        i->appendPhoton(r, ph.scatterings);
        ++i;
    }

    return i - buff.points.begin();
}

void DiffMCApp::output()
{
    std::stringstream ss;
    ss << options_.workDir << "out.txt";

    std::fstream stream(ss.str().c_str(), std::fstream::out | std::fstream::trunc);
    stream << m_dataBuff;
    stream.close();
}

void DiffMCApp::printHelp()
{
    using namespace std;

    cerr << "Usage: " << options_.execName << "[options]" << endl;
    cerr << "Available options:" << endl;
    cerr << "--workdir [path]\t\t\toutput path" << endl;
    cerr << "--seed [seed]\t\t\tseed for random numbers generator" << endl;
    cerr << "--photons [photons]\t\t\tnumber of photons to scatter" << endl;
    cerr << "--scatterings [scatterings]\t\t\tmax scatterings for each photon" << endl;
    cerr << "--maxtime [maxtime]\t\t\tmax scattering time for each photon" << endl;
    cerr << "--points [points]\t\t\tnumber of sampling points" << endl;

    cerr << "--loadoeprofile [filename]\t\tload o-e profile from file"  << endl;
    cerr << "--loadeoprofile [filename]\t\tload e-o profile from file"  << endl;
    cerr << "--loadeeprofile [filename]\t\tload e-e profile from file"  << endl;
    cerr << "--loadoefreepath [filename]\t\tload o-beam free path from file"  << endl;
    cerr << "--loadeofreepath [filename]\t\tload e-beam free path from file"  << endl;
    cerr << "--loadeeprobability [filename]\t\tload e-e probability from file"  << endl;

    cerr << "--saveoeprofile [filename]\t\tsave o-e profile to file"    << endl;
    cerr << "--saveeoprofile [filename]\t\tsave e-o profile to file"    << endl;
    cerr << "--saveeeprofile [filename]\t\tsave e-e profile to file" << endl;
    cerr << "--saveoefreepath [filename]\t\tsave o-beam free path to file"  << endl;
    cerr << "--saveeofreepath [filename]\t\tsave e-beam free path to file"  << endl;
    cerr << "--saveeeprobability [filename]\t\tsave e-e probability to file"  << endl;
}

bool DiffMCApp::prepareOFreePath(LinearInterpolation& l)
{
    using namespace std;

    if (options_.oFreePathOptions == Load) {

        cerr << "loading o-beam free path file..." << endl;

        if (!l.load(options_.oFreePathName)) {

            cerr << "can't load o-beam free path data" << endl;
            return false;
        }
    }

    if (options_.oFreePathOptions == Create || options_.oFreePathOptions == Save) {

        cerr << "calculating o-beam free path data..." << endl;
        createFreePath<Optics::OBeam>(l);
    }

    if (options_.oFreePathOptions == Save) {

        cerr << "saving o-beam free path data to file..." << endl;

        if (!l.save(options_.oFreePathName)) {

            cerr << "can't save o-beam free path data" << endl;
            return false;
        }
    }

    cerr << "done" << endl;

    return true;
}


bool DiffMCApp::prepareEFreePath(LinearInterpolation& l)
{
    using namespace std;

    if (options_.eFreePathOptions == Load) {

        cerr << "loading e-beam free path file..." << endl;

        if (!l.load(options_.eFreePathName)) {

            cerr << "can't load e-beam free path data" << endl;
            return false;
        }
    }

    if (options_.eFreePathOptions == Create || options_.eFreePathOptions == Save) {

        cerr << "calculating e-beam free path data..." << endl;
        createFreePath<Optics::EBeam>(l);
    }

    if (options_.eFreePathOptions == Save) {

        cerr << "saving e-beam free path data to file..." << endl;

        if (!l.save(options_.eFreePathName)) {

            cerr << "can't save e-beam free path data" << endl;
            return false;
        }
    }

    cerr << "done" << endl;

    return true;
}


bool DiffMCApp::prepareEChannelProb(LinearInterpolation& l)
{
    using namespace std;

    if (options_.eeProbabilityOptions == Load) {

        cerr << "loading e-e probability file..." << endl;

        if (!l.load(options_.eeProbabilityName)) {

            cerr << "can't load e-e probability data" << endl;
            return false;
        }
    }

    if (options_.eeProbabilityOptions == Create || options_.eeProbabilityOptions == Save) {

        cerr << "calculating e-e probability data..." << endl;
        createEChannelProb<Optics::EBeam>(l);
    }

    if (options_.eeProbabilityOptions == Save) {

        cerr << "saving e-e probability data to file..." << endl;

        if (!l.save(options_.eeProbabilityName)) {

            cerr << "can't save e-e probability data" << endl;
            return false;
        }
    }

    cerr << "done" << endl;

    return true;
}


bool DiffMCApp::prepareEEPartition(Partition& p)
{
    using namespace std;

    if (options_.eeProfileOptions == Load) {

        cerr << "loading e-e profile file..." << endl;

        if (!p.load(options_.eeProfileName)) {

            cerr << "can't load e-e profile data" << endl;
            return false;
        }
    }

    if (options_.eeProfileOptions == Create || options_.eeProfileOptions == Save) {

        cerr << "calculating e-e profile data..." << endl;

        const int chunksNum = m_chunkParams.size();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < chunksNum; ++i) {

            p.addChunk<IndicatrixEE>(i == 0 ? 0. : m_chunkParams[i-1].first,
                    m_chunkParams[i].first, m_chunkParams[i].second);
        }
    }

    if (options_.eeProfileOptions == Save) {

        cerr << "saving e-e profile data to file..." << endl;

        if (!p.save(options_.eeProfileName)) {

            cerr << "can't save e-e profile data" << endl;
            return false;
        }
    }

    cerr << "done" << endl;

    return true;
}

bool DiffMCApp::prepareOEPartition(Partition& p)
{
    using namespace std;

    if (options_.oeProfileOptions == Load) {

        cerr << "loading o-e profile file..." << endl;

        if (!p.load(options_.oeProfileName)) {

            cerr << "can't load o-e profile data" << endl;
            return false;
        }
    }

    if (options_.oeProfileOptions == Create || options_.oeProfileOptions == Save) {

        cerr << "calculating o-e profile data..." << endl;

        const int chunksNum = m_chunkParams.size();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < chunksNum; ++i) {

            p.addChunk<IndicatrixOE>(i == 0 ? 0. : m_chunkParams[i-1].first,
                    m_chunkParams[i].first, m_chunkParams[i].second);
        }
    }

    if (options_.oeProfileOptions == Save) {

        cerr << "saving o-e profile data to file..." << endl;

        if (!p.save(options_.oeProfileName)) {

            cerr << "can't save o-e profile data" << endl;
            return false;
        }
    }

    cerr << "done" << endl;

    return true;
}


bool DiffMCApp::prepareEOPartition(Partition& p)
{
    using namespace std;

    if (options_.eoProfileOptions == Load) {

        cerr << "loading e-o profile file..." << endl;

        if (!p.load(options_.eoProfileName)) {

            cerr << "can't load e-o profile data" << endl;
            return false;
        }
    }

    if (options_.eoProfileOptions == Create || options_.eoProfileOptions == Save) {

        cerr << "calculating e-o profile data..." << endl;

        const int chunksNum = m_chunkParams.size();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < chunksNum; ++i) {

            p.addChunk<IndicatrixEO>(i == 0 ? 0. : m_chunkParams[i-1].first,
                    m_chunkParams[i].first, m_chunkParams[i].second);
        }
    }

    if (options_.eoProfileOptions == Save) {

        cerr << "saving e-o profile data to file..." << endl;

        if (!p.save(options_.eoProfileName)) {

            cerr << "can't save e-o profile data" << endl;
            return false;
        }
    }

    cerr << "done" << endl;

    return true;
}

void DiffMCApp::flushBuffers(const int scatteredCount, const DataBuff& buff)
{
    #pragma omp critical
    {
        m_dataBuff += buff;
        m_photonCnt += scatteredCount;

        fprintf(stderr, "Photons: %d\n", m_photonCnt);

        if (0 == m_photonCnt % m_saveRate)
            output();
    }
}

#include "mathcompat.h"

#include <cstdio>
#include <memory.h>
#include <string.h>
#include <sstream>
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


DiffMCApp::DiffMCApp() :
    m_workDir(),
    m_execFileName(),
    m_oFreePathFileName(),
    m_eFreePathFileName(),
    m_eChannelProbFileName(),
    m_oePartitionFileName(),
    m_eoPartitionFileName(),
    m_eePartitionFileName(),
    m_loadOFreePath(false),
    m_saveOFreePath(false),
    m_loadEFreePath(false),
    m_saveEFreePath(false),
    m_loadEChannelProb(false),
    m_saveEChannelProb(false),
    m_loadOEPartition(false),
    m_saveOEPartition(false),
    m_loadEOPartition(false),
    m_saveEOPartition(false),
    m_loadEEPartition(false),
    m_saveEEPartition(false),
    m_eLength(),
    m_oLength(),
    m_eChannelProb(),
    m_seed(1000),
    m_maxPhotons(1000),
    m_maxScatterings(1000),
    m_maxTime(100.),
    m_points(500),
    m_photonCnt(0),
    m_saveRate(0),
    m_chunkParams(),
    m_dataBuff(),
    m_positions()
    //, m_meanCos(0.)
{
    m_positions.reserve(m_maxPhotons);
}

bool DiffMCApp::getOpts(int argc, char ** argv)
{
    m_execFileName = argv[0];

    for (int i = 1; i < argc; ++i) {

        if (!strcmp(argv[i], "--seed")) {

            if (++i == argc)
                return false;

            m_seed = atoi(argv[i]);
        }
        else if(!strcmp(argv[i], "--workdir")) {

            if (++i == argc)
                return false;

            m_workDir = argv[i];

            if (!m_workDir.empty())
                m_workDir = m_workDir + '/';
            }
            else if(!strcmp(argv[i], "--loadofreepath")) {

            if (m_saveOFreePath)
                return false;

            if (++i == argc)
                return false;

            m_loadOFreePath     = true;
            m_oFreePathFileName = argv[i];
        }
        else if(!strcmp(argv[i], "--loadefreepath")) {

            if (m_saveEFreePath)
                return false;

            if (++i == argc)
                return false;

            m_loadEFreePath     = true;
            m_eFreePathFileName = argv[i];
        }
        else if(!strcmp(argv[i], "--loadechannelprob")) {

            if (m_saveEChannelProb)
                return false;

            if (++i == argc)
                return false;

            m_loadEChannelProb     = true;
            m_eChannelProbFileName = argv[i];
        }
        else if (!strcmp(argv[i], "--loadoepartition")) {

            if (m_saveOEPartition)
                return false;

            if (++i == argc)
                return false;

            m_loadOEPartition     = true;
            m_oePartitionFileName = argv[i];
        }
        else if (!strcmp(argv[i], "--loadeopartition")) {

            if (m_saveEOPartition)
                return false;

            if (++i == argc)
                return false;

            m_loadEOPartition     = true;
            m_eoPartitionFileName = argv[i];
        }
        else if (!strcmp(argv[i], "--loadeepartition")) {

            if (m_saveEEPartition)
                return false;

            if (++i == argc)
                return false;

            m_loadEEPartition     = true;
            m_eePartitionFileName = argv[i];
        }
        else if (!strcmp(argv[i], "--saveofreepath")) {

            if (m_loadOFreePath)
                return false;

            if (++i == argc)
                return false;

            m_saveOFreePath     = true;
            m_oFreePathFileName = argv[i];

        }
        else if (!strcmp(argv[i], "--saveefreepath")) {

            if (m_loadEFreePath)
                return false;

            if (++i == argc)
                return false;

            m_saveEFreePath     = true;
            m_eFreePathFileName = argv[i];

        }
        else if (!strcmp(argv[i], "--saveechannelprob")) {

            if (m_loadEChannelProb)
                return false;

            if (++i == argc)
                return false;

            m_saveEChannelProb     = true;
            m_eChannelProbFileName = argv[i];

        }
        else if (!strcmp(argv[i], "--saveoepartition")) {

            if (m_loadOEPartition)
                return false;

            if (++i == argc)
                return false;

            m_saveOEPartition     = true;
            m_oePartitionFileName = argv[i];
        }
        else if (!strcmp(argv[i], "--saveeopartition")) {

            if (m_loadEOPartition)
                return false;

            if (++i == argc)
                return false;

                m_saveEOPartition     = true;
                m_eoPartitionFileName = argv[i];
        }
        else if (!strcmp(argv[i], "--saveeepartition")) {

            if (m_loadEEPartition)
                return false;

            if (++i == argc)
                return false;

            m_saveEEPartition     = true;
            m_eePartitionFileName = argv[i];
        }
        else if (!strcmp(argv[i], "--photons")) {

            if (++i == argc)
                return false;

            m_maxPhotons = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--scatterings")) {

            if (++i == argc)
                return false;

            m_maxScatterings = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--maxtime")) {

            if (++i == argc)
                return false;

            m_maxTime = atof(argv[i]);
        }
        else if (!strcmp(argv[i], "--points")) {

            if (++i == argc)
                return false;

            m_points = atoi(argv[i]);
        }
/*        else if (!strcmp(argv[i], "--H")) {

            if (++i == argc)
                return false;

            Optics::H = atof(argv[i]);
        }*/
        else
            return false;
    }

    return true;
}




int DiffMCApp::run()
{
    fprintf(stderr, "# seed = %d\n", getSeed());
    fprintf(stderr, "# maxtime = %.17e\n", m_maxTime);
    fprintf(stderr, "# H = %.17e\n", Optics::H);
    fprintf(stderr, "# lambda = %.17e\n", Optics::lambda);

    //Optics::init();


    m_maxTime *= Optics::c;


    prepareDataBuff(m_dataBuff);


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


    fprintf(stderr, "scattering...\n");
    Photon::init(&m_oLength, &m_eLength, &pOE, &pEO, &pEE, &m_eChannelProb);

    const int flushRate = 1;
    m_saveRate  = omp_get_max_threads()*flushRate;

    const Float t = 0.5*M_PI; //angle with director
    const Vector3 initVector = Vector3(cos(t), 0, sin(t)).normalize();
    fprintf(stderr, "initial angle: %.17e\n", t);


    //main loop

    #pragma omp parallel
    {
        RngEngine rng_engine;
        rng_engine.seed(m_seed + kSeedIncrement*omp_get_thread_num());

        bool flush = false;
        int  scatteredCount = 0;
        DataBuff  buff;

        #pragma omp for schedule (dynamic)
        for (int i = 0; i < m_maxPhotons; ++i) {

            Photon ph(rng_engine, initVector, Optics::ECHANNEL);
            size_t timeIdx = 0;
            prepareDataBuff(buff);

            if (flush) {

                flush = false;
                scatteredCount = 0;
                buff.clear();
            }


            while ((ph.scatterings < m_maxScatterings) && (ph.time <= m_maxTime)) {

                ph.move();

                timeIdx = processScattering(ph, buff, timeIdx);

 //               Vector3 s_i = ph.s_i;
                ph.scatter();

 //               m_meanCos += s_i*ph.s_i;

            }

            processLastScattering(ph);

            if (++scatteredCount == flushRate) {

                flush = true;
                flushBuffers(scatteredCount, buff);
            }
        }

        if (m_maxPhotons % flushRate)
            flushBuffers(scatteredCount, buff);
    }

//    fprintf(stderr, "# <cos> = %f\n", m_meanCos / m_maxScatterings);

    output();

    outputPositions();

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

        Vector3 pos = ph.pos - ph.s_i * (ph.time - i->time)/nn;

        Float x2 = pos.x() * pos.x();
        Float y2 = pos.y() * pos.y();
        Float z2 = pos.z() * pos.z();

        i->x2 += x2;
        i->y2 += y2;
        i->z2 += z2;

        Float x4 = x2*x2;
        Float y4 = y2*y2;
        Float z4 = z2*z2;

        i->x4 += x4;
        i->y4 += y4;
        i->z4 += z4;

        i->x6 += x2*x4;
        i->y6 += y2*y4;
        i->z6 += z2*z4;

        i->scatterings += ph.scatterings;
        i->measurements++;

        //fprintf(stderr, "%.17e\t%.17e\n", ph.time, i->time);

        ++i;
    }

    return i - buff.points.begin();
}

void DiffMCApp::processLastScattering(const Photon& ph)
{
    if (ph.time < m_maxTime)
        return;

    const Angle a = Angle(ph.s_i, Optics::director);
    const Float nn = (Optics::OCHANNEL == ph.channel) ? Optics::OBeam::n(a) : Optics::EBeam::n(a);

    Vector3 pos = ph.pos - ph.s_i * (ph.time - m_maxTime)/nn;
    m_positions.push_back(pos);
}

void DiffMCApp::output()
{
    std::stringstream ss;
    ss << m_workDir << "out.txt";
    FILE* file = fopen(ss.str().c_str(), "w");


    DataBuff::Points::const_iterator i = m_dataBuff.points.begin();
    for (; i != m_dataBuff.points.end(); ++i) {

        if (i->measurements) {

            Float x2 = i->x2 / i->measurements;
            Float y2 = i->y2 / i->measurements;
            Float z2 = i->z2 / i->measurements;

            Float x4 = i->x4 / i->measurements;
            Float y4 = i->y4 / i->measurements;
            Float z4 = i->z4 / i->measurements;

            Float x6 = i->x6 / i->measurements;
            Float y6 = i->y6 / i->measurements;
            Float z6 = i->z6 / i->measurements;

            fprintf(file, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%d\n",
                    i->time / Optics::c,
                    x2, y2, z2,
                    x4, y4, z4,
                    x6, y6, z6,
                    ((Float)i->scatterings) / i->measurements,
                    i->measurements);
        }
    }


    fclose(file);
}

void DiffMCApp::outputPositions()
{
    std::stringstream ss;
    ss << m_workDir << "pos.txt";
    FILE* file = fopen(ss.str().c_str(), "w");

    std::vector<Vector3>::iterator i = m_positions.begin();
    for (; i != m_positions.end(); ++i) {

        fprintf(file, "%.17e\t%.17e\t%.17e\n", i->x(), i->y(), i->z());
    }

    fclose(file);
}


void DiffMCApp::printHelp()
{
    fprintf(stderr, "Usage: %s [options]", m_execFileName.c_str());
    fprintf(stderr, "\n\nAvailable options:");
    fprintf(stderr, "\n--seed [seed]\t\t\t\tseed for random numbers generator");
    fprintf(stderr, "\n--workdir [path]\t\t\toutput path");
    fprintf(stderr, "\n--loadofreepath [filename]\t\tload o-beam free path from file");
    fprintf(stderr, "\n--loadefreepath [filename]\t\tload e-beam free path from file");
    fprintf(stderr, "\n--loadechannelprob [filename]\t\tload e-e probability from file");
    fprintf(stderr, "\n--loadoepartition [filename]\t\tload o-e partition from file");
    fprintf(stderr, "\n--loadeopartition [filename]\t\tload e-o partition from file");
    fprintf(stderr, "\n--loadeepartition [filename]\t\tload e-e partition from file");
    fprintf(stderr, "\n--saveofreepath [filename]\t\tsave o-beam free path to file");
    fprintf(stderr, "\n--saveefreepath [filename]\t\tsave e-beam free path to file");
    fprintf(stderr, "\n--saveechannelprob [filename]\t\tsave e-e probability to file");
    fprintf(stderr, "\n--saveoepartition [filename]\t\tsave o-e partition to file");
    fprintf(stderr, "\n--saveeopartition [filename]\t\tsave e-o partition to file");
    fprintf(stderr, "\n--saveeepartition [filename]\t\tsave e-e partition to file");
    fprintf(stderr, "\n--photons [photons]\t\t\tnumber of photons to scatter");
    fprintf(stderr, "\n--scatterings [scatterings]\t\tmax scatterings for each photon");
    fprintf(stderr, "\n--maxtime [time]\t\tmax scattering time for each photon");
    fprintf(stderr, "\n--points [points]\t\tnumber of sampling points");
    fprintf(stderr, "\n");
}

bool DiffMCApp::prepareOFreePath(LinearInterpolation& l)
{
    if (isLoadOFreePath()) {

        fprintf(stderr, "loading o-beam free path file...");

        if (!l.load(getOFreePathFileName())) {

            fprintf(stderr, "can't load o-beam free path data\n");
            return false;
        }
        else {

            fprintf(stderr, "\tdone\n");
        }
    }
    else {

        fprintf(stderr, "calculating o-beam free path data...\n");
        createFreePath<Optics::OBeam>(l);
    }

    if (isSaveOFreePath()) {

        fprintf(stderr, "saving o-beam free path data to file...");

        if (!l.save(getOFreePathFileName())) {

            fprintf(stderr, "can't save o-beam free path data\n");
            return false;
        }
        else {

            fprintf(stderr, "\tdone\n");
        }
    }

    return true;
}


bool DiffMCApp::prepareEFreePath(LinearInterpolation& l)
{
    if (isLoadEFreePath()) {

        fprintf(stderr, "loading e-beam free path file...");

        if (!l.load(getEFreePathFileName())) {

            fprintf(stderr, "can't load e-beam free path data\n");
            return false;
        }
        else {

            fprintf(stderr, "\tdone\n");
        }
    }
    else {

        fprintf(stderr, "calculating e-beam free path data...\n");
        createFreePath<Optics::EBeam>(l);
    }

    if (isSaveEFreePath()) {

        fprintf(stderr, "saving e-beam free path data to file...");

        if (!l.save(getEFreePathFileName())) {

            fprintf(stderr, "can't save e-beam free path data\n");
            return false;
        }
        else {

            fprintf(stderr, "\tdone\n");
        }
    }

    return true;
}



bool DiffMCApp::prepareEChannelProb(LinearInterpolation& l)
{
    if (isLoadEChannelProb()) {

        fprintf(stderr, "loading e-e probability file...");

        if (!l.load(getEChannelProbFileName())) {

            fprintf(stderr, "can't load e-e probability data\n");
                return false;
        }
        else {

            fprintf(stderr, "\tdone\n");
        }
    }
    else {

        fprintf(stderr, "calculating e-e probability data...\n");
        createEChannelProb<Optics::EBeam>(l);
    }

    if (isSaveEChannelProb()) {

        fprintf(stderr, "saving e-e probability data to file...");

        if (!l.save(getEChannelProbFileName())) {

            fprintf(stderr, "can't save e-e probability data\n");
            return false;
        }
        else {

            fprintf(stderr, "\tdone\n");
        }
    }

    return true;
}


bool DiffMCApp::prepareEEPartition(Partition& p)
{
    if (isLoadEEPartition()) {

        fprintf(stderr, "loading e-e partition...");

        if (!p.load(getEEPartitionFileName())) {

            fprintf(stderr, "can't load e-e partition\n");
            return false;
        }
        else {

            fprintf(stderr, "\tdone\n");
        }
    }
    else {

        fprintf(stderr, "creating e-e partition...\n");

        const int chunksNum = m_chunkParams.size();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < chunksNum; ++i) {

            p.addChunk<IndicatrixEE>(i == 0 ? 0. : m_chunkParams[i-1].first,
                    m_chunkParams[i].first, m_chunkParams[i].second);
        }

    }

    if (isSaveEEPartition()) {

        fprintf(stderr, "saving e-e partition...");

        if (!p.save(getEEPartitionFileName())) {

            fprintf(stderr, "can't save e-e partition\n");
            return false;
        }
        else {

            fprintf(stderr, "\tdone\n");
        }
    }

    return true;
}

bool DiffMCApp::prepareOEPartition(Partition& p)
{
    if (isLoadOEPartition()) {

        fprintf(stderr, "loading o-e partition...");

        if (!p.load(getOEPartitionFileName())) {

            fprintf(stderr, "can't load o-e partition\n");
            return false;
        }
        else {

            fprintf(stderr, "\tdone\n");
        }
    }
    else {

        fprintf(stderr, "creating o-e partition...\n");

        const int chunksNum = m_chunkParams.size();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < chunksNum; ++i) {

            p.addChunk<IndicatrixOE>(i == 0 ? 0. : m_chunkParams[i-1].first,
                    m_chunkParams[i].first, m_chunkParams[i].second);
        }

    }

    if (isSaveOEPartition()) {

        fprintf(stderr, "saving o-e partition...");

        if (!p.save(getOEPartitionFileName())) {

            fprintf(stderr, "can't save o-e partition\n");
            return false;
        }
        else {

            fprintf(stderr, "\tdone\n");
        }
    }

    return true;
}


bool DiffMCApp::prepareEOPartition(Partition& p)
{
    if (isLoadEOPartition()) {

        fprintf(stderr, "loading e-o partition...");

        if (!p.load(getEOPartitionFileName())) {

            fprintf(stderr, "can't load e-o partition\n");
            return false;
        }
        else {

            fprintf(stderr, "\tdone\n");
        }
    }
    else {

        fprintf(stderr, "creating e-o partition...\n");

        const int chunksNum = m_chunkParams.size();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < chunksNum; ++i) {

            p.addChunk<IndicatrixEO>(i == 0 ? 0. : m_chunkParams[i-1].first,
                    m_chunkParams[i].first, m_chunkParams[i].second);
        }

    }

    if (isSaveEOPartition()) {

        fprintf(stderr, "saving e-o partition...");

        if (!p.save(getEOPartitionFileName())) {

            fprintf(stderr, "can't save e-o partition\n");
                return false;
        }
        else {

            fprintf(stderr, "\tdone\n");
        }
    }

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

void DiffMCApp::prepareDataBuff(DataBuff& buff)
{
    for (int j = 0; j < m_points; ++j) {

        DataBuff::Point pt;
        pt.measurements = 0;
        pt.x2 = 0.;
        pt.y2 = 0.;
        pt.z2 = 0.;
        pt.x4 = 0.;
        pt.y4 = 0.;
        pt.z4 = 0.;
        pt.x6 = 0.;
        pt.y6 = 0.;
        pt.z6 = 0.;
        pt.time = j * (m_maxTime / m_points);
        buff.points.push_back(pt);
    }
}

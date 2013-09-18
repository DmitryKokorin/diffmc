#ifndef _PHOTON_H_
#define _PHOTON_H_

#ifdef WIN32
#include <random>
#else
#include <tr1/random>
#endif

#include <limits.h>

#include "vector3.h"
#include "matrix3.h"
#include <stdint.h>


class LinearInterpolation;
class Partition;
class PartitionChunk;
class Angle;


typedef std::tr1::mt19937 RngEngine;
typedef std::tr1::uniform_real<Float> RngDist;
typedef std::tr1::variate_generator <RngEngine, RngDist > Rng;


class Photon

{
public:

    Photon(RngEngine& rng_engine_, const Vector3& s = Vector3(0., 0., 1.), const int channel_ = Optics::ECHANNEL);

    static void init(   LinearInterpolation*    oLength,
                        LinearInterpolation*    eLength,
                        Partition*              oePartition,
                        Partition*              eoPartition,
                        Partition*              eePartition,
                        LinearInterpolation*    eChannelProb);

    void move();
    void scatter();

    Vector3   pos;  //position
    Vector3   s_i;  //normalized wave vector

    int64_t  scatterings;
    int      channel;
    Float    time;

private:

    void createTransformToPartitionCoords(Matrix3& mtx, Vector3& nn, Angle& a_i);

    void choosePointInRect(Float& x, Float& y, const int rectIdx, const Float randX, const Float randY);

    RngEngine& rng_engine;

    inline Float random() { return Float(Photon::rng_engine()) / (uint64_t(UINT_MAX) + 1); }

    static LinearInterpolation* s_oLength;
    static LinearInterpolation* s_eLength;

    static Partition*           s_oePartition;
    static Partition*           s_eoPartition;
    static Partition*           s_eePartition;

    static LinearInterpolation* s_eChannelProb;

    //these are to simulate static behavior for a reference (without ugly pointer syntax)
    LinearInterpolation&    oLength;
    LinearInterpolation&    eLength;

    Partition&              oePartition;
    Partition&              eoPartition;
    Partition&              eePartition;

    LinearInterpolation&    eChannelProb;

    PartitionChunk *m_chunk; //current chunk

    //disable copying
    Photon(const Photon&);
    Photon& operator=(const Photon&);
};


#endif /* _PHOTON_H_ */

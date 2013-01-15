#include "mathcompat.h"

#include <functional>
#include <iostream>
#include <omp.h>

#include "partition.h"
#include "partchunk.h"
#include "optics.h"
#include "indicatrix.h"
#include "freepath.h"
#include "coords.h"
#include "photon.h"


LinearInterpolation*    Photon::s_oLength       = NULL;
LinearInterpolation*    Photon::s_eLength       = NULL;

Partition*              Photon::s_oePartition   = NULL;
Partition*              Photon::s_eoPartition   = NULL;
Partition*              Photon::s_eePartition   = NULL;

LinearInterpolation*    Photon::s_eChannelProb  = NULL;



void Photon::init(  LinearInterpolation*    oLength_,
                    LinearInterpolation*    eLength_,
                    Partition*              oePartition_,
                    Partition*              eoPartition_,
                    Partition*              eePartition_,
                    LinearInterpolation*    eChannelProb_)
{
    s_oLength     = oLength_;
    s_eLength     = eLength_;

    s_oePartition = oePartition_;
    s_eoPartition = eoPartition_;
    s_eePartition = eePartition_;

    s_eChannelProb = eChannelProb_;
}


Photon::Photon(RngEngine& rng_engine_, const Vector3& s, const int channel_) :
    pos(0.,0.,0.),
    s_i(s),
    scatterings(0),
    channel(channel_),
    time(0),
    rng_engine(rng_engine_),
    oLength(*s_oLength),
    eLength(*s_eLength),
    oePartition(*s_oePartition),
    eoPartition(*s_eoPartition),
    eePartition(*s_eePartition),
    eChannelProb(*s_eChannelProb),
    m_chunk(NULL)
{
}


void Photon::move()
{
    Float rnd;
    const Angle a = Angle(s_i, Optics::director);
    const Float theta = symmetrizeTheta(a.theta);

    const Float meanFreePath = (Optics::OCHANNEL == channel) ? oLength(theta) : eLength(theta);
    const Float nn = (Optics::OCHANNEL == channel) ? Optics::OBeam::n(a) : Optics::EBeam::n(a);

    rnd = random();

    const Float d = -log1p(-rnd)*meanFreePath;

    pos  += d*s_i;
    time += d*nn;
}

void Photon::scatter()
{
    //prepare random numbers
    Float randChannel = random();
    Float randRect    = random();
    Float randX       = random();
    Float randY       = random();
    Float randPhi     = random();

    Vector3 nn;
    Matrix3 mtx;
    Angle   a_i;

    createTransformToPartitionCoords(mtx, nn, a_i);


    //select channel
    int newChannel;

    if (Optics::OCHANNEL == channel) {

        newChannel = Optics::ECHANNEL;
        m_chunk = eoPartition.getChunk(a_i.theta);
    }
    else {

        newChannel = randChannel > eChannelProb(symmetrizeTheta(Angle(s_i, Optics::director).theta)) ?
            Optics::OCHANNEL : Optics::ECHANNEL;

        m_chunk = Optics::OCHANNEL == newChannel ? eoPartition.getChunk(a_i.theta) : eePartition.getChunk(a_i.theta);
    }

    typedef PartitionChunk::ValuesVector ValuesVector;
    const ValuesVector &values = m_chunk->values();
    size_t rectIdx = 0;
    ValuesVector::const_iterator i = std::lower_bound(values.begin(), values.end(), randRect);

    if (values.end() == i) {

        std::cerr << "rect not found!" << randRect << ' ' << values.front() << ' '<< values.back() << std::endl;
    }
    else
        rectIdx = i - values.begin();

    //adjust point
    Float p, t;
    choosePointInRect(t, p, rectIdx, randX, randY);

    if (randPhi > 0.5)   //choose one of symmetrical cases
        p = 2*M_PI - p;

    Float sintheta = sin(t);
    Vector3 s_s =  Vector3(sintheta*cos(p), sintheta*sin(p), cos(t));

    s_i = invert(mtx)*s_s;
    s_i.normalize();  //to be sure

    channel = newChannel;

    scatterings++;
}

void Photon::createTransformToPartitionCoords(Matrix3& mtx, Vector3& nn, Angle& a_i)
{
    Vector3 v2;
    nn = Optics::director;

    if (Angle(s_i, nn).costheta < 0)
        nn = -nn;

    a_i = Angle(s_i, nn);

    if (fabs(a_i.sintheta) > kMachineEpsilon)
        v2 = crossProduct(s_i, nn).normalize();
    else
        v2 = createSomePerpendicular(s_i).normalize();

    Vector3 v3 = crossProduct(s_i, v2).normalize();

    mtx = createTransformMatrix(v2, v3, s_i);

    nn = mtx*nn;                            //director in s_i -based coordinate system
}

void Photon::choosePointInRect(Float& x, Float& y, const int rectNum, const Float randX, const Float randY)
{
    Rect& rect = m_chunk->rects()[rectNum];
    KnotsVector &knots = m_chunk->knots();

    Float b1 = knots[rect.tl].value;
    Float b2 = knots[rect.tr].value - b1;
    Float b3 = knots[rect.bl].value - b1;
    Float b4 = b1 - knots[rect.tr].value - knots[rect.bl].value + knots[rect.br].value;

    int roots;

    {
        Float x1, x2;
        roots = solveQuadric(b2 + 0.5*b4, 0.5*(b1 + 0.5*b3), -randX*(b2+0.5*b4+0.5*b1+0.25*b3), x1, x2);

        if (roots == 1) {

            x = x1;
        }
        else if (roots == 2) {

            if (x1 >= 0. && x1 < 1.) {

                x = x1;
            }
            else if (x2 >= 0. && x2 < 1.) {

                x = x2;
            }
            else {

                x = 0.5;
                std::cerr << "x out of range, " << x1 << '\t' << x2 << std::endl;
            }
        }
    }

    {
        Float y1, y2;
        roots = solveQuadric(0.5*(b3 + b4*x), b1 + b2*x, -randY*(b1+b2*x + 0.5*(b3+b4*x)), y1, y2);

        if (roots == 1) {

            y = y1;
        }
        else if (roots == 2) {

            if (y1 >= 0. && y1 < 1.) {

                y = y1;
            }
            else if (y2 >= 0. && y2 < 1.) {

                y = y2;
            }
            else {

                y = 0.5;
                std::cerr << "y out of range, " << y1 << '\t' << y2 << std::endl;
            }
        }
    }

    x *= rect.width;
    y *= rect.height;

    x += m_chunk->knots()[rect.tl].x;
    y += m_chunk->knots()[rect.tl].y;
}

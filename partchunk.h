#pragma once

#include <vector>
#include <string>
#include <map>

#include "common.h"
#include "node.h"
#include "rect.h"
#include "optics.h"
#include "indicatrix.h"
#include "matrix3.h"
#include "coords.h"


typedef std::vector<Rect> RectsVector;
typedef std::vector<Knot> KnotsVector;


class PartitionChunk
{
public:

    PartitionChunk();
    virtual ~PartitionChunk();

    template <class T>
    bool create(const Float minAngle, const Float maxAngle);

    bool load(FILE *file);
    bool save(FILE *file);

    void setData(Float** const data, const Float& cellSquare);
    void refine();

    inline size_t getRectsCount() const { return m_rects.size(); }
    inline size_t getKnotsCount() const { return m_knots.size(); }
    inline bool  isAngleInRange(const Float angle) const { return angle > m_minAngle && angle <= m_maxAngle;}

    static const Float kEpsilon;
    static const int   kThetaDegree = 15;
    static const int   kThetaSize   = (1 << kThetaDegree) + 1;
    static const int   kPhiDegree   = 7;
    static const int   kPhiSize     = (1 << kPhiDegree) + 1;

    static const Float kThetaResolution;
    static const Float kPhiResolution;


    RectsVector m_rects;
    KnotsVector m_knots;


private:

    Float integral(const GreedRect& rect);
    Float approxIntegral(const GreedRect& rect);
    Float rectError(const GreedRect& rect);

    template <class T>
    void createPartitionTree();

    void createRectsList();

    void cleanUp();
    void normalize();

    Float m_minAngle;
    Float m_maxAngle;

    Node*    m_root;
    Float**  m_data;             //knots
    Float**  m_cellIntegrals;    //integrals of elementary cells

    std::map<int, int> m_knotsMap;

    //used on tree creation
    void refineNode(Node* node);

    //used on list creation
    void processTreeNode(Node* node);

    Float m_cellSquare;
    Float m_fullIntegral;

    //disable copying
    PartitionChunk(const PartitionChunk&);
    PartitionChunk& operator=(const PartitionChunk&);
};

template <class T>
bool PartitionChunk::create(const Float minAngle, const Float maxAngle)
{
    m_minAngle   = minAngle;
    m_maxAngle   = maxAngle;

    m_root = new Node(GreedRect(0,0, kThetaSize-1, kPhiSize-1));

    m_cellIntegrals = allocate2dArray<Float>(kPhiSize-1, kThetaSize-1);

    createPartitionTree<T>();
    createRectsList();

    cleanUp();

    return true;
}

template <class T>
void PartitionChunk::createPartitionTree()
{
    Float **data = allocate2dArray<Float>(kPhiSize, kThetaSize);

    Float theta_i = 0.5*(m_minAngle + m_maxAngle);
    T ind = createIndicatrix<T>(theta_i);

    //calculate array values

    Float t, p;

    for (int i = 0; i < kThetaSize; ++i) {

        t = i*kThetaResolution;

        for (int j = 0; j < kPhiSize; ++j) {

            p = j*kPhiResolution;

            Vector3 ss_s = Vector3(sin(t)*cos(p), sin(t)*sin(p), cos(t));
            data[j][i]  = ind(ss_s)*sin(t);
        }
    }

    setData(data, kThetaResolution * kPhiResolution);
    refine();
    normalize();

    //free2dArray(data);
}

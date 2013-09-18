#include "mathcompat.h"

#include <cstdio>
#include <iostream>

#include "mathcompat.h"
#include "common.h"

#include "partchunk.h"



////////////////////////////////////////////////////////////////////////////////
//
// theta - x
// phi   - y
//

//const Float PartitionChunk::kEpsilon = 0.01;
//const Float PartitionChunk::kEpsilon = 0.01;
const Float PartitionChunk::kEpsilon = 0.01;

const Float PartitionChunk::kMinX = -1.;
const Float PartitionChunk::kMaxX = 1.;
const Float PartitionChunk::kMinY = 0.;
const Float PartitionChunk::kMaxY = M_PI;

const Float PartitionChunk::kIntegralEpsilon = 1e-10;

PartitionChunk::PartitionChunk() :
    m_rects(),
    m_knots(),
    m_values(),
    m_minAngle(0.),
    m_maxAngle(0.),
    m_root(NULL),
    m_knotsMap(),
    m_infos(),
    m_fullIntegral(0.),
    m_rectsCount(0)

{
}


PartitionChunk::~PartitionChunk()
{
    cleanUp();
}

void PartitionChunk::cleanUp()
{
    delete m_root; m_root = NULL;
    m_knotsMap.clear();

    m_infos.clear();
}

void PartitionChunk::processNode(Node* node)
{
    if (node->isLeaf()) {

        int  indeces[4];
        Knot knots[4];

        Rect &rect = node->interpolation.rect;

        knots[0] = Knot(rect.x1, rect.y1, node->interpolation.f_x1y1);
        knots[1] = Knot(rect.x2, rect.y1, node->interpolation.f_x2y1);
        knots[2] = Knot(rect.x1, rect.y2, node->interpolation.f_x1y2);
        knots[3] = Knot(rect.x2, rect.y2, node->interpolation.f_x2y2);

        for (int i = 0; i < 4; ++i) {

            if (m_knotsMap.find(knots[i]) == m_knotsMap.end()) {

                indeces[i] = m_knots.size();
                m_knots.push_back(knots[i]);
                m_knotsMap[knots[i]] = indeces[i];
            }
            else {

                indeces[i] = m_knotsMap[knots[i]];
            }
        }

        m_rects.push_back(PartitionRect(indeces[0], indeces[1], indeces[2], indeces[3], &m_knots));
    }
    else {

        processNode(node->pChild1);
        processNode(node->pChild2);
    }
}


void PartitionChunk::process()
{
    InfosContainer::const_iterator i = m_infos.begin();

    for (; i != m_infos.end(); ++i) {

        int  indeces[4];
        Knot knots[4];

        const Rect &rect = i->interpolation.rect;

        knots[0] = Knot(rect.x1, rect.y1, i->interpolation.f_x1y1);
        knots[1] = Knot(rect.x2, rect.y1, i->interpolation.f_x2y1);
        knots[2] = Knot(rect.x1, rect.y2, i->interpolation.f_x1y2);
        knots[3] = Knot(rect.x2, rect.y2, i->interpolation.f_x2y2);

        for (int j = 0; j < 4; ++j) {

            if (m_knotsMap.find(knots[j]) == m_knotsMap.end()) {

                indeces[j] = m_knots.size();
                m_knots.push_back(knots[j]);
                m_knotsMap[knots[j]] = indeces[j];
            }
            else {

                indeces[j] = m_knotsMap[knots[j]];
            }
        }

        m_rects.push_back(PartitionRect(indeces[0], indeces[1], indeces[2], indeces[3], &m_knots));
    }
}

Float PartitionChunk::currentError() const
{
    Float error = 0.;

    InfosContainer::const_iterator i = m_infos.begin();
    for (; i != m_infos.end(); ++i)
        error += i->error;

    return error;
}



bool PartitionChunk::load(FILE *file)
{
    m_knots.clear();
    m_rects.clear();

    if (   !fscanf(file, "%le", &m_minAngle)
        || !fscanf(file, "%le", &m_maxAngle))
        return false;

    int knotsNumber;
    if (!fscanf(file, "%d", &knotsNumber))
        return false;

    m_knots.reserve(knotsNumber);

    for (int i = 0; i < knotsNumber; ++i) {

        Float x, y, value;
        if (!fscanf(file, "%le\t%le\t%le", &x, &y, &value))
            return false;

        m_knots.push_back(Knot(x, y, value));
    }

    unsigned long int rectsNumber, i;

    if (!fscanf(file, "%lu", &rectsNumber))
        return false;

    m_rects.reserve(rectsNumber);

    for (i = 0; i < rectsNumber; ++i) {

        int tl, tr, bl, br;
        if (!fscanf(file, "%d\t%d\t%d\t%d", &tl, &tr, &bl, &br))
            return false;

        m_rects.push_back(PartitionRect(tl, tr, bl, br, &m_knots));
    }

    fillValues();

    fprintf(stderr, "chunk loaded: %e -- %e\t%lu\n", m_minAngle, m_maxAngle, (unsigned long int)m_rects.size());
    return true;
}

bool PartitionChunk::save(FILE *file)
{
    //angles
    fprintf(file, "%.17e\n", m_minAngle);
    fprintf(file, "%.17e\n", m_maxAngle);

    //knots
    {
        KnotsVector::iterator i;

        fprintf(file, "%lu\n", (long unsigned int)m_knots.size());

        for (i = m_knots.begin(); i != m_knots.end(); ++i) {

            fprintf(file, "%.17e\t%.17e\t%.17e\n", i->x, i->y, i->value);
        }
    }

    //rects
    {
        PartitionRectsVector::iterator i;

        fprintf(file, "%lu\n", (long unsigned int)m_rects.size());

        for (i = m_rects.begin(); i != m_rects.end(); ++i) {

            fprintf(file, "%d\t%d\t%d\t%d\n", i->tl, i->tr, i->bl, i->br);
        }
    }

    return true;
}

void PartitionChunk::normalize()
{
    //calculate full integral
    Float integral = 0.;

    PartitionRectsVector::iterator i = m_rects.begin();
    for (; i != m_rects.end(); ++i) {

        Float rectIntegral = 0.25* (m_knots[i->tl].value +
                                    m_knots[i->tr].value +
                                    m_knots[i->bl].value +
                                    m_knots[i->br].value) *
                            i->square;

        integral += rectIntegral;
    }

    KnotsVector::iterator j = m_knots.begin();
    for (; j != m_knots.end(); ++j) {

        j->value /= integral;
    }
}

void PartitionChunk::fillValues()
{
    PartitionRectsVector::iterator i;
    Float value = 0.;

    m_values.clear();
    m_values.reserve(m_rects.size());

    for (i = m_rects.begin(); i != m_rects.end(); ++i) {

        Float rectIntegral = 0.25* (m_knots[i->tl].value +
                                    m_knots[i->tr].value +
                                    m_knots[i->bl].value +
                                    m_knots[i->br].value) *
                            i->square;

        value += rectIntegral;

        m_values.push_back(value);
    }

#if defined TEST
    std::cerr << "full integral: "<<  value << std::endl;
#endif
}

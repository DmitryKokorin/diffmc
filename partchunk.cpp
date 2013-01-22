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

const Float PartitionChunk::kEpsilon = 0.01;
const Float PartitionChunk::kThetaResolution = M_PI / PartitionChunk::kThetaSize;
const Float PartitionChunk::kPhiResolution   = M_PI / PartitionChunk::kPhiSize;


PartitionChunk::PartitionChunk() :
    m_rects(),
    m_knots(),
    m_values(),
    m_minAngle(0.),
    m_maxAngle(0.),
    m_root(NULL),
    m_data(NULL),
    m_cellIntegrals(NULL),
    m_knotsMap(),
    m_cellSquare(0.),
    m_fullIntegral(0.)
{
}


PartitionChunk::~PartitionChunk()
{
    cleanUp();
}

void PartitionChunk::cleanUp()
{
    delete m_root; m_root = NULL;

    if (m_cellIntegrals) {

        free2dArray(m_cellIntegrals);
        m_cellIntegrals = NULL;
    }

    if (m_data) {

        free2dArray(m_data);
        m_data = NULL;
    }
}

void PartitionChunk::setData(Float** const data, const Float& cellSquare)
{
    m_data = data;
    m_cellSquare = cellSquare;
    m_fullIntegral = 0.;

    for (int i = 0; i < kThetaSize-1; ++i)
        for (int j = 0; j < kPhiSize-1; ++j) {

            m_cellIntegrals[j][i] = approxIntegral(GreedRect(i, j, i+1, j+1));
            m_fullIntegral += m_cellIntegrals[j][i];
        }
}

void PartitionChunk::refine()
{
    refineNode(m_root);
}

void PartitionChunk::refineNode(Node* node)
{
    if (node->isLeaf()) {

        Float nodeIntegral    = integral(node->rect);
        Float rectMaxError    = nodeIntegral*kEpsilon;

        //std::cerr << "nodeIntegral " << nodeIntegral << std::endl;

        if (node->rect.canSplitX() && node->rect.canSplitY()) {

            Float xSplitError = rectError(node->rect.leftHalf()) +
                                rectError(node->rect.rightHalf());

            Float ySplitError = rectError(node->rect.topHalf()) +
                                rectError(node->rect.bottomHalf());

            //std::cerr << "xSplitError " << xSplitError << "ySplitError " << ySplitError << std::endl;

            if ((xSplitError > rectMaxError) && (xSplitError >= ySplitError)) {

                node->splitX();
            }
            else if ((ySplitError > rectMaxError) && (ySplitError > xSplitError)) {

                node->splitY();
            }
        }
        else if (node->rect.canSplitX()) {

            Float xSplitError = rectError(node->rect.leftHalf()) +
                                rectError(node->rect.rightHalf());


            if (xSplitError > rectMaxError) {

                node->splitX();
            }
        }
        else if (node->rect.canSplitY()) {

            Float ySplitError = rectError(node->rect.topHalf()) +
                                rectError(node->rect.bottomHalf());


            if (ySplitError > rectMaxError) {

                node->splitY();
            }
        }
    }

    if (!node->isLeaf()) {

        refineNode(node->pChild1);
        refineNode(node->pChild2);
    }
}

Float PartitionChunk::integral(const GreedRect& rect)
{
    Float res = 0.;

    for (int i = rect.x1; i < rect.x2; ++i)
        for (int j = rect.y1; j < rect.y2; ++j)
            res += m_cellIntegrals[j][i];

    return res;
}

Float PartitionChunk::approxIntegral(const GreedRect& rect)
{
    return 0.25 * ( m_data[rect.y1][rect.x1] +
                    m_data[rect.y2][rect.x1] +
                    m_data[rect.y1][rect.x2] +
                    m_data[rect.y2][rect.x2]) *
    (rect.square * m_cellSquare);
}

Float PartitionChunk::rectError(const GreedRect& rect)
{
    Float res = 0.;
    Float a, b, c, d;

    a = m_data[rect.y1][rect.x1] / rect.square;
    b = m_data[rect.y1][rect.x2] / rect.square;
    c = m_data[rect.y2][rect.x1] / rect.square;
    d = m_data[rect.y2][rect.x2] / rect.square;

    for (int i = rect.x1; i < rect.x2; ++i)
        for (int j = rect.y1; j < rect.y2; ++j) {

            res +=  std::abs(m_data[j][i]  -  (a*(rect.x2 - i)*(rect.y2 - j)
                                           + b*(i - rect.x1)*(rect.y2 - j)
                                           + c*(rect.x2 - i)*(j - rect.y1)
                                           + d*(i - rect.x1)*(j - rect.y1)));
        }

    return res*m_cellSquare;
}


void PartitionChunk::createRectsList()
{
    processTreeNode(m_root);
}

void PartitionChunk::processTreeNode(Node* node)
{
    if (node->isLeaf()) {

        int  keys[4];
        int  indeces[4];
        Knot knots[4];
        Float values[4];

        GreedRect &r = node->rect;

        keys[0] = r.x1 + r.y1*kThetaSize; //tl
        keys[1] = r.x2 + r.y1*kThetaSize; //tr
        keys[2] = r.x1 + r.y2*kThetaSize; //bl
        keys[3] = r.x2 + r.y2*kThetaSize; //br

        Float x1 = r.x1 * kThetaResolution;
        Float x2 = r.x2 * kThetaResolution;
        Float y1 = r.y1 * kPhiResolution;
        Float y2 = r.y2 * kPhiResolution;

        values[0] = m_data[r.y1][r.x1];
        values[1] = m_data[r.y1][r.x2];
        values[2] = m_data[r.y2][r.x1];
        values[3] = m_data[r.y2][r.x2];

        knots[0] = Knot(x1, y1, values[0]);
        knots[1] = Knot(x2, y1, values[1]);
        knots[2] = Knot(x1, y2, values[2]);
        knots[3] = Knot(x2, y2, values[3]);

        for (int i = 0; i < 4; ++i) {

            if (m_knotsMap.find(keys[i]) == m_knotsMap.end()) {

                indeces[i] = m_knots.size();
                m_knots.push_back(knots[i]);
                m_knotsMap[keys[i]] = indeces[i];
            }
            else {

                indeces[i] = m_knotsMap[keys[i]];
            }
        }

        m_rects.push_back(Rect(indeces[0], indeces[1], indeces[2], indeces[3], &m_knots));
    }
    else {

        processTreeNode(node->pChild1);
        processTreeNode(node->pChild2);
    }
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

        m_rects.push_back(Rect(tl, tr, bl, br, &m_knots));
    }

    normalize();
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
        RectsVector::iterator i;

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

    RectsVector::iterator i = m_rects.begin();
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
    RectsVector::iterator i;
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
}

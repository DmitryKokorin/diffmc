#pragma once

#include <cstdlib>
#include <vector>
#include <string>
#include <tr1/unordered_map>
#include <iostream>

#include <list>
#include <algorithm>

#include "common.h"
#include "node.h"
#include "bilinear_interpolation.h"
#include "partition_rect.h"
#include "optics.h"
#include "indicatrix.h"
#include "matrix3.h"
#include "coords.h"
#include "rect_integrate.h"


typedef std::vector<PartitionRect> PartitionRectsVector;
typedef std::vector<Knot>          KnotsVector;


class PartitionChunk
{
public:

    typedef std::vector<Float> ValuesVector;

    PartitionChunk();
    virtual ~PartitionChunk();

    template <class T>
    bool create(const Float minAngle, const Float maxAngle);

    bool load(FILE *file);
    bool save(FILE *file);

    inline size_t getRectsCount() const { return m_rects.size(); }
    inline size_t getKnotsCount() const { return m_knots.size(); }
    inline bool  isAngleInRange(const Float angle) const { return angle >= m_minAngle && angle <= m_maxAngle;}

    static const Float kEpsilon;

    inline PartitionRectsVector &rects() {return m_rects;}
    inline KnotsVector &knots() {return m_knots;}
    inline ValuesVector &values() {return m_values;}


private:

    template <typename T>
    inline Float rectError(const T& func, const BilinearInterpolation &interpolation) const;

    void cleanUp();
    void normalize();
    void fillValues();

    PartitionRectsVector m_rects;
    KnotsVector m_knots;

    std::vector<Float> m_values;  //integrals upp to the n-th rect

    Float m_minAngle;
    Float m_maxAngle;

    Node*    m_root;
    std::tr1::unordered_map<Knot, int, Knot::hash> m_knotsMap;

    //used on tree creation
    template <typename T>
    void refineNode(Node *node, const T &func, const Float minIntegral, const Float currentError);

    template <typename T>
    void refineNode2(Node *node, const T &func, const Float minIntegral);

    template <typename T>
    void refine(const T &func, const Float maxError);

    void process();

    struct InterpolationInfo
    {
        InterpolationInfo()
            : interpolation()
            , error(0.)
        {}

        BilinearInterpolation interpolation;
        Float                 error;

        bool operator<(const InterpolationInfo &other)
        {
            return error > other.error;
        }
    };

    typedef std::list<InterpolationInfo> InfosContainer;

    Float currentError() const;

    InfosContainer m_infos;

    //used on list creation
    void processNode(Node* node);

    Float m_fullIntegral;

    //disable copying
    PartitionChunk(const PartitionChunk&);
    PartitionChunk& operator=(const PartitionChunk&);

    int m_rectsCount;

    static const Float kMinX, kMaxX;
    static const Float kMinY, kMaxY;
    static const Float kIntegralEpsilon;
};

template <class T>
struct PartitionFunctor
{
    PartitionFunctor(const T& _ind)
        : ind(_ind)
    {}

    inline Float operator()(const Float ct, const Float phi) const
    {
        Float st  = sqrt(1-ct*ct);
        Vector3 s = Vector3(st*cos(phi), st*sin(phi), ct);

        return ind(s);
    }

    const T &ind;
};

template <class T>
struct ErrorFunctor
{
    ErrorFunctor(const T& _func, const BilinearInterpolation &_interpolation)
        : func(_func)
        , interpolation(_interpolation)
    {}

    inline Float operator()(const Float x, const Float y) const
    {
        return std::abs(func(x, y) - interpolation(x, y));
    }

    const T &func;
    const BilinearInterpolation &interpolation;
};

template <class T>
bool PartitionChunk::create(const Float minAngle, const Float maxAngle)
{
    m_minAngle   = minAngle;
    m_maxAngle   = maxAngle;

    //create direction
    Float theta_i = 0.5*(m_minAngle + m_maxAngle);
    Angle   a_i = Angle(theta_i);
    Vector3 s_i = createSomeDeviantVector(Optics::director, a_i).normalize();

    //create coordinate system
    Vector3 v2 = crossProduct(s_i, Optics::director).normalize();
    Vector3 v3 = crossProduct(s_i, v2).normalize();

    Matrix3 mtx = createTransformMatrix(v2, v3, s_i);
    Vector3 nn  = mtx*Optics::director;

    T ind(Vector3(0., 0., 1.), nn);
    PartitionFunctor<T> func(ind);
//    m_root = new Node(BilinearInterpolation::create(func, Rect(-1., 0., 1., M_PI)));

//    ErrorFunctor<PartitionFunctor<T> > errorFunctor(func, m_root->interpolation);


//    Float rootError = rect_integrate<ErrorFunctor<PartitionFunctor<T> > >(errorFunctor,
//                                               m_root->interpolation.rect, kEpsilon);

//    Float minIntegral = rect_integrate<PartitionFunctor<T> >(func, m_root->interpolation.rect, 1e-15);

//    minIntegral *= kIntegralEpsilon;
//    refineNode(m_root, func, minIntegral, rootError);
//    processNode(m_root);
//

    m_root = new Node(BilinearInterpolation::create(func, Rect(-1., 0., 1., M_PI)));

    Float minIntegral = 1e-15*rect_integrate<PartitionFunctor<T> >(func, m_root->interpolation.rect, 1e-15);

    refineNode2(m_root, func, minIntegral);
    processNode(m_root);


/*    InterpolationInfo   root;
    root.interpolation = BilinearInterpolation::create(func, Rect(-1., 0., 1., M_PI));

    ErrorFunctor<PartitionFunctor<T> > errorFunctor(func, root.interpolation);
    root.error = rect_integrate(errorFunctor, root.interpolation.rect, 1e-15);

    m_infos.push_back(root);

    Float integral = rect_integrate(func, root.interpolation.rect, 1e-15);

    refine(func, kEpsilon*integral);
    process()*/



    normalize();

    cleanUp();
    fillValues();

    return true;
}

template <typename T>
void PartitionChunk::refineNode(Node *node, const T &func, const Float minIntegral, const Float currentError)
{
    Float err1 = 0.;
    Float err2 = 0.;

    if (node->isLeaf()) {

        Rect &r = node->interpolation.rect;
        //std::cerr << "x1="<< r.x1 << " x2=" << r.x2 << " y1=" << r.y1 << " y2=" << r.y2 << std::endl;

        const Float eps2 = 0.001;

        bool canSplitX = r.width  > (kMaxX - kMinX) / (1e+6);
        bool canSplitY = r.height > (kMaxY - kMinY) / (1e+6);

        //std::cerr << "minIntegral: " << minIntegral << std::endl;

        //std::cerr << "currentError: " << currentError << std::endl;
        Float integral = rect_integrate(func, node->interpolation.rect, eps2);
        //std::cerr << "Integral: " << integral << std::endl;
        Float maxError = kEpsilon * integral*(1. - 2.*eps2/(eps2+1.));
        //std::cerr << "maxError: " << maxError << std::endl;

        if (currentError > maxError && integral > minIntegral) {

            if (canSplitX && canSplitY) {

                BilinearInterpolation top, bottom, left, right;
                node->interpolation.splitX(func, left, right);
                node->interpolation.splitY(func, top, bottom);

                Float leftError  = rect_integrate(ErrorFunctor<T>(func, left), left.rect, eps2);
                Float rightError = rect_integrate(ErrorFunctor<T>(func, right), right.rect, eps2);
                Float xError     = leftError + rightError;
                //std::cerr << "xError: " << xError << std::endl;

                Float topError    = rect_integrate(ErrorFunctor<T>(func, top), top.rect, eps2);
                Float bottomError = rect_integrate(ErrorFunctor<T>(func, bottom), bottom.rect, eps2);
                Float yError      = topError + bottomError;
                //std::cerr << "yError: " << yError << std::endl;

                if (xError < yError) {

                        node->split(left, right);
                //        std::cerr << "split x" << std::endl;
                        err1 = leftError;
                        err2 = rightError;
                }
                else {

                        node->split(top, bottom);
                //        std::cerr << "split y" << std::endl;
                        err1 = topError;
                        err2 = bottomError;
                }
            }
            else if (canSplitX) {

                BilinearInterpolation left, right;
                node->interpolation.splitX(func, left, right);

                Float leftError  = rect_integrate(ErrorFunctor<T>(func, left), left.rect, eps2);
                Float rightError = rect_integrate(ErrorFunctor<T>(func, right), right.rect, eps2);
                //Float xError     = leftError + rightError;
                //std::cerr << "xError: " << xError << std::endl;

                node->split(left, right);
                //std::cerr << "split x" << std::endl;
                err1 = leftError;
                err2 = rightError;

            }
            else if (canSplitY) {

                BilinearInterpolation top, bottom;
                node->interpolation.splitY(func, top, bottom);

                Float topError    = rect_integrate(ErrorFunctor<T>(func, top), top.rect, eps2);
                Float bottomError = rect_integrate(ErrorFunctor<T>(func, bottom), bottom.rect, eps2);
                //Float yError      = topError + bottomError;
                //std::cerr << "yError: " << yError << std::endl;

                node->split(top, bottom);
                //std::cerr << "split y" << std::endl;
                err1 = topError;
                err2 = bottomError;
            }
        }
    }

    if (!node->isLeaf()) {

        ++m_rectsCount;

        //std::cerr << ++m_rectsCount << std::endl;

        refineNode(node->pChild1, func, minIntegral, err1);
        refineNode(node->pChild2, func, minIntegral, err2);
    }
}


template <typename T>
void PartitionChunk::refineNode2(Node *node, const T &func, const Float minIntegral)
{
    if (node->isLeaf()) {

        BilinearInterpolation &i = node->interpolation;

        if (i.integral() < minIntegral)
            return;

        Rect &r = i.rect;
        //std::cerr << "x1="<< r.x1 << " x2=" << r.x2 << " y1=" << r.y1 << " y2=" << r.y2 << std::endl;

        Float maxError = kEpsilon * std::max(std::max(std::max(i.f_x1y1, i.f_x1y2), i.f_x2y1), i.f_x2y2);
        //std::cerr << "maxError: " << maxError << std::endl;

        BilinearInterpolation left, right, top, bottom;
        node->interpolation.splitX(func, left, right);
        node->interpolation.splitY(func, top, bottom);

        Float xError     = std::max(std::abs(i(r.midx, r.y1) - left.f_x2y1), std::abs(i(r.midx, r.y2) - left.f_x2y2));
        //std::cerr << "xError: " << xError << std::endl;

        Float yError     = std::max(std::abs(i(r.x1, r.midy) - top.f_x1y2), std::abs(i(r.x2, r.midy) - top.f_x2y2));
        //std::cerr << "yError: " << yError << std::endl;

        const bool splitX = xError > maxError && (xError <  yError || yError < maxError);
        const bool splitY = yError > maxError && (yError <= xError || xError < maxError);

        if (splitX) {

            node->split(left, right);
            //std::cerr << "split x" << std::endl;
        }
        else if (splitY) {

            node->split(top, bottom);
            //std::cerr << "split y" << std::endl;
        }
    }

    if (!node->isLeaf()) {

        ++m_rectsCount;

        //std::cerr << ++m_rectsCount << std::endl;

        refineNode2(node->pChild1, func, minIntegral);
        refineNode2(node->pChild2, func, minIntegral);
    }
}

template <typename T>
void PartitionChunk::refine(const T &func, const Float maxError)
{
    const Float eps2 = 0.01;

    //int rects = 1;
    Float err;

    while ((err = currentError()) > maxError) {

        //std::cerr << "currentError : " << err << std::endl;
        //std::cerr << "maxError : " << maxError << std::endl;

        InterpolationInfo info = m_infos.front();

        InterpolationInfo top, bottom, left, right;
        info.interpolation.splitX(func, left.interpolation, right.interpolation);
        info.interpolation.splitY(func, top.interpolation, bottom.interpolation);

        left.error  = rect_integrate(ErrorFunctor<T>(func, left.interpolation), left.interpolation.rect, eps2);
        right.error = rect_integrate(ErrorFunctor<T>(func, right.interpolation), right.interpolation.rect, eps2);
        Float xError     = left.error + right.error;
        //std::cerr << "xError: " << xError << std::endl;

        top.error    = rect_integrate(ErrorFunctor<T>(func, top.interpolation), top.interpolation.rect, eps2);
        bottom.error = rect_integrate(ErrorFunctor<T>(func, bottom.interpolation), bottom.interpolation.rect, eps2);
        Float yError      = top.error + bottom.error;
        //std::cerr << "yError: " << yError << std::endl;

        if (xError < yError) {

            //std::cerr << "split x" << std::endl;
            m_infos.push_back(left);
            m_infos.push_back(right);
        }
        else {

            //std::cerr << "split y" << std::endl;
            m_infos.push_back(top);
            m_infos.push_back(bottom);
        }

        m_infos.pop_front();

        //std::cerr << ++rects << std::endl;

        m_infos.sort();
    }
}


template <typename T>
inline Float PartitionChunk::rectError(const T& func, const BilinearInterpolation &interpolation) const
{
    ErrorFunctor<T> errorFunctor(func, interpolation);
    return rect_integrate(errorFunctor, interpolation.rect, 0.1*kEpsilon, 1);
}


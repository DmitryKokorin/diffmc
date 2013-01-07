#pragma once

#include "mathcompat.h"


#include <omp.h>

#include "common.h"
#include "indicatrix.h"
#include "coords.h"
#include "linterpol.h"
#include "spherical.h"


Float symmetrizeTheta(const Float theta);


namespace freepath {

template <class T>
class Functor
{
public:

    Functor(Indicatrix<T, Optics::OBeam>& indO_, Indicatrix<T, Optics::EBeam>& indE_,
                const Vector3& nn_, const Angle& a_i_)
        : indO(indO_)
        , indE(indE_)
        , nn(nn_)
        , a_i(a_i_)
    {}

    inline Float operator()(const Vector3 &s)
    {
        Angle a_s   = Angle(s, nn);

        return  (indO(s)*Optics::OBeam::cosd(a_s)/Optics::OBeam::f2(a_s) +
                 indE(s)*Optics::EBeam::cosd(a_s)/Optics::EBeam::f2(a_s) ) / T::cosd(a_i);
    }

protected:

    Indicatrix<T, Optics::OBeam>& indO;
    Indicatrix<T, Optics::EBeam>& indE;

    const Vector3& nn;
    const Angle& a_i;
};

} //namespace


template <class T>
void createFreePath(LinearInterpolation& li, const int kPoints = 400)
{

    li.resize(kPoints);
    li.setRange(0., 0.5*M_PI);


    const Float kResolution = li.resolution();

    #pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < kPoints; ++i) {

        Float theta_i = i*kResolution;
        Angle   a_i   = Angle(theta_i);
        Vector3 s_i   = createSomeDeviantVector(Optics::director, a_i).normalize();

        //create coordinate system

        Vector3 v2 = crossProduct(s_i, Optics::director).normalize();
        Vector3 v3 = crossProduct(v2, s_i);

        Matrix3 mtx = invert(createTransformMatrix(s_i, v2, v3));
        Vector3 nn = mtx*Optics::director;

        Indicatrix<T, Optics::OBeam> indO = Indicatrix<T, Optics::OBeam>(Vector3(1., 0., 0.), nn);
        Indicatrix<T, Optics::EBeam> indE = Indicatrix<T, Optics::EBeam>(Vector3(1., 0., 0.), nn);

        freepath::Functor<T> functor(indO, indE, nn, a_i);
        Float integral = spherical::integral(functor);

        li[i] = 1./(integral);
    }

    li[0] = li[1];
}

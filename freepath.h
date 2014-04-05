#pragma once

#include "mathcompat.h"

#include <omp.h>

#include "common.h"
#include "indicatrix.h"
#include "coords.h"
#include "linterpol.h"
#include "spherical.h"

#include "rect_integrate.h"


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

#if !defined TEST
    inline Float operator()(const Vector3 &s) const
    {
        Angle a_s   = Angle(s, nn);

#if !defined EE_ONLY
        return  (indO(s)*Optics::OBeam::cosd(a_s)/Optics::OBeam::f2(a_s) +
                 indE(s)*Optics::EBeam::cosd(a_s)/Optics::EBeam::f2(a_s) ) / T::cosd(a_i);
#else
        return  (indE(s)*Optics::EBeam::cosd(a_s)/Optics::EBeam::f2(a_s) ) / T::cosd(a_i);
#endif
    }
#else
    inline Float operator()(const Vector3 &s) const
    {
        return  indE(s);
    }
#endif

protected:

    Indicatrix<T, Optics::OBeam>& indO;
    Indicatrix<T, Optics::EBeam>& indE;

    const Vector3& nn;
    const Angle& a_i;
};


template <class T>
struct RectFunctor
{
    RectFunctor(const T& _func)
        : func(_func)
    {}

    inline Float operator()(const Float ct, const Float phi) const
    {
        Float st  = sqrt(1-ct*ct);
        Vector3 s = Vector3(st*cos(phi), st*sin(phi), ct);

        return func(s);
    }

    const T &func;
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
        Angle a_i;
        Vector3 nn;

        Indicatrix<T, Optics::OBeam> indO = createIndicatrix<Indicatrix<T, Optics::OBeam> >(theta_i, a_i, nn);
        Indicatrix<T, Optics::EBeam> indE = createIndicatrix<Indicatrix<T, Optics::EBeam> >(theta_i);

        typedef freepath::Functor<T> Functor;
        Functor functor(indO, indE, nn, a_i);
//        Float integral = spherical::integral(functor);
        Float integral = rect_integrate(freepath::RectFunctor<Functor>(functor), Rect(-1., 0., 1., 2*M_PI), 1e-8, 20);

        li[i] = 1./(integral);
    }

    li[0] = li[1];
}

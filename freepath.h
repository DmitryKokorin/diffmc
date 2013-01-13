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
        Angle a_i;
        Vector3 nn;

        Indicatrix<T, Optics::OBeam> indO = createIndicatrix<Indicatrix<T, Optics::OBeam> >(theta_i, a_i, nn);
        Indicatrix<T, Optics::EBeam> indE = createIndicatrix<Indicatrix<T, Optics::EBeam> >(theta_i);

        freepath::Functor<T> functor(indO, indE, nn, a_i);
        Float integral = spherical::integral(functor);

        li[i] = 1./(integral);
    }

    li[0] = li[1];
}

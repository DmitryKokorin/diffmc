#pragma once

#include "common.h"
#include "linterpol.h"
#include "matrix3.h"
#include "indicatrix.h"
#include "spherical.h"

template <class T>
void createEChannelProb(LinearInterpolation& li, const int kPoints = 400)
{
    li.resize(kPoints);
    li.setRange(0., 0.5*M_PI);

    const Float kResolution = li.resolution();

    #pragma omp parallel for
    for (int i = 1; i < kPoints; ++i) {

        Float theta_i = i*kResolution;

        Indicatrix<T, Optics::OBeam> indO = createIndicatrix<Indicatrix<T, Optics::OBeam> >(theta_i);
        Indicatrix<T, Optics::EBeam> indE = createIndicatrix<Indicatrix<T, Optics::EBeam> >(theta_i);

        Float oIntegral = spherical::integral(indO);
        Float eIntegral = spherical::integral(indE);

        li[i] = eIntegral /(eIntegral + oIntegral);
    }

    li[0] = li[1];
}


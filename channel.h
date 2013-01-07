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

        Angle   a_i = Angle(theta_i);
        Vector3 s_i = createSomeDeviantVector(Optics::director, a_i).normalize();

        //create coordinate system

        Vector3 v2 = crossProduct(s_i, Optics::director).normalize();
        Vector3 v3 = crossProduct(v2, s_i);

        Matrix3 mtx = invert(createTransformMatrix(s_i, v2, v3));
        Vector3 nn = mtx*Optics::director;

        Indicatrix<T, Optics::OBeam> indO = Indicatrix<T, Optics::OBeam>(Vector3(1., 0., 0.), nn);
        Indicatrix<T, Optics::EBeam> indE = Indicatrix<T, Optics::EBeam>(Vector3(1., 0., 0.), nn);

        Float oIntegral = spherical::integral(indO);
        Float eIntegral = spherical::integral(indE);

        li[i] = eIntegral /(eIntegral + oIntegral);
    }

    li[0] = li[1];
}


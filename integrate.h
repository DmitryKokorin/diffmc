#pragma once

#include <algorithm>

#include "mathcompat.h"
#include "common.h"


class Adapt
{
public:

    Adapt(const Float eps = 100*kMachineEpsilon):
        m_failed(false),
        m_eps(eps),
        m_estimation()
    {}

    template <typename T>
    Float integrate(T& func, const Float l, const Float r);

    bool failed() const { return m_failed; }

private:

    template <typename T>
    Float step(T& func, const Float l,  const Float r,
                        const Float fl, const Float fr);

    bool  m_failed;
    Float m_eps;
    Float m_estimation;

    static const Float kAlpha;
    static const Float kBeta;
};

template <typename T>
Float Adapt::integrate(T& func, const Float l, const Float r)
{
    m_failed = false;

    //15-point Kronrod for integral estimation
    Float h  = 0.5*(r - l);
    Float x0 = 0.5*(r + l);

    Float f[15] = { func(x0 - h*0.991455371120813),
                    func(x0 - h*0.949107912342759),
                    func(x0 - h*0.864864423359769),
                    func(x0 - h*0.741531185599394),
                    func(x0 - h*0.586087235467691),
                    func(x0 - h*0.405845151377397),
                    func(x0 - h*0.207784955007898),
                    func(x0),
                    func(x0 + h*0.207784955007898),
                    func(x0 + h*0.405845151377397),
                    func(x0 + h*0.586087235467691),
                    func(x0 + h*0.741531185599394),
                    func(x0 + h*0.864864423359769),
                    func(x0 + h*0.949107912342759),
                    func(x0 + h*0.991455371120813) };


    Float w[8] = { 0.022935322010529,
                   0.063092092629979,
                   0.104790010322250,
                   0.140653259715525,
                   0.169004726639267,
                   0.190350578064785,
                   0.204432940075298,
                   0.209482141084728 };

    m_estimation = h*(  w[0]*(f[0] + f[14]) +
                        w[1]*(f[1] + f[13]) +
                        w[2]*(f[2] + f[12]) +
                        w[3]*(f[3] + f[11]) +
                        w[4]*(f[4] + f[10]) +
                        w[5]*(f[5] + f[9])  +
                        w[6]*(f[6] + f[8])  +
                        w[7]*f[7]
                    );

    m_estimation = fabs(m_estimation);

    return step(func, l, r, func(l), func(r));
}

template <typename T>
Float Adapt::step(T& func, const Float l,  const Float r,
                           const Float fl, const Float fr)
{
    const Float h   = 0.5*(r - l);
    const Float x0  = 0.5*(r + l);

    const Float xl2 = x0 - kAlpha*h;
    const Float xl1 = x0 - kBeta*h;
    const Float xr1 = x0 + kBeta*h;
    const Float xr2 = x0 + kAlpha*h;

    if (xl2 <= l || xr2 >= r) {

        m_failed = true;
        return 0.;  //zero-dimensional set
    }

    const Float fl2 = func(xl2);
    const Float fl1 = func(xl1);
    const Float f0  = func(x0);
    const Float fr1 = func(xr1);
    const Float fr2 = func(xr2);

    //4-point Gauss-Lobatto
    const Float int1 = h * (1./6.*(fl + fr) + 5./6.*(fl1 + fr1));

    //7-point Kronrod
    const Float int2 = h * (   11./210.*(fl  + fr)
                            +  72./245.*(fl2 + fr2)
                            + 125./294.*(fl1 + fr1)
                            +  16./35. *f0 );

    if (fabs(int1 - int2) < m_eps*h/**std::max(m_estimation, 
    fabs(int2))*/)
        return int2;

    return    step(func, l,   xl2, fl,  fl2)
            + step(func, xl2, xl1, fl2, fl1)
            + step(func, xl1, x0,  fl1, f0 )
            + step(func, x0,  xr1, f0,  fr1)
            + step(func, xr1, xr2, fr1, fr2)
            + step(func, xr2, r,   fr2, fr );
}

#pragma once

#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "mathcompat.h"
#include "common.h"


namespace integrate {

template <typename T>
Float adaptiveStep(T &func, const Float a, const Float b,
                            const Float fa, const Float fb,
                            const Float tolerance, const Float estimate, const int depth, const int maxdepth)
{
    static const Float alpha = std::sqrt(2.0/3.0);
    static const Float beta  = 1.0/std::sqrt(5.0);

    const Float m = 0.5*(a + b);
    const Float h = 0.5*(b - a);

    const Float mll = m - alpha*h;
    const Float ml  = m - beta*h;
    const Float mr  = m + beta*h;
    const Float mrr = m + alpha*h;

    const Float fmll = func(mll);
    const Float fml  = func(ml);
    const Float fm   = func(m);
    const Float fmr  = func(mr);
    const Float fmrr = func(mrr);

    const Float i2 = h/6.0*(fa + fb + 5.0*(fml + fmr));

    const Float i1 = h/1470.0*(77.0*(fa   + fb)  +
                        432.0*(fmll + fmrr)+
                        625.0*(fml  + fmr) +
                        672.0*fm);

    if (std::abs(i2 - i1) < tolerance*estimate || depth > maxdepth || mll <= a || b <= mrr) {

        if (mll <= a || b <= mrr)
            throw std::runtime_error("no more numbers");

        //std::cerr << "i1 = " << i1 << std::endl;
        return i1;
    }
    else {

        int d = depth + 1;

        return  adaptiveStep(func, a,   mll, fa,   fmll, tolerance, estimate, d, maxdepth) +
                adaptiveStep(func, mll, ml,  fmll, fml,  tolerance, estimate, d, maxdepth) +
                adaptiveStep(func, ml,  m,   fml,  fm,   tolerance, estimate, d, maxdepth) +
                adaptiveStep(func, m,   mr,  fm,   fmr,  tolerance, estimate, d, maxdepth) +
                adaptiveStep(func, mr,  mrr, fmr,  fmrr, tolerance, estimate, d, maxdepth) +
                adaptiveStep(func, mrr, b,   fmrr, fb,   tolerance, estimate, d, maxdepth);
    }
}

template <typename T>
Float adaptive(T &func, const Float a, const Float b,
                        const Float tolerance = 1.0e-10, const int maxdepth = 10)
{
    Float y[13];
    Float m = 0.5*(a + b);
    Float h = 0.5*(b - a);

    static const Float alpha = std::sqrt(2.0/3.0);
    static const Float beta  = 1.0/std::sqrt(5.0);
    static const Float x1    = 0.942882415695480;
    static const Float x2    = 0.641853342345781;
    static const Float x3    = 0.236383199662150;

    static const Float x[] = {0, -x1, -alpha, -x2, -beta, -x3, 0.0, x3, beta, x2, alpha, x1};

    y[0]  = func(a);
    y[12] = func(b);
    for (int i = 1; i < 12; ++i)
        y[i]=func(m+x[i]*h);

    Float estimate = h*
        (0.0158271919734802*(y[0]+y[12])
        +0.0942738402188500*(y[1]+y[11])
        +0.155071987336585 *(y[2]+y[10])
        +0.188821573960182 *(y[3]+y[9])
        +0.199773405226859 *(y[4]+y[8])
        +0.224926465333340 *(y[5]+y[7])
        +0.242611071901408 *y[6]);

    estimate = std::abs(estimate);
    if (0.0 == estimate)
        estimate = b - a;

    //std::cerr << estimate << std::endl;


    return adaptiveStep(func, a, b, y[0], y[12], tolerance, estimate, 0, maxdepth);
}


}

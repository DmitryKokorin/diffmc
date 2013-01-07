#pragma once

#include <algorithm>
#include <stdexcept>

#include "mathcompat.h"
#include "common.h"


namespace integrate {

template <typename T>
Float adaptiveStep(T &func, const Float a, const Float b,
                            const Float fa, const Float fb,
                            const Float absTolerance, const int depth, const int maxdepth)
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

    if (std::abs(i2 - i1) < absTolerance || depth > maxdepth || mll <= a || b <= mrr) {

        if (mll <= a || b <= mrr)
            throw std::runtime_error("no more numbers");

        return i1;
    }
    else {

        Float ih = 1.0/h;
        int d = depth + 1;

        return  adaptiveStep(func, a,   mll, fa,   fmll, absTolerance*(mll-a)*ih,    d, maxdepth) +
                adaptiveStep(func, mll, ml,  fmll, fml,  absTolerance*(ml - mll)*ih, d, maxdepth) +
                adaptiveStep(func, ml,  m,   fml,  fm,   absTolerance*(m - ml)*ih,   d, maxdepth) +
                adaptiveStep(func, m,   mr,  fm,   fmr,  absTolerance*(mr-m)*ih,     d, maxdepth) +
                adaptiveStep(func, mr,  mrr, fmr,  fmrr, absTolerance*(mrr - mr)*ih, d, maxdepth) +
                adaptiveStep(func, mrr, b,   fmrr, fb,   absTolerance*(b-mrr)*ih,    d, maxdepth);
    }
}

template <typename T>
Float adaptive(T &func, const Float a, const Float b,
                        const Float absTolerance, const int maxdepth = 10)
{
    return adaptiveStep(func, a, b, func(a), func(b), absTolerance, 0, maxdepth);
}


}

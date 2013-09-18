#pragma once

#include "integrate.h"
#include "rect.h"

namespace detail {


template <class T>
class YFunctor
{
public:

    YFunctor(const T &_func, const Float _x)
      : func(_func)
      , x(_x)
    {}

    inline Float operator()(const Float y) const
    {
        return func(x, y);
    }

protected:

    const T &func;
    const Float x;
};


template <class T>
class XFunctor
{
public:

    XFunctor(const T &_func, const Rect &_rect, const double _tolerance = 1.0e-10, const int _max_depth = 10)
      : func(_func)
      , rect(_rect)
      , tolerance(_tolerance)
      , max_depth(_max_depth)
    {}

    inline Float operator()(const Float x) const
    {
        YFunctor<T> yFunctor(func, x);
        return integrate::adaptive(yFunctor, rect.y1, rect.y2, tolerance, max_depth);
    }

protected:

    const T &func;
    const Rect &rect;
    const double tolerance;
    const int max_depth;
};

}

template <class T>
Float rect_integrate(const T &_func, const Rect &_rect, const double _tolerance = 1.0e-10, const int _max_depth = 10)
{
    using namespace detail;

    const XFunctor<T> xFunctor(_func, _rect, _tolerance, _max_depth);
    return integrate::adaptive(xFunctor, _rect.x1, _rect.x2, _tolerance, _max_depth);
}

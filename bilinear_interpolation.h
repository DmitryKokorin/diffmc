#pragma once

#include "common.h"
#include "rect.h"

struct BilinearInterpolation
{
    BilinearInterpolation()
        : rect()
        , f_x1y1()
        , f_x1y2()
        , f_x2y1()
        , f_x2y2()
    {
    }

    BilinearInterpolation(const Rect &_rect,
                          const Float& _f_x1y1, const Float& _f_x1y2, const Float& _f_x2y1, const Float& _f_x2y2)
        : rect(_rect)
        , f_x1y1(_f_x1y1)
        , f_x1y2(_f_x1y2)
        , f_x2y1(_f_x2y1)
        , f_x2y2(_f_x2y2)
    {
    }

    Float operator()(const Float x, const Float y) const
    {
        return (f_x1y1*(rect.x2 - x)*(rect.y2 - y)
             + f_x2y1*(x - rect.x1)*(rect.y2 - y)
             + f_x1y2*(rect.x2 - x)*(y - rect.y1)
             + f_x2y2*(x - rect.x1)*(y - rect.y1)) / rect.square;
    }

    Float integral() const
    {
        return 0.25 * rect.square * (f_x1y1 + f_x2y1 + f_x1y2 + f_x2y2);
    }

    template <typename T>
    static BilinearInterpolation create(const T &_func, const Rect &_rect)
    {
        return BilinearInterpolation(_rect, _func(_rect.x1, _rect.y1), _func(_rect.x1, _rect.y2),
                                            _func(_rect.x2, _rect.y1), _func(_rect.x2, _rect.y2));
    }

    template <typename T>
    inline void splitX(const T &_func, BilinearInterpolation &_leftHalf, BilinearInterpolation &_rightHalf)
    {
        Float f_midxy1 = _func(rect.midx, rect.y1);
        Float f_midxy2 = _func(rect.midx, rect.y2);

        _leftHalf = BilinearInterpolation(rect.leftHalf(), f_x1y1,   f_x1y2,
                                                           f_midxy1, f_midxy2);

        _rightHalf = BilinearInterpolation(rect.rightHalf(), f_midxy1, f_midxy2,
                                                             f_x2y1,   f_x2y2);
     }

    template <typename T>
    inline void splitY(const T &_func, BilinearInterpolation &_topHalf, BilinearInterpolation &_bottomHalf)
    {
        Float f_x1midy = _func(rect.x1, rect.midy);
        Float f_x2midy = _func(rect.x2, rect.midy);

        _topHalf = BilinearInterpolation(rect.topHalf(), f_x1y1, f_x1midy,
                                                         f_x2y1, f_x2midy);

        _bottomHalf = BilinearInterpolation(rect.bottomHalf(), f_x1midy, f_x1y2,
                                                               f_x2midy, f_x2y2);
     }

    Rect rect;
    Float f_x1y1, f_x1y2, f_x2y1, f_x2y2;
};

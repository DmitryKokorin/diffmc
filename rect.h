#pragma once

#include <vector>
#include "common.h"

struct Knot
{
    Knot(const Float _x, const Float _y)
        : x(_x)
        , y(_y)
    {};

    Knot()
        : x(0.)
        , y(0.)
    {};

    Float x;
    Float y;
};

struct Rect
{

    typedef std::vector<Knot> KnotsVector;

    Rect(const int tl, const int tr,
         const int bl, const int br, KnotsVector *_knots);

    int tl, tr,
        bl, br;

    Float x1, x2, y1, y2;

    Float   width,
            height;

    Float   square;
};

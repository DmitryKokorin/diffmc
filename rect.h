#pragma once

#include <vector>
#include "common.h"

struct Knot
{
    Knot(const Float _x, const Float _y, const Float _value)
        : x(_x)
        , y(_y)
        , value(_value)
    {};

    Knot()
        : x()
        , y()
        , value()
    {};

    Float x;
    Float y;
    Float value;
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

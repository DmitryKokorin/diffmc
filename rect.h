#pragma once

#include "common.h"


class Rect
{
public:
    Rect()
      : x1()
      , y1()
      , x2()
      , y2()
      , midx()
      , midy()
      , width()
      , height()
      , square()
    {}

    Rect(const Float& x1_, const Float& y1_, const Float& x2_, const Float& y2_) :
    x1(x1_),
    y1(y1_),
    x2(x2_),
    y2(y2_),
    midx((x1+x2)/2),
    midy((y1+y2)/2),
    width(x2-x1),
    height(y2-y1),
    square(width*height)
    {
    }

    Rect topHalf()    const { return Rect(x1, y1, x2, midy); }
    Rect bottomHalf() const { return Rect(x1, midy, x2, y2); }
    Rect leftHalf()   const { return Rect(x1, y1, midx, y2); }
    Rect rightHalf()  const { return Rect(midx, y1, x2, y2); }

    Float x1, y1, x2, y2;
    Float midx, midy;
    Float width, height;
    Float square;
};

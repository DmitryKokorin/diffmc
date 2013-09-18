#include "partition_rect.h"


PartitionRect::PartitionRect(const int _tl,
           const int _tr,
           const int _bl,
           const int _br,
           KnotsVector *knots)
    : tl(_tl)
    , tr(_tr)
    , bl(_bl)
    , br(_br)
    , x1(0.), x2(0.), y1(0.), y2(0.), width(0.), height(0.), square(0.)
{
    KnotsVector &k = *knots;

    x1 = k[tl].x;
    x2 = k[tr].x;
    y1 = k[tl].y;
    y2 = k[bl].y;

    width  = x2 - x1;
    height = y2 - y1;

    square = width*height;
}

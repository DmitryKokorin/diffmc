#pragma once

#include <cstdlib>
#include <vector>
#include <tr1/functional>
#include "common.h"

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::tr1::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}


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

    bool operator==(const Knot &other) const
    {
        return x == other.x && y == other.y && value == other.value;
    }

    struct hash
    {
        inline size_t operator()(const Knot &k) const
        {
          size_t seed = 0;
          ::hash_combine(seed, k.x);
          ::hash_combine(seed, k.y);
          return seed;
        }
    };

    Float x;
    Float y;
    Float value;
};

struct PartitionRect
{

    typedef std::vector<Knot> KnotsVector;

    PartitionRect(const int tl, const int tr,
         const int bl, const int br, KnotsVector *_knots);

    int tl, tr,
        bl, br;

    Float x1, x2, y1, y2;

    Float   width,
            height;

    Float   square;
};




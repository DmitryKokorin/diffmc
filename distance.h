#pragma once

#include "spherical.h"

namespace detail {

template <class T>
class DistanceFunctor
{
    Float f1_norm, f2_norm;
    T &func1;
    T &func2;

public:

    DistanceFunctor(T &func1_, const Float f1_norm_,
                    T &func2_, const Float f2_norm_)
      : f1_norm(f1_norm_)
      , f2_norm(f2_norm_)
      , func1(func1_)
      , func2(func2_)
    {}

    Float operator()(const Vector3 &s)
    {
        return std::abs(func1(s)/f1_norm - func2(s)/f2_norm);
    }
};

} //namespace detail

//calculates distance between two functions defined on a sphere
template <class T>
Float distance(T &func1, T &func2)
{
    using namespace detail;
    using namespace spherical;

    Float f1_norm = integral(func1);
    Float f2_norm = integral(func2);

    DistanceFunctor<T> distanceFunctor = DistanceFunctor<T>(func1, f1_norm, func2, f2_norm);

    return integral(distanceFunctor);
}


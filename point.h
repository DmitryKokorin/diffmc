#pragma once

#include <ostream>
#include "vector3.h"


struct Point
{
    Point()
      : r(), r2(), r3(), r4(), r5(), r6(), time(), scatterings(), measurements()
    {
    }

    Vector3 r, r2, r3, r4, r5, r6;

    Float time;

    Float   scatterings;
    int     measurements;

    void clear()
    {
        r.clear();
        r2.clear();
        r3.clear();
        r4.clear();
        r5.clear();
        r6.clear();

        scatterings  = 0.;
        measurements = 0;
    }

    Point& operator+=(const Point& rhv)
    {
        r  += rhv.r;
        r2 += rhv.r2;
        r3 += rhv.r3;
        r4 += rhv.r4;
        r5 += rhv.r5;
        r6 += rhv.r6;

        scatterings  += rhv.scatterings;
        measurements += rhv.measurements;

        return *this;
    }

    inline void appendPhoton(const Vector3 &_r, const int _scatterings)
    {
        Vector3 _r2(_r.x()*_r.x(),   _r.y()*_r.y(),   _r.z()*_r.z());
        Vector3 _r3(_r2.x()*_r.x(),  _r2.y()*_r.y(),  _r2.z()*_r.z());
        Vector3 _r4(_r2.x()*_r2.x(), _r2.y()*_r2.y(), _r2.z()*_r2.z());
        Vector3 _r5(_r3.x()*_r2.x(), _r3.y()*_r2.y(), _r3.z()*_r2.z());
        Vector3 _r6(_r3.x()*_r3.x(), _r3.y()*_r3.y(), _r3.z()*_r3.z());

        r  += _r;
        r2 += _r2;
        r3 += _r3;
        r4 += _r4;
        r5 += _r5;
        r6 += _r6;

        scatterings += _scatterings;
        ++measurements;
    }

    void average()
    {
        r  /= measurements;
        r2 /= measurements;
        r3 /= measurements;
        r4 /= measurements;
        r5 /= measurements;
        r6 /= measurements;

        scatterings /= measurements;
    }

    friend std::ostream &operator<<(std::ostream &stream, const Point &p)
    {
        return stream   << p.r << "\t"
                        << p.r2 << "\t"
                        << p.r3 << "\t"
                        << p.r4 << "\t"
                        << p.r5 << "\t"
                        << p.r6 << "\t"
                        << p.scatterings << "\t"
                        << p.measurements;
    }
};


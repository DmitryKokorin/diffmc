#pragma once

#include <ostream>

#include "point.h"
#include "optics.h"

class DataBuff
{
public:

    typedef ::Point Point;
    typedef std::vector<Point> Points;

    DataBuff()
      : points()
    {}

    DataBuff(const int pointsNum, const Float maxTime)
      : points()
    {
        prepare(pointsNum, maxTime);
    }

    void prepare(const int pointsNum, const Float maxTime)
    {
        points.clear();
        points.reserve(pointsNum);

        for (int i = 0; i < pointsNum; ++i) {

            Point pt;
            pt.time = i * (maxTime / pointsNum);
            points.push_back(pt);
        }
    }

    void clear()
    {
        Points::iterator i = points.begin();
        for (; i != points.end(); ++i) {

            i->clear();
        }
    }

    DataBuff& operator+=(const DataBuff& rhv)
    {
        Points::iterator i = points.begin();
        Points::const_iterator j = rhv.points.begin();

        for (; i != points.end(); ++i, ++j) {

            *i += *j;
        }

        return *this;
    }

    void average()
    {
        Points::iterator i = points.begin();
        for (; i != points.end(); ++i) {

            if (i->measurements)
                i->average();
        }
    }


    friend std::ostream &operator<< (std::ostream &stream, const DataBuff &buff)
    {
        DataBuff::Points::const_iterator i = buff.points.begin();
        for (; i != buff.points.end(); ++i) {

            if (i->measurements)
                stream << i->time / Optics::c << "\t" << *i << std::endl;
        }

        return stream;
    }

    Points points;
};

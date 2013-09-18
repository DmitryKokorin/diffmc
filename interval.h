#pragma once

struct Interval
{
    Interval(const Float _left, const Float _right)
      : left(_left), right(_right)
    {}

    Interval()
      : left(), right()
    {}

    Float left;
    Float right;

    bool operator<(const Interval &other) const
    {
        return left < other.left;
    }

    inline Float mid() const {return 0.5*(left + right);}
};

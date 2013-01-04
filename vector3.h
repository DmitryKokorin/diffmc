#pragma once

#include <ostream>
#include "common.h"



class Vector3
{
public:

    Vector3() :
        x_(0.0), y_(0.0), z_(0.0)
    {}

    Vector3(const Vector3& rhv) :
        x_(rhv.x_), y_(rhv.y_), z_(rhv.z_)
    {}

    Vector3(const Float x, const Float y, const Float z) :
        x_(x), y_(y), z_(z)
    {}

    inline Vector3& operator= (const Vector3& rhv);
    inline Vector3& operator+=(const Vector3& rhv);
    inline Vector3  operator+ (const Vector3& rhv) const;
    inline Vector3& operator-=(const Vector3& rhv);
    inline Vector3  operator- (const Vector3& rhv) const;

    inline Vector3& operator*=(const Float& rhv);
    inline Vector3  operator* (const Float& rhv) const;
    inline Vector3& operator/=(const Float& rhv);
    inline Vector3  operator/ (const Float& rhv) const;

    inline Float    operator*(const Vector3& rhv) const;

    inline Float x() const;
    inline Float y() const;
    inline Float z() const;

    inline Float norm() const;
    inline void clear();
    inline Vector3& normalize();

    friend std::ostream &operator<<(std::ostream &stream, const Vector3 &v)
    {
        return stream << v.x() << "\t" << v.y() << "\t" << v.z();
    }


private:

    Float x_, y_, z_;
};

inline Vector3 operator*(const Float& lhv, const Vector3& rhv);
inline Vector3 operator-(const Vector3& v);
inline Vector3 crossProduct(const Vector3& lhv, const Vector3& rhv);


#include "vector3_inl.h"

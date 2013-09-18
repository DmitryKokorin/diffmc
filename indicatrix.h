#pragma once

#include "mathcompat.h"
#include "vector3.h"
#include "angle.h"
#include "coords.h"

#include "optics.h"

template <class T1, class T2>
class Indicatrix
{
public:

    Indicatrix(const Vector3& s_i, const Vector3& director_i);
    Float operator()(const Vector3& s_s) const;

protected:

    Vector3   s_i;
    Vector3   k_i;
    Vector3   director_i;   //director in k_i based coordinate system
    Angle     a_i;          //angle between k_i and n vectors 
    Vector3   e_i;

    Float     ei_n;
    Float     factor1;

    static const Float kCalculationEpsilon;
};

template <class T1, class T2>
const Float Indicatrix<T1, T2>::kCalculationEpsilon = kMachineEpsilon;


template <class T1, class T2>
Indicatrix<T1, T2>::Indicatrix(const Vector3& s_i_, const Vector3& director_i_) 
    : s_i(s_i_)
    , k_i()
    , director_i(director_i_)
    , a_i(s_i, director_i_)
    , e_i()
    , ei_n()
    , factor1()
{
    k_i  = T1::k(s_i, a_i);
    e_i  = T1::e(s_i, director_i, a_i);

    ei_n = e_i*director_i;
    factor1 = Optics::s0/(T1::n(a_i)*T1::cosd(a_i));
}

#if !defined TEST

template <class T1, class T2>
Float Indicatrix<T1, T2>::operator()(const Vector3& s_s) const
{
    const Angle a_s     = Angle(s_s, director_i);
    const Vector3 k_s   = T2::k(s_s, a_s);
    const Float cosd_s  = T2::cosd(a_s);

    Float res = factor1 * T2::f2(a_s) * T2::n(a_s)/(cosd_s*cosd_s*cosd_s);

    const Vector3 q       = k_s - k_i;
    const Vector3 q_par   = director_i*(q*director_i);
    const Vector3 q_perp  = q - q_par;

    Vector3 a1 = q_perp;

    if (a1.norm() < q.norm()*kCalculationEpsilon) { //along the optical axis

        //we need some a1 vector here, any unit vector that is perpendicular to n
        a1 = createSomePerpendicular(director_i);
    }

    a1.normalize();

    Vector3 a2 = crossProduct(director_i, a1);
    a2.normalize();   //to be sure

    const Vector3 e_s  = T2::e(s_s, director_i, a_s);
    const Float   es_n = e_s*director_i;

    const Float es_a1  = e_s*a1;
    const Float es_a2  = e_s*a2;
    const Float ei_a1  = e_i*a1;
    const Float ei_a2  = e_i*a2;

    const Float q_perp_q_perp   = q_perp*q_perp;
    const Float q_par_q_par_add = q_par*q_par + Optics::add;


    res *= (ei_a1*ei_a1*es_n*es_n + 2.*ei_a1*es_n*ei_n*es_a1 + es_a1*es_a1*ei_n*ei_n) /
                   (Optics::t1*q_perp_q_perp + q_par_q_par_add) +
           (ei_a2*ei_a2*es_n*es_n + 2.*ei_a2*es_n*ei_n*es_a2 + es_a2*es_a2*ei_n*ei_n) /
                   (Optics::t2*q_perp_q_perp + q_par_q_par_add);

    return res;
}

#else

template <class T1, class T2>
Float Indicatrix<T1, T2>::operator()(const Vector3& s_s) const
{
    Angle a_s     = Angle(s_s, director_i);
    Vector3 k_s   = T2::k(s_s, a_s);

    Float res = factor1*T2::n(a_s);

    Vector3 q       = k_s - k_i;

    res /= (q*q + Optics::add);

    return res;
}


#endif


typedef Indicatrix<Optics::EBeam, Optics::EBeam> IndicatrixEE;
typedef Indicatrix<Optics::EBeam, Optics::OBeam> IndicatrixEO;
typedef Indicatrix<Optics::OBeam, Optics::EBeam> IndicatrixOE;
typedef Indicatrix<Optics::OBeam, Optics::OBeam> IndicatrixOO;  //should always be equal to 0.


//creates some indicatrix that has theta_i angle to director
template<typename T>
T createIndicatrix(const Float theta_i, Angle &a_i, Vector3 &nn)
{
    a_i = Angle(theta_i);
    Vector3 s_i = createSomeDeviantVector(Optics::director, a_i).normalize();

    //create coordinate system

    Vector3 v2 = crossProduct(s_i, Optics::director).normalize();
    Vector3 v3 = crossProduct(v2, s_i);

    Matrix3 mtx = createTransformMatrix(s_i, v2, v3);
    nn  = mtx*Optics::director;

    return T(Vector3(1., 0., 0.), nn);

//    Matrix3 mtx = invert(createTransformMatrix(v2, v3, s_i));
//   nn  = mtx*Optics::director;

//    return T(Vector3(0., 0., 1.), nn);
}


template<typename T>
T createIndicatrix(const Float theta_i)
{
    Angle a_i;
    Vector3 nn;

    return createIndicatrix<T>(theta_i, a_i, nn);
}


#include "optics.h"


namespace Optics {

//physical constants

const Vector3 director = Vector3(1., 0., 0.);

#define STARK_LUBENSKY_PRE_1997 0
#define STARK_VAN_TIGGELEN_RMP_2000_1 0
#define STARK_VAN_TIGGELEN_RMP_2000_2 0
#define VAN_TIGGELEN_MCLC_1997_PAA 0
#define VAN_TIGGELEN_MCLC_1997_MBBA 0
#define HEIDERICH_1997_1 0
#define HEIDERICH_1997_2 0
#define KAO_1997 0
#define WIERSMA_1999 0
#define FIGURES 0
#define FIGURES_V3 1

#if STARK_LUBENSKY_PRE_1997

const Float eps_par  = 2.923;
const Float eps_perp = 2.381;
const Float eps_a = eps_par - eps_perp;

const Float _no = sqrt(eps_perp);

const Float K3 = 5.3e-7;
const Float K1 = 4.187e-7;
const Float K2 = 2.279e-7;

const Float t1 = K1/K3;
const Float t2 = K2/K3;

const Float lambda = 5.464e-5; //stark97, calculation
const Float k0 = 2.*M_PI / lambda;

const Float T = 300.;
const Float Kb = 1.38e-16;

const Float H = 0.5e+4;
const Float Hi_alpha = 0.95e-7;
const Float xi = sqrt(K3/Hi_alpha)/H;


#elif STARK_VAN_TIGGELEN_RMP_2000_1

const Float eps_par  = 2.923;
const Float eps_perp = 2.381;
const Float eps_a = eps_par - eps_perp;

const Float _no = sqrt(eps_perp);

const Float K3 = 5.3e-7;
const Float K1 = 2.3e-7;
const Float K2 = 2.3e-7;

const Float t1 = K1/K3;
const Float t2 = K2/K3;

const Float lambda = 5.464e-5; //stark97, calculation
const Float k0 = 2.*M_PI / lambda;

const Float T = 300.;
const Float Kb = 1.38e-16;

const Float H = 0.0;
const Float xi = 1.8e-4;

#elif STARK_VAN_TIGGELEN_RMP_2000_2

const Float eps_par  = 2.923;
const Float eps_perp = 2.381;
const Float eps_a = eps_par - eps_perp;

const Float _no = sqrt(eps_perp);

const Float K3 = 5.3e-7;
const Float K1 = 4.187e-7;
const Float K2 = 2.279e-7;

const Float t1 = K1/K3;
const Float t2 = K2/K3;

const Float lambda = 5.464e-5; //stark97, calculation
const Float k0 = 2.*M_PI / lambda;

const Float T = 300.;
const Float Kb = 1.38e-16;

const Float H = 0.0;
const Float xi = 1.8e-4;

#elif VAN_TIGGELEN_MCLC_1997_PAA

const Float eps_par  = 3.35;
const Float eps_perp = 2.47;
const Float eps_a = eps_par - eps_perp;

const Float _no = sqrt(eps_perp);

const Float K3 = 9.51e-7;
const Float K1 = 3.0e-7;
const Float K2 = 3.0e-7;

const Float t1 = K1/K3;
const Float t2 = K2/K3;

const Float lambda = 5.0e-5; 
const Float k0 = 2.*M_PI / lambda;

const Float T = 400.;
const Float Kb = 1.38e-16;

const Float H = 0.0;
const Float xi = 2.2e-4;

#elif VAN_TIGGELEN_MCLC_1997_MBBA

const Float eps_par  = 4.7;
const Float eps_perp = 5.4;
const Float eps_a = eps_par - eps_perp;

const Float _no = sqrt(eps_perp);

const Float K3 = 7.45e-7;
const Float K1 = 3.7e-7;
const Float K2 = 3.7e-7;

const Float t1 = K1/K3;
const Float t2 = K2/K3;

const Float lambda = 5.0e-5; 
const Float k0 = 2.*M_PI / lambda;

const Float T = 300.;
const Float Kb = 1.38e-16;

const Float H = 0.0;
const Float xi = 2.2e-4;

#elif HEIDERICH_1997_1

const Float eps_perp = 4.7;
const Float eps_a    = 0.3 * eps_perp;
const Float eps_par  = eps_perp + eps_a;

const Float _no = sqrt(eps_perp);

const Float K3 = 8.4e-7;
const Float K1 = 6.7e-7;
const Float K2 = 6.7e-7;

const Float t1 = K1/K3;
const Float t2 = K2/K3;

const Float lambda = 2.31e-5; 
const Float k0 = 2.*M_PI / lambda;

const Float T = 300.;
const Float Kb = 1.38e-16;

const Float H = 0.0;
const Float xi = 1.23e-4;

#elif HEIDERICH_1997_2

const Float eps_perp = 4.7;
const Float eps_a    = -0.3 * eps_perp;
const Float eps_par  = eps_perp + eps_a;

const Float _no = sqrt(eps_perp);

const Float K3 = 8.4e-7;
const Float K1 = 6.7e-7;
const Float K2 = 6.7e-7;

const Float t1 = K1/K3;
const Float t2 = K2/K3;

const Float lambda = 2.31e-5; 
const Float k0 = 2.*M_PI / lambda;

const Float T = 300.;
const Float Kb = 1.38e-16;

const Float H = 0.0;
const Float xi = 1.23e-4;

#elif KAO_1997

const Float eps_par  = 2.923;
const Float eps_perp = 2.381;
const Float eps_a = eps_par - eps_perp;

const Float _no = sqrt(eps_perp);

const Float K3 = 5.3e-7;
const Float K1 = 4.187e-7;
const Float K2 = 2.279e-7;

const Float t1 = K1/K3;
const Float t2 = K2/K3;

const Float lambda = 5.145e-5; 
const Float k0 = 2.*M_PI / lambda;

const Float T = 303.15;
const Float Kb = 1.38e-16;

const Float H = 0.2e+4;
const Float Hi_alpha = 0.95e-7;
const Float xi = sqrt(K3/Hi_alpha)/H;

#elif WIERSMA_1999

const Float eps_par  = 2.923;
const Float eps_perp = 2.381;
const Float eps_a = eps_par - eps_perp;

const Float _no = sqrt(eps_perp);

const Float K3 = 7.5e-7;
const Float K1 = 5.93e-7;
const Float K2 = 3.23e-7;

const Float t1 = K1/K3;
const Float t2 = K2/K3;

const Float lambda = 4.05e-5;
const Float k0 = 2.*M_PI / lambda;

const Float T = 300;
const Float Kb = 1.38e-16;

const Float H = 0.5e+4;
const Float Hi_alpha = 1.1e-7;
const Float xi = sqrt(K3/Hi_alpha)/H;

#elif FIGURES

const Float eps_par  = 2.923;
const Float eps_perp = 2.381;
const Float eps_a = eps_par - eps_perp;

const Float _no = sqrt(eps_perp);

const Float K3 = 5.3e-7;
const Float K1 = 4.187e-7;
const Float K2 = 2.279e-7;

const Float t1 = K1/K3;
const Float t2 = K2/K3;

const Float lambda = 5.145e-5; 
const Float k0 = 2.*M_PI / lambda;

const Float T = 303.15;
const Float Kb = 1.38e-16;

const Float H = 9.0e+4;
const Float Hi_alpha = 0.95e-7;
const Float xi = sqrt(K3/Hi_alpha)/H;

#elif FIGURES_V3

const Float eps_par = 3.0;
const Float eps_perp = 2.2;
const Float eps_a = eps_par - eps_perp;

const Float _no = sqrt(eps_perp);

const Float K3 = 7.5e-7;
const Float t1 = 0.79;
const Float t2 = 0.43;
const Float K1 = t1*K3;
const Float K2 = t2*K3;

const Float lambda = 5.145e-5; 
const Float k0 = 2.*M_PI / lambda;

const Float T = 301.;
const Float Kb = 1.38e-16;

const Float H = 0.2e+4;
const Float Hi_alpha = 1.1e-7;
const Float xi = sqrt(K3/Hi_alpha)/H;


#if 0

const Float eps_par = 3.0;
const Float eps_perp = 2.2;
const Float eps_a = eps_par - eps_perp;

const Float _no = sqrt(eps_perp);


//const Float K3 = 6.1e-7;
/*const Float K3 = 7.5e-7;
const Float t1 = 0.79;
const Float t2 = 0.43;
const Float K1 = t1*K3;
const Float K2 = t2*K3;
*/
/*
const Float sum = 6.1e-7*(1. + 0.79 + 0.43);
const Float K3 = sum*0.7;
const Float K1 = sum*0.2;
const Float K2 = sum*0.1;
*/
const Float sum = 7.5e-7*(1. + 0.79 + 0.43);
const Float K3 = sum/3.0;
const Float K1 = sum/3.0;
const Float K2 = sum/3.0;


const Float t1 = K1/K3;
const Float t2 = K2/K3;


//const Float lambda = 4.88e-5; //cm  scatmc
const Float lambda = 5.145e-5; //stark97, experimental
//const Float lambda = 5.464e-5; //stark97, calculation
//const Float lambda = 4.05e-5;  //wiersma99
//const Float lambda = 5.50e-5;
const Float k0 = 2.*M_PI / lambda;
//const Float xi = 4.2e-4;

const Float T = 301.;
const Float Kb = 1.38e-16;

const Float H = 35.0e+4;
//const Float Hi_alpha = 5e-6;
const Float Hi_alpha = 1.1e-7;
const Float xi = sqrt(K3/Hi_alpha)/H;
#endif
#endif

const Float c = 29979245800; //cm/s 299 792 458

#if TEST
const Float la = 0.25;
#endif

//precalculated constants

const Float eps_perp2 = eps_perp*eps_perp;
const Float eps_par2 = eps_par*eps_par; 
const Float eps_perp_eps_par = eps_perp*eps_par;
const Float eps_par_eps_perp2 = eps_par*eps_perp2;


const Float s0 = (0.25*Kb*T*eps_a*eps_a)/(lambda*lambda*K3);
const Float add = 0.25*lambda*lambda/(M_PI*M_PI*xi*xi);

}  //namespace Optics

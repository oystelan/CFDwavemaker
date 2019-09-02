#ifndef Stokes5_H
#define Stokes5_H

#define eta Stokes5_eta // rename functions
#define u Stokes5_u
#define v Stokes5_v

#include<stdio.h>
#include<math.h>

typedef struct{
/* Coefficients depending on wave and environmental properties */
  double A11, A22, A31, A33, A42, A44, A51, A53, A55,
  B22, B31, B42, B44, B53, B55,
  C0, C2, C4,
  D2, D4,
  E2, E4,
  S, eps, wave_length, wave_height, current, gravity, depth, k, kd, kH, c, cs, u_mean, Q, R,
  x0, y0;
} Stokes5;


extern void set_stokes5_properties(Stokes5 *stokes5, double wave_length, double wave_height, double current, double depth, double gravity, double x0, double y0);
extern double eta(Stokes5 *stokes5, double t, double X);
extern double u(Stokes5 *stokes5, double t, double X, double Y);
extern double v(Stokes5 *stokes5, double t, double X, double Y);

#endif



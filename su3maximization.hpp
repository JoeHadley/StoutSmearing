#ifndef _SU3MAXIMIZATION_HPP
#define _SU3MAXIMIZATION_HPP

#include <complex>
#include <lapacke.h>
#include "color_algebra.hpp"
#include <numbers>

// Maximizes ReTr(S^+*U) over all S in SU(3).
// If S is a sum of SU(3) matrices U_i, this is equivalent to finding the point U in SU(3) 
//that minimizes the sum of geodesic distances to all U_i.
int MaximizeSU3Overlap(t_complex* S, t_complex* U, int max_iter=1000, double precision=1e-14);

#endif // _SU3MAXIMIZATION_HPP
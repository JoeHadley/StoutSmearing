#ifndef _SPINOR_ALGEBRA_HPP
#define _SPINOR_ALGEBRA_HPP

using namespace std;

#include <complex>

typedef complex<double>       t_complex;

constexpr int NSPINOR = 4;

const t_complex* gamma_matrix(int mu);

inline void s4x4_times_s4(t_complex* res, const t_complex* A, const t_complex* v) 
{
    res[0] = A[0]*v[0] + A[1]*v[1] + A[2]*v[2] + A[3]*v[3];
    res[1] = A[4]*v[0] + A[5]*v[1] + A[6]*v[2] + A[7]*v[3];
    res[2] = A[8]*v[0] + A[9]*v[1] + A[10]*v[2] + A[11]*v[3];
    res[3] = A[12]*v[0] + A[13]*v[1] + A[14]*v[2] + A[15]*v[3];
}

inline void s4x4_times_s4x4(t_complex* res, const t_complex* A, const t_complex* B) 
{
    res[0] = A[0]*B[0] + A[1]*B[4] + A[2]*B[8] + A[3]*B[12];
    res[1] = A[0]*B[1] + A[1]*B[5] + A[2]*B[9] + A[3]*B[13];
    res[2] = A[0]*B[2] + A[1]*B[6] + A[2]*B[10] + A[3]*B[14];
    res[3] = A[0]*B[3] + A[1]*B[7] + A[2]*B[11] + A[3]*B[15];

    res[4] = A[4]*B[0] + A[5]*B[4] + A[6]*B[8] + A[7]*B[12];
    res[5] = A[4]*B[1] + A[5]*B[5] + A[6]*B[9] + A[7]*B[13];
    res[6] = A[4]*B[2] + A[5]*B[6] + A[6]*B[10] + A[7]*B[14];
    res[7] = A[4]*B[3] + A[5]*B[7] + A[6]*B[11] + A[7]*B[15];

    res[8] = A[8]*B[0] + A[9]*B[4] + A[10]*B[8] + A[11]*B[12];
    res[9] = A[8]*B[1] + A[9]*B[5] + A[10]*B[9] + A[11]*B[13];
    res[10] = A[8]*B[2] + A[9]*B[6] + A[10]*B[10] + A[11]*B[14];
    res[11] = A[8]*B[3] + A[9]*B[7] + A[10]*B[11] + A[11]*B[15];

    res[12] = A[12]*B[0] + A[13]*B[4] + A[14]*B[8] + A[15]*B[12];
    res[13] = A[12]*B[1] + A[13]*B[5] + A[14]*B[9] + A[15]*B[13];
    res[14] = A[12]*B[2] + A[13]*B[6] + A[14]*B[10] + A[15]*B[14];
    res[15] = A[12]*B[3] + A[13]*B[7] + A[14]*B[11] + A[15]*B[15];
}


inline void gamma_times_psi(t_complex* res, const t_complex* psi, const int mu)
{
    switch (mu) {
        case 0:
            res[0] = psi[2];
            res[1] = psi[3];
            res[2] = psi[0];
            res[3] = psi[1];
            break;
        case 1:
            res[0] = {-psi[3].imag(), psi[3].real()};
            res[1] = {-psi[2].imag(), psi[2].real()};
            res[2] = {psi[1].imag(), -psi[1].real()};
            res[3] = {psi[0].imag(), -psi[0].real()};
            break;
        case 2:
            res[0] = psi[3];
            res[1] = -psi[2];
            res[2] = -psi[1];
            res[3] = psi[0];
            break;
        case 3:
            res[0] = {-psi[2].imag(), psi[2].real()};
            res[1] = {psi[3].imag(), -psi[3].real()};
            res[2] = {psi[0].imag(), -psi[0].real()};
            res[3] = {-psi[1].imag(), psi[1].real()};
            break;
    }
}

inline void sigma_times_psi(t_complex* res, const t_complex* psi, const int mu, const int nu)
{
    switch (mu) {
        case 0:
            switch (nu) {
                case 1:
                    res[0] = psi[1];
                    res[1] = psi[0];
                    res[2] = -psi[3];
                    res[3] = -psi[2];
                    break;
                case 2:
                    res[0] = {psi[1].imag(), -psi[1].real()};
                    res[1] = {-psi[0].imag(), psi[0].real()};
                    res[2] = {-psi[3].imag(), psi[3].real()};
                    res[3] = {psi[2].imag(), -psi[2].real()};
                    break;
                case 3:
                    res[0] = psi[0];
                    res[1] = -psi[1];
                    res[2] = -psi[2];
                    res[3] = psi[3];
                    break;
            }
            break;
        case 1:
            switch (nu) {
                case 2:
                    res[0] = -psi[0];
                    res[1] = psi[1];
                    res[2] = -psi[2];
                    res[3] = psi[3];
                    break;
                case 3:
                    res[0] = {psi[1].imag(), -psi[1].real()};
                    res[1] = {-psi[0].imag(), psi[0].real()};
                    res[2] = {psi[3].imag(), -psi[3].real()};
                    res[3] = {-psi[2].imag(), psi[2].real()};
                    break;
            }
            break;
        case 2:
            switch (nu) {
                case 3:
                    res[0] = -psi[1];
                    res[1] = -psi[0];
                    res[2] = -psi[3];
                    res[3] = -psi[2];
                    break;
            }
            break;
    }
}

inline void i_sigma_times_psi(t_complex* res, const t_complex* psi, const int mu, const int nu)
{
    switch (mu) {
        case 0:
            switch (nu) {
                case 1:
                    res[0] = {-psi[1].imag(), psi[1].real()};
                    res[1] = {-psi[0].imag(), psi[0].real()};
                    res[2] = {psi[3].imag(), -psi[3].real()};
                    res[3] = {psi[2].imag(), -psi[2].real()};
                    break;
                case 2:
                    res[0] = psi[1];
                    res[1] = -psi[0];
                    res[2] = -psi[3];
                    res[3] = psi[2];
                    break;
                case 3:
                    res[0] = {-psi[0].imag(), psi[0].real()};
                    res[1] = {psi[1].imag(), -psi[1].real()};
                    res[2] = {psi[2].imag(), -psi[2].real()};
                    res[3] = {-psi[3].imag(), psi[3].real()};
                    break;
            }
            break;
        case 1:
            switch (nu) {
                case 2:
                    res[0] = {psi[0].imag(), -psi[0].real()};
                    res[1] = {-psi[1].imag(), psi[1].real()};
                    res[2] = {psi[2].imag(), -psi[2].real()};
                    res[3] = {-psi[3].imag(), psi[3].real()};
                    break;
                case 3:
                    res[0] = psi[1];
                    res[1] = -psi[0];
                    res[2] = psi[3];
                    res[3] = -psi[2];
                    break;
            }
            break;
        case 2:
            switch (nu) {
                case 3:
                    res[0] = {psi[1].imag(), -psi[1].real()};
                    res[1] = {psi[0].imag(), -psi[0].real()};
                    res[2] = {psi[3].imag(), -psi[3].real()};
                    res[3] = {psi[2].imag(), -psi[2].real()};
                    break;
            }
            break;
    }
}

inline void projection_plus_times_psi(t_complex* res, const t_complex* psi) 
{
    res[0] = 0.5 * (psi[0] + psi[2]);
    res[1] = 0.5 * (psi[1] + psi[3]);
    res[2] = 0.5 * (psi[2] + psi[0]);
    res[3] = 0.5 * (psi[3] + psi[1]);
}

inline void projection_minus_times_psi(t_complex* res, const t_complex* psi) 
{
    res[0] = 0.5 * (psi[0] - psi[2]);
    res[1] = 0.5 * (psi[1] - psi[3]);
    res[2] = 0.5 * (psi[2] - psi[0]);
    res[3] = 0.5 * (psi[3] - psi[1]);
}



#endif // _SPINOR_ALGEBRA_HPP
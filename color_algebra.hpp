#ifndef _COLOR_ALGEBRA_HPP
#define _COLOR_ALGEBRA_HPP

using namespace std;

#include <complex>

typedef complex<double>       t_complex;

constexpr int NCOLOR = 3;
constexpr int NCOLOR2 = NCOLOR*NCOLOR;

inline void c3x3_set_unity(t_complex* A) 
{
    A[0] = t_complex(1.0, 0.0);
    A[1] = t_complex(0.0, 0.0);
    A[2] = t_complex(0.0, 0.0);

    A[3] = t_complex(0.0, 0.0);
    A[4] = t_complex(1.0, 0.0);
    A[5] = t_complex(0.0, 0.0);

    A[6] = t_complex(0.0, 0.0);
    A[7] = t_complex(0.0, 0.0);
    A[8] = t_complex(1.0, 0.0);
}


inline void c3x3_set_zero(t_complex* A) 
{
    A[0] = t_complex(0.0, 0.0);
    A[1] = t_complex(0.0, 0.0);
    A[2] = t_complex(0.0, 0.0);

    A[3] = t_complex(0.0, 0.0);
    A[4] = t_complex(0.0, 0.0);
    A[5] = t_complex(0.0, 0.0);

    A[6] = t_complex(0.0, 0.0);
    A[7] = t_complex(0.0, 0.0);
    A[8] = t_complex(0.0, 0.0);
}

inline void c3x3_setter(t_complex* A, const t_complex* B)
{
    for (int i = 0; i < NCOLOR*NCOLOR; i++) {
        A[i] = B[i];
    }
}

inline void c3x3_times_c3x3(t_complex* res, const t_complex* A, const t_complex* B) 
{
    res[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
    res[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
    res[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8];

    res[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
    res[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
    res[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8];

    res[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6];
    res[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7];
    res[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
}

inline t_complex tr_c3x3_times_c3x3(const t_complex* A, const t_complex* B) 
{
    return A[0]*B[0] + A[1]*B[3] + A[2]*B[6] +
           A[3]*B[1] + A[4]*B[4] + A[5]*B[7] + 
           A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
}

inline t_complex tr_c3x3_times_c3x3_conj(const t_complex* A, const t_complex* B) 
{
    return A[0]*conj(B[0]) + A[1]*conj(B[1]) + A[2]*conj(B[2]) + 
           A[3]*conj(B[3]) + A[4]*conj(B[4]) + A[5]*conj(B[5]) + 
           A[6]*conj(B[6]) + A[7]*conj(B[7]) + A[8]*conj(B[8]);
}

inline void c3x3_times_c3x3_conj(t_complex* res, const t_complex* A, const t_complex* B) 
{
    res[0] = A[0]*conj(B[0]) + A[1]*conj(B[1]) + A[2]*conj(B[2]);
    res[1] = A[0]*conj(B[3]) + A[1]*conj(B[4]) + A[2]*conj(B[5]);
    res[2] = A[0]*conj(B[6]) + A[1]*conj(B[7]) + A[2]*conj(B[8]);

    res[3] = A[3]*conj(B[0]) + A[4]*conj(B[1]) + A[5]*conj(B[2]);
    res[4] = A[3]*conj(B[3]) + A[4]*conj(B[4]) + A[5]*conj(B[5]);
    res[5] = A[3]*conj(B[6]) + A[4]*conj(B[7]) + A[5]*conj(B[8]);

    res[6] = A[6]*conj(B[0]) + A[7]*conj(B[1]) + A[8]*conj(B[2]);
    res[7] = A[6]*conj(B[3]) + A[7]*conj(B[4]) + A[8]*conj(B[5]);
    res[8] = A[6]*conj(B[6]) + A[7]*conj(B[7]) + A[8]*conj(B[8]);
}

inline void c3x3_conj_times_c3x3(t_complex* res, const t_complex* A, const t_complex* B) 
{
    res[0] = conj(A[0])*B[0] + conj(A[3])*B[3] + conj(A[6])*B[6];
    res[1] = conj(A[0])*B[1] + conj(A[3])*B[4] + conj(A[6])*B[7];
    res[2] = conj(A[0])*B[2] + conj(A[3])*B[5] + conj(A[6])*B[8];

    res[3] = conj(A[1])*B[0] + conj(A[4])*B[3] + conj(A[7])*B[6];
    res[4] = conj(A[1])*B[1] + conj(A[4])*B[4] + conj(A[7])*B[7];
    res[5] = conj(A[1])*B[2] + conj(A[4])*B[5] + conj(A[7])*B[8];

    res[6] = conj(A[2])*B[0] + conj(A[5])*B[3] + conj(A[8])*B[6];
    res[7] = conj(A[2])*B[1] + conj(A[5])*B[4] + conj(A[8])*B[7];
    res[8] = conj(A[2])*B[2] + conj(A[5])*B[5] + conj(A[8])*B[8];
}

inline void c3x3_conj_times_c3x3_conj(t_complex* res, const t_complex* A, const t_complex* B) 
{
    res[0] = conj(A[0])*conj(B[0]) + conj(A[3])*conj(B[1]) + conj(A[6])*conj(B[2]);
    res[1] = conj(A[0])*conj(B[3]) + conj(A[3])*conj(B[4]) + conj(A[6])*conj(B[5]);
    res[2] = conj(A[0])*conj(B[6]) + conj(A[3])*conj(B[7]) + conj(A[6])*conj(B[8]);

    res[3] = conj(A[1])*conj(B[0]) + conj(A[4])*conj(B[1]) + conj(A[7])*conj(B[2]);
    res[4] = conj(A[1])*conj(B[3]) + conj(A[4])*conj(B[4]) + conj(A[7])*conj(B[5]);
    res[5] = conj(A[1])*conj(B[6]) + conj(A[4])*conj(B[7]) + conj(A[7])*conj(B[8]);

    res[6] = conj(A[2])*conj(B[0]) + conj(A[5])*conj(B[1]) + conj(A[8])*conj(B[2]);
    res[7] = conj(A[2])*conj(B[3]) + conj(A[5])*conj(B[4]) + conj(A[8])*conj(B[5]);
    res[8] = conj(A[2])*conj(B[6]) + conj(A[5])*conj(B[7]) + conj(A[8])*conj(B[8]);
}

inline t_complex tr_c3x3(const t_complex* A) 
{
    return A[0] + A[4] + A[8];
}

inline void c3x3_times_c3(t_complex* res, const t_complex* A, const t_complex* q, const int stride) 
{
    res[0*stride] = A[0]*q[0*stride] + A[1]*q[1*stride] + A[2]*q[2*stride];
    res[1*stride] = A[3]*q[0*stride] + A[4]*q[1*stride] + A[5]*q[2*stride];
    res[2*stride] = A[6]*q[0*stride] + A[7]*q[1*stride] + A[8]*q[2*stride];
}

inline void c3x3_conj(t_complex* res, const t_complex* A) 
{
    res[0] = conj(A[0]);
    res[1] = conj(A[3]);
    res[2] = conj(A[6]);

    res[3] = conj(A[1]);
    res[4] = conj(A[4]);
    res[5] = conj(A[7]);

    res[6] = conj(A[2]);
    res[7] = conj(A[5]);
    res[8] = conj(A[8]);
}

inline double norm_diff(t_complex* a, t_complex* b, int n)
{
    double res = 0.0;
    for(int i=0; i<n; i++)
        res += norm(a[i] - b[i]);
    return sqrt(res);
}

inline double c3x3_unitarity_norm(const t_complex* u) 
{
    t_complex uud[9], id[9];
    c3x3_set_unity(id);
    c3x3_conj_times_c3x3(uud, u, u);
    return norm_diff(uud, id, 9);
}

inline t_complex c3x3_det(const t_complex* u) 
{
    return u[0]*(u[4]*u[8] - u[5]*u[7]) - 
           u[1]*(u[3]*u[8] - u[5]*u[6]) + 
           u[2]*(u[3]*u[7] - u[4]*u[6]);
}

#endif // _COLOR_ALGEBRA_HPP

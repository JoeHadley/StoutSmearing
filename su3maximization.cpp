#include "su3maximization.hpp"

constexpr double pi = numbers::pi;

double asin_branches(double x, int m){ return (m%2==0? asin(x) : pi - asin(x) ); }

double FN(double xi, double phi, double* sigma, int m0, int m1)
{
    return xi + sigma[2]*sin(phi + asin_branches(xi/sigma[0], m0) + asin_branches(xi/sigma[1], m1));
}

int MaximizeSU3Overlap(t_complex* S, t_complex* U, int max_iter, double precision)
{
    t_complex tmp[NCOLOR2], A[NCOLOR];
    std::copy(S, S + NCOLOR2, tmp);

    double sigma[NCOLOR];
    t_complex w[NCOLOR2], vt[NCOLOR2];
    double superb[NCOLOR];

    //lapack_int info = LAPACKE_zgesvd(LAPACK_ROW_MAJOR, 'A', 'A', NCOLOR, NCOLOR, tmp, NCOLOR, sigma, w, NCOLOR, vt, NCOLOR, superb);
    lapack_int info = LAPACKE_zgesvd(
    LAPACK_ROW_MAJOR,
    'A', 'A',
    NCOLOR, NCOLOR,
    reinterpret_cast<lapack_complex_double*>(tmp),  // <--- cast
    NCOLOR,
    sigma,
    reinterpret_cast<lapack_complex_double*>(w),    // cast if complex
    NCOLOR,
    reinterpret_cast<lapack_complex_double*>(vt),   // cast if complex
    NCOLOR,
    superb
    );
    
    if (info > 0) { return -info; } //SVD failed
    
    c3x3_times_c3x3(tmp, w, vt);
    t_complex wvt_det = c3x3_det(tmp);
    double phi = std::arg(wvt_det);
    double max_overlap = -1e10; double alpha_optimal[NCOLOR];  int m0optimal = 0, m1optimal = 0;

    int actual_max_iter = 0;
    
    for(int im=0; im<3; im++)
    {
        int m0 = im%2; int m1 = (im/2)%2; //Scanning over several branches of asin
        //Finding the root of FN(xi) = 0 using the bisection method
        double xiL = -sigma[1], xiR = sigma[1], xi = 0.0; int iter = 0;
        double fL = FN(xiL, phi, sigma, m0, m1), fR = FN(xiR, phi, sigma, m0, m1); if(fL*fR > 0.0) return -523; //No root in this interval

        do
        {
            xi = 0.5*(xiL + xiR);
            double f = FN(xi, phi, sigma, m0, m1);
            if(f*fL < 0.0) { xiR = xi; fR = f; } else { xiL = xi; fL = f; };
            iter++;
        }
        while(xiR - xiL > precision && iter < max_iter);
        if(abs(xi)>sigma[2]) return -524; //Root is out of range, something went wrong

        actual_max_iter = std::max(actual_max_iter, iter);

        double alpha[3];
        alpha[0] = asin_branches(xi/sigma[0], m0);
        alpha[1] = asin_branches(xi/sigma[1], m1);
        double alpha_tot = alpha[0] + alpha[1] + phi;
        double alpha2A = asin_branches(xi/sigma[2], 0);
        double alpha2B = asin_branches(xi/sigma[2], 1);
        //Which of the two branches is closer to giving alpha0 + alpha1 + alpha2 +phi = 2*pi*m, m is Integer?
        alpha[2] = (abs(sin(0.5*(alpha_tot + alpha2A))) < abs(sin(0.5*(alpha_tot + alpha2B))) ? alpha2A : alpha2B);
        //Calculating the overlap = re(Tr(S^+*U)) for this set of alphas
        double overlap = sigma[0]*cos(alpha[0]) + sigma[1]*cos(alpha[1]) + sigma[2]*cos(alpha[2]);
        if(overlap > max_overlap) { max_overlap = overlap; std::copy(alpha, alpha + NCOLOR, alpha_optimal); m0optimal = m0; m1optimal = m1; };
    };

    if(m0optimal != 0 && m1optimal != 0) return -525; //Something went wrong, we should have found the maximum on the first branch

    //Finally, reconstruct U = W * D * V^+ where D = diag(exp(i*alpha_optimal))
    for(int i=0; i<NCOLOR; i++)  A[i] = std::exp(t_complex(0.0, alpha_optimal[i]));

    for(int i=0; i<NCOLOR; i++)
        for(int j=0; j<NCOLOR; j++)
            tmp[i*NCOLOR + j] = w[i*NCOLOR + j]*A[j];

    c3x3_times_c3x3(U, tmp, vt);

    return actual_max_iter;
}

//This was a part of MaximizeSU3Overlap function earlier for testing purposes (to test SVD accuracy)
//    // Multiply singular values into w (w = w * sigma)
//     for (int i = 0; i < NCOLOR; ++i)
//         for (int j = 0; j < NCOLOR; ++j)
//             tmp[i*NCOLOR + j] = w[i*NCOLOR + j]*sigma[j];

//     // Reconstruct: recon = w * vt using color_algebra functions
//     c3x3_times_c3x3(tmp1, tmp, vt);
//     double err = norm_diff(tmp1, S, NCOLOR2);
//     stat[0] = err;
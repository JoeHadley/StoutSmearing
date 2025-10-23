#include "linalg.hpp"

ranlux48 rng_engine;
uniform_real_distribution<double> rng_uniform_dist;
normal_distribution<double> rng_normal_dist{0.0, 1.0};
uniform_int_distribution<int> z4_dist(0, 3);

void init_rng(int rng_engine_seed)
{
//Initialize the random number generator
	if(rng_engine_seed==0)
	{
		auto t1 = chrono::high_resolution_clock::now();
		int nanos = chrono::duration_cast<chrono::nanoseconds>(t1.time_since_epoch()).count();
		rng_engine.seed(nanos);
	}
	else
		rng_engine.seed(rng_engine_seed);
}

t_complex rand_complex()
{
    double re = rng_uniform_dist(rng_engine);
    double im = rng_uniform_dist(rng_engine);
    return t_complex(re, im);
}

void rand_complex(t_complex* z, int n)
{
    for(int i=0; i<n; i++)
        z[i] = rand_complex();
}

void rand_complex_z4(t_complex* z, int n)
{
    for(int i=0; i<n; i++)
    {
        int r = z4_dist(rng_engine);
        switch(r)
        {
            case 0: z[i] = t_complex( 1.0,  0.0); break;
            case 1: z[i] = t_complex( 0.0,  1.0); break;
            case 2: z[i] = t_complex(-1.0,  0.0); break;
            case 3: z[i] = t_complex( 0.0, -1.0); break;
        }
    }
}

void split_array(const t_complex* arr, int n, vector<t_complex>& pos, vector<t_complex>& neg, vector<t_complex>& zero)
{
    double tol = 1e-8;

    pos.clear();
    neg.clear();
    zero.clear();

    for (int i = 0; i < n; i++) {
        double val = arr[i].real();
        if (abs(val) <= tol) {
            zero.push_back(val);
        } else if (val > tol) {
            pos.push_back(val);
        } else {
            neg.push_back(val);
        }
    }

    sort(pos.begin(), pos.end(), [](const t_complex& a, const t_complex& b) {
        return a.real() < b.real();}
    );
    sort(neg.begin(), neg.end(), [](const t_complex& a, const t_complex& b) {
        return a.real() > b.real();}
    );
}

bool match_arrays(vector<t_complex> zero1, vector<t_complex> zero2, vector<t_complex> pos1, vector<t_complex> pos2, vector<t_complex> neg1, vector<t_complex> neg2)
{
    double tol = 1e-8;
    bool matching = true;

    if ((pos2.size() > pos1.size()) || (neg2.size() > neg1.size()) || (zero2.size() > zero1.size())) {
        matching = false;
        return matching;
    }
    if ((zero1.size() != zero2.size()) && ((pos2.size() > 0) || (neg2.size() > 0))) {
        matching = false;
        return matching;
    }
    for (int i = 0; i < pos2.size(); i++) {
        if ((pos1[i].real() - pos2[i].real()) > tol) {
            matching = false;
            return matching;
        }
    }
    for (int i = 0; i < neg2.size(); i++) {
        if ((neg1[i].real() - neg2[i].real() > tol)) {
            matching = false;
            return matching;
        }
    }
    return matching;
}
//Note that A will be overwritten!
int diagonalize(t_complex* A, t_complex* evals, t_complex* revecs, t_complex* levecs, int n)
{
    //int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'V', 'V', n, A, n, evals, levecs, n, revecs, n);
    lapack_complex_double* A_ptr    = reinterpret_cast<lapack_complex_double*>(A);
    lapack_complex_double* evals_ptr= reinterpret_cast<lapack_complex_double*>(evals);
    lapack_complex_double* levecs_ptr = reinterpret_cast<lapack_complex_double*>(levecs);
    lapack_complex_double* revecs_ptr = reinterpret_cast<lapack_complex_double*>(revecs);

    int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'V', 'V', n,
                          A_ptr, n,
                          evals_ptr,
                          levecs_ptr, n,
                          revecs_ptr, n);
//
    
    if(info != 0) return info;
    conjugate_transpose(levecs, n);

    //Restoring the proper normalisation of revecs and levecs
	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
		{
			A[i*n + j] = 0.0;
			for(int k=0; k<n; k++)
				A[i*n + j] += levecs[i*n + k] * revecs[k*n + j];
		};

    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            revecs[i*n + j] /= A[j*n + j];
	return 0;
}

void conjugate_transpose(t_complex* A, int n)
{
	for(int i=0; i<n; i++)
	{
		A[i*n + i] = conj(A[i*n + i]);
		for(int j=i+1; j<n; j++)
		{
			t_complex tmp = conj(A[i*n + j]);
			A[i*n + j] = conj(A[j*n + i]);
			A[j*n + i] = tmp;
		}
	}
}

void  A_pluseq_bB(t_complex* A, const t_complex b, const t_complex *B, const uint n)
{
	cblas_zaxpy(n,&b,B,1,A,1);
}

void   aA_plus_bB(t_complex* out, const t_complex a, const t_complex* A, const t_complex b, const t_complex* B, const uint n)
{
    #pragma omp parallel for
    for(uint i=0; i<n; i++)
        out[i] = a*A[i] + b*B[i];
};

t_complex dot_product(const t_complex* x, const t_complex* y, const int n) {
    t_complex res = 0;
    for (int i = 0; i < n; ++i) {
        res += conj(x[i]) * y[i];
    }
    return res;
}

double complex_norm(const t_complex* x, const int n) {
    return sqrt(real(dot_product(x, x, n)));
}

// Givens rotation
void apply_givens(t_complex& h1, t_complex& h2,
                 t_complex& c, t_complex& s) {
    if (abs(h2) < 1e-20) {
        c = 1.0;
        s = 0.0;
    }
    else if (abs(h2) > abs(h1)) {
        t_complex temp = h1 / h2;
        s = 1.0 / sqrt(1.0 + norm(temp));
        c = temp * s;
    } else {
        t_complex temp = h2 / h1;
        c = 1.0 / sqrt(1.0 + norm(temp));
        s = temp * c;
    }
    t_complex temp = conj(c) * h1 + conj(s) * h2;
    h2 = -s * h1 + c * h2;
    h1 = temp;
}

void GMRES(t_complex* res,
           const int n, const int maxiter, const int restart, const double tol,
           void (*Alpha)(t_complex*, const t_complex*),
           void (*Beta)(t_complex*, const t_complex*),
           const t_complex* input)
{
    // TO DO: Check against wiki
    // GMRES linear solver for x: A*x = b = B*y
    // See [https://en.wikipedia.org/wiki/Generalized_minimal_residual_method]

    // Allocate work arrays
    t_complex* r = new t_complex[n];  // Residual (error) vector: r = b - A*x
    t_complex* w = new t_complex[n];  // Work vector: temporary storage
    t_complex** V = new t_complex*[restart + 1];  // Krylov basis vectors
    for (int i = 0; i < restart + 1; ++i) {
        V[i] = new t_complex[n];  // Arnoldi basis vectors
    }
    t_complex** H = new t_complex*[restart + 1];  // Upper Hessenberg matrix H
    for (int i = 0; i < restart + 1; ++i) {
        H[i] = new t_complex[restart];
    }
    t_complex* cs = new t_complex[restart];  // Givens rotation cosines
    t_complex* sn = new t_complex[restart];  // Givens rotation sines
    t_complex* e1 = new t_complex[restart + 1];  // RHS vector in Krylov space
    t_complex* y = new t_complex[restart];  // Coefficients of solution in Krylov basis

    // Initialization
    for (int i = 0; i < n; ++i) {
        res[i] = 0;
    }
    Beta(r, input);  // Initial residual: r = b - A*x = b = B*y
    double bnorm = complex_norm(r, n);  // Initial bnorm = rnorm
    double error = 1.0;  // Initial error: rnorm/bnorm = 1.0
    int iter = 0;

    // Main loop
    while ((iter < maxiter) && (error > tol)) {
         // First Krylov basis vector
        double rnorm = complex_norm(r, n);
        t_complex inv_rnorm = 1.0 / rnorm;
        for (int j = 0; j < n; ++j) {
            V[0][j] = r[j] * inv_rnorm;
        }

        // Initialize RHS for least squares system (Arnoldi coefficients)
        for (int i = 0; i < restart + 1; ++i) {
            e1[i] = 0;
        }
        e1[0] = rnorm;

        int j;
        for (j = 0; j < restart && iter < maxiter; ++j, ++iter) {
            // Arnoldi function and Modified Gram-Schmidt
            Alpha(w, V[j]);

            for (int i = 0; i < j + 1; ++i) {
                H[i][j] = dot_product(V[i], w, n);
                for (int k = 0; k < n; ++k) w[k] -= H[i][j] * V[i][k];
            }
            H[j+1][j] = complex_norm(w, n);

            // If nearly zero, Krylov subspace has stagnated
            if (abs(H[j+1][j]) < 1e-15) {                      
                j++;
                break;
            }

            // Normalize Krylov vector
            t_complex inv_h = 1.0 / H[j+1][j];                 
            for (int k = 0; k < n; ++k) {
                V[j+1][k] = w[k] * inv_h;
            }

            // Apply previous Givens rotations to new column of H
            for (int i = 0; i < j; ++i) {
                t_complex temp = conj(cs[i]) * H[i][j] + conj(sn[i]) * H[i+1][j];
                H[i+1][j] = -sn[i] * H[i][j] + cs[i] * H[i+1][j];
                H[i][j] = temp;
            }

            // Apply next Givens rotation
            apply_givens(H[j][j], H[j+1][j], cs[j], sn[j]);

            // Update residual vectors
            e1[j+1] = -sn[j] * e1[j];
            e1[j] = conj(cs[j]) * e1[j];
            error = abs(e1[j+1]) / bnorm;
            
            // Early exit if convergence reached
            if (error < tol) {
                j++;
                break;
            }
        }

        // Solve least-squares problem H*y = e1
        for (int i = j-1; i >= 0; --i) {
            y[i] = e1[i];
            for (int k = i+1; k < j; ++k) {
                y[i] -= H[i][k] * y[k];
            }
            y[i] /= H[i][i];
        }

        // Update result
        for (int i = 0; i < j; ++i)
            for (int k = 0; k < n; ++k)
                res[k] += V[i][k] * y[i];

        // Compute new residual explicitly if restart needed
        if (error > tol) {
            Beta(r, input);
            Alpha(w, res);
            for (int i = 0; i < n; ++i) {
                r[i] -= w[i];
            }
        }
    }


    if (error > tol) {
        cerr << "GMRES did not converge with error: " << error << "\n";
    }

    // Free memory
    delete[] r;
    delete[] w;
    delete[] cs;
    delete[] sn;
    delete[] e1;
    delete[] y;
    for (int i = 0; i <= restart; ++i) {
        delete[] V[i];
    }
    delete[] V;
    for (int i = 0; i <= restart; ++i) {
        delete[] H[i];
    }
    delete[] H;
}


t_complex element_product(const t_complex* A, const int n)
{
    t_complex res = {1.0, 0.0};

    for (int i = 0; i < n; i++)
    {
        res *= A[i];
    }

    return res;
}

t_complex element_log_sum(const t_complex* A, const int n)
{
    t_complex res = {0.0, 0.0};

    for (int i = 0; i < n; i++)
    {
        res += log(abs(A[i]));
    }

    return res;
}

void matrix_multiplication(t_complex* res, const t_complex* A, const t_complex* B, const int n)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res[i*n + j] = 0.0;
            for (int k = 0; k < n; k++) {
                res[i*n + j] += A[i*n + k] * B[k*n + j];
            }
        }
    }
}

void matrix_vector_mult(t_complex* res, const t_complex* A, const t_complex* v, const int n)
{
    for (int i = 0; i < n; i++) {
        res[i] = 0.0;
        for (int j = 0; j < n; j++) {
            res[i] += A[i*n + j] * v[j];
        }
    }
}

void rescale(t_complex* A, t_complex a, uint n)
{
	cblas_zscal(n, &a, A, 1);
}

double   norm(t_complex* psi, uint n)
{
	return cblas_dznrm2(n, psi, 1);
}

double   norm_diff(t_complex* psi1, t_complex* psi2, uint n)
{
    double res = 0.0;
    #pragma omp parallel for
    for(int i=0; i<n; i++)
        res += std::norm(psi1[i] - psi2[i]);
    return res;
}
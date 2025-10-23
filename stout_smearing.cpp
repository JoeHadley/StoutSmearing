#include "stout_smearing.hpp"

// Computes the staples for the link at position x in direction mu
void get_staples(t_complex* C,
                Lattice &lat,
                GaugeLikeField &gf,
                int x, int mu,
                double* rho){

    for (int i = 0; i < NCOLOR2; i++)
        C[i] = t_complex(0.0, 0.0);

    t_complex tmp_f[NCOLOR2], fwd_staple[NCOLOR2];
    t_complex tmp_b[NCOLOR2], bwd_staple[NCOLOR2];

    double rho_val;

    // For each direction nu, 2 staples
    for(uint nu=0; nu<4; nu++){
        if(nu==mu) continue; // Skip the principle direction
        

        // Set up isotropic rho. Assume rho[0] is time, rho[1], [2], [3] are all space
        if (mu ==0 || nu==0 ){
            rho_val = rho[0];
        }
        else {
            rho_val = rho[1];
        }

        // Forward staple 
        t_complex* F1 = gf.Link(x, nu);
        t_complex* F2 = gf.Link(lat.shift_fwd(x, nu), mu);
        t_complex* F3 = gf.Link(lat.shift_fwd(x, mu), nu); // To conjugate

        c3x3_times_c3x3(tmp_f, F1, F2);
        c3x3_times_c3x3_conj(fwd_staple, tmp_f, F3);

        // Backward staple
        t_complex* B1 = gf.Link(lat.shift_bwd(x,nu), nu); // To conjugate
        t_complex* B2 = gf.Link(lat.shift_bwd(x,nu), mu);
        t_complex* B3 = gf.Link(lat.shift_fwd(lat.shift_bwd(x, nu), mu), nu);

        c3x3_conj_times_c3x3(tmp_b, B1, B2);
        c3x3_times_c3x3(bwd_staple, tmp_b, B3);


        // Accumulate staples
        for (int i=0; i<NCOLOR2; i++){
            C[i] += rho_val*(fwd_staple[i] + bwd_staple[i]);
        };

    };


}


void project_to_su3(t_complex* res, t_complex* M ){

    t_complex Mdagger[NCOLOR2];
    t_complex diff[NCOLOR2];

    c3x3_conj(Mdagger, M);

    // diff = M - Mdagger
    for(int i=0; i<NCOLOR2; i++)
        diff[i] = .5*(M[i] - Mdagger[i]);
    
    t_complex tr = tr_c3x3(diff)/3.0;
    

    t_complex trace_Identity[NCOLOR2] = {
        tr, t_complex(0.0, 0.0), t_complex(0.0, 0.0),
        t_complex(0.0, 0.0), tr, t_complex(0.0, 0.0),
        t_complex(0.0, 0.0), t_complex(0.0, 0.0), tr
    };

    for (int i=0; i<NCOLOR2; i++)

        res[i] = diff[i] - trace_Identity[i];

};




// signature: projection_field should be a caller-allocated t_complex array of length 4*gf.vol*NCOLOR2
// rho is pointer to two values {rho_time, rho_space} (your isotropic convention)
void compute_antihermitian_matrix(t_complex* projection_field,
                                Lattice &lat,
                                GaugeLikeField &gf,
                                double* rho)
{
    // temporaries (reuse inside loops to avoid allocations)
    t_complex C[NCOLOR2];    // staple C_mu(x)
    t_complex Omega[NCOLOR2]; // Omega = C * U^\dagger
    t_complex P[NCOLOR2];    // projection (anti-Hermitian traceless result)

    // sanity: ensure projection_field not null
    if (!projection_field) return;

    // Loop over lattice sites and directions
    for (int x = 0; x < gf.vol; ++x) {
        for (uint mu = 0; mu < 4; ++mu) {
            // 1) compute staple C_mu(x) into C
            get_staples(C, lat, gf, x, mu, rho);

            // 2) compute Omega = C * U_mu(x)^\dagger
            //    We assume you have c3x3_times_c3x3_conj(dest, A, B) which performs dest = A * B^\dagger
            c3x3_times_c3x3_conj(Omega, C, gf.Link(x, mu));

            // 3) project Omega to traceless anti-Hermitian (su3 element) into P
            //    project_to_su3(res, M) should write 9 elements to res
            project_to_su3(P, Omega);

            // 4) store P into the projection_field array with same layout as Links:
            //    index = (x*4 + mu)*NCOLOR2 + i
            int base = (x * 4 + mu) * NCOLOR2;
            for (int i = 0; i < NCOLOR2; ++i)
                projection_field[base + i] = P[i];
        }
    }
}

void exp_iQ(t_complex* res,  t_complex* Q, t_complex* evals_out,t_complex* eigenvectors_out ) {
    const int N = 3;
    int info;

    // Step 1: Build iQ (since Q is anti-Hermitian, iQ is Hermitian)
    t_complex iQ[NCOLOR2];
    for (int i = 0; i < N*N; ++i)
        iQ[i] = t_complex(0.0, 1.0) * Q[i];

    // Step 2: Compute eigen-decomposition: iQ = V * diag(λ) * V†
    double eigvals[N];             // real eigenvalues
    t_complex V[NCOLOR2];          // will hold eigenvectors

    // Copy iQ because LAPACKE overwrites input
    for (int i = 0; i < N*N; ++i)
        V[i] = iQ[i];

    info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U',
                         N,
                         reinterpret_cast<lapack_complex_double*>(V),
                         N,
                         eigvals);

    if (info != 0) {
        std::cerr << "LAPACKE_zheev failed with info = " << info << std::endl;
        return;
    }

    // Step 3: Compute exp(iQ) = V * exp(i*λ) * V†
    t_complex expD[NCOLOR2] = {0};
    for (int i = 0; i < N; ++i)
        expD[i*N + i] = std::exp(t_complex(0.0, eigvals[i]));

    t_complex temp[NCOLOR2] = {0};
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k)
                temp[i*N + j] += V[i*N + k] * expD[k*N + k] * std::conj(V[j*N + k]);

    for (int i = 0; i < N*N; ++i)
        res[i] = temp[i];

    // Optional: output eigenvalues
    if (evals_out) {
        for (int i = 0; i < N; ++i)
            evals_out[i] = t_complex(eigvals[i], 0.0);
    }
    // Optional: output eigenvectors
    if (eigenvectors_out) {
        for (int i = 0; i < N*N; ++i)
            eigenvectors_out[i] = V[i];
    }
}



void smear_field( GaugeLikeField &output,
                  GaugeLikeField &gf,
                  Lattice &lat,
                  int n_steps,
                  double* rho) {

    // Allocate projection field
    std::vector<t_complex> projection_field(4 * gf.vol * NCOLOR2, t_complex(0.0, 0.0));

    // Temporary arrays for link updates
    t_complex exp_iQ_matrix[NCOLOR2];
    t_complex updated_link[NCOLOR2];

    // Start with input field
    GaugeLikeField current = gf;

    for (int step = 0; step < n_steps; ++step) {
        // 1) Compute anti-Hermitian traceless matrices
        compute_antihermitian_matrix(projection_field.data(), lat, current, rho);

        // 2) Update gauge links: U_mu(x) -> exp(i * Q_mu(x)) * U_mu(x)
        for (int x = 0; x < current.vol; ++x) {
            for (uint mu = 0; mu < 4; ++mu) {
                // Get Q_mu(x)
                t_complex* Q = projection_field.data() + (x * 4 + mu) * NCOLOR2;

                // Compute exp(iQ)
                exp_iQ(exp_iQ_matrix, Q);

                // Update link: U_mu(x) = exp(iQ) * U_mu(x)
                t_complex* U_link = current.Link(x, mu);
                c3x3_times_c3x3(updated_link, exp_iQ_matrix, U_link);

                // Store updated link back
                for (int i = 0; i < NCOLOR2; ++i)
                    U_link[i] = updated_link[i];
            }
        }
    }
        // Copy result to output
    for (int x = 0; x < output.vol; ++x)
        for (uint mu = 0; mu < 4; ++mu) {
            t_complex* U_src = current.Link(x, mu);
            t_complex* U_dst = output.Link(x, mu);
            for (int i = 0; i < NCOLOR2; ++i)
                U_dst[i] = U_src[i];
        }
}
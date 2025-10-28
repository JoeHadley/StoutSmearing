#include <iostream>
#include <fstream>
#include <complex>
#include <random>
#include <string>



using namespace std;

#include "color_algebra.hpp"
#include "lattice.hpp"
#include "gauge_field.hpp"
#include "stout_smearing.hpp"





void printMatrix(const t_complex* A, const string& name = "") 
{


  string title;
  if (!name.empty()) {
    title = "Matrix " + name + ":\n";
  } else {
    title = "Matrix:\n";
  }

  

  std::cout << title;
    
    for(int i = 0; i < NCOLOR; ++i){
        for(int j = 0; j < NCOLOR; ++j){
            std::cout << A[i*NCOLOR + j] << "  ";
        }
        std::cout << "\n";
    }
}

void randomlyFillMatrix(t_complex* A) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < NCOLOR2; ++i) {
        double realPart = dis(gen);
        double imagPart = dis(gen);
        A[i] = t_complex(realPart, imagPart);
    }
}






void exp_iQ(t_complex* U, const t_complex* Q) {
    constexpr int N = 3;
    t_complex A[N*N];
    double w[N];  // eigenvalues

    // Copy input (LAPACKE_zheev overwrites it)
    for (int i = 0; i < N*N; ++i)
        A[i] = Q[i];

    // LAPACKE expects lapack_complex_double*
    lapack_int info = LAPACKE_zheev(
        LAPACK_ROW_MAJOR, 'V', 'U', N,
        reinterpret_cast<lapack_complex_double*>(A),
        N, w
    );

    if (info != 0) {
        std::cerr << "Error in LAPACKE_zheev: info=" << info << std::endl;
        return;
    }

    // Compute exp(i * lambda)
    t_complex expD[N*N] = {};
    for (int i = 0; i < N; ++i)
        expD[i*N + i] = std::exp(t_complex(0, w[i])); // e^{i λ_i}

    // Reconstruct exp(iQ) = V * expD * V†
    t_complex tmp[N*N];
    c3x3_times_c3x3(tmp, A, expD);        // tmp = V * expD
    c3x3_times_c3x3_conj(U, tmp, A);      // U = tmp * V†
}



void exp_iQ2(t_complex* res, t_complex* Q ){

    // Use LAPACK to compute the matrix exponential of iQ
    // Q is anti-Hermitian and traceless, so iQ is Hermitian

    // Prepare data for LAPACK
    t_complex A[NCOLOR2];
    for (int i = 0; i < NCOLOR2; ++i) {
        A[i] = Q[i];
    }

    // Compute the matrix exponential using LAPACK
    int n = 3;
    int lda = 3;
    int info;

    // Create a copy of A to hold the result
    t_complex expA[NCOLOR2];

    // Use LAPACK's zgeev to compute eigenvalues and eigenvectors
    t_complex w[3]; // eigenvalues
    t_complex vl[NCOLOR2]; // left eigenvectors
    t_complex vr[NCOLOR2]; // right eigenvectors

    info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U',
                         N,
                         reinterpret_cast<lapack_complex_double*>(V),
                         N,
                         eigvals);

    if (info != 0) {
        // Handle error
        return;
    }

    // Construct exp(iQ) from eigenvalues and eigenvectors
    for (int i = 0; i < n; ++i) {
        t_complex lambda = w[i];
        t_complex exp_lambda = exp(lambda);
        for (int j = 0; j < n; ++j) {
            expA[i * n + j] = exp_lambda * vr[i * n + j];
        }
    }

    // Copy result to res
    for (int i = 0; i < NCOLOR2; ++i) {
        res[i] = expA[i];
    }
};










int main() {


  std::string cfg_file = "unit_LT8_LS4.bin";
  GaugeField U(cfg_file);

  std::cout << U << std::endl;


  t_complex A[NCOLOR2];
  t_complex B[NCOLOR2];
  t_complex C[NCOLOR2];


  randomlyFillMatrix(A);
  randomlyFillMatrix(B);

  printMatrix(A, "A");
  printMatrix(B, "B");



  // Multiply A * B
  c3x3_times_c3x3(C, A, B);

  // Print the result
  printMatrix(C, "C = A * B");







  GaugeField gf("unit_LT8_LS4.bin");

  std::cout << "Before randomization, mean plaquette: " << gf.MeanPlaquette() << std::endl;

  gf.RandomizeLinks();

  std::cout << "After randomization, mean plaquette: " << gf.MeanPlaquette() << std::endl;

  return 0;

  // Reconstruct B from the eigen-decomposition of exp(iB)
  // Vdagger * D * V


  
  //t_complex conj_link[9];
  //c3x3_conj(conj_link, first_link);

  //std::cout << "Conjugate of first link matrix:\n";
  //for (int i = 0; i < 3; ++i) {
  //    for (int j = 0; j < 3; ++j) {
  //        std::cout << conj_link[i*3 + j] << " ";
  //    }
  //    std::cout << "\n";
  //}
}










  // --- Reconstruct B from eigen-decomposition ---
  t_complex Lambda[NCOLOR2] = {0};
  for (int i = 0; i < NCOLOR; ++i)
      Lambda[i*NCOLOR + i] = evals[i]; // evals = real eigenvalues of iB

  // temp = V * diag(Lambda)
  //t_complex temp[NCOLOR2];
  //c3x3_times_c3x3(temp, eigenvectors, Lambda);

  // iB_reconstructed = V * diag(Lambda) * V†
  t_complex iB_reconstructed[NCOLOR2];
  
  VtimesDiagLtimesVdagger(iB_reconstructed, eigenvectors, Lambda);
  
  //c3x3_conj_times_c3x3(iB_reconstructed, eigenvectors, temp);

  // B_reconstructed = -i * iB_reconstructed
  t_complex B_reconstructed[NCOLOR2];
  for (int i = 0; i < NCOLOR2; ++i)
      B_reconstructed[i] = t_complex(0.0, -1.0) * iB_reconstructed[i];

  // Compute reconstruction error
  double recon_error = norm_diff(B, B_reconstructed, NCOLOR2);
  std::cout << "Reconstruction error ||B - (-i V diag(evals) V†)|| = "
            << recon_error << std::endl;

  return 0;
}
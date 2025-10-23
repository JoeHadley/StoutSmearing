#include <iostream>
#include <random>
#include "lattice.hpp"
#include "gauge_field.hpp"
#include "color_algebra.hpp"
#include "stout_smearing.hpp"

#include <string>
 



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

int main() {

  t_complex A[NCOLOR2];
  t_complex B[NCOLOR2];

  // Fill A with random complex numbers
  randomlyFillMatrix(A);
  printMatrix(A, "A");

  // Project A onto anti-Hermitian, traceless matrix B
  project_to_su3(B, A);
  printMatrix(B, "B = proj_SU3(A)");

  std::cout << "Trace(B) = " << tr_c3x3(B) << std::endl;

  // Check anti-Hermiticity: B + B† ≈ 0
  t_complex Bdagger[NCOLOR2];
  c3x3_conj(Bdagger, B);

  double antiHermNorm = 0.0;
  for (int i = 0; i < NCOLOR2; i++)
      antiHermNorm += norm(B[i] + Bdagger[i]);
  antiHermNorm = sqrt(antiHermNorm);
  std::cout << "||B + B†|| = " << antiHermNorm << std::endl;

  // --- Compute exp(iB) and eigen-decomposition ---
  t_complex expA[NCOLOR2];
  t_complex evals[NCOLOR];           // real eigenvalues of iB
  t_complex eigenvectors[NCOLOR2];   // eigenvectors of iB

  exp_iQ(expA, B, evals, eigenvectors);
  printMatrix(expA, "exp(iB)");

  // --- Sanity check: unitarity of exp(iB) ---
  double unitarity = 0.0;
  for (int i = 0; i < NCOLOR; ++i) {
      for (int j = 0; j < NCOLOR; ++j) {
          t_complex sum = 0.0;
          for (int k = 0; k < NCOLOR; ++k)
              sum += std::conj(expA[k*NCOLOR + i]) * expA[k*NCOLOR + j];
          if (i == j) sum -= t_complex(1.0, 0.0);
          unitarity += std::norm(sum);
      }
  }
  std::cout << "||U†U - I|| = " << sqrt(unitarity) << std::endl;

  // --- Check determinant ---
  t_complex det =
      expA[0]*(expA[4]*expA[8] - expA[5]*expA[7])
    - expA[1]*(expA[3]*expA[8] - expA[5]*expA[6])
    + expA[2]*(expA[3]*expA[7] - expA[4]*expA[6]);
  std::cout << "det(U) = " << det
            << "  |det| = " << std::abs(det) << std::endl;

  // --- Reconstruct B from eigen-decomposition ---
  t_complex Lambda[NCOLOR2] = {0};
  for (int i = 0; i < NCOLOR; ++i)
      Lambda[i*NCOLOR + i] = evals[i]; // evals = real eigenvalues of iB

  // temp = V * diag(Lambda)
  t_complex temp[NCOLOR2];
  c3x3_times_c3x3(temp, eigenvectors, Lambda);

  // iB_reconstructed = V * diag(Lambda) * V†
  t_complex iB_reconstructed[NCOLOR2];
  c3x3_conj_times_c3x3(iB_reconstructed, eigenvectors, temp);

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



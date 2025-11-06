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


void traceTest(const t_complex* A) {
    t_complex trace = tr_c3x3(A);
    std::cout << "Trace: " << trace << std::endl;
}   

void antiHermitianTest(const t_complex* A, string name = "B") {
    t_complex Adagger[NCOLOR2];
    c3x3_conj(Adagger, A);

    double sum = 0.0;
    for (int i = 0; i < NCOLOR2; i++)
        sum += std::norm(A[i] + Adagger[i]);  // <-- std::norm(), not variable name
    double antiHermNorm = std::sqrt(sum);

    string message;
    message = "||" + name + " + " + name + "†|| = ";
    std::cout << message << antiHermNorm << std::endl;
}


void unitarity_test(const t_complex* U, std::string name = "U") {
    double sum = 0.0;

    for (int i = 0; i < NCOLOR; ++i) {
        for (int j = 0; j < NCOLOR; ++j) {
            t_complex inner = 0.0;
            for (int k = 0; k < NCOLOR; ++k)
                inner += std::conj(U[k*NCOLOR + i]) * U[k*NCOLOR + j];
            if (i == j)
                inner -= t_complex(1.0, 0.0);  // subtract identity
            sum += std::norm(inner);
        }
    }

    double norm = std::sqrt(sum);
    std::cout << "||" << name << "†" << name << " - I|| = " << norm << std::endl;
}


void determinant_test(const t_complex* U, std::string name = "U") {
    t_complex det = U[0]*(U[4]*U[8] - U[5]*U[7])
                  - U[1]*(U[3]*U[8] - U[5]*U[6])
                  + U[2]*(U[3]*U[7] - U[4]*U[6]);
    std::cout << "det(" << name << ") = " << det
              << "   |det| = " << std::abs(det) << std::endl;
}

void reconstruction_test(t_complex* B, t_complex* eigenvectors, t_complex* evals, std::string name = "B") {
  // Make diagonal matrix of eigenvalues
  t_complex Lambda[NCOLOR2] = {0};
  for (int i = 0; i < NCOLOR; ++i)
      Lambda[i*NCOLOR + i] = evals[i]; // evals = real eigenvalues of iB


  t_complex iB_reconstructed[NCOLOR2];
  // Reconstruct iB: V * diag(evals) * V†  
  VtimesDiagLtimesVdagger(iB_reconstructed, eigenvectors, Lambda);

  t_complex B_reconstructed[NCOLOR2];
  for (int i = 0; i < NCOLOR2; ++i)
      B_reconstructed[i] = t_complex(0.0, -1.0) * iB_reconstructed[i];

  // Compute reconstruction error
  double recon_error = norm_diff(B, B_reconstructed, NCOLOR2);
  std::cout << "Reconstruction error ||B - (-i V diag(evals) V†)|| = "
            << recon_error << std::endl;
}


void run_tests(){
  t_complex A[NCOLOR2];
  t_complex B[NCOLOR2];

  // Fill A with random complex numbers
  randomlyFillMatrix(A);
  //printMatrix(A, "A");

  // Project A onto anti-Hermitian, traceless matrix B
  project_to_su3(B, A);
  //printMatrix(B, "B = proj_SU3(A)");

  // Verify B is traceless
  traceTest(B);

  // Verify B is anti-Hermitian
  antiHermitianTest(B,"B");



  // Compute C = exp(iB) and eigen-decomposition
  t_complex C[NCOLOR2];
  t_complex evals[NCOLOR];           // real eigenvalues of iB
  t_complex eigenvectors[NCOLOR2];   // eigenvectors of iB

  exp_iQ(C, B, evals, eigenvectors);
  printMatrix(C, "exp(iB)");

  // Verify C is unitary
  unitarity_test(C);


  // Verify det(C) = 1
  determinant_test(C);


  // Verify reconstruction of B from eigen-decomposition
  reconstruction_test(B, eigenvectors, evals, "B");

}

int main() {

  // load a gauge field
  GaugeField gf("Gen2_16x24n9820");

  std::cout << "Mean plaquette before smearing: " << gf.MeanPlaquette() << std::endl;
  // set smearing parameters
  int n_steps = 1;
  double rho[4] = {0.0, 0.14, 0.14, 0.14}; // isotropic smearing
  


  Lattice lat_copy(gf.LT, gf.LS);
  GaugeField smeared_gf(lat_copy);
  //GaugeField smeared_gf("unit_LT8_LS4.bin");

  //GaugeField smeared_gf(gf.vol);   // create a gauge field with the same lattice size
  for(int x=0; x<gf.vol; x++)
    for(uint mu=0; mu<4; mu++)
      std::copy(gf.Link(x, mu), gf.Link(x, mu)+NCOLOR2, smeared_gf.Link(x, mu));


  smear_field(smeared_gf, gf, n_steps, rho);
  std::cout << "Mean plaquette after smearing: " << smeared_gf.MeanPlaquette() << std::endl;


  return 0;
}



#ifndef _STOUT_SMEARING_HPP
#define _STOUT_SMEARING_HPP

using namespace std;

#include <complex>
#include "color_algebra.hpp"
#include <vector>
#include "lattice.hpp"
#include <lapacke.h>
#include "gauge_field.hpp"

struct GaugeLikeField {
    std::vector<t_complex> data;
    int vol;

    GaugeLikeField(int vol_) : vol(vol_) {
        data.resize(4 * vol * 9, t_complex(0.0, 0.0));
    }

    inline t_complex* Link(int x, uint mu) {
        return data.data() + (x * 4 + mu) * 9;
    }
    inline const t_complex* Link(int x, uint mu) const {
    return data.data() + (x * 4 + mu) * 9;
    }
};


void get_staples(t_complex* C,
                GaugeField &gf,
                int x, int mu,
                double* rho);

void project_to_su3(t_complex* res, t_complex* M );

void compute_antihermitian_matrix(t_complex* projection_field,
                                GaugeField &gf,
                                double* rho);

void VtimesDiagLtimesVdagger(t_complex* res, const t_complex* V, const t_complex* L);

void exp_iQ(t_complex* res,  t_complex* Q, t_complex* evals_out= nullptr,t_complex* eigenvectors_out= nullptr );

void smear_field( GaugeField &output_gf, GaugeField &input_gf, int n_steps, double* rho);
#endif // _STOUT_SMEARING_HPP
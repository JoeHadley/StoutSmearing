#ifndef _COLOR_SPINOR_HPP
#define _COLOR_SPINOR_HPP

#include <random>
#include <iostream>
#include <iomanip>
#include "color_algebra.hpp"

using namespace std;

constexpr int print_width = 18;

inline void print_c3x3(t_complex* A)
{
    for (int i = 0; i < NCOLOR; i++) {
        cout << "[";
        for (int j = 0; j < NCOLOR; j++)
            cout << std::setw(print_width) << A[NCOLOR*i + j];
        cout << "]" << endl;
    }
}

inline void print_s4x4(t_complex* A)
{
    for (int i = 0; i < NSPINOR; i++) {
        cout << "[";
        for (int j = 0; j < NSPINOR; j++)
            cout << std::setw(print_width) << A[NSPINOR*i + j];
        cout << "]" << endl;
    }
}

#endif // _COLOR_SPINOR_HPP

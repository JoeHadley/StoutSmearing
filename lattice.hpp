#ifndef _LATTICE_HPP
#define _LATTICE_HPP


#include <fstream>
#include <iostream>  // if you use std::cout or std::cerr

#include <boost/program_options.hpp>

#include "ansi_io.hpp"
#include "spinor_algebra.hpp"

namespace po = boost::program_options;

typedef unsigned int uint;

class Lattice 
{
private:
    void init()
    {
        size[0] = LT;
        size[1] = LS;
        size[2] = LS;
        size[3] = LS;

        vol   = LT*LS*LS*LS;
        vol3D =    LS*LS*LS;

        block_size[4] = 1;
        for(int mu=3; mu>=0; mu--)
            block_size[mu] = size[mu]*block_size[mu + 1];

        beta = (double)LT;
    }
public:
    int size[4];  // size = {LT, LS, LS, LS}
    int block_size[5];  // block_size = {LT*LS^3, LS^3, LS^2, LS, 1}
    int vol;
    int vol3D;
    int LT, LS;
    double beta;
    po::options_description lattice_params{"Lattice parameters"};

    //Constructor from file
    Lattice(std::string fname)
    {
        std::ifstream fin(fname, std::ios::binary);
        int lat_size[4];
        fin.read(reinterpret_cast<char*>(lat_size), sizeof(lat_size));
        LT = lat_size[0];
        LS = lat_size[1];
        fin.close();
        init();
    }
    //Constructor from values
    Lattice(int LT, int LS) : LT(LT), LS(LS) { init(); }
    //Constructor from command line arguments
    Lattice(int argc, char* argv[])
    {
        lattice_params.add_options()
            ("LT", po::value<int>(&LT)->default_value(3), "Temporal lattice size")
            ("LS", po::value<int>(&LS)->default_value(24),  "Spatial lattice size");

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(lattice_params).allow_unregistered().run(), vm);
        po::notify(vm);

        // Check if help is needed
        if(vm.count("help")) { std::cout << lattice_params << std::endl; }

        init();
    }
    //Shift functions
    int shift_fwd(int x, uint mu)
    {
        int nx = x + block_size[mu+1];
        int tx = x/block_size[mu+1];
        if((tx%size[mu])+1==size[mu])  // Anti-periodic boundary condition
            nx -= block_size[mu];
        return nx;
    };
    int shift_bwd(int x, uint mu)
    {
        int nx = x - block_size[mu+1];
        int tx = x/block_size[mu+1];
        if(tx%size[mu]==0)  // Anti-periodic boundary condition
            nx += block_size[mu];
        return nx;
    };
    //coordinate functions
    void get_coords(int x, int* coords)
    {
        for(uint mu=0; mu<4; mu++)
            coords[mu] = (x/block_size[mu+1])%size[mu];
    };
    friend ostream& operator<<(ostream& os, const Lattice& lattice)
    {
        os << "\tLattice size: " << ansi::yellow << " [" << lattice.size[0] << ", " << lattice.size[1] << ", " << lattice.size[2] << ", " << lattice.size[3] << "]" << ansi::reset << endl;
        os << "\tVolume:       " << ansi::yellow << lattice.vol << ansi::reset << endl;
        os << "\t3D Volume:    " << ansi::yellow << lattice.vol3D << ansi::reset <<  endl;
        os << "\tBlock sizes:  " << ansi::yellow << " [" << lattice.block_size[0] << ", " << lattice.block_size[1] << ", " << lattice.block_size[2] << ", " << lattice.block_size[3] << ", " << lattice.block_size[4] << "]" << ansi::reset << endl;
        os << "\tBeta = 1/T:   " << ansi::yellow << lattice.beta << ansi::reset << endl;
        return os;
    }
};

#endif // _LATTICE_HPP

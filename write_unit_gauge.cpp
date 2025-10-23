#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <string>

typedef std::complex<double> t_complex;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " LT LS" << std::endl;
        return 1;
    }

    int LT = std::atoi(argv[1]);
    int LS = std::atoi(argv[2]);
    int size[4] = {LT, LS, LS, LS};
    int vol = LT * LS * LS * LS;

    double plaq0 = 3.0;

    // Generate filename based on lattice size
    std::string filename = "unit_LT" + std::to_string(LT) + "_LS" + std::to_string(LS) + ".bin";

    std::ofstream fout(filename, std::ios::binary);
    if(!fout) {
        std::cerr << "Error opening file " << filename << " for writing!" << std::endl;
        return 1;
    }

    // Write lattice size and initial plaquette
    fout.write(reinterpret_cast<char*>(size), sizeof(size));
    fout.write(reinterpret_cast<char*>(&plaq0), sizeof(plaq0));

    // Fill gauge links with identity matrices for odd sites
    for(int x = 0; x < vol; ++x) {
        int xc[4];
        xc[0] = x / (LS*LS*LS);
        xc[1] = (x / (LS*LS)) % LS;
        xc[2] = (x / LS) % LS;
        xc[3] = x % LS;
        int sum = xc[0] + xc[1] + xc[2] + xc[3];

        if(sum % 2 == 1) {  // odd sites
            for(int mu = 0; mu < 4; ++mu) {
                t_complex identity[9];
                for(int i = 0; i < 9; ++i)
                    identity[i] = (i%4==0) ? t_complex(1.0, 0.0) : t_complex(0.0, 0.0);

                fout.write(reinterpret_cast<char*>(identity), sizeof(identity));
                fout.write(reinterpret_cast<char*>(identity), sizeof(identity));
            }
        }
    }

    fout.close();
    std::cout << "Identity gauge configuration written to " << filename << std::endl;
    return 0;
}

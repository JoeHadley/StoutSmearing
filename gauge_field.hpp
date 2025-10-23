#ifndef _GAUGE_FIELD_HPP
#define _GAUGE_FIELD_HPP

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#include "lattice.hpp"
#include "color_spinor.hpp"
#include "ansi_io.hpp"
#include "linalg.hpp"
#include "su3maximization.hpp"

class GaugeField : public Lattice {
private:
    void init();
public:
    GaugeField(const Lattice& lat) : Lattice(lat) {init();}
    GaugeField(string fname);
    GaugeField(int argc, char* argv[]) : Lattice(argc, argv) {init();}
    ~GaugeField() {delete[] Links;}
    string ConfigFileName;

    t_complex* Links;
    int        LinksSize;
    t_complex* Link(int x, uint mu) {return Links + x*4*9 + mu*9;};
    t_complex* Link(int t, int xs, uint mu) {return Links + (t*vol3D + xs)*4*9 + mu*9;};
    void       Plaquette(t_complex* res, int x, uint mu, uint nu);
    double     ReTrPlaquette(int x, uint mu, uint nu);
    double     MeanSpatialPlaquette();
    double     MeanTemporalPlaquette();
    double     MeanPlaquette();
    void       PolyakovLoop(t_complex* res, int xs);
    t_complex  MeanPolyakovLoop();
    void       WriteToFile(string fname);
    void       StaticGauge();
    t_complex  StaticProjection(int ProjectionMode=0, bool diagnostic_output=true);
    void       RandomizeLinks();
};

/*std::ostream& operator<<(std::ostream& os, const GaugeField& U) 
{
    for (int x=0; x<U.lat->vol; x++) 
        for (uint mu=0; mu<4; mu++) //TODO: Link is a pointer, print out the content - create a function for this
            os << "\t x=" << x << ", mu=" << mu << ": " << ansi::cyan << U.Link(x, mu) << ansi::reset << std::endl;
    return os;
}*/

#endif // _GAUGE_FIELD_HPP
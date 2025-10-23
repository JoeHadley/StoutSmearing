#include "gauge_field.hpp"

void GaugeField::init()
{
    LinksSize = vol*4*9;
    Links = new t_complex[LinksSize];
    for(int x=0; x<vol; x++)
        for(uint mu=0; mu<4; mu++)
            c3x3_set_unity(Link(x, mu));
};

void GaugeField::Plaquette(t_complex* res, int x, uint mu, uint nu)
{
    t_complex U1[9], U2[9];
    c3x3_times_c3x3(U1, Link(x, mu), Link(shift_fwd(x, mu), nu));
    c3x3_times_c3x3(U2, Link(x, nu), Link(shift_fwd(x, nu), mu));
    c3x3_times_c3x3_conj(res, U1, U2);
}

double GaugeField::ReTrPlaquette(int x, uint mu, uint nu)
{
    t_complex res[9];
    Plaquette(res, x, mu, nu);
    return real(tr_c3x3(res));
}

double GaugeField::MeanSpatialPlaquette()
{
    double res = 0.0;
    for(int x=0; x<vol; x++)
        for(uint mu=1; mu<4; mu++)
            for(uint nu=mu+1; nu<4; nu++)
                res += ReTrPlaquette(x, mu, nu);
    return res/(3.0*vol);
}

double GaugeField::MeanTemporalPlaquette()
{
    double res = 0.0;
    for(int x=0; x<vol; x++)
        for(uint mu=1; mu<4; mu++)
                res += ReTrPlaquette(x, mu, 0);
    return res/(3.0*vol);
}

double GaugeField::MeanPlaquette()
{
    return 0.5*(MeanSpatialPlaquette() + MeanTemporalPlaquette());
}

void GaugeField::PolyakovLoop(t_complex* res, int xs)
{
    t_complex nP[9];
    c3x3_set_unity(res);
    for(int t=0; t<LT; t++)
    {
        int x = t*vol3D + xs;
        c3x3_times_c3x3(nP, res, Link(x, 0));
        std::copy(nP, nP+9, res);
    };
}

t_complex GaugeField::MeanPolyakovLoop()
{
    t_complex res = 0.0;
    for(int xs=0; xs<vol3D; xs++)
    {
        t_complex P[9];
        PolyakovLoop(P, xs);
        res += tr_c3x3(P)/3.0;
    };
    return res/(double)(vol3D);
}

GaugeField::GaugeField(string fname) : Lattice(fname)
{
    // Getting just the configuration file name
    size_t pos = fname.find_last_of("/\\");
    ConfigFileName = (pos == std::string::npos) ? fname : fname.substr(pos + 1);

    std::ifstream fin(fname, ios::binary);
    int lat_size[4];
    double plaq0;
    fin.read((char*)lat_size, sizeof(lat_size));
    fin.read((char*)&plaq0, sizeof(plaq0));
    cout << "lat_size: " << lat_size[0] << " " << lat_size[1] << " " << lat_size[2] << " " << lat_size[3] << endl;

    init();

    for(int x=0; x<vol; x++)
    {
        int xc[4];
        get_coords(x, xc);
        if((xc[0] + xc[1] + xc[2] + xc[3])%2==1) //This is a definition of an odd lattice site
            for(uint mu=0; mu<4; mu++)
            {
                fin.read((char*)Link(               x, mu), sizeof(t_complex)*9);
                fin.read((char*)Link(shift_bwd(x, mu), mu), sizeof(t_complex)*9);
            };
    };
    fin.close();
        
    double max_unitarity_err  = 0.0;
    for(int x=0; x<vol; x++)
        for(uint mu=0; mu<4; mu++)
            max_unitarity_err = max(max_unitarity_err, c3x3_unitarity_norm(Link(x, mu)));

    double plaq1 = MeanPlaquette();

    cout << setprecision(16) << "plaq1 = " << plaq1 << ", plaq0 = " << plaq0 << endl;   

    if(abs(plaq1 - plaq0)>1.0E-8 || max_unitarity_err>1.0E-6)
    {
        cerr << ansi::red << "Error reading gauge field from file " << ansi::cyan << fname << ":" << ansi::reset << endl;
        cerr << "\t" << ansi::red << "plaq0 = " << ansi::yellow << plaq0 << ansi::red << ", plaq1 = " << ansi::yellow << plaq1 << ansi::reset << endl;
        cerr << "\t" << ansi::red << "Unitarity error: " << ansi::yellow << max_unitarity_err << ansi::reset << endl;
        cerr << endl;
    }
    else
    {
        cout << ansi::green  << "Gauge field read from file " << ansi::cyan << fname << ansi::green << " successfully, ConfigFileName = " << ansi::cyan << ConfigFileName << ansi::reset << endl;
        cout << ansi::green  << "\tMean spatial and temporal plaquettes: ";
        cout << ansi::yellow << MeanSpatialPlaquette() << ansi::reset << ", ";
        cout << ansi::yellow << MeanTemporalPlaquette() << ansi::reset << endl;
        cout << ansi::green  << "\tMean Polyakov loop:                   ";
        cout << ansi::yellow << MeanPolyakovLoop() << ansi::reset << endl;
        cout << endl;
    };
}

void GaugeField::WriteToFile(string fname)
{
    std::ofstream fout(fname, ios::binary);
    if (!fout.is_open()) {
        cerr << ansi::red << "Error: Cannot open file " << ansi::cyan << fname << ansi::red << " for writing." << ansi::reset << endl;
        return;
    }
    
    // Write header: lattice size and mean plaquette
    int lat_size[4] = {LT, LS, LS, LS};
    double plaq0 = MeanPlaquette();
    fout.write((char*)lat_size, sizeof(lat_size));
    fout.write((char*)&plaq0, sizeof(plaq0));
    
    // Write gauge links in the same format as read
    for(int x=0; x<vol; x++)
    {
        int xc[4];
        get_coords(x, xc);
        if((xc[0] + xc[1] + xc[2] + xc[3])%2==1) // Only odd lattice sites
            for(uint mu=0; mu<4; mu++)
            {
                fout.write((char*)Link(               x, mu), sizeof(t_complex)*9);
                fout.write((char*)Link(shift_bwd(x, mu), mu), sizeof(t_complex)*9);
            };
    };
    
    fout.close();
    
    cout << ansi::green << "Gauge field written to file " << ansi::cyan << fname << ansi::green << " successfully." << ansi::reset << endl;
}

void GaugeField::StaticGauge()
{
    double max_gauge_err = 0.0;
    for(int xs=0; xs<vol3D; xs++)
    {
        t_complex P[NCOLOR2], tmp[NCOLOR2], evals[NCOLOR], levecs[NCOLOR2], revecs[NCOLOR2], U0A[NCOLOR2], evals_d[NCOLOR2];
        PolyakovLoop(P, xs);
        //Calculating the LT's root of P - this will be the value of the time-like links in the static gauge
        copy(P, P+NCOLOR2, tmp);
        diagonalize(tmp, evals, revecs, levecs, NCOLOR);
        fill(evals_d, evals_d+NCOLOR2, t_complex(0.0, 0.0));
        for(int i=0; i<NCOLOR; i++) evals_d[i*NCOLOR + i] = pow(evals[i], 1.0/((double)LT));
        c3x3_times_c3x3(tmp, evals_d, levecs); //tmp = D*L
        c3x3_times_c3x3(U0A, revecs, tmp); //Now U0A is the LT's root of P

        //Now gauge transform all time-like links at (t, xs) so that they become U0A
        //Space-like links are also transformed accordingly
        t_complex omega[NCOLOR2];
        c3x3_set_unity(omega); //We assume that the gauge transformation at t=0 is unity
        for(int t=0; t<LT; t++)
        {
            //Apply the current gauge transformation to all links at (t, xs)
            for(int i=1; i<=3; i++)
            {
                //Gauge-transforming the forward link attached to site xs
                std::copy(Link(t, xs, i), Link(t, xs, i)+NCOLOR2, tmp); //tmp = Link(t, xs, i)
                c3x3_times_c3x3(Link(t, xs, i), omega, tmp); //Link(t, xs, i) -> omega * Link(t, xs, i)
                int xsb = shift_bwd(t*vol3D + xs, i);
                //Gauge-transforming the backward link attached to site xs
                std::copy(Link(xsb, i), Link(xsb, i)+NCOLOR2, tmp); //tmp = Link(xsb, i)
                c3x3_times_c3x3_conj(Link(xsb, i), tmp, omega); //Link(xsb, i) -> Link(xsb, i) * omega^+
            };

            //Compute the gauge transformation matrix for the next time slice
            c3x3_times_c3x3(tmp, omega, Link(t, xs, 0)); //Link(t, xs, 0) is the time-like link at (t, xs)
            c3x3_conj_times_c3x3(omega, U0A, tmp);  //omega = U0A^+ * tmp
            std::copy(U0A, U0A+NCOLOR2, Link(t, xs, 0)); //Set the time-like link to U0A
        };
        c3x3_set_unity(tmp);
        double err = norm_diff(omega, tmp, NCOLOR2);
        max_gauge_err = max(max_gauge_err, err);
    };
    cout << ansi::green << "Max error of omega after static gauge fixing = " << ansi::magenta << max_gauge_err << ansi::reset << endl;
}

t_complex GaugeField::StaticProjection(int ProjectionMode, bool diagnostic_output)
{
    // First apply static gauge fixing
    StaticGauge();

    //Determine the Z3 center sector and replace the time-like links at the last time slice with the corresponding center element
    t_complex PolyakovLoopAvg = MeanPolyakovLoop();
    double angle = std::arg(PolyakovLoopAvg);
    angle = (abs(angle)<pi/3 ? 0 : 2*pi/3*angle/abs(angle)); //angle/abs(angle) is just the sign of angle
    t_complex Z = std::exp(t_complex(0.0, angle));

    std::cout << ansi::green << "Polyakov loop average:    " << ansi::magenta << PolyakovLoopAvg << ansi::reset << endl;
    std::cout << ansi::green << "Center sector projection: " << ansi::magenta << Z << ansi::reset << endl;

    for(int xs=0; xs<vol3D; xs++)
        for(int it=0; it<LT; it++)
        {
            t_complex* link = Link(it, xs, 0);
            std::fill(link, link+NCOLOR2, t_complex(0.0, 0.0));
            for(int i=0; i<NCOLOR; i++)
                link[i*NCOLOR + i] = (it==LT-1? Z : t_complex(1.0, 0.0));
        };

    if(ProjectionMode==0)
    {
        //All links are set to be equal to the links at t=0
        for(int it=1; it<LT; it++)
            for(int xs=0; xs<vol3D; xs++)
                for(int i=1; i<=3; i++) // Spatial directions only
                    std::copy(Link(0, xs, i), Link(0, xs, i)+NCOLOR2, Link(it, xs, i));
    };

    if(ProjectionMode==1)
    {
        double max_det_err = 0.0;
        double max_unitarity_err = 0.0;
        double max_projection_err = 0.0;
        double min_overlap = 1.0E10, max_overlap = -1.0E10, mean_overlap = 0.0;
        int max_iter = 0; double mean_iter = 0.0;
        
        // For each spatial point and spatial direction, average over time
        for(int xs=0; xs<vol3D; xs++)
            for(int i=1; i<=3; i++) // Spatial directions only
            {
                t_complex avg_link[NCOLOR2]; fill(avg_link, avg_link+NCOLOR2, t_complex(0.0, 0.0));
                
                // Sum all time slices
                for(int t=0; t<LT; t++)
                {
                    t_complex* current_link = Link(t, xs, i);
                    for(int a=0; a<NCOLOR2; a++) avg_link[a] += current_link[a];
                };
                for(int a=0; a<NCOLOR2; a++) avg_link[a] /= (double)LT;
                
                t_complex U[NCOLOR2];
                // Project back to SU(3)
                int iter = MaximizeSU3Overlap(avg_link, U);
                if(iter<0)
                {
                    cerr << ansi::red << "Error: MaximizeSU3Overlap failed at xs=" << xs << ", i=" << i << ", error code " << iter << ansi::reset << endl;
                    continue;
                };

                if(diagnostic_output)
                {
                    max_iter = max(max_iter, iter);
                    mean_iter += (double)iter;

                    //Check the projection quality
                    t_complex tmp1[NCOLOR2], tmp2[NCOLOR2];
                    c3x3_conj_times_c3x3(tmp1, U, avg_link);
                    c3x3_conj_times_c3x3(tmp2, avg_link, U);
                    t_complex diag = (tr_c3x3(tmp1) - tr_c3x3(tmp2))/3.0;
                    double projection_err  = 0.0;
                    for(int i=0; i<NCOLOR; i++)
                        for(int j=0; j<NCOLOR; j++)
                        {
                            t_complex d = tmp1[i*NCOLOR + j] - tmp2[i*NCOLOR + j];
                            if(i==j) d -= diag;
                            projection_err += norm(d);
                        }
                    max_projection_err = max(max_projection_err, sqrt(projection_err));

                    //And finally the actual value of the overlap
                    double overlap = real(tr_c3x3(tmp1))/3.0;
                    min_overlap = min(min_overlap, overlap);
                    max_overlap = max(max_overlap, overlap);
                    mean_overlap += overlap;

                    t_complex det = c3x3_det(U);
                    max_det_err = max(max_det_err, abs(1.0 - det));
                    max_unitarity_err = max(max_unitarity_err, c3x3_unitarity_norm(U));
                };

                // Replace all time slices with the projected averaged link
                for(int t=0; t<LT; t++) copy(U, U+NCOLOR2, Link(t, xs, i));
            };

        if(diagnostic_output)
        {
            mean_overlap /= (double)(vol3D*3);
            mean_iter /= (double)(vol3D*3);

            std::cout << ansi::green << "Max. determinant error =  " << ansi::magenta << max_det_err << ansi::reset << endl;
            std::cout << ansi::green << "Max. unitarity error   =  " << ansi::magenta << max_unitarity_err << ansi::reset << endl;
            std::cout << ansi::green << "Max. projection error  =  " << ansi::magenta << max_projection_err << ansi::reset << endl;
            std::cout << ansi::green << "Range of overlap values:  " << ansi::magenta << min_overlap << " to " << max_overlap << " (mean = " << mean_overlap << ")" << ansi::reset << endl;
            std::cout << ansi::green << "Max. iterations:          " << ansi::magenta << max_iter << ansi::reset << endl;
            std::cout << ansi::green << "Mean iterations:          " << ansi::magenta << mean_iter << ansi::reset << endl;
        };
    }; //End of ProjectionMode==1

    return Z;
}

void GaugeField::RandomizeLinks()
{
    for(int x=0; x<vol; x++)
        for(uint mu=0; mu<4; mu++)
        {
            t_complex random_matrix[NCOLOR2];
            for(int i=0; i<NCOLOR2; i++)
                random_matrix[i] = t_complex(drand48(), drand48());
            int iter = MaximizeSU3Overlap(random_matrix, Link(x, mu));
            if(iter<0) { cerr << ansi::red << "Error: MaximizeSU3Overlap failed at x=" << x << ", mu=" << mu << ", error code " << iter << ansi::reset << endl; continue; }; 
        };
}

//These checks were a part of the StaticGauge function
//Checking if the diagonalization was successful
        // t_complex check[9], evals_d[9];
        // std::fill(evals_d, evals_d+9, t_complex(0.0, 0.0));
        // for(int i=0; i<3; i++)
        //     evals_d[3*i + i] = evals[i];
        
        // c3x3_times_c3x3(tmp, evals_d, levecs); //tmp = D*L
        // c3x3_times_c3x3(check, revecs, tmp);   //check = R*D*L
        // double err = norm_diff(P, check, 9);
        // max_diag_err = max(max_diag_err, err);

        //Checking the root calculation
        // t_complex check[NCOLOR2]; c3x3_set_unity(check);
        // for(int i=0; i<LT; i++)
        // {
        //     c3x3_times_c3x3(tmp, check, U0A); //tmp = check*U0A
        //     copy(tmp, tmp+NCOLOR2, check);
        // };

        // double err = norm_diff(P, check, NCOLOR2);
        // max_root_err = max(max_root_err, err);
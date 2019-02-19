//
//  main.cpp
//  opt_pricing_binomial_american_put
//
//  Created by Changjie Ma on 2/16/19.
//  Copyright © 2019 Changjie Ma. All rights reserved.
//  credit to http://www.goddardconsulting.ca/option-pricing-binomial-alts.html#crrdrift
#include "opt_cmpt.h"
#include <iostream>
#include <fstream>
#include <math.h>
//#include <algorithm>
using namespace std;

int main(int argc, const char *argv[])
{
    // Comparing different parametrization methods
    ofstream write_a;
    write_a.open ("results_a.csv");
    for (int i=25;i<=N_max;i++){
        set_step(i);
        write_a << delta_t << ",";

        crr();
        initialize();
        write_a << f(0,0) << ",";

        jarrow_rudd();
        initialize();
        write_a << f(0,0) << ",";

        tian();
        initialize();
        write_a << f(0,0) << "," << "\n";
    }
    write_a.close();
    cout << "done for 1.a!" << endl;


    // binomial Black-Scholes and Richardson extrapolation
    ofstream write_b;
    write_b.open("results_b.csv");
    for (int i = 25; i <= N_max; i++) {
        write_b << i << ",";

        set_step(i);
        crr();
        initialize();
        write_b << f(0, 0) << ",";

        set_step(i);
        crr();
        initialize();
        write_b << f_BS(0, 0) << ",";

        write_b << f_BBSR(i) << ",";

        write_b << f_BAM(i) << ",\n";

    }
    write_b.close();
    cout << "done for 1.b!" << endl;

    // Trinomial model
    ofstream write_c;
    write_c.open("results_c.csv");
    for (int i = 25; i <= N_max; i++) {
        set_step(i);
        write_c << i << ",";

        tri();
        initialize();
        write_c << f_tri(0, 0) << ",\n";

    }
    write_c.close();
    cout << "done for 1.c!" << endl;
    
// Binomial pyramid for pricing basket option
    
    K = 100;
    int step = 60;
    _initialize();
    set_step(step);
    clock_t t = clock();
    float price = f_pyramid(0,0,0,0);
    t = clock() - t;
    free_map();
    double cpu_time = (float)t / CLOCKS_PER_SEC;
    cout << step << " steps, price = " << price << ", cpu time = " << cpu_time << "s" << endl;
    return 0;
}

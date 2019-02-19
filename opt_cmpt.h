//
//  opt_cmpt.h
//  opt_pricing_binomial_american_put
//
//  Created by Changjie Ma on 2/16/19.
//  Copyright © 2019 Changjie Ma. All rights reserved.
//

#ifndef opt_cmpt_h
#define opt_cmpt_h
#define N_max 500
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <numeric>
using namespace std;

// define features of options, euro put price calculate from B-S: 9.024846
float r = 0.03;
float sigma = 0.2;
float S_0 = 100;
float K = 105;
float T = 1;

// define features of binomial lattice model
// features doesn't change with different parametrization methods
int N;// number of steps
float delta_t;

// features change with different parametrization methods
float u;
float d;
float p;
float _p;

// Normal CDF
float Normal(const float& z)
{
    if (z > 6.0) {
        return 1.0;
    }                             // this guards against overflow
    if (z < -6.0) {
        return 0.0;
    }
    float b1 = 0.31938153;
    float b2 = -0.356563782;
    float b3 = 1.781477937;
    float b4 = -1.821255978;
    float b5 = 1.330274429;
    float p = 0.2316419;
    float c2 = 0.3989423;
    float a = fabs(z);
    float t = 1.0 / (1.0 + a * p);
    float b = c2 * exp((-z) * (z / 2.0));
    float n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
    n = 1.0 - b * n;
    if (z < 0.0) n = 1.0 - n;
    return n;
}

// Black-Scholes Function
float option_price_put_black_scholes(const float& S,      // spot price
                                     const float& Strike,       // Strike (exercise) price,
                                     const float& rate,       // interest rate
                                     const float& volatility,                                         // volatility
                                     const float& time)
{
    float time_sqrt = sqrt(time);
    float d1 = (log(S / Strike) + rate * time) / (volatility * time_sqrt) + 0.5 * volatility * time_sqrt;
    float d2 = d1 - (volatility * time_sqrt);
    return Strike * exp(-rate * time) * Normal(-d2) - S * Normal(-d1);
}

// This is the stock price for node(i,j)
float s(int i, int j)
{
    return pow(u, (i + j) / 2) * pow(d, (i - j) / 2) * S_0;
}

// This is the stock price for node(i,j) in binomial model
float _s(int i, int j)
{
    return pow(u, j) * S_0;
}

// This is the stock price for node(i,j) in trinomial model
float map[N_max][2 * N_max + 1]; // set up a 2d array for storage

// set up memorize
void initialize()
{
    for (int i = 0; i < N_max; i++) {
        for (int j = 0; j < 2 * N_max + 1; j++) {
            map[i][j] = -1;
        }
    }
}

// set steps
void set_step(int n)
{
    N = n;
    delta_t = T / N;
}

// set parameters
void crr()
{
    u = exp(sigma * sqrt(delta_t));
    d = exp(-sigma * sqrt(delta_t));
    p = (exp(r * delta_t) - d) / (u - d);
}
void jarrow_rudd()
{
    u = exp((r - sigma * sigma / 2) * delta_t + sigma * sqrt(delta_t));
    d = exp((r - sigma * sigma / 2) * delta_t - sigma * sqrt(delta_t));
    p = 0.5;
}
void tian()
{
    float niu = exp(sigma * sigma * delta_t);
    u = 0.5 * exp(r * delta_t) * niu * (niu + 1 + sqrt(niu * niu + 2 * niu - 3));
    d = 0.5 * exp(r * delta_t) * niu * (niu + 1 - sqrt(niu * niu + 2 * niu - 3));
    p = (exp(r * delta_t) - d) / (u - d);
}
void tri()
{
    float niu = r - 0.5 * sigma * sigma;
    float lambda = sqrt(3 / 2);
    u = exp(lambda * sigma * sqrt(delta_t));
    d = 1 / u;
    p = 1 / (2 * lambda * lambda) + niu * sqrt(delta_t) / (2 * lambda * sigma);
    _p = 1 - 1 / (lambda * lambda);
}

// using recurssion method to price options
float f(int i, int j)
{
    if (i < N) {
        if (map[i][i + j] != -1) return map[i][i + j];
        else {
            float result = p * (f(i + 1, j + 1)) + (1 - p) * f(i + 1, j - 1);
            map[i][i + j] = result;
            return max(K - s(i, j), result);
        }
    } else {
        return std::max<float>(0, K - s(i, j));
    }
}

// using BS+binomial model
float f_BS(int i, int j)
{
    if (i < N - 1) {
        if (map[i][i + j] != -1) return map[i][i + j];
        else {
            float result = p * (f_BS(i + 1, j + 1)) + (1 - p) * f_BS(i + 1, j - 1);
            map[i][i + j] = result;
            return max(K - s(i, j), result);
        }
    } else {
        return option_price_put_black_scholes(s(i, j), K, r, sigma, delta_t);
    }
}

// using BBSR+binomial model
float f_BBSR(int n)
{
    set_step(n);
    crr();
    initialize();
    float _v1 = f_BS(0, 0);

    set_step(n / 2);
    crr();
    initialize();
    float _v2 = f_BS(0, 0);

    return 4/3 * _v1 - 1/3*_v2;
}

// using BAM+binomial model
float f_BAM(int n)
{
    set_step(n);
    crr();
    initialize();
    float _v1 = f(0, 0);

    set_step(n - 1);
    crr();
    initialize();
    float _v2 = f(0, 0);

    return (_v1 + _v2) / 2;
}

// using trinomial lattice to conduct recurssion
float f_tri(int i, int j)
{
    if (i < (N)) {
        if (map[i][i + j] != -1) return map[i][i + j];
        else {
            float result = p * (f_tri(i + 1, j + 1)) + _p * (f_tri(i + 1, j)) + (1 - p - _p) * f_tri(i + 1, j - 1);
            map[i][i + j] = result;
            return max(K - _s(i, j), result);
        }
    } else return std::max<float>(0, K - _s(i, j));
}







// Question 2

// The memorization array for Q2
#define size 130
float ****_map;
void _initialize()
{
    _map = new float ***[size];
    for (int i = 0; i < size; i++) {
        _map[i] = new float **[size];
        for (int j = 0; j < size; j++) {
            _map[i][j] = new float *[size];
            for (int k = 0; k < size; k++) {
                _map[i][j][k] = new float[size];
            }
        }
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                for (int l = 0; l < size; l++) {
                    _map[i][j][k][l] = -1;
                }
            }
        }
    }
}
// Free memory
void free_map()
{
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                free(_map[i][j][k]);
            }
            free(_map[i][j]);
        }
        free(_map[i]);
    }
    free(_map);
}

// using binomial pyramid to price basket options
// http://www.haas.berkeley.edu/MFE/download/student_papers/mfe02_wan-pricing_basket_options.pdf
float s_pyramid_1(int i, int j, int k, int l)
{
    return S_0 * exp(0.03 * delta_t + sqrt(delta_t) * (0.2000000 * j - 0.1789975 * k - 0.1962642 * l));
}

float s_pyramid_2(int i, int j, int k, int l)
{
    return S_0 * exp(0.03 * delta_t + sqrt(delta_t) * (-0.2000000 * j + 0.1789975 * k - 0.2000828 * l));
}

float s_pyramid_3(int i, int j, int k, int l)
{
    return S_0 * exp(0.03 * delta_t + sqrt(delta_t) * (-0.2000000 * j - 0.2189975 * k + 0.1600828 * l));
}
float f_pyramid(int i, int j, int k, int l)
{
    if (i < N) {
        if (_map[i][j + i][k + i][l + i] != -1) {
            return _map[i][j + i][k + i][l + i];
        } else {
            float result =
                (f_pyramid(i + 1, j + 1, k + 1, l + 1) + f_pyramid(i + 1, j + 1, k + 1, l - 1) +
                 f_pyramid(i + 1, j + 1, k - 1, l + 1) + f_pyramid(i + 1, j + 1, k - 1, l) +
                 f_pyramid(i + 1, j - 1, k + 1, l + 1) + f_pyramid(i + 1, j - 1, k + 1, l) +
                 f_pyramid(i + 1, j - 1, k - 1, l + 1) + f_pyramid(i + 1, j - 1, k - 1, l - 1)) / 8;
            _map[i][j + i][k + i][l + i] = result;
            return max(max(max(s_pyramid_1(i, j, k, l), s_pyramid_2(i, j, k, l)), s_pyramid_3(i, j, k, l)) - K, result);
        }
    } else
//        return std::max<float>(0,K-s(i,j));
        return max<float>(max(max(s_pyramid_1(i, j, k, l), s_pyramid_2(i, j, k, l)), s_pyramid_3(i, j, k, l)) - K, 0);
}





// MC for Q2
//float get_uniform()
//{
//    return (((float)random()) / (pow(2.0, 31.0) - 1.0));
//}
//int no_of_trials = 100;
//
//float monte_carlo(int no_of_divisions){
//
//    float strike_price = 100;
//
//    float delta_T = T / ((float)no_of_divisions);
//    float delta_R = (r - 0.5 * pow(sigma, 2)) * delta_T;
//    float delta_SD = sigma * sqrt(delta_T);
//
//    float payoff=0;
//
//    for (int i = 0; i < no_of_trials; i++) {
//        // by sharing random variables we create 4 paths
//        float current_stock_price1 = S_0;
//        float current_stock_price2 = S_0;
//        float current_stock_price3 = S_0;
//        //        float current_stock_price4 = S_0;
//
//        for (int j = 0; j < no_of_divisions; j++) {
//            // create the unit normal variates using the Box-Muller Transform
//            float x = get_uniform();
//            float y = get_uniform();
//            float a =  sqrt(-2.0 * log(x)) * cos(6.283185307999998 * y);
//            float b =  sqrt(-2.0 * log(x)) * sin(6.283185307999998 * y);
//
//            float _x = get_uniform();
//            float _y = get_uniform();
//            float c =  sqrt(-2.0 * log(_x)) * cos(6.283185307999998 * _y);
//
//            float e1 = a;
//            float e2 = 0.1*a+0.9949874*b;
//            float e3 = 0.1*a+0.0904534*b+0.9908674*c;
//
//            current_stock_price1 = current_stock_price1 * exp(delta_R - delta_SD * e1);
//            current_stock_price2 = current_stock_price1 * exp(delta_R - delta_SD * e2);
//            current_stock_price3 = current_stock_price1 * exp(delta_R - delta_SD * e3);
//        }
//        payoff=payoff+max<float>(max(max(current_stock_price1,
//                                         current_stock_price2),
//                                     current_stock_price3)-strike_price,0);
//    }
//    float average = payoff/no_of_trials;
//    return average;
//}

#endif /* opt_cmpt_h */

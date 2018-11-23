//
//  main.cpp
//  black_scholes_formula
//  Routine that prices European call and put options using the Black-Scholes formula
//  Input: present time, spot price, strike, maturity date, volatility, constant interest rate, continuous dividend rate
//  Reference Stefanica, Primer for the Maths of Fin Eng, Section 3.8
//
//  Created by carlos on 23/11/2018.
//  Copyright © 2018 carlos. All rights reserved.
//

#include <iostream>
#include <cmath>
using namespace std;

#include <iomanip> // to use setprecision() so that std output gives me more decimal places

//  function that computes the cumulative distribution Z up to t
double cum_dist_normal(double t) {
    double z = abs(t);
    double y = 1/(1+0.2316419*z);
    double a[] = {0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429};
    double m = 1-exp(-t*t/2)*(a[0]*y+a[1]*y*y+a[2]*y*y*y+a[3]*y*y*y*y+a[4]*y*y*y*y*y)/(sqrt(2*M_PI));
    double nn;
    
    if (t>0) {
        nn = m;
    }
    else {
        nn = 1-m;
    }
    return nn;
}

int main() {
    double t, S, K, T, sigma, r, q;
    
    cout << "present time t (often equal to 0): ";
    cin >> t;
    cout << "spot price of the underlying asset (at time t) S: ";
    cin >> S;
    cout << "option strike K: ";
    cin >> K;
    cout << "maturity date T (time to maturity is T-t): ";
    cin >> T;
    cout << "volatility of the underlying asset sigma: ";
    cin >> sigma;
    cout << "constant interest rate r: ";
    cin >> r;
    cout << "continuous dividend rate of the underlying asset q: ";
    cin >> q;

    double d1 = (log(S/K)+(r-q+sigma*sigma/2)*(T-t))/(sigma*sqrt(T-t));
    double d2 = d1-sigma*sqrt(T-t);
    
    cout << "\nC = " << setprecision (8) << S*exp(-q*(T-t))*cum_dist_normal(d1) - K*exp(-r*(T-t))*cum_dist_normal(d2) << endl;
    cout << "P = " << setprecision (8) << K*exp(-r*(T-t))*cum_dist_normal(-d2) - S*exp(-q*(T-t))*cum_dist_normal(-d1) << endl;
 //   cout << cum_dist_normal(d1) << endl;
    return 0;
}

#include "itensor/all.h"
#include "header.h"

namespace lv = liouvillian;

double liouvillian::calc_rel_current(const ITensor& rho, const int& site, const double& J)
{
    // get indices
    IndexSet s = rho.inds();
    int N = s.r();
    int i = site;
    
    if ((N-1) < i) {printf("WARNING: In function calc_rel_current: given site is larger than total number of sites");}
    
    // initialize magnetization
    Complex rel_current = 0.0;
    
    // get coefficient for current on site i:
    // only relevant coefficient for the trace is the one with all sites in state 0.
    // so all sites must be in state 0, only site i in state 2 and site i+1 in state 3. Then applying of sigma_z,i on site i will change these states to 0 as well. 
    // Multiplication with 2 (trace of the 2x2 unitary matrix of the ith subspace) gets us tr(sigma_z,i*rho).
    
    int k = (N-2+i) - ((N-2+i) >= N)*N;
    int l = (N-3+i) - ((N-3+i) >= N)*N;
    int m = (N-4+i) - ((N-4+i) >= N)*N;
    int n = (N-5+i) - ((N-5+i) >= N)*N;
    int o = (N-6+i) - ((N-6+i) >= N)*N; 
    int p = (N-7+i) - ((N-7+i) >= N)*N;
    int q = (N-8+i) - ((N-8+i) >= N)*N;
    int r = (N-9+i) - ((N-9+i) >= N)*N;
    int ss = (N-10+i) - ((N-10+i) >= N)*N;
    int t = (N-11+i) - ((N-11+i) >= N)*N; 
    int u = (N-12+i) - ((N-12+i) >= N)*N;
    int v = (N-13+i) - ((N-13+i) >= N)*N;
    int w = (N-14+i) - ((N-14+i) >= N)*N;
    int x = (N-15+i) - ((N-15+i) >= N)*N;
    int y = (N-16+i) - ((N-16+i) >= N)*N;
    int z = (N-17+i) - ((N-17+i) >= N)*N; 
    int a = (N-18+i) - ((N-18+i) >= N)*N;
    int b = (N-19+i) - ((N-19+i) >= N)*N;
    
    if (N == 2) {rel_current = (rho.cplx(s[0](2),s[1](3)) - rho.cplx(s[0](3),s[1](2)));}
    else if (N == 3) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1)));}
    else if (N == 4) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1)));}
    else if (N == 5) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1)));}
    else if (N == 6) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1)));}
    else if (N == 7) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1)));}
    else if (N == 8) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1)));}
    else if (N == 9) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1)));}
    else if (N == 10) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1)));}
    else if (N == 11) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1)));}
    else if (N == 12) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1)));}
    else if (N == 13) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1)));}
    else if (N == 14) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1)));}
    else if (N == 15) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1)));}
    else if (N == 16) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1)));}
    else if (N == 17) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1),s[y](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1),s[y](1)));}
    else if (N == 18) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1),s[y](1),s[z](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1),s[y](1),s[z](1)));}
    else if (N == 19) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1),s[y](1),s[z](1),s[a](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1),s[y](1),s[z](1),s[a](1)));}
    else if (N == 20) {rel_current = (rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1),s[y](1),s[z](1),s[a](1),s[b](1)) - rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1),s[y](1),s[z](1),s[a](1),s[b](1)));}
    
    else std::cout << "WARNING, current calculation only implemented for N<=20" << std::endl;
    
        
//     PrintData(rho.cplx(s[i-1](2),s[i](3),s[k](1),s[l](1)));
//     PrintData(rho.cplx(s[i-1](3),s[i](2),s[k](1),s[l](1)));
//     std::cout << " looking at coefficient i-1iklmn: " << (i-1) << i << k  << l << m << n << std::endl;
    
    //rescale magnetization for compatibility with spin 1/2 picture
    //rel_current *= J*0.25*pow(2.0,N);
    
    // calculate relative current without spin-1/2 factor
    rel_current *= pow(2.0,N)*0.5;
    
    return (-1.0)*real(rel_current); 
}

        

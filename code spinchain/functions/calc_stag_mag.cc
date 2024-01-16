#include "itensor/all.h"
#include "header.h"

namespace lv = liouvillian;

double liouvillian::calc_stag_mag(ITensor& rho)
{
    // get indices
    IndexSet s = rho.inds();
    int N = s.r();
    
    // initialize magnetization
    Complex stag_mag = 0.0;
    Complex site_mag = 0.0;
    
    // get coefficient for magnetization on site i:
    // only relevant coefficient for the trace is the one with all sites in state 0.
    // so all sites must be in state 0, only site i in state 3. Then applying of sigma_z,i on site i will change this 
    // coefficient to 0 as well. Multiplication with 2 (trace of the 2x2 unitary matrix of the ith subspace) gets us tr(sigma_z,i*rho).
    for(int i = 0; i < N; i++)
    {
        int j = (N-1+i) - ((N-1+i) >= N)*N; 
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
        
        if (N == 2) {site_mag = rho.real(s[i](4),s[j](1));}
        else if (N == 3) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1));}
        else if (N == 4) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1));}
        else if (N == 5) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1));}
        else if (N == 6) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1));}
        
        else if (N == 7) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1));}
        else if (N == 8) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1));}
        else if (N == 9) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1));}
        else if (N == 10) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1));}
        else if (N == 11) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1));}
        else if (N == 12) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1));}
        else if (N == 13) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1));}
        else if (N == 14) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1));}
        else if (N == 15) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1));}
        else if (N == 16) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1));}
        else if (N == 17) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1),s[y](1));}
        else if (N == 18) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1),s[y](1),s[z](1));}
        else if (N == 19) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1),s[y](1),s[z](1),s[a](1));}
        else if (N == 20) {site_mag = rho.cplx(s[i](4),s[j](1),s[k](1),s[l](1),s[m](1),s[n](1),s[o](1),s[p](1),s[q](1),s[r](1),s[ss](1),s[t](1),s[u](1),s[v](1),s[w](1),s[x](1),s[y](1),s[z](1),s[a](1),s[b](1));}
        
        else std::cout << "Warning, mag calculation only implemented for N<=20" << std::endl;
        
        //rescale magnetization for compatibility with spin 1/2 picture
        site_mag *= 0.5*pow(2.0,N);
        
        if (i%2 == 0) {stag_mag += site_mag/N;}
        else {stag_mag -= site_mag/N;}
//         std::cout << " looking at coefficient ijklm: " << i << j  << k  << l << m << n << std::endl;
//         std::cout << "site_mag: " << site_mag << " stag_mag: " << stag_mag << std::endl;
    }
    
    return (-1.0)*real(stag_mag); 
}

        

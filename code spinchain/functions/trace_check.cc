#include "itensor/all.h"
#include "header.h"

double liouvillian::trace_check(const ITensor& rho)
{
    // initial value for trace is zero
    Complex trace = 0.0;
    
    // get indices
    IndexSet s = rho.inds();
    int N = s.r();
    
    // this part gets the coefficient c_0000 only for varying sizes of the Tensor rho:
    if (N == 2) {trace+=rho.cplx(s[0](1),s[1](1));}
    else if (N == 3) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1));}
    else if (N == 4) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1));}
    else if (N == 5) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1));}
    else if (N == 6) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1));}
    else if (N == 7) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1));}
    else if (N == 8) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1));}
    else if (N == 9) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1),s[8](1));}
    else if (N == 10) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1),s[8](1),s[9](1));}
    else if (N == 11) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1),s[8](1),s[9](1),s[10](1));}
    else if (N == 12) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1),s[8](1),s[9](1),s[10](1),s[11](1));}
    else if (N == 13) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1),s[8](1),s[9](1),s[10](1),s[11](1),s[12](1));}
    else if (N == 14) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1),s[8](1),s[9](1),s[10](1),s[11](1),s[12](1),s[13](1));}
    else if (N == 15) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1),s[8](1),s[9](1),s[10](1),s[11](1),s[12](1),s[13](1),s[14](1));}
    else if (N == 16) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1),s[8](1),s[9](1),s[10](1),s[11](1),s[12](1),s[13](1),s[14](1),s[15](1));}
    else if (N == 17) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1),s[8](1),s[9](1),s[10](1),s[11](1),s[12](1),s[13](1),s[14](1),s[15](1),s[16](1));}
    else if (N == 18) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1),s[8](1),s[9](1),s[10](1),s[11](1),s[12](1),s[13](1),s[14](1),s[15](1),s[16](1),s[17](1));}
    else if (N == 19) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1),s[8](1),s[9](1),s[10](1),s[11](1),s[12](1),s[13](1),s[14](1),s[15](1),s[16](1),s[17](1),s[18](1));}
    else if (N == 20) {trace+=rho.cplx(s[0](1),s[1](1),s[2](1),s[3](1),s[4](1),s[5](1),s[6](1),s[7](1),s[8](1),s[9](1),s[10](1),s[11](1),s[12](1),s[13](1),s[14](1),s[15](1),s[16](1),s[17](1),s[18](1),s[19](1));}


    else std::cout << "Warning, trace_check only possible for N<=20" << std::endl;
    return trace.real();
}


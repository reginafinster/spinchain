#include "itensor/all.h"
#include "header.h"
#include <cmath>

bool liouvillian::steady_state_check(std::vector<double>& cont, const double& mag, const double& tol)
{

    // loop through vector and compare last value with all others. Only if all of them are equal (within a certain tolerance), we have reached a steady state.
    unsigned int N = cont.size();
    int reached_steady_state = 1;
    if (cont[N-1] > 0)
    {
        for (unsigned int k = 2; k<=N; k++)
        {
            if (!(( cont[N-k] > (cont[N-1] - (cont[N-1]*tol)) ) && ( cont[N-k] < (cont[N-1] + (cont[N-1]*tol)) ))) {reached_steady_state = 0;}
        }
    }
    else
    { 
        for (unsigned int k = 2; k<=N; k++)
        {
            if (!(( cont[N-k] < (cont[N-1] - (cont[N-1]*tol)) ) && ( cont[N-k] > (cont[N-1] + (cont[N-1]*tol)) ))) {reached_steady_state = 0;}
        }
    }
    return reached_steady_state;
}

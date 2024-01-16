#include "itensor/all.h"
#include "header.h"
#include <complex>

using namespace itensor;

void liouvillian::set_hamilton(ITensor& l_tensor, const double& h_one, const double& h_two, const double& J, const double& Delta)
{

    //get indices
    IndexSet s = l_tensor.inds();
    
    //set J-coupling
    l_tensor.set(s[0](1),s[1](2),s[2](3),s[3](4), -2.0*J);
    l_tensor.set(s[0](1),s[1](3),s[2](2),s[3](4), +2.0*J);
    l_tensor.set(s[0](1),s[1](4),s[2](2),s[3](3), -2.0*J);
    l_tensor.set(s[0](1),s[1](4),s[2](3),s[3](2), +2.0*J);
    l_tensor.set(s[0](2),s[1](1),s[2](4),s[3](3), -2.0*J);
    l_tensor.set(s[0](2),s[1](3),s[2](1),s[3](4), +2.0*J);
    l_tensor.set(s[0](2),s[1](3),s[2](4),s[3](1), -2.0*J);
    l_tensor.set(s[0](2),s[1](4),s[2](1),s[3](3), -2.0*J);
    l_tensor.set(s[0](3),s[1](1),s[2](4),s[3](2), +2.0*J);
    l_tensor.set(s[0](3),s[1](2),s[2](1),s[3](4), -2.0*J);
    l_tensor.set(s[0](3),s[1](2),s[2](4),s[3](1), +2.0*J);
    l_tensor.set(s[0](3),s[1](4),s[2](1),s[3](2), +2.0*J);
    l_tensor.set(s[0](4),s[1](1),s[2](2),s[3](3), +2.0*J);
    l_tensor.set(s[0](4),s[1](1),s[2](3),s[3](2), -2.0*J);
    l_tensor.set(s[0](4),s[1](2),s[2](3),s[3](1), -2.0*J);
    l_tensor.set(s[0](4),s[1](3),s[2](2),s[3](1), +2.0*J);
    
    //set Delta-coupling 
    l_tensor.set(s[0](1),s[1](2),s[2](4),s[3](3), +2.0*Delta);
    l_tensor.set(s[0](1),s[1](3),s[2](4),s[3](2), -2.0*Delta);
    l_tensor.set(s[0](2),s[1](1),s[2](3),s[3](4), +2.0*Delta);
    l_tensor.set(s[0](2),s[1](4),s[2](3),s[3](1), +2.0*Delta);
    l_tensor.set(s[0](3),s[1](1),s[2](2),s[3](4), -2.0*Delta);
    l_tensor.set(s[0](3),s[1](4),s[2](2),s[3](1), -2.0*Delta);
    l_tensor.set(s[0](4),s[1](2),s[2](1),s[3](3), +2.0*Delta);
    l_tensor.set(s[0](4),s[1](3),s[2](1),s[3](2), -2.0*Delta);
    
    //set disorder
    for (int j = 1; j<=4; j++)
    {
        l_tensor.set(s[0](2),s[1](j),s[2](3),s[3](j), +4.0*h_one);
        l_tensor.set(s[0](3),s[1](j),s[2](2),s[3](j), -4.0*h_one);
        l_tensor.set(s[0](j),s[1](2),s[2](j),s[3](3), +4.0*h_two);
        l_tensor.set(s[0](j),s[1](3),s[2](j),s[3](2), -4.0*h_two);
    }
    
    return;
}
    

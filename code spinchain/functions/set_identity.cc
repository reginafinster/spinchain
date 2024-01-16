#include "itensor/all.h"
#include "header.h"

using namespace itensor;

void liouvillian::set_identity(ITensor& my_identity)
{
    // check if s1 and s2 have the same dimension
    if (my_identity.r() != 2) {std::cout << "Tensor is not a matrix!" << std::endl;}
    IndexSet inds = my_identity.inds();
    // check if matrix is square
    if (inds[0].m() != inds[1].m()){std::cout << "Matrix is not square!" << std::endl;}
    
    for(int i=0; i<(inds[0].m()); i++)
    {
        my_identity.set(inds[0](i+1), inds[1](i+1), 1.0);   
    }
    
    return;
}

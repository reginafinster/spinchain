#include "itensor/all.h"
#include "header.h"

void liouvillian::MultiplyTwoMPO(MPO const& Aorig, MPO const& Borig, MPO& res, Args args)
    {


    if(!args.defined("Cutoff")) args.add("Cutoff",1E-14);

    if(Aorig.N() != Borig.N()) Error("MultplyTwoMPO: Mismatched Number of sites");
    const int mposize = Borig.N();

    auto A = Aorig;
    A.position(1);

    MPO B;
    if(&Borig == &Aorig)
        {
        B = A;
        }
    else
        {
        B = Borig;
        B.position(1);
        }

    B.primeall();

    res=A;
    res.primelinks(0,4);
    res.mapprime(1,2,Site);

    ITensor clust,nfork;
    //function to get the right index
    auto siteistypeSiteandnoprime = [](Index const& i) { return (i.type() == Site && i.primeLevel()==0); };
    for(int i = 1; i < mposize; ++i)
        {
        if(i == 1)
            {
            clust = A.A(i) * B.A(i);
            }
        else
            {
            clust = nfork * A.A(i) * B.A(i);
            }

        if(i == mposize-1) break;

        nfork = ITensor(linkInd(A,i),linkInd(B,i),linkInd(res,i));

        denmatDecomp(clust,res.Anc(i),nfork,Fromleft,args);

        auto mid = commonIndex(res.A(i),nfork,Link);
        mid.dag();


        //here we included bin[i+1] instead of res.sites()(i+1) such that we dont need a siteset and can use our own indexcontainer  we need to be careful if we multiply MPOs which has bins which are not of type Site and not at least one with primelevel 0 than this function will produce an error
         res.Anc(i+1) = ITensor(mid,dag(findindex(res.A(i+1),siteistypeSiteandnoprime)),prime(findindex(res.A(i+1),siteistypeSiteandnoprime),2),rightLinkInd(res,i+1));
        }

    nfork = clust * A.A(mposize) * B.A(mposize);

    res.svdBond(mposize-1,nfork,Fromright);
    res.noprimelink();
    res.mapprime(2,1,Site);
    res.orthogonalize();

    }

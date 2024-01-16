#include "itensor/all.h"
#include "header.h"
#include <time.h>
#include <unistd.h>

void liouvillian::set_disorder(std::vector<double>& h, std::vector<double>& check_average, const double& disorder)
{
    int N = h.size();
    
    if (disorder == 0.0)
    {
        for (int i=0; i<N; i++)
        {
            h.at(i) = 0.0;
        }
    }
    else
    {
        printf("set random disorder");
        
        //seed for random number generator
        sleep(3);
        double temp;
        //srand48 (time(NULL));
        std::cout << time(NULL) << std::endl;
        srand48((long int)time(NULL));

        
        for(int i=0; i<N; i++)
        {
            temp = drand48(); // creates random number in [0.0,1.0)
            // sets h[i] in [-disorder, disorder], as drand does not cover the negetive part of the interval [-1,1]
            // we need to map [0,1) to the interval [-1,1).
            h[i]  =   ((2.0 * disorder * temp) - disorder); 
            printf("disorder at site i: %2f\n", h[i]);
            check_average[i] += h[i]; // stores all random values for each site for checking their average values
        }
    }    
    return;
}





//     else if (random==0)
//     {
//         printf("set predefinied disorder");
//         double temp[20] = {0.945103, 0.044628, 0.142511, 0.368976, 0.012486, 0.246508, 0.976066, 0.903735, 0.984476, 0.907800, 0.326360, 0.278685, 0.040283, 0.700396, 0.640707, 0.373804, 0.500995, 0.033604, 0.916084, 0.277847};
//         for(int i=0; i<N; i++)
//         {
//             h[i]  =   ((2.0 * disorder * temp) - disorder); 
//             printf("disorder at site i: %2f\n", h[i]);
//             check_average[i] += h[i]; // stores all random values for each site for checking their average values
//         }
//     }

//     else
//     {
//         for (int i=0; i<N; i++)
//         {
//             if (i%2 == 0) {h.at(i) = disorder;}
//             else {h.at(i) = (-1.0)*disorder;}
//             
//         }
//     }

#include "Riostream.h"

void errors_calc(){

    Double_t X[5][5];
    Int_t states[5] = {4274, 4685, 4630, 4500, 4700};

    ifstream in;
    in.open("results_short.txt");
    for(int i=0; i<25; i++){
        in >> X[i/5][i%5];
    }

    for(int i=0; i<5; i++){
        Double_t err = 0;
        for(int j=1; j<5; j++)
            err+= pow(X[i][j]/X[i][0]-1,2);
        printf("Extra contribution into the systematic error for X(%i) = %.4f%\n", states[i], 100*err);
        // printf("Updated systematic error for X(%i) = %.3f%\n", states[i], 100*sqrt(pow(0.04087,2)+pow(err,2)));
        
    }
}

// {
//     Int_t states[5] = {4274, 4685, 4630, 4500, 4700};
//     Double_t UL_mean[5] = {0.771, 0.772, 0.798, 0.778, 0.795};
//     Double_t UL_variation[5][4] = {{0.771, 0.771, 0.771, 0.771},
//                                    {0.772, 0.772, 0.773, 0.772},
//                                    {0.800, 0.797, 0.801, 0.795},
//                                    {0.778, 0.778, 0.779, 0.778},
//                                    {0.795, 0.794, 0.796, 0.793}};

//     for(int i=0; i<5; i++){
//         Double_t err =0;
//         for(int j=0; j<4; j++)
//             err+= pow(UL_variation[i][j]/UL_mean[i]-1,2);
//         err = sqrt(err);
//         // printf("Extra contribution into the systematic error for X(%i) = %.3f%\n", states[i], 100*err);
//         printf("Updated systematic error for X(%i) = %.3f%\n", states[i], 100*sqrt(pow(0.04087,2)+pow(err,2)));
//     }

// }
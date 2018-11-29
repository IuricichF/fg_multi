#include "forman/formangradient.h"

using namespace std;

int main(int argc, char* argv[])
{

    Timer time;
    MemoryUsage mem;


    //reading the input
    FormanGradient grad = FormanGradient(argc,argv);
    mem.getValue_in_MB(true);

    time.start();
    //computing the Forman gradient
    grad.computeFormanGradient(true);
//    grad.computeFormanGradientAllili(true);
    mem.getValue_in_MB(true);

    time.stop();
    cout << "Computation ended " << time.getElapsedTime() << endl;


   SUMap indexSimplices;
   // grad.printBoundaryMatricesSimplicialComplex(indexSimplices);
   grad.printBoundaryMatrices(indexSimplices);



    return 0;
}

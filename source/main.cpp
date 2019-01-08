#include "forman/formangradient.h"

using namespace std;

int main(int argc, char* argv[])
{

    Timer time;
    MemoryUsage mem;


    //reading the input
    time.start();
    FormanGradient grad = FormanGradient(argc,argv);
    time.stop();
    cout << "Simplicial complex loaded in " << time.getElapsedTime() << " seconds " << endl;
    mem.getValue_in_MB(true);


    time.start();
    //computing the Forman gradient
    grad.computeFormanGradient(true);
    time.stop();
    cout << "Forman gradient computed in " << time.getElapsedTime() << " seconds " << endl;
    mem.getValue_in_MB(true);


   grad.printBoundaryMatrices();



    return 0;
}

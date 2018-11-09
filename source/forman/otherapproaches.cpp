#include "formangradient.h"



void FormanGradient::computeFormanGradientAllili(bool computeTopsByBatch){


    Timer time;
    MemoryUsage mem;

    time.start();
    if(computeTopsByBatch)
        sc.storeFullStar();
    time.stop();

    cout << "Tops computed " << time.getElapsedTime() << endl;

    mem.getValue_in_MB(true);

    auto foo = bind(&FormanGradient::topologicalSorting, this,_1,_2);
    //insert all simplices in a ordered set
    SSet orderedSimplices = SSet(foo);

    for(int v=0; v<sc.getVerticesNum(); v++){
        orderedSimplices.insert(implicitS(v));

        vector<explicitS>* top = sc.topStar(explicitS(0,v));

        for(int dim=1; dim<=2; dim++){

            for(uint i=0; i<top->size(); i++){

                if((*top)[i].getDim() == dim){
                    orderedSimplices.insert(sc.toImplicit((*top)[i]));
                }
                else{

                    vector<implicitS>* bd = sc.boundaryk(sc.toImplicit((*top)[i]),dim);

                    for(vector<implicitS>::iterator it=bd->begin(); it!=bd->end(); it++){
                        orderedSimplices.insert(*it);
                    }

                    delete bd;
                }
            }
        }
    }

    mem.getValue_in_MB(true);

    SBMap classified = SBMap(foo);
    for(implicitS simpl : orderedSimplices){
        classified[simpl]=false;
    }

    mem.getValue_in_MB(true);

    time.start();
    auto simplSetFunction = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
    for(implicitS simpl : orderedSimplices){

        if(classified[simpl])
            continue;

        SSet lowerStar = SSet(simplSetFunction);
        simplexLowerStarAllili(simpl,lowerStar,classified);

        lowerStar.insert(simpl);

        homotopy_expansion(lowerStar);

        for(auto s : lowerStar){
            classified[s]=true;
        }
    }
    time.stop();


    cout << "Effective computation with Allili " << time.getElapsedTime() << endl;
    mem.getValue_in_MB(true);
 }


void FormanGradient::simplexLowerStarAllili(implicitS simpl,SSet& lwStars,SBMap& classified){

    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);

    vector<uint> simplFiltr = simplexFiltration(simpl);
    vector<int> vertices = simpl.getConstVertices();

    vector<SSet> lwPerVertex;
    for(int v : vertices)
    {
        SSet lw = SSet(foo);
        vertexLowerStarAllili(v,lw);
        lwPerVertex.push_back(lw);
    }

    for(implicitS s : lwPerVertex[0]){

        if(s.getDim() > simpl.getDim()){
            int i=1;
            for(; i<lwPerVertex.size(); i++){
                if(lwPerVertex[i].find(s) == lwPerVertex[i].end()){
                    break;
                }
            }

            if(i == lwPerVertex.size()){

                vector<uint> starFiltr = simplexFiltration(s);

                int f=0;
                for(; f<simplFiltr.size(); f++){
                    if(simplFiltr[f] < starFiltr[f])
                        break;
                }

                if(f == simplFiltr.size() && !classified[s]){
                    lwStars.insert(s);
                }
            }
        }
    }

}


void FormanGradient::vertexLowerStarAllili(int vertex,SSet& lwStars){

    vector<explicitS>* top = sc.topStar(explicitS(0,vertex));
    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
    lwStars = SSet(foo);

    for(uint i=0; i<top->size(); i++){
        explicitS sTop = (*top)[i];
        for(uint d=1; d<sTop.getDim();d++)
        {
            vector<implicitS>* bd = sc.boundaryk(sc.toImplicit((*top)[i]),d);

            for(vector<implicitS>::iterator it=bd->begin(); it!=bd->end(); it++){

                vector<int> vertices = it->getConstVertices();
                set<int> setV(vertices.begin(), vertices.end());

                if(setV.find(vertex) != setV.end()){
                    lwStars.insert(*it);
                }
            }

            delete bd;
        }

        lwStars.insert(sc.toImplicit(sTop));
    }

    delete top;

}

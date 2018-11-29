#include "formangradient.h"


void FormanGradient::printBoundaryMatrices(SUMap& indexSimplexes){

    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
    SUMap criticalCells = SUMap(foo);
    uint index=0;
    vector<map<float,uint> > fieldValues = vector<map<float,uint> >(nFields());

    for(auto lvl : criticalS){
        for(auto s : lvl.second){

            vector<float> fVal = simplexScalarValue(s);
            for(int i=0; i<fieldValues.size(); i++){
                if(fieldValues[i].find(fVal[i]) == fieldValues[i].end()){
                    fieldValues[i][fVal[i]]=0;
                }
            }
        }
    }

    ofstream out("boundaryForman.txt");

    out << fieldValues.size() << endl;
    //indices of the filtration steps
    for(int i=0; i<fieldValues.size(); i++){
        int index=0;
        for(auto it : fieldValues[i]){
            fieldValues[i][it.first]=index++;
            out << it.first << " ";
        }
        out << endl;
    }

    //indices of the cells
    index=0;
    for(auto lvl : criticalS){
        for(auto s : lvl.second){
            criticalCells[s]=index++;
        }
    }

    Timer time;

    time.start();
    cout << "Critical cells: " << index << endl;
    for(auto lvl : criticalS){
        for(auto s : lvl.second){
            SSet cells;

            if(lvl.first != 0)
                computeBoundaryCell(s, cells);


            // if(indexSimplexes.size() != 0){
            //     assert(indexSimplexes.find(s) != indexSimplexes.end());
            //     out << "[" << indexSimplexes[s] << "] ";
            // }
            // else
            //     out << "[-1] ";
            //
            // out << lvl.first << " ";
            // for(auto c : cells)
            //     out << criticalCells[c] << " ";
            // out << ": ";
            //
            // vector<float> fVal = simplexScalarValue(s);
            // for(int i=0; i<fVal.size(); i++){
            //     out << fieldValues[i][fVal[i]] << " ";
            // }
            // out << endl;
        }
    }

    time.stop();
    cout << time.getElapsedTime() << endl;

    out.close();
}


void FormanGradient::printBoundaryMatricesSimplicialComplex(SUMap& criticalCells){

    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
    criticalCells = SUMap(foo);
    uint index=0;
    vector<map<float,uint> > fieldValues = vector<map<float,uint> >(nFields());

    map<uint, SSet > incidenceGraph;

    for(int i=0; i<=sc.getComplexDim(); i++){
        incidenceGraph[i]=SSet(foo);

        if(i==0){
            for(int j=0; j<sc.getVerticesNum(); j++){
                implicitS vert = sc.toImplicit(explicitS(0,j));
                incidenceGraph[i].insert(vert);
            }
        }
        else{

            for(int index=0; index<sc.getTopSimplexesNum(i); index++){
                implicitS simpl = sc.toImplicit(explicitS(i,index));
                incidenceGraph[i].insert(simpl);
            }

            for(int k=i+1; k<=sc.getComplexDim(); k++){

                for(int index=0; index<sc.getTopSimplexesNum(k); index++){
                    implicitS simpl = sc.toImplicit(explicitS(k,index));
                    vector<implicitS>* bds = sc.boundaryk(simpl,i);
                    for(auto s : *bds){
                        incidenceGraph[i].insert(s);
                    }
                    delete bds;
                }
            }

        }

    }

    for(auto lvl : incidenceGraph){
        for(auto s : lvl.second){

            vector<float> fVal = simplexScalarValue(s);
            for(int i=0; i<fieldValues.size(); i++){
                if(fieldValues[i].find(fVal[i]) == fieldValues[i].end()){
                    fieldValues[i][fVal[i]]=0;
                }
            }
        }
    }

    ofstream out("boundarySimplex.txt");

    out << fieldValues.size() << endl;
    //indices of the filtration steps
    for(int i=0; i<fieldValues.size(); i++){
        int index=0;
        for(auto it : fieldValues[i]){
            fieldValues[i][it.first]=index++;
            out << it.first << " ";
        }
        out << endl;
    }

    //indices of the cells
    index=0;
    for(auto lvl : incidenceGraph){
        for(auto s : lvl.second){
            criticalCells[s]=index++;
        }
    }

    for(auto lvl : incidenceGraph){
        for(auto s : lvl.second){
            vector<implicitS>* bds = new vector<implicitS>;

            if(lvl.first != 0)
                bds = sc.boundaryk(s, lvl.first-1);

            out << lvl.first << " ";
            for(auto c : *bds)
                out << criticalCells[c] << " ";
            out << ": ";

            delete bds;

            vector<float> fVal = simplexScalarValue(s);
            for(int i=0; i<fVal.size(); i++){
                out << fieldValues[i][fVal[i]] << " ";
            }
            out << endl;
        }
    }

    out.close();
}

#include "formangradient.h"

bool FormanGradient::cmpSimplexesFiltr(const implicitS& lhs, const implicitS& rhs){

    if(lhs.getDim() == rhs.getDim()){
        vector<uint> fValuesL;
        vector<uint> fValuesR;

        for(int i=0; i<=lhs.getDim();i++){
            fValuesL.push_back(filtration[lhs.getConstVertices()[i]]);
            fValuesR.push_back(filtration[rhs.getConstVertices()[i]]);
        }

        sort(fValuesL.begin(),fValuesL.end(), std::greater<uint>());
        sort(fValuesR.begin(),fValuesR.end(), std::greater<uint>());

        for(uint i=0; i<fValuesL.size(); i++){
            if(fValuesL[i]==fValuesR[i])
                continue;

            return fValuesL[i] < fValuesR[i];
        }
    }

    return lhs.getDim() < rhs.getDim();
}



bool FormanGradient::sortVerticesFiltration(const int& v1,const int& v2){
    return filtration[v1] > filtration[v2];
}


bool FormanGradient::filtrComparer(const vector<float>& v1, const vector<float>& v2) const{

    for( int i=0; i<v1.size(); i++){
        if (v1[i] != v2[i]){
          return v1[i] < v2[i];
        }
    }

    return false;
}

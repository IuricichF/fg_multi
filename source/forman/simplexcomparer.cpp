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

bool FormanGradient::cmpInjectiveFiltr(const implicitS& lhs, const implicitS& rhs){

    if(lhs.getDim() == rhs.getDim()){
        vector<uint> fValuesL = simplexFiltration(lhs);
        vector<uint> fValuesR = simplexFiltration(rhs);

        int greater=0;
        int smaller=0;

        for(uint i=0; i<fValuesL.size(); i++){
            if(fValuesL[i]>fValuesR[i])
                greater++;
            if(fValuesL[i]<fValuesR[i])
                smaller++;

        }

        if(greater == 0)
            return true;
        else if(greater == smaller){
            for(uint i=0; i<fValuesL.size(); i++){
                return fValuesL[i]<fValuesR[i];
            }
        }
        else{
            return false;
        }
    }

    return lhs.getDim() < rhs.getDim();
}

bool FormanGradient::topologicalSorting(const implicitS& lhs, const implicitS& rhs){

    if(sc.theSame(lhs,rhs))
        return false;

    vector<uint> fValuesL = simplexFiltration(lhs);
    vector<uint> fValuesR = simplexFiltration(rhs);

    for(int i=0; i<fValuesL.size(); i++){
        if(fValuesL[i] != fValuesR[i]){
            return fValuesL[i] < fValuesR[i];
        }
    }

    if(lhs.getDim() != rhs.getDim()){
        return lhs.getDim() < rhs.getDim();
    }
    else{
        vector<int> vLhs = lhs.getConstVertices();
        vector<int> vRhs = rhs.getConstVertices();

        for(int i=0; i<vLhs.size(); i++){
            if(vLhs[i] != vRhs[i])
                return vLhs[i] < vRhs[i];
        }
    }
}

bool FormanGradient::sortVerticesFiltration(const int& v1,const int& v2){
    return filtration[v1] > filtration[v2];
}


bool FormanGradient::filtrComparer(const pair<float,uint>& v1, const pair<float,uint>& v2) const{
    if(v1.first == v2.first){
        return v1.second < v2.second;
    }

    return v1.first < v2.first;
}


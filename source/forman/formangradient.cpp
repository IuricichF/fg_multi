#include "formangradient.h"

FormanGradient::FormanGradient(int argc, char** argv)
{

    sc = SimplicialComplex();

    string fileMesh(argv[1]);

    if(fileMesh.find(".off") != string::npos){
        sc.readOFF(fileMesh.c_str());
    }
    else if(fileMesh.find(".bin") != string::npos || fileMesh.find(".raw") != string::npos){

        vector<int> res;
        for(int i=0; i<argc-2; i++){
            res.push_back(atoi(argv[i+2]));
        }

        if(res.size() == 2){
            sc.readSquareGrid(fileMesh.c_str(),res[0],res[1]);
        }
        else if(res.size() == 3){
            sc.readCubicalGrid(fileMesh.c_str(),res[0],res[1],res[2]);
        }
        else{
           cout << "Library works only for 2D or 3D datasets" << endl;
        }
    }
    else{
        cout << "Unknown format" << endl;
    }
    MemoryUsage mem;
    mem.getValue_in_MB(true);

    cout << "Computing fields" << endl;

    int nField = sc.getVertex(0).getCoordinates().size()-3;
    cout << "Simplicial complex dimension " << sc.getComplexDim() << endl;
    cout << "Number of fields " << nField << endl;

    //create gradient encoding
    gradient = GradientEncoding(sc);

    filtration = vector<uint>(sc.getVerticesNum(),0); //final filtration
    //cout << "Start" << endl;

    vector<vector<float> > injectiveF(sc.getVerticesNum());

    for(int i=0; i<nField; i++){
        for(int j=0;j<sc.getVerticesNum(); j++){
            float val=sc.getVertex(j).getCoordinate(i+3);

            vector<float> values = sc.getVertex(j).getCoordinates();
            values.erase(values.begin()+(i+3));
            values.erase(values.begin(), values.begin()+3);
            values.insert(values.begin(),val);
            values.push_back(j);

            //cout << val << endl;
            injectiveF[j]=vector<float>(values);
        }
        //cout << "Stop" << endl;
    }

    sort(injectiveF.begin(),injectiveF.end(),bind(&FormanGradient::filtrComparer, this,_1,_2));

    //filtration created
    int ind=0;
    for(auto vec : injectiveF){
        filtration[vec.back()]=ind++;
    }

}

FormanGradient::~FormanGradient()
{

}


void FormanGradient::computeFormanGradient(bool computeTopsByBatch){

    if(computeTopsByBatch)
        sc.storeFullStar();

    // map<pair<int,int>, SSet > filtrationAll;
    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);

    #pragma omp parallel for
    for(uint i=0; i<filtration.size(); i++){
        vector<SSet> lwStars;
        splitVertexLowerStar(i,lwStars);

        for(auto lw : lwStars){
            //uncomment here for old software
            homotopy_expansion(lw);

        }

    }

}

void FormanGradient::splitVertexLowerStar(int v,vector<SSet>& lwStars){

    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
    map<vector<float>, SSet> lws;
    implicitS vert = implicitS(v);

    vector<float> filtrLvl = simplexScalarValue(vert);


    if(lws.find(filtrLvl) == lws.end()){
        lws[filtrLvl]= SSet(foo);
    }
    lws[filtrLvl].insert(implicitS(v));


    //cout << sc.getComplexDim() << endl;
    for(int i=1; i<=sc.getComplexDim(); i++){

        SSet* lw = vertexLowerStar(v,i);
        //cout << lw->size() << endl;

        for(auto s : *lw){
            //cout << s << endl;
            filtrLvl = simplexScalarValue(s);

//            for(auto v : filtrLvl)
//                cout << v << " ";
//            cout << endl;

            if(lws.find(filtrLvl) == lws.end()){
                lws[filtrLvl]=SSet(foo);
            }
            lws[filtrLvl].insert(s);
        }
        delete lw;
    }

    for(auto l : lws)
        lwStars.push_back(l.second);
}

void FormanGradient::homotopy_expansion(SSet& simplexes){

    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
    vector<SSet> sdiv(4,SSet(foo));
    uint alive=0;

    for(auto s : simplexes){
        sdiv[s.getDim()].insert(s);
        alive++;
    }

    int d=1;

    while(alive != 0){

        if(d < sdiv.size()){

            while(sdiv[d-1].size() > 0){
                //cout << "Here again " << d << endl;

                list<implicitS> toRemove;
                for(auto s : sdiv[d]){
                    //cout << "Next simplex " << s << endl;
                    implicitS nextPair;

                    bool nmPairable = numPairableLowerStar(s,sdiv[d-1],nextPair) == 1;

                    if(nmPairable){

                        //cout << "prima di pair con " << nextPair << endl;
                        setPair(s,nextPair);

                        sdiv[d-1].erase(nextPair);
                        toRemove.push_back(s);
                        alive-=2;
                    }
                }

                if(!toRemove.empty()){
                    //cout << "Prima di rimuovere " << toRemove.size() << " " << sdiv[d].size() << endl;
//                    for(auto s : sdiv[d]){
//                        //cout << s << endl;
//                    }

                    for(auto s : toRemove){
                        //cout << s << endl;
                        sdiv[d].erase(s);
                    }
                    //cout << sdiv[d].size() << endl;
                }
                else{
                    implicitS critical=*sdiv[d-1].begin();
                    sdiv[d-1].erase(critical);
                    alive--;

                    #pragma omp critical
                    {
                        if(criticalS.find(critical.getDim()) == criticalS.end()){
                            SSet crit = SSet(foo);
                            criticalS[critical.getDim()]=crit;
                        }
                        criticalS[critical.getDim()].insert(critical);
                    }
                }
            }

            d++;

        }
        else{
            implicitS critical=*sdiv[d-1].begin();
            sdiv[d-1].erase(critical);
            alive--;

            #pragma omp critical
            {
                if(criticalS.find(critical.getDim()) == criticalS.end()){
                    SSet crit = SSet(foo);
                    criticalS[critical.getDim()]=crit;
                }
                criticalS[critical.getDim()].insert(critical);
            }
        }
    }
}


 vector<float> FormanGradient::simplexScalarValue(const implicitS &simpl){

     vector<int> vertis = simpl.getConstVertices();
     vector<float> filtr = sc.getScalarFields(vertis[0]);

     for(uint i=0; i<filtr.size(); i++){
         for(auto v : vertis){
             if(sc.getScalarField(v,i) > filtr[i])
                 filtr[i]=sc.getScalarField(v,i);
         }
     }


     return filtr;
 }

int FormanGradient::numPairableLowerStar(const implicitS &next, const SSet& sset, implicitS& pair){

    vector<implicitS>* boundary = sc.boundaryk(next,next.getDim()-1);

    int num=0;
    for(auto s : *boundary){
        if(sset.find(s) != sset.end()){
            num++;
            pair=s;
        }
    }

    delete boundary;

    return num;
}

bool FormanGradient::isPaired(const implicitS &simpl){

    vector<int> vertices = simpl.getConstVertices();
    map<uint, vector<explicitS> > tops;

    //retrieve the fan of top incident into the vertices of next
    for(auto v : vertices){
        vector<explicitS>* topSimpl = sc.topStar(implicitS(v));
        tops[v]=*topSimpl;
        delete topSimpl;
    }

    vector<explicitS> cobBig = tops[vertices[0]];
    for(int i=1; i<vertices.size(); i++){

        vector<explicitS> other = tops[vertices[i]];

        vector<explicitS> result(cobBig.size());
        vector<explicitS>::iterator it = set_intersection(cobBig.begin(), cobBig.end(),
                         other.begin(), other.end(), result.begin());
        result.resize(it-result.begin());
        cobBig=result;
    }


    if(cobBig.size() == 0){
        cout << "the coboundary was empty " << endl;
        cout << simpl.getDim() << endl;

        vector<int> verts = simpl.getConstVertices();

        for(auto v : verts)
            cout << v << " ";
        cout << endl;

        explicitS top = sc.toExplicit(simpl);
        return gradient.isPaired(simpl,top,sc);
    }

    for(auto s : cobBig){
        if(gradient.isPaired(simpl, s, sc))
            return true;
    }

    return false;
}



void FormanGradient::setPair(const implicitS &next, const implicitS &pair){
    //next has to be bigger than pair
    assert(next.getDim() > pair.getDim());

    vector<explicitS>* tops = sc.topStar(next);

    #pragma omp critical
    {
        for(auto s : *tops){
            gradient.pair(next,pair,s,sc);
        }
    }
}

void FormanGradient::freePair(const implicitS &next, const implicitS &pair){
    //next is bigger than pair

    vector<explicitS>* tops = sc.topStar(next);
    for(auto s : *tops){
        gradient.free(next,pair,s,sc);
    }
}


bool FormanGradient::getPair(const implicitS& next, implicitS& pair){

    vector<explicitS>* tops = sc.topStar(next);
    for(auto s : *tops){
        if(gradient.getPair(next,pair,s,sc)){
            return true;
        }
    }

    return false;
}


//compute the lower star of a vertex already knowing the top simplexes
SSet* FormanGradient::vertexLowerStar(uint vert, uint dimension){

    vector<explicitS>* top = sc.topStar(explicitS(0,vert));
//    cout << "Vertex " << vert << endl;
//    cout << "Size top " << top->size() << " " << dimension <<endl;
    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
    SSet* ret = new SSet(foo);

    for(uint i=0; i<top->size(); i++){
        if((*top)[i].getDim() < dimension)
            continue;

        if((*top)[i].getDim() == dimension){
            vector<int> vertices = sc.getTopSimplex((*top)[i]).getVertices();
            sort(vertices.begin(), vertices.end(), bind(&FormanGradient::sortVerticesFiltration, this,_1,_2));


            if(vertices[0]==vert){
//                for(auto v : vertices)
//                    cout << v << " ";
//                cout << endl;

                ret->insert(sc.toImplicit((*top)[i]));
            }
            continue;
        }

        vector<implicitS>* bd = sc.boundaryk(sc.toImplicit((*top)[i]),dimension);
        //cout << bd->size() << endl;

        for(vector<implicitS>::iterator it=bd->begin(); it!=bd->end(); it++){
            vector<int> vertices = it->getConstVertices();
            sort(vertices.begin(), vertices.end(), bind(&FormanGradient::sortVerticesFiltration, this,_1,_2));

            if(vertices[0]==vert){
//                for(auto v : vertices)
//                    cout << v << " ";
//                cout << endl;

                ret->insert(*it);
            }
        }

        delete bd;
    }

    delete top;

    return ret;
}

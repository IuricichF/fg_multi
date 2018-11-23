#ifndef FORMANGRADIENT_H
#define FORMANGRADIENT_H

#include <stack>
#include <math.h>

#include "gradientencoding.h"


typedef set<implicitS, boost::function<bool(const implicitS &, const implicitS &)>> SSet;
typedef map<implicitS, unsigned, boost::function<bool(const implicitS &, const implicitS &)>> SUMap;
typedef map<implicitS, bool, boost::function<bool(const implicitS &, const implicitS &)>> SBMap;


class FormanGradient
{
private:
    vector<uint> filtration; //for each vertex its filtration value
    vector<vector<float> > scalarValues; //for each vertex its field value [Vertices x number of fields]
    vector<vector<uint> > componentBasedFiltration; //injective function for each component [number of fields x Vertices ]

    GradientEncoding gradient;
    map<uint, SSet > criticalS;

    SimplicialComplex sc;


public:
    FormanGradient(int,char**);
    ~FormanGradient();

    //compute the Forman gradient on the input dataset, if batchTop=true the relation among vertices and incident top simplices
    //is computed once and stored explicitly (higher memory consumption but lower timings)
    void computeFormanGradient(bool batchTop);

    //compute the Forman gradient on the input dataset, if batchTop=true the relation among vertices and incident top simplices
    //is computed once and stored explicitly (higher memory consumption but lower timings)
    void computeFormanGradientAllili(bool);

    //print critical points and 1-skeleton
    void print1skeleton();
    void printBoundaryMatrices(SUMap& criticalCells);
    void printBoundaryMatricesSimplicialComplex(SUMap &criticalCells);

private:
    SSet* vertexLowerStar(uint vert, uint d); //compute the lower star of vert (only simplices of dimension=d)
    void splitVertexLowerStar(int v,vector<SSet>& lwStars); //split the lower star of a vertex v according to the sublevelsets of the function
    void homotopy_expansion(SSet&); //apply homotopy expansion on a set of simplices
    int numPairableLowerStar(const implicitS &next, const SSet& sset, implicitS& pair); //return the number of simplices pairable with next in the lower star, pair is one of these
    bool isPaired(const implicitS &simpl); //true if simpl is paired with another simplex
    void setPair(const implicitS &next, const implicitS &pair); //set the new gradient pair between next and pair (NOTE: next has to be bigger than pair)
    bool getPair(const implicitS &simpl, implicitS& next); //next is the simplex paired with simpl
    void freePair(const implicitS &next, const implicitS &pair); //remove pair (next,pair) from the gradient (NOTE: next has to be bigger than pair)

    vector<uint> simplexFiltration(const implicitS &simpl); //return the vector-valued filtration for a simplex. Each component is obtained as the maximum of the filtrations of its vertices
    vector<float> simplexScalarValue(const implicitS &simpl); //return the vector-valued function for a simplex. Each component is obtained as the maximum of the function values of its vertices

    void computeBoundaryCell(implicitS const& cell, SSet& descCells); //given a critical i-simplex return the critical (i-1)-cells connected to it


    //Simplex comparer - to organize simplices in maps and sets based on their function/filtration values
    bool cmpSimplexesFiltr(const implicitS& lhs, const implicitS& rhs); //true if the indexing of lhs is less than rhs (or if lhs is smaller than rhs). Used for homotopy expansion when working on a single sublevelset
    bool cmpInjectiveFiltr(const implicitS& lhs, const implicitS& rhs); //true if the injective filtration of lhs is less that rhs (or if lhs is smaller than rhs)
    bool sortVerticesFiltration(const int& v1,const int& v2); //true if the filtration of v1 is less than v2
    bool filtrComparer(const vector<float>& v1, const vector<float>& v2) const; //true if v1 has a smaller value scalar value (float) or a smaller index (uint)
    bool topologicalSorting(const implicitS& lhs, const implicitS& rhs); //true if the filtration of lhs is lower than rhs, or if lhs is smaller than rhs

    //comparison with the algorithm by [Allili et al. 2016]
    void vertexLowerStarAllili(int vertex,SSet& lwStars);
    void simplexLowerStarAllili(implicitS simpl,SSet& lwStars,SBMap& classified);

};

#endif // FORMANGRADIENT_H

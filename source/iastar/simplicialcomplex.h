﻿#ifndef SIMPLICIALCOMPLEX_H
#define SIMPLICIALCOMPLEX_H

#include <vector>
#include <forward_list>
#include <queue>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <list>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include "topsimplex.h"
#include "vertex.h"



class SimplicialComplex
{


protected:

    map<int,int> realIndex;
    vector<vector<TopSimplex> > topSimplexes;
    vector<forward_list<int> > adjRelations;

    vector<Vertex> vertices;
    vector< vector<explicitS> > topPerVertex; //empty by default. Call storeFullStart() to encode the VTop relations and boost computation speed.


    //EmbeddingSpace == number of coordinates for each vertex
    //SimplicialComplexSpace == dimension of the biggest top simplex

public:
    SimplicialComplex();

    //------------------------------------------------------------------------------------------
    //functions for reding input files. Supported formats [OFF,GMV]
    void readOFF(const char* file); //read .OFF file
    void readTS(const char* file); //read .TS  file
    void readGMV(const char* file); //read .GMV file
    void readSimplexes(char* file); //TODO //set of simplexes. Vertices have no coordinates
    void readPoints(char* name_in_file, float threshold, uint dimension); //set of points on which a Vietoris-Rips complex is computed
    void readGraph(char* name_in_file); //set of arcs of a graph based on which clicks are computed
    void readIA(char* file); //read IA file format
    void readSquareGrid(const char* file, int xRes, int yRes);
    void readCubicalGrid(const char* file, int xRes, int yRes, int zRes);

    void saveIA(char* file); // write IA file
    void saveVTK(char* file); // write VTK file

    void storeFullStar();
    void emptyFullStar();
    //explicitly store/remove the complete set of top-simplexes for each vertex (speedup computation bur requires more memory)

    int getVerticesNum();
    //return: the number of vertices

    int getTopSimplexesNum(int d);
    //return: the number of top-simplexes of dimension d

    int getComplexDim();
    //return: the simplicial complex dimension (dimension of its maximal simplex)

    vector<int> getTopSimplexesSet();
    //return: return the number of top-simplexes of each dimension (dimensions having 0 top-simplexes are ignored)

    Vertex& getVertex(int i);
    //return: vertex of index i

    vector<float> getScalarFields(int i);
    //return: scalar values for vertex i

    float getScalarField(int i, int j);
    //return scalar value j for vertex i

    vector<TopSimplex>& getTopSimplices(int d);
    //return: set of top-simplexes of dimendion d

    TopSimplex& getTopSimplex(const explicitS);
    //return: top-simplex corresponding to the input simplex encoded in the explicit way

    //------------------------------------------------------------------------------------------
    //functions for simplexes explicitly encoded in the IA* data structure (explicitS)

    implicitS toImplicit(const explicitS&);
    //return: implicit representation of a simplex computed from its explicit representation

    vector<implicitS>* link(const explicitS&);
    //return: the set of simplexes in the link of s.
    //NOTE: only the simplexes of highets dimension are returned, i.e. if edge e is in the link of s, e is returned by the function but its vertices are not.

    vector<implicitS>* boundaryk(const explicitS& s,uint k);
    //return: simplexes of dimension k on the boundary of s

    vector<implicitS>* coboundaryk(const explicitS& s,uint k);
    //return: simplexes of dimension k on the coboundary of s

    vector<explicitS>* topAdjacent(const explicitS& s, uint i);
    //return: simplexes adjacent to s and sharing the i-th face of s

    vector<explicitS>* topStar(const explicitS& s);
    //return: top-simplexes on the coboundary of s (or, in other words, incident to s)
    //NOTE: valid only for vertices (a top-simplex cannot have simplexes on their coboundary)


    vector<explicitS>* topStar(const explicitS& s, int d);
    //return: top-simplexes of dimension d on the coboundary of s (or, in other words, incident to s)
    //NOTE: valid only for vertices (top-simplex cannot have simplexes on their coboundary)


    //------------------------------------------------------------------------------------------
    //functions for simplexes not encoded in the IA* and implicitly represented (implicitS) [see implicitS.h for details]

    explicitS toExplicit(const implicitS&);
    //return: explicit representation of a simplex computed from its implicit representation
    //NOTE: only a vertex or a top simplex have a valid explicit representation;

    bool theSame(const implicitS& s1,const implicitS& s2);
    //return: true if s1 and s2 represent the same simplex
    //NOTE: use theSame() when you do not know if s1 and s2 are represented based on the same top-simplex. Otherwise operator== is faster. [see implicitS.h for details]

    vector<implicitS>* link(const implicitS& s);
    //return: the set of simplexes in the link of s.
    //NOTE: only the simplexes of highets dimension are returned, i.e. if edge e is in the link of s, e is returned by the function but its vertices are not.

    vector<implicitS>* boundaryk(const implicitS& s,uint k);
    //return: simplexes of dimension k on the boundary of s

    vector<implicitS>* coboundaryk(const implicitS&,uint);
    //return: simplexes of dimension k on the coboundary of s

    vector<implicitS>* adjacents(const implicitS& s);
    //return: simplexes adjacent to s

    vector<explicitS>* topStar(const implicitS& s);
    //return: top-simplexes on the coboundary of s (or, in other words, incident to s)

    bool areIncident(const implicitS& small, const implicitS& big);
    //return: true if "small" is on the boundary of "big"

    //------------------------------------------------------------------------------------------
    //functions for geometry, vectors, etc
    Vertex barycenter(const implicitS& s);
    Vertex barycenter(const explicitS& s);
    vector<float> vectorSum(const vector<float> &v1, const vector<float> &v2);
    vector<float> vectorSubtr(const vector<float> &v1, const vector<float> &v2);
    vector<float> scalarProd(const vector<float> &v1, const vector<float> &v2);
    vector<float> prod(const vector<float> &v1, float prod);
    float dotProd(const vector<float> &v1, const vector<float> &v2);
    float dist(const vector<float> &v1, const vector<float> &v2);
    float norm(const vector<float> &v1);

protected:
    void buildDataStructure();
    //function for initializing the IA*

    implicitS toImplicitInsideTop(const explicitS vertex, const explicitS simpl);
    //return: implicit representation of a vertex (inside top simplex simpl) computed from its explicit representation

    vector<implicitS>* inTop(const explicitS& topSimpl, uint dim);

    uint32_t nextSubset(uint32_t v);
    vector<int> filter(const vector<int>& v, uint32_t mask);
    vector<implicitS>* getSubsets(const vector<int>& arr, uint32_t k);


    forward_list<explicitS>* incidentCluster(const explicitS v, const explicitS s);
    //(given a vertex v, and a top simplex of dimension k incident in v, retrieve all the k-1 connected simplexes still incident in v)

    void recursive_insert(implicitS simplex, const implicitS& original,uint pos, uint dim, forward_list<implicitS>* ret);
    //recursive function used for computing boundary and coboundary simplexes.

    void build_top_simplex(set<uint> setR, set<uint> setP, set<uint> setX, const vector<set<uint>* >& arcs, map<int, list<TopSimplex>* >* top_simplexes_local);
    //Function for computing the top simplexes given a graph


};


#endif // SIMPLICIALCOMPLEX_H

/* Written by : Humayun Kabir
 *          humayun.k1@gmail.com
 */

#ifndef _GRAPH_COARSENING_H
#define _GRAPH_COARSENING_H

#include<iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <vector>
#include <cassert>
#include <set>
#include <map>

#include "graph.cpp"

using namespace std;


struct neighborWeight {

        unsigned int neighbor;
        unsigned int weight;

        neighborWeight(unsigned int _neighbor, unsigned int _weight) : neighbor(_neighbor), weight(_weight) {
        }

};

template<typename T1, typename T2>
struct sortBySecond {
        typedef pair<T1, T2> typePair;
        bool operator () (typePair const& a, typePair const& b) const {
                return a.second < b.second;
        }
};
 


//Find a maximal matching by visiting vertices in random order and matching a vertex to one of neighbors randomly
void randomMatching(graph_t g, unsigned int *perm, graph_t & coarsenedG, int target);

//Find a maximal matching by visiting vertices in random order and matching a vertex the neighbor with maximum weight
void heavyEdgeMatching(graph_t g, unsigned int *perm, graph_t & coarsenedG, int target);

//Find a maximal matching by visiting vertices in random order and matching a vertex the neighbor with least weight
void lightEdgeMatching(graph_t g, unsigned int *perm, graph_t & coarsenedG);

//Find final mapping from the coarsest vertices to finest vertices
void findFinalMapping(int numSupRows, vector<unsigned int> & vectorSizes, vector<unsigned int *> & mappingVectors, unsigned int* & numEdgesSupRowsToRows, unsigned int* & mapSupRowstoRows);

//Coarsen a graph g and finds a graph coarsendGraph such that each vertex of coarsendGraph has supRowSize vertices of g
void coarsenGraph(graph_t & g, int k, unsigned int* & numEdgesSupRowsToRows, unsigned int* & mapSupRowstoRows, graph_t & coarsendGraph);



#endif


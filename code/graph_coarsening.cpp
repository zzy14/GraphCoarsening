/* Written by : Humayun Kabir
 *          humayun.k1@gmail.com
 */

#include "graph_coarsening.h"


/* Coarsen a graph and find a coarsend graph until a vertex in coarsened graph has supRowSize vertices 
 * of the input graph
 * Input: g - input graph
 * 	  supRowSize - super row size
 * Output: coarsendGraph - coarsened graph
 * 	   numEdgesSupRowsToRows - points to mapSupRowstoRows to indicate where each coarsened vertex starts
 * 	   mapSupRowstoRows - contains fine graph vertices
 */
void coarsenGraph(graph_t & g, int supRowSize, unsigned int* & numEdgesSupRowsToRows, unsigned int* & mapSupRowstoRows, graph_t & coarsendGraph) {
	
	vector<unsigned int *> mappingVectors;	
	vector<unsigned int> vectorSizes;	

	graph_t localG = g;

	int finalNumVtxs = localG.n / supRowSize;

	//default value is g.n
	//unsigned int maxSupRowSize = supRowSize * alpha;

	unsigned int *perm;
	while(  finalNumVtxs < localG.n ){ 
	
		perm = (unsigned int *) malloc(localG.n * sizeof(unsigned int));

		randomMatching(localG, perm, coarsendGraph);
		//heavyEdgeMatching(localG, perm, coarsendGraph);
		//lightEdgeMatching(localG, perm, coarsendGraph);

		vectorSizes.push_back( localG.n);
		mappingVectors.push_back(perm);
		localG = coarsendGraph;
	}

	//Find the final mapping from the coarsest graph to the fine graph 
	findFinalMapping(localG.n, vectorSizes, mappingVectors, numEdgesSupRowsToRows, mapSupRowstoRows);

        for(int i = 0; i < mappingVectors.size(); i++)
                if( mappingVectors[i] != NULL )
                        free( mappingVectors[i] );

}

/* Given intermediate mappinf of vertices find the final mapping from the coarsest vertices to the finest vertices
 * Input: numSupRows - number of coarsest vertices
 * 	  vectorSizes - vector containing number of vertices in intermediate coarsened graphs
 *	  mappingVectors - vector containing the mapping of fine vertices to coarse vertices 
 * Output: numEdgesSupRowsToRows - points to mapSupRowstoRows to indicate where each coarsest row starts
 * 	   mapSupRowstoRows - contains the vertices of the finest graph
 */

void findFinalMapping(int numSupRows, vector<unsigned int> & vectorSizes, vector<unsigned int *> & mappingVectors, unsigned int* & numEdgesSupRowsToRows, unsigned int* & mapSupRowstoRows) {

	vector< vector<unsigned int> > finalAdjLists;
	vector< vector<unsigned int> > tempVectAdjLists;

	for(int i = vectorSizes.size() -1; i >= 0; i--) {
	
		if( i == vectorSizes.size() -1 ) {
			finalAdjLists.resize( numSupRows );

			for(int j = 0; j < vectorSizes[i]; j++) 
				finalAdjLists[ mappingVectors[i][j] ].push_back( j );
		}
		else {
			tempVectAdjLists.resize(vectorSizes[i+1]);
			for(int j = 0; j < vectorSizes[i+1]; j++) 
				tempVectAdjLists[ j ].resize(0);
			

			for(int j = 0; j < vectorSizes[i]; j++) 
				tempVectAdjLists[ mappingVectors[i][j] ].push_back( j );
		
			for(int j = 0; j < numSupRows; j++) {
				vector<unsigned int> curAdjList;
				for(int k = 0; k < finalAdjLists[j].size(); k++) {
					int finerRowIndex = finalAdjLists[j][k];
					for(int l = 0; l < tempVectAdjLists[finerRowIndex].size(); l++) 
						curAdjList.push_back( tempVectAdjLists[finerRowIndex][l] );
				}
				finalAdjLists[j] = curAdjList;
			}
		}

	}


	numEdgesSupRowsToRows = (unsigned int *)malloc( (numSupRows + 1) * sizeof(unsigned int));

	numEdgesSupRowsToRows[0] = 0;
	int mapIndex = 0;
	for(int i = 0; i < numSupRows; i++) {
		numEdgesSupRowsToRows[i+1] = numEdgesSupRowsToRows[i] + finalAdjLists[i].size();
		for(int j = 0; j < finalAdjLists[i].size(); j++) {
			mapSupRowstoRows[mapIndex] = finalAdjLists[i][j];
			mapIndex++;
		}
	}	

}


/* randomMatching - Find a maximal matching of the vertcies. Visit a vertex u in random
 * order and match u randomly with one of its unmatched neighbors
 * Input: g - the graph with edge and vertex weight
 * Output: perm - mapping of vertices of g to the vertices of coarsendG
 *         coarsenedG - graph of the coarsened vertices */


void randomMatching(graph_t g, unsigned int *perm, graph_t & coarsenedG) {

	//Matched edges   
        set<Edge> matchEdges; 

	//Initially vertices are not matched 
        bool isMatched[g.n];

	srand ( unsigned ( std::time(0) ) );
	vector<int> myvector;

	//Initialize these values
	for (int i = 0; i < g.n; ++i) {
                isMatched[i] = false; 
		myvector.push_back(i);
        }
 
	//Randomly shuffles the vertices
	std::random_shuffle ( myvector.begin(), myvector.end() );

	//Visit a vertex at random and randomly match it with one of its neighbors
	for (int i = 0; i < g.n; ++i) {
                int currVtx = myvector[i];

                if( !isMatched[currVtx] ) {
			//Make a vector of unmatched neighbors	
			vector<unsigned int> neighbors;
			
			for(int j = g.num_edges[currVtx]; j < g.num_edges[currVtx+1]; j++) {
				if( (g.adj[j] != currVtx) && !isMatched[ g.adj[j] ] ) {
					neighbors.push_back( g.adj[j] );
				}
			}

			//Randomly select one of its unmatched neighbors	
			int numVtxs = neighbors.size();
			if( numVtxs == 1) {	
				isMatched[ currVtx ] = true;
                                isMatched [ neighbors[0] ] = true;
                                matchEdges.insert( Edge(currVtx, neighbors[0]) );
			}
			else if( numVtxs > 1) {
				int randIndex = rand() % numVtxs;	
				isMatched[ currVtx ] = true;
                                isMatched [  neighbors[ randIndex] ] = true;         
                                matchEdges.insert( Edge(currVtx, neighbors[ randIndex]) ); 
			}	
                }
        }

        //Count the vertices in the coarsened graph 
        //Store a mapping from coarsend vertex to the vertices in the finer graph        
        int numConargenedVertices = 0;
        vector< vector<unsigned int> > superRowsToRows;


	/*	
        //Match edges
        cout<<"Matched edges\n";
        for (set<Edge>::iterator it = matchEdges.begin(); it != matchEdges.end(); ++it) {
                    unsigned int u = it->either();
                    unsigned int v = it->other(u);

                    cout<<"Edge("<<u<<","<<v<<") "<<"\n";
	}*/	

	for (set<Edge>::iterator it = matchEdges.begin(); it != matchEdges.end(); ++it) {
		    unsigned int u = it->either();
		    unsigned int v = it->other(u);
	
                    perm[u] = numConargenedVertices;
                    perm[v] = numConargenedVertices;

                    vector<unsigned int> supRowContains(2,0);
		    supRowContains[0] = u;
		    supRowContains[1] = v;
		    superRowsToRows.push_back(supRowContains);
                    numConargenedVertices++;
        }
	
	/*	
	cout<<"Unmatched vertices\n";	
        for(int i = 0; i < g.n; ++i) {
               if( !isMatched[i] )
			cout<<i<<"\n"; 
	} */   


        for(int i = 0; i < g.n; ++i) {
               if( !isMatched[i] ) {
			perm[i] = numConargenedVertices;
			vector<unsigned int> supRowContains(1,i);
	                superRowsToRows.push_back(supRowContains);
			numConargenedVertices++;
               }
        }

        coarsenedG.n = numConargenedVertices;

	//Calculate the vertex weights of the coarsened graph
	coarsenedG.vertexWeight = (unsigned int *)malloc( numConargenedVertices * sizeof(unsigned int));
	for(int i = 0; i < coarsenedG.n; i++) { 
		unsigned int currVtxWeight = 0;
		for(int j = 0; j < superRowsToRows[i].size(); j++) {
			 currVtxWeight = currVtxWeight + g.vertexWeight[ superRowsToRows[i][j] ];
		}
		coarsenedG.vertexWeight[ i ] = currVtxWeight;
	}

        coarsenedG.num_edges = (unsigned int *)malloc( (coarsenedG.n+1) * sizeof(unsigned int) );

	//vectors to compute adjacency of coarsened graph and edge weights
	vector< unsigned int > adjVector;
	vector< unsigned int > weightVector;
        
	coarsenedG.num_edges[0] = 0;
	
	for(int i = 0; i < coarsenedG.n; i++) {

		//Maps neighbor to weight
		map<unsigned int, unsigned int> adjNeighborWeightMap;

		for(int j = 0; j < superRowsToRows[i].size(); j++) {
                      unsigned int vtxInCurSuperRow = superRowsToRows[i][j];
		      for(int j = g.num_edges[vtxInCurSuperRow]; j < g.num_edges[vtxInCurSuperRow+1]; j ++) {
				unsigned int adjVtx = g.adj[j];
				unsigned int connectedSupRow = perm[ adjVtx ];
			
				map<unsigned int, unsigned int>::iterator it = adjNeighborWeightMap.find( connectedSupRow );
				
				if( it != adjNeighborWeightMap.end() ) {	
					int totalWeight = g.edgeWeight[j] +  it->second;
					it->second = totalWeight;
				}
				else {
					adjNeighborWeightMap.insert( pair<unsigned int, unsigned int>(connectedSupRow, g.edgeWeight[j]) );
				}
                      }
                }
		
		//Sort the neighbors according to weight in increasing order
		vector< pair<unsigned int, unsigned int> > adjWeightVector( adjNeighborWeightMap.begin(), adjNeighborWeightMap.end());
		sort(adjWeightVector.begin(), adjWeightVector.end(), sortBySecond<unsigned int, unsigned int>()); 

                //Put the neighbors in the graph_t struct
                coarsenedG.num_edges[ i + 1 ] =  coarsenedG.num_edges[i] + adjWeightVector.size();

                for(int j = 0; j < adjWeightVector.size(); j++) {
			adjVector.push_back( adjWeightVector[j].first );
			weightVector.push_back( adjWeightVector[j].second ); 
		}
	}

	coarsenedG.m = adjVector.size();
	coarsenedG.orgM = coarsenedG.m;

	assert( coarsenedG.m == coarsenedG.num_edges[coarsenedG.n] );

        coarsenedG.adj = (unsigned int *)malloc( coarsenedG.m  * sizeof(unsigned int) );
        coarsenedG.edgeWeight = (unsigned int *)malloc( coarsenedG.m * sizeof(unsigned int) );

	for(int j = 0; j < coarsenedG.m; j++) {
		coarsenedG.adj [j] = adjVector[j];
		coarsenedG.edgeWeight [j] = weightVector[j];
	}
}

/* headyEdgeMatching - Find a maximal matching of the vertcies. Visit a vertex u in random
 * order and match u with v such that w(u,v) - weight of edge (u,v) is the maximum among the
 * valid neighbors of u.
 * Input: g - the graph with edge and vertex weight
 * Output: perm - mapping of vertices of g to the vertices of coarsendG
 *         coarsenedG - graph of the coarsened vertices */

void heavyEdgeMatching(graph_t g, unsigned int *perm, graph_t & coarsenedG) {
	
	//Matched edges   
        set<Edge> matchEdges; 

	//Initially vertices are not matched 
        bool isMatched[g.n];

	srand ( unsigned ( std::time(0) ) );
	vector<int> myvector;

	//Initalize thse values
	for (int i = 0; i < g.n; i++) {
                isMatched[i] = false; 
		myvector.push_back(i);
        }
	
	//Randomly shuffles the vertices
	std::random_shuffle ( myvector.begin(), myvector.end() );
	
	//Visit a vertex at random and match it with the heaviest edge
	for (int i = 0; i < g.n; i++) {
                int currVtx = myvector[i];

                if( !isMatched[currVtx] ) {
			for(int j = (g.num_edges[currVtx+1] - 1); j >= (int)g.num_edges[currVtx]; j--) {
				
				if( (g.adj[j] != currVtx) && !isMatched[ g.adj[j] ] ) {
                                	isMatched[ currVtx ] = true;
                                	isMatched [ g.adj[j] ] = true;
                                	matchEdges.insert( Edge(currVtx, g.adj[j]) );
					break;
				}
			}
                }
        }
      
 
	//Count the vertices in the coarsened graph 
	//Store a mapping from coarsend vertex to the vertices in the finer graph
        int numConargenedVertices = 0;
        vector< vector<unsigned int> > superRowsToRows;
	
        /*//Match edges
        cout<<"Matched edges\n";
        for (set<Edge>::iterator it = matchEdges.begin(); it != matchEdges.end(); ++it) {
                    unsigned int u = it->either();
                    unsigned int v = it->other(u);

                    cout<<"Edge("<<u<<","<<v<<") "<<"\n";
	}*/	

	for (set<Edge>::iterator it = matchEdges.begin(); it != matchEdges.end(); ++it) {
		    unsigned int u = it->either();
		    unsigned int v = it->other(u);
	
                    perm[u] = numConargenedVertices;
                    perm[v] = numConargenedVertices;

                    vector<unsigned int> supRowContains(2,0);
		    supRowContains[0] = u;
		    supRowContains[1] = v;
		    superRowsToRows.push_back(supRowContains);
                    numConargenedVertices++;
        }
	
/*	
	cout<<"Unmatched vertices\n";	
        for(int i = 0; i < g.n; ++i) {
               if( !isMatched[i] )
			cout<<i<<"\n"; 
	}    
	*/

        for(int i = 0; i < g.n; ++i) {
               if( !isMatched[i] ) {
			perm[i] = numConargenedVertices;
			vector<unsigned int> supRowContains(1,i);
	                superRowsToRows.push_back(supRowContains);
			numConargenedVertices++;
               }
        }

        coarsenedG.n = numConargenedVertices;

	//Calculate the vertex weights of the coarsened graph
	coarsenedG.vertexWeight = (unsigned int *)malloc( numConargenedVertices * sizeof(unsigned int));
	for(int i = 0; i < coarsenedG.n; i++) { 
		unsigned int currVtxWeight = 0;
		for(int j = 0; j < superRowsToRows[i].size(); j++) {
			 currVtxWeight = currVtxWeight + g.vertexWeight[ superRowsToRows[i][j] ];
		}
		coarsenedG.vertexWeight[ i ] = currVtxWeight;
	}

        coarsenedG.num_edges = (unsigned int *)malloc( (coarsenedG.n+1) * sizeof(unsigned int) );

	//vectors to compute adjacency and edge weights of coarsened graph
	vector< unsigned int > adjVector;
	vector< unsigned int > weightVector;
        
	coarsenedG.num_edges[0] = 0;
	
	for(int i = 0; i < coarsenedG.n; i++) {

		//Maps neighbor to weight
		map<unsigned int, unsigned int> adjNeighborWeightMap;

		for(int j = 0; j < superRowsToRows[i].size(); j++) {
                      unsigned int vtxInCurSuperRow = superRowsToRows[i][j];
		      for(int j = g.num_edges[vtxInCurSuperRow]; j < g.num_edges[vtxInCurSuperRow+1]; j ++) {
				unsigned int adjVtx = g.adj[j];
				unsigned int connectedSupRow = perm[ adjVtx ];
			
				map<unsigned int, unsigned int>::iterator it = adjNeighborWeightMap.find( connectedSupRow );
				
				if( it != adjNeighborWeightMap.end() ) {	
					int totalWeight = g.edgeWeight[j] +  it->second;
					it->second = totalWeight;
				}
				else {
					adjNeighborWeightMap.insert( pair<unsigned int, unsigned int>(connectedSupRow, g.edgeWeight[j]) );
				}
                      }
                }
		
		//Sort the neighbors according to weights in increasing order
		vector< pair<unsigned int, unsigned int> > adjWeightVector( adjNeighborWeightMap.begin(), adjNeighborWeightMap.end());
		sort(adjWeightVector.begin(), adjWeightVector.end(), sortBySecond<unsigned int, unsigned int>()); 

                //Put the neighbors in the graph_t struct
                coarsenedG.num_edges[ i + 1 ] =  coarsenedG.num_edges[i] + adjWeightVector.size();

                for(int j = 0; j < adjWeightVector.size(); j++) {
			adjVector.push_back( adjWeightVector[j].first );
			weightVector.push_back( adjWeightVector[j].second ); 
		}
	}

	coarsenedG.m = adjVector.size();
	coarsenedG.orgM = coarsenedG.m;

	assert( coarsenedG.m == coarsenedG.num_edges[coarsenedG.n] );

        coarsenedG.adj = (unsigned int *)malloc( coarsenedG.m  * sizeof(unsigned int) );
        coarsenedG.edgeWeight = (unsigned int *)malloc( coarsenedG.m * sizeof(unsigned int) );

	for(int j = 0; j < coarsenedG.m; j++) {
		coarsenedG.adj [j] = adjVector[j];
		coarsenedG.edgeWeight [j] = weightVector[j];
	}
}

/* lightEdgeMatching - Finds a maximal matching of the vertcies. Visit a vertex u in random
 * order and match u with v such that w(u,v) - weight of edge (u,v) is the minimum among the
 * valid neighbors of u.
 * Input: g - the graph with edge and vertex weight
 * Output: perm - mapping of vertices of g to the vertices of coarsendG
 *         coarsenedG - graph of the coarsened vertices */

void lightEdgeMatching(graph_t g, unsigned int *perm, graph_t & coarsenedG) {
	
	//Matched edges   
        set<Edge> matchEdges; 

	//Initially vertices are not matched 
        bool isMatched[g.n];

	srand ( unsigned ( std::time(0) ) );
	vector<int> myvector;

	//Initialize these values
	for (int i = 0; i < g.n; i++) {
                isMatched[i] = false; 
		myvector.push_back(i);
        }
	
	//Randomly shuffles the vertices
	std::random_shuffle ( myvector.begin(), myvector.end() );
	
	//Visit a vertex at random and match it with the heaviest edge
	for (int i = 0; i < g.n; i++) {
                int currVtx = myvector[i];

                if( !isMatched[currVtx] ) {
			for(int j = g.num_edges[currVtx]; j < g.num_edges[currVtx+1]; j++) {
				
				if( (g.adj[j] != currVtx) && !isMatched[ g.adj[j] ] ) {
                                	isMatched[ currVtx ] = true;
                                	isMatched [ g.adj[j] ] = true;
                                	matchEdges.insert( Edge(currVtx, g.adj[j]) );
					break;
				}
			}
                }
        }
      
 
	//Count the vertices in the coarsened graph 
	//Store a mapping from coarsend vertex to the vertices in the finer graph
        int numConargenedVertices = 0;
        vector< vector<unsigned int> > superRowsToRows;
	
        /*//Match edges
        cout<<"Matched edges\n";
        for (set<Edge>::iterator it = matchEdges.begin(); it != matchEdges.end(); ++it) {
                    unsigned int u = it->either();
                    unsigned int v = it->other(u);

                    cout<<"Edge("<<u<<","<<v<<") "<<"\n";
	}*/	

	for (set<Edge>::iterator it = matchEdges.begin(); it != matchEdges.end(); ++it) {
		    unsigned int u = it->either();
		    unsigned int v = it->other(u);
	
                    perm[u] = numConargenedVertices;
                    perm[v] = numConargenedVertices;

                    vector<unsigned int> supRowContains(2,0);
		    supRowContains[0] = u;
		    supRowContains[1] = v;
		    superRowsToRows.push_back(supRowContains);
                    numConargenedVertices++;
        }
	
/*	
	cout<<"Unmatched vertices\n";	
        for(int i = 0; i < g.n; ++i) {
               if( !isMatched[i] )
			cout<<i<<"\n"; 
	}    
	*/

        for(int i = 0; i < g.n; ++i) {
               if( !isMatched[i] ) {
			perm[i] = numConargenedVertices;
			vector<unsigned int> supRowContains(1,i);
	                superRowsToRows.push_back(supRowContains);
			numConargenedVertices++;
               }
        }

        coarsenedG.n = numConargenedVertices;

	//Calculate the vertex weights of the coarsened graph
	coarsenedG.vertexWeight = (unsigned int *)malloc( numConargenedVertices * sizeof(unsigned int));
	for(int i = 0; i < coarsenedG.n; i++) { 
		unsigned int currVtxWeight = 0;
		for(int j = 0; j < superRowsToRows[i].size(); j++) {
			 currVtxWeight = currVtxWeight + g.vertexWeight[ superRowsToRows[i][j] ];
		}
		coarsenedG.vertexWeight[ i ] = currVtxWeight;
	}

        coarsenedG.num_edges = (unsigned int *)malloc( (coarsenedG.n+1) * sizeof(unsigned int) );

	//vectors to compute adjacency and edge weights of coarsened graph
	vector< unsigned int > adjVector;
	vector< unsigned int > weightVector;
        
	coarsenedG.num_edges[0] = 0;
	
	for(int i = 0; i < coarsenedG.n; i++) {

		//Maps neighbor to weight
		map<unsigned int, unsigned int> adjNeighborWeightMap;

		for(int j = 0; j < superRowsToRows[i].size(); j++) {
                      unsigned int vtxInCurSuperRow = superRowsToRows[i][j];
		      for(int j = g.num_edges[vtxInCurSuperRow]; j < g.num_edges[vtxInCurSuperRow+1]; j ++) {
				unsigned int adjVtx = g.adj[j];
				unsigned int connectedSupRow = perm[ adjVtx ];
			
				map<unsigned int, unsigned int>::iterator it = adjNeighborWeightMap.find( connectedSupRow );
				
				if( it != adjNeighborWeightMap.end() ) {	
					int totalWeight = g.edgeWeight[j] +  it->second;
					it->second = totalWeight;
				}
				else {
					adjNeighborWeightMap.insert( pair<unsigned int, unsigned int>(connectedSupRow, g.edgeWeight[j]) );
				}
                      }
                }
		
		//Sort the neighbors according to weights in increasing order
		vector< pair<unsigned int, unsigned int> > adjWeightVector( adjNeighborWeightMap.begin(), adjNeighborWeightMap.end());
		sort(adjWeightVector.begin(), adjWeightVector.end(), sortBySecond<unsigned int, unsigned int>()); 

                //Put the neighbors in the graph_t struct
                coarsenedG.num_edges[ i + 1 ] =  coarsenedG.num_edges[i] + adjWeightVector.size();

                for(int j = 0; j < adjWeightVector.size(); j++) {
			adjVector.push_back( adjWeightVector[j].first );
			weightVector.push_back( adjWeightVector[j].second ); 
		}
	}

	coarsenedG.m = adjVector.size();
	coarsenedG.orgM = coarsenedG.m;

	assert( coarsenedG.m == coarsenedG.num_edges[coarsenedG.n] );

        coarsenedG.adj = (unsigned int *)malloc( coarsenedG.m  * sizeof(unsigned int) );
        coarsenedG.edgeWeight = (unsigned int *)malloc( coarsenedG.m * sizeof(unsigned int) );

	for(int j = 0; j < coarsenedG.m; j++) {
		coarsenedG.adj [j] = adjVector[j];
		coarsenedG.edgeWeight [j] = weightVector[j];
	}
}

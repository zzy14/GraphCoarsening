/* Written by : Humayun Kabir
 *          humayun.k1@gmail.com
 */

#include<iostream>
using namespace std;

#include "graph.h"


//Print a graph
void print_graph(graph_t & g) {

	cout<<"N: "<< g.n <<" M: "<<g.m;
	
	cout<<"\nVertex Weight:";
	for(int i = 0; i < g.n; i++) {
		cout<<"\nVertex: "<<i <<" has weight:" << g.vertexWeight[i];
	}
	cout<<"\n";
	
	for(int i = 0; i < g.n; i++) {
		cout<<"\nVertex: "<<i <<" has neighbors:";
		for(int j = g.num_edges[i]; j < g.num_edges[i+1]; j++) {
			cout<<"\nNeighbor: "<<g.adj[j]<<" Weight: "<<g.edgeWeight[j];
		}
	cout<<"\n";
	}
}


void load_graph_from_file(char *inFileName, graph_t *g) {

        FILE *infp;

        infp = fopen(inFileName, "rb");
        if (infp == NULL) {
                fprintf(stderr, "Error: could not open file to read graph: %s.\n", inFileName);
                exit(1);
        }

        //Read number of vertices and edges
        fread(&(g->n), sizeof(long), 1, infp);
        fread(&(g->m), sizeof(long), 1, infp);
        fread(&(g->orgM), sizeof(long), 1, infp);
        long notNeededVal;
        fread(&notNeededVal, sizeof(long), 1, infp);
        fread(&notNeededVal, sizeof(long), 1, infp);
        fread(&notNeededVal, sizeof(long), 1, infp);
        fread(&notNeededVal, sizeof(long), 1, infp);

        printf("N: %ld, M: %ld, orgM: %ld \n", g->n, g->m, g->orgM);

       long m = g->m;

        //Allocate space
        g->num_edges = (unsigned int *) malloc((g->n + 1) * sizeof(unsigned int));
        assert(g->num_edges != NULL);
        g->adj = (unsigned int *) malloc(m * sizeof(unsigned int));
        assert(g->adj != NULL);
        double * val = (double *) malloc(m * sizeof(double));
        assert(val != NULL);

        fread(g->num_edges, sizeof(unsigned int), (g->n) + 1, infp);
        fread(g->adj, sizeof(unsigned int), m, infp);
        fread(val, sizeof(double), m, infp);

        fclose(infp);

        //Free memory
        free(val);

        g->edgeWeight = (unsigned int *) malloc(m * sizeof(unsigned int));
        assert(g->edgeWeight != NULL);

        for(int i = 0; i < m; i++)
                g->edgeWeight[i] = 1;

        g->vertexWeight = (unsigned int *) malloc(g->n * sizeof(unsigned int));
        assert(g->vertexWeight != NULL);

        for(int i = 0; i < g->n; i++)
                g->vertexWeight[i] = 1;

        fprintf(stdout, " Finished reading the original graph from disk. \n");
}        
        




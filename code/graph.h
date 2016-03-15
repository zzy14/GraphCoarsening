/* Written by : Humayun Kabir
 * 	    humayun.k1@gmail.com
 */


#ifndef _GRAPH_H
#define _GRAPH_H



//Structure to represent a graph
typedef struct {
        long n;
        long m;
        long orgM;

        unsigned int *num_edges;
        unsigned int *adj;
        unsigned int *edgeWeight;
        unsigned int *vertexWeight;
} graph_t;


//Structure to represent an edge
struct Edge {
   unsigned int u;
   unsigned int v;

   Edge(unsigned int inU, unsigned int inV) : u(inU), v(inV) {}

   unsigned int either() const { return u; }
   unsigned int other(unsigned int inU) const {
        if (inU == u )
                return v;
        else
                return u;
   }

   bool operator < (const Edge right) const {
        if( u < right.u )
                return true;
        else if( v < right.v)
                return true;
        return false;
   }
};

//Print the graph
void print_graph(graph_t & g);

//Load a graph from file
void load_graph_from_file(char *filename, graph_t *g);

#endif

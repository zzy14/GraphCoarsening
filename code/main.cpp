/* Written by : Humayun Kabir
 *          humayun.k1@gmail.com
 */


#include "graph_coarsening.cpp"

int main(int argc, char *argv[]) {

	//Declare a graph
        graph_t g;

	//The program needs three parameters
        if(argc < 3) {
                cout<<"Usage: "<<argv[0]<<" graph-name "<< "super-row-size"<<endl;
                exit(1);
        }

	//locad the graph
        load_graph_from_file(argv[1], &g);
	
	// Declare the coarsened graph
        graph_t coarsenedG;

        //Maximum number of rows in a super-row
        int superRowSize = atoi( argv[2] );
	
	//Use these two arrays to find mapping from coarsest to finest vertices
	unsigned int *numEdgesSupRowsToRows;
	unsigned int *mapSupRowstoRows = (unsigned int *)malloc(g.n * sizeof(unsigned int));

	//coarsen g and get coarsenedG where each vertex contains superRowSize of vertices of g
	coarsenGraph(g, superRowSize, numEdgesSupRowsToRows, mapSupRowstoRows, coarsenedG);

	// cout<<"number of coarsest vertices: "<<coarsenedG.n<<endl;

	// for(int i = 0; i < coarsenedG.n; i++) {
	// 	cout<<"Coarsest vertex: "<< i <<" contains vertices: \n"; 
	// 	for(int j = numEdgesSupRowsToRows[i]; j < numEdgesSupRowsToRows[i+1]; j++)
	// 		cout<<mapSupRowstoRows[j]<<"\t";
	// 	cout<<endl;
	// }

    // write look up
    FILE* fout;
    fout = fopen("lookup.txt", "wb");
    if (fout == NULL) {
        fprintf(stderr, "Error: could not open file to write lookup.\n");
        exit(1);
    }

    for(int i = 0; i < coarsenedG.n; i++) {
        for(int j = numEdgesSupRowsToRows[i]; j < numEdgesSupRowsToRows[i+1]; j++)
        {
            fprintf(fout, "%u", mapSupRowstoRows[j]);
            for(int k = coarsenedG.num_edges[i]; k < coarsenedG.num_edges[i+1]; k++) {
                fprintf(fout, " %u", coarsenedG.adj[k]);
            }
            fprintf(fout, "\n");
        }
    }
    fclose(fout);

	//Free memory
	free( numEdgesSupRowsToRows );
	free( mapSupRowstoRows );

        return 0;
}



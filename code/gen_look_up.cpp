#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <queue>
#include <algorithm>  
#include <vector>  
#include <functional> 
#include "graph_coarsening.cpp"
using namespace std;

void load_coarse(char *filename, unsigned int * coarse, long n, unsigned int* coarse_num)
{
    FILE *infp;

    infp = fopen(filename, "rb");
    if (infp == NULL) {
            fprintf(stderr, "Error: could not open file to read coarse: %s.\n", filename);
            exit(1);
    }

    unsigned int i, j;
    fscanf(infp, "%u %u", &i, coarse_num);

    for (long k = 0; k < n; k++)
    {
        fscanf(infp, "%u %u", &i, &j);
        coarse[i] = j;
    }

    fclose( infp );
}

int main(int argc, char *argv[]) {
    int support_num = 5;
    int max_depth = 5;
    int walk_num = 40;
    unsigned int coarse_num;

    //Declare a graph
    graph_t g;

    //The program needs three parameters
    if(argc < 3) {
            cout<<"Usage: "<<argv[0]<<" graph-name "<< " coarse-name "<<endl;
            exit(1);
    }

    //load the graph
    load_graph_from_file(argv[1], &g);
    
    //load coarse dict
    unsigned int* coarse = (unsigned int *) malloc(g.n * sizeof(unsigned int));
    load_coarse(argv[2], coarse, g.n, &coarse_num);

    unsigned int* look_up = (unsigned int *) malloc(g.n * sizeof(unsigned int) * support_num);
    unsigned int* addr = look_up;

    int* table = (int *) malloc(coarse_num * sizeof(int));

    unsigned int cur_n; //当前节点
    unsigned int neigh_num;
    for (int i = 0; i < g.n; i++)
    {
        memset(table, 0, coarse_num * sizeof(unsigned int));
        for (int j = 0; j < walk_num; j++)
        {
            cur_n = i;
            for (int k = 0; k < max_depth; k++)
            {
                neigh_num = g.num_edges[cur_n+1] - g.num_edges[cur_n];
                if (neigh_num == 0)
                    break;
                cur_n = g.adj[g.num_edges[cur_n] + rand() % neigh_num];
                table[coarse[cur_n]] += 1;
            }
        }
        priority_queue <int, vector<int>, greater<int> > pq; // 小的在首
        for (int j = 0; j < support_num; j++)
        {
            pq.push(table[j]);
        }
        for (int j = support_num; j < coarse_num; j++)
        {
            pq.push(table[j]);
            pq.pop();
        }
        int min_count = pq.top();
        int pos = 0;
        for (int j = 0; j < coarse_num; j++)
        {
            if (table[j] >= min_count)
            {
                addr[pos] = j;
                pos++;
            }
            if (pos == support_num)
                break;
        }
        if (i % 1000 == 0)
        {
            printf("Progress: %.3f%%", (double)i / (double)(g.n + 1) * 100);
            fflush(stdout);
        }
        addr += support_num;
    }
    
    // write look up
    FILE* fout;
    fout = fopen("look_up_raw.txt", "wb");
    if (fout == NULL) {
        fprintf(stderr, "Error: could not open file to write lookup.\n");
        exit(1);
    }
    fprintf(fout, "%ld %u\n", g.n, coarse_num);
    addr = look_up;
    for(int i = 0; i < g.n; i++) {
        fprintf(fout, "%d", i);
        for (int j = 0; j < support_num; j++)
        {
            fprintf(fout, " %u", addr[j]);
        }
        fprintf(fout, "\n");
        addr += support_num;
    }
    fclose(fout);

    //Free memory
    free( coarse );
    free( look_up );

    return 0;
}



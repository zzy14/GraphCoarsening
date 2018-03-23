#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <queue>
#include <algorithm>  
#include <vector>  
#include <functional> 
#include "graph_coarsening.cpp"
using namespace std;

int MAX_DEPTH = 5;
int WALK_NUM = 40;
int SUPPORT_NUM = 5;
int total_count = 0;
unsigned int coarse_num;
graph_t g;
unsigned int* coarse;
unsigned int* look_up;

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

void* look_up_thread(void* params)
{
    int begin = ((int*)params)[0];
    int end = ((int*)params)[1];
    int count = 0;
    // printf("%d %d\n", begin, end);

    int* table = (int *) malloc(coarse_num * sizeof(int));
    unsigned int cur_n; //当前节点
    unsigned int neigh_num;

    unsigned int* addr = look_up + begin*SUPPORT_NUM;

    for (int i = begin; i < end; i++)
    {
        memset(table, 0, coarse_num * sizeof(unsigned int));
        for (int j = 0; j < WALK_NUM; j++)
        {
            cur_n = i;
            for (int k = 0; k < MAX_DEPTH; k++)
            {
                neigh_num = g.num_edges[cur_n+1] - g.num_edges[cur_n];
                if (neigh_num == 0)
                    break;
                cur_n = g.adj[g.num_edges[cur_n] + rand() % neigh_num];
                table[coarse[cur_n]] += 1;
            }
        }
        priority_queue <int, vector<int>, greater<int> > pq; // 小的在首
        for (int j = 0; j < SUPPORT_NUM; j++)
        {
            pq.push(table[j]);
        }
        for (int j = SUPPORT_NUM; j < coarse_num; j++)
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
            if (pos == SUPPORT_NUM)
                break;
        }
        addr += SUPPORT_NUM;
        count ++;
        if (count % 1000 == 0)
        {
            total_count += 1000;
            printf("%cProgress: %.3lf%%", 13, (double)total_count / (double)g.n);
            fflush(stdout);
        }
    }

    free(table);
    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    int num_threads = 1;

    //The program needs three parameters
    if(argc < 4) {
            cout<<"Usage: "<<argv[0]<<" graph-name "<< " coarse-name " << "num-threads" <<endl;
            exit(1);
    }

    num_threads = atoi(argv[3]);
    //load the graph
    load_graph_from_file(argv[1], &g);
    
    //load coarse dict
    coarse = (unsigned int *) malloc(g.n * sizeof(unsigned int));
    load_coarse(argv[2], coarse, g.n, &coarse_num);

    look_up = (unsigned int *) malloc(g.n * sizeof(unsigned int) * SUPPORT_NUM);
    unsigned int* addr = look_up;

    pthread_t *pt = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    int* params = (int*)malloc(2*num_threads*sizeof(int));

    int base_size = (int)g.n / num_threads + 1;
    for (int a = 0; a < num_threads; a++)
    {
        params[2*a] = base_size * a;
        if (base_size*(a+1) > g.n)
            params[2*a+1] = g.n;
        else
            params[2*a+1] = base_size*(a+1);
    }

    for (int a = 0; a < num_threads; a++) 
    {
        pthread_create(&pt[a], NULL, look_up_thread, (void *)(params + 2*a));
    }
    for (int a = 0; a < num_threads; a++) pthread_join(pt[a], NULL);
    
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
        for (int j = 0; j < SUPPORT_NUM; j++)
        {
            fprintf(fout, " %u", addr[j]);
        }
        fprintf(fout, "\n");
        addr += SUPPORT_NUM;
    }
    fclose(fout);

    //Free memory
    free( coarse );
    free( look_up );
    free ( params );

    return 0;
}



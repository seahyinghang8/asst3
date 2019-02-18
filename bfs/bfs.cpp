#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>
#include <iostream>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
}

void vertex_set_clear(vertex_set_complex* list) {
    list->count = 0;
    list->vertices.clear();
}

void vertex_set_init(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

void vertex_set_init(vertex_set_complex* list, int count) {
    list->max_vertices = count;
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{

    for (int i=0; i<frontier->count; i++) {

        int node = frontier->vertices[i];

        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->outgoing_starts[node + 1];

        // attempt to add all neighbors to the new frontier
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor];

            if (distances[outgoing] == NOT_VISITED_MARKER) {
                distances[outgoing] = distances[node] + 1;
                int index = new_frontier->count++;
                new_frontier->vertices[index] = outgoing;
            }
        }
    }
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step_parallel(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{

    #pragma omp parallel for schedule(dynamic, 32)
    for (int i=0; i<frontier->count; i++) {

        int node = frontier->vertices[i];

        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->outgoing_starts[node + 1];

        // create an array for all the neighbours that might
        // potentially be added to the new frontier
        int num_neighbours = end_edge - start_edge;
        int to_be_added[num_neighbours];
        int counter = 0;

        // attempt to add all neighbors to the new frontier
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor];
            
            if (distances[outgoing] == NOT_VISITED_MARKER) {
                if (__sync_bool_compare_and_swap(&distances[outgoing], NOT_VISITED_MARKER, distances[node] + 1)) {
                    to_be_added[counter] = outgoing;
                    counter++;
                }
            }
        }
        // add all of the newly found neighbour to the new frontier in one shot
        if (counter > 0) {
            int index = __sync_fetch_and_add(&new_frontier->count, counter);
            for (int j=0; j<counter; j++) {
                new_frontier->vertices[index+j] = to_be_added[j];
            }
        }
    }
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution* sol) {

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        if (frontier->count >= 32 * 4){
            top_down_step_parallel(graph, frontier, new_frontier, sol->distances);
        } else {
            top_down_step(graph, frontier, new_frontier, sol->distances);
        }

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}

// Take one step of "bottom-up" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void bottom_up_step(
    Graph g,
    vertex_set_complex* frontier,
    vertex_set_complex* new_frontier,
    int* distances)
{
    for (int i=0; i<g->num_nodes; i++) {
        if (distances[i] == NOT_VISITED_MARKER) {
            int start_edge = g->incoming_starts[i];
            int end_edge = (i == g->num_nodes - 1)
                            ? g->num_edges
                            : g->incoming_starts[i + 1];
            for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
                int node = g->incoming_edges[neighbor];
                if (frontier->vertices.count(node)){
                    int index = new_frontier->count++;
                    new_frontier->vertices.insert(i);
                    distances[i] = distances[node] + 1;
                    break;
                }
            }
        }

    }
}

// Take one step of "bottom-up" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void bottom_up_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    for (int i=0; i<g->num_nodes; i++) {
        if (distances[i] == NOT_VISITED_MARKER) {
            int start_edge = g->incoming_starts[i];
            int end_edge = (i == g->num_nodes - 1)
                            ? g->num_edges
                            : g->incoming_starts[i + 1];
            for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
                int node = g->incoming_edges[neighbor];
                if (frontier->vertices[node]){
                    int index = new_frontier->count++;
                    new_frontier->vertices[i] = 1;
                    distances[i] = distances[node] + 1;
                    break;
                }
            }
        }

    }
}

// Take one step of "bottom-up" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void bottom_up_step_parallel(
    Graph g,
    vertex_set_complex* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    //int size = g->num_nodes / omp_get_num_threads();
    int size = 1;
    #pragma omp parallel for schedule(auto)
    for (int j=0; j<g->num_nodes; j+=size) {
        int counter = 0;
        int temp[size];
        int d[size];
        for (int i=j; i <g->num_nodes && i<j+size; i++) {
            if (distances[i] == NOT_VISITED_MARKER) {
                int start_edge = g->incoming_starts[i];
                int end_edge = (i == g->num_nodes - 1)
                    ? g->num_edges
                    : g->incoming_starts[i + 1];
                for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
                    int node = g->incoming_edges[neighbor];
                    if (frontier->vertices.count(node)){
                        /*temp[counter] = i;
                        d[counter] = distances[node] + 1;
                        counter ++; */
                        int index = __sync_fetch_and_add(&new_frontier->count, 1);
                        new_frontier->vertices[index] = i;
                        distances[i] = distances[node] + 1;
                        break;
                    }
                }
            }
        }
        /*int index = __sync_fetch_and_add(&new_frontier->count, counter);
        for (int i=0; i<counter; i++){
            new_frontier->vertices[index+i] = temp[i];
            distances[temp[i]] = d[i];
        }*/

    }
}

// Take one step of "bottom-up" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void bottom_up_step_parallel(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    int size = g->num_nodes / 32 / 2; 
    #pragma omp parallel for schedule(auto)
    for (int j=0; j<g->num_nodes; j+=size) {
        int counter = 0;
        int temp[size];
        int d[size];
        for (int i=j; i <g->num_nodes && i<j+size; i++) {
            if (distances[i] == NOT_VISITED_MARKER) {
                int start_edge = g->incoming_starts[i];
                int end_edge = (i == g->num_nodes - 1)
                    ? g->num_edges
                    : g->incoming_starts[i + 1];
                for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
                    int node = g->incoming_edges[neighbor];
                    if (frontier->vertices[node]){
                        temp[counter] = i;
                        d[counter] = distances[node] + 1;
                        counter ++;
                        /*int index = __sync_fetch_and_add(&new_frontier->count, 1);
                        new_frontier->vertices[i] = 1;
                        distances[i] = distances[node] + 1;*/
                        break;
                    }
                }
            }
        }
        int index = __sync_fetch_and_add(&new_frontier->count, counter);
        for (int i=0; i<counter; i++){
            new_frontier->vertices[temp[i]] = 1;
            distances[temp[i]] = d[i];
        }

    }
}
/*
// Take one step of "bottom-up" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void bottom_up_step_parallel_vec(
    Graph g,
    std::vector<int>* to_check,
    vertex_set_complex* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    #pragma omp parallel for schedule(auto)
    for (int j=0; j<to_check->size(); j++) {
        int i = (*to_check)[j];
            int start_edge = g->incoming_starts[i];
            int end_edge = (i == g->num_nodes - 1)
                            ? g->num_edges
                            : g->incoming_starts[i + 1];
            for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
                int node = g->incoming_edges[neighbor];
                if (frontier->vertices.count(node)){
                    int index = __sync_fetch_and_add(&new_frontier->count, 1);
                    new_frontier->vertices[index] = i;
                    distances[i] = distances[node] + 1;
                    break;
                }
            }

    }
}
*/

void bfs_bottom_up(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier= &list1;
    vertex_set* new_frontier= &list2;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier with the root node
    frontier->count++;
    memset(frontier->vertices, 0, graph->num_nodes*sizeof(int));
    frontier->vertices[ROOT_NODE_ID] = 1;
    sol->distances[ROOT_NODE_ID] = 0;



    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);
        memset(new_frontier->vertices, 0, graph->num_nodes*sizeof(int));

        bottom_up_step_parallel(graph, frontier, new_frontier, sol->distances);

        // int count = 0;
        // for (int i=0; i<graph->num_nodes; i++)
        //     if (sol->distances[i] != NOT_VISITED_MARKER) {
        //         count += 1;
        //     }
        // std::cout<< count << std::endl;

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier= new_frontier;
        new_frontier= tmp;
    }
}

void bfs_bottom_up_complex(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.

    vertex_set_complex list1;
    vertex_set_complex list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set_complex* frontier_complex = &list1;
    vertex_set_complex* new_frontier_complex = &list2;

    vertex_set list4;
    vertex_set_init(&list4, graph->num_nodes);

    vertex_set* new_frontier = &list4;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier with the root node
    frontier_complex->count++;
    frontier_complex->vertices.insert(ROOT_NODE_ID);
    sol->distances[ROOT_NODE_ID] = 0;



    while (frontier_complex->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        bottom_up_step_parallel(graph, frontier_complex, new_frontier, sol->distances);

        vertex_set_clear(new_frontier_complex);
        new_frontier_complex->count = new_frontier->count;
        for (int i=0; i < new_frontier->count; i++) {
            new_frontier_complex->vertices.insert(new_frontier->vertices[i]);
        }
        // int count = 0;
        // for (int i=0; i<graph->num_nodes; i++)
        //     if (sol->distances[i] != NOT_VISITED_MARKER) {
        //         count += 1;
        //     }
        // std::cout<< count << std::endl;

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set_complex* tmp = frontier_complex;
        frontier_complex = new_frontier_complex;
        new_frontier_complex = tmp;
    }
}

void bfs_hybrid(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
    
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;
    
    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    int* tmp_vertices = (int*)malloc(sizeof(int) * graph->num_nodes);

    bool bottom_up = false;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);
        if (bottom_up) {
            memset(new_frontier->vertices, 0, graph->num_nodes*sizeof(int));
            bottom_up_step_parallel(graph, frontier, new_frontier, sol->distances);
            if (new_frontier->count > graph->num_nodes/2) {
            } else {
                bottom_up = false;
                int j = 0;
                for (int i = 0; i < graph->num_nodes; i++) {
                    if (new_frontier->vertices[i] == 1){
                        tmp_vertices[j] = i;
                        j++;
                    }
                }
                for (int i = 0; i < new_frontier->count; i++) {
                    new_frontier->vertices[i] = tmp_vertices[i];
                }
            }
        }
        else {
            top_down_step_parallel(graph, frontier, new_frontier, sol->distances);
            if (new_frontier->count > graph->num_nodes/2) {
                bottom_up = true;
                for (int i = 0; i < new_frontier->count; i++) {
                    tmp_vertices[i] = new_frontier->vertices[i];
                }
                memset(new_frontier->vertices, 0, graph->num_nodes*sizeof(int));
                for (int i = 0; i < new_frontier->count; i++) {
                    new_frontier->vertices[tmp_vertices[i]] = 1;
                }
            } else {
            }
        }

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}

void bfs_hybrid_complex(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
    
    vertex_set_complex list1;
    vertex_set_complex list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set_complex* frontier_complex = &list1;
    vertex_set_complex* new_frontier_complex = &list2;
    
    vertex_set list3;
    vertex_set list4;
    vertex_set_init(&list3, graph->num_nodes);
    vertex_set_init(&list4, graph->num_nodes);

    vertex_set* frontier = &list3;
    vertex_set* new_frontier = &list4;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    bool complex_set = false;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);
        if (complex_set) {
            bottom_up_step_parallel(graph, frontier_complex, new_frontier, sol->distances);
        }
        else {
            top_down_step_parallel(graph, frontier, new_frontier, sol->distances);
        }

        if (new_frontier->count > graph->num_nodes/2) {
            complex_set = true;

            vertex_set_clear(frontier_complex);
            frontier_complex->count = new_frontier->count;
            for (int i=0; i < new_frontier->count; i++) {
                frontier_complex->vertices.insert(new_frontier->vertices[i]);
            }
        } else {
            complex_set = false;
        }

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        if (!complex_set) {
            vertex_set* tmp = frontier;
            frontier = new_frontier;
            new_frontier = tmp;
        }
    }
}

void bfs_hybrid_old(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
    
    vertex_set_complex list1;
    vertex_set_complex list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set_complex* frontier_complex = &list1;
    vertex_set_complex* new_frontier_complex = &list2;
    
    vertex_set list3;
    vertex_set list4;
    vertex_set_init(&list3, graph->num_nodes);
    vertex_set_init(&list4, graph->num_nodes);

    vertex_set* frontier = &list3;
    vertex_set* new_frontier = &list4;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    bool complex_set = false;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        if (frontier->count > graph->num_nodes/2) {
            vertex_set_clear(new_frontier_complex);
            if (complex_set) {
                bottom_up_step(graph, frontier_complex, new_frontier_complex, sol->distances);
                frontier->count = new_frontier_complex->count;
            } else {
                complex_set = true;
                vertex_set_clear(frontier_complex);
                frontier_complex->count = frontier->count;
                for (int i=0; i < frontier->count; i++) {
                    frontier_complex->vertices.insert(frontier->vertices[i]);
                }

                vertex_set_clear(new_frontier_complex);

                bottom_up_step(graph, frontier_complex, new_frontier_complex, sol->distances);
                frontier->count = new_frontier_complex->count;
            }
        } else {
            vertex_set_clear(new_frontier);
            if (complex_set) {
                complex_set = false;
                vertex_set_clear(frontier);
                vertex_set_clear(new_frontier);

                frontier->count = frontier_complex->count;
                int counter = 0;
                for (int node: frontier_complex->vertices) {
                    frontier->vertices[counter] = node;
                    counter++;
                }

                top_down_step_parallel(graph, frontier, new_frontier, sol->distances);
            } else {
                top_down_step_parallel(graph, frontier, new_frontier, sol->distances);
            }
        }

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        if (complex_set) {
            vertex_set_complex* tmp = frontier_complex;
            frontier_complex = new_frontier_complex;
            new_frontier_complex = tmp;
        } else {
            vertex_set* tmp = frontier;
            frontier = new_frontier;
            new_frontier = tmp;
        }
    }
}

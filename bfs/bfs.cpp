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
//#define VERBOSE true

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
}

void vertex_set_init(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}


// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighbouring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{   
    // obtain the distance for the current frontier
    int new_frontier_distance = distances[frontier->vertices[0]] + 1;

    for (int i=0; i<frontier->count; i++) {

        int node = frontier->vertices[i];

        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->outgoing_starts[node + 1];

        // attempt to add all neighbours to the new frontier
        for (int neighbour=start_edge; neighbour<end_edge; neighbour++) {
            int outgoing = g->outgoing_edges[neighbour];

            if (distances[outgoing] == NOT_VISITED_MARKER) {
                distances[outgoing] = new_frontier_distance;
                int index = new_frontier->count++;
                new_frontier->vertices[index] = outgoing;
            }
        }
    }
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighbouring vertices to the
// new_frontier.
void top_down_step_parallel(
    Graph g,
    bool* visited,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    // obtain the distance for the current frontier
    int new_frontier_distance = distances[frontier->vertices[0]] + 1;

    #pragma omp parallel for schedule(static, 64)
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

        // attempt to add all neighbours to the new frontier
        for (int neighbour=start_edge; neighbour<end_edge; neighbour++) {
            int outgoing = g->outgoing_edges[neighbour];
            
            if (visited[outgoing]) continue;
            if (distances[outgoing] != NOT_VISITED_MARKER) continue;
            if (__sync_bool_compare_and_swap(&distances[outgoing], NOT_VISITED_MARKER, new_frontier_distance)) {
                to_be_added[counter] = outgoing;
                counter++;
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

void update_visited(vertex_set* frontier, bool* visited) {
    for (size_t i = 0; i < frontier->count; i++) {
        int index = frontier->vertices[i];
        visited[index] = true;
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

    // create a new array that shows if a node has been visited
    // which is only updated at the end of each iteration
    bool* visited = (bool*) malloc(graph->num_nodes * sizeof(bool));

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
        visited[i] = false;
    }

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        if (frontier->count >= 32 * 32){
            top_down_step_parallel(graph, visited, frontier, new_frontier, sol->distances);
        } else {
            top_down_step(graph, frontier, new_frontier, sol->distances);
        }

        update_visited(new_frontier, visited);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }

    free(visited);
}

// Take one step of "bottom-up" BFS.  For each vertex in the graph,
// check if the edge is not visited and its neighbour is a frontier, 
// add it to the new_frontier.
void bottom_up_step(
    Graph g,
    bool* visited,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{   
    // obtain the distance for the current frontier
    int new_frontier_distance = distances[frontier->vertices[0]] + 1;

    // loop through all the vertices in the graph
    for (int i = 0; i < g->num_nodes; i++) {
        if (visited[i]) continue;
        int start_edge = g->incoming_starts[i];
        int end_edge = (i == g->num_nodes - 1)
                        ? g->num_edges
                        : g->incoming_starts[i + 1];
        // loop through all the neighbours for the given vertex
        for (int neighbour = start_edge; neighbour < end_edge; neighbour++) {
            int incoming = g->incoming_edges[neighbour];
            // check if any neighbour has already been visited
            // if true, this means that the neighbour is on the frontier
            // and i should be added to the new_frontier
            if (!visited[incoming]) continue;
            int index = new_frontier->count++;
            new_frontier->vertices[index] = i;
            distances[i] = new_frontier_distance;
            break;
        }
    }
}

// Take one step of "bottom-up" BFS.  For each vertex in the graph,
// check if the edge is not visited and its neighbour is a frontier, 
// add it to the new_frontier.
void bottom_up_step_parallel(
    Graph g,
    bool* visited,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    // obtain the distance for the current frontier
    int new_frontier_distance = distances[frontier->vertices[0]] + 1;

    int nodes_per_thread = 256; 

    #pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < g->num_nodes; j += nodes_per_thread) {
        // for the last iteration, if the number of nodes assigned
        // exceeds the maximum, set a ceiling to the size
        int size = (nodes_per_thread + j >= g->num_nodes) ? g->num_nodes - j : nodes_per_thread;
        int counter = 0;
        int to_be_added[size];

        for (int i = j; i < j + size; i++) {
            if (visited[i]) continue;
            int start_edge = g->incoming_starts[i];
            int end_edge = (i == g->num_nodes - 1)
                ? g->num_edges
                : g->incoming_starts[i + 1];
            for (int neighbour=start_edge; neighbour<end_edge; neighbour++) {
                int incoming = g->incoming_edges[neighbour];
                // check if any neighbour has already been visited
                // if true, this means that the neighbour is on the frontier
                // and i should be added to the new_frontier
                if (!visited[incoming]) continue;
                to_be_added[counter] = i;
                counter++;
                break;
            }
        }

        int index = __sync_fetch_and_add(&new_frontier->count, counter);
        for (int i = 0; i < counter; i++){
            new_frontier->vertices[index + i] = to_be_added[i];
            distances[to_be_added[i]] = new_frontier_distance;
        }
    }
}

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

    // create a new array that shows if a node has been visited
    // which is only updated at the end of each iteration
    bool* visited = (bool*) malloc(graph->num_nodes * sizeof(bool));

    // initialize all nodes to NOT_VISITED
    for (int i = 0; i < graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
        visited[i] = false;
    }

    // setup frontier with the root node
    frontier->vertices[0] = ROOT_NODE_ID;
    frontier->count++;
    visited[ROOT_NODE_ID] = true;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif
        vertex_set_clear(new_frontier);

        //bottom_up_step(graph, visited, frontier, new_frontier, sol->distances);
        bottom_up_step_parallel(graph, visited, frontier, new_frontier, sol->distances);
        
        update_visited(new_frontier, visited);

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();

        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier= new_frontier;
        new_frontier= tmp;
    }

    free(visited);
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

    vertex_set* frontier= &list1;
    vertex_set* new_frontier= &list2;

    // create a new array that shows if a node has been visited
    // which is only updated at the end of each iteration
    bool* visited = (bool*) malloc(graph->num_nodes * sizeof(bool));

    // initialize all nodes to NOT_VISITED
    for (int i = 0; i < graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
        visited[i] = false;
    }

    // setup frontier with the root node
    frontier->vertices[0] = ROOT_NODE_ID;
    frontier->count++;
    visited[ROOT_NODE_ID] = true;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif
        vertex_set_clear(new_frontier);

        if ((double)frontier->count / graph->num_nodes > 0.03) {
            bottom_up_step_parallel(graph, visited, frontier, new_frontier, sol->distances);
        } else if (frontier->count >= 32 * 32){
            top_down_step_parallel(graph, visited, frontier, new_frontier, sol->distances);
        } else {
            top_down_step(graph, frontier, new_frontier, sol->distances);
        }
        
        update_visited(new_frontier, visited);

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();

        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier= new_frontier;
        new_frontier= tmp;
    }

    free(visited);
}

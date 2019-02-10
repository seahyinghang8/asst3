#include "page_rank.h"

#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <utility>

#include "../common/CycleTimer.h"
#include "../common/graph.h"


#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>

// pageRank --
//
// g:           graph to process (see common/graph.h)
// solution:    array of per-vertex vertex scores (length of array is numNodes(g))
// damping:     page-rank algorithm's damping parameter
// convergence: page-rank algorithm's convergence threshold
//
void pageRank(Graph g, double* solution, double damping, double convergence)
{


  // initialize vertex weights to uniform probability. Double
  // precision scores are used to avoid underflow for large graphs
  int numNodes = num_nodes(g);
  double equal_prob = 1.0 / numNodes;

  double* score_old = (double*) aligned_alloc(64, numNodes * sizeof(double));
  double* score_new = (double*) aligned_alloc(64, numNodes * sizeof(double));

  // track the sink nodes that will redistribute
  // its probability to other nodes
  Vertex* sink_nodes = (Vertex*) aligned_alloc(64, numNodes * sizeof(double));
  size_t num_sink_nodes = 0;

  for (int i = 0; i < numNodes; ++i) {
    score_old[i] = equal_prob;
    // add the sink node if it has no outgoing node
    if (outgoing_size(g, i) == 0) {
      sink_nodes[num_sink_nodes] = i;
      num_sink_nodes++;
    }
  }
  
  
  /*
     CS149 students: Implement the page rank algorithm here.  You
     are expected to parallelize the algorithm using openMP.  Your
     solution may need to allocate (and free) temporary arrays.

     Basic page rank pseudocode is provided below to get you started:

     // initialization: see example code above
     score_old[vi] = 1/numNodes;

     while (!converged) {

       // compute score_new[vi] for all nodes vi:
       score_new[vi] = sum over all nodes vj reachable from incoming edges
                          { score_old[vj] / number of edges leaving vj  }
       score_new[vi] = (damping * score_new[vi]) + (1.0-damping) / numNodes;

       score_new[vi] += sum over all nodes vj with no outgoing edges
                          { damping * score_old[vj] / numNodes }

       // compute how much per-node scores have changed
       // quit once algorithm has converged

       global_diff = sum over all nodes vi { abs(score_new[vi] - score_old[vi]) };
       converged = (global_diff < convergence)
     }

   */

  bool converged = false;
  while (!converged) {
    double partial_prob[omp_get_max_threads()];
    double partial_diff[omp_get_max_threads()];
    // set partial prob and partial diff to be 0
    for (size_t i = 0; i < omp_get_max_threads(); ++i) {
      partial_prob[i] = 0.;
      partial_diff[i] = 0.;
    }
    // calculate the partial increment in probability due to sink nodes
    #pragma omp parallel for schedule(auto)
    for (size_t i = 0; i < num_sink_nodes; ++i) {
      partial_prob[omp_get_thread_num()] += score_old[sink_nodes[i]] * damping / numNodes;
    }

    double sink_prob = 0.;
    for (size_t i = 0; i < omp_get_max_threads(); i++) {
      sink_prob += partial_prob[i];
    }

    // each chunk should be set to obtain the same cache line to prevent false sharing
    #pragma omp parallel for schedule(auto)
    for (size_t i = 0; i < numNodes; ++i) {
      // loop through all the incoming nodes
      score_new[i] = 0.;
      int num_incoming = incoming_size(g, i);
      const Vertex* v_arr = incoming_begin(g, i);
      for (size_t j = 0; j < num_incoming; ++j) {
        int incoming_vertex = v_arr[j];
        score_new[i] += score_old[incoming_vertex] / outgoing_size(g, incoming_vertex);
      }
      // account for damping
      score_new[i] = (damping * score_new[i]) + (1.0 - damping) / numNodes;
      score_new[i] += sink_prob;
      partial_diff[omp_get_thread_num()] += fabs(score_new[i] - score_old[i]);
    }


    double global_diff = 0.;
    for (size_t i = 0; i < omp_get_max_threads(); i++) {
      global_diff += partial_diff[i];
    }
    
    converged = (global_diff < convergence);

    // swap score_old and score_new
    double* temp = score_new;
    score_new = score_old;
    score_old = temp;
  }
  // update solution with score_new
  for (size_t i = 0; i < numNodes; ++i) {
    solution[i] = score_old[i];
  }

  free(score_old);
  free(score_new);
  free(sink_nodes);
}

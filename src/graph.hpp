/**
 *    @file  graph.hpp
 *   @brief  Interface function definitions for the sequence graph class.
 *
 *  This header file defines interface functions for the sequence graph class.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Fri Nov 11, 2016  01:08
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_GRAPH_HPP__
#define PSI_GRAPH_HPP__

#include <algorithm>
#include <stdexcept>
#include <functional>
#include <random>
#include <cstdint>

#include <gum/seqgraph.hpp>
#include <vg/io/stream.hpp>


namespace psi {
  namespace util {
    template< class TGraph >
    inline typename TGraph::offset_type
    max_node_len( TGraph const& graph )
    {
      typedef typename TGraph::id_type id_type;
      typedef typename TGraph::rank_type rank_type;
      typedef typename TGraph::offset_type offset_type;

      offset_type max = 1;
      graph.for_each_node(
          [&graph, &max]( rank_type rank , id_type id ) {
            auto len = graph.node_length( id );
            if ( max < len ) max = len;
            return true;
          } );
      return max;
    }

    template< class TGraph >
    inline typename TGraph::offset_type
    total_nof_loci( TGraph const& graph )
    {
      offset_type total = 0;
      graph.for_each_node(
          [&graph, &total]( rank_type rank, id_type id ) {
            total += graph.node_length( id );
            return true;
          } );
      return total;
    }

    template< class TGraph, typename TNodeIter, typename TEdgeIter >
    inline void
    induced_graph( TGraph const& graph,
                   TNodeIter nbegin, TNodeIter nend,
                   TEdgeIter ebegin, TEdgeIter eend,
                   vg::Graph* induced )
    {
      for ( ; nbegin != nend; ++nbegin ) {
        vg::Node* new_node = induced->add_node();
        new_node->set_id( *nbegin );
        new_node->set_sequence( graph.node_sequence( *nbegin ) );
      }

      for ( ; ebegin != eend; ++ebegin ) {
        vg::Edge* new_edge = induced->add_edge();
        new_edge->set_from( graph.from_id( *ebegin ) );
        new_edge->set_to( graph.to_id( *ebegin ) );
        new_edge->set_from_start( graph.is_from_start( *ebegin ) );
        new_edge->set_to_end( graph.is_to_end( *ebegin ) );
        new_edge->set_overlap( graph.edge_overlap( *ebegin ) );
      }
    }

    /**
     *  @brief  Get the induced graph of a set of nodes and edges in `vg::Graph` objects.
     *
     *  @param  nbegin The begin iterator of the nodes set.
     *  @param  nend The end iterator of the nodes set.
     *  @param  ebegin The begin iterator of the edges set.
     *  @param  eend The end iterator of the edges set.
     *  @param  callback The callback function passing each generated `vg::Graph` messages.
     *  @param  chunk_size The maximum number of nodes allowed in a `vg::Graph` message.
     *
     *  :TODO:Sat Apr 20 12:29:\@cartoonist: add ability of including path in the
     *  `vg::Graph`.
     */
    template< class TGraph, typename TNodeIter, typename TEdgeIter >
    inline void
    induced_graph( TGraph const& graph,
                   TNodeIter nbegin, TNodeIter nend,
                   TEdgeIter ebegin, TEdgeIter eend,
                   std::function< void( vg::Graph& ) > callback,
                   std::ptrdiff_t chunk_size ) const
    {
      assert( chunk_size < std::PTRDIFF_MAX );
      auto nodes_l = nbegin;
      auto nodes_r = nodes_l + 1;
      auto edges_l = ebegin;
      auto edges_r = edges_l + 1;
      auto nofmsg = ( nend - nbegin ) / chunk_size + 1;
      for ( unsigned int i = 0; i < nofmsg && nodes_r != nend; ++i ) {
        nodes_r = nodes_l + std::min( chunk_size, nend - nodes_l );
        auto maxid = *( nodes_r - 1 );
        while ( edges_r != eend && graph.from_id( *edges_r ) <= maxid ) ++edges_r;
        vg::Graph g;
        induced_graph( nodes_l, nodes_r, edges_l, edges_r, &g );
        callback( g );
        nodes_l = nodes_r;
        edges_l = edges_r;
      }
    }

    /**
     *  @brief  Get the node ID of an ajacent node randomly.
     *
     *  @param  vargraph The variation graph.
     *  @param  node_id The ID of the node whose an adjacent node should be returned.
     *  @return an adjacent node ID if available any; otherwise 0 -- an invalid node ID.
     *
     *  Picking one of the adjacent nodes by generating a pseudo-random number with
     *  uniform distribution in [0, out-degree(v)].
     */
    template< class TGraph >
    inline typename TGraph::id_type
    random_adjacent( TGraph const& graph, typename TGraph::id_type node_id,
                         unsigned int seed=0 )
    {
      auto odeg = graph.outdegree( node_id );
      if ( odeg == 0 ) return 0;

      std::random_device rd;  // Will be used to obtain a seed for the random no. engine
      if ( seed == 0 ) seed = rd(); // use random_device to generate a seed if seed is not provided
      std::mt19937 gen( seed );  // Standard mersenne_twister_engine seeded with seed
      std::uniform_int_distribution<> dis(0, odeg - 1);
      auto idx = dis(gen);

      typename TGraph::id_type candidate = 0;
      graph.for_each_edges_out(
          node_id,
          [&candidate, &idx]( id_type to, linktype_type type ) {
            if ( idx == 0 ) {
              candidate = to;
              return false;
            }
            --idx;
            return true;
          } );
      assert( candidate != 0 );
      return candidate;
    }

    /**
     *  @brief  Get the ID of an adjacent node with least coverage.
     *
     *  @param  vargraph The variation graph.
     *  @param  node_id The ID ot the node whose an adjacent node should be returned.
     *  @param  paths_set A set of paths as a container of `Path`.
     *  @return the ID of the one of adjacent nodes with least coverage. If there are
     *          multiple nodes with minimum coverage, it returns one of them. If all nodes
     *          are covered equally or no forward edge exists, it returns zero.
     *
     *  Calculate coverage for all adjacent nodes and find the smallest one.
     */
    template< class TGraph, typename TContainer >
    inline typename TGraph::id_type
    least_covered_adjacent( TGraph const& graph, typename TGraph::id_type node_id,
                            TContainer const& paths_set )
    {
      typedef TGraph graph_type;
      typedef typename graph_type::id_type id_type;

      id_type lc_id = 0;
      int lc_value = -1;
      bool equally_covered = true;
      auto cov =
          [&graph, &paths_set]( id_type id ) {
            return path_coverage( id, paths_set );
          };

      graph.for_each_edges_out(
          node_id,
          [&lc_id, &lc_value, &cov, &equally_covered]( id_type to, linktype_type type ) {
            auto value = cov( to );
            if ( lc_value != -1 && equally_covered && lc_value != value ) {
              equally_covered = false;
            }
            if ( lc_value == -1 || value < lc_value ) {
              lc_id = to;
              lc_value = value;
            }
            return true;
          } );

      return equally_covered ? 0 : lc_id;
    }

    template< class TGraph, typename TPath, typename TContainer >
    inline typename TGraph::id_type
    least_covered_adjacent( TGraph const& vargraph, TPath& tail,
                            TContainer const& paths_set )
    {
      typedef TGraph graph_type;
      typedef typename graph_type::id_type id_type;

      id_type lc_id = 0;
      int lc_value = -1;
      bool equally_covered = true;
      auto cov =
          [&graph, &paths_set]( TPath const& tail ) {
            return path_coverage( tail.begin(), tail.end(), paths_set );
          };

      if ( !tail.empty() ) {
        graph.for_each_edges_out(
            tail.back(),
            [&lc_id, &lc_value, &cov, &equally_covered, &tail]
            ( id_type to, linktype_type type ) {
              tail.push_back( to );  /* XXX: should be popped later! */
              auto value = cov( tail );
              tail.pop_back();
              if ( lc_value != -1 && equally_covered && lc_value != value ) {
                equally_covered = false;
              }
              if ( lc_value == -1 || value < lc_value ) {
                lc_id = to;
                lc_value = value;
              }
              return true;
            } );
      }

      return equally_covered ? 0 : lc_id;
    }
  }  /* --- end of namespace util --- */
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_GRAPH_HPP__ --- */

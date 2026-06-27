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
#include <functional>
#include <random>
#include <cstdint>
#include <vector>

#include <gum/graph.hpp>

#include "utils.hpp"


namespace psi {
  template< typename TId, typename TOffset >
  class PositionBase {
    public:
      /* === TYPE MEMBERS === */
      using id_type = TId;
      using offset_type = TOffset;
      /* === LIFECYCLE === */
      PositionBase( id_type id=0, offset_type offset=0 )
        : m_id( id ), m_offset( offset )
      { }

      PositionBase( const PositionBase& ) = default;
      PositionBase( PositionBase&& ) = default;
      ~PositionBase( ) = default;
      /* === OPERATORS === */
      PositionBase& operator=( const PositionBase& ) = default;
      PositionBase& operator=( PositionBase&& ) = default;
      /* === ACCESSORS === */
      inline id_type
      node_id( ) const
      {
        return this->m_id;
      }

      inline offset_type
      offset( ) const
      {
        return this->m_offset;
      }
      /* === MUTATORS === */
      inline void
      set_node_id( id_type e_id )
      {
        this->m_id = e_id;
      }

      inline void
      set_offset( offset_type e_offset )
      {
        this->m_offset = e_offset;
      }
    private:
      /* === DATA MEMBERS === */
      TId m_id;
      TOffset m_offset;
  };

  template< typename TId = typename gum::GraphBaseTrait< gum::Dynamic >::id_type,
            typename TOffset = typename gum::GraphBaseTrait< gum::Dynamic >::offset_type >
  using Position = psi::PositionBase< TId, TOffset >;

  namespace util {
    template< class TGraph, typename TNodeIter, typename TEdgeIter, typename TVGGraph,
              typename TCoordinate = gum::CoordinateType< TGraph, gum::coordinate::Identity,
                                                          decltype( TVGGraph().node(0).id() ) > >
    inline void
    induced_graph( TGraph const& graph,
                   TNodeIter nbegin, TNodeIter nend,
                   TEdgeIter ebegin, TEdgeIter eend,
                   TVGGraph* induced, TCoordinate&& coord={} )
    {
      for ( ; nbegin != nend; ++nbegin ) {
        auto new_node = induced->add_node();
        new_node->set_id( coord( *nbegin ) );
        std::string label = graph.node_sequence( *nbegin );
        new_node->set_sequence( label );
      }

      for ( ; ebegin != eend; ++ebegin ) {
        auto new_edge = induced->add_edge();
        new_edge->set_from( coord( graph.from_id( *ebegin ) ) );
        new_edge->set_to( coord( graph.to_id( *ebegin ) ) );
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
    template< class TGraph, typename TNodeIter, typename TEdgeIter, typename TVGGraph,
              typename TCoordinate = gum::CoordinateType< TGraph, gum::coordinate::Identity,
                                                          decltype( TVGGraph().node(0).id() ) > >
    inline void
    induced_graph( TGraph const& graph,
                   TNodeIter nbegin, TNodeIter nend,
                   TEdgeIter ebegin, TEdgeIter eend,
                   std::function< void( TVGGraph& ) > callback,
                   std::ptrdiff_t chunk_size, TCoordinate&& coord={} )
    {
      assert( chunk_size < PTRDIFF_MAX );
      auto nodes_l = nbegin;
      auto nodes_r = nodes_l + 1;
      auto edges_l = ebegin;
      auto edges_r = edges_l + 1;
      auto nofmsg = ( nend - nbegin ) / chunk_size + 1;
      for ( unsigned int i = 0; i < nofmsg && nodes_r != nend; ++i ) {
        nodes_r = nodes_l + std::min( chunk_size, nend - nodes_l );
        auto maxid = *( nodes_r - 1 );
        while ( edges_r != eend && graph.from_id( *edges_r ) <= maxid ) ++edges_r;
        TVGGraph g;
        induced_graph( graph, nodes_l, nodes_r, edges_l, edges_r, &g, coord );
        callback( g );
        nodes_l = nodes_r;
        edges_l = edges_r;
      }
    }

    /**
     *  @brief  Get the node ID of an ajacent node randomly.
     *
     *  @param  graph The graph.
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
      typedef TGraph graph_type;
      typedef typename graph_type::id_type id_type;
      typedef typename graph_type::rank_type rank_type;
      typedef typename graph_type::linktype_type linktype_type;

      thread_local static std::mt19937 lgen;
      thread_local static unsigned int lseed = std::mt19937::default_seed;

      auto odeg = graph.outdegree( node_id );
      if ( odeg == 0 ) return 0;

      rank_type idx;
      if ( seed == 0 ) idx = random::random_index( odeg );
      else {
        if ( seed != lseed ) {
          lseed = seed;
          lgen.seed( seed );
        }
        std::uniform_int_distribution< rank_type > dis( 0, odeg - 1 );
        idx = dis( lgen );
      }

      id_type candidate = 0;
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
     *  @param  graph The graph.
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
      typedef typename graph_type::linktype_type linktype_type;

      id_type lc_id = 0;
      std::size_t lc_value = UINTMAX_MAX;
      bool equally_covered = true;
      auto cov =
          [&graph, &paths_set]( id_type id ) {
            return path_coverage( id, paths_set );
          };

      graph.for_each_edges_out(
          node_id,
          [&lc_id, &lc_value, &cov, &equally_covered]( id_type to, linktype_type type ) {
            auto value = cov( to );
            if ( equally_covered && lc_value != UINTMAX_MAX && lc_value != value ) {
              equally_covered = false;
            }
            if ( value < lc_value ) {
              lc_id = to;
              lc_value = value;
            }
            return true;
          } );

      return equally_covered ? 0 : lc_id;
    }

    template< class TGraph, typename TPath, typename TContainer >
    inline typename TGraph::id_type
    least_covered_adjacent( TGraph const& graph, TPath& tail,
                            TContainer const& paths_set )
    {
      typedef TGraph graph_type;
      typedef typename graph_type::id_type id_type;
      typedef typename graph_type::linktype_type linktype_type;

      id_type lc_id = 0;
      std::size_t lc_value = UINTMAX_MAX;
      bool equally_covered = true;
      auto cov =
          [&paths_set]( TPath const& tail ) {
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
              if ( equally_covered && lc_value != UINTMAX_MAX && lc_value != value ) {
                equally_covered = false;
              }
              if ( value < lc_value ) {
                lc_id = to;
                lc_value = value;
              }
              return true;
            } );
      }

      return equally_covered ? 0 : lc_id;
    }

    /**
     *  @brief Get component node rank boundaries of the given graph.
     *
     *  @param[in]  graph The graph.
     *  @return A sorted vector containing the smallest node rank in each component.
     *
     *  NOTE: This method assumes that the input graph is sorted such that node rank
     *  ranges of components are disjoint.
     *
     *  NOTE: This function assumes that the graph is augmented by one path per region
     *        and nothing more.
     */
    template< class TGraph >
    inline std::vector< typename TGraph::rank_type >
    components_ranks( TGraph const& graph )
    {
      std::vector< typename TGraph::rank_type > result;
      graph.for_each_path(
          [&graph, &result]( auto path_rank, auto path_id ) {
            auto sid = *graph.path( path_id ).begin();
            result.push_back( graph.id_to_rank( sid ) );
            return true;
          } );
      std::sort( result.begin(), result.end() );
      return result;
    }
  }  /* --- end of namespace util --- */
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_GRAPH_HPP__ --- */

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
#include <vector>

#include <gum/graph.hpp>
#include <vg/vg.pb.h>
#include <vg/io/stream.hpp>

#include "utils.hpp"


namespace psi {
  namespace util {
    template< class TGraph >
    inline typename TGraph::offset_type
    max_node_len( TGraph const& graph )
    {
      typedef TGraph graph_type;
      typedef typename graph_type::id_type id_type;
      typedef typename graph_type::rank_type rank_type;
      typedef typename graph_type::offset_type offset_type;

      offset_type max = 1;
      graph.for_each_node(
          [&graph, &max]( rank_type rank , id_type id ) {
            auto len = graph.node_length( id );
            if ( max < len ) max = len;
            return true;
          } );
      return max;
    }

    /**
     *  @brief  Compute total number of loci in the subgraph indicated by node ranks
     *          [lower, upper).
     *
     *  @param[in]  graph The input graph.
     *  @param[in]  lower The minimum node rank in the subgraph.
     *  @param[in]  upper The maximum node rank in the subgraph.
     */
    template< class TGraph >
    inline typename TGraph::offset_type
    total_nof_loci( TGraph const& graph, typename TGraph::rank_type lower,
                    typename TGraph::rank_type upper=0 )
    {
      typedef TGraph graph_type;
      typedef typename graph_type::id_type id_type;
      typedef typename graph_type::rank_type rank_type;
      typedef typename graph_type::offset_type offset_type;

      offset_type total = 0;
      graph.for_each_node(
          [&graph, &total, &upper]( rank_type rank, id_type id ) {
            total += graph.node_length( id );
            if ( rank + 1 == upper ) return false;
            return true;
          },
          lower );
      return total;
    }

    template< class TGraph >
    inline typename TGraph::offset_type
    total_nof_loci( TGraph const& graph, gum::Dynamic )
    {
      return total_nof_loci( graph, 1 /* from the first node to the end */ );
    }

    template< class TGraph >
    inline typename TGraph::offset_type
    total_nof_loci( TGraph const& graph, gum::Succinct )
    {
      return gum::util::length_sum( graph.get_node_prop().sequences() );
    }

    template< class TGraph >
    inline typename TGraph::offset_type
    total_nof_loci( TGraph const& graph )
    {
      return total_nof_loci( graph, typename TGraph::spec_type() );
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
                   std::ptrdiff_t chunk_size )
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
     *  @brief  Count the nodes of a subgraph in the graph or of the whole graph.
     *
     *  The subgraph is indicated by the node range [lower, upper).
     */
    template< class TGraph >
    inline typename TGraph::rank_type
    node_count( TGraph const& graph, typename TGraph::rank_type lower=1,
                typename TGraph::rank_type upper=0 )
    {
      if ( upper == 0 ) upper = graph.get_node_count() + 1;
      assert( lower <= graph.get_node_count() && lower > 0 );
      assert( upper <= graph.get_node_count()+1 && upper > lower );
      return upper - lower;
    }

    /**
     *  @brief  Count the edges of a graph component or of the whole graph.
     *
     *  The component is indicated by the node range [lower, upper).
     */
    template< class TGraph >
    inline typename TGraph::rank_type
    edge_count( TGraph const& graph, typename TGraph::rank_type lower=1,
                typename TGraph::rank_type upper=0 )
    {
      typedef typename TGraph::id_type id_type;
      typedef typename TGraph::rank_type rank_type;

      rank_type edge_count = 0;
      graph.for_each_node(
          [&graph, &edge_count, &upper]( rank_type rank, id_type id ){
            edge_count += graph.outdegree( id );
            if ( rank + 1 == upper ) return false;
            return true;
          },
          lower );
      return edge_count;
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

    /**
     *  @brief  Get the adjacency matrix of the graph in CRS format.
     *
     *  @param[in]  graph The graph.
     *  @param[in]  tag  The CRS trait tag.
     *  @param[in]  lower  The lower node rank (inclusive).
     *  @param[in]  upper  The upper node rank (exclusive).
     *  @return  The adjacency matrix.
     *
     *  Compute adjacency matrix of a component in the given `graph` or of the whole
     *  graph. The component is indicated by nodes whose ranks are in the range [lower,
     *  upper). The resulting adjacency matrix is stored in CRS format.
     */
    template< class TGraph, typename TCrsTraits >
    inline typename TCrsTraits::crsMat_t
    adjacency_matrix( TGraph const& graph, TCrsTraits /* tag */,
                      typename TGraph::rank_type lower=1,
                      typename TGraph::rank_type upper=0 )
    {
      typedef TGraph graph_type;
      typedef typename graph_type::id_type id_type;
      typedef typename graph_type::offset_type offset_type;
      typedef typename graph_type::rank_type rank_type;
      typedef typename graph_type::linktype_type linktype_type;

      typedef typename TCrsTraits::lno_t lno_t;
      typedef typename TCrsTraits::size_type size_type;
      typedef typename TCrsTraits::lno_nnz_view_t lno_nnz_view_t;
      typedef typename TCrsTraits::scalar_view_t scalar_view_t;
      typedef typename TCrsTraits::lno_view_t lno_view_t;
      typedef typename TCrsTraits::crsMat_t crsMat_t;

      if ( upper == 0 ) upper = graph.get_node_count() + 1;
      lno_t nrows = total_nof_loci( graph, lower, upper );
      size_type nnz = nrows - node_count( graph, lower, upper ) +
          edge_count( graph, lower, upper );

      lno_nnz_view_t entries( "entries", nnz );
      scalar_view_t values( "values", nnz );
      lno_view_t rowmap( "rowmap", nrows + 1 );

      for ( size_type i = 0; i < nnz; ++i ) values( i ) = 1;  // boolean

      offset_type cursor = 0;
      offset_type start = gum::util::id_to_charorder( graph, graph.rank_to_id( lower ) );
      size_type i = 0;
      size_type irow = 0;
      rowmap( irow++ ) = i;
      graph.for_each_node(
          [&]( rank_type rank, id_type id ) {
            assert( gum::util::id_to_charorder( graph, id ) == cursor + start );
            for ( offset_type offset = 1; offset < graph.node_length( id ); ++offset ) {
              entries( i++ ) = ++cursor;
              rowmap( irow++ ) = i;
            }
            ++cursor;
            graph.for_each_edges_out(
                id,
                [&graph, &entries, &i, start]( id_type to, linktype_type ) {
                  entries( i++ ) = gum::util::id_to_charorder( graph, to ) - start;
                  return true;
                } );
            rowmap( irow++ ) = i;
            if ( rank + 1 == upper ) return false;
            return true;
          },
          lower );
      assert( i == nnz );
      assert( irow == static_cast< unsigned int >( nrows + 1 ) );

      return crsMat_t( "adjacency matrix", nrows, nrows, nnz, values, rowmap, entries );
    }

    /**
     *  @brief Compress a distance index by removing intra-node loci pairs
     *
     *  @param  dindex input distance index
     *  @param  graph underlying graph
     *  @return a mutable compressed distance index of type `TMutableCRSMatrix`
     *
     *  NOTE: The resulting mutable matrix can be assigned to a immutable compressed
     *        matrix afterwards.
     *
     *  NOTE: The input uncompressed distance index is passed by non-const reference,
     *        since containers in const Buffered specialisations cannot be iterated.
     */
    template< typename TMutableCRSMatrix,
              typename TCRSMatrix,
              typename TGraph >
    inline TMutableCRSMatrix
    compress_distance_index( TCRSMatrix& dindex, TGraph const& graph )
    {
      typedef TGraph graph_type;
      typedef typename graph_type::rank_type rank_type;

      typedef TMutableCRSMatrix crsmat_mutable_type;
      typedef TCRSMatrix crsmat_type;
      typedef typename crsmat_type::ordinal_type ordinal_type;

      typename crsmat_mutable_type::entries_type entries;
      typename crsmat_mutable_type::rowmap_type rowmap;
      crsmat_mutable_type::base_type::traits_type::init( entries );
      crsmat_mutable_type::base_type::traits_type::init( rowmap );
      rank_type cnode_rank = 0;  // current node rank
      ordinal_type start = 0;    // row start index
      ordinal_type end;          // row end index
      ordinal_type cloc = 0;     // current node loci index
      ordinal_type nloc = 0;     // next node loci index
      for ( ordinal_type nrow = 0; nrow < dindex.numRows(); ++nrow ) {
        rowmap.push_back( entries.size() );
        if ( nrow == nloc ) {
          ++cnode_rank;
          cloc = nloc;
          nloc += graph.node_length( graph.rank_to_id( cnode_rank ) );
        }
        assert( nrow < nloc );
        end = dindex.rowMap( nrow + 1 );
        for ( ; start < end; ++start ) {
          if ( cloc <= dindex.entry( start ) && dindex.entry( start ) < nloc ) continue;
          else entries.push_back( dindex.entry( start ) );
        }
      }
      rowmap.push_back( entries.size() );
      assert( start == dindex.nnz() );

      return crsmat_mutable_type( dindex.numCols(), std::move( entries ), std::move( rowmap ) );
    }
  }  /* --- end of namespace util --- */
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_GRAPH_HPP__ --- */

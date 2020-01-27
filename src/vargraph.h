/**
 *    @file  vargraph.h
 *   @brief  VarGraph class definition.
 *
 *  This header file contains VarGraph, Paths, and VarGraph iterators class definitions
 *  and interface functions.
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

#ifndef VARGRAPH_H__
#define VARGRAPH_H__

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <unordered_set>
#include <deque>
#include <utility>
#include <functional>
#include <random>

#include "path.h"
#include "graph_iter.h"

#ifdef HAVE_MEMCPY
#  define MACRO_STACK
#  pragma push_macro("HAVE_MEMCPY")
#  undef HAVE_MEMCPY
#endif

#include <xg.hpp>

#ifdef MACRO_STACK
#  pragma pop_macro("HAVE_MEMCPY")
#  undef MACRO_STACK
#endif


namespace grem
{
  class VarGraph : public xg::XG
  {
    public:
      // typedefs
      typedef vg::Node node_type;                                               /**< @brief Node type. */
      typedef std::make_unsigned_t< decltype( vg::Node().id() ) > nodeid_type;  /**< @brief Node ID type. */
      typedef std::size_t rank_type;                                            /**< @brief Node ID type. */
      typedef decltype( vg::Position().offset() ) offset_type;                  /**< @brief Node offset type. */

      using xg::XG::XG;

      // Public methods
      using xg::XG::id_to_rank;

        inline rank_type
      id_to_rank( vg::Position pos ) const
      {
        return this->id_to_rank( pos.node_id() );
      }

        inline bool
      is_branch ( nodeid_type node_id ) const
      {
        if ( this->edges_from_count( node_id ) > 1 ) {
          return true;
        }
        return false;
      }  /* -----  end of method is_branch  ----- */

        inline bool
      is_merge ( nodeid_type node_id ) const
      {
        if ( this->edges_to_count( node_id ) > 1 ) {
          return true;
        }
        return false;
      }  /* -----  end of method is_merge  ----- */

        inline bool
      has_edges_from( nodeid_type node_id ) const
      {
        return this->edges_from_count( node_id ) != 0;
      }

        inline bool
      has_edges_to( nodeid_type node_id ) const
      {
        return this->edges_to_count( node_id ) != 0;
      }

        inline vg::Edge
      get_edge( nodeid_type from, bool from_start, nodeid_type to, bool to_end ) const
      {
        auto type = this->edge_type( from_start, to_end );
        return this->edge_from_encoding( from, to, type );
      }

        static inline constexpr int
      get_overlap( nodeid_type from, bool from_start, nodeid_type to, bool to_end )
      {
        // :TODO:Fri Apr 19 13:57:\@cartoonist: Need overlap?
        return 0;
      }

        static inline constexpr int
      get_overlap( nodeid_type from, nodeid_type to, int type )
      {
        // :TODO:Fri Apr 19 13:57:\@cartoonist: Need overlap?
        return 0;
      }

        inline int
      edges_count_by_id( nodeid_type id, bool is_reverse=false, bool go_left=false ) const
      {
        int counter = 0;
        auto handle = this->get_handle( id, is_reverse );
        this->follow_edges( handle, go_left, [&counter]( vg::handle_t const& h ) {
            ++counter;
            return true;
          });
        return counter;
      }

        inline int
      edges_from_count( nodeid_type id, bool is_reverse=false ) const
      {
        return edges_count_by_id( id, is_reverse, false );
      }

        inline int
      edges_to_count( nodeid_type id, bool is_reverse=false ) const
      {
        return edges_count_by_id( id, is_reverse, true );
      }

        inline std::make_unsigned_t< offset_type >
      get_max_node_len( ) const
      {
        std::make_unsigned_t< offset_type > max = 1;
        for ( std::size_t rank = 1; rank <= this->max_node_rank(); ++rank ) {
          auto&& id = this->rank_to_id( rank );
          if ( max < this->node_length( id ) ) max = this->node_length( id );
        }
        return max;
      }

        inline std::make_unsigned_t< offset_type >
      get_total_nof_loci( ) const
      {
        std::make_unsigned_t< offset_type > total = 0;
        for ( std::size_t rank = 1; rank <= this->max_node_rank(); ++rank ) {
          auto id = this->rank_to_id( rank );
          total += this->node_length( id );
        }
        return total;
      }

      template< typename TNodesIter, typename TEdgesIter >
          inline void
        induced_graph( TNodesIter nbegin, TNodesIter nend,
            TEdgesIter ebegin, TEdgesIter eend,
            vg::Graph* graph ) const
        {
          for ( ; nbegin != nend; ++nbegin ) {
            vg::Node* new_node = graph->add_node();
            new_node->set_id( *nbegin );
            new_node->set_sequence( this->node_sequence( *nbegin ) );
          }

          for ( ; ebegin != eend; ++ebegin ) {
            vg::Edge* new_edge = graph->add_edge();
            new_edge->set_from( std::get< 0 >( *ebegin ) );
            new_edge->set_to( std::get< 1 >( *ebegin ) );
            new_edge->set_from_start( std::get< 2 >( *ebegin ) > 2 );
            new_edge->set_to_end( !( std::get< 2 >( *ebegin ) % 2 ) );
            new_edge->set_overlap( this->get_overlap( std::get< 0 >( *ebegin ),
                  std::get< 1 >( *ebegin ), std::get< 2 >( *ebegin ) ) );
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
      template< typename TNodesIter, typename TEdgesIter >
          inline void
        induced_graph( TNodesIter nbegin, TNodesIter nend,
            TEdgesIter ebegin, TEdgesIter eend,
            std::function< void( vg::Graph& ) > callback, std::ptrdiff_t chunk_size ) const
        {
          auto nodes_l = nbegin;
          auto nodes_r = nodes_l + 1;
          auto edges_l = ebegin;
          auto edges_r = edges_l + 1;
          auto nofmsg = ( nend - nbegin ) / chunk_size + 1;
          for ( unsigned int i = 0; i < nofmsg && nodes_r != nend; ++i ) {
            nodes_r = nodes_l + std::min( chunk_size, nend - nodes_l );
            auto maxid = *( nodes_r - 1 );
            while ( edges_r != eend && std::get< 0 >( *edges_r ) <= maxid ) ++edges_r;
            vg::Graph g;
            this->induced_graph( nodes_l, nodes_r, edges_l, edges_r, &g );
            callback( g );
            nodes_l = nodes_r;
            edges_l = edges_r;
          }
        }
  };

  /* Graph interface functions  ------------------------------------------------ */

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
    inline VarGraph::nodeid_type
  get_random_adjacent( const VarGraph &vargraph, VarGraph::nodeid_type node_id,
     unsigned int seed=0 )
  {
    if ( !vargraph.has_edges_from( node_id ) ) {
      return 0;
    }

    auto fwd_edges = vargraph.edges_from(node_id);

    std::random_device rd;  // Will be used to obtain a seed for the random no. engine
    if ( seed == 0 ) seed = rd(); // use random_device to generate a seed if seed is not provided
    std::mt19937 gen( seed );  // Standard mersenne_twister_engine seeded with seed
    std::uniform_int_distribution<> dis(0, fwd_edges.size() - 1);

    return fwd_edges[ dis(gen) ].to();
  }  /* -----  end of function get_random_adjacent  ----- */

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
  template< typename TContainer >
      inline VarGraph::nodeid_type
    least_covered_adjacent( const VarGraph& vargraph, VarGraph::nodeid_type node_id,
        const TContainer& paths_set )
    {
      if ( vargraph.has_edges_from( node_id ) ) {
        auto fwd_edges = vargraph.edges_from( node_id );

        VarGraph::nodeid_type first_adj_id = ( *fwd_edges.begin() ).to();

        VarGraph::nodeid_type lc_node_id = first_adj_id;
        unsigned int lc = get_path_coverage( first_adj_id, paths_set );
        bool equally_covered = true;

        for ( auto e_itr = ++fwd_edges.begin(); e_itr != fwd_edges.end(); ++e_itr ) {
          const VarGraph::nodeid_type& next_node = ( *e_itr ).to();
          unsigned int next_node_cov = get_path_coverage( next_node, paths_set );

          if ( equally_covered && lc != next_node_cov ) {
            equally_covered = false;
          }

          if ( next_node_cov < lc ) {
            lc = next_node_cov;
            lc_node_id = next_node;
          }
        }

        if ( !equally_covered ) {
          return lc_node_id;
        }
      }

      return 0;
    }  /* -----  end of template function least_covered_adjacent  ----- */

  template< typename TPath, typename TContainer >
      inline VarGraph::nodeid_type
    least_covered_adjacent( const VarGraph& vargraph, TPath& tail,
        const TContainer& paths_set )
    {
      if ( !tail.empty() && vargraph.has_edges_from( tail.back() ) ) {
        auto fwd_edges = vargraph.edges_from( tail.back() );

        VarGraph::nodeid_type first_adj_id = ( *fwd_edges.begin() ).to();
        VarGraph::nodeid_type lc_node_id = first_adj_id;
        tail.push_back( first_adj_id );         /* XXX: should be popped later! */
        unsigned int lc = get_path_coverage( tail.begin(), tail.end(), paths_set );
        tail.pop_back();

        bool equally_covered = true;
        for ( auto e_itr = ++fwd_edges.begin(); e_itr != fwd_edges.end(); ++e_itr ) {
          const VarGraph::nodeid_type& next_node = ( *e_itr ).to();
          tail.push_back( next_node );         /* XXX: should be popped later! */
          unsigned int next_node_cov = get_path_coverage( tail.begin(), tail.end(), paths_set );
          tail.pop_back();

          if ( equally_covered && lc != next_node_cov ) equally_covered = false;

          if ( next_node_cov < lc ) {
            lc = next_node_cov;
            lc_node_id = next_node;
          }
        }

        if ( !equally_covered ) {
          return lc_node_id;
        }
      }

      return 0;
    }

  /* END OF graph interface functions  ----------------------------------------- */

  /* GRAPH ITERATORS  ============================================================ */

  /* Traits template specialization  ------------------------------------------- */

  /**
   *  @brief  Breadth-first search graph iterator tag.
   *
   *  Specialization of generic graph iterator tag BFS for VarGraph.
   */
  template< >
    struct GraphIterTraits < VarGraph, BFS >
    {
      typedef VarGraph::nodeid_type Value;
      typedef Value Level;
      typedef std::deque< std::pair< Value, Level > > TContainer;

      struct pair_hash
      {
        inline std::size_t operator()(const std::pair< Value, Level > & v) const
        {
          return std::hash< Value >()(v.first);
        }
      };

      struct pair_pred
      {
        inline bool operator()
          (const std::pair< Value, Level > & v,
           const std::pair< Value, Level > & u)
          const noexcept
        {
          if (v.first == u.first) return true;
          else return false;
        }
      };

      typedef std::unordered_set< TContainer::value_type, pair_hash, pair_pred > TSet;
      typedef struct {
        std::size_t lb_visited_rank;  /**< @brief lower-bound for rank of visited nodes. */
      } TState;
      typedef void* TParameter;
      constexpr static const TParameter param_default = nullptr;
    };

  /**
   *  @brief  Backtracker graph iterator tag.
   *
   *  Specialization of generic graph iterator tag Backtracker for VarGraph.
   */
  template< >
    struct GraphIterTraits< VarGraph, Backtracker > {
      typedef VarGraph::nodeid_type Value;
      typedef Value Level;
      typedef std::deque< std::pair< Value, Value > > TContainer;
      typedef void* TSet;
      typedef struct {
        Value start;                            /**< @brief Start node ID. */
        Value buffer;                           /**< @brief Buffer node ID. 0=nothing */
        bool end;                               /**< @brief End flag. */
      } TState;
      typedef void* TParameter;
      constexpr static const TParameter param_default = nullptr;
    };  /* ----------  end of struct Backtracker  ---------- */

  /**
   *  @brief  Haplotyper graph iterator tag.
   *
   *  Specialization of generic graph iterator tag Haplotyper for VarGraph.
   */
  template< typename TSpec >
    struct GraphIterTraits< VarGraph, Haplotyper< TSpec > > {
      typedef VarGraph::nodeid_type Value;
      typedef std::unique_ptr< Path< VarGraph, Dynamic > > TContainer;
      /**< @brief Set of visited paths. */
      typedef std::vector< Path< VarGraph, Haplotype > > TSet;
      typedef TSet::size_type Level;        /**< @brief No. of selected path so far. */
      typedef struct {
        Value start;                        /**< @brief Start node ID. */
        bool end;                           /**< @brief End flag. */
        std::unique_ptr< Path< VarGraph, Haplotype > > current_path;
        unsigned int setback;
        unsigned int entropy;
      } TState;
      typedef unsigned int TParameter;
      static const TParameter param_default = 0;
    };  /* ----------  end of struct HaplotyperIter  ---------- */

  /**
   *  @brief  Haplotyper graph iterator tag.
   *
   *  Specialization of generic graph iterator tag Haplotyper< Random > for VarGraph.
   */
  template< >
    struct GraphIterTraits< VarGraph, Haplotyper< Random > > {
      typedef VarGraph::nodeid_type Value;
      typedef Value Level;
      typedef void* TContainer;
      typedef void* TSet;
      typedef struct {
        Value start;                            /**< @brief Start node ID. */
        Level level;                            /**< @brief Level. */
        bool end;                               /**< @brief End flag. */
      } TState;
      typedef unsigned int TParameter;
      constexpr static const TParameter param_default = 0;
    };

  /* END OF tags template specialization  -------------------------------------- */

}  /* -----  end of namespace grem  ----- */

namespace seqan {
  /**
   *  @brief  Iterator template class specialization for VarGraph.
   *
   *  This class extends existing Iterator class in seqan namespace.
   */
  template< typename TSpec >
    class Iterator< grem::VarGraph, TSpec >
    {
      public:
        typedef grem::GraphIter< grem::VarGraph, TSpec > Type;
        /* ====================  TYPEDEFS      ======================================= */
    };  /* ----------  end of template class Iterator  ---------- */

  template< typename TSpec >
    struct Value< grem::GraphIter< grem::VarGraph, TSpec > > {
      using Type = typename grem::GraphIter< grem::VarGraph, TSpec >::TTraits::Value;
    };

  template< typename T >
    struct Level;

  template< typename TSpec >
    struct Level< grem::GraphIter< grem::VarGraph, TSpec > > {
      using Type = typename grem::GraphIter< grem::VarGraph, TSpec >::TTraits::Level;
    };
}  /* -----  end of namespace seqan  ----- */

namespace grem {
  /* VarGraph iterators template specialization  --------------------------------- */

  /* BFS template specialization  ------------------------------------------------ */

  /* Internal functions specialization. */

  /**
   *  @brief  Search the graph for next unvisited node.
   *
   *  @return the node ID of the next visited node or `0` if all are visited.
   *
   *  The lower-bound for visited rank is also updated so that any node with smaller
   *  rank is visited. So in order to find the next unvisited node, we should search
   *  among nodes with higher ranks.
   */
  template< >
      inline typename seqan::Value< GraphIter< VarGraph, BFS > >::Type
    GraphIter< VarGraph, BFS >::next_unvisited( )
    {
      TTraits::Value next_node = 0;
      while ( this->state.lb_visited_rank <= this->vargraph_ptr->max_node_rank() ) {
        next_node = this->vargraph_ptr->rank_to_id( this->state.lb_visited_rank );
        if ( !(*this)[ next_node ] ) break;
        ++this->state.lb_visited_rank;
      }
      if ( this->state.lb_visited_rank <= this->vargraph_ptr->max_node_rank() ) {
        return next_node;
      }
      else {
        return 0;
      }
    }  /* -----  end of template function get_next_unvisited  ----- */

  /* Interface functions specialization. */

  template< >
      inline bool
    at_end( GraphIter< VarGraph, BFS >& it )
    {
      return it.visiting_buffer.empty();
    }

  template< >
      inline void
    begin( GraphIter< VarGraph, BFS >& it, const VarGraph* g,
        typename seqan::Value< GraphIter< VarGraph, BFS > >::Type start,
        typename GraphIter< VarGraph, BFS >::parameter_type )
    {
      if ( start == 0 ) start = g->rank_to_id( 1 );

      it.state.lb_visited_rank = 1;
      if ( g->id_to_rank( start ) == 1 ) {
        ++it.state.lb_visited_rank;
      }

      it.vargraph_ptr = g;
      it.visiting_buffer.push_back( std::make_pair( start, 0 ) );
      it.visited.insert( std::make_pair( start, 0 ) );
      it.itr_value = it.visiting_buffer.front().first;
    }

  template< >
      inline void
    go_begin( GraphIter< VarGraph, BFS >& it,
        typename seqan::Value< GraphIter< VarGraph, BFS > >::Type start,
        typename GraphIter< VarGraph, BFS >::parameter_type )
    {
      if ( start == 0 ) start = it.vargraph_ptr->rank_to_id( 1 );

      it.state.lb_visited_rank = 1;
      if ( it.vargraph_ptr->id_to_rank( start ) == 1 ) {
        ++it.state.lb_visited_rank;
      }

      it.visiting_buffer.clear();
      it.visiting_buffer.push_back( std::make_pair( start, 0 ) );
      it.visited.clear();
      it.visited.insert( std::make_pair( start, 0 ) );
      it.itr_value = it.visiting_buffer.front().first;
    }  /* -----  end of template function go_begin  ----- */

  template< >
      inline typename seqan::Level< GraphIter< VarGraph, BFS > >::Type
    level( GraphIter< VarGraph, BFS >& it )
    {
      if ( !it.visiting_buffer.empty() ) {
        return it.visiting_buffer.front().second;
      }

      throw std::runtime_error( "invalid level query on the end of iterator." );
    }

  /* Member functions specialization. */

  template< >
      inline GraphIter< VarGraph, BFS >&
    GraphIter< VarGraph, BFS >::operator++( )
    {
      typedef typename seqan::Level< GraphIter< VarGraph, BFS > >::Type TLevel;
      typedef typename seqan::Value< GraphIter< VarGraph, BFS > >::Type TValue;

      if ( at_end( *this ) && this->raise_on_end ) {
        throw std::range_error( "The iterator has reached the end." );
      }

      if ( this->visiting_buffer.empty() ) return *this;

      TLevel plevel = level( *this );
      auto edges = this->vargraph_ptr->edges_from( this->itr_value );
      for ( auto it = edges.begin(); it != edges.end(); ++it ) {
        TValue adj_node = (*it).to();
        if ( !(*this)[ adj_node ] ) {
          this->visiting_buffer.push_back( std::make_pair( adj_node, plevel + 1 ) );
          this->visited.insert( std::make_pair( adj_node, plevel + 1 ) );
        }
      }
      this->visiting_buffer.pop_front();

      if ( !this->visiting_buffer.empty() ) {
        this->itr_value = this->visiting_buffer.front().first;
      }
      else {
        this->itr_value = this->next_unvisited();
        if ( this->itr_value != 0 ) {
          this->visiting_buffer.push_back( std::make_pair( this->itr_value, 0 ) );
          this->visited.insert( std::make_pair( this->itr_value, 0 ) );
        }
      }

      if ( this->itr_value != 0 &&
          this->state.lb_visited_rank == this->vargraph_ptr->id_to_rank( this->itr_value ) ) {
        ++this->state.lb_visited_rank;
      }

      return *this;
    }

  /**
   *  @brief  Check whether a node ID is visited by the BFS iterator or not.
   *
   *  @param  id The ID of the query node.
   *  @return `true` if node is visited by BFS; otherwise `false`.
   *
   *  It queries visited set for the node ID. Level doesn't matter so it is set to zero
   *  (see method `pair_pred`).
   */
  template< >
  template< typename TId >
      inline bool
    GraphIter< VarGraph, BFS >::operator[]( const TId& id )
    {
      return this->visited.find( std::make_pair( id, 0 ) ) != this->visited.end();
    }  /* -----  end of method GraphIter< VarGraph, BFS >::operator[]  ----- */

  /* END OF BFS template specialization  ----------------------------------------- */

  /* Backtracker template specialization  ---------------------------------------- */

  /* Interface functions specialization. */
  template< >
      inline bool
    at_end( GraphIter< VarGraph, Backtracker >& it )
    {
      return it.state.end;
    }  /* -----  end of template function at_end  ----- */

  template< >
      inline void
    begin( GraphIter< VarGraph, Backtracker >& it, const VarGraph* g,
        typename seqan::Value< GraphIter< VarGraph, Backtracker > >::Type start,
        typename GraphIter< VarGraph, Backtracker >::parameter_type )
    {
      if ( start == 0 ) start = g->rank_to_id( 1 );

      it.vargraph_ptr = g;
      it.itr_value = start;
      it.state.buffer = 0;
      it.state.end = false;
      it.state.start = start;
    }  /* -----  end of template function begin  ----- */

  template< >
      inline void
    go_begin( GraphIter< VarGraph, Backtracker >& it,
        typename seqan::Value< GraphIter< VarGraph, Backtracker > >::Type start,
        typename GraphIter< VarGraph, Backtracker >::parameter_type )
    {
      if ( start == 0 ) start = it.state.start;  // Re-use stored start node.

      it.itr_value = start;
      it.state.buffer = 0;  // Re-set buffer.
      it.state.end = false;  // Re-set at-end flag.
      it.visiting_buffer.clear();
    }  /* -----  end of template function go_begin  ----- */

  /* Member functions specialization. */

  template< >
      inline GraphIter< VarGraph, Backtracker >&
    GraphIter< VarGraph, Backtracker >::operator++( )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Backtracker > >::Type TValue;

      if ( at_end( *this ) && this->raise_on_end ) {
        throw std::range_error( "The iterator has reached the end." );
      }

      if ( this->state.buffer != 0 ) {                             // Any node buffered?
        this->itr_value = this->state.buffer;                      // Use it.
        this->state.buffer = 0;                                    // Clear up buffer.
      }
      else {                                                  // else
        TValue cnode_id = this->itr_value;
        if ( this->vargraph_ptr->has_edges_from( cnode_id ) ) {  // Any forward edge?
          // Go forward.
          this->itr_value = this->vargraph_ptr->edges_from( cnode_id )[ 0 ].to();
          // On each branch nodes enqueue other branches for traversing later.
          auto edges = this->vargraph_ptr->edges_from( cnode_id );
          for ( int i = edges.size() - 1; i >= 1; --i ) {
            this->visiting_buffer.push_back( std::make_pair( cnode_id, edges[i].to() ) );
          }
        }
        else {
          this->state.end = true;
        }
      }

      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Backtracker >::operator++  ----- */

  template< >
      inline GraphIter< VarGraph, Backtracker >&
    GraphIter< VarGraph, Backtracker >::operator--( )
    {
      if ( this->state.buffer != 0 ) {                             // Any node buffered?
        while (                // Remove all buffered branches of the current node.
            !this->visiting_buffer.empty() &&
            this->visiting_buffer.back().first == this->itr_value ) {
          this->visiting_buffer.pop_back();
        }
        this->state.buffer = 0;
      }

      if ( !this->visiting_buffer.empty() ) {                 // Go back in buffer.
        this->itr_value = this->visiting_buffer.back().first;
        this->state.buffer = this->visiting_buffer.back().second;
        this->visiting_buffer.pop_back();
        this->state.end = false;  // Reset at-end flag.
      }
      else {
        this->state.end = true;  // Set at-end flag.
      }

      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Backtracker >::operator--  ----- */

  /* END OF Backtracker template specialization  --------------------------------- */

  /* Haplotyper template specialization  ----------------------------------------- */

  /* Interface functions specialization. */
  template< >
      inline bool
    at_end( GraphIter< VarGraph, Haplotyper<> >& it )
    {
      return it.state.end;
    }  /* -----  end of template function at_end  ----- */

  template< >
      inline void
    begin( GraphIter< VarGraph, Haplotyper<> >& it, const VarGraph* g,
        typename seqan::Value< GraphIter< VarGraph, Haplotyper<> > >::Type start,
        typename GraphIter< VarGraph, Haplotyper<> >::parameter_type )
    {
      if ( start == 0 ) start = g->rank_to_id( 1 );

      it.vargraph_ptr = g;
      it.itr_value = start;
      it.visiting_buffer = std::make_unique< Path< VarGraph, Dynamic > >( g );
      it.state.start = start;
      it.state.end = false;
      it.state.current_path = std::make_unique< Path< VarGraph, Haplotype > >( g );
      it.state.current_path->push_back( it.itr_value );
      it.state.setback = 0;
      it.state.entropy = 1;
    }  /* -----  end of template function begin  ----- */

  template< >
      inline void
    go_begin( GraphIter< VarGraph, Haplotyper<> >& it,
        typename seqan::Value< GraphIter< VarGraph, Haplotyper<> > >::Type start,
        typename GraphIter< VarGraph, Haplotyper<> >::parameter_type )
    {
      if ( start == 0 ) start = it.state.start;  // Re-use start node.

      it.itr_value = start;
      it.state.start = start;
      it.visiting_buffer->clear();
      it.state.end = false;  // Re-set at-end flag.
      it.visited.clear();
      it.state.current_path->clear();
      it.state.current_path->push_back( it.itr_value );
      it.state.setback = 0;
      it.state.entropy = 1;
    }  /* -----  end of template function go_begin  ----- */

  template< >
      inline typename seqan::Level< GraphIter< VarGraph, Haplotyper<> > >::Type
    level( GraphIter< VarGraph, Haplotyper<> >& it )
    {
      return it.visited.size();
    }

  /* Member functions specialization. */

  template< >
      inline void
    GraphIter< VarGraph, Haplotyper<> >::set_setback( )
    {
      this->state.setback = this->visited.size();
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper<> >::set_setback  ----- */

  /**
   *  A setback path is a sequence of last 's' nodes of currently generating haplotype.
   *  An unvisited setback path is a path that does not occur as a subset of any
   *  previously generated haplotypes. This function search in adjacent nodes set for
   *  a node that together with 's-1' previously selected nodes forms an unvisited
   *  setback in order to cover more k-mers from all paths in the graph; i.e. generating
   *  more diverse haplotypes.
   */
  template< >
      inline GraphIter< VarGraph, Haplotyper<> >&
    GraphIter< VarGraph, Haplotyper<> >::operator++( )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Haplotyper<> > >::Type TValue;

      if ( at_end( *this ) && this->raise_on_end ) {
        throw std::range_error( "The iterator has reached the end." );
      }

      if ( !this->vargraph_ptr->has_edges_from( this->itr_value ) ) {
        this->state.end = true;
        return *this;
      }

      if ( this->state.setback > 1 ) {
        while ( this->visiting_buffer->size() != 0 &&
            ( this->state.entropy > this->state.setback ) ) {
          this->state.entropy /= this->vargraph_ptr->edges_from_count( this->visiting_buffer->front() );
          this->visiting_buffer->pop_front();
        }
      }

      TValue next_candidate = 0;
      auto fwd_edges = this->vargraph_ptr->edges_from( this->itr_value );
      if ( this->state.setback == 0 || fwd_edges.size() == 1 ) {
        next_candidate = fwd_edges[0].to();
      }
      else {
        // Search for a forward node such that the setback path is not in previous paths.
        do {
          for ( auto e : fwd_edges ) {
            this->visiting_buffer->push_back( e.to() );
            if ( (*this)[ *this->visiting_buffer ] ) {  // Visited?
              this->visiting_buffer->pop_back();
              continue;                             // Next edge.
            }
            this->visiting_buffer->pop_back();       // No change to the iterator state.
            next_candidate = e.to();                // Found!
            break;
          }
        } while ( this->state.setback == 1 &&
            next_candidate == 0 &&
            this->visiting_buffer->empty() &&
            ( this->visiting_buffer->push_back( this->itr_value ), true ) );
        if ( this->state.setback == 1 && !this->visiting_buffer->empty() ) this->visiting_buffer->pop_back();
      }
      // If no unvisited setback path found, use a node with least path coverage.
      if ( next_candidate == 0 ) {
        next_candidate = least_covered_adjacent( *this->vargraph_ptr,
            *this->visiting_buffer, this->visited );
      }
      // If all forward edges are visited, pick one randomly with uniform distribution.
      if ( next_candidate == 0 ) {
        next_candidate =
          get_random_adjacent( ( *this->vargraph_ptr ),  this->itr_value );
      }

      this->itr_value = next_candidate;
      if ( this->state.setback > 1 ) {
        this->visiting_buffer->push_back( this->itr_value );
        this->state.entropy *= this->vargraph_ptr->edges_from_count( this->itr_value );
      }
      this->state.current_path->push_back( this->itr_value );

      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper<> >::operator++  ----- */

  template< >
      inline GraphIter< VarGraph, Haplotyper<> >&
    GraphIter< VarGraph, Haplotyper<> >::operator--( int )
    {
      this->itr_value = this->state.start;    // Reset the iterator to the start node.
      this->visiting_buffer->clear();
      this->state.entropy = 1;
      if ( this->state.setback > 1 ) {
        this->visiting_buffer->push_back( this->itr_value );
        this->state.entropy *= this->vargraph_ptr->edges_from_count( this->itr_value );
      }
      this->state.end = false;                // Reset at-end flag.
      this->state.current_path->clear();
      this->state.current_path->push_back( this->itr_value );
      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper<> >::operator--  ----- */

  template< >
      inline GraphIter< VarGraph, Haplotyper<> >&
    GraphIter< VarGraph, Haplotyper<> >::operator--( )
    {
      this->visited.push_back( std::move( *this->state.current_path ) );
      this->set_setback();
      (*this)--;
      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper<> >::operator--  ----- */

  /**
   *  @brief  Check if the given path is present in paths generated so far.
   *
   *  @param  path A container of node IDs indicating nodes in a path.
   *  @return `true` if the path is present; `false` otherwise.
   *
   *  Check whether the given path is generated before or not.
   */
  template< >
  template< typename TContainer >
      inline bool
    GraphIter< VarGraph, Haplotyper<> >::operator[]( const TContainer& path )
    {
      // :TODO:Sun Apr 01 00:01:\@cartoonist: Unordered check!
      // ... `PathSet< Path< TGraph, Compact >, InMemory >` should be used for
      // ... `this->visited`.
      return covered_by( path.begin(), path.end(), this->visited );
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper<> >::operator[]  ----- */

  /* Interface functions specialization. */

  template< >
      inline bool
    at_end( GraphIter< VarGraph, Haplotyper< Local > >& it )
    {
      return it.state.end;
    }  /* -----  end of template function at_end  ----- */

  template< >
      inline void
    begin( GraphIter< VarGraph, Haplotyper< Local > >& it, const VarGraph* g,
        typename seqan::Value< GraphIter< VarGraph, Haplotyper< Local > > >::Type start,
        typename GraphIter< VarGraph, Haplotyper< Local > >::parameter_type p )
    {
      if ( start == 0 ) start = g->rank_to_id( 1 );
      if ( p == 0 ) {
        throw std::runtime_error( "Parameter value of Local Haplotyper cannot be zero" );
      }

      it.vargraph_ptr = g;
      it.itr_value = start;
      it.visiting_buffer = std::make_unique< Path< VarGraph, Dynamic > >( g );
      it.state.start = start;
      it.state.end = false;
      it.state.current_path = std::make_unique< Path< VarGraph, Haplotype > >( g );
      it.state.current_path->push_back( it.itr_value );
      it.state.setback = 0;
      it.param = p;
    }  /* -----  end of template function begin  ----- */

  template< >
      inline void
    go_begin( GraphIter< VarGraph, Haplotyper< Local > >& it,
        typename seqan::Value< GraphIter< VarGraph, Haplotyper< Local > > >::Type start,
        typename GraphIter< VarGraph, Haplotyper< Local > >::parameter_type p )
    {
      if ( start == 0 ) start = it.state.start;  // Re-use start node.
      if ( p == 0 ) p = it.param;  // Re-use parameter value.

      it.itr_value = start;
      it.state.start = start;
      it.visiting_buffer->clear();
      it.state.end = false;  // Re-set at-end flag.
      it.visited.clear();
      it.state.current_path->clear();
      it.state.current_path->push_back( it.itr_value );
      it.state.setback = 0;
      it.param = p;
    }  /* -----  end of template function go_begin  ----- */

  template< >
      inline typename seqan::Level< GraphIter< VarGraph, Haplotyper< Local > > >::Type
    level( GraphIter< VarGraph, Haplotyper< Local > >& it )
    {
      return it.visited.size();
    }

  /* Member functions specialization. */

  template< >
      inline void
    GraphIter< VarGraph, Haplotyper< Local > >::set_setback( )
    {
      if ( this->visited.size() == 0 ) {
        this->state.setback = 0;
      }
      else {
        this->state.setback = 2 * ceil( log2( this->visited.size() + 1 ) ) - 1;
      }
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper<> >::set_setback  ----- */

  /**
   *  A setback path is a sequence of last 's' nodes of currently generating haplotype
   *  encompassing last k basepair. An unvisited setback path is a path that does
   *  not occur as a subset of any previously generated haplotypes. This function search
   *  in adjacent nodes set for a node that together with 's-1' previously selected nodes
   *  forms an unvisited setback in order to cover more k-mers from all paths in the graph;
   *  i.e. generating more diverse haplotypes.
   */
  template< >
      inline GraphIter< VarGraph, Haplotyper< Local > >&
    GraphIter< VarGraph, Haplotyper< Local > >::operator++( )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Haplotyper< Local > > >::Type TValue;

      if ( at_end( *this ) && this->raise_on_end ) {
        throw std::range_error( "The iterator has reached the end." );
      }

      if ( !this->vargraph_ptr->has_edges_from( this->itr_value ) ) {
        this->state.end = true;
        return *this;
      }

      if ( this->state.setback != 0 ) {
        rtrim_front_by_len( *this->visiting_buffer, this->param - 1 );
      }

      TValue next_candidate = 0;
      auto fwd_edges = this->vargraph_ptr->edges_from( this->itr_value );
      if ( this->state.setback == 0 || fwd_edges.size() == 1 ) {
        next_candidate = fwd_edges[0].to();
      }
      else {
        // Search for a forward node such that the setback path is not in previous paths.
        for ( auto e : fwd_edges ) {
          this->visiting_buffer->push_back( e.to() );
          if ( (*this)[ *this->visiting_buffer ] ) {  // Visited?
            this->visiting_buffer->pop_back();
            continue;                             // Next edge.
          }
          this->visiting_buffer->pop_back();      // No change to the iterator state.
          next_candidate = e.to();                // Found!
          break;
        }
      }
      // If no unvisited setback path found, use a node with least path coverage.
      if ( next_candidate == 0 ) {
        next_candidate = least_covered_adjacent( *this->vargraph_ptr,
            *this->visiting_buffer, this->visited );
      }
      // If all forward edges are visited, pick one randomly with uniform distribution.
      if ( next_candidate == 0 ) {
        next_candidate =
          get_random_adjacent( ( *this->vargraph_ptr ),  this->itr_value );
      }

      this->itr_value = next_candidate;
      if ( this->state.setback != 0 ) {
        this->visiting_buffer->push_back( this->itr_value );
      }
      this->state.current_path->push_back( this->itr_value );

      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper< Local > >::operator++  ----- */

  template< >
      inline GraphIter< VarGraph, Haplotyper< Local > >&
    GraphIter< VarGraph, Haplotyper< Local > >::operator--( int )
    {
      this->itr_value = this->state.start;    // Reset the iterator to the start node.
      this->visiting_buffer->clear();
      if ( this->state.setback != 0 ) {
        this->visiting_buffer->push_back( this->itr_value );
      }
      this->state.end = false;                // Reset at-end flag.
      this->state.current_path->clear();
      this->state.current_path->push_back( this->itr_value );
      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper< Local > >::operator--  ----- */

  template< >
      inline GraphIter< VarGraph, Haplotyper< Local > >&
    GraphIter< VarGraph, Haplotyper< Local > >::operator--( )
    {
      this->visited.push_back( std::move( *this->state.current_path ) );
      this->set_setback();
      (*this)--;
      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper< Local > >::operator--  ----- */

  /**
   *  @brief  Check if the given path is present in paths generated so far.
   *
   *  @param  path A container of node IDs indicating nodes in a path.
   *  @return `true` if the path is present; `false` otherwise.
   *
   *  Check whether the given path is generated before or not.
   */
  template< >
  template< typename TContainer >
      inline bool
    GraphIter< VarGraph, Haplotyper< Local > >::operator[]( const TContainer& path )
    {
      // :TODO:Sun Apr 01 00:01:\@cartoonist: Unordered check!
      // ... `PathSet< Path< TGraph, Compact >, InMemory >` should be used for
      // ... `this->visited`.
      return covered_by( path.begin(), path.end(), this->visited );
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper< Local > >::operator[]  ----- */

  /* Interface functions specialization. */
  template< >
      inline bool
    at_end( GraphIter< VarGraph, Haplotyper< Random > >& it )
    {
      return it.state.end;
    }  /* -----  end of template function at_end  ----- */

  template< >
      inline void
    begin( GraphIter< VarGraph, Haplotyper< Random > >& it, const VarGraph* g,
        typename seqan::Value< GraphIter< VarGraph, Haplotyper< Random > > >::Type start,
        typename GraphIter< VarGraph, Haplotyper< Random > >::parameter_type p )
    {
      if ( start == 0 ) start = g->rank_to_id( 1 );

      it.vargraph_ptr = g;
      it.itr_value = start;
      it.state.start = start;
      it.state.level = 1;
      it.state.end = false;
      it.param = p;
    }  /* -----  end of template function begin  ----- */

  template< >
      inline void
    go_begin( GraphIter< VarGraph, Haplotyper< Random > >& it,
        typename seqan::Value< GraphIter< VarGraph, Haplotyper< Random > > >::Type start,
        typename GraphIter< VarGraph, Haplotyper< Random > >::parameter_type p )
    {
      if ( start == 0 ) start = it.state.start;  // Re-use start node.

      it.itr_value = start;
      it.state.start = start;
      it.state.level = 1;
      it.state.end = false;  // Re-set at-end flag.
      it.param = p;
    }  /* -----  end of template function go_begin  ----- */

  template< >
      inline typename seqan::Level< GraphIter< VarGraph, Haplotyper< Random > > >::Type
    level( GraphIter< VarGraph, Haplotyper< Random > >& it )
    {
      return it.state.level;
    }

  /* Member functions specialization. */

  template< >
      inline GraphIter< VarGraph, Haplotyper< Random > >&
    GraphIter< VarGraph, Haplotyper< Random > >::operator++( )
    {
      if ( at_end( *this ) && this->raise_on_end ) {
        throw std::range_error( "The iterator has reached the end." );
      }

      if ( !this->vargraph_ptr->has_edges_from( this->itr_value ) ) {
        this->state.end = true;
        return *this;
      }

      this->itr_value =
        get_random_adjacent( *this->vargraph_ptr,  this->itr_value, this->param );
      ++this->state.level;

      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper< Random > >::operator++  ----- */

  template< >
      inline GraphIter< VarGraph, Haplotyper< Random > >&
    GraphIter< VarGraph, Haplotyper< Random > >::operator--( )
    {
      this->itr_value = this->state.start;
      this->state.level = 1;
      this->state.end = false;
      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper< Local > >::operator--  ----- */

  /* END OF Haplotyper template specialization  ---------------------------------- */

  /* Haplotyper iterator interface functions  ------------------------------------ */

  /**
   *  @brief  Extend a path to length k using a graph iterator.
   *
   *  @param  path The extended path.
   *  @param  iter The graph iterator.
   *  @param  k The extension length.
   *
   *  Add nodes using graph iterator until the length of path reach k or longer.
   */
  template< typename TPathSpec, typename TIterSpec >
      inline void
    extend_to_k( Path< VarGraph, TPathSpec >& path,
        grem::GraphIter< VarGraph, TIterSpec >& iter,
        unsigned int k )
    {
      while ( !at_end( iter ) && path.get_sequence_len() < k ) {
        add_node( path, *iter );
        ++iter;
      }
      iter.operator*();  /**< @brief Trigger the exception if it's at end and raise_on_end is true. */
    }

  /**
   *  @overload Prevent using BFS for extension.
   */
  template< typename TPathSpec >
      inline void
    extend_to_k( Path< VarGraph, TPathSpec >&,
        grem::GraphIter< VarGraph, BFS >&,
        unsigned int )
    {
      throw std::runtime_error( "Cannot be used by BFS iterator." );
    }

  /**
   *  @brief  Simulate a unique haplotype.
   *
   *  @param[out]  haplotype The simulated haplotype as a Path.
   *  @param[in,out]  iter Haplotyper graph iterator.
   *  @param[in]  tries Number of tries if the generated haplotype is not unique.
   *
   *  This function gets a Haplotyper graph iterator and generate a unique haplotype
   *  if available. The input Haplotyper iterator stores required information of the
   *  previous simulated haplotypes for which the iterator is used. So, in order to
   *  simulate multiple unique haplotypes use the same iterator as the input. It tries
   *  `tries` times to generated a unique haplotype.
   */
  template< typename TSpec >
      inline void
    get_uniq_full_haplotype( Path< VarGraph >& haplotype,
        GraphIter< VarGraph, Haplotyper< TSpec > >& iter,
        int tries=0 )
    {
      do {
        while ( !at_end( iter ) ) {
          add_node( haplotype, *iter );
          ++iter;
        }
        if ( tries-- && iter[ haplotype.get_nodes() ] ) {
          iter--;  // discard the traversed path and reset the Haplotyper iterator.
        }
        else {
          --iter;  // save the traversed path and reset the Haplotyper iterator.
          break;
        }
        /* trying again */
      } while ( true );
    }

  /**
   *  @brief  Add a simulated haplotype to the given pathset.
   *
   *  @param  paths The pathset.
   *  @param[in,out]  iter Haplotyper graph iterator.
   *  @param[in]  tries Number of tries if the generated haplotype is not unique.
   *
   *  @overload for `PathSet`.
   */
  template< typename TPathSet, typename TSpec >
      inline void
    get_uniq_full_haplotype( TPathSet& paths,
        GraphIter< VarGraph, Haplotyper< TSpec > >& iter,
        int tries=0 )
    {
      Path< VarGraph > haplotype( iter.get_vargraph() );
      get_uniq_full_haplotype( haplotype, iter, tries );
      if ( length( haplotype ) != 0 ){
        paths.push_back( std::move( haplotype ) );
      }
    }

  template< typename TPathSet, typename TSpec >
      inline void
    get_uniq_patches( TPathSet& paths,
        GraphIter< VarGraph, Haplotyper< TSpec > >& iter,
        unsigned int k )
    { // :TODO:Fri Dec 01 21:26:\@cartoonist: the length of pre-context sequence won't be k all the time.
      iter.raise_on_end = true;
      Path< VarGraph > patch( iter.get_vargraph() );
      Path< VarGraph, Dynamic > frontier( iter.get_vargraph() );
      typename Path< VarGraph, Dynamic >::nodes_type::value_type marked;
      try {
        while ( true ) {
          marked = 0;
          if ( !frontier.empty() ) marked = frontier.get_nodes().back();
          // Bootstrap.
          if ( !marked ) extend_to_k( frontier, iter, k );
          else extend_to_k( frontier, iter,
              2*k + frontier.get_sequence_len() - frontier.get_seqlen_tail() );
          // Check the next patch distance to merge with previous patch if is less than k.
          if ( !patch.empty() && iter[ frontier.get_nodes() ] ) {
            patch.set_right_by_len( k - 1 );
            paths.push_back( std::move( patch ) );
            clear( patch );
            rtrim_front_by_len( frontier, k, true );
          }
          else if ( !patch.empty() ) {
            // Nodes from first to the `marked` are already added.
            trim_front( frontier, marked );
            marked = 0;
            extend_to_k( frontier, iter, k );
          }
          if ( patch.empty() ) {
            // Search for a patch of length k that is not covered by visited paths of `iter`.
            while ( iter[ frontier.get_nodes() ] ) {
              add_node( frontier, *iter );
              ltrim_front_by_len( frontier, k, true );
              ++iter;
            }
          }
          // Extend the patch.
          patch += frontier;
          rtrim_front_by_len( frontier, k );
          while ( !iter[ frontier.get_nodes() ] ) {
            add_node( frontier, *iter );
            add_node( patch, *iter );
            rtrim_front_by_len( frontier, k );
            ++iter;
          }
        }
      }
      catch( const std::range_error& ) {
        if ( length( patch ) > 0 ) {
          if ( !iter[ frontier.get_nodes() ] &&
              !rcontains( patch, frontier.get_nodes().rbegin(), frontier.get_nodes().rend() ) )
          {
            if ( marked != 0 ) trim_front( frontier, marked );
            patch += frontier;
          }
          paths.push_back( std::move( patch ) );
        }
        --iter;  // save the traversed path and reset the Haplotyper iterator.
      }
      iter.raise_on_end = false;
    }

  template< typename TPathSet, typename TSpec >
      inline bool
    get_uniq_patched_haplotype( TPathSet& paths,
        GraphIter< VarGraph, Haplotyper< TSpec > >& iter,
        unsigned int context_len )
    {
      assert( context_len != 0 );
      if ( level( iter ) == 0 ) {
        get_uniq_full_haplotype( paths, iter );
        return true;
      }
      unsigned int paths_no = paths.size();
      get_uniq_patches( paths, iter, context_len );
      if ( paths_no == paths.size() ) return false;
      return true;
    }

  template< typename TPath >
      inline void
    get_rnd_full_haplotype( TPath& haplotype,
        GraphIter< VarGraph, Haplotyper< Random > >& iter )
    {
      while ( !at_end( iter ) ) {
        haplotype.push_back( *iter );
        ++iter;
      }
      --iter;
    }

  template< typename TGraph >
    inline std::size_t
    count_kmers( TGraph const& vargraph, unsigned int k, bool forward=false )
    {
      typedef typename TGraph::nodeid_type rank_type;
      typedef typename TGraph::rank_type nodeid_type;
      typedef typename seqan::Iterator< TGraph, Backtracker >::Type bter_type;

      if ( !k ) return 0;

      if ( !forward ) {
        throw std::runtime_error( "Counting k-mers on both strands is not implemented" );
      }

      bter_type bt_itr( vargraph );
      Path< TGraph > trav_path( &vargraph );

      std::size_t counter = 0;
      for ( rank_type rank = 1; rank <= vargraph.max_node_rank(); ++rank ) {
        nodeid_type id = vargraph.rank_to_id( rank );
        auto label_len = vargraph.node_length( id );

        // Count k-mers in the node
        std::size_t precontext = k - 1;  // precontext = max( label_len, k - 1 )
        if ( label_len >= k ) {
          counter += label_len - k + 1;
          precontext = label_len;
        }

        // Count k-mers across the node
        go_begin( bt_itr, id );
        while ( !at_end( bt_itr ) ) {
          extend_to_k( trav_path, bt_itr, label_len - 1 + k );
          if ( trav_path.get_sequence_len() >= k ) {
            // Number of k-mers across the node: min( |trav_path| - precontext, k - 1 )
            counter += std::min( trav_path.get_sequence_len() - precontext,
                static_cast< std::size_t >( k - 1 ) );
          }
          --bt_itr;
          trim_back( trav_path, *bt_itr );
        }
      }
      return counter;
    }

  /* END OF Haplotyper iterator interface functions  ----------------------------- */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef VARGRAPH_H__  ----- */

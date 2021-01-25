/**
 *    @file  graph_iter.hpp
 *   @brief  Generic graph iterator template classes.
 *
 *  This is a header file for generic graph iterators. All definitions/declarations are
 *  generic and they should be specialized for a specific graph type.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Fri Mar 03, 2017  11:28
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  PSI_GRAPH_ITER_HPP__
#define  PSI_GRAPH_ITER_HPP__

#include <map>
#include <deque>
#include <vector>
#include <utility>
#include <unordered_set>
#include <stdexcept>

#include "path.hpp"


namespace psi {
  /* GraphIter strategies. */
  struct BFS { };
  struct DFS { };
  struct Backtracker { };

  struct Global { };
  struct Local { };
  struct Random { };
  template< typename TSpec = Global >
  struct Haplotyper { };

  /* Graph iterator traits. */
  template< class TGraph, typename TSpec >
  struct GraphIterTraits;

  /**
   *  @brief  Breadth-first search graph iterator trait.
   */
  template< class TGraph >
  struct GraphIterTraits< TGraph, BFS > {
    typedef TGraph graph_type;
    typedef typename graph_type::id_type value_type;
    typedef value_type level_type;
    typedef std::deque< std::pair< value_type, level_type > > container_type;
    typedef std::unordered_set< value_type > set_type;
    typedef struct {
      std::size_t lb_rank;  /**< @brief lower-bound for rank of visited nodes. */
    } state_type;
    typedef void* param_type;

    constexpr static const param_type param_default = nullptr;
  };  /* --- end of struct GraphIterTraits (BFS specialisation) --- */

  /**
   *  @brief  Backtracker graph iterator trait.
   */
  template< class TGraph >
  struct GraphIterTraits< TGraph, Backtracker > {
    typedef TGraph graph_type;
    typedef typename graph_type::id_type value_type;
    typedef value_type level_type;
    typedef std::deque< std::pair< value_type, value_type > > container_type;
    typedef void* set_type;
    typedef struct {
      value_type start;   /**< @brief Start node ID. */
      value_type buffer;  /**< @brief Buffer node ID. 0=nothing */
    } state_type;
    typedef void* param_type;

    constexpr static const param_type param_default = nullptr;
  };  /* --- end of struct GraphIterTraits (Backtracker specialisation) --- */

  /**
   *  @brief  Haplotyper graph iterator trait.
   */
  template< class TGraph, typename TSpec >
  struct GraphIterTraits< TGraph, Haplotyper< TSpec > > {
    typedef TGraph graph_type;
    typedef typename graph_type::id_type value_type;
    typedef Path< graph_type, Dynamic > path_type;
    typedef std::unique_ptr< path_type > container_type;
    /**< @brief Set of visited paths. */
    typedef Path< graph_type, Haplotype > happath_type;
    typedef std::vector< happath_type > set_type;
    typedef typename set_type::size_type level_type;  /**< @brief No. of selected path so far. */
    typedef struct {
      value_type start;  /**< @brief Start node ID. */
      std::unique_ptr< happath_type > current_path;
      unsigned int setback;
      unsigned int entropy;
    } state_type;
    typedef unsigned int param_type;

    static const param_type param_default = 0;
  };  /* --- end of struct GraphIterTraits (Haplotyper general specialisation) --- */

  /**
   *  @brief  Haplotyper graph iterator trait (Random specialisation).
   */
  template< class TGraph >
  struct GraphIterTraits< TGraph, Haplotyper< Random > > {
    typedef TGraph graph_type;
    typedef typename graph_type::id_type value_type;
    typedef value_type level_type;
    typedef void* container_type;
    typedef void* set_type;
    typedef struct {
      value_type start;  /**< @brief Start node ID. */
      level_type level;  /**< @brief Level. */
    } state_type;
    typedef unsigned int param_type;

    constexpr static const param_type param_default = 0;
  };  /* --- end of struct GraphIterTraits (Random Haplotyper specialisation) --- */

  /**
   *  @brief  Graph iterator template base class.
   *
   *  @tparam  TGraph The graph type.
   *  @tparam  TSpec  Specify the traversal algorithm.
   *
   *  Graph iterator template class for traversing a sequence graph. The iterator will
   *  go forward by operator++ and backward by operator-- (if supported by iterator
   *  type). They can be initialised and check if it reached to the end by `begin` and
   *  `end` interface functions, respectively.
   */
  template< class TGraph, typename TSpec >
  class GraphIterBase {
  public:
    /* === TYPEDEFS === */
    typedef TGraph graph_type;
    typedef TSpec spec_type;
    typedef GraphIterTraits< graph_type, spec_type > traits_type;
    typedef typename traits_type::value_type value_type;
    typedef typename traits_type::level_type level_type;
    typedef typename traits_type::container_type container_type;
    typedef typename traits_type::set_type set_type;
    typedef typename traits_type::state_type state_type;
    typedef typename traits_type::param_type param_type;

    /* === CONSTANTS === */
    constexpr static const int end_value = 0;

    /* === TYPEDEFS === */
    struct EndIterator {
      constexpr static const value_type value = GraphIterBase::end_value;
    };
    typedef EndIterator end_type;

    /* === LIFECYCLE === */
    GraphIterBase( graph_type const* ptr=nullptr )
      : graph_ptr( ptr ), value( GraphIterBase::end_value ), raise_on_end( false )
    { }

    GraphIterBase( GraphIterBase const& ) = default;
    GraphIterBase( GraphIterBase&& ) = default;
    ~GraphIterBase( ) = default;

    /* === ACCESSORS === */
    inline graph_type const*
    get_graph_ptr( ) const
    {
      return this->graph_ptr;
    }

    static constexpr inline param_type
    get_default_param( )
    {
      return traits_type::param_default;
    }

    inline bool
    get_raise_on_end( ) const
    {
      return this->raise_on_end;
    }

    /* === ACCESSORS === */
    inline void
    set_raise_on_end( bool value )
    {
      this->raise_on_end = value;
    }

    /* === OPERATORS === */
    GraphIterBase& operator=( GraphIterBase const& ) = default;
    GraphIterBase& operator=( GraphIterBase&& ) = default;

    inline value_type
    operator*( ) const
    {
      return this->value;
    }

  protected:
    /* === DATA MEMBERS === */
    graph_type const* graph_ptr;  /**< @brief Pointer to graph. */
    value_type value;             /**< @brief Iter. current value. */
    container_type visiting;      /**< @brief Visiting buffer. */
    set_type visited;             /**< @brief Visited set. */
    state_type state;             /**< @brief Special-purpose vars. */
    param_type param;             /**< @brief Input parameter. */
    bool raise_on_end;            /**< @brief Throw an exception if it hits the end. */
  };  /* --- end of template class GraphIterBase --- */

  template< class TGraph, typename TSpec >
  class GraphIter;

  template< class TGraph, typename TSpec >
  inline GraphIter< TGraph, TSpec >
  begin( TGraph const& graph,
         TSpec /* tag */,
         typename GraphIter< TGraph, TSpec >::value_type start=0,
         typename GraphIter< TGraph, TSpec >::param_type p=
           GraphIter< TGraph, TSpec >::get_default_param() )
  {
    return GraphIter< TGraph, TSpec >( graph, start, p );
  }

  template< class TGraph, typename TSpec >
  inline typename GraphIter< TGraph, TSpec >::end_type
  end( TGraph const& graph, TSpec /* tag */ )
  {
    return typename GraphIter< TGraph, TSpec >::end_type();
  }

  template< class TGraph, typename TSpec >
  inline bool
  operator==( GraphIterBase< TGraph, TSpec > const& it,
              typename GraphIterBase< TGraph, TSpec >::end_type& end )
  {
    typedef typename GraphIterBase< TGraph, TSpec >::end_type end_type;
    return *it == end_type::value;
  }

  template< class TGraph, typename TSpec >
  inline bool
  operator==( typename GraphIterBase< TGraph, TSpec >::end_type& end,
              GraphIterBase< TGraph, TSpec > const& it )
  {
    return it == end;
  }

  template< class TGraph, typename TSpec >
  inline bool
  operator!=( GraphIterBase< TGraph, TSpec > const& it,
              typename GraphIterBase< TGraph, TSpec >::end_type& end )
  {
    return !( it == end );
  }

  template< class TGraph, typename TSpec >
  inline bool
  operator!=( typename GraphIterBase< TGraph, TSpec >::end_type& end,
              GraphIterBase< TGraph, TSpec > const& it )
  {
    return !( it == end );
  }

  template< class TGraph >
  class GraphIter< TGraph, BFS >
    : public GraphIterBase< TGraph, BFS >
  {
  public:
    /* === TYPEDEFS === */
    typedef TGraph graph_type;
    typedef BFS spec_type;
    typedef GraphIterBase< graph_type, spec_type > base_type;
    using typename base_type::value_type;
    using typename base_type::level_type;
    using typename base_type::container_type;
    using typename base_type::set_type;
    using typename base_type::state_type;
    using typename base_type::param_type;

    /* === LIFECYCLE === */
    GraphIter( graph_type const& g,
               value_type start=0,
               param_type=GraphIter::get_default_param() )
      : base_type( &g )
    {
      if ( start == 0 ) start = g.rank_to_id( 1 );

      this->state.lb_rank = 1;
      if ( g.id_to_rank( start ) == 1 ) {
        ++this->state.lb_rank;
      }

      this->value = start;
      this->visiting.push_back( { this->value, 0 } );
      this->visited.insert( this->value );
    }

    GraphIter( ) : base_type( ) { }

    GraphIter( GraphIter const& ) = default;
    GraphIter( GraphIter&& ) = default;
    ~GraphIter( ) = default;

    /* === OPERATORS === */
    GraphIter& operator=( GraphIter const& ) = default;
    GraphIter& operator=( GraphIter&& ) = default;

    inline GraphIter&
    operator++( )
    {
      typedef typename graph_type::id_type id_type;
      typedef typename graph_type::linktype_type linktype_type;

      if ( this->visiting.empty() ) {
        assert( this->value == base_type::end_value );
        return *this;
      }

      level_type plevel = this->level();
      this->graph_ptr->for_each_edges_out(
          this->value,
          [this, &plevel]( id_type to, linktype_type type ) {
            if ( !( *this )[ to ] ) {
              this->visiting.push_back( { to, plevel + 1 } );
              this->visited.insert( to );
            }
            return true;
          } );
      this->visiting.pop_front();

      if ( !this->visiting.empty() ) {
        this->value = this->visiting.front().first;
      }
      else {
        this->value = this->next_unvisited();
        if ( this->value != base_type::end_value ) {
          this->visiting.push_back( { this->value, 0 } );
          this->visited.insert( this->value );
        }
        else if ( this->raise_on_end ) {
          throw std::range_error( "end of iteration" );
        }
      }

      if ( this->value != base_type::end_value &&
          this->state.lb_rank == this->graph_ptr->id_to_rank( this->value ) ) {
        ++this->state.lb_rank;
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
    inline bool
    operator[]( value_type id )
    {
      return this->visited.find( id ) != this->visited.end();
    }

    /* === METHODS === */
    inline void
    reset( value_type start=0 )
    {
      if ( start == 0 ) start = this->graph_ptr->rank_to_id( 1 );

      this->state.lb_rank = 1;
      if ( this->graph_ptr->id_to_rank( start ) == 1 ) {
        ++this->state.lb_rank;
      }

      this->value = start;
      this->visiting.clear();
      this->visiting.push_back( { this->value, 0 } );
      this->visited.clear();
      this->visited.insert( this->value );
    }

    inline level_type
    level( ) const
    {
      if ( !this->visiting.empty() ) return this->visiting.front().second;
      throw std::runtime_error( "invalid level query on the end of iterator" );
    }

  private:
    /* === METHODS === */
    /**
     *  @brief  Search the graph for next unvisited node.
     *
     *  @return the node ID of the next visited node or `0` if all are visited.
     *
     *  The lower-bound for visited rank is also updated so that any node with smaller
     *  rank is visited. So in order to find the next unvisited node, we should search
     *  among nodes with higher ranks.
     */
    inline value_type
    next_unvisited( )
    {
      typedef typename graph_type::id_type id_type;
      typedef typename graph_type::rank_type rank_type;

      value_type unvisited = base_type::end_value;
      this->graph_ptr->for_each_node(
          [this, &unvisited]( rank_type rank, id_type id ) {
            if ( !( *this )[ id ] ) {
              unvisited = id;
              this->state.lb_rank = rank;
              return false;
            }
            return true;
          },
          this->state.lb_rank
        );
      return unvisited;
    }
  };  /* --- end of template class GraphIter (BFS specialisation) --- */

  template< class TGraph >
  class GraphIter< TGraph, Backtracker >
    : public GraphIterBase< TGraph, Backtracker >
  {
  public:
    /* === TYPEDEFS === */
    typedef TGraph graph_type;
    typedef Backtracker spec_type;
    typedef GraphIterBase< TGraph, Backtracker > base_type;
    typedef typename graph_type::id_type id_type;
    typedef typename graph_type::rank_type rank_type;
    typedef typename graph_type::linktype_type linktype_type;
    using typename base_type::value_type;
    using typename base_type::level_type;
    using typename base_type::container_type;
    using typename base_type::set_type;
    using typename base_type::state_type;
    using typename base_type::param_type;

    /* === LIFECYCLE === */
    GraphIter( graph_type const& g,
               value_type start=0,
               param_type=GraphIter::get_default_param() )
      : base_type( &g )
    {
      if ( start == 0 ) start = g.rank_to_id( 1 );

      this->value = start;
      this->state.buffer = base_type::end_value;
      this->state.start = start;
    }

    GraphIter( ) : base_type( ) { }

    GraphIter( GraphIter const& ) = default;
    GraphIter( GraphIter&& ) = default;
    ~GraphIter( ) = default;

    /* === OPERATORS === */
    GraphIter& operator=( GraphIter const& ) = default;
    GraphIter& operator=( GraphIter&& ) = default;

    inline GraphIter&
    operator++( )
    {
      if ( this->state.buffer != base_type::end_value ) {  // Any node buffered?
        this->value = this->state.buffer;                  // Use it.
        this->state.buffer = base_type::end_value;         // Clear up the buffer.
      }
      else {                               // else
        value_type cnode = this->value;
        this->value = base_type::end_value;
        this->graph_ptr->for_each_edges_out(
            cnode,
            [this, &cnode]( id_type to, linktype_type type ) {
              if ( this->value == base_type::end_value ) {
                this->value = to;
                return true;
              }
              this->visiting.push_back( { cnode, to } );
              return true;
            } );
        if ( this->value == base_type::end_value && this->raise_on_end ) {
          throw std::range_error( "end of iteration" );
        }
      }
      return *this;
    }

    inline GraphIter&
    operator--( )
    {
      if ( this->state.buffer != base_type::end_value ) {  // Any node buffered?
        // Remove all buffered branches of the current node.
        while ( !this->visiting.empty() && this->visiting.back().first == this->value ) {
          this->visiting.pop_back();
        }
        this->state.buffer = base_type::end_value;
      }

      this->value = base_type::end_value;

      if ( !this->visiting.empty() ) {                 // Go back in buffer.
        this->value = this->visiting.back().first;
        this->state.buffer = this->visiting.back().second;
        this->visiting.pop_back();
      }

      return *this;
    }

    /* === METHODS === */
    inline void
    reset( value_type start=0 )
    {
      if ( start == 0 ) start = this->state.start;  // Re-use stored start node.

      this->value = start;
      this->state.buffer = base_type::end_value;  // Re-set buffer.
      this->visiting.clear();
    }
  };  /* --- end of template class GraphIter (Backtracker specialisation) --- */

  template< class TGraph >
  class GraphIter< TGraph, Haplotyper<> >
    : public GraphIterBase< TGraph, Haplotyper<> >
  {
  public:
    /* === TYPEDEFS === */
    typedef TGraph graph_type;
    typedef Haplotyper<> spec_type;
    typedef GraphIterBase< graph_type, spec_type > base_type;
    typedef typename graph_type::id_type id_type;
    typedef typename graph_type::rank_type rank_type;
    typedef typename graph_type::linktype_type linktype_type;
    using typename base_type::value_type;
    using typename base_type::level_type;
    using typename base_type::container_type;
    typedef typename container_type::element_type path_type;
    using typename base_type::set_type;
    typedef typename set_type::value_type happath_type;
    using typename base_type::state_type;
    using typename base_type::param_type;

    /* === LIFECYCLE === */
    GraphIter( graph_type const& g,
               value_type start=0,
               param_type=GraphIter::get_default_param() )
      : base_type( &g )
    {
      if ( start == 0 ) start = g.rank_to_id( 1 );

      this->value = start;
      this->state.start = start;
      this->visiting = std::make_unique< path_type >( &g );
      this->state.current_path = std::make_unique< happath_type >( &g );
      this->state.current_path->push_back( this->value );
      this->state.setback = 0;
      this->state.entropy = 1;
    }

    GraphIter( ) : base_type( ) { }

    GraphIter( GraphIter const& ) = default;
    GraphIter( GraphIter&& ) = default;
    ~GraphIter( ) = default;

    /* === OPERATORS === */
    GraphIter& operator=( GraphIter const& ) = default;
    GraphIter& operator=( GraphIter&& ) = default;

    /**
     *  A setback path is a sequence of last 's' nodes of currently generating haplotype.
     *  An unvisited setback path is a path that does not occur as a subset of any
     *  previously generated haplotypes. This function search in adjacent nodes set for
     *  a node that together with 's-1' previously selected nodes forms an unvisited
     *  setback in order to cover more k-mers from all paths in the graph; i.e. generating
     *  more diverse haplotypes.
     */
    inline GraphIter&
    operator++( )
    {
      if ( !this->graph_ptr->has_edges_out( this->value ) ) {
        this->value = base_type::end_value;
        if ( this->raise_on_end ) throw std::range_error( "end of iteration" );
        return *this;
      }

      if ( this->state.setback > 1 ) {
        while ( this->visiting->size() != 0 &&
                this->state.entropy > this->state.setback ) {
          this->state.entropy /= this->graph_ptr->outdegree( this->visiting->front() );
          this->visiting->pop_front();
        }
      }

      value_type candidate = base_type::end_value;
      if ( this->state.setback == 0 || this->graph_ptr->outdegree( this->value ) == 1 ) {
        this->graph_ptr->for_each_edges_out(
            this->value,
            [&candidate]( id_type id, linktype_type ) {
              candidate = id;
              return false;
            } );
      }
      else {
        // Search for a forward node such that the setback path is not in previous paths.
        do {
          this->graph_ptr->for_each_edges_out(
              this->value,
              [this, &candidate]( id_type to, linktype_type ) {
                this->visiting->push_back( to );
                if ( ( *this )[ *this->visiting ] ) {  // Visited?
                  this->visiting->pop_back();
                  return true;
                }
                this->visiting->pop_back();            // No change to the iterator state.
                candidate = to;                        // Found!
                return false;
              } );
        } while ( this->state.setback == 1 &&
                  candidate == base_type::end_value &&
                  this->visiting->empty() &&
                  ( this->visiting->push_back( this->value ), true ) );
        if ( this->state.setback == 1 && !this->visiting->empty() ) this->visiting->pop_back();
      }
      // If no unvisited setback path found, use a node with least path coverage.
      if ( candidate == base_type::end_value ) {
        candidate = util::least_covered_adjacent( *this->graph_ptr,
                                                  *this->visiting,
                                                  this->visited );
      }
      // If all forward edges are visited, pick one randomly with uniform distribution.
      if ( candidate == base_type::end_value ) {
        candidate = util::random_adjacent( *this->graph_ptr, this->value );
      }

      this->value = candidate;
      if ( this->state.setback > 1 ) {
        this->visiting->push_back( this->value );
        this->state.entropy *= this->graph_ptr->outdegree( this->value );
      }
      this->state.current_path->push_back( this->value );
      return *this;
    }

    inline GraphIter&
    operator--( int )
    {
      this->value = this->state.start;  // Reset the iterator to the start node.
      this->visiting->clear();
      this->state.entropy = 1;
      if ( this->state.setback > 1 ) {
        this->visiting->push_back( this->value );
        this->state.entropy *= this->graph_ptr->outdegree( this->value );
      }
      this->state.current_path->clear();
      this->state.current_path->push_back( this->value );
      return *this;
    }

    inline GraphIter&
    operator--( )
    {
      this->visited.push_back( std::move( *this->state.current_path ) );
      this->set_setback();
      ( *this )--;
      return *this;
    }

    /**
     *  @brief  Check if the given path is present in paths generated so far.
     *
     *  @param  path A container of node IDs indicating nodes in a path.
     *  @return `true` if the path is present; `false` otherwise.
     *
     *  Check whether the given path is generated before or not.
     */
    template< class TContainer >
    inline bool
    operator[]( TContainer const& path )
    {
      // :TODO:Sun Apr 01 00:01:\@cartoonist: Unordered check!
      // ... `PathSet< Path< TGraph, Compact >, InMemory >` should be used for
      // ... `this->visited`.
      return covered_by( path.begin(), path.end(), this->visited );
    }

    /* === METHODS === */
    inline void
    reset( value_type start=0 )
    {
      if ( start == 0 ) start = this->state.start;  // Re-use start node.

      this->value = start;
      this->state.start = start;
      this->visiting->clear();
      this->visited.clear();
      this->state.current_path->clear();
      this->state.current_path->push_back( this->value );
      this->state.setback = 0;
      this->state.entropy = 1;
    }

    inline level_type
    level( ) const
    {
      return this->visited.size();
    }

  private:
    /* === METHODS === */
    inline void
    set_setback( )
    {
      this->state.setback = this->visited.size();
    }
  };  /* --- end of template class GraphIter (Haplotyper<> specialisation) --- */

  template< class TGraph >
  class GraphIter< TGraph, Haplotyper< Local > >
    : public GraphIterBase< TGraph, Haplotyper< Local > >
  {
  public:
    /* === TYPEDEFS === */
    typedef TGraph graph_type;
    typedef Haplotyper< Local > spec_type;
    typedef GraphIterBase< graph_type, spec_type > base_type;
    typedef typename graph_type::id_type id_type;
    typedef typename graph_type::rank_type rank_type;
    typedef typename graph_type::linktype_type linktype_type;
    using typename base_type::value_type;
    using typename base_type::level_type;
    using typename base_type::container_type;
    typedef typename container_type::element_type path_type;
    using typename base_type::set_type;
    typedef typename set_type::value_type happath_type;
    using typename base_type::state_type;
    using typename base_type::param_type;

    /* === LIFECYCLE === */
    GraphIter( graph_type const& g,
               value_type start=0,
               param_type p=GraphIter::get_default_param() )
      : base_type( &g )
    {
      if ( start == 0 ) start = g.rank_to_id( 1 );
      if ( p == 0 ) {
        throw std::runtime_error( "Parameter value of Local Haplotyper cannot be zero" );
      }

      this->value = start;
      this->visiting = std::make_unique< path_type >( &g );
      this->state.start = start;
      this->state.current_path = std::make_unique< happath_type >( &g );
      this->state.current_path->push_back( this->value );
      this->state.setback = 0;
      this->param = p;
    }

    GraphIter( ) : base_type( ) { }

    GraphIter( GraphIter const& ) = default;
    GraphIter( GraphIter&& ) = default;
    ~GraphIter( ) = default;

    /* === OPERATORS === */
    GraphIter& operator=( GraphIter const& ) = default;
    GraphIter& operator=( GraphIter&& ) = default;

    /**
     *  A setback path is a sequence of last 's' nodes of currently generating haplotype
     *  encompassing last k basepair. An unvisited setback path is a path that does
     *  not occur as a subset of any previously generated haplotypes. This function search
     *  in adjacent nodes set for a node that together with 's-1' previously selected nodes
     *  forms an unvisited setback in order to cover more k-mers from all paths in the graph;
     *  i.e. generating more diverse haplotypes.
     */
    inline GraphIter&
    operator++( )
    {
      if ( !this->graph_ptr->has_edges_out( this->value ) ) {
        this->value = base_type::end_value;
        if ( this->raise_on_end ) throw std::range_error( "end of iteration" );
        return *this;
      }

      if ( this->state.setback != 0 ) {
        rtrim_front_by_len( *this->visiting, this->param - 1 );
      }

      value_type candidate = base_type::end_value;
      if ( this->state.setback == 0 || this->graph_ptr->outdegree( this->value ) == 1 ) {
        this->graph_ptr->for_each_edges_out(
            this->value,
            [&candidate]( id_type id, linktype_type ) {
              candidate = id;
              return false;
            } );
      }
      else {
        // Search for a forward node such that the setback path is not in previous paths.
        this->graph_ptr->for_each_edges_out(
            this->value,
            [this, &candidate]( id_type to, linktype_type ) {
              this->visiting->push_back( to );
              if ( ( *this )[ *this->visiting ] ) {  // Visited?
                this->visiting->pop_back();
                return true;
              }
              this->visiting->pop_back();            // No change to the iterator state.
              candidate = to;                        // Found!
              return false;
            } );
      }
      // If no unvisited setback path found, use a node with least path coverage.
      if ( candidate == base_type::end_value ) {
        candidate = util::least_covered_adjacent( *this->graph_ptr,
                                                  *this->visiting,
                                                  this->visited );
      }
      // If all forward edges are visited, pick one randomly with uniform distribution.
      if ( candidate == base_type::end_value ) {
        candidate = util::random_adjacent( *this->graph_ptr,  this->value );
      }

      this->value = candidate;
      if ( this->state.setback != 0 ) {
        this->visiting->push_back( this->value );
      }
      this->state.current_path->push_back( this->value );
      return *this;
    }

    inline GraphIter&
    operator--( int )
    {
      this->value = this->state.start;    // Reset the iterator to the start node.
      this->visiting->clear();
      if ( this->state.setback != 0 ) {
        this->visiting->push_back( this->value );
      }
      this->state.current_path->clear();
      this->state.current_path->push_back( this->value );
      return *this;
    }

    inline GraphIter&
    operator--( )
    {
      this->visited.push_back( std::move( *this->state.current_path ) );
      this->set_setback();
      ( *this )--;
      return *this;
    }

    /**
     *  @brief  Check if the given path is present in paths generated so far.
     *
     *  @param  path A container of node IDs indicating nodes in a path.
     *  @return `true` if the path is present; `false` otherwise.
     *
     *  Check whether the given path is generated before or not.
     */
    template< class TContainer >
    inline bool
    operator[]( TContainer const& path )
    {
      // :TODO:Sun Apr 01 00:01:\@cartoonist: Unordered check!
      // ... `PathSet< Path< TGraph, Compact >, InMemory >` should be used for
      // ... `this->visited`.
      return covered_by( path.begin(), path.end(), this->visited );
    }

    /* === METHODS === */
    inline void
    reset( value_type start=0,
           param_type p=GraphIter::get_default_param() )
    {
      if ( start == 0 ) start = this->state.start;  // Re-use start node.
      if ( p == 0 ) p = this->param;                // Re-use parameter value.

      this->value = start;
      this->state.start = start;
      this->visiting->clear();
      this->visited.clear();
      this->state.current_path->clear();
      this->state.current_path->push_back( this->value );
      this->state.setback = 0;
      this->param = p;
    }

    inline level_type
    level( ) const
    {
      return this->visited.size();
    }

  private:
    /* === METHODS === */
    inline void
    set_setback( )
    {
      if ( this->visited.size() == 0 ) {
        this->state.setback = 0;
      }
      else {
        this->state.setback = 2 * ceil( log2( this->visited.size() + 1 ) ) - 1;
      }
    }
  };  /* --- end of template class GraphIter (Haplotyper< Local > specialisation) --- */

  template< class TGraph >
  class GraphIter< TGraph, Haplotyper< Random > >
    : public GraphIterBase< TGraph, Haplotyper< Random > >
  {
  public:
    /* === TYPEDEFS === */
    typedef TGraph graph_type;
    typedef Haplotyper< Random > spec_type;
    typedef GraphIterBase< graph_type, spec_type > base_type;
    using typename base_type::value_type;
    using typename base_type::level_type;
    using typename base_type::container_type;
    using typename base_type::set_type;
    using typename base_type::state_type;
    using typename base_type::param_type;

    /* === LIFECYCLE === */
    GraphIter( graph_type const& g,
               value_type start=0,
               param_type p=GraphIter::get_default_param() )
      : base_type( &g )
    {
      if ( start == 0 ) start = g.rank_to_id( 1 );

      this->value = start;
      this->state.start = start;
      this->state.level = 1;
      this->param = p;
    }

    GraphIter( ) : base_type( ) { }

    GraphIter( GraphIter const& ) = default;
    GraphIter( GraphIter&& ) = default;
    ~GraphIter( ) = default;

    /* === OPERATORS === */
    GraphIter& operator=( GraphIter const& ) = default;
    GraphIter& operator=( GraphIter&& ) = default;

    inline GraphIter&
    operator++( )
    {
      if ( !this->graph_ptr->has_edges_out( this->value ) ) {
        this->value = base_type::end_value;
        if ( this->raise_on_end ) throw std::range_error( "end of iteration" );
        return *this;
      }
      this->value = util::random_adjacent( *this->graph_ptr, this->value, this->param );
      ++this->state.level;
      return *this;
    }

    inline GraphIter&
    operator--( )
    {
      this->reset( );
      return *this;
    }

    /* === METHODS === */
    inline void
    reset( value_type start=0,
           param_type p=GraphIter::get_default_param() )
    {
      if ( start == 0 ) start = this->state.start;  // Re-use start node.
      if ( p == 0 ) p = this->param;                // Re-use parameter value.

      this->value = start;
      this->state.start = start;
      this->state.level = 1;
      this->param = p;
    }

    inline level_type
    level( ) const
    {
      return this->state.level;
    }
  };  /* --- end of template class GraphIter (Haplotyper< Random > specialisation) --- */

  namespace util {
    /**
     *  @brief  Extend a path to length k using a graph iterator.
     *
     *  @param  path The extended path.
     *  @param  iter The graph iterator.
     *  @param  k The extension length.
     *
     *  Add nodes using graph iterator until the length of path reach k or longer.
     */
    template< class TGraph, typename TPathSpec, typename TIterSpec >
    inline void
    extend_to_k( Path< TGraph, TPathSpec >& path,
                 GraphIter< TGraph, TIterSpec >& iter,
                 typename GraphIter< TGraph, TIterSpec >::end_type& end,
                 unsigned int k )
    {
      while ( iter != end && path.get_sequence_len() < k ) {
        add_node( path, *iter );
        ++iter;
      }
    }

    /**
     *  @overload Prevent using BFS for extension.
     */
    template< class TGraph, typename TPathSpec >
    inline void
    extend_to_k( Path< TGraph, TPathSpec >&,
                 GraphIter< TGraph, BFS >&,
                 typename GraphIter< TGraph, BFS >::end_type&,
                 unsigned int )
    {
      throw std::runtime_error( "a path cannot be extended by a BFS iterator" );
    }

    template< class TGraph >
    inline std::size_t
    count_kmers( TGraph const& graph, unsigned int k, bool forward=false )
    {
      typedef typename TGraph::id_type id_type;
      typedef typename TGraph::rank_type rank_type;

      if ( !k ) return 0;

      if ( !forward ) {
        throw std::runtime_error( "counting k-mers on both strands is not implemented" );
      }

      auto bt_itr = begin( graph, Backtracker() );
      auto bt_end = end( graph, Backtracker() );
      Path< TGraph > trav_path( &graph );

      std::size_t counter = 0;
      graph.for_each_node(
          [&]( rank_type rank, id_type id ) {
            auto label_len = graph.node_length( id );

            // Count k-mers in the node
            std::size_t precontext = k - 1;  // precontext = max( label_len, k - 1 )
            if ( label_len >= k ) {
              counter += label_len - k + 1;
              precontext = label_len;
            }

            // Count k-mers across the node
            bt_itr.reset( id );
            while ( bt_itr != bt_end ) {
              extend_to_k( trav_path, bt_itr, bt_end, label_len - 1 + k );
              if ( trav_path.get_sequence_len() >= k ) {
                // Number of k-mers across the node: min( |trav_path| - precontext, k - 1 )
                counter += std::min( trav_path.get_sequence_len() - precontext,
                                     static_cast< std::size_t >( k - 1 ) );
              }
              --bt_itr;
              trim_back( trav_path, *bt_itr );
            }
          } );
      return counter;
    }
  }  /* --- end of namespace util --- */
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_GRAPH_ITER_HPP__  --- */

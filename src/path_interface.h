/**
 *    @file  path_interface.h
 *   @brief  Interface function definitions for Path template class.
 *
 *  This header file contains interface functions for Path template class.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Fri Mar 30, 2018  11:58
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  PATH_INTERFACE_H__
#define  PATH_INTERFACE_H__

#include <functional>

#include "path_base.h"


namespace grem {
  /* Default Path interface functions  ----------------------------------------- */

  /**
   *  @brief  Initialize the internal data structure of the path.
   *
   *  @param  path The path.
   *
   *  See Path `initialize` member function.
   */
  template< typename TGraph, typename TSpec >
      inline void
    initialize( Path< TGraph, TSpec >& path )
    {
      path.initialize();
    }

  /**
   *  @brief  Extend the path forward.
   *
   *  @param  path The path.
   *  @param  node_id The new node ID.
   *
   *  Internal function to add the `node_id` to the end of the current path.
   */
  template< typename TGraph, typename TSpec >
      inline void
    add_node( Path< TGraph, TSpec >& path, typename TGraph::nodeid_type const& node_id )
    {
      path.push_back( node_id );
    }  /* -----  end of template function add_node  ----- */

  /**
   *  @brief  Extend a path with another one (interface function).
   *
   *  @param  path The path to be extended.
   *  @param  other The other path to use for extension.
   *
   *  @note It does not prevent self-extension. Use `operator+=` if this functionality
   *        is desired.
   */
  template< typename TGraph, typename TSpec1, typename TSpec2 >
      inline void
    extend( Path< TGraph, TSpec1 >& path, const Path< TGraph, TSpec2 >& other )
    {
      for ( const auto& node_id : other.get_nodes() ) {
        add_node( path, node_id );
      }
    }

  /**
   *  @brief  Extend the path with another path of the same type.
   *
   *  @param  path The path to be extended.
   *  @param  other The other path to use for extension.
   *  @return The reference to this instance after extension.
   */
  template< typename TGraph, typename TSpec >
      inline Path< TGraph, TSpec >&
    operator+=( Path< TGraph, TSpec >& path, const Path< TGraph, TSpec >& other )
    {
      if ( &path != &other ) {
        extend( path, other );
      }
      return path;
    }

  /**
   *  @brief  Extend the path with another path of different type.
   *
   *  @param  path The path to be extended.
   *  @param  other The other path to use for extension.
   *  @return The reference to this instance after extension.
   */
  template< typename TGraph, typename TSpec1, typename TSpec2 >
      inline Path< TGraph, TSpec1 >&
    operator+=( Path< TGraph, TSpec1 >& path, const Path< TGraph, TSpec2 >& other )
    {
      extend( path, other );
      return path;
    }

  template< typename TGraph, typename TSpec >
      inline typename Path< TGraph, TSpec >::size_type
    rank( const Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::seqsize_type pos )
    {
      return path.rank( pos );
    }

  template< typename TGraph, typename TSpec >
      inline typename Path< TGraph, TSpec >::seqsize_type
    select( const Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::size_type rank )
    {
      return path.select( rank );
    }

  /**
   *  @brief  Map a position in the path sequence to corresponding node ID.
   *
   *  @param  path The path.
   *  @param  pos The position on the path sequence.
   *  @return The node ID on whose label the given position lies.
   *
   *  It gets the node ID at rank `rank(pos)` in the node queue.
   */
  template< typename TGraph, typename TSpec >
      inline typename TGraph::nodeid_type
    position_to_id( const Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::seqsize_type pos )
    {
      return path.get_nodes()[ rank( path, pos ) ];
    }

  /**
   *  @brief  Map a position in the path sequence to corresponding node offset.
   *
   *  @param  path The path.
   *  @param  pos The position on the path sequence.
   *  @return The offset from the first base in the node on whose label the given
   *          position lies.
   */
  template< typename TGraph, typename TSpec >
      inline typename TGraph::offset_type
    position_to_offset( const Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::seqsize_type pos )
    {
      auto&& sel = select( path, rank( path, pos ) );
      assert( pos >= sel );
      return pos - sel;
    }

  /**
   *  @brief  Compute forward sequence from the nodes queue.
   *
   *  @param  path The path.
   *  @return Forward sequence represented by the path.
   *
   *  Compute the forward sequence of the path from nodes queue.
   */
  template< typename TGraph, typename TSpec >
      inline typename Path< TGraph, TSpec >::string_type
    sequence( const Path< TGraph, TSpec >& path, Forward, unsigned int context=0 )
    {
      assert( path.get_vargraph() != nullptr );

      if ( length( path ) == 0 ) return "";

      typename Path< TGraph, TSpec >::string_type repr_str;
      repr_str.reserve( path.get_sequence_len() );

      /* If context is set (not zero) the first node is trimmed by length `context`. */
      auto iter = path.get_nodes().begin();
      if ( context != 0 ){
        unsigned int off = ( path.get_vargraph()->node_length( *iter ) + 1 > context ?
            path.get_vargraph()->node_length( *iter ) - context + 1 :
            0 );
        repr_str += path.get_vargraph()->node_sequence( *iter )
          .substr( off );
        ++iter;
      }

      /* Add sequence of intermediate nodes. */
      while ( ( context != 0 && iter != path.get_nodes().end() - 1 ) ||
          ( context == 0 && iter != path.get_nodes().end() ) ) {
        repr_str += path.get_vargraph()->node_sequence( *iter );
        ++iter;
      }

      /* If context is set (not zero) the last node is trimmed by length `context`. */
      if ( context != 0 ) {
        repr_str += path.get_vargraph()->node_sequence( *iter ).substr( 0, context - 1 );
      }

      return repr_str;
    }  /* -----  end of template function sequence  ----- */

  /**
   *  @brief  Compute reversed sequence from the nodes queue.
   *
   *  @param  path The path.
   *  @return Reversed sequence represented by the path.
   *
   *  Compute the reversed sequence of the path from nodes queue.
   */
  template< typename TGraph, typename TSpec >
      inline typename Path< TGraph, TSpec >::string_type
    sequence( const Path< TGraph, TSpec >& path, Reversed, unsigned int context=0 )
    {
      typename Path< TGraph, TSpec >::string_type repr_str =
        sequence( path, Forward(), context );
      std::reverse( repr_str.begin(), repr_str.end() );
      return repr_str;
    }  /* -----  end of template function sequence  ----- */

  template< typename TGraph, typename TSpec >
      inline typename Path< TGraph, TSpec >::string_type
    sequence( const Path< TGraph, TSpec >& path, unsigned int context=0 )
    {
      return sequence( path, Forward(), context );
    }

  /**
   *  @brief  Clear the path.
   *
   *  @param  path The path.
   *
   *  Internal function to clear nodes queue and nodes set.
   */
  template< typename TGraph, typename TSpec >
      inline void
    clear( Path< TGraph, TSpec >& path )
    {
      path.clear();
    }  /* -----  end of template function clear  ----- */

  /**
   *  @brief  Reserve memory for the path.
   *
   *  @param  path The path.
   *  @param  size The target capacity.
   *
   *  Internal function to reserve memory for nodes and nodes set.
   */
  template< typename TGraph, typename TSpec >
      inline void
    reserve( Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::size_type size )
    {
      path.reserve( size );
    }  /* -----  end of template function reserve  ----- */

  /**
   *  @brief  Get path length (number of nodes).
   *
   *  @param  path The path.
   *  @return the number of nodes on the path.
   *
   *  It gets the length of the path which is the number of nodes in the path.
   */
  template< typename TGraph, typename TSpec >
      inline typename Path< TGraph, TSpec >::size_type
    length( const Path< TGraph, TSpec >& path )
    {
      return path.size();
    }  /* -----  end of template function length  ----- */

  /**
   * @brief  pop the last node from the path.
   */
  template< typename TGraph, typename TSpec >
      inline void
    pop_back( Path< TGraph, TSpec >& path )
    {
      path.pop_back();
    }  /* -----  end of method pop_back  ----- */

  /**
   * @brief  pop the first node from the path.
   *
   * NOTE: only implemented for Dynamic Path.
   */
  template< typename TGraph >
      inline void
    pop_front( Path< TGraph, Dynamic >& path )
    {
      path.pop_front();
    }

  /**
   *  @brief  Trim the path from the node whose ID matches the given ID from back.
   *
   *  @param  path The path to be trimmed.
   *  @param  node_id The node ID from where the path should be trimmed.
   *
   *  It pops the last nodes in the nodes list until it reaches to the given node ID
   *  inclusive. It will yield an empty path if `node_id` is not among the path nodes.
   *  If `node_id` is equal to zero, an invalid node ID, it pops only the last node from
   *  the path.
   */
  template< typename TGraph, typename TSpec >
      inline void
    trim_back( Path< TGraph, TSpec >& path, typename TGraph::nodeid_type node_id=0 )
    {
      bool found = false;
      while ( !found && length( path ) != 0 ) {
        auto&& last_node = path.get_nodes().back();
        if ( node_id == 0 || last_node == node_id ) found = true;
        pop_back( path );
      }
    }  /* -----  end of template function trim_back  ----- */

  /**
   *  @brief  Trim the path from the back until further trimming would lead length of < k.
   *
   *  @param  path The path to be trimmed.
   *  @param  k The length of k.
   *
   *  The nodes from back will be dropped from the path until to the point that further
   *  trimming leads its length to be less than k.
   */
  template< typename TGraph, typename TSpec >
      inline void
    trim_back_by_len( Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::seqsize_type k )
    {
      while ( length( path ) != 0 && path.get_sequence_len() -
          path.get_vargraph()->node_length( path.get_nodes().back() ) >= k ) {
        pop_back( path );
      }
    }  /* -----  end of template function trim_back_by_len  ----- */

  /**
   *  @brief  Trim the path from the node whose ID matches the given ID from front.
   *
   *  @param  path The path to be trimmed.
   *  @param  node_id The node ID from where the path should be trimmed.
   *
   *  It pops the first nodes in the nodes list until it reaches to the given node ID
   *  inclusive. It will yield an empty path if `node_id` is not among the path nodes.
   *  If `node_id` is equal to zero, an invalid node ID, it pops only the first node
   *  from the path.
   *
   *  @note specific for Dynamic Paths.
   */
  template< typename TGraph >
      inline void
    trim_front( Path< TGraph, Dynamic >& path, typename TGraph::nodeid_type node_id=0 )
    {
      bool found = false;
      while ( !found && length( path ) != 0 ) {
        auto&& first_node = path.get_nodes().front();
        if ( node_id == 0 || first_node == node_id ) found = true;
        pop_front( path );
      }
    }  /* -----  end of template function trim_front  ----- */

  /**
   *  @brief  Trim the path from the front until further trimming would lead length of < k.
   *
   *  @param  path The path to be trimmed.
   *  @param  k The length of k.
   *
   *  The nodes from front will be dropped from the path until to the point that further
   *  trimming leads its length to be less than k.
   */
  template< typename TGraph >
      inline void
    trim_front_by_len( Path< TGraph, Dynamic >& path,
        typename Path< TGraph, Dynamic >::seqsize_type k )
    {
      while ( length( path ) != 0 && path.get_sequence_len() -
          path.get_vargraph()->node_length( path.get_nodes().front() ) >= k ) {
        pop_front( path );
      }
    }  /* -----  end of template function trim_front_by_len  ----- */

  /* END OF Default Path interface functions  ---------------------------------- */

  /* Path interface functions  ------------------------------------------------- */

  /**
   *  @brief  Check whether the given node ID is on the path or not.
   *
   *  @param  path The path.
   *  @param  node_id The ID of the node.
   *  @return `true` if the node ID is on the path.
   *
   *  Check if the given node ID presents in the path.
   */
  template< typename TGraph, typename TSpec >
      inline bool
    contains( const Path< TGraph, TSpec >& path, typename TGraph::nodeid_type node_id )
    {
      return path.contains( node_id );
    }  /* -----  end of template function contains  ----- */

  /**
   *  @brief  Check whether the given node ID is on the path or not.
   *
   *  @param  path The path.
   *  @param  node_id The ID of the node.
   *  @return `true` if the node ID is on the path.
   *
   *  @overload for Haplotype path.
   */
  template< typename TGraph >
      inline bool
    contains( const Path< TGraph, Haplotype >& path, typename TGraph::nodeid_type node_id )
    {
      if ( node_id < 1 || node_id > path.nodes.size() ) return false;
      return path.nodes[ node_id - 1 ] == 1;
    }

  /**
   *  @brief  Check whether this path contains a set of node IDs.
   *
   *  @param  path The path.
   *  @param  begin The begin iterator of node IDs set to be checked.
   *  @param  end The end iterator of node IDs set to be checked.
   *  @return `true` if this path contains all node IDs in the range `[begin, end)`;
   *          otherwise `false` -- including the case that the range `[begin, end)` is
   *          empty; i.e. `begin == end`.
   *
   *  It checks if the nodes of the given path is present in the node set. It does not
   *  check the order of the nodes.
   */
  template< typename TGraph, typename TIter >
      inline bool
    contains( const Path< TGraph, Micro >& path, TIter begin, TIter end )
    {
      auto on_path = [&path]( typename TGraph::nodeid_type i ) {
        return contains( path, i );
      };

      if ( begin != end && std::all_of( begin, end, on_path ) ) return true;
      return false;
    }

  /**
   *  @brief  Check whether this path contains a set of node IDs.
   *
   *  @param  path The path.
   *  @param  begin The begin iterator of node IDs set to be checked.
   *  @param  end The end iterator of node IDs set to be checked.
   *  @return `true` if this path contains all node IDs in the range `[begin, end)`;
   *          otherwise `false` -- including the case that the range `[begin, end)` is
   *          empty; i.e. `begin == end`.
   *
   *  It checks if the nodes of the given path is present in the node set. It DOES check
   *  the order of the nodes.
   */
  template< typename TGraph, typename TIter >
      inline bool
    contains( const Path< TGraph, Haplotype >& path, TIter begin, TIter end )
    {
      if ( begin == end ) return false;

      typename TGraph::nodeid_type prev = 0;
      for ( ; begin != end; ++begin ) {
        if ( *begin <= prev || !contains( path, *begin ) ) return false;
        prev = *begin;
      }
      return true;
    }

  /**
   *  @brief  Check whether a shorter list of node IDs is a subsequent of bigger one.
   *
   *  @param  pn_begin The begin iterator of the node IDs in bigger list.
   *  @param  pn_end The end iterator of the node IDs in bigger list.
   *  @param  begin The begin iterator of the node IDs in smaller list.
   *  @param  end The end iterator of the node IDs in smaller list.
   *  @return `true` if the smaller list is a subsequence of the bigger one; otherwise
   *          `false` -- including the case that the lists are empty; i.e.
   *          `end == begin && pn_begin == pn_end`.
   *
   *  It checks if the nodes of the smaller list is present in the bigger one. It DOES
   *  check the order of the nodes and they should be present without gap.
   */
  template< typename TIter1, typename TIter2 >
      inline bool
    _contains( TIter1 pn_begin, TIter1 pn_end, TIter2 begin, TIter2 end )
    {
      if ( begin != end && pn_begin != pn_end && pn_end - pn_begin >= end - begin ) {
        auto lc = std::find( pn_begin, pn_end, *begin );
        if ( lc == pn_end || pn_end - lc < end - begin ) return false;
        if ( std::equal( begin, end, lc ) ) return true;
      }

      return false;
    }  /* -----  end of template function _contains  ----- */

  /**
   *  @brief  Check whether this path contains another path.
   *
   *  @param  path The path.
   *  @param  begin The begin iterator of node IDs set of the path to be checked.
   *  @param  end The end iterator of node IDs set of the path to be checked.
   *  @return `true` if this path is a superset of the given path; otherwise `false`
   *          -- including the case that the given path is empty; i.e.
   *          `end == begin`.
   *
   *  The input path should be smaller than the path. It checks if the nodes of the
   *  given path is present in the node set. It DOES check the order of the nodes.
   */
  template< typename TGraph, typename TSpec, typename TIter >
      inline bool
    contains( const Path< TGraph, TSpec >& path, TIter begin, TIter end )
    {
      return _contains( path.get_nodes().begin(), path.get_nodes().end(), begin, end );
    }  /* -----  end of template function contains  ----- */

  template< typename TGraph, typename TSpec, typename TIter >
      inline bool
    rcontains( const Path< TGraph, TSpec >& path, TIter rbegin, TIter rend )
    {
      return
        _contains( path.get_nodes().rbegin(), path.get_nodes().rend(), rbegin, rend );
    }  /* -----  end of template function rcontains  ----- */

  template< typename TGraph, typename TIter >
      inline bool
    rcontains( const Path< TGraph, Compact >& path, TIter rbegin, TIter rend )
    {
      const auto& nodes = path.get_nodes();
      std::size_t qlen = rend - rbegin;
      if ( rbegin != rend && nodes.size() >= qlen ) {
        auto lc = rfind( nodes, *rbegin );
        if ( lc == nodes.begin() || lc - nodes.begin() < rend - rbegin ) return false;
        if ( requal( rbegin, rend, lc, nodes.begin() ) ) return true;
      }

      return false;
    }  /* -----  end of template function rcontains  ----- */

  template< typename TGraph, typename TIter >
      inline bool
    rcontains( const Path< TGraph, Micro >& path, TIter rbegin, TIter rend )
    {
      return contains( path, rbegin, rend );
    }

  template< typename TGraph, typename TIter >
      inline bool
    rcontains( const Path< TGraph, Haplotype >& path, TIter rbegin, TIter rend )
    {
      if ( rbegin == rend ) return false;

      typename TGraph::nodeid_type prev = path.get_vargraph()->max_node_rank() + 1;
      for ( ; rbegin != rend; ++rbegin ) {
        if ( *rbegin >= prev || !contains( path, *rbegin ) ) return false;
        prev = *rbegin;
      }
      return true;
    }

  /**
   *  @brief  Check whether this path contains another path.
   *
   *  @param  path The path.
   *  @param  begin The begin iterator of node IDs set of the path to be checked.
   *  @param  end The end iterator of node IDs set of the path to be checked.
   *  @param  idx_begin The begin iterator of the index list of candidates.
   *  @param  idx_end The end iterator of the index list of candidates.
   *  @return `true` if this path is a superset of the given path; otherwise `false`
   *          -- including the case that the given path is empty; i.e.
   *          `end == begin`.
   *
   *  The input path should be smaller than the path. It checks if the nodes of the
   *  given path is present in the node set. It DOES check the order of the nodes.
   *
   *  @overload it uses some candidate indices in the path node vector as hint to find
   *            the location of the first node in the path nodes list.
   */
  template< typename TGraph, typename TSpec, typename TIter1, typename TIter2 >
      inline bool
    contains( const Path< TGraph, TSpec >& path, TIter1 begin, TIter1 end,
        TIter2 idx_begin, TIter2 idx_end )
    {
      if ( begin != end ) {
        for ( ; idx_begin != idx_end; ++idx_begin ) {
          auto lc = path.get_nodes().begin() + *idx_begin;
          if ( std::equal( begin, end, lc ) ) return true;
        }
      }

      return false;
    }  /* -----  end of template function contains  ----- */

  /**
   *  @brief  Check whether a path is covered by a set of paths.
   *
   *  @param  begin The begin iterator of the path as set of node IDs.
   *  @param  end The end iterator of the path as set of node IDs.
   *  @param  path_set_begin The begin iterator of the paths set.
   *  @param  path_set_end The end iterator of the paths set.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  This function checks if all nodes of the given path is covered by at least one
   *  path in the given set of paths. The path set should be a container of the Path
   *  class.
   */
  template< class TIter1, class TIter2 >
      inline bool
    covered_by( TIter1 begin, TIter1 end, TIter2 paths_set_begin, TIter2 paths_set_end )
    {
      for ( ; paths_set_begin != paths_set_end; ++paths_set_begin ) {
        if ( contains( (*paths_set_begin), begin, end ) ) return true;
      }
      return false;
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Check whether a path is covered by a set of paths.
   *
   *  @param  begin The begin iterator of the path as set of node IDs.
   *  @param  end The end iterator of the path as set of node IDs.
   *  @param  paths_set A set of paths as a container of `Path`.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  Overloaded. See `covered_by( TIter1, TIter1, TIter2, TIter2 )`.
   */
  template< class TIter1, class TIter2, class TContainer >
      inline bool
    covered_by( TIter1 begin, TIter2 end, const TContainer& paths_set )
    {
      return covered_by( begin, end, paths_set.begin(), paths_set.end() );
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Check whether a path is covered by a set of paths.
   *
   *  @param  path_nodes The given path to be checked as a vector of node IDs.
   *  @param  paths_set A set of paths as a container of `Path`.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  Overloaded. See `covered_by( TIter1, TIter1, TIter2, TIter2 )`.
   */
  template< typename TNodeID, class TContainer >
      inline bool
    covered_by( const std::vector< TNodeID >& path_nodes, const TContainer& paths_set )
    {
      return covered_by( path_nodes.begin(), path_nodes.end(),
          paths_set.begin(), paths_set.end() );
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Check whether a path is covered by a set of paths.
   *
   *  @param  path The given path as Path class instance.
   *  @param  paths_set A set of paths as a container of `Path`.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  Overloaded. See `covered_by( TIter1, TIter1, TIter2, TIter2 )`.
   */
  template< typename TGraph, typename TSpec, class TContainer >
      inline bool
    covered_by( const Path< TGraph, TSpec >& path, const TContainer& paths_set )
    {
      return covered_by( path.get_nodes().begin(), path.get_nodes().end(),
          paths_set.begin(), paths_set.end() );
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Check whether a node is covered by a set of paths.
   *
   *  @param  node_id The ID of the given node to be checked.
   *  @param  paths_set A set of paths as a container of `Path`.
   *  @return `true` if the node is on one of the path in the paths set, otherwise
   *          `false`.
   *
   *  This function simply checks if a node is on any of the path in the given set of
   *  paths.
   */
  template< typename TNodeID, typename TContainer >
      inline bool
    covered_by( TNodeID node_id, const TContainer& paths_set )
    {
      for ( const auto& path : paths_set ) {
        if ( contains( path, node_id ) ) return true;
      }

      return false;
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Return the path coverage of a given node by its ID.
   *
   *  @param  node_id The node ID of the node to be checked.
   *  @param  paths_set A set of paths as a container of `Path`.
   *  @return the number of paths that cover the given node.
   *
   *  This function get a node ID and return the number of paths that cover the node.
   */
  template< typename TNodeID, typename TContainer >
      inline std::size_t
    get_path_coverage( TNodeID const& node_id, const TContainer& paths_set )
    {
      std::size_t coverage = 0;
      for ( const auto & path : paths_set ) {
        if ( contains( path, node_id ) ) {
          ++coverage;
        }
      }
      return coverage;
    }  /* -----  end of template function get_path_coverage  ----- */
  /* END OF Path interface functions  ------------------------------------------ */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef PATH_INTERFACE_H__  ----- */

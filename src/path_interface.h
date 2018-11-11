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

#include "cpp/vg.pb.h"

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
   *  @brief  Extend the path forward.
   *
   *  @param  path The path.
   *  @param  node_id The new node ID.
   *
   *  Internal function to add the `node_id` to the end of the current path.
   */
  template< typename TGraph, typename TSpec >
      inline void
    add_node( Path< TGraph, TSpec >& path, typename TGraph::nodeid_type const& node_id,
        typename TGraph::offset_type const& node_offset )
    {
      path.push_back( node_id, node_offset );
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
  template< typename TGraph, typename TSpec1, typename TSpec2,
    typename=std::enable_if_t< !is_generic_path< TSpec1 >::value || !is_generic_path< TSpec2 >::value > >
      inline void
    extend( Path< TGraph, TSpec1 >& path, const Path< TGraph, TSpec2 >& other )
    {
      for ( const auto& node_id : other.get_nodes() ) {
        add_node( path, node_id );
      }
    }

  template< typename TGraph, typename TSpec1, typename TSpec2,
    typename=std::enable_if_t< is_generic_path< TSpec1 >::value >,
    typename=std::enable_if_t< is_generic_path< TSpec2 >::value > >
      inline void
    extend( Path< TGraph, TSpec1 >& path, const Path< TGraph, TSpec2 >& other )
    {
      if ( other.empty() ) return;
      if ( path.empty() ) {
        add_node( path, other.front(), other.get_head_offset() );
      }
      else {
        add_node( path, other.front() );
      }
      if ( other.size() > 1 ) {
        for ( auto it = other.begin()+1; it != other.end()-1; ++it ) {
          add_node( path, *it );
        }
        add_node( path, other.back(), other.get_seqlen_tail() );
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
      auto rnk = rank( path, pos );
      auto sel = select( path, rnk );
      assert( pos >= sel );
      return pos - sel + ( rnk ? 0 : path.get_head_offset() );
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
    sequence( const Path< TGraph, TSpec >& path, Forward )
    {
      TGraph const* vargraph = path.get_vargraph();
      assert( vargraph != nullptr );

      if ( path.empty() ) return "";

      typename Path< TGraph, TSpec >::string_type repr;
      repr.reserve( path.get_sequence_len() );

      /* Add sequence of the head node. */
      auto off = path.get_head_offset();
      repr += vargraph->node_sequence( path.front() ).substr( off, path.get_sequence_len() );
      if ( path.size() > 1 ) {
        /* Add sequence of intermediate nodes. */
        for ( auto itr = path.begin()+1; itr != path.end()-1; ++itr ) {
          repr += vargraph->node_sequence( *itr );
        }
        /* Add sequence of the tail node. */
        auto len = path.get_seqlen_tail();
        repr += vargraph->node_sequence( path.back() ).substr( 0, len );
      }

      return repr;
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
    sequence( const Path< TGraph, TSpec >& path, Reversed )
    {
      typename Path< TGraph, TSpec >::string_type repr_str =
        sequence( path, Forward() );
      std::reverse( repr_str.begin(), repr_str.end() );
      return repr_str;
    }  /* -----  end of template function sequence  ----- */

  template< typename TGraph, typename TSpec >
      inline typename Path< TGraph, TSpec >::string_type
    sequence( const Path< TGraph, TSpec >& path )
    {
      return sequence( path, Forward() );
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
    trim_back( Path< TGraph, TSpec >& path, typename TGraph::nodeid_type node_id=0,
        bool exclusive=false )
    {
      bool found = false;
      while ( !found && !path.empty() ) {
        if ( node_id == 0 || path.back() == node_id ) found = true;
        if ( exclusive && found ) break;
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
   *
   *  We have called it "ltrim" because this trimming preserve inclusion of the left-most
   *  k-mer in the resulting path.
   */
  template< typename TGraph, typename TSpec >
      inline void
    ltrim_back_by_len( Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::seqsize_type k, bool hard=false )
    {
      while ( !path.empty() && path.get_sequence_len() - path.get_seqlen_tail() >= k ) {
        pop_back( path );
      }
      if ( !path.empty() && hard ) {
        path.set_right_by_len( k + path.get_seqlen_tail() - path.get_sequence_len() );
      }
    }  /* -----  end of template function ltrim_back_by_len  ----- */

  /**
   *  @brief  Trim the path from the back until further trimming would lead length of
   *          less than `k + <length of first node> - 1`.
   *
   *  @param  path The path to be trimmed.
   *  @param  k The value of k.
   *
   *  The nodes from back will be dropped from the path until to the point that further
   *  trimming leads its length to be less than `k + <length of first node> - 1`.
   *
   *  We have called it "rtrim" because this trimming preserve inclusion of the right-most
   *  k-mer of the first node in the resulting path.
   */
  template< typename TGraph, typename TSpec >
      inline void
    rtrim_back_by_len( Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::seqsize_type k, bool hard=false )
    {
      if ( path.size() < 2 ) return;
      while ( !path.empty() && path.get_sequence_len() -
          path.get_seqlen_head() - path.get_seqlen_tail() >= k-1 ) {
        pop_back( path );
      }
      if ( !path.empty() && hard ) {
        path.set_right_by_len( k - 1 + path.get_seqlen_tail() + path.get_seqlen_head()
            - path.get_sequence_len() );
      }
    }  /* -----  end of template function rtrim_back_by_len  ----- */

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
    trim_front( Path< TGraph, Dynamic >& path, typename TGraph::nodeid_type node_id=0,
        bool exclusive=false )
    {
      bool found = false;
      while ( !found && !path.empty() ) {
        if ( node_id == 0 || path.front() == node_id ) found = true;
        if ( exclusive && found ) break;
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
   *
   *  We have called it "ltrim" because this trimming preserve inclusion of the left-most
   *  k-mer of the last node in the resulting path.
   */
  template< typename TGraph >
      inline void
    ltrim_front_by_len( Path< TGraph, Dynamic >& path,
        typename Path< TGraph, Dynamic >::seqsize_type k, bool hard=false )
    {
      if ( path.size() < 2 ) return;
      while ( !path.empty() && path.get_sequence_len() -
          path.get_seqlen_head() - path.get_seqlen_tail() >= k-1 ) {
        pop_front( path );
      }
      if ( !path.empty() && hard ) {
        path.set_left_by_len( k - 1 + path.get_seqlen_tail() + path.get_seqlen_head()
            - path.get_sequence_len() );
      }
    }  /* -----  end of template function ltrim_front_by_len  ----- */

  /**
   *  @brief  Trim the path from the front until further trimming would lead length of < k.
   *
   *  @param  path The path to be trimmed.
   *  @param  k The length of k.
   *
   *  The nodes from front will be dropped from the path until to the point that further
   *  trimming leads its length to be less than k.
   *
   *  We have called it "rtrim" because this trimming preserve inclusion of the right-most
   *  k-mer in the resulting path.
   */
  template< typename TGraph >
      inline void
    rtrim_front_by_len( Path< TGraph, Dynamic >& path,
        typename Path< TGraph, Dynamic >::seqsize_type k, bool hard=false )
    {
      while ( !path.empty() && path.get_sequence_len() - path.get_seqlen_head() >= k ) {
        pop_front( path );
      }
      if ( !path.empty() && hard ) {
        path.set_left_by_len( k + path.get_seqlen_head() - path.get_sequence_len() );
      }
    }  /* -----  end of template function rtrim_front_by_len  ----- */

  template< typename TGraph, typename TSpec >
      inline YaPair< typename TGraph::nodeid_type, typename TGraph::offset_type >
    leftmost_kmer_pos( Path< TGraph, TSpec > const& path,
        typename Path< TGraph, TSpec >::seqsize_type k )
    {
      YaPair< typename TGraph::nodeid_type, typename TGraph::offset_type > endpos;
      TGraph const* vargraph = path.get_vargraph();
      typename Path< TGraph, TSpec >::seqsize_type len = 0;
      for ( auto it = path.begin(); it != path.end(); ++it ) {
        len += vargraph->node_length( *it );
        if ( len >= k ) {
          endpos.first = *it;
          endpos.second = k + vargraph->node_length( *it ) - len - 1;
          break;
        }
      }
      return endpos;
    }

  template< typename TGraph, typename TSpec >
      inline YaPair< typename TGraph::nodeid_type, typename TGraph::offset_type >
    rightmost_kmer_pos( Path< TGraph, TSpec > const& path,
        typename Path< TGraph, TSpec >::seqsize_type k )
    {
      YaPair< typename TGraph::nodeid_type, typename TGraph::offset_type > startpos;
      TGraph const* vargraph = path.get_vargraph();
      typename Path< TGraph, TSpec >::seqsize_type len = 0;
      for ( auto it = path.end(); it != path.begin(); --it ) {
        len += vargraph->node_length( *(it - 1) );
        if ( len >= k ) {
          startpos.first = *(it - 1);
          startpos.second = len - k;
          break;
        }
      }
      return startpos;
    }

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
      return path.contains( begin, end );
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
      return path.rcontains( rbegin, rend );
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
   *  @param  path The given path to be checked as a set of node IDs.
   *  @param  paths_set A set of paths as a container of `Path`.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  Overloaded. See `covered_by( TIter1, TIter1, TIter2, TIter2 )`.
   */
  template< typename TPath, class TContainer >
      inline bool
    covered_by( const TPath& path, const TContainer& paths_set )
    {
      return covered_by( path.begin(), path.end(), paths_set.begin(), paths_set.end() );
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
  template< typename TNodeID, typename TContainer, typename = std::enable_if_t< std::is_scalar< TNodeID >::value, void > >
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

  /**
   *  @brief  Return the path coverage of a given path.
   *
   *  @param  begin The begin iterator of the nodes list whose coverage should be computed.
   *  @param  end The end iterator of the nodes list whose coverage should be computed.
   *  @param  paths_set A set of paths as a container of `Path` instances.
   *  @return the number of paths that cover the given path.
   *
   *  This function get a path and return the number of paths that cover the path.
   */
  template< typename TIter, typename TContainer >
      inline std::size_t
    get_path_coverage( TIter begin, TIter end, TContainer const& paths_set )
    {
      std::size_t coverage = 0;
      for ( const auto& path : paths_set ) {
        if ( contains( path, begin, end ) ) {
          ++coverage;
        }
      }
      return coverage;
    }  /* -----  end of template function get_path_coverage  ----- */

  /**
   *  @brief  Convert a grem::Path to vg::Path.
   *
   *  @param  path The native path.
   *  @param  vgpath The pointer to an allocated vg path object.
   *
   *  A vg path is a protobuf message which represents a path aligned to the graph
   *  containing a set of Mappings ordered by their ranks. A Mapping is defined for a
   *  loci in the graph (here a first location of each node) defining a set of edits to
   *  transform the node label to the path sequence aligned with that node which is a
   *  full-match here.
   */
  template< typename TGraph, typename TSpec >
      inline void
    convert( Path< TGraph, TSpec > const& path, vg::Path* vgpath )
    {
      TGraph const* vargraph = path.get_vargraph();
      typename TGraph::rank_type rank = 1;
      for ( auto it = path.begin(); it != path.end(); ++it ) {
        typename TGraph::offset_type label_len;
        typename TGraph::offset_type noff = 0;
        if ( it == path.begin() ) {
          label_len = path.get_seqlen_head();
          noff = path.get_head_offset();
        }
        else if ( it == path.end()-1 ) label_len = path.get_seqlen_tail();
        else label_len = vargraph->node_length( *it );
        vg::Mapping* mapping = vgpath->add_mapping();
        mapping->mutable_position()->set_node_id( *it );
        mapping->mutable_position()->set_offset( noff );
        vg::Edit* edit = mapping->add_edit();
        edit->set_from_length( label_len );
        edit->set_to_length( label_len );
        mapping->set_rank( rank++ );
      }
    }

  template< typename TGraph, typename TSpec >
      inline void
    convert( Path< TGraph, TSpec > const& path, vg::Path* vgpath,
        std::vector< vg::Position > const& loci )
    {
      TGraph const* vargraph = path.get_vargraph();
      typename TGraph::rank_type rank = 1;

      auto comp_id =
        [vargraph]( vg::Position const& elem, vg::Position const& value ) {
          return vargraph->id_to_rank( elem ) < vargraph->id_to_rank( value );
        };

      auto comp_both =
        [vargraph]( vg::Position const& elem, vg::Position const& value ) {
          return vargraph->id_to_rank( elem ) < vargraph->id_to_rank( value ) ||
            ( vargraph->id_to_rank( elem ) == vargraph->id_to_rank( value ) &&
              elem.offset() < value.offset() );
        };

      for ( auto it = path.begin(); it != path.end(); ++it ) {
        typename TGraph::offset_type label_len = vargraph->node_length( *it );
        typename TGraph::offset_type coffset = 0;
        if ( it == path.begin() ) coffset = path.get_head_offset();
        else if ( it == path.end()-1 ) label_len = path.get_seqlen_tail();
        vg::Mapping* mapping = vgpath->add_mapping();
        mapping->mutable_position()->set_node_id( *it );
        mapping->mutable_position()->set_offset( coffset );

        vg::Position cpos;
        cpos.set_node_id( *it );
        cpos.set_offset( coffset );
        std::vector< vg::Position>::const_iterator nextedit, lastedit;
        if ( it == path.end()-1 ) {
          nextedit = std::lower_bound( loci.begin(), loci.end(), cpos, comp_id );
          lastedit = std::upper_bound( loci.begin(), loci.end(), cpos, comp_both );
        }
        else {
          nextedit = std::lower_bound( loci.begin(), loci.end(), cpos, comp_both );
          lastedit = std::upper_bound( loci.begin(), loci.end(), cpos, comp_id );
        }
        typename TGraph::offset_type toffset =
              nextedit != lastedit ?
              (*nextedit).offset() :
              label_len;
        do {
          vg::Edit* edit = mapping->add_edit();
          if ( coffset > toffset ) {
            ++nextedit;
            toffset =
              nextedit != lastedit ?
              (*nextedit).offset() :
              label_len;
          }
          if ( coffset == toffset ) {
            edit->set_from_length( 1 );
            edit->set_to_length( 1 );
            edit->set_sequence( "S" );
            ++coffset;
          }
          else {
            edit->set_from_length( toffset - coffset );
            edit->set_to_length( toffset - coffset );
            coffset = toffset;
          }
        } while ( coffset < label_len );
        mapping->set_rank( rank++ );
      }
    }
  /* END OF Path interface functions  ------------------------------------------ */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef PATH_INTERFACE_H__  ----- */

/**
 *    @file  path.h
 *   @brief  Path template class definitions.
 *
 *  This header file defines Path template class in variation graph, its specializations
 *  for different purposes, and some interface functions to work with these paths.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sun Nov 19, 2017  14:38
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PATH_H__
#define PATH_H__

#include <fstream>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <string>
#include <unordered_set>
#include <deque>
#include <algorithm>
#include <utility>

#include <seqan/basic.h>
#include <sdsl/bit_vectors.hpp>

#include "utils.h"


namespace grem{
  /* Path specialization tags. */
  struct DefaultStrategy;
  struct DynamicStrategy;
  struct CompactStrategy;
  typedef seqan::Tag< DefaultStrategy > Default;
  typedef seqan::Tag< DynamicStrategy > Dynamic;
  typedef seqan::Tag< CompactStrategy > Compact;

  template< typename TGraph, typename TSpec >
    struct PathTraits;

  template< typename TGraph >
    struct PathTraits< TGraph, Default > {
      typedef std::vector< typename TGraph::nodeid_type > TNodeSequence;
    };

  template< typename TGraph >
    struct PathTraits< TGraph, Dynamic > {
      typedef std::deque< typename TGraph::nodeid_type > TNodeSequence;
    };

  /* Path interface functions forwards  ---------------------------------------- */

  template< typename TGraph, typename TSpec >
    class Path;

  template< typename TGraph, typename TSpec >
      void
    initialize( Path< TGraph, TSpec >& path );
  template< typename TGraph, typename TSpec >
      void
    save( const Path< TGraph, TSpec >& path, std::ostream& out );
  template< typename TGraph, typename TSpec >
      void
    load( Path< TGraph, TSpec >& path, std::istream& in );
  template< typename TGraph, typename TSpec >
      void
    add_node( Path< TGraph, TSpec >& path,
        typename TGraph::nodeid_type const& node_id );
  template< typename TGraph, typename TSpec >
      typename Path< TGraph, TSpec >::size_type
    rank( const Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::seqsize_type pos );
  template< typename TGraph, typename TSpec >
      typename Path< TGraph, TSpec >::seqsize_type
    select( const Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::size_type rank );
  template< typename TGraph, typename TSpec >
      typename Path< TGraph, TSpec >::string_type
    sequence( const Path< TGraph, TSpec >& path );
  template< typename TGraph, typename TSpec >
      void
    clear( Path< TGraph, TSpec >& path );
  template< typename TGraph, typename TSpec >
      void
    reserve( Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::size_type size );
  template< typename TGraph, typename TSpec >
      typename Path< TGraph, TSpec >::size_type
    length( const Path< TGraph, TSpec >& path );
  template< typename TGraph, typename TSpec >
      bool
    contains( const Path< TGraph, TSpec >& path,
        typename TGraph::nodeid_type node_id );
  template< typename TGraph, typename TSpec >
      void
    pop_back( Path< TGraph, TSpec >& path );
  template< typename TGraph >
      void
    pop_front( Path< TGraph, Dynamic >& path );

  /* END OF path interface functions forwards  --------------------------------- */

  /**
   *  @brief  Path template class.
   *
   *  Represent a path in the variation graph with efficient node ID query.
   */
  template< typename TGraph, typename TSpec = Default >
    class Path
    {
      private:
        /* ====================  TYPEDEFS      ======================================= */
        typedef PathTraits< TGraph, TSpec > TTraits;
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef std::string string_type;
        typedef string_type::size_type seqsize_type;
        typedef typename TTraits::TNodeSequence nodes_type;
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        const TGraph* vargraph;
        nodes_type nodes;
        std::unordered_set< typename TGraph::nodeid_type > nodes_set;
        seqsize_type seqlen;
        /* Loaded on demand. */
        string_type seq;
        /* Loaded after calling `initialize`. */
        bool initialized;
        sdsl::bit_vector bv_node_breaks;
        sdsl::rank_support_v<> rs_node_breaks;
        sdsl::select_support_mcl<> ss_node_breaks;
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename decltype( nodes )::size_type size_type;
        /* ====================  LIFECYCLE     ======================================= */
        Path( const TGraph* g )
          : vargraph( g ), seqlen( 0 ), initialized( false )
        { }

        Path( const TGraph* g, nodes_type&& p )
          : Path( g )
        {
          this->set_nodes( std::move( p ) );
        }

        Path( const TGraph* g, const nodes_type& p )
          : Path( g, nodes_type( p ) )
        { }

        Path( const Path& ) = default;
        Path( Path&& ) = default;
        Path& operator=( const Path& ) = default;
        Path& operator=( Path&& ) = default;
        ~Path() = default;
        /* ====================  OPERATORS     ======================================= */
        /**
         *  @brief  Assign the path with another type of path.
         *
         *  @param  other The other path to use for assignment.
         *  @return The reference to this instance after assignment.
         *
         *  @note If the type of the `other` path is the same as this path, the copy
         *        assignment operator would be used.
         */
        template< typename TSpec2 >
            inline Path&
          operator=( const Path< TGraph, TSpec2 >& other ) {
            if ( this->vargraph != other.get_vargraph() ) {
              throw std::runtime_error( "Mismatching variation graphs" );
            }
            nodes_type d;
            std::copy( other.get_nodes().begin(), other.get_nodes().end(), std::back_inserter( d ) );
            this->set_nodes( std::move( d ) );
            return *this;
          }
        /* ====================  ACCESSORS     ======================================= */
        /**
         *  @brief  getter function for vargraph.
         */
          inline const TGraph*
        get_vargraph( ) const
        {
          return this->vargraph;
        }  /* -----  end of method get_vargraph  ----- */

        /**
         *  @brief  getter function for nodes.
         */
          inline const nodes_type&
        get_nodes( ) const
        {
          return this->nodes;
        }  /* -----  end of method get_nodes  ----- */

        /**
         *  @brief  getter function for seqlen.
         */
          inline seqsize_type
        get_sequence_len( ) const
        {
          return this->seqlen;
        }  /* -----  end of method get_sequence_len  ----- */

        /**
         *  @brief  getter function for seq.
         *
         *  @note The sequence is constructed on demand by calling this function.
         */
          inline const string_type&
        get_sequence( )
        {
          if ( this->seq.length() == 0 ) this->seq = sequence( *this );
          return this->seq;
        }  /* -----  end of method get_sequence  ----- */

        /**
         *  @brief  Is path initialized?
         *
         *  @return `true` if the path is initialized; `false` otherwise.
         *
         *  Initialization constructs node breaks bit vector, rank, and select supports.
         */
          inline bool
        is_initialized( ) const
        {
          return this->initialized;
        }
        /* ====================  MUTATORS      ======================================= */
        /**
         *  @brief  setter function for vargraph.
         */
          inline void
        set_vargraph( const TGraph* value )
        {
          this->vargraph = value;
        }  /* -----  end of method set_vargraph  ----- */

        /**
         *  @brief  setter function for nodes.
         */
          inline void
        set_nodes( nodes_type&& value )
        {
          assert( path.vargraph != nullptr );

          this->nodes = std::move( value );
          this->nodes_set.clear();
          this->nodes_set.reserve( this->nodes.size() );
          this->seqlen = 0;
          for ( const auto& node_id : this->nodes ) {
            this->nodes_set.insert( node_id );
            this->seqlen += this->vargraph->node_length( node_id );
          }
          this->seq.clear();
          this->initialized = false;
        }  /* -----  end of method set_nodes  ----- */

        /**
         *  @brief  setter function for nodes.
         */
          inline void
        set_nodes( const nodes_type& value )
        {
          this->set_nodes( nodes_type( value ) );
        }  /* -----  end of method set_nodes  ----- */
        /* ====================  INTERFACE FUNCTIONS  ================================ */
          friend void
        initialize< TGraph, TSpec >( Path< TGraph, TSpec >& path );
          friend void
        save< TGraph, TSpec >( const Path< TGraph, TSpec >& path, std::ostream& out );
          friend void
        load< TGraph, TSpec >( Path< TGraph, TSpec >& path, std::istream& in );
          friend void
        add_node< TGraph, TSpec >( Path< TGraph, TSpec >& path,
            typename TGraph::nodeid_type const& node_id );
          friend size_type
        rank< TGraph, TSpec >( const Path< TGraph, TSpec >& path, seqsize_type i );
          friend seqsize_type
        select< TGraph, TSpec >( const Path< TGraph, TSpec >& path, size_type rank );
          friend void
        clear< TGraph, TSpec >( Path< TGraph, TSpec >& path );
          friend void
        reserve< TGraph, TSpec >( Path< TGraph, TSpec >& path, size_type size );
          friend bool
        contains< TGraph, TSpec >( const Path< TGraph, TSpec >& path,
            typename TGraph::nodeid_type node_id );
          friend void
        pop_back< TGraph, TSpec >( Path< TGraph, TSpec >& path );
          friend void
        pop_front< TGraph >( Path< TGraph, Dynamic >& path );
    };  /* -----  end of template class Path  ----- */

  /* Normal Path interface functions  ------------------------------------------ */

  /**
   *  @brief  Initialize data structure for efficient rank and select queries.
   *
   *  @param  path The path.
   *
   *  It constructs node breaks bit vector and corresponding rank and select supports.
   *  The node breaks bit vector is a bit vector of sequence length. A bit at position
   *  `i` is set if a node starts from that position in the sequence except for the
   *  first node. For example for this path:
   *
   *  (GCAAT) -> (A) -> (TTAGCC) -> (GCA)
   *
   *  the corresponding path is:
   *
   *  GCAATATTAGCCGCA
   *
   *  and its bit vector is:
   *
   *  000001100000100
   *
   *  which has a set bit at the first position of each node in the path.
   */
  template< typename TGraph, typename TSpec >
      inline void
    initialize( Path< TGraph, TSpec >& path )
    {
      if ( path.is_initialized() ) return;

      assert( path.vargraph != nullptr );

      sdsl::util::assign( path.bv_node_breaks,
          sdsl::bit_vector( path.get_sequence_len(), 0 ) );
      typename Path< TGraph, TSpec >::seqsize_type cursor = 0;
      for ( const auto& node_id : path.nodes ) {
        cursor += path.vargraph->node_length( node_id );
        path.bv_node_breaks[ cursor - 1 ] = 1;
      }
      sdsl::util::init_support( path.rs_node_breaks, &path.bv_node_breaks );
      sdsl::util::init_support( path.ss_node_breaks, &path.bv_node_breaks );
      path.initialized = true;
    }

  /**
   *  @brief  Save the path to an output stream.
   *
   *  @param  path The path.
   *  @param  out The output stream.
   *
   *  It saves the sequence of the node IDs into the given output stream.
   */
  template< typename TGraph, typename TSpec >
      inline void
    save( const Path< TGraph, TSpec >& path, std::ostream& out )
    {
      assert( path.is_initialized() );
      serialize( out, path.nodes, path.nodes.begin(), path.nodes.end() );
      serialize( path.bv_node_breaks, out );
    }  /* -----  end of function save  ----- */

  /**
   *  @brief  Save the path to an output stream.
   *
   *  @param  path The path.
   *  @param  out The output stream.
   *
   *  It saves the sequence of the node IDs into the given output stream.
   */
  template< typename TGraph, typename TSpec >
      inline void
    save( Path< TGraph, TSpec >& path, std::ostream& out )
    {
      initialize( path );
      save( std::add_const_t< Path< TGraph, TSpec >& >( path ), out );
    }  /* -----  end of function save  ----- */

  /**
   *  @brief  Save the path to file.
   *
   *  @param  path The path.
   *  @param  file_name The name of the file to be written.
   *
   *  It saves the sequence of the node IDs into the file.
   *
   *  @overload for `std::string` as file name.
   */
  template< typename TGraph, typename TSpec >
      inline void
    save( const Path< TGraph, TSpec >& path, const std::string& file_name )
    {
      std::ofstream ofs( file_name, std::ofstream::out | std::ofstream::binary );
      if( !ofs ) {
        throw std::runtime_error( "cannot open file '" + file_name + "'" );
      }
      save( path, ofs );
    }  /* -----  end of function save  ----- */

  /**
   *  @brief  Save the path to file.
   *
   *  @param  path The path.
   *  @param  file_name The name of the file to be written.
   *
   *  It saves the sequence of the node IDs into the file.
   *
   *  @overload for non-const Path instance.
   */
  template< typename TGraph, typename TSpec >
      inline void
    save( Path< TGraph, TSpec >& path, const std::string& file_name )
    {
      std::ofstream ofs( file_name, std::ofstream::out | std::ofstream::binary );
      if( !ofs ) {
        throw std::runtime_error( "cannot open file '" + file_name + "'" );
      }
      save( path, ofs );
    }  /* -----  end of function save  ----- */

  /**
   *  @brief  Load the path from an input stream.
   *
   *  @param  path The path.
   *  @param  in The input stream.
   *
   *  It loads the sequence of the node IDs from the given input stream.
   */
  template< typename TGraph, typename TSpec >
      inline void
    load( Path< TGraph, TSpec >& path, std::istream& in )
    {
      typename Path< TGraph, TSpec >::nodes_type nodes;
      deserialize( in, nodes, std::back_inserter( nodes ) );
      path.set_nodes( std::move( nodes ) );
      load( path.bv_node_breaks, in );
      sdsl::util::init_support( path.rs_node_breaks, &path.bv_node_breaks );
      sdsl::util::init_support( path.ss_node_breaks, &path.bv_node_breaks );
      path.initialized = true;
    }  /* -----  end of template function load  ----- */

  /**
   *  @brief  Load the path from file.
   *
   *  @param  path The path.
   *  @param  file_name The name of the file to be read.
   *
   *  It loads the sequence of the node IDs from the file.
   *
   *  @overload for `std::string` as file name.
   */
  template< typename TGraph, typename TSpec >
      inline void
    load( Path< TGraph, TSpec >& path, const std::string& file_name )
    {
      std::ifstream ifs( file_name, std::ifstream::in | std::ifstream::binary );
      if( !ifs ) {
        throw std::runtime_error( "cannot open file '" + file_name + "'" );
      }
      load( path, ifs );
    }  /* -----  end of template function load  ----- */

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
      assert( path.vargraph != nullptr );
      path.nodes.push_back( node_id );
      path.nodes_set.insert( node_id );
      path.seqlen += path.vargraph->node_length( node_id );
      if ( path.seq.length() != 0 ) {
        path.seq += path.vargraph->node_sequence( node_id );
      }
      path.initialized = false;
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
      if ( path.get_vargraph() != other.get_vargraph() ) {
        throw std::runtime_error( "Mismatching variation graphs" );
      }
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

  /**
   *  @brief  Get the node rank (0-based) in the path by a position in its sequence.
   *
   *  @param  path The path.
   *  @param  pos The position in the path (0-based).
   *  @return The node rank in the nodes queue (0-based) on whose label the position
   *          `pos` in the path sequence relies.
   *
   *  The value of `rank(pos)` is the node rank in the path's nodes queue.
   */
  template< typename TGraph, typename TSpec >
      inline typename Path< TGraph, TSpec >::size_type
    rank( const Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::seqsize_type pos )
    {
      assert( path.is_initialized() );
      if ( pos < 0 || pos >= path.get_sequence_len() ) {
        throw std::runtime_error( "Position out of range." );
      }
      return path.rs_node_breaks( pos );
    }

  /**
   *  @brief  Get the position in the sequence from which the node with rank `rank` starts.
   *
   *  @param  path The path.
   *  @return The position in the path sequence from which the node with rank `rank`
   *          starts.
   *
   *  Since the position corresponding to the last base in the node label is set in the
   *  bit vector, the value of `select(rank)` is one base before the actual position to
   *  be returned. So, `select(rank) + 1` would be the desired position.
   */
  template< typename TGraph, typename TSpec >
      inline typename Path< TGraph, TSpec >::seqsize_type
    select( const Path< TGraph, TSpec >& path,
        typename Path< TGraph, TSpec >::size_type rank )
    {
      assert( path.is_initialized() );
      if ( rank < 0 || rank >= length( path ) ) {
        throw std::runtime_error( "Rank out of range." );
      }
      if ( rank == 0 ) return 0;
      return path.ss_node_breaks( rank ) + 1;
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
      return path.get_nodes().at( rank( path, pos ) );
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
    sequence( const Path< TGraph, TSpec >& path, Forward )
    {
      assert( path.get_vargraph() != nullptr );

      typename Path< TGraph, TSpec >::string_type repr_str;
      repr_str.reserve( path.get_sequence_len() );
      for ( const auto& node_id : path.get_nodes() ) {
        repr_str += path.get_vargraph()->node_sequence( node_id );
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
    sequence( const Path< TGraph, TSpec >& path, Reversed )
    {
      typename Path< TGraph, TSpec >::string_type repr_str = sequence( path, Forward() );
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
      path.nodes.clear();
      path.nodes_set.clear();
      path.seqlen = 0;
      path.seq.clear();
      sdsl::util::clear( path.bv_node_breaks );
      sdsl::util::clear( path.rs_node_breaks );
      sdsl::util::clear( path.ss_node_breaks );
      path.initialized = false;
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
      grem::reserve( path.nodes, size );        /* defined in `utils.h` */
      path.nodes_set.reserve( size );
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
      return path.get_nodes().size();
    }  /* -----  end of template function length  ----- */

  /**
   * @brief  pop the last node from the path.
   */
  template< typename TGraph, typename TSpec >
      inline void
    pop_back( Path< TGraph, TSpec >& path )
    {
      assert( path.vargraph != nullptr );

      if ( path.nodes.empty() ) return;

      auto&& last_node = path.nodes.back();
      auto&& last_node_len = path.vargraph->node_length( last_node );
      assert( path.seqlen >= last_node_len );
      path.seqlen -= last_node_len;
      path.nodes_set.erase( path.nodes_set.find( last_node ) );
      path.nodes.pop_back();
      path.initialized = false;

      if ( path.seq.length() != 0 ) path.seq.resize( path.get_sequence_len() );
    }  /* -----  end of method pop_back  ----- */

  /**
   * @brief  pop the last node from the path.
   *
   * NOTE: only implemented for Dynamic Path.
   */
  template< typename TGraph >
      inline void
    pop_front( Path< TGraph, Dynamic >& path )
    {
      assert( path.vargraph != nullptr );

      if ( path.nodes.empty() ) return;

      auto&& first_node = path.nodes.front();
      auto&& first_node_len = path.vargraph->node_length( first_node );
      assert( path.seqlen >= first_node_len );
      path.seqlen -= first_node_len;
      path.nodes_set.erase( path.nodes_set.find( first_node ) );
      path.nodes.pop_front();
      path.initialized = false;

      if ( path.seq.length() != 0 ) path.seq = path.seq.substr( first_node_len );
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
      while ( path.get_sequence_len() -
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
      while ( path.get_sequence_len() -
          path.get_vargraph()->node_length( path.get_nodes().front() ) >= k ) {
        pop_front( path );
      }
    }  /* -----  end of template function trim_front_by_len  ----- */

  /* END OF Normal Path interface functions  ----------------------------------- */

  /* Compact Path interface functions forwards  -------------------------------- */

  template< typename TGraph >
      inline void
    save( const Path< TGraph, Compact >& path, std::ostream& out );
  template< typename TGraph >
      inline void
    add_node( Path< TGraph, Compact >& path, typename TGraph::nodeid_type const& node_id );
  template< typename TGraph >
      inline void
    clear( Path< TGraph, Compact >& path );
  template< typename TGraph >
      inline void
    reserve( Path< TGraph, Compact >& path,
        typename Path< TGraph, Compact >::size_type size );
  template< typename TGraph >
      inline typename Path< TGraph, Compact >::size_type
    length( const Path< TGraph, Compact >& path );

  /* END OF Compact Path interface functions forwards  ------------------------- */

  /**
   *  @brief  Path template class Compact specialization.
   *
   *  Represent a path in the variation graph with efficient node ID query.
   */
  template< typename TGraph >
    class Path< TGraph, Compact >
    {
      public:
        /* ====================  DATA MEMBERS  ======================================= */
        std::unordered_set< typename TGraph::nodeid_type > nodes_set;
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename decltype( nodes_set )::size_type size_type;
        /* ====================  LIFECYCLE     ======================================= */
        Path( ) = default;
        Path( const std::vector< typename TGraph::nodeid_type >& p )
        {
          this->set_nodes( p );
        }
        /* ====================  MUTATORS      ======================================= */
        /**
         *  @brief  Set the nodes in the path (clear the previous state).
         */
          inline void
        set_nodes( const std::vector< typename TGraph::nodeid_type >& value )
        {
          this->nodes_set.clear();
          this->nodes_set.reserve( value.size() );
          std::copy( value.begin(), value.end(),
              std::inserter( this->nodes_set, this->nodes_set.end() ) );
        }  /* -----  end of method set_nodes  ----- */
        /* ====================  INTERFACE FUNCTIONS  ================================ */
          friend void
        save< TGraph >( const Path< TGraph, Compact >& path, std::ostream& out );
          friend void
        add_node< TGraph >( Path< TGraph, Compact >& path,
            typename TGraph::nodeid_type const& node_id );
          friend void
        clear< TGraph >( Path< TGraph, Compact >& path );
          friend void
        reserve< TGraph >( Path< TGraph, Compact >& path, size_type size );
          friend size_type
        length< TGraph >( const Path< TGraph, Compact >& path );
          friend bool
        contains< TGraph, Compact >( const Path< TGraph, Compact >& path,
            typename TGraph::nodeid_type node_id );
    };  /* -----  end of specialized template class Path  ----- */

  /* Compact Path interface functions definitions  ----------------------------- */

  /**
   *  @brief  Save the path to an output stream -- Compact specialization.
   *
   *  @param  path The path.
   *  @param  out The output stream.
   *
   *  It saves the sequence of the node IDs into the given output stream.
   */
  template< typename TGraph >
      inline void
    save( const Path< TGraph, Compact >& path, std::ostream& out )
    {
      serialize( out, path.nodes_set, path.nodes_set.begin(), path.nodes_set.end() );
    }  /* -----  end of function save  ----- */

  /**
   *  @brief  Save the path to file -- Compact specialization.
   *
   *  @param  path The path.
   *  @param  file_name The name of the file to be written.
   *
   *  It saves the sequence of the node IDs into the file.
   */
  template< typename TGraph >
      inline void
    save( const Path< TGraph, Compact >& path, const std::string& file_name )
    {
      std::ofstream ofs( file_name, std::ofstream::out | std::ofstream::binary );
      if( !ofs ) {
        throw std::runtime_error( "cannot open file '" + file_name + "'" );
      }
      save( path, ofs );
    }  /* -----  end of function save  ----- */

  /**
   *  @brief  Extend the path.
   *
   *  @param  path The path.
   *  @param  node_id The new node ID to be added into the path node set.
   */
  template< typename TGraph >
      inline void
    add_node( Path< TGraph, Compact >& path,
        typename TGraph::nodeid_type const& node_id )
    {
      path.nodes_set.insert( node_id );
    }

  /**
   *  @brief  Clear the path.
   *
   *  @param  path The path.
   */
  template< typename TGraph >
      inline void
    clear( Path< TGraph, Compact >& path )
    {
      path.nodes_set.clear();
    }

  /**
   *  @brief  Reserve memory for the path.
   *
   *  @param  path The path.
   *  @param  size The target capacity.
   *
   *  It reserves memory for nodes set.
   */
  template< typename TGraph >
      inline void
    reserve( Path< TGraph, Compact >& path,
        typename Path< TGraph, Compact >::size_type size )
    {
      path.nodes_set.reserve( size );
    }

  /**
   *  @brief  Get path length.
   *
   *  @param  path The path.
   *  @return The path length.
   *
   *  It gets the length of the path.
   */
  template< typename TGraph >
      inline typename Path< TGraph, Compact >::size_type
    length( const Path< TGraph, Compact >& path )
    {
      return path.nodes_set.size();
    }  /* -----  end of function length  ----- */

  /**
   *  @brief  Load a Compact Path from an input stream.
   *
   *  @param  path The path.
   *  @param  in The input stream.
   *
   *  It loads the sequence of the node IDs from the given input stream.
   *
   *  Specialized for Compact Path.
   */
  template< typename TGraph, typename TSpec >
      inline void
    load( Path< TGraph, Compact >& path, std::istream& in )
    {
      std::vector< typename TGraph::nodeid_type > nodes;
      deserialize( in, nodes, std::back_inserter( nodes ) );
      path.set_nodes( std::move( nodes ) );
    }

  /* END OF Compact Path interface functions  ---------------------------------- */

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
      return path.nodes_set.find( node_id ) != path.nodes_set.end();
    }  /* -----  end of template function contains  ----- */

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
   *  given path is present in the node set. It does not check the order of the
   *  nodes.
   */
  template< typename TGraph, typename TSpec, typename TIter >
      inline bool
    contains( const Path< TGraph, TSpec >& path, TIter begin, TIter end )
    {
      if ( end - begin > 0 ) {
        auto&& on_path = [&path]( typename TGraph::nodeid_type const& i ) {
          return contains( path, i );
        };

        if ( std::all_of( begin, end, on_path ) ) return true;
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
      for ( auto itr = paths_set_begin; itr != paths_set_end; ++itr ) {
        if ( contains( (*itr), begin, end ) )
        {
          return true;
        }
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

#endif  /* ----- #ifndef PATH_H__  ----- */

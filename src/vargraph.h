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
#include <ios>
#include <stdexcept>
#include <exception>
#include <string>
#include <vector>
#include <unordered_set>
#include <deque>
#include <utility>
#include <cstdint>
#include <functional>
#include <random>

#ifdef HAVE_MEMCPY
#  define MACRO_STACK
#  pragma push_macro("HAVE_MEMCPY")
#  undef HAVE_MEMCPY
#endif

#include <xg/xg.hpp>

#ifdef MACRO_STACK
#  pragma pop_macro("HAVE_MEMCPY")
#  undef MACRO_STACK
#endif

#include <seqan/basic.h>
#include <sdsl/bit_vectors.hpp>
#include "stream/stream.h"

#include "graph_iter.h"
#include "utils.h"

// :TODO:Tue Sep 26 17:03:\@cartoonist: split the header file into multiple ones?


namespace grem
{
  class VarGraph : public xg::XG
  {
    public:
      // typedefs
      typedef vg::Node node_type;                               /**< @brief Node type. */
      typedef decltype( vg::Node().id() ) nodeid_type;          /**< @brief Node ID type. */
      typedef std::size_t rank_type;                            /**< @brief Node ID type. */
      typedef decltype( vg::Position().offset() ) offset_type;  /**< @brief Node offset type. */

      // Constructors
      VarGraph( void ) : XG( ) { }

      VarGraph( std::ifstream &ifs, bool fmt_xg = true )
      {
        if ( fmt_xg ) {
          load( ifs );
        }
        else {
          from_stream( ifs );
        }
      }

      VarGraph( vg::Graph &vg_graph ) : XG( vg_graph )
      { }

      // Public methods
        inline bool
      is_branch ( nodeid_type node_id ) const
      {
        if ( this->edges_from( node_id ).size() > 1 ) {
          return true;
        }
        return false;
      }  /* -----  end of method is_branch  ----- */

        inline bool
      is_merge ( nodeid_type node_id ) const
      {
        if ( this->edges_to( node_id ).size() > 1 ) {
          return true;
        }
        return false;
      }  /* -----  end of method is_merge  ----- */

        inline bool
      has_edges_from( nodeid_type node_id ) const
      {
        return this->edges_from( node_id ).size() != 0;
      }

        inline bool
      has_edges_to( nodeid_type node_id ) const
      {
        return this->edges_to( node_id ).size() != 0;
      }
  };

  /* Path interface functions forwards  ---------------------------------------- */

  template< typename TSpec >
    class Path;

  template< typename TSpec >
      void
    initialize( Path< TSpec >& path );
  template< typename TSpec >
      void
    save( const Path< TSpec >& path, std::ostream& out );
  template< typename TSpec >
      void
    load( Path< TSpec >& path, std::istream& in );
  template< typename TSpec >
      void
    add_node( Path< TSpec >& path, const VarGraph::nodeid_type& node_id );
  template< typename TSpec >
      typename Path< TSpec >::size_type
    rank( const Path< TSpec >& path, typename Path< TSpec >::seqsize_type pos );
  template< typename TSpec >
      typename Path< TSpec >::seqsize_type
    select( const Path< TSpec >& path, typename Path< TSpec >::size_type rank );
  template< typename TSpec >
      typename Path< TSpec >::string_type
    sequence( const Path< TSpec >& path );
  template< typename TSpec >
      void
    clear( Path< TSpec >& path );
  template< typename TSpec >
      void
    reserve( Path< TSpec >& path, typename Path< TSpec >::size_type size );
  template< typename TSpec >
      typename Path< TSpec >::size_type
    length( const Path< TSpec >& path );
  template< typename TSpec >
      bool
    contains( const Path< TSpec >& path, VarGraph::nodeid_type node_id );
  template< typename TSpec >
      void
    pop_back( Path< TSpec >& path );
  template< typename TSpec >
      void
    pop_front( Path< TSpec >& path );

  /* END OF path interface functions forwards  --------------------------------- */

  /* Path specialization tags. */
  struct DefaultStrategy;
  struct DynamicStrategy;
  struct CompactStrategy;
  typedef seqan::Tag< DefaultStrategy > Default;
  typedef seqan::Tag< DynamicStrategy > Dynamic;
  typedef seqan::Tag< CompactStrategy > Compact;

  template< typename TSpec >
    struct PathTraits;

  template< >
    struct PathTraits< Default > {
      typedef std::vector< VarGraph::nodeid_type > TNodeSequence;
    };

  template< >
    struct PathTraits< Dynamic > {
      typedef std::deque< VarGraph::nodeid_type > TNodeSequence;
    };

  /**
   *  @brief  Path template class.
   *
   *  Represent a path in the variation graph with efficient node ID query.
   */
  template< typename TSpec = Default >
    class Path
    {
      private:
        /* ====================  TYPEDEFS      ======================================= */
        typedef PathTraits< TSpec > TTraits;
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef std::string string_type;
        typedef string_type::size_type seqsize_type;
        typedef typename TTraits::TNodeSequence nodes_type;
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        const VarGraph* vargraph;
        nodes_type nodes;
        std::unordered_set< VarGraph::nodeid_type > nodes_set;
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
        Path( const VarGraph* g )
          : vargraph( g ), seqlen( 0 ), initialized( false )
        { }

        Path( const VarGraph* g, nodes_type&& p )
          : Path( g )
        {
          this->set_nodes( std::move( p ) );
        }

        Path( const VarGraph* g, const nodes_type& p )
          : Path( g, nodes_type( p ) )
        { }
        /* ====================  ACCESSORS     ======================================= */
        /**
         *  @brief  getter function for vargraph.
         */
          inline const VarGraph*
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
        set_vargraph( const VarGraph* value )
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
        initialize< TSpec >( Path< TSpec >& path );
          friend void
        save< TSpec >( const Path< TSpec >& path, std::ostream& out );
          friend void
        load< TSpec >( Path< TSpec >& path, std::istream& in );
          friend void
        add_node< TSpec >( Path< TSpec >& path, const VarGraph::nodeid_type& node_id );
          friend size_type
        rank< TSpec >( const Path< TSpec >& path, seqsize_type i );
          friend seqsize_type
        select< TSpec >( const Path< TSpec >& path, size_type rank );
          friend void
        clear< TSpec >( Path< TSpec >& path );
          friend void
        reserve< TSpec >( Path< TSpec >& path, size_type size );
          friend bool
        contains< TSpec >( const Path< TSpec >& path, VarGraph::nodeid_type node_id );
          friend void
        pop_back< TSpec >( Path< TSpec >& path );
          friend void
        pop_front< TSpec >( Path< TSpec >& path );
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
  template< typename TSpec >
      inline void
    initialize( Path< TSpec >& path )
    {
      if ( path.is_initialized() ) return;

      assert( path.vargraph != nullptr );

      sdsl::util::assign( path.bv_node_breaks, bit_vector( path.get_sequence_len(), 0 ) );
      typename Path< TSpec >::seqsize_type cursor = 0;
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
  template< typename TSpec >
      inline void
    save( const Path< TSpec >& path, std::ostream& out )
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
  template< typename TSpec >
      inline void
    save( Path< TSpec >& path, std::ostream& out )
    {
      initialize( path );
      save( std::add_const_t< Path< TSpec >& >( path ), out );
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
  template< typename TSpec >
      inline void
    save( const Path< TSpec >& path, const std::string& file_name )
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
  template< typename TSpec >
      inline void
    save( Path< TSpec >& path, const std::string& file_name )
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
  template< typename TSpec >
      inline void
    load( Path< TSpec >& path, std::istream& in )
    {
      typename Path< TSpec >::nodes_type  nodes;
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
  template< typename TSpec >
      inline void
    load( Path< TSpec >& path, const std::string& file_name )
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
  template< typename TSpec >
      inline void
    add_node( Path< TSpec >& path, const VarGraph::nodeid_type& node_id )
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
   *  @brief  Get the node rank (0-based) in the path by a position in its sequence.
   *
   *  @param  path The path.
   *  @param  pos The position in the path (0-based).
   *  @return The node rank in the nodes queue (0-based) on whose label the position
   *          `pos` in the path sequence relies.
   *
   *  The value of `rank(pos)` is the node rank in the path's nodes queue.
   */
  template< typename TSpec >
      inline typename Path< TSpec >::size_type
    rank( const Path< TSpec >& path, typename Path< TSpec >::seqsize_type pos )
    {
      assert( path.is_initialized() );
      assert( 0 <= pos && pos < path.get_sequence_len() );
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
  template< typename TSpec >
      inline typename Path< TSpec >::seqsize_type
    select( const Path< TSpec >& path, typename Path< TSpec >::size_type rank )
    {
      assert( path.is_initialized() );
      assert( 0 <= rank && rank < length( path ) );
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
  template< typename TSpec >
      inline VarGraph::nodeid_type
    position_to_id( const Path< TSpec >& path, typename Path< TSpec >::seqsize_type pos )
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
  template< typename TSpec >
      inline VarGraph::offset_type
    position_to_offset( const Path< TSpec >& path, typename Path< TSpec >::seqsize_type pos )
    {
      auto&& sel = select( path, rank( path, pos ) );
      assert( pos >= sel );
      return pos - sel;
    }

  /**
   *  @brief  Compute sequence from the nodes queue.
   *
   *  @param  path The path.
   *  @return Sequence represented by the path.
   *
   *  Compute the sequence of the path from nodes queue.
   */
  template< typename TSpec >
      inline typename Path< TSpec >::string_type
    sequence( const Path< TSpec >& path )
    {
      assert( path.get_vargraph() != nullptr );

      typename Path< TSpec >::string_type repr_str;
      repr_str.reserve( path.get_sequence_len() );
      for ( const auto& node_id : path.get_nodes() ) {
        repr_str += path.get_vargraph()->node_sequence( node_id );
      }
      return repr_str;
    }  /* -----  end of template function sequence  ----- */

  /**
   *  @brief  Clear the path.
   *
   *  @param  path The path.
   *
   *  Internal function to clear nodes queue and nodes set.
   */
  template< typename TSpec >
      inline void
    clear( Path< TSpec >& path )
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
  template< typename TSpec >
      inline void
    reserve( Path< TSpec >& path, typename Path< TSpec >::size_type size )
    {
      path.nodes.reserve( size );
      path.nodes_set.reserve( size );
    }  /* -----  end of template function reserve  ----- */

  /**
   *  @brief  Reserve memory for the path.
   *
   *  @param  path The path.
   *  @param  size The target capacity.
   *
   *  Internal function to reserve memory for nodes and nodes set.
   *
   *  @overload for Dynamic Path.
   */
  template< >
      inline void
    reserve( Path< Dynamic >& path, Path< Dynamic >::size_type size )
    {
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
  template< typename TSpec >
      inline typename Path< TSpec >::size_type
    length( const Path< TSpec >& path )
    {
      return path.get_nodes().size();
    }  /* -----  end of template function length  ----- */

  /**
   * @brief  pop the last node from the path.
   */
  template< typename TSpec >
      inline void
    pop_back( Path< TSpec >& path )
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
  template< >
      inline void
    pop_front( Path< Dynamic >& path )
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
   *  @brief  setter function for nodes (overloaded for `std::vector`).
   */
    inline void
  set_nodes( Path< Dynamic >& path, const std::vector< VarGraph::nodeid_type >& value )
  {
    Path< Dynamic >::nodes_type d;
    std::copy( value.begin(), value.end(), std::back_inserter( d ) );
    path.set_nodes( std::move( d ) );
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
  template< typename TSpec >
      inline void
    trim_back( Path< TSpec >& path, VarGraph::nodeid_type node_id=0 )
    {
      bool found = false;
      while ( !found && length( path ) != 0 ) {
        auto&& last_node = path.get_nodes().back();
        if ( node_id == 0 || last_node == node_id ) found = true;
        pop_back( path );
      }
    }  /* -----  end of template function trim_back  ----- */

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
    inline void
  trim_front( Path< Dynamic >& path, VarGraph::nodeid_type node_id=0 )
  {
    bool found = false;
    while ( !found && length( path ) != 0 ) {
      auto&& first_node = path.get_nodes().front();
      if ( node_id == 0 || first_node == node_id ) found = true;
      pop_front( path );
    }
  }  /* -----  end of template function trim_front  ----- */

  /* END OF Normal Path interface functions  ----------------------------------- */

  /**
   *  @brief  Path template class Compact specialization.
   *
   *  Represent a path in the variation graph with efficient node ID query.
   */
  template< >
    class Path< Compact >
    {
      public:
        /* ====================  DATA MEMBERS  ======================================= */
        std::unordered_set< VarGraph::nodeid_type > nodes_set;
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename decltype( nodes_set )::size_type size_type;
        /* ====================  LIFECYCLE     ======================================= */
        Path( ) = default;
        Path( const std::vector< VarGraph::nodeid_type >& p )
        {
          this->set_nodes( p );
        }
        /* ====================  MUTATORS      ======================================= */
        /**
         *  @brief  Set the nodes in the path (clear the previous state).
         */
          inline void
        set_nodes( const std::vector< VarGraph::nodeid_type >& value )
        {
          this->nodes_set.clear();
          this->nodes_set.reserve( value.size() );
          std::copy( value.begin(), value.end(),
              std::inserter( this->nodes_set, this->nodes_set.end() ) );
        }  /* -----  end of method set_nodes  ----- */
        /* ====================  INTERFACE FUNCTIONS  ================================ */
          friend void
        save< Compact >( const Path< Compact >& path, std::ostream& out );
          friend void
        add_node< Compact >( Path< Compact >& path,
            const VarGraph::nodeid_type& node_id );
          friend void
        clear< Compact >( Path< Compact >& path );
          friend void
        reserve< Compact >( Path< Compact >& path, size_type size );
          friend size_type
        length< Compact >( const Path< Compact >& path );
          friend bool
        contains< Compact >( const Path< Compact >& path,
            VarGraph::nodeid_type node_id );
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
  template< >
      inline void
    save( const Path< Compact >& path, std::ostream& out )
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
  template< >
      inline void
    save( const Path< Compact >& path, const std::string& file_name )
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
  template< >
      inline void
    add_node( Path< Compact >& path, const VarGraph::nodeid_type& node_id )
    {
      path.nodes_set.insert( node_id );
    }

  /**
   *  @brief  Clear the path.
   *
   *  @param  path The path.
   */
  template< >
      inline void
    clear( Path< Compact >& path )
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
  template< >
      inline void
    reserve( Path< Compact >& path, Path< Compact >::size_type size )
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
  template< >
      inline Path< Compact >::size_type
    length( const Path< Compact >& path )
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
  template< typename TSpec >
      inline void
    load( Path< Compact >& path, std::istream& in )
    {
      std::vector< VarGraph::nodeid_type > nodes;
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
  template< typename TSpec >
      inline bool
    contains( const Path< TSpec >& path, VarGraph::nodeid_type node_id )
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
  template< typename TSpec, typename TIter >
      inline bool
    contains( const Path< TSpec >& path, TIter begin, TIter end )
    {
      if ( end - begin > 0 ) {
        auto&& on_path = [&path]( const VarGraph::nodeid_type& i ) {
          return contains( path, i );
        };

        if ( std::all_of( begin, end, on_path ) ) return true;
      }

      return false;
    }  /* -----  end of template function contains  ----- */

  /**
   *  @brief  Check whether this path contains another path.
   *
   *  @param  path The path.
   *  @param  path_nodes The node IDs of the given path to be checked.
   *  @return `true` if this path is a superset of the given path.
   *
   *  Overloaded. See `contains( const Path< TSpec >&, TIter, TIter )`.
   */
  template< typename TSpec, typename TContainer >
      inline bool
    contains( const Path< TSpec >& path, const TContainer& path_nodes )
    {
      return contains( path, path_nodes.begin(), path_nodes.end() );
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
  template< class Iter1, class Iter2 >
      inline bool
    covered_by( Iter1 begin, Iter1 end, Iter2 paths_set_begin, Iter2 paths_set_end )
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
   *  @param  path_nodes The given path to be checked as a container of node IDs
   *  @param  paths_set A set of paths as a container of `Path`.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  Overloaded. See `covered_by( TIter1, TIter1, TIter2, TIter2 )`.
   */
  template< class TContainer1, class TContainer2 >
      inline bool
    covered_by( const TContainer1& path_nodes, const TContainer2& paths_set )
    {
      return covered_by( path_nodes.begin(), path_nodes.end(),
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
  template< typename TContainer >
      inline bool
    covered_by( VarGraph::nodeid_type node_id, const TContainer& paths_set )
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
  template< typename TContainer >
      inline std::size_t
    get_path_coverage( const VarGraph::nodeid_type& node_id, const TContainer& paths_set )
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
  get_random_adjacent ( const VarGraph &vargraph, VarGraph::nodeid_type node_id )
  {
    if ( !vargraph.has_edges_from( node_id ) ) {
      return 0;
    }

    auto fwd_edges = vargraph.edges_from(node_id);

    std::random_device rd;  // Will be used to obtain a seed for the random no. engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, fwd_edges.size() - 1);

    return fwd_edges.at(dis(gen)).to();
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
    };

  /**
   *  @brief  Backtracker graph iterator tag.
   *
   *  Specialization of generic graph iterator tag BacktrackerIter for VarGraph.
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
    };  /* ----------  end of struct Backtracker  ---------- */

  /**
   *  @brief  Haplotyper graph iterator tag.
   *
   *  Specialization of generic graph iterator tag HaplotyperIter for VarGraph.
   */
  template< >
    struct GraphIterTraits< VarGraph, Haplotyper > {
      typedef VarGraph::nodeid_type Value;
      typedef std::deque< Value > TContainer;
      /**< @brief Set of visited paths. */
      typedef std::vector< Path< Compact > > TSet;
      typedef TSet::size_type Level;        /**< @brief No. of selected path so far. */
      typedef struct {
        Value start;                        /**< @brief Start node ID. */
        bool end;                           /**< @brief End flag. */
        Path< Compact > current_path;
        unsigned int setback;
      } TState;
    };  /* ----------  end of struct HaplotyperIter  ---------- */

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
      inline GraphIter< VarGraph, BFS >
    begin( const VarGraph& g,
        typename seqan::Value< GraphIter< VarGraph, BFS > >::Type start )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, BFS > >::Type TValue;

      GraphIter< VarGraph, BFS > begin_itr;

      TValue start_node_id;
      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = g.rank_to_id( 1 );
      }

      begin_itr.state.lb_visited_rank = 1;
      if ( g.id_to_rank( start_node_id ) == 1 ) {
        ++begin_itr.state.lb_visited_rank;
      }

      begin_itr.vargraph_ptr = &g;
      begin_itr.visiting_buffer.push_back( std::make_pair( start_node_id, 0 ) );
      begin_itr.visited.insert( std::make_pair( start_node_id, 0 ) );
      begin_itr.itr_value = begin_itr.visiting_buffer.front().first;

      return begin_itr;
    }

  template< >
      inline void
    go_begin( GraphIter< VarGraph, BFS >& it,
        typename seqan::Value< GraphIter< VarGraph, BFS > >::Type start )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, BFS > >::Type TValue;

      TValue start_node_id;
      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = it.vargraph_ptr->rank_to_id( 1 );
      }

      it.state.lb_visited_rank = 1;
      if ( it.vargraph_ptr->id_to_rank( start_node_id ) == 1 ) {
        ++it.state.lb_visited_rank;
      }

      it.visiting_buffer.clear();
      it.visiting_buffer.push_back( std::make_pair( start_node_id, 0 ) );
      it.visited.clear();
      it.visited.insert( std::make_pair( start_node_id, 0 ) );
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
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper >::operator[]  ----- */

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
      inline GraphIter< VarGraph, Backtracker >
    begin( const VarGraph& g,
        typename seqan::Value< GraphIter< VarGraph, Backtracker > >::Type start )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Backtracker > >::Type TValue;

      GraphIter< VarGraph, Backtracker > begin_itr;
      TValue start_node_id;

      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = g.rank_to_id( 1 );
      }

      begin_itr.vargraph_ptr = &g;
      begin_itr.itr_value = start_node_id;
      begin_itr.state.buffer = 0;
      begin_itr.state.end = false;
      begin_itr.state.start = start_node_id;

      return begin_itr;
    }  /* -----  end of template function begin  ----- */

  template< >
      inline void
    go_begin( GraphIter< VarGraph, Backtracker >& it,
        typename seqan::Value< GraphIter< VarGraph, Backtracker > >::Type start )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Backtracker > >::Type TValue;

      TValue start_node_id;

      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = it.state.start;  // Re-use stored start node.
      }

      it.itr_value = start_node_id;
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

      if ( this->state.end && this->raise_on_end ) {
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
          this->itr_value = this->vargraph_ptr->edges_from( cnode_id ).at( 0 ).to();
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
    at_end( GraphIter< VarGraph, Haplotyper >& it )
    {
      return it.state.end;
    }  /* -----  end of template function at_end  ----- */

  template< >
      inline GraphIter< VarGraph, Haplotyper >
    begin( const VarGraph& g,
        typename seqan::Value< GraphIter< VarGraph, Haplotyper > >::Type start )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Haplotyper > >::Type TValue;

      GraphIter< VarGraph, Haplotyper > begin_itr;

      TValue start_node_id;
      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = g.rank_to_id( 1 );
      }

      begin_itr.vargraph_ptr = &g;
      begin_itr.itr_value = start_node_id;
      begin_itr.state.start = start_node_id;
      begin_itr.state.end = false;
      add_node( begin_itr.state.current_path, begin_itr.itr_value );
      begin_itr.state.setback = 0;

      return begin_itr;
    }  /* -----  end of template function begin  ----- */

  template< >
      inline void
    go_begin( GraphIter< VarGraph, Haplotyper >& it,
        typename seqan::Value< GraphIter< VarGraph, Haplotyper > >::Type start )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Haplotyper > >::Type TValue;

      TValue start_node_id;
      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = it.state.start;  // Re-use start node.
      }

      it.itr_value = start_node_id;
      it.state.start = start_node_id;
      it.visiting_buffer.clear();
      it.state.end = false;  // Re-set at-end flag.
      it.visited.clear();
      clear( it.state.current_path );
      add_node( it.state.current_path, it.itr_value );
      it.state.setback = 0;
    }  /* -----  end of template function go_begin  ----- */

  template< >
      inline typename seqan::Level< GraphIter< VarGraph, Haplotyper > >::Type
    level( GraphIter< VarGraph, Haplotyper >& it )
    {
      return it.visited.size();
    }

  /* Member functions specialization. */

  template< >
      inline void
    GraphIter< VarGraph, Haplotyper >::set_setback( )
    {
      this->state.setback = (( this->visited.size() == 0 /* first path */||
                               this->visited.size() % 2 /* odd */) ?
                             this->visited.size() : this->visited.size() + 1 );
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper >::set_setback  ----- */

  /**
   *  A setback path is a sequence of last 's' nodes of currently generating haplotype.
   *  An unvisited setback path is a path that does not occur as a subset of any
   *  previously generated haplotypes. This function search in adjacent nodes set for
   *  a node that together with 's-1' previously selected nodes forms an unvisited
   *  setback in order to cover more k-mers from all paths in the graph; i.e. generating
   *  more diverse haplotypes.
   */
  template< >
      inline GraphIter< VarGraph, Haplotyper >&
    GraphIter< VarGraph, Haplotyper >::operator++( )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Haplotyper > >::Type TValue;

      if ( this->state.end && this->raise_on_end ) {
        throw std::range_error( "The iterator has reached the end." );
      }

      if ( !this->vargraph_ptr->has_edges_from( this->itr_value ) ) {
        this->state.end = true;
        return *this;
      }

      if ( this->state.setback != 0 &&
          this->visiting_buffer.size() >= this->state.setback ) {
        this->visiting_buffer.pop_front();
      }

      TValue next_candidate = 0;
      auto fwd_edges = this->vargraph_ptr->edges_from( this->itr_value );
      if ( this->state.setback == 0 || fwd_edges.size() == 1 ) {
        next_candidate = fwd_edges[0].to();
      }
      else {
        // Search for a forward node such that the setback path is not in previous paths.
        for ( auto e : fwd_edges ) {
          this->visiting_buffer.push_back( e.to() );
          if ( covered_by( this->visiting_buffer, this->visited ) ) {  // Visited?
            this->visiting_buffer.pop_back();
            continue;                             // Next edge.
          }
          this->visiting_buffer.pop_back();       // No change to the iterator state.
          next_candidate = e.to();                // Found!
        }
      }
      // If no unvisited setback path found, use a node with least path coverage.
      if ( next_candidate == 0 ) {
        next_candidate = least_covered_adjacent( *this->vargraph_ptr,
            this->itr_value, this->visited );
      }
      // If all forward edges are visited, pick one randomly with uniform distribution.
      if ( next_candidate == 0 ) {
        next_candidate =
          get_random_adjacent( ( *this->vargraph_ptr ),  this->itr_value );
      }

      this->itr_value = next_candidate;
      if ( this->state.setback != 0 ) {
        this->visiting_buffer.push_back( this->itr_value );
      }
      add_node( this->state.current_path, this->itr_value );

      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper >::operator++  ----- */

  template< >
      inline GraphIter< VarGraph, Haplotyper >&
    GraphIter< VarGraph, Haplotyper >::operator--( int )
    {
      this->itr_value = this->state.start;    // Reset the iterator to the start node.
      this->visiting_buffer.clear();
      if ( this->state.setback != 0 ) {
        this->visiting_buffer.push_back( this->itr_value );
      }
      this->state.end = false;                // Reset at-end flag.
      clear( this->state.current_path );
      add_node( this->state.current_path, this->itr_value );
      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper >::operator--  ----- */

  template< >
      inline GraphIter< VarGraph, Haplotyper >&
    GraphIter< VarGraph, Haplotyper >::operator--( )
    {
      this->visited.push_back( this->state.current_path );
      this->set_setback();
      (*this)--;
      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper >::operator--  ----- */

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
    GraphIter< VarGraph, Haplotyper >::operator[]( const TContainer& path )
    {
      return covered_by( path, this->visited );
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper >::operator[]  ----- */

  /* END OF Haplotyper template specialization  ---------------------------------- */

  /* Haplotyper iterator interface function  ------------------------------------- */

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
    void
  get_uniq_haplotype( Path<>& haplotype,
      typename seqan::Iterator< VarGraph, Haplotyper >::Type& iter,
      int tries=0 )
  {
    do {
      clear( haplotype );
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

  /* END OF Haplotyper iterator interface function  ------------------------------ */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef VARGRAPH_H__  ----- */

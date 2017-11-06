/*
 * =====================================================================================
 *
 * Filename: vargraph.h
 *
 * Created: Fri Nov 11, 2016  01:08
 *
 * Description: VarGraph class definition.
 *
 * Copyright (c) 2016, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari, <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#ifndef VARGRAPH_H__
#define VARGRAPH_H__

#include <algorithm>
#include <fstream>
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

#include "graph_iter.h"
#include "utils.h"


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
  struct CompactStrategy;
  typedef seqan::Tag< CompactStrategy > Compact;

  /**
   *  @brief  Path template class.
   *
   *  Represent a path in the variation graph with efficient node ID query.
   */
  template< typename TSpec = void >
    class Path
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef std::string string_type;
        typedef string_type::size_type seqsize_type;
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        const VarGraph* vargraph;
        std::deque< VarGraph::nodeid_type > nodes;
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

        Path( const VarGraph* g, std::deque< VarGraph::nodeid_type >&& p )
          : Path( g )
        {
          this->set_nodes( std::move( p ) );
        }

        Path( const VarGraph* g, const std::deque< VarGraph::nodeid_type >& p )
          : Path( g, std::deque< VarGraph::nodeid_type >( p ) )
        { }

        Path( const VarGraph* g, const std::vector< VarGraph::nodeid_type >& p )
          : Path( g )
        {
          this->set_nodes( p );
        }
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
          inline const std::deque< VarGraph::nodeid_type >&
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
        set_nodes( std::deque< VarGraph::nodeid_type >&& value )
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
        set_nodes( const std::deque< VarGraph::nodeid_type >& value )
        {
          this->set_nodes( std::deque< VarGraph::nodeid_type >( value ) );
        }  /* -----  end of method set_nodes  ----- */

        /**
         *  @brief  setter function for nodes (overloaded for `std::vector`).
         */
          inline void
        set_nodes( const std::vector< VarGraph::nodeid_type >& value )
        {
          std::deque< VarGraph::nodeid_type > d;
          d.reserve( value.size() );
          std::copy( value.begin(), value.end(), std::back_inserter( d ) );
          this->set_nodes( std::move( d ) );
        }
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
      Path< TSpec >::seqsize_type cursor = 0;
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
      std::deque< VarGraph::nodeid_type > nodes;
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
    position_to_id( const Path< TSpec >& path, Path< TSpec >::seqsize_type pos )
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
    position_to_offset( const Path< TSpec >& path, Path< TSpec >::seqsize_type pos )
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

  template< typename TSpec >
      inline void
    pop_front( Path< TSpec >& path )
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
   *  @brief  Trim the path from the node whose ID matches the given ID.
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
    trim( Path< TSpec >& path, VarGraph::nodeid_type node_id )
    {
      bool found = false;
      while ( !found && length( path ) != 0 ) {
        auto&& last_node = path.get_nodes().back();
        if ( node_id == 0 || last_node == node_id ) found = true;
        pop_back( path );
      }
    }  /* -----  end of template function trim  ----- */

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
  /* Haplotyper iterator interface function declaration  ------------------------- */

    void
  get_uniq_haplotype( Path<>& haplotype,
      typename seqan::Iterator< VarGraph, Haplotyper >::Type& iter,
      int tries=0 );

  /* END OF Haplotyper iterator interface function declaration  ------------------ */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef VARGRAPH_H__  ----- */

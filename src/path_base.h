/**
 *    @file  path_base.h
 *   @brief  Path template class definitions.
 *
 *  This header file defines Path template class in variation graph, its specializations
 *  for different purposes.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sat Mar 31, 2018  22:17
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  PATH_BASE_H__
#define  PATH_BASE_H__

#include <fstream>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <string>
#include <set>
#include <deque>
#include <algorithm>
#include <utility>

#include <seqan/basic.h>
#include <sdsl/bit_vectors.hpp>

#include "sequence.h"
#include "utils.h"


namespace grem{
  /* Path specialization tags. */
  struct DefaultStrategy;
  struct DynamicStrategy;
  struct MicroStrategy;
  struct CompactStrategy;
  struct HaplotypeStrategy;
  typedef seqan::Tag< DefaultStrategy > Default;
  typedef seqan::Tag< DynamicStrategy > Dynamic;
  typedef seqan::Tag< CompactStrategy > Compact;
  typedef seqan::Tag< MicroStrategy > Micro;
  typedef seqan::Tag< HaplotypeStrategy > Haplotype;

  /* Path node existence query strategies */
  struct OrderedStrategy;
  struct UnorderedStrategy;
  typedef seqan::Tag< OrderedStrategy > Ordered;
  typedef seqan::Tag< UnorderedStrategy > Unordered;

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

  template< typename TGraph >
    struct PathTraits< TGraph, Compact > {
      typedef sdsl::enc_vector< sdsl::coder::elias_delta > TNodeSequence;
    };

  template< typename TGraph >
    struct PathTraits< TGraph, Haplotype > {
      typedef sdsl::bit_vector TNodeSequence;
    };

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
        typedef TGraph graph_type;
        typedef TSpec spec_type;
        typedef std::string string_type;
        typedef string_type::size_type seqsize_type;
        typedef typename TTraits::TNodeSequence nodes_type;
        typedef typename graph_type::nodeid_type value_type;
        typedef typename nodes_type::size_type size_type;
        typedef ptrdiff_t difference_type;
        typedef typename nodes_type::const_iterator const_iterator;
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        const graph_type* vargraph;
        nodes_type nodes;
        seqsize_type seqlen;
        /* Loaded on demand. */
        string_type seq;
        /* Loaded after calling `initialize`. */
        bool initialized;
        sdsl::bit_vector bv_node_breaks;
        sdsl::bit_vector::rank_1_type rs_node_breaks;
        sdsl::bit_vector::select_1_type ss_node_breaks;
      public:
        /* ====================  LIFECYCLE     ======================================= */
        Path( const graph_type* g )
          : vargraph( g ), seqlen( 0 ), initialized( false )
        { }

        Path( const graph_type* g, nodes_type&& p )
          : Path( g )
        {
          this->set_nodes( std::move( p ) );
        }

        Path( const graph_type* g, const nodes_type& p )
          : Path( g, nodes_type( p ) )
        { }

        Path( const Path& other )
        {
          this->vargraph = other.vargraph;
          this->nodes = other.nodes;
          this->seqlen = other.seqlen;
          this->seq = other.seq;
          this->initialized = other.initialized;
          this->bv_node_breaks = other.bv_node_breaks;
          sdsl::util::init_support( this->rs_node_breaks, &this->bv_node_breaks );
          sdsl::util::init_support( this->ss_node_breaks, &this->bv_node_breaks );
        }

        Path( Path&& other )
        {
          this->vargraph = other.vargraph;
          this->nodes = std::move( other.nodes );
          this->seqlen = other.seqlen;
          this->seq = std::move( other.seq );
          this->initialized = other.initialized;
          this->bv_node_breaks.swap( other.bv_node_breaks );
          sdsl::util::init_support( this->rs_node_breaks, &this->bv_node_breaks );
          sdsl::util::init_support( this->ss_node_breaks, &this->bv_node_breaks );
          // Free up memory first of `other`.
          sdsl::util::clear( other.rs_node_breaks );
          sdsl::util::clear( other.ss_node_breaks );
        }

        Path& operator=( const Path& other )
        {
          this->vargraph = other.vargraph;
          this->nodes = other.nodes;
          this->seqlen = other.seqlen;
          this->seq = other.seq;
          this->initialized = other.initialized;
          this->bv_node_breaks = other.bv_node_breaks;
          sdsl::util::init_support( this->rs_node_breaks, &this->bv_node_breaks );
          sdsl::util::init_support( this->ss_node_breaks, &this->bv_node_breaks );
          return *this;
        }

        Path& operator=( Path&& other )
        {
          this->vargraph = other.vargraph;
          this->nodes = std::move( other.nodes );
          this->seqlen = other.seqlen;
          this->seq = std::move( other.seq );
          this->initialized = other.initialized;
          // Free up memory of `this` first.
          sdsl::util::clear( this->bv_node_breaks );
          sdsl::util::clear( this->rs_node_breaks );
          sdsl::util::clear( this->ss_node_breaks );

          this->bv_node_breaks.swap( other.bv_node_breaks );
          sdsl::util::init_support( this->rs_node_breaks, &this->bv_node_breaks );
          sdsl::util::init_support( this->ss_node_breaks, &this->bv_node_breaks );
          // Free up memory first of `other`.
          sdsl::util::clear( other.rs_node_breaks );
          sdsl::util::clear( other.ss_node_breaks );
          return *this;
        }

        ~Path() = default;
        /* ====================  OPERATORS     ======================================= */
        template< typename TSpec2 >
            inline Path& operator=( Path< graph_type, TSpec2 >&& other );
        template< typename TSpec2 >
            inline Path& operator=( const Path< graph_type, TSpec2 >& other );
        /* ====================  ACCESSORS     ======================================= */
        /**
         *  @brief  getter function for vargraph.
         */
          inline const graph_type*
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

          inline value_type
        operator[]( size_t idx ) const
        {
          return this->get_nodes()[ idx ];
        }

          inline value_type
        at( size_t idx ) const
        {
          return this->get_nodes().at( idx );
        }

          inline size_type
        size( ) const
        {
          return this->get_nodes().size();
        }

          inline bool
        empty( ) const
        {
          return this->size() == 0;
        }

          inline const_iterator
        begin( ) const
        {
          return this->get_nodes().begin();
        }

          inline const_iterator
        end( ) const
        {
          return this->get_nodes().end();
        }

          inline value_type
        front( ) const
        {
          return *this->begin();
        }

          inline value_type
        back( ) const
        {
          return *( this->end() - 1 );
        }
        /* ====================  MUTATORS      ======================================= */
        /**
         *  @brief  setter function for vargraph.
         */
          inline void
        set_vargraph( const graph_type* value )
        {
          this->vargraph = value;
        }  /* -----  end of method set_vargraph  ----- */

          inline void
        set_nodes( nodes_type&& value );

        /**
         *  @brief  setter function for nodes.
         */
          inline void
        set_nodes( const nodes_type& value )
        {
          this->set_nodes( nodes_type( value ) );
        }  /* -----  end of method set_nodes  ----- */

        template< typename TIter >
            inline void
          set_nodes( TIter begin, TIter end )
          {
            nodes_type nd( end - begin );
            std::copy( begin, end, nd.begin() );
            this->set_nodes( std::move( nd ) );
          }
        /* ====================  METHODS       ======================================= */
        /**
         *  @brief  Initialize data structures for efficient rank and select queries.
         *
         *  @param  path The path.
         *
         *  It constructs node breaks bit vector and corresponding rank and select supports.
         */
          inline void
        initialize( )
        {
          if ( this->is_initialized() || this->size() == 0 ) return;

          assert( this->vargraph != nullptr );

          this->init_bv_node_breaks();
          sdsl::util::init_support( this->rs_node_breaks, &this->bv_node_breaks );
          sdsl::util::init_support( this->ss_node_breaks, &this->bv_node_breaks );
          this->initialized = true;
        }

          inline void
        push_back( value_type const& nid )
        {
          assert( this->vargraph != nullptr );
          this->nodes.push_back( nid );
          this->seqlen += this->vargraph->node_length( nid );
          if ( this->seq.length() != 0 ) {
            this->seq += this->vargraph->node_sequence( nid );
          }
          this->initialized = false;
        }

          inline void
        pop_back( );

          inline void
        pop_front( );

          inline void
        clear( )
        {
          grem::clear( this->nodes );
          this->seqlen = 0;
          this->seq.clear();
          sdsl::util::clear( this->bv_node_breaks );
          sdsl::util::clear( this->rs_node_breaks );
          sdsl::util::clear( this->ss_node_breaks );
          this->initialized = false;
        }

          inline void
        reserve( size_type size )
        {
          grem::reserve( this->nodes, size );  /* defined in `utils.h` */
        }

          inline void
        serialize( std::ostream& out ) const
        {
          if ( ! this->is_initialized() )
            throw std::runtime_error( "Path must be initialized to be serialized" );
          grem::serialize( out, this->nodes );
          this->bv_node_breaks.serialize( out );
        }

          inline void
        serialize( std::ostream& out )
        {
          this->initialize();
          const Path* const_this = this;
          const_this->serialize( out );
        }

          inline void
        load( std::istream& in )
        {
          nodes_type tmp;
          deserialize( in, tmp );
          this->set_nodes( std::move( tmp ) );
          this->bv_node_breaks.load( in );
          sdsl::util::init_support( this->rs_node_breaks, &this->bv_node_breaks );
          sdsl::util::init_support( this->ss_node_breaks, &this->bv_node_breaks );
          this->initialized = true;
        }

        /**
         *  @brief  Get the node rank (0-based) in the path by a position in its sequence.
         *
         *  @param  pos The position in the path (0-based).
         *  @return The node rank in the nodes queue (0-based) on whose label the position
         *          `pos` in the path sequence relies.
         *
         *  The value of `rank(pos)` is the node rank in the path's nodes queue.
         */
          inline size_type
        rank( seqsize_type pos ) const
        {
          assert( this->is_initialized() );
          if ( pos < 0 || pos >= this->seqlen ) {
            throw std::runtime_error( "Position out of range." );
          }
          return this->rs_node_breaks( pos );
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
          inline seqsize_type
        select( size_type rank ) const
        {
          assert( this->is_initialized() );
          if ( rank < 0 || rank >= this->size() ) {
            throw std::runtime_error( "Rank out of range." );
          }
          if ( rank == 0 ) return 0;
          return this->ss_node_breaks( rank ) + 1;
        }

          inline bool
        contains( value_type nid ) const
        {
          if ( std::find( this->begin(), this->end(), nid ) != this->end() ) return true;
          return false;
        }
        /* ====================  FRIENDSHIPS   ======================================= */
        friend class Path< graph_type, Default >;
        friend class Path< graph_type, Dynamic >;
        friend class Path< graph_type, Compact >;
        friend class Path< graph_type, Haplotype >;
      private:
        /**
         *  @brief  Initialize node breaks bitvector.
         *
         *  @param  path The path.
         *
         *  The node breaks bit vector is a bit vector of sequence length. A bit at position
         *  `i` is set if a node ends at that position in the sequence. For example for this
         *  path:
         *
         *  (GCAAT) -> (A) -> (TTAGCC) -> (GCA)
         *
         *  the corresponding path is:
         *
         *  GCAATATTAGCCGCA
         *
         *  and its bit vector is:
         *
         *  000011000001001
         *
         *  which has a set bit at the first position of each node in the path.
         */
          inline void
        init_bv_node_breaks( )
        {
          assert( this->size() != 0 );
          sdsl::util::assign( this->bv_node_breaks, sdsl::bit_vector( this->seqlen, 0 ) );
          seqsize_type cursor = 0;
          for ( auto it = this->begin(); it != this->end(); ++it ) {
            cursor += this->vargraph->node_length( *it );
            this->bv_node_breaks[ cursor - 1 ] = 1;
          }
        }
    };  /* -----  end of template class Path  ----- */

  /* Path member functions  ---------------------------------------------------- */

  /**
   *  @brief  Move assignment from other type of path.
   *
   *  @param  other The other path.
   *  @return The reference to this instance after move.
   */
  template< typename TGraph, typename TSpec1 >
    template< typename TSpec2 >
        inline Path< TGraph, TSpec1 >&
      Path< TGraph, TSpec1 >::operator=( Path< TGraph, TSpec2 >&& other )
      {
        using namespace sdsl::util;

        this->vargraph = other.vargraph;

        assign( this->nodes, other.nodes );
        grem::clear( other.nodes );

        this->seqlen = other.seqlen;
        this->seq = std::move( other.seq );
        this->initialized = other.initialized;
        sdsl::util::clear( this->bv_node_breaks );
        this->bv_node_breaks.swap( other.bv_node_breaks );
        init_support( this->rs_node_breaks, &this->bv_node_breaks );
        init_support( this->ss_node_breaks, &this->bv_node_breaks );
        sdsl::util::clear( other.rs_node_breaks );
        sdsl::util::clear( other.ss_node_breaks );
        return *this;
      }

  /**
   *  @brief  Assign a path with another type of path.
   *
   *  @param  other The Dynamic path.
   *  @return The reference to this instance after assignment.
   */
  template< typename TGraph, typename TSpec1 >
    template< typename TSpec2 >
        inline Path< TGraph, TSpec1 >&
      Path< TGraph, TSpec1 >::operator=( const Path< TGraph, TSpec2 >& other )
      {
        using namespace sdsl::util;

        this->vargraph = other.vargraph;
        assign( this->nodes, other.nodes );
        this->seqlen = other.seqlen;
        this->seq = other.seq;
        this->initialized = other.initialized;
        this->bv_node_breaks = other.bv_node_breaks;
        init_support( this->rs_node_breaks, &this->bv_node_breaks );
        init_support( this->ss_node_breaks, &this->bv_node_breaks );
        return *this;
      }

  /**
   *  @brief  setter function for nodes.
   */
  template< typename TGraph, typename TSpec >
      inline void
    Path< TGraph, TSpec >::set_nodes( nodes_type&& value )
    {
      this->clear();
      if ( value.empty() ) return;

      assert( this->vargraph != nullptr );

      this->nodes = std::move( value );
      for ( auto&& n : this->nodes ) this->seqlen += this->vargraph->node_length( n );
    }  /* -----  end of method set_nodes  ----- */

  template< typename TGraph, typename TSpec >
      inline void
    Path< TGraph, TSpec >::pop_back( )
    {
      assert( this->vargraph != nullptr );

      if ( this->nodes.empty() ) return;

      auto&& last_node_len = this->vargraph->node_length( this->nodes.back() );
      assert( this->seqlen >= last_node_len );
      this->seqlen -= last_node_len;
      this->nodes.pop_back();
      this->initialized = false;

      if ( this->seq.length() != 0 ) this->seq.resize( this->seqlen );
    }

  template< typename TGraph, typename TSpec >
      inline void
    Path< TGraph, TSpec >::pop_front( )
    {
      assert( this->vargraph != nullptr );

      if ( this->nodes.empty() ) return;

      auto&& first_node_len = this->vargraph->node_length( this->nodes.front() );
      assert( this->seqlen >= first_node_len );
      this->seqlen -= first_node_len;
      this->nodes.pop_front();
      this->initialized = false;

      if ( this->seq.length() != 0 ) this->seq = this->seq.substr( first_node_len );
    }

  /* END OF Path member functions  --------------------------------------------- */

  /**
   *  @brief  Path template class Micro specialization.
   *
   *  Represent a path in the variation graph with efficient node ID query.
   */
  template< typename TGraph >
    class Path< TGraph, Micro >
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TGraph graph_type;
        typedef Micro spec_type;
        typedef std::set< typename TGraph::nodeid_type > nodes_set_type;
        typedef typename graph_type::nodeid_type value_type;
        typedef typename nodes_set_type::size_type size_type;
        typedef typename nodes_set_type::difference_type difference_type;
        typedef typename nodes_set_type::const_iterator const_iterator;
        /* ====================  LIFECYCLE     ======================================= */
        Path( ) = default;
        Path( const std::vector< typename graph_type::nodeid_type >& p )
        {
          this->set_nodes( p );
        }

        template< typename TSpec >
          Path( const Path< graph_type, TSpec >& other )
          {
            this->set_nodes( other.get_nodes() );
          }

        template< typename TSpec >
          Path& operator=( const Path< graph_type, TSpec >& other )
          {
            this->set_nodes( other.get_nodes() );
          }

        Path( const Path& ) = default;
        Path( Path&& ) = default;
        Path& operator=( const Path& ) = default;
        Path& operator=( Path&& ) = default;
        ~Path() = default;
        /* ====================  ACCESSORS     ======================================= */
        /**
         *  @brief  getter function for nodes_set.
         */
          inline const nodes_set_type&
        get_nodes_set( ) const
        {
          return this->nodes_set;
        }  /* -----  end of method get_nodes_set  ----- */

          inline size_type
        size( ) const
        {
          return this->nodes_set.size();
        }

          inline bool
        empty( ) const
        {
          return this->size() == 0;
        }

          inline const_iterator
        begin( ) const
        {
          return nodes_set.begin();
        }

          inline const_iterator
        end( ) const
        {
          return nodes_set.end();
        }
        /* ====================  MUTATORS      ======================================= */
        /**
         *  @brief  Set the nodes in the path (clear the previous state).
         */
          inline void
        set_nodes( const std::vector< value_type >& value )
        {
          this->clear( );
          grem::reserve( this->nodes_set, value.size() );
          std::copy( value.begin(), value.end(),
              std::inserter( this->nodes_set, this->nodes_set.end() ) );
        }  /* -----  end of method set_nodes  ----- */
        /* ====================  METHODS       ======================================= */
          inline void
        push_back( value_type const& node_id )
        {
          this->nodes_set.insert( node_id );
        }

          inline void
        clear( )
        {
          this->nodes_set.clear();
        }

          inline void
        reserve( size_type size )
        {
          grem::reserve( this->nodes_set, size );
        }

          inline void
        serialize( std::ostream& out ) const
        {
          grem::serialize( out, this->nodes_set, this->nodes_set.begin(), this->nodes_set.end() );
        }

          inline void
        load( std::istream& in )
        {
          deserialize( in, this->nodes_set, std::inserter( this->nodes_set, this->nodes_set.end() ) );
        }

          inline bool
        contains( value_type nid ) const
        {
          return this->nodes_set.find( nid ) != this->nodes_set.end();
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        nodes_set_type nodes_set;
    };  /* -----  end of specialized template class Path  ----- */

  /**
   *  @brief  Path specialized template class for representing a Haplotype in DAG graphs.
   *
   *  NOTE: Apart from assuming that the underlying variation graph is DAG, it assumes
   *        that the node ranks are a topological sort of nodes.
   */
  template< typename TGraph >
    class Path< TGraph, Haplotype > {
      private:
        /* ====================  TYPEDEFS      ======================================= */
        typedef PathTraits< TGraph, Haplotype > TTraits;
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TGraph graph_type;
        typedef Haplotype spec_type;
        typedef typename TTraits::TNodeSequence nodes_type;
        typedef typename graph_type::nodeid_type value_type;
        typedef typename nodes_type::size_type size_type;
        typedef ptrdiff_t difference_type;
        typedef sdsl::random_access_const_iterator< Path > const_iterator;
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        const graph_type* vargraph;
        nodes_type nodes;
        typename graph_type::rank_type last_node_rank;
        bool initialized;
        typename nodes_type::rank_1_type rs_nodes;
        typename nodes_type::select_1_type ss_nodes;
      public:
        /* ====================  LIFECYCLE     ======================================= */
        Path( const graph_type* g ) :
          vargraph( g ),
          nodes( g->max_node_rank(), 0 ),
          last_node_rank( 0 ),
          initialized( false )
        { }

        template< typename TContainer >
          Path( const graph_type* g, const TContainer& node_ids )
            : Path( g )
          {
            this->set_nodes( node_ids );
          }

        Path( const Path& other )
        {
          this->vargraph = other.vargraph;
          this->nodes = other.nodes;
          this->last_node_rank = other.last_node_rank;
          this->initialize();
        }

        Path( Path&& other )
        {
          this->vargraph = other.vargraph;
          this->nodes = std::move( other.nodes );
          this->last_node_rank = other.last_node_rank;
          this->initialize();
          sdsl::util::clear( other.rs_nodes );
          sdsl::util::clear( other.ss_nodes );
        }

          Path&
        operator=( const Path& other )
        {
          this->vargraph = other.vargraph;
          this->nodes = other.nodes;
          this->last_node_rank = other.last_node_rank;
          this->initialize();
          return *this;
        }

          Path&
        operator=( Path&& other )
        {
          this->vargraph = other.vargraph;
          this->nodes = std::move( other.nodes );
          this->last_node_rank = other.last_node_rank;
          this->initialize();
          sdsl::util::clear( other.rs_nodes );
          sdsl::util::clear( other.ss_nodes );
          return *this;
        }

        ~Path() = default;
        /* ====================  OPERATORS     ======================================= */
        template< typename TSpec2 >
            inline Path&
          operator=( const Path< graph_type, TSpec2 >& other )
          {
            this->vargraph = other.vargraph;
            this->set_nodes( other.nodes );
            return *this;
          }
        /* ====================  ACCESSORS     ======================================= */
        /**
         *  @brief  getter function for vargraph.
         */
          inline const graph_type*
        get_vargraph( ) const
        {
          return this->vargraph;
        }  /* -----  end of method get_vargraph  ----- */

        /**
         *  @brief  getter function for nodes.
         */
          inline const Path&
        get_nodes( ) const
        {
          return *this;
        }  /* -----  end of method get_nodes  ----- */

        /**
         *  @brief  Is path initialized?
         *
         *  @return `true` if the path is initialized; `false` otherwise.
         *
         *  This function initializes nodes rank and select support data structures.
         */
          inline bool
        is_initialized( ) const
        {
          return this->initialized;
        }

          inline value_type
        operator[]( size_t idx ) const
        {
          assert( this->is_initialized() );
          return this->vargraph->rank_to_id( this->ss_nodes( idx + 1 ) + 1 );
        }

          inline value_type
        at( size_t idx ) const
        {
          if ( ! this->is_initialized() )
            throw std::runtime_error( "Path must be initialized to be operational" );
          if ( idx < 0 || idx >= this->size() )
            throw std::runtime_error( "Index out of range" );
          return (*this)[ idx ];
        }

          inline size_type
        size( ) const
        {
          assert( this->is_initialized() );
          return this->rs_nodes( this->nodes.size() );
        }

          inline bool
        empty( ) const
        {
          return this->size() == 0;
        }

          inline const_iterator
        begin( ) const
        {
          return const_iterator( this, 0 );
        }

          inline const_iterator
        end( ) const
        {
          return const_iterator( this, this->size() );
        }

          inline value_type
        front( ) const
        {
          return *this->begin();
        }

          inline value_type
        back( ) const
        {
          return *( this->end() - 1 );
        }
        /* ====================  MUTATORS      ======================================= */
        /**
         *  @brief  setter function for vargraph.
         */
          inline void
        set_vargraph( const graph_type* value )
        {
          this->vargraph = value;
        }  /* -----  end of method set_vargraph  ----- */

        template< typename TIter >
            inline void
          set_nodes( TIter begin, TIter end )
          {
            this->clear();
            for ( ; begin != end; ++begin ) this->push_back( *begin );
            this->initialize();
          }

        template< typename TContainer >
            inline void
          set_nodes( const TContainer& node_ids )
          {
            this->set_nodes( node_ids.begin(), node_ids.end() );
          }  /* -----  end of template method set_nodes  ----- */
        /* ====================  METHODS       ======================================= */
          inline void
        push_back( value_type const& nid )
        {
          assert( this->vargraph != nullptr );

          auto nrank = this->vargraph->id_to_rank( nid );
          if ( nrank <= this->last_node_rank )
            throw std::runtime_error( "Path IDs sequence must be non-decreasing" );

          this->nodes[ nrank - 1 ] = 1;
          this->initialized = false;
          this->last_node_rank = nrank;
        }

        /**
         *  @brief  Remove the last node from the path.
         *
         *  XXX: This method requires re-initialization of the path which is costly.
         */
          inline void
        pop_back( )
        {
          if ( this->empty() ) return;
          if ( this->size() == 1 ) this->last_node_rank = 0;
          else this->last_node_rank = this->ss_nodes( this->size() - 1 ) + 1;
          this->nodes[ this->ss_nodes( this->size() ) ] = 0;
          this->initialize();
        }


        /**
         *  @brief  Remove the first node from the path.
         *
         *  XXX: This method requires re-initialization of the path which is costly.
         */
          inline void
        pop_front( )
        {
          if ( this->empty() ) return;
          if ( this->size() == 1 ) this->last_node_rank = 0;
          this->nodes[ this->ss_nodes( 1 ) ] = 0;
          this->initialize();
        }

        /**
         *  @brief  Initialize data structures for efficient rank and select queries.
         *
         *  @param  path The path.
         *
         *  It initializes the rank and select supports for nodes bit vector.
         *
         *  @overload for Haplotype paths.
         */
          inline void
        initialize( )
        {
          sdsl::util::init_support( this->rs_nodes, &this->nodes );
          sdsl::util::init_support( this->ss_nodes, &this->nodes );
          this->initialized = true;
        }

          inline void
        serialize( std::ostream& out ) const
        {
          this->nodes.serialize( out );
        }

          inline void
        load( std::istream& in )
        {
          this->nodes.load( in );
          this->initialize();
          this->last_node_rank = this->ss_nodes( this->get_nodes().size() ) + 1;
        }

          inline void
        clear( )
        {
          sdsl::util::assign( this->nodes,
              sdsl::bit_vector( this->vargraph->max_node_rank(), 0 ) );
          this->initialize();
          this->last_node_rank = 0;
        }  /* -----  end of method clear  ----- */

          inline bool
        contains( value_type nid ) const
        {
          if ( ! this->vargraph->has_node( nid ) ) return false;

          auto nrank = this->vargraph->id_to_rank( nid );
          return this->nodes[ nrank - 1 ] == 1;
        }

        template< typename TIter >
            inline bool
          contains( TIter begin, TIter end ) const
          {
            if ( begin == end ) return false;

            auto brank = this->vargraph->id_to_rank( *begin );
            auto erank = this->vargraph->id_to_rank( *( end - 1 ) );
            if ( erank < brank || brank == 0 || erank == 0 ) return false;

            long int plen = this->rs_nodes( erank ) - this->rs_nodes( brank - 1 );
            long int qlen = end - begin;
            if ( plen != qlen || plen < 0 ) return false;

            typename graph_type::rank_type prev = 0;
            for ( ; begin != end; ++begin ) {
              auto curr = this->vargraph->id_to_rank( *begin );
              if ( curr <= prev || !this->contains( *begin ) ) return false;
              prev = curr;
            }

            return true;
          }

        template< typename TIter >
            inline bool
          rcontains( TIter rbegin, TIter rend ) const
          {
            if ( rbegin == rend ) return false;

            auto rbrank = this->vargraph->id_to_rank( *rbegin );
            auto rerank = this->vargraph->id_to_rank( *( rend - 1 ) );
            if ( rbrank < rerank || rbrank == 0 || rerank == 0 ) return false;

            long int plen = this->rs_nodes( rbrank ) - this->rs_nodes( rerank - 1 );
            long int qlen = rend - rbegin;
            if ( plen != qlen || plen < 0 ) return false;

            auto prev = this->vargraph->id_to_rank( *rbegin ) + 1;
            for ( ; rbegin != rend; ++rbegin ) {
              auto curr = this->vargraph->id_to_rank( *rbegin );
              if ( curr >= prev || !this->contains( *rbegin ) ) return false;
              prev = curr;
            }

            return true;
          }
    };  /* -----  end of specialized template class Path  ----- */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef PATH_BASE_H__  ----- */

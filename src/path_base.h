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
  typedef seqan::Tag< DefaultStrategy > Default;
  typedef seqan::Tag< DynamicStrategy > Dynamic;
  typedef seqan::Tag< CompactStrategy > Compact;
  typedef seqan::Tag< MicroStrategy > Micro;

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

  /* Path interface functions forwards  ---------------------------------------- */

  template< typename TGraph, typename TSpec >
    class Path;

  template< typename TGraph, typename TSpec >
      void
    init_bv_node_breaks( Path< TGraph, TSpec >& path );
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
        typedef TSpec spec_type;
        typedef TGraph graph_type;
        typedef std::string string_type;
        typedef string_type::size_type seqsize_type;
        typedef typename TTraits::TNodeSequence nodes_type;
        typedef typename nodes_type::size_type size_type;
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        const TGraph* vargraph;
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
            inline Path& operator=( Path< TGraph, TSpec2 >&& other );
        template< typename TSpec2 >
            inline Path& operator=( const Path< TGraph, TSpec2 >& other );
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
        /* ====================  FRIENDSHIPS   ======================================= */
        friend class Path< TGraph, Default >;
        friend class Path< TGraph, Dynamic >;
        friend class Path< TGraph, Compact >;
        /* ====================  INTERFACE FUNCTIONS  ================================ */
          friend void
        init_bv_node_breaks< TGraph, TSpec >( Path< TGraph, TSpec >& path );
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
        clear( other.nodes );

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
      clear( *this );
      if ( value.empty() ) return;

      assert( this->vargraph != nullptr );

      this->nodes = std::move( value );
      for ( auto&& n : this->nodes ) this->seqlen += this->vargraph->node_length( n );
    }  /* -----  end of method set_nodes  ----- */

  /* END OF Path member functions  --------------------------------------------- */

  /* Micro Path interface functions forwards  ---------------------------------- */

  template< typename TGraph >
      void
    save( const Path< TGraph, Micro >& path, std::ostream& out );
  template< typename TGraph >
      void
    add_node( Path< TGraph, Micro >& path, typename TGraph::nodeid_type const& node_id );
  template< typename TGraph >
      void
    clear( Path< TGraph, Micro >& path );
  template< typename TGraph >
      void
    reserve( Path< TGraph, Micro >& path,
        typename Path< TGraph, Micro >::size_type size );
  template< typename TGraph >
      typename Path< TGraph, Micro >::size_type
    length( const Path< TGraph, Micro >& path );
  template< typename TGraph >
      bool
    contains( const Path< TGraph, Micro >& path,
        typename TGraph::nodeid_type node_id );

  /* END OF Micro Path interface functions forwards  --------------------------- */

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
        typedef Micro spec_type;
        typedef std::set< typename TGraph::nodeid_type > nodes_set_type;
        typedef typename nodes_set_type::size_type size_type;
        /* ====================  LIFECYCLE     ======================================= */
        Path( ) = default;
        Path( const std::vector< typename TGraph::nodeid_type >& p )
        {
          this->set_nodes( p );
        }

        template< typename TSpec >
          Path( const Path< TGraph, TSpec >& other )
          {
            this->set_nodes( other.get_nodes() );
          }

        template< typename TSpec >
          Path& operator=( const Path< TGraph, TSpec >& other )
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
        /* ====================  MUTATORS      ======================================= */
        /**
         *  @brief  Set the nodes in the path (clear the previous state).
         */
          inline void
        set_nodes( const std::vector< typename TGraph::nodeid_type >& value )
        {
          this->nodes_set.clear();
          grem::reserve( this->nodes_set, value.size() );
          std::copy( value.begin(), value.end(),
              std::inserter( this->nodes_set, this->nodes_set.end() ) );
        }  /* -----  end of method set_nodes  ----- */
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        nodes_set_type nodes_set;
        /* ====================  INTERFACE FUNCTIONS  ================================ */
          friend void
        save< TGraph >( const Path< TGraph, Micro >& path, std::ostream& out );
          friend void
        add_node< TGraph >( Path< TGraph, Micro >& path,
            typename TGraph::nodeid_type const& node_id );
          friend void
        clear< TGraph >( Path< TGraph, Micro >& path );
          friend void
        reserve< TGraph >( Path< TGraph, Micro >& path, size_type size );
          friend size_type
        length< TGraph >( const Path< TGraph, Micro >& path );
          friend bool
        contains< TGraph >( const Path< TGraph, Micro >& path,
            typename TGraph::nodeid_type node_id );
    };  /* -----  end of specialized template class Path  ----- */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef PATH_BASE_H__  ----- */

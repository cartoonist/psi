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
        typedef typename graph_type::offset_type offset_type;
        typedef typename nodes_type::size_type size_type;
        typedef ptrdiff_t difference_type;
        typedef typename nodes_type::const_iterator const_iterator;
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        const graph_type* vargraph;
        nodes_type nodes;
        offset_type left;   /**< @brief The length of the first node sequence */
        offset_type right;  /**< @brief The length of the last node sequence */
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
          : vargraph( g ), left( 0 ), right( 0 ), seqlen( 0 ), initialized( false )
        { }

        /**
         *  XXX NOTE: Be careful of assigning `left` and `right` offset. The zero value
         *            means all sequence of first node or last one is included in the
         *            path. Any value < len(node) indicates the lenght of the suffix or
         *            the prefix of the sequence of first or last node respectively.
         *            They are NOT sequence offsets of the first or last node.
         */
        Path( const graph_type* g, nodes_type p, offset_type l=0, offset_type r=0 )
          : Path( g )
        {
          this->set_nodes( std::move( p ), l, r );
        }

        Path( const Path& other )
        {
          this->vargraph = other.vargraph;
          this->nodes = other.nodes;
          this->left = other.left;
          this->right = other.right;
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
          this->left = other.left;
          this->right = other.right;
          this->seqlen = other.seqlen;
          this->seq = std::move( other.seq );
          this->initialized = other.initialized;
          this->bv_node_breaks.swap( other.bv_node_breaks );
          sdsl::util::init_support( this->rs_node_breaks, &this->bv_node_breaks );
          sdsl::util::init_support( this->ss_node_breaks, &this->bv_node_breaks );
          // Free up memory first of `other`.
          // :TODO:Sun Oct 21 23:33:\@cartoonist: Clear does not free up memory!
          sdsl::util::clear( other.rs_node_breaks );
          sdsl::util::clear( other.ss_node_breaks );
        }

        Path& operator=( const Path& other )
        {
          this->vargraph = other.vargraph;
          this->nodes = other.nodes;
          this->left = other.left;
          this->right = other.right;
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
          this->left = other.left;
          this->right = other.right;
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
         *  @brief  get the node offset of the head node (front node).
         */
          inline offset_type
        get_head_offset( ) const
        {
          if ( this->left == 0 ) return 0;
          assert( !this->empty() );  /* Left offset of an empty path should be always 0 */
          return this->vargraph->node_length( this->front() ) - this->left;
        }

        /**
         *  @brief  getter function for seqlen.
         */
          inline seqsize_type
        get_sequence_len( ) const
        {
          return this->seqlen;
        }  /* -----  end of method get_sequence_len  ----- */

        /**
         *  @brief  get left offset value.
         */
          inline offset_type
        get_left_len( ) const
        {
          assert( !this->empty() );
          return this->left ? this->left : this->vargraph->node_length( this->front() );
        }

        /**
         *  @brief  get right offset value.
         */
          inline offset_type
        get_right_len( ) const
        {
          assert( !this->empty() );
          return this->right ? this->right : this->vargraph->node_length( this->back() );
        }

        /**
         *  @brief  get the sequence length of the head node (front node).
         */
          inline seqsize_type
        get_seqlen_head( ) const
        {
          if ( this->empty() ) return 0;
          if ( this->size() == 1 ) return this->seqlen;
          return this->get_left_len();
        }

        /**
         *  @brief  get the sequence length of the tail node (back node).
         */
          inline seqsize_type
        get_seqlen_tail( ) const
        {
          if ( this->empty() ) return 0;
          if ( this->size() == 1 ) return this->seqlen;
          return this->get_right_len();
        }

        /**
         *  @brief  getter function for seq.
         *
         *  @note The sequence is constructed on demand by calling this function.
         */
          inline const string_type&
        get_sequence( )
        {
          if ( this->seq.empty() ) this->seq = sequence( *this );
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
        set_left_by_len( offset_type value )
        {
          if ( value < 0 ) throw std::runtime_error( "Invalid left offset value" );
          if ( this->empty() ) {
            throw std::runtime_error( "Cannot set offset for an empty path" );
          }

          offset_type front_len = this->vargraph->node_length( this->front() );
          if ( value > front_len || value == 0 ) value = front_len;
          if ( this->size() == 1 && front_len - value >= this->get_right_len() ) {
            throw std::runtime_error( "left exceeds right on the one-node path" );
          }
          std::ptrdiff_t diff = value - this->get_left_len();
          if ( diff == 0 ) return;
          this->seqlen += diff;
          if ( !this->seq.empty() ) {
            if ( diff < 0 ) {
              this->seq.erase( 0, -diff );
            }
            else {
              auto nstr = this->vargraph->node_sequence( this->front() );
              this->seq.insert( 0, nstr.substr( front_len - value, diff ) );
            }
          }
          this->left = ( value == front_len ) ? 0 : value;
          this->initialized = false;
        }

          inline void
        set_right_by_len( offset_type value )
        {
          if ( value < 0 ) throw std::runtime_error( "Invalid right offset value" );
          if ( this->empty() ) {
            throw std::runtime_error( "Cannot set offset for an empty path" );
          }

          offset_type back_len = this->vargraph->node_length( this->back() );
          if ( value > back_len || value == 0 ) value = back_len;
          if ( this->size() == 1 && value <= this->get_head_offset() ) {
            throw std::runtime_error( "right exceeds left on the one-node path" );
          }
          std::ptrdiff_t diff = value - this->get_right_len();
          if ( diff == 0 ) return;
          this->seqlen += diff;
          if ( !this->seq.empty() ) {
            if ( diff < 0 ) {
              this->seq.resize( this->seqlen );
            }
            else {
              auto nstr = this->vargraph->node_sequence( this->back() );
              this->seq += nstr.substr( this->right, diff );  /**< @brief `right` is always non-zero here */
            }
          }
          this->right = ( value == back_len ) ? 0 : value;
          this->initialized = false;
        }

          inline void
        set_nodes( nodes_type value, offset_type l=0, offset_type r=0 );

        template< typename TIter >
            inline void
          set_nodes( TIter begin, TIter end, offset_type l=0, offset_type r=0 )
          {
            nodes_type nd( end - begin );
            std::copy( begin, end, nd.begin() );
            this->set_nodes( std::move( nd ), l, r );
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
          if ( this->right != 0 ) this->set_right_by_len( 0 );
          this->nodes.push_back( nid );
          this->seqlen += this->vargraph->node_length( nid );
          if ( !this->seq.empty() ) {
            this->seq += this->vargraph->node_sequence( nid );
          }
          this->initialized = false;
        }

        /**
         *  @brief  Append a new node to the back of the path.
         *
         *  @param  nid Node ID.
         *  @param  noff Node offset.
         *
         *  XXX NOTE: When appending the first node `noff` indicates the first locus
         *  should be included in the path. So, `noff == 0` will append all the node
         *  sequence to the path (i.e. `left == 0`). However, for further appending
         *  `noff` points to "one locus after" the locus need to be appended. In this
         *  case, `noff == 0` or `noff == len` will also append the whole node to the
         *  path (i.e. `right == 0`); where `len` is the sequence length of node `nid`.
         */
          inline void
        push_back( value_type const& nid, offset_type noff )
        {
          assert( this->vargraph != nullptr );
          bool first = this->empty();
          offset_type nlen = this->vargraph->node_length( nid );
          if ( noff < 0 ) noff = 0;
          this->initialized = false;
          if ( first ) {
            if ( noff >= nlen ) noff = nlen - 1;
            this->nodes.push_back( nid );
            this->seqlen += nlen - noff;
            this->left = noff ? seqlen : 0;
            assert( this->seq.empty() );
          }
          else {
            if ( this->right != 0 ) this->set_right_by_len( 0 );
            if ( noff > nlen || noff == 0 ) noff = nlen;
            this->nodes.push_back( nid );
            this->seqlen += noff;
            this->right = ( noff == nlen ) ? 0 : noff;
            if ( !this->seq.empty() ) {
              this->seq += this->vargraph->node_sequence( nid ).substr( 0, noff );
            }
          }
        }

          inline void
        pop_back( );

          inline void
        pop_front( );

          inline void
        clear( )
        {
          grem::clear( this->nodes );
          this->left = 0;
          this->right = 0;
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
          ASSERT( this->is_initialized() );  // Path must be initialized to be serialized
          grem::serialize( out, this->nodes );
          grem::serialize( out, static_cast< uint64_t >( this->left ) );
          grem::serialize( out, static_cast< uint64_t >( this->right ) );
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
          offset_type l, r;
          uint64_t offset;
          deserialize( in, tmp );
          deserialize( in, offset );
          l = offset;
          deserialize( in, offset );
          r = offset;
          this->set_nodes( std::move( tmp ), l, r );
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
          seqsize_type cursor = this->get_seqlen_head();
          this->bv_node_breaks[ cursor - 1 ] = 1;
          if ( this->size() > 1 ) {
            for ( auto it = this->begin()+1; it != this->end()-1; ++it ) {
                cursor += this->vargraph->node_length( *it );
                this->bv_node_breaks[ cursor - 1 ] = 1;
              }
            cursor += this->get_seqlen_tail();
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

        this->left = other.left;
        this->right = other.right;
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
        this->left = other.left;
        this->right = other.right;
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
    Path< TGraph, TSpec >::set_nodes( nodes_type value, offset_type l, offset_type r )
    {
      this->clear();
      if ( value.empty() ) return;

      assert( this->vargraph != nullptr );

      this->nodes = std::move( value );
      for ( auto n : this->nodes ) this->seqlen += this->vargraph->node_length( n );
      this->set_left_by_len( l );
      this->set_right_by_len( r );
    }  /* -----  end of method set_nodes  ----- */

  template< typename TGraph, typename TSpec >
      inline void
    Path< TGraph, TSpec >::pop_back( )
    {
      assert( this->vargraph != nullptr );

      if ( this->empty() ) return;

      this->seqlen -= this->get_seqlen_tail();
      this->nodes.pop_back();
      this->initialized = false;
      if ( !this->seq.empty() ) this->seq.resize( this->seqlen );
      this->right = 0;
      if ( this->empty() ) this->left = 0;
    }

  template< typename TGraph, typename TSpec >
      inline void
    Path< TGraph, TSpec >::pop_front( )
    {
      assert( this->vargraph != nullptr );

      if ( this->nodes.empty() ) return;

      offset_type diff = this->get_seqlen_head();
      this->seqlen -= diff;
      this->nodes.pop_front();
      this->initialized = false;
      if ( !this->seq.empty() ) this->seq = this->seq.substr( diff );
      this->left = 0;
      if ( this->empty() ) this->right = 0;
    }

  template< typename TPathSpec >
    class is_generic_path : public std::true_type {};

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

  template< >
    class is_generic_path< Micro > : public std::false_type {};

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
          this->last_node_rank = this->ss_nodes( this->size() ) + 1;
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
          return this->contains_by_rank( nrank );
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
              if ( curr <= prev || !this->contains_by_rank( curr ) ) return false;
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
              if ( curr >= prev || !this->contains_by_rank( curr ) ) return false;
              prev = curr;
            }

            return true;
          }
      private:
          inline bool
        contains_by_rank( typename graph_type::rank_type rank ) const
        {
          if ( rank == 0 ) return false;
          return this->nodes[ rank - 1 ] == 1;
        }
    };  /* -----  end of specialized template class Path  ----- */

  template< >
    class is_generic_path< Haplotype > : public std::false_type {};
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef PATH_BASE_H__  ----- */

/**
 *    @file  fmindex.hpp
 *   @brief  FM-Index wrapper on `sdsl::csa_wt`.
 *
 *  This is a wrapper module for `sdsl::csa_wt`.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Feb 13, 2018  10:10
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  PSI_FMINDEX_HPP__
#define  PSI_FMINDEX_HPP__

#include <fstream>
#include <string>

#include <seqan/index.h>
#include <sdsl/suffix_arrays.hpp>

#include "sequence.hpp"
#include "utils.hpp"


namespace grem {
  template< class TWT = sdsl::wt_huff<>, uint32_t TDens = 32, uint32_t TInvDens = 64 >
    struct FMIndex;

  template< typename TSpec >
    class is_fmindex : public std::false_type {
    };

  template< class TWT, uint32_t TDens, uint32_t TInvDens >
    class is_fmindex< FMIndex< TWT, TDens, TInvDens > > : public std::true_type {
    };
}  /* -----  end of namespace grem  ----- */

namespace seqan {
  /* Forwards */
  template< class TWT, uint32_t TDens, uint32_t TInvDens >
      void
    indexRequire(
        Index< grem::YaString< grem::DiskBased >, grem::FMIndex< TWT, TDens, TInvDens > >& index,
        FibreSALF );

  template< class TWT, uint32_t TDens, uint32_t TInvDens >
      void
    indexRequire(
        Index< grem::YaString< grem::InMemory >, grem::FMIndex< TWT, TDens, TInvDens > >& index,
        FibreSALF );

  template< class TWT, uint32_t TDens, uint32_t TInvDens >
      void
    indexRequire(
        Index< StringSet< grem::YaString< grem::DiskBased > >, grem::FMIndex< TWT, TDens, TInvDens > >& index,
        FibreSALF );

  template< class TWT, uint32_t TDens, uint32_t TInvDens >
      void
    indexRequire(
        Index< StringSet< grem::YaString< grem::InMemory > >, grem::FMIndex< TWT, TDens, TInvDens > >& index,
        FibreSALF );

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      typename Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >::text_type&
    getFibre( Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >& index, FibreText );

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      typename Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >::text_type const&
    getFibre( Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > const& index, FibreText );

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
    class Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TText text_type;                               /**< @brief input text type */
        typedef typename text_type::pos_type pos_type;
        typedef grem::FMIndex< TWT, TDens, TInvDens > spec_type;
        typedef sdsl::csa_wt< TWT, TDens, TInvDens > value_type;
        typedef typename value_type::string_type string_type;  /**< @brief output string type */
        typedef typename value_type::size_type savalue_type;
        typedef typename value_type::index_category index_category;
        typedef char char_type;
        typedef typename value_type::comp_char_type comp_char_type;
        /* ====================  LIFECYCLE     ======================================= */
        Index ( )
          : text_p( nullptr ), owner( true ) { }
        Index( text_type* t )
          : text_p( t ), owner( false ) { }
        Index( text_type& t )
          : Index( &t ) { }

        Index( Index const& ) = delete;
        Index( Index&& ) = default;
        Index& operator=( Index const& ) = delete;
        Index& operator=( Index&& ) = default;
        ~Index( ) noexcept
        {
          this->clear();
        }
        /* ====================  ACCESSORS     ======================================= */
          inline bool
        owns_text( ) const
        {
          return this->owner;
        }
        /* ====================  METHODS       ======================================= */
          inline savalue_type
        size( ) const
        {
          return this->fm.size();
        }

          inline bool
        empty( ) const
        {
          return this->size() == 0;
        }

          inline void
        clear_fibres( )
        {
          sdsl::util::clear( this->fm );
        }

          inline void
        clear( )
        {
          this->clear_fibres();
          if ( this->owner && this->text_p != nullptr ) {
            delete this->text_p;
            this->text_p = nullptr;
          }
        }

          inline bool
        constructible( ) const
        {
          return ( this->fm.size() == 0 ) && ( this->text_p != nullptr );
        }

          inline void
        serialize( std::ostream& out )
        {
          this->fm.serialize( out );
          this->text_p->serialize( out );
        }

          inline void
        load( std::istream& in )
        {
          this->clear();
          this->fm.load( in );
          assert( this->text_p == nullptr );
          this->text_p = new text_type();
          this->text_p->load( in );
          this->owner = true;
        }

        // :TODO:Wed Apr 04 13:17:\@cartoonist: FIXME: a Holder class should be
        //     responsible for the text fibre not client class.
        //     === PROBLEM ===
        //     The index has a pointer to the text fibre which it might not own it. It
        //     is valid until the text object lives. When a class aggregates an index
        //     (particularly, FMIndex) and its text fibre when an object of that class
        //     is moved to another one. The text fibre of new object's index is no
        //     longer valid since it points to the text fibre copy of the moved object
        //     not the new one.
        //     This is the problem for all classes that aggregates the index and its
        //     text fibre. The text fibre should be owned by a Holder class not the
        //     client class (e.g. PathSet, or PathIndex here).
        //     === SOLUTION ===
        //     This function is a work-around for this issue. The true solution would be
        //     wrapping the text fibre in a Holder class which is responsible for the
        //     memory management of the text fibre not the class aggregating the index
        //     and its text. In other words, the text fibre object will be handed to
        //     objects on the fly and finaly the last holder would free the memory.
        //     This Holder class can be developed from scratch or use smart pointers in
        //     C++ Standard Library.
          inline void
        set_text_fibre( text_type* ext_p, bool update=true )
        {
          if ( update ) this->clear_fibres();
          this->text_p = ext_p;
          this->owner = false;
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        value_type fm;
        text_type* text_p;
        bool owner;
        /* ====================  INTERFACE FUNCTIONS  ================================ */
          friend void
        indexRequire< TWT, TDens, TInvDens >( Index& index, FibreSALF );
          friend text_type&
        getFibre< TText, TWT, TDens, TInvDens >( Index&, FibreText );
          friend text_type const&
        getFibre< TText, TWT, TDens, TInvDens >( Index const&, FibreText );
        /* ====================  FRIENDSHIPS   ======================================= */
        friend class Finder< Index >;
        friend class Iter< Index, seqan::TopDown<> >;
        friend class Iter< Index, seqan::TopDown< seqan::ParentLinks<> > >;
        template< typename TPath, typename TSpec >
          friend class PathSet;
    };

  template< class TWT, uint32_t TDens, uint32_t TInvDens >
      inline void
    indexRequire(
        Index< grem::YaString< grem::DiskBased >, grem::FMIndex< TWT, TDens, TInvDens > >& index,
        FibreSALF )
    {
      if ( index.constructible() ) {
        std::string tmpdir = grem::get_tmpdir_env();
        sdsl::cache_config config;
        if ( tmpdir.size() != 0 ) {
          config.dir = std::move( tmpdir );
        }
        construct( index.fm, index.text_p->get_file_path(), config, 1 );
      }
    }

  template< class TWT, uint32_t TDens, uint32_t TInvDens >
      inline void
    indexRequire(
        Index< grem::YaString< grem::InMemory >, grem::FMIndex< TWT, TDens, TInvDens > >& index,
        FibreSALF )
    {
      if ( index.constructible() ) construct_im( index.fm, index.text_p->c_str(), 1 );
    }

  template< class TWT, uint32_t TDens, uint32_t TInvDens >
      inline void
    indexRequire(
        Index< StringSet< grem::YaString< grem::DiskBased > >, grem::FMIndex< TWT, TDens, TInvDens > >& index,
        FibreSALF )
    {
      if ( index.constructible() ) {
        std::string tmpdir = grem::get_tmpdir_env();
        sdsl::cache_config config;
        if ( tmpdir.size() != 0 ) {
          config.dir = std::move( tmpdir );
        }
        construct( index.fm, index.text_p->get_file_path(), config, 1 );
      }
    }

  template< class TWT, uint32_t TDens, uint32_t TInvDens >
      inline void
    indexRequire(
        Index< StringSet< grem::YaString< grem::InMemory > >, grem::FMIndex< TWT, TDens, TInvDens > >& index,
        FibreSALF )
    {
      if ( index.constructible() ) construct_im( index.fm, index.text_p->c_str(), 1 );
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline typename Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >::text_type&
    getFibre( Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >& index, FibreText )
    {
      return *index.text_p;
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline typename Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >::text_type const&
    getFibre( Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > const& index, FibreText )
    {
      return *index.text_p;
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline bool
    open( Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >& index,
        const std::string& file_name )
    {
      std::ifstream ifs( file_name, std::ifstream::in | std::ifstream::binary );
      if( !ifs ) return false;
      open( index, ifs );
      return true;
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline bool
    save( Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >& index,
        const std::string& file_name )
    {
      std::ofstream ofs( file_name, std::ofstream::out | std::ofstream::binary );
      if( !ofs ) return false;
      save( index, ofs );
      return true;
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline void
    clear( Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >& index )
    {
      index.clear();
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
    struct Value< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > > {
      typedef typename Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >::char_type Type;
    };

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
    struct Fibre< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, FibreText > {
      typedef typename Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >::text_type Type;
    };

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
    struct Fibre< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > const, FibreText > {
      typedef typename Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >::text_type const Type;
    };

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
    struct SAValue< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > > {
      typedef typename Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >::pos_type Type;
    };

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
    struct SAValue< const Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > > {
      typedef typename Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >::pos_type Type;
    };

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
    class Finder< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TText text_type;                               /**< @brief input text type */
        typedef typename text_type::pos_type pos_type;
        typedef Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > index_type;
        typedef typename index_type::string_type string_type;  /**< @brief output string type */
        typedef typename index_type::savalue_type savalue_type;
        typedef typename index_type::index_category index_category;
        typedef typename index_type::char_type char_type;
        typedef typename index_type::comp_char_type comp_char_type;
        /* ====================  ASSERTS       ======================================= */
        static_assert( std::is_same< sdsl::csa_tag, index_category >::value, "index category should be `csa`" );
        /* ====================  LIFECYCLE     ======================================= */
        Finder( index_type* i_p )
          : index_p( i_p ), occ_cur( 0 ), occ_end( 0 ), initiated( false ) {}
        Finder( index_type& i )
          : Finder( &i ) {}
        /* ====================  OPERATORS     ======================================= */
          inline Finder&
        operator++( )
        {
          if ( !this->at_end( ) ) ++this->occ_cur;
          return *this;
        }

          inline Finder&
        operator++( int )
        {
          auto tmp = *this;
          ++(*this);
          return tmp;
        }
        /* ====================  METHODS       ======================================= */
          inline bool
        empty( ) const
        {
          return !this->initiated;
        }

          inline bool
        at_end( ) const
        {
          return this->empty( ) || ( this->occ_cur > this->occ_end );
        }

          inline void
        clear( )
        {
          this->occ_cur = 0;
          this->occ_end = 0;
          this->initiated = false;
        }

          inline pos_type
        get_position( ) const
        {
          assert( !this->at_end() );
          return this->index_p->text_p->get_position( this->index_p->fm[ this->occ_cur ] );
        }

          inline savalue_type
        count( ) const
        {
          if ( this->at_end() ) return 0;
          return this->occ_end + 1 - this->occ_cur;
        }

        template< typename TIter >
            inline void
          backward_search( TIter pt_begin, TIter pt_end )
          {
            indexRequire( *(this->index_p), FibreSALF() );
            this->occ_cur = 0;
            this->occ_end = this->index_p->size()-1;
            while ( pt_begin < pt_end && this->occ_cur <= this->occ_end ) {
              --pt_end;
              sdsl::backward_search( this->index_p->fm, this->occ_cur, this->occ_end,
                  (char)*pt_end, this->occ_cur, this->occ_end);
            }
            this->initiated = true;
          }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        index_type* index_p;
        savalue_type occ_cur;
        savalue_type occ_end;
        bool initiated;
    };

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline bool
    empty( Finder< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > >& finder )
    {
      return finder.empty();
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline bool
    atEnd( Finder< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > >& finder )
    {
      return finder.at_end();
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline typename Finder< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > >::pos_type
    beginPosition( Finder< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > >& finder )
    {
      return finder.get_position();
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TPattern >
      inline bool
    find( Finder< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > >& finder,
        TPattern const& pattern )
    {
      using std::begin;
      using std::end;
      if ( empty( finder ) ) finder.backward_search( begin( pattern ), end( pattern ) );
      else ++finder;

      return !atEnd( finder );
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline bool
    find( Finder< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > >& finder,
        const char* pattern )
    {
      std::string p( pattern );
      return find( finder, p );
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline void
    clear( Finder< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > >& finder )
    {
      finder.clear();
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
    class Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TText text_type;                               /**< @brief input text type */
        typedef typename text_type::pos_type pos_type;
        typedef Index< TText, grem::FMIndex< TWT, TDens, TInvDens > > index_type;
        typedef typename index_type::string_type string_type;  /**< @brief output string type */
        typedef typename index_type::savalue_type savalue_type;
        typedef typename index_type::index_category index_category;
        typedef typename index_type::char_type char_type;
        typedef typename index_type::comp_char_type comp_char_type;
        typedef typename std::vector< pos_type > occs_type;
        typedef typename std::pair< savalue_type, savalue_type > range_type;
        /* ====================  ASSERTS       ======================================= */
        static_assert( std::is_same< sdsl::csa_tag, index_category >::value, "index category should be `csa`" );
        /* ====================  LIFECYCLE     ======================================= */
        Iter( index_type* i_p )
          : index_p( i_p ), occ_cur( 0 ), occ_end( 0 ), initialized( false ) { }
        Iter( index_type& i )
          : Iter( &i ) { }
        /* ====================  METHODS       ======================================= */
          inline bool
        is_initialized( ) const
        {
          return this->initialized;
        }

          inline bool
        at_end( ) const
        {
          return ( !this->is_initialized() ) || ( this->occ_cur > this->occ_end );
        }

          inline bool
        is_root( ) const
        {
          return ( !this->is_initialized() ) || ( this->count() == this->index_size() );
        }

          inline void
        clear( )
        {
          this->occ_cur = 0;
          this->occ_end = 0;
          this->initialized = false;
          this->history.clear();
        }

          inline savalue_type
        get_raw_position( savalue_type i ) const
        {
          ASSERT( this->is_initialized() );
          assert( this->occ_cur + i <= this->occ_end );
          return this->index_p->fm[ this->occ_cur + i ];
        }

          inline savalue_type
        get_raw_position( savalue_type i )
        {
          if ( !this->is_initialized() ) this->init();
          const Iter* _this = this;
          return _this->get_raw_position( i );
        }

          inline pos_type
        get_position( savalue_type i ) const
        {
          ASSERT( this->is_initialized() );
          return this->index_p->text_p->get_position( this->get_raw_position( i ) );
        }

          inline pos_type
        get_position( savalue_type i )
        {
          if ( !this->is_initialized() ) this->init();
          const Iter* _this = this;
          return _this->get_position( i );
        }

          inline occs_type
        get_occurrences( ) const
        {
          savalue_type occs_no;
          if ( this->is_root() ) occs_no = 0;
          else occs_no = this->count();
          occs_type occs( occs_no );

          for ( savalue_type i = 0; i < occs_no; ++i ) {
            occs[i] = this->get_position( i );
          }
          return occs;
        }

          inline void
        reserve_history( savalue_type size )
        {
          this->history.reserve( size );
        }

          inline savalue_type
        history_size( ) const
        {
          return this->history.size();
        }

          inline void
        history_push( )
        {
          this->history.push_back( std::make_pair( this->occ_cur, this->occ_end ) );
        }

          inline void
        history_pop( )
        {
          this->occ_cur = this->history.back().first;
          this->occ_end = this->history.back().second;
          this->history.pop_back();
        }

          inline savalue_type
        index_size( ) const
        {
          return this->index_p->size();
        }

          inline savalue_type
        count( ) const
        {
          ASSERT( this->is_initialized() );
          if ( this->at_end() ) return 0;
          return this->occ_end + 1 - this->occ_cur;
        }

          inline savalue_type
        count( )
        {
          if ( !this->is_initialized() ) this->init();
          const Iter* _this = this;
          return _this->count();
        }

          inline void
        init( )
        {
          indexRequire( *(this->index_p), FibreSALF() );
          this->occ_cur = 0;
          this->occ_end = this->index_size() - 1;
          this->initialized = true;
        }

          inline char_type
        last_char( )
        {
          assert( !this->at_end() );
          return first_row_symbol( this->occ_cur, this->index_p->fm );
        }

          inline savalue_type
        go_down( char_type c )
        {
          if ( !this->is_initialized() ) this->init();
          this->history_push();
          savalue_type no = sdsl::backward_search( this->index_p->fm,
              this->occ_cur, this->occ_end, (char)c,
              this->occ_cur, this->occ_end );

          if ( no == 0 ) this->history_pop();
          return no;
        }

          inline savalue_type
        go_down_gt( char_type c )  /**< @brief go down with any character larger than c */
        {
          comp_char_type cc = this->index_p->fm.char2comp[c];
          for ( ++cc; cc < this->index_p->fm.sigma; ++cc ) {
            auto no = this->go_down( this->index_p->fm.comp2char[cc] );
            if ( no != 0 ) return no;
          }
          return 0;
        }

          inline savalue_type
        rep_length( ) const
        {
          return this->history_size();
        }

          inline savalue_type
        parent_edge_length( ) const
        {
          if ( this->is_root() ) return 0;
          return 1;
        }

          inline string_type
        parent_edge_label( ) const
        {
          if ( this->is_root() ) return string_type("");
          auto a_loc = this->get_raw_position( 0 );
          return extract( this->index_p->fm, a_loc, a_loc );
        }

          inline string_type
        representative( ) const
        {
          auto a_loc = this->get_raw_position( 0 );
          return extract( this->index_p->fm, a_loc, a_loc + this->rep_length() - 1 );
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        index_type* index_p;
        savalue_type occ_cur;
        savalue_type occ_end;
        bool initialized;
        std::vector< range_type > history;
    };

  template< typename TIterator, typename TSAValue >
      inline void
    reserveHistory( TIterator& iter, TSAValue size )
    { /* NO-OP */ }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline void
    reserveHistory( Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >& iter,
        typename Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >::savalue_type size )
    {
      iter.reserve_history( size );
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline bool
    isRoot( Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >& iter )
    {
      return iter.is_root();
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline void
    goRoot( Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >& iter )
    {
      if ( isRoot( iter ) ) return;
      iter.init();
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline bool
    goDown( Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >& iter,
        typename Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >::char_type c )
    {
      return iter.go_down( c ) != 0;
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline bool
    goUp( Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >& iter )
    {
      if ( isRoot( iter ) ) return false;
      iter.history_pop();
      return true;
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline bool
    goRight( Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >& iter )
    {
      if ( isRoot( iter ) ) return false;
      auto c = iter.last_char();
      goUp( iter );
      if ( iter.go_down_gt( c ) != 0 ) return true;
      goDown( iter, c );  /**< @brief cannot go right, undo goUp. */
      return false;
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline typename Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >::savalue_type
    parentEdgeLength( Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > > const& iter )
    {
      return iter.parent_edge_length();
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline typename Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >::string_type
    parentEdgeLabel( Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > > const& iter )
    {
      return iter.parent_edge_label();
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline typename Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >::savalue_type
    repLength( Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > > const& iter )
    {
      return iter.rep_length();
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline typename Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >::string_type
    representative( Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > > const& iter )
    {
      return iter.representative();
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline typename Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > >::occs_type
    getOccurrences( Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > > const& iter )
    {
      return iter.get_occurrences();
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
    struct Iterator< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > > {
      typedef Iter< Index< TText, grem::FMIndex< TWT, TDens, TInvDens > >, TopDown< TSpec > > Type;
    };
}  /* -----  end of namespace seqan  ----- */

#endif  /* --- #ifndef PSI_FMINDEX_HPP__ --- */

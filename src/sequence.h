/**
 *    @file  sequence.h
 *   @brief  Sequence abstract data types.
 *
 *  Sequence header file contains type definitions, abstract data types, and helper
 *  interface functions (mostly implemented at the top of SeqAn sequence library; i.e.
 *  `seqan/sequence.h`) to work with sequences and sequence sets.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Mar 07, 2017  12:50
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  SEQUENCE_H__
#define  SEQUENCE_H__

#include <fstream>
#include <stdexcept>
#include <memory>

#include <seqan/seq_io.h>
#include <kseq++/seqio.h>
#include <seqan/sequence.h>
#include <sdsl/bit_vectors.hpp>

#include "utils.h"
#include "logger.h"


#define SEQUENCE_DEFAULT_SENTINEL_CHAR '$'


namespace grem {
  // :TODO:Tue Sep 05 09:36:\@cartoonist: clear the code up from direct usage of seqan::Dna5QString.
  /* Typedefs  ------------------------------------------------------------------- */
  template< typename TSpec = seqan::Owner<> >
    using CharStringSet = seqan::StringSet< seqan::CharString, TSpec >;
  template< typename TSpec = seqan::Owner<> >
    using Dna5QStringSet = seqan::StringSet< seqan::Dna5QString, TSpec >;
  typedef seqan::Dependent< seqan::Generous > Dependent;
  /* END OF Typedefs  ------------------------------------------------------------ */

  /* Meta-functions  ------------------------------------------------------------- */

  template< typename TStringSet >
    class MakeOwner;

  template< typename TStringSet >
    class MakeDependent;

  template< typename TText >
    class MakeOwner< seqan::StringSet< TText, grem::Dependent > > {
      public:
        typedef seqan::StringSet< TText, seqan::Owner<> > Type;
    };

  template< typename TText >
    class MakeOwner< seqan::StringSet< TText, seqan::Owner<> > > {
      public:
        typedef seqan::StringSet< TText, seqan::Owner<> > Type;
    };

  template< typename TText >
    class MakeDependent< seqan::StringSet< TText, seqan::Owner<> > > {
      public:
        typedef seqan::StringSet< TText, grem::Dependent > Type;
    };

  template< typename TText >
    class MakeDependent< seqan::StringSet< TText, grem::Dependent > > {
      public:
        typedef seqan::StringSet< TText, grem::Dependent > Type;
    };

  template< typename TContainer >
    class Ownership;

  template< typename TText, typename TSpec >
    class Ownership< seqan::StringSet< TText, TSpec > > {
      public:
        typedef grem::Dependent Type;
    };

  template< typename TText >
    class Ownership< seqan::StringSet< TText, seqan::Owner<> > > {
      public:
        typedef seqan::Owner<> Type;
    };

  /* END OF Meta-functions  ------------------------------------------------------ */

  /* Data structures  ------------------------------------------------------------ */

  /* Forwards */
  template< typename TSpec >
    class YaString;

  struct DiskBasedStrategy;
  struct InMemoryStrategy;
  typedef seqan::Tag< DiskBasedStrategy > DiskBased;
  typedef seqan::Tag< InMemoryStrategy > InMemory;

  typedef YaString< DiskBased > DiskString;
  typedef YaString< InMemory > MemString;

  template< >
    class YaString< InMemory > : public std::string {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef std::string string_type;
        typedef InMemory spec_type;
        typedef spec_type device_type;
        typedef string_type::size_type size_type;
        typedef size_type pos_type;
        /* ====================  LIFECYCLE     ======================================= */
        using std::string::string;
        /* ====================  METHODS       ======================================= */
          inline void
        serialize( std::ostream& out )
        {
          grem::serialize( out, this->begin(), this->end() );
        }

          inline void
        load( std::istream& in )
        {
          this->clear();
          grem::deserialize( in, *this, std::back_inserter( *this ) );
        }

          inline pos_type
        get_position( pos_type p )
        {
          return p;
        }

          inline size_type
        raw_length( ) const
        {
          return std::string::size();  /* call base class size function */
        }
    };

  template< >
    class YaString< DiskBased > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef DiskBased spec_type;
        typedef spec_type device_type;
        typedef std::string string_type;
        typedef string_type::size_type size_type;
        typedef uint64_t sint_type;       /**< @brief Type for serializing the length. */
        typedef size_type pos_type;
        /* ====================  LIFECYCLE     ======================================= */
        YaString( const string_type& data, std::string _fpath )
          : fpath( std::move( _fpath ) ), out( this->fpath ), len( 0 )
        {
          this->append( data );
        }

        YaString( string_type&& data, std::string _fpath )
          : YaString( data, std::move( _fpath ) ) { }

        YaString( const string_type& data )
          : YaString( data, get_tmpfile() ) { }

        YaString( string_type&& data )
          : YaString( data ) { }

        YaString( )
          : YaString( string_type() ) { }

        YaString( const YaString& ) = delete;
        YaString( YaString&& ) = default;
        YaString& operator=( const YaString& ) = delete;
        YaString& operator=( YaString&& ) = default;
        ~YaString( ) = default;
        /* ====================  ACCESSORS     ======================================= */
          inline std::string
        get_file_path()
        {
          this->close();
          return this->fpath;
        }
        /* ====================  OPERATORS     ======================================= */
        YaString& operator=( const string_type& data );

          inline YaString&
        operator+=( const string_type& str )
        {
          this->append( str );
          return *this;
        }

          inline YaString&
        operator+=( string_type&& str )
        {
          *this += str;  // re-use above overload of `operator+=`.
          return *this;
        }
        /* ====================  METHODS       ======================================= */
          inline bool
        is_open( )
        {
          return this->out.is_open();
        }

          inline size_type
        length( ) const
        {
          return this->len;
        }

          inline size_type
        raw_length( ) const
        {
          /* to avoid ambiguity of calling `length()` which will be overriden by StringSet class. */
          return YaString::length();
        }

          inline pos_type
        get_position( pos_type p )
        {
          return p;
        }

          inline void
        clear( )
        {
          this->close();
          this->fpath = get_tmpfile();
          this->out.open( this->fpath );
          this->len = 0;
        }

        template< typename TSize, typename TTag = seqan::Exact >
          inline void reserve( TSize size, TTag = TTag() ) { /* NO-OP */ }

          inline void
        serialize( std::ostream& out )
        {
          out.flush();
          grem::serialize( out, this->fpath.begin(), this->fpath.end() );
          grem::serialize( out, static_cast< sint_type >( len ) );
        }

          inline void
        load( std::istream& in )
        {
          this->close();
          this->fpath.clear();
          grem::deserialize( in, this->fpath, std::back_inserter( this->fpath ) );
          if ( !readable( this->fpath ) ) {
            get_logger( "main" )->warn(
                "File '{}' does not exist: disk-based string content cannot be read.",
                this->fpath );
          }

          if ( appendable( this->fpath ) ) {
            this->out.open( this->fpath, std::ofstream::app );
          }

          sint_type l;
          grem::deserialize( in, l );
          this->len = l;
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        std::string fpath;
        std::ofstream out;
        size_type len;
        /* ====================  METHODS       ======================================= */
          inline void
        close( )
        {
          if ( this->is_open() ) this->out.close();
        }

          inline void
        append( const string_type& data )
        {
          if ( data.length() == 0 ) return;
          if ( !this->is_open() )
            throw std::runtime_error( "attempting to write to a closed disk-based string." );
          this->len += data.length();
          this->out << data;
        }
    };

    inline YaString< DiskBased >::size_type
  length( YaString< DiskBased >& dstr )
  {
    return dstr.length();
  }

    inline void
  clear( YaString< DiskBased >& dstr )
  {
    dstr.clear();
  }

    inline YaString< DiskBased >&
  YaString< DiskBased >::operator=( const YaString< DiskBased >::string_type& data )
  {
    this->clear();
    this->append( data );
    return *this;
  }

  template< typename TSize >
      inline void
    reserve( YaString< DiskBased >& dstr, TSize const size )
    {
      dstr.reserve( size );
    }

  template< typename TString >
    class YaInfix
    : public std::pair< typename TString::size_type, typename TString::size_type > { };

  template< typename TString >
      inline typename TString::size_type
    length( YaInfix< TString > const& inf )
    {
      assert( inf.second >= inf.first );
      return inf.second - inf.first;
    }

  template< typename T1, typename T2 >
    class YaPair : public std::pair< T1, T2 > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef std::pair< T1, T2 > base_type;
        /* ====================  LIFECYCLE     ======================================= */
        YaPair( T1 t1, T2 t2 )
          : base_type( t1, t2 ), i1( this->first ), i2( this->second )
        { }

        YaPair( )
          : base_type( ), i1( this->first ), i2( this->second )
        { }

        YaPair( YaPair const& other )
          : base_type( other ), i1( this->first ), i2( this->second )
        { }

        YaPair( YaPair&& other )
          : base_type( std::move( other ) ), i1( this->first ), i2( this->second )
        { }

        template< typename TSpec >
          YaPair( seqan::Pair< T1, T2, TSpec > const& other )
            : base_type( other.i1, other.i2 ), i1( this->first ), i2( this->second )
          { }

        template< typename TSpec >
          YaPair( seqan::Pair< T1, T2, TSpec >&& other )
            : base_type( other.i1, other.i2 ), i1( this->first ), i2( this->second )
          { }

        ~YaPair( ) = default;

          inline YaPair&
        operator=( YaPair const& other )
        {
          this->first = other.first;
          this->second = other.second;
          return *this;
        }

          inline YaPair&
        operator=( YaPair&& other )
        {
          this->first = other.first;
          this->second = other.second;
          return *this;
        }

        template< typename TSpec >
            inline YaPair&
          operator=( seqan::Pair< T1, T2, TSpec > const& other )
          {
            this->first = other.i1;
            this->second = other.i2;
            return *this;
          }

        template< typename TSpec >
            inline YaPair&
          operator=( seqan::Pair< T1, T2, TSpec >&& other )
          {
            this->first = other.i1;
            this->second = other.i2;
            return *this;
          }
        /* ====================  DATA MEMBERS  ======================================= */
        T1& i1;
        T2& i2;
    };
}  /* -----  end of namespace grem  ----- */

namespace seqan {
  template< >
    class StringSet< grem::DiskString, Owner<> > : public grem::DiskString {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef grem::DiskString value_type;
        typedef grem::DiskString::string_type string_type;
        typedef string_type::size_type stringsize_type;
        typedef uint64_t size_type;       /**< @brief Type for serializing the length. */
        typedef grem::DiskBased device_type;
        typedef grem::YaPair< size_type, stringsize_type > pos_type;
        /* ====================  LIFECYCLE     ======================================= */
        StringSet( )
          : grem::DiskString( ), count( 0 ), initialized( false ), bv_str_breaks( )
        { }

        StringSet( std::string _fpath )
          : grem::DiskString( string_type(), std::move( _fpath ) ),
          count( 0 ), initialized( false ), bv_str_breaks( )
        { }

        StringSet( const StringSet& ) = delete;
        StringSet& operator=( const StringSet& ) = delete;

        StringSet( StringSet&& other )
          : grem::DiskString( std::move( other ) )
        {
          this->count = other.count;
          this->bv_str_breaks.swap( other.bv_str_breaks );
          this->initialize();
          sdsl::util::clear( other.rs_str_breaks );
          sdsl::util::clear( other.ss_str_breaks );
        }

        StringSet& operator=( StringSet&& other )
        {
          grem::DiskString::operator=( std::move( other ) );
          this->count = other.count;
          sdsl::util::clear( this->bv_str_breaks );
          this->bv_str_breaks.swap( other.bv_str_breaks );
          sdsl::util::clear( this->rs_str_breaks );
          sdsl::util::clear( this->ss_str_breaks );
          this->initialize();
          sdsl::util::clear( other.rs_str_breaks );
          sdsl::util::clear( other.ss_str_breaks );
          return *this;
        }

        ~StringSet( ) = default;
        /* ====================  CONST MEMBERS ======================================= */
        const char SENTINEL = SEQUENCE_DEFAULT_SENTINEL_CHAR;
        /* ====================  OPERATORS     ======================================= */
          inline grem::YaInfix< StringSet >
        operator[]( size_type idx ) const
        {
          ASSERT( this->is_initialized() );
          grem::YaInfix< StringSet > retval;
          retval.first = this->select( idx );
          retval.second = this->select( idx + 1 ) - 1;
          return retval;
        }

          inline grem::YaInfix< StringSet >
        operator[]( size_type idx )
        {
          if ( !this->is_initialized() ) this->initialize();
          const StringSet* _this = this;
          return _this->operator[]( idx );
        }
        /* ====================  METHODS       ======================================= */
        void push_back( const string_type& str );

          inline void
        push_back( string_type&& str )
        {
          this->push_back( str );
        }

          inline size_type
        get_id( stringsize_type strpos )
        {
          if ( !this->is_initialized() ) this->initialize();
          return this->rank( strpos );
        }

          inline stringsize_type
        get_offset( stringsize_type strpos )
        {
          if ( !this->is_initialized() ) this->initialize();
          return strpos - this->select( this->rs_str_breaks( strpos ) );
        }

          inline pos_type
        get_position( stringsize_type strpos )
        {
          pos_type pos( 0, 0 );
          pos.first = this->get_id( strpos );
          pos.second = this->get_offset( strpos );
          return pos;
        }

          inline void
        initialize( )
        {
          this->shrink_bv_str_breaks();
          sdsl::util::init_support( this->rs_str_breaks, &this->bv_str_breaks );
          sdsl::util::init_support( this->ss_str_breaks, &this->bv_str_breaks );
          this->initialized = true;
        }

          inline bool
        is_initialized( ) const
        {
          return this->initialized;
        }

          inline size_type
        length( ) const
        {
          return this->count;
        }

          inline void
        clear( )
        {
          grem::DiskString::clear();
          this->count = 0;
          sdsl::util::clear( this->bv_str_breaks );
          sdsl::util::clear( this->rs_str_breaks );
          sdsl::util::clear( this->ss_str_breaks );
          initialized = false;
        }

        template< typename TSize, typename TTag = seqan::Exact >
          inline void reserve( TSize, TTag = TTag() ) { /* NO-OP */ }

          inline void
        serialize( std::ostream& out )
        {
          this->shrink_bv_str_breaks();
          grem::DiskString::serialize( out );
          grem::serialize( out, this->count );
          this->bv_str_breaks.serialize( out );
        }

          inline void
        load( std::istream& in )
        {
          this->clear();
          grem::DiskString::load( in );
          grem::deserialize( in, this->count );
          this->bv_str_breaks.load( in );
          this->initialize();
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        size_type count;
        bool initialized;
        sdsl::bit_vector bv_str_breaks;
        sdsl::bit_vector::rank_1_type rs_str_breaks;
        sdsl::bit_vector::select_1_type ss_str_breaks;
        /* ====================  METHODS       ======================================= */
          inline size_type
        rank( stringsize_type strpos ) const
        {
          assert( this->is_initialized() );
          assert( this->bv_str_breaks[ strpos ] != 1 );
          return this->rs_str_breaks( strpos );
        }

          inline stringsize_type
        select( size_type r ) const
        {
          assert( this->is_initialized() );
          if ( r == 0 ) return 0;
          return this->ss_str_breaks( r ) + 1;
        }

          inline void
        shrink_bv_str_breaks( )
        {
          this->bv_str_breaks.resize( this->raw_length() + 1 );  // shrink
        }
    };

    inline void
  StringSet< grem::DiskString, Owner<> >::push_back(
      typename StringSet< grem::DiskString, Owner<> >::string_type const& str )
  {
    if ( this->length() != 0 ) *this += std::string( 1, SENTINEL );
    ++this->count;
    *this += str;
    stringsize_type breakpoint = this->raw_length();
    if ( breakpoint >= this->bv_str_breaks.size() ) {
      sdsl::bit_vector new_bv( grem::roundup64( breakpoint + 1 ), 0 );
      grem::bv_icopy( this->bv_str_breaks, new_bv, 0, breakpoint - str.size() );
      sdsl::util::assign( this->bv_str_breaks, std::move( new_bv ) );
    }
    this->bv_str_breaks[ breakpoint ] = 1;
    this->initialized = false;
  }

    inline void
  push_back( StringSet< grem::DiskString, Owner<> >& dstr,
      typename StringSet< grem::DiskString, Owner<> >::string_type const& str )
  {
    dstr.push_back( str );
  }

    inline void
  push_back( StringSet< grem::DiskString, Owner<> >& dstr,
      typename StringSet< grem::DiskString, Owner<> >::string_type&& str )
  {
    dstr.push_back( str );
  }

    inline void
  appendValue( StringSet< grem::DiskString, Owner<> >& dstr,
      StringSet< grem::DiskString, Owner<> >::string_type& str )
  {
    push_back( dstr, str );
  }

    inline void
  appendValue( StringSet< grem::DiskString, Owner<> >& dstr,
      StringSet< grem::DiskString, Owner<> >::string_type&& str )
  {
    push_back( dstr, str );
  }

    inline StringSet< grem::DiskString, Owner<> >::size_type
  length( const StringSet< grem::DiskString, Owner<> >& dstr )
  {
    return dstr.length();
  }

    inline void
  clear( StringSet< grem::DiskString, Owner<> >& dstr )
  {
    dstr.clear();
  }

  template< typename TSize >
      inline void
    reserve( StringSet< grem::DiskString, Owner<> >& dstr,
        TSize const size )
    {
      dstr.reserve( size );
    }

  template< >
    class StringSet< grem::MemString, Owner<> > : public grem::MemString {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef grem::MemString value_type;
        typedef grem::MemString::string_type string_type;
        typedef string_type::size_type stringsize_type;
        typedef uint64_t size_type;
        typedef grem::InMemory device_type;
        typedef grem::YaPair< size_type, stringsize_type > pos_type;
        /* ====================  LIFECYCLE     ======================================= */
        StringSet( )
          : grem::MemString( ), count( 0 ), initialized( false ), bv_str_breaks( )
        { }

        StringSet( const StringSet& other )
          : grem::MemString( other )
        {
          this->count = other.count;
          this->bv_str_breaks = other.bv_str_breaks;
          this->initialize();
        }

        StringSet( StringSet&& other )
          : grem::MemString( std::move( other ) )
        {
          this->count = other.count;
          this->bv_str_breaks.swap( other.bv_str_breaks );
          this->initialize();
          sdsl::util::clear( other.rs_str_breaks );
          sdsl::util::clear( other.ss_str_breaks );
        }

        StringSet& operator=( const StringSet& other )
        {
          grem::MemString::operator=( other );
          this->count = other.count;
          this->bv_str_breaks = other.bv_str_breaks;
          this->initialize();
          return *this;
        }

        StringSet& operator=( StringSet&& other )
        {
          grem::MemString::operator=( std::move( other ) );
          this->count = other.count;
          sdsl::util::clear( this->bv_str_breaks );
          this->bv_str_breaks.swap( other.bv_str_breaks );
          sdsl::util::clear( this->rs_str_breaks );
          sdsl::util::clear( this->ss_str_breaks );
          this->initialize();
          sdsl::util::clear( other.rs_str_breaks );
          sdsl::util::clear( other.ss_str_breaks );
          return *this;
        }

        ~StringSet( ) = default;
        /* ====================  CONST MEMBERS ======================================= */
        const char SENTINEL = SEQUENCE_DEFAULT_SENTINEL_CHAR;
        /* ====================  OPERATORS     ======================================= */
          inline grem::YaInfix< StringSet >
        operator[]( size_type idx ) const
        {
          ASSERT( this->is_initialized() );
          grem::YaInfix< StringSet > retval;
          retval.first = this->select( idx );
          retval.second = this->select( idx + 1 ) - 1;
          return retval;
        }

          inline grem::YaInfix< StringSet >
        operator[]( size_type idx )
        {
          if ( !this->is_initialized() ) this->initialize();
          const StringSet* _this = this;
          return _this->operator[]( idx );
        }
        /* ====================  METHODS       ======================================= */
        void push_back( const string_type& str );

          inline void
        push_back( string_type&& str )
        {
          this->push_back( str );
        }

          inline size_type
        get_id( stringsize_type strpos )
        {
          if ( !this->is_initialized() ) this->initialize();
          return this->rank( strpos );
        }

          inline stringsize_type
        get_offset( stringsize_type strpos )
        {
          if ( !this->is_initialized() ) this->initialize();
          return strpos - this->select( this->rs_str_breaks( strpos ) );
        }

          inline pos_type
        get_position( stringsize_type strpos )
        {
          pos_type pos( 0, 0 );
          pos.first = this->get_id( strpos );
          pos.second = this->get_offset( strpos );
          return pos;
        }

          inline void
        initialize( )
        {
          this->shrink_bv_str_breaks( );
          sdsl::util::init_support( this->rs_str_breaks, &this->bv_str_breaks );
          sdsl::util::init_support( this->ss_str_breaks, &this->bv_str_breaks );
          this->initialized = true;
        }

          inline bool
        is_initialized( ) const
        {
          return this->initialized;
        }

          inline size_type
        length( ) const
        {
          return this->count;
        }

          inline void
        clear( )
        {
          grem::MemString::clear();
          this->count = 0;
          sdsl::util::clear( this->bv_str_breaks );
          sdsl::util::clear( this->rs_str_breaks );
          sdsl::util::clear( this->ss_str_breaks );
          initialized = false;
        }

        template< typename TSize, typename TTag = seqan::Exact >
          inline void reserve( TSize, TTag = TTag() ) { /* NO-OP */ }

          inline void
        serialize( std::ostream& out )
        {
          this->shrink_bv_str_breaks( );
          grem::MemString::serialize( out );
          grem::serialize( out, this->count );
          this->bv_str_breaks.serialize( out );
        }

          inline void
        load( std::istream& in )
        {
          this->clear();
          grem::MemString::load( in );
          grem::deserialize( in, this->count );
          this->bv_str_breaks.load( in );
          this->initialize();
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        size_type count;
        bool initialized;
        sdsl::bit_vector bv_str_breaks;
        sdsl::bit_vector::rank_1_type rs_str_breaks;
        sdsl::bit_vector::select_1_type ss_str_breaks;
        /* ====================  METHODS       ======================================= */
          inline size_type
        rank( stringsize_type strpos ) const
        {
          assert( this->is_initialized() );
          assert( this->bv_str_breaks[ strpos ] != 1 );
          return this->rs_str_breaks( strpos );
        }

          inline stringsize_type
        select( size_type r ) const
        {
          assert( this->is_initialized() );
          if ( r == 0 ) return 0;
          return this->ss_str_breaks( r ) + 1;
        }

          inline void
        shrink_bv_str_breaks( )
        {
          this->bv_str_breaks.resize( this->raw_length() + 1 );  // shrink
        }
    };

    inline void
  StringSet< grem::MemString, Owner<> >::push_back(
      typename StringSet< grem::MemString, Owner<> >::string_type const& str )
  {
    if ( this->length() != 0 ) *this += std::string( 1, SENTINEL );
    ++this->count;
    *this += str;
    stringsize_type breakpoint = this->raw_length();
    if ( breakpoint >= this->bv_str_breaks.size() ) {
      sdsl::bit_vector new_bv( grem::roundup64( breakpoint + 1 ), 0 );
      grem::bv_icopy( this->bv_str_breaks, new_bv, 0, breakpoint - str.size() );
      sdsl::util::assign( this->bv_str_breaks, std::move( new_bv ) );
    }
    this->bv_str_breaks[ breakpoint ] = 1;
    this->initialized = false;
  }

    inline void
  push_back( StringSet< grem::MemString, Owner<> >& dstr,
      typename StringSet< grem::MemString, Owner<> >::string_type const& str )
  {
    dstr.push_back( str );
  }

    inline void
  push_back( StringSet< grem::MemString, Owner<> >& dstr,
      typename StringSet< grem::MemString, Owner<> >::string_type&& str )
  {
    dstr.push_back( str );
  }

    inline void
  appendValue( StringSet< grem::MemString, Owner<> >& dstr,
      StringSet< grem::MemString, Owner<> >::string_type& str )
  {
    push_back( dstr, str );
  }

    inline void
  appendValue( StringSet< grem::MemString, Owner<> >& dstr,
      StringSet< grem::MemString, Owner<> >::string_type&& str )
  {
    push_back( dstr, str );
  }

    inline StringSet< grem::MemString, Owner<> >::size_type
  length( const StringSet< grem::MemString, Owner<> >& dstr )
  {
    return dstr.length();
  }

    inline void
  clear( StringSet< grem::MemString, Owner<> >& dstr )
  {
    dstr.clear();
  }

  template< typename TSize >
      inline void
    reserve( StringSet< grem::MemString, Owner<> >& dstr,
        TSize const size )
    {
      dstr.reserve( size );
    }
}  /* -----  end of namespace seqan  ----- */

namespace grem {
  /* StringSets interface functions */
  template< typename TText >
      inline typename seqan::Id< seqan::StringSet< TText, seqan::Owner<> > >::Type
    position_to_id( const seqan::StringSet< TText, seqan::Owner<> >& strset,
        typename seqan::Id< seqan::StringSet< TText, seqan::Owner<> > >::Type rel_id )
    {
      if ( rel_id >= length( strset ) || rel_id < 0 ) {
        throw std::runtime_error( "position out of range" );
      }
      return rel_id;  /* Relative id is the absolute one in Owner string set */
    }

  template< typename TText >
      inline typename seqan::Id< seqan::StringSet< TText, seqan::Owner<> > >::Type
    position_to_id( const seqan::StringSet< TText, seqan::Owner<> >& strset,
        typename seqan::StringSetPosition< seqan::StringSet< TText, seqan::Owner<> > >::Type const& pos )
    {
      return position_to_id( strset, pos.i1 );
    }

  template< typename TText >
      inline typename seqan::Position< seqan::StringSet< TText, seqan::Owner<> > >::Type
    position_to_offset( const seqan::StringSet< TText, seqan::Owner<> >& strset,
        typename seqan::StringSetPosition< seqan::StringSet< TText, seqan::Owner<> > >::Type const& pos )
    {
      using seqan::length;
      if ( pos.i2 >= length( strset[pos.i1] ) || pos.i2 < 0 ) {
        throw std::runtime_error( "position out of range" );
      }
      return pos.i2;
    }

  /* Forwards */
  template< typename TStringSet >
    class Records;

  template< typename TText, typename TSpec >
    class Ownership< Records< seqan::StringSet< TText, TSpec > > > {
      public:
        typedef grem::Dependent Type;
    };

  template< typename TText >
    class Ownership< Records< seqan::StringSet< TText, seqan::Owner<> > > > {
      public:
        typedef seqan::Owner<> Type;
    };

  /* Records interface functions */
  template< typename TText >
      inline typename Records< seqan::StringSet< TText, seqan::Owner<> > >::TId
    position_to_id( const Records< seqan::StringSet< TText, seqan::Owner<> > >& records,
        typename Records< seqan::StringSet< TText, seqan::Owner<> > >::TId rec_id )
    {
      if ( rec_id >= length( records.str ) || rec_id < 0 ) {
        throw std::runtime_error( "position out of range" );
      }
      return records.position_to_id( rec_id );
    }

  template< typename TText >
      inline typename Records< seqan::StringSet< TText, grem::Dependent > >::TId
    position_to_id( const Records< seqan::StringSet< TText, grem::Dependent > >& records,
        typename Records< seqan::StringSet< TText, grem::Dependent > >::TId rec_id )
    {
      if ( rec_id >= length( records.str ) || rec_id < 0 ) {
        throw std::runtime_error( "position out of range" );
      }
      assert( records.o_str != nullptr );
      return records.rec_offset + rec_id;
    }

  template< typename TText, typename TSpec >
      inline typename Records< seqan::StringSet< TText, TSpec > >::TId
    position_to_id( const Records< seqan::StringSet< TText, TSpec > >& records,
        typename Records< seqan::StringSet< TText, TSpec > >::TStringSetPosition const& pos )
    {
      return position_to_id( records, pos.i1 );
    }

  template< typename TText >
      inline typename Records< seqan::StringSet< TText, grem::Dependent > >::TPosition
    position_to_offset( const Records< seqan::StringSet< TText, grem::Dependent > >& records,
        typename Records< seqan::StringSet< TText, grem::Dependent > >::TStringSetPosition const& pos )
    {
      using seqan::length;
      if ( pos.i2 >= length( records.str[pos.i1] ) || pos.i2 < 0 ) {
        throw std::runtime_error( "position out of range" );
      }
      return pos.i2;
    }

  template< typename TText >
      inline typename Records< seqan::StringSet< TText, seqan::Owner<> > >::TPosition
    position_to_offset( const Records< seqan::StringSet< TText, seqan::Owner<> > >& records,
        typename Records< seqan::StringSet< TText, seqan::Owner<> > >::TStringSetPosition const& pos )
    {
      using seqan::length;
      if ( pos.i2 >= length( records.str[pos.i1] ) || pos.i2 < 0 ) {
        throw std::runtime_error( "position out of range" );
      }
      return records.position_to_offset( pos );
    }

  template< typename TText >
      inline void
    clear( Records< seqan::StringSet< TText, seqan::Owner<> > >& records )
    {
      records.clear();
    }

  template< typename TText >
      inline void
    clear( Records< seqan::StringSet< TText, grem::Dependent > >& records )
    {
      clear( records.str );
      records.rec_offset = 0;
      records.o_str = nullptr;
    }

  template< typename TText, typename TStringSetSpec >
      inline typename Records< seqan::StringSet< TText, TStringSetSpec > >::TSize
    length( const Records< seqan::StringSet< TText, TStringSetSpec > >& records )
    {
      return length( records.str );
    }

  template< typename TText, typename TStringSetSpec, typename TPosition >
      inline typename seqan::Reference< seqan::StringSet< TText, TStringSetSpec > const >::Type
    get_value( const Records< seqan::StringSet< TText, TStringSetSpec > >& records,
        TPosition pos )
    {
      return records.str[pos];
    }

  template< typename TText, typename TStringSetSpec, typename TPosition >
      inline typename seqan::Reference< seqan::StringSet< TText, TStringSetSpec > >::Type
    get_value( Records< seqan::StringSet< TText, TStringSetSpec > >& records,
        TPosition pos )
    {
      return records.str[pos];
    }

  template< typename TText >
      inline bool
    load_chunk( Records< seqan::StringSet< TText, grem::Dependent > >& records,
        const Records< seqan::StringSet< TText, seqan::Owner<> > >& ref,
        typename Records< seqan::StringSet< TText, seqan::Owner<> > >::TPosition n,
        typename Records< seqan::StringSet< TText, seqan::Owner<> > >::TPosition start_pos )
    {
      if ( start_pos >= length( ref.str ) || start_pos < 0 ) return false;

      clear( records.str );
      records.rec_offset = start_pos;
      records.o_str = &ref.str;
      auto i = start_pos;
      for ( ; i < n + start_pos && i < length( ref.str ); ++i ) {
        appendValue( records.str, ref.str[ i ] );
      }
      return true;
    }

  template< typename TText >
    class Records< seqan::StringSet< TText, seqan::Owner<> > > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef seqan::Owner<> TSpec;
        typedef seqan::StringSet< TText, TSpec > TStringSet;
        typedef typename MakeOwner< TStringSet >::Type TRefStringSet;
        typedef typename seqan::StringSetPosition< TStringSet >::Type TStringSetPosition;
        typedef typename seqan::Position< TStringSet >::Type TPosition;
        typedef typename seqan::Id< TStringSet >::Type TId;
        typedef typename seqan::Size< TStringSet >::Type TSize;
        /* ====================  CLASSES       ======================================= */
        /**
         *  @brief  Map seeds to their location in the reads set.
         *
         *  This class contains some data structure to map a position in the seeds set
         *  to its original position in the reads set.
         */
        class SeedMap {
          public:
            /* ====================  TYPEDEFS      =================================== */
            using bv_type = sdsl::bit_vector;
            using rank_type = bv_type::rank_1_type;
            using select_type = bv_type::select_1_type;
            using id_type = Records::TId;
            using offset_type = Records::TPosition;
            using pos_type = Records::TStringSetPosition;
            /* ====================  LIFECYCLE     =================================== */
            SeedMap( bv_type _bv, unsigned int st )
              : step( st )
            {
              this->bv.swap( _bv );
              this->initialize();
            }

            SeedMap( const SeedMap& other )
              : bv( other.bv ), step( other.step )
            {
              this->initialize();
            }

            SeedMap( SeedMap&& other ) noexcept
              : step( other.step )
            {
              this->bv.swap( other.bv );
              this->initialize();
            }

            SeedMap& operator=( const SeedMap& other )
            {
              this->bv = other.bv;
              this->step = other.step;
              this->initialize();
            }

            SeedMap& operator=( SeedMap&& other ) noexcept
            {
              this->bv.swap( other.bv );
              this->step = other.step;
              this->initialize();
            }

            ~SeedMap( ) = default;
            /* ====================  METHODS       =================================== */
              inline void
            initialize( )
            {
              sdsl::util::init_support( this->rs, &this->bv );
              sdsl::util::init_support( this->ss, &this->bv );
            }

              inline id_type
            get_reads_id( id_type seeds_id ) const
            {
              return this->rs( seeds_id );
            }

              inline offset_type
            get_reads_offset( pos_type seeds_pos ) const
            {
              id_type rid = this->get_reads_id( seeds_pos.i1 );
              id_type first_seed_id = rid ? this->ss( rid )+1 : 0;
              return ( seeds_pos.i1 - first_seed_id ) * this->step + seeds_pos.i2;
            }
          private:
            /* ====================  DATA MEMBERS  =================================== */
            bv_type bv;
            rank_type rs;
            select_type ss;
            unsigned int step;
        };
        /* ====================  DATA MEMBERS  ======================================= */
        CharStringSet<> name;
        //TStringSet2 comment;
        TStringSet str;
        //TStringSet2 qual;
        /* ====================  LIFECYCLE     ======================================= */
        Records( TId roff=0 ) : rec_offset( roff ) { }
        /* ====================  OPERATORS     ======================================= */
          inline typename seqan::Reference< TStringSet const >::Type
        operator[]( TPosition pos ) const { return get_value( *this, pos ); }
          inline typename seqan::Reference< TStringSet >::Type
        operator[]( TPosition pos ) { return get_value( *this, pos ); }
        /* ====================  ACCESSORS     ======================================= */
          inline bool
        has_seedmap( ) const
        {
          return this->sm_ptr != nullptr;
        }

          inline TId
        get_record_offset( ) const
        {
          return this->rec_offset;
        }
        /* ====================  MUTATORS      ======================================= */
          inline void
        set_record_offset( TId value )
        {
          this->rec_offset = value;
        }

          inline void
        add_record_offset( TId value )
        {
          this->rec_offset += value;
        }

          inline void
        set_seedmap( typename SeedMap::bv_type bv, unsigned int step )
        {
          this->sm_ptr = std::make_unique< SeedMap >( std::move( bv ), step );
        }
        /* ====================  METHODS       ======================================= */
          inline void
        clear( )
        {
          using seqan::clear;
          using grem::clear;
          clear( this->name );
          //clear( records.comment );
          clear( this->str );
          //clear( records.qual );
          this->set_record_offset( 0 );
          this->sm_ptr.reset( nullptr );
        }

          inline TId
        position_to_id( TId rec_id ) const
        {
          if ( this->has_seedmap( ) ) rec_id = this->sm_ptr->get_reads_id( rec_id );
          return this->rec_offset + rec_id;
        }

          inline TPosition
        position_to_offset( TStringSetPosition pos ) const
        {
          if ( this->has_seedmap( ) ) pos.i2 = this->sm_ptr->get_reads_offset( pos );
          return pos.i2;
        }
      protected:
        /* ====================  DATA MEMBERS  ======================================= */
        TId rec_offset;
        std::unique_ptr< SeedMap > sm_ptr;
    };

  template< typename TText >
    class Records< seqan::StringSet< TText, grem::Dependent > > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef grem::Dependent TSpec;
        typedef seqan::StringSet< TText, TSpec > TStringSet;
        typedef typename MakeOwner< TStringSet >::Type TRefStringSet;
        typedef typename seqan::StringSetPosition< TStringSet >::Type TStringSetPosition;
        typedef typename seqan::Position< TStringSet >::Type TPosition;
        typedef typename seqan::Id< TStringSet >::Type TId;
        typedef typename seqan::Size< TStringSet >::Type TSize;
        /* ====================  DATA MEMBERS  ======================================= */
        TStringSet str;
        /* ====================  LIFECYCLE     ======================================= */
        Records( TId roff=0 ) : rec_offset( roff ), o_str( nullptr ) { }
        /* ====================  METHODS       ======================================= */
          inline typename seqan::Reference< TStringSet const >::Type
        operator[]( TPosition pos ) const { return get_value( *this, pos ); }
          inline typename seqan::Reference< TStringSet >::Type
        operator[]( TPosition pos ) { return get_value( *this, pos ); }
        /* ====================  INTERFACE FUNCTIONS  ================================ */
          friend TId
        position_to_id< TText >( const Records& records, TId pos );

          friend void
        clear< TText >( Records& records );

          friend bool
        load_chunk< TText >( Records& records,
            const Records< TRefStringSet >& ref,
            typename Records< TRefStringSet >::TPosition n,
            typename Records< TRefStringSet >::TPosition start_pos );
      protected:
        /* ====================  DATA MEMBERS  ======================================= */
        TId rec_offset;  /**< @brief First string ID: id(i) = offset + pos(i). */
        const TRefStringSet* o_str;  /**< @brief original string set. */
    };

  /* Sequence directions strategies */
  struct ForwardStrategy;
  struct ReversedStrategy;

  /* Sequence direction tags */
  typedef seqan::Tag< ForwardStrategy > Forward;
  typedef seqan::Tag< ReversedStrategy > Reversed;

  /**
   *  @brief  Meta-function getting direction of a set of paths.
   */
  template< typename TSpec >
    struct Direction;

  /**
   *  @brief  Meta-function getting direction of paths in a `StringSet`.
   */
  template< typename TText, typename TSpec >
    struct Direction< seqan::StringSet< TText, TSpec > > {
      typedef Forward Type;
    };

  /* Seeding strategies */
  struct OverlapStrategy;
  struct GreedyOverlapStrategy;
  struct NonOverlapStrategy;
  struct GreedyNonOverlapStrategy;

  /* Seeding strategy tags */
  typedef seqan::Tag< OverlapStrategy > Overlapping;
  typedef seqan::Tag< GreedyOverlapStrategy > GreedyOverlapping;
  typedef seqan::Tag< NonOverlapStrategy > NonOverlapping;
  typedef seqan::Tag< GreedyNonOverlapStrategy > GreedyNonOverlapping;

  /* _RecordsIterBase forward declaration  --------------------------------------- */
  template< typename TRecords, typename TSpec >
    class _RecordsIterBase;
  /* END OF _RecordsIterBase forward declaration  -------------------------------- */

  /* _RecordsIterBase interface functions  --------------------------------------- */
  template< typename TRecords, typename TSpec >
      inline bool
    at_end( const _RecordsIterBase< TRecords, TSpec >& iter );
  template< typename TRecords, typename TSpec >
      inline typename _RecordsIterBase< TRecords, TSpec >::TId const&
    get_id( const _RecordsIterBase< TRecords, TSpec >& iter );
  template< typename TRecords, typename TSpec >
      inline typename _RecordsIterBase< TRecords, TSpec >::TTextPosition const&
    get_offset( const _RecordsIterBase< TRecords, TSpec >& iter );
  template< typename TRecords, typename TSpec >
      inline typename _RecordsIterBase< TRecords, TSpec >::TStringSetPosition const&
    get_position( const _RecordsIterBase< TRecords, TSpec >& iter );
  /* END OF _RecordsIterBase interface functions  -------------------------------- */

  template< typename TStringSet, typename TSpec >
    class _RecordsIterBase< Records< TStringSet >, TSpec > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef Records< TStringSet > TRecords;
        typedef typename seqan::StringSetPosition< TStringSet >::Type TStringSetPosition;
        typedef typename seqan::Value< TStringSet >::Type TText;
        typedef typename seqan::Id< TStringSet >::Type TId;
        typedef typename seqan::Position< TText >::Type TTextPosition;
        typedef typename seqan::Infix< TText const >::Type TInfix;
        /* ====================  LIFECYCLE     ======================================= */
        _RecordsIterBase( const TRecords* recs, TTextPosition len )
          : records( recs ), current_pos( { 0, 0 } ), infix_len( len ) { }

        _RecordsIterBase( const _RecordsIterBase& ) = default;
        _RecordsIterBase( _RecordsIterBase&& ) = default;
        _RecordsIterBase& operator=( const _RecordsIterBase& ) = default;
        _RecordsIterBase& operator=( _RecordsIterBase&& ) = default;
        ~_RecordsIterBase() = default;
        /* ====================  ACCESSORS     ======================================= */
          inline const TRecords*
        get_records_ptr( ) const
        {
          return this->records;
        }
        /* ====================  OPERATORS     ======================================= */
          inline TInfix
        operator*() const {
#ifndef NDEBUG
          if ( at_end( *this ) ) {
            throw std::range_error( "Iterator has already reached at the end." );
          }
#endif  /* ----- #ifndef NDEBUG  ----- */
          return infixWithLength(
              ( *this->records )[ this->current_pos.i1 ],
              this->current_pos.i2,
              this->infix_len );
        }
        /* ====================  INTERFACE FUNCTIONS  ================================ */
        friend bool
          at_end< TRecords, TSpec >( const _RecordsIterBase< TRecords, TSpec >& );
        friend const TId&
          get_id< TRecords, TSpec >( const _RecordsIterBase< TRecords, TSpec >& );
        friend const TTextPosition&
          get_offset< TRecords, TSpec >( const _RecordsIterBase< TRecords, TSpec >& );
        friend const TStringSetPosition&
          get_position< TRecords, TSpec >( const _RecordsIterBase< TRecords, TSpec >& );
      protected:
        const TRecords* records;
        TStringSetPosition current_pos;
        TTextPosition infix_len;
    };

  template< typename TRecords, typename TSpec >
      inline bool
    at_end( const _RecordsIterBase< TRecords, TSpec >& iter )
    {
      if ( iter.current_pos.i1 >= length( *iter.records ) ) return true;
      return false;
    }

  template< typename TRecords, typename TSpec >
      inline typename _RecordsIterBase< TRecords, TSpec >::TId const&
    get_id( const _RecordsIterBase< TRecords, TSpec >& iter )
    {
      return iter.current_pos.i1;
    }

  template< typename TRecords, typename TSpec >
      inline typename _RecordsIterBase< TRecords, TSpec >::TTextPosition const&
    get_offset( const _RecordsIterBase< TRecords, TSpec >& iter )
    {
      return iter.current_pos.i2;
    }

  template< typename TRecords, typename TSpec >
      inline typename _RecordsIterBase< TRecords, TSpec >::TStringSetPosition const&
    get_position( const _RecordsIterBase< TRecords, TSpec >& iter )
    {
      return iter.current_pos;
    }

  template< typename TRecords, typename TSpec >
    class RecordsIter;

  template< typename TRecords >
    class RecordsIter< TRecords, Overlapping >
    : public _RecordsIterBase< TRecords, Overlapping > {
      private:
        /* ====================  TYPEDEFS      ======================================= */
        typedef _RecordsIterBase< TRecords, Overlapping > TBase;
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename TRecords::TStringSet TStringSet;
        typedef typename TBase::TStringSetPosition TStringSetPosition;
        typedef typename TBase::TText TText;
        typedef typename TBase::TId TId;
        typedef typename TBase::TTextPosition TTextPosition;
        typedef typename TBase::TInfix TInfix;
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        TTextPosition step;
      public:
        /* ====================  LIFECYCLE     ======================================= */
        RecordsIter( const TRecords* recs, TTextPosition len, TTextPosition stp )
          : TBase( recs, len ), step( stp ) { }
        /* ====================  OPERATORS     ======================================= */
          inline RecordsIter&
        operator++()
        {
#ifndef NDEBUG
          if ( at_end( *this ) ) {
            throw std::range_error( "Iterator has already reached at the end." );
          }
#endif  /* ----- #ifndef NDEBUG  ----- */
          const auto& current_strlen = length( ( *this->records )[ this->current_pos.i1 ] );
          if ( this->current_pos.i2 + this->infix_len + this->step <= current_strlen ) {
            this->current_pos.i2 += this->step;
          }
          else {
            this->current_pos.i2 = 0;
            ++this->current_pos.i1;
          }
          return *this;
        }
    };

  /**
   *  @brief  NonOverlapping Records iterator.
   *
   *  It is defined as "Overlapping" with step size equals to infix length.
   */
  template< typename TRecords >
    class RecordsIter< TRecords, NonOverlapping >
    : public RecordsIter< TRecords, Overlapping > {
      private:
        /* ====================  TYPEDEFS      ======================================= */
        typedef RecordsIter< TRecords, Overlapping > TBase;
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename TRecords::TStringSet TStringSet;
        typedef typename TBase::TStringSetPosition TStringSetPosition;
        typedef typename TBase::TText TText;
        typedef typename TBase::TId TId;
        typedef typename TBase::TTextPosition TTextPosition;
        typedef typename TBase::TInfix TInfix;
        /* ====================  LIFECYCLE     ======================================= */
        RecordsIter( const TRecords* recs, TTextPosition len )
          : TBase( recs, len, len ) { }
    };

  /**
   *  @brief  GreedyOverlapping Records iterator.
   *
   *  It is defined as "Overlapping" with step size equals to 1.
   */
  template< typename TRecords >
    class RecordsIter< TRecords, GreedyOverlapping >
    : public RecordsIter< TRecords, Overlapping > {
      private:
        /* ====================  TYPEDEFS      ======================================= */
        typedef RecordsIter< TRecords, Overlapping > TBase;
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename TRecords::TStringSet TStringSet;
        typedef typename TBase::TStringSetPosition TStringSetPosition;
        typedef typename TBase::TText TText;
        typedef typename TBase::TId TId;
        typedef typename TBase::TTextPosition TTextPosition;
        typedef typename TBase::TInfix TInfix;
        /* ====================  LIFECYCLE     ======================================= */
        RecordsIter( const TRecords* recs, TTextPosition len )
          : TBase( recs, len, 1 ) { }
    };

  /**
   *  @overload postfix increment operator.
   */
  template< typename TRecords, typename TSpec >
      inline RecordsIter< TRecords, TSpec >
    operator++( RecordsIter< TRecords, TSpec >& iter, int )
    {
      RecordsIter< TRecords, TSpec > tmp = iter;
      iter.operator++();
      return tmp;
    }

  /* END OF Data structures  ----------------------------------------------------- */

  /* Interface functions  -------------------------------------------------------- */

  /**
   *  @brief  Read records from the input file into a sequence record set.
   *
   *  @param[out]  records Sequence record set to store records in the input file.
   *  @param[in,out]  infile The input file.
   *  @param[in]  num_record Read this number of record from the input file.
   *
   *  A wrapper function for `seqan::readRecords` method to read the records into
   *  sequence record set. If `num_record` is equal to zero, it reads all recrods.
   */
  template< typename TText >
      inline void
    readRecords( Records< seqan::StringSet< TText, seqan::Owner<> > >& records,
        seqan::SeqFileIn& infile,
        unsigned int num_record = 0 )
    {
      CharStringSet<> quals;
      if ( num_record != 0 ) {
        seqan::readRecords( records.name, records.str, quals, infile, num_record );
      }
      else {
        seqan::readRecords( records.name, records.str, quals, infile );
      }
      assignQualities( records.str, quals );
      return;
    }  /* -----  end of template function readRecords  ----- */

  template< typename TText >
      inline void
    readRecords( Records< seqan::StringSet< TText, seqan::Owner<> > >& records,
        klibpp::SeqStreamIn& iss,
        unsigned int num_record=0 )
    {
      klibpp::KSeq rec;
      clear( records );
      records.set_record_offset( iss.counts() );
      unsigned int i = 0;
      while ( iss >> rec ) {
        appendValue( records.name, rec.name );
        appendValue( records.str, rec.seq );
        if ( ++i == num_record ) break;
      }
    }

  /**
   *  @brief  Get next lexicographical k-mer in a specific position in the string.
   *
   *  @param  str The input k-mer.
   *  @param  rank The rank of the character which should be incremented (1-based).
   *  @return The lowest rank at which the character is not modified. In case that
   *          there is no k-mer it returns negative value of -1.
   *
   *  The value of input string would be changed to next lexicographical k-mer by
   *  incrementing the character at position '`rank` - 1', the carry over may change
   *  the other characters at lower ranks. The lowest rank at which the character is not
   *  modified will be returned.
   */
  template< typename TText >
      inline int
    increment_kmer( TText& str, unsigned long int rank = 0 )
    {
      static const unsigned int max_value = ordValue( maxValue( str[0] ) );
      static const unsigned int min_value = ordValue( minValue( str[0] ) );

      int _rank = std::min( rank - 1, length( str ) - 1 );

      while ( 0 <= _rank && ordValue( str[_rank] ) == max_value ) {
        str[_rank] = min_value;
        _rank--;
      }
      if ( _rank >= 0 ) {
        str[_rank] = ordValue( str[_rank] ) + 1;
        return _rank;
      }
      for( unsigned int i = 0; i < length( str ); ++i ) str[i] = max_value;
      return -1;
    }

  /**
   *  @brief  Add any k-mers from the given string set with `step` distance to seed set.
   *
   *  @param  seeds The seed set.
   *  @param  string_set The string set from which seeds are extracted.
   *  @param  k The length of the seeds.
   *  @param  step The step size.
   *
   *  For each string in string set, it add all substring of length `k` starting from 0
   *  to end of string with `step` distance with each other. If `step` is equal to `k`,
   *  it gets non-overlapping substrings of length k.
   */
  template< typename TText, typename TStringSetSpec >
      inline void
    seeding( seqan::StringSet< TText, seqan::Owner<> >& seeds,
        const seqan::StringSet< TText, TStringSetSpec >& string_set,
        unsigned int k,
        unsigned int step,
        sdsl::bit_vector* bv_ptr=nullptr )
    {
      typedef typename seqan::Size< seqan::StringSet< TText, TStringSetSpec > >::Type size_type;
      typedef typename seqan::Position< TText >::Type pos_type;

      using seqan::length;
      clear( seeds );
      // The total number of seeds is always less than: (len(R) - |R|k)/s + |R|;
      // where len(R) is the total sequence length of reads set R, and |R| is the number
      // of reads in R, k is seed length, and s is step size.
      auto lensum = lengthSum( string_set );
      auto nofreads = length( string_set );
      assert( lensum >= nofreads*k );
      auto est_nofseeds = static_cast<int>( ( lensum - nofreads*k ) / step ) + nofreads;
      reserve( seeds, est_nofseeds );
      if ( bv_ptr ) sdsl::util::assign( *bv_ptr, sdsl::bit_vector( est_nofseeds, 0 ) );

      for ( size_type idx = 0; idx < length( string_set ); ++idx ) {
        for ( pos_type i = 0; i < length( string_set[idx] ) - k + 1; i += step ) {
          appendValue( seeds, seqan::infixWithLength( string_set[idx], i, k ) );
        }
        if ( bv_ptr ) ( *bv_ptr )[ length( seeds ) - 1 ] = 1;
      }
      if ( bv_ptr ) ( *bv_ptr ).resize( length( seeds ) );
    }  /* -----  end of template function seeding  ----- */

  /**
   *  @brief  Add any k-mers from the given records with `step` distance to seeds record.
   *
   *  @param  seeds The seeds record.
   *  @param  reads The string set from which seeds are extracted.
   *  @param  k The length of the seeds.
   *  @param  step The step size.
   *
   *  For each string in reads record, it add all substring of length `k` starting from
   *  0 to end of string with `step` distance with each other. If `step` is equal to
   *  `k`, it gets non-overlapping substrings of length k.
   */
  template< typename TRecords1, typename TRecords2,
    typename = std::enable_if_t< std::is_same< typename TRecords1::TSpec, seqan::Owner<> >::value, void > >
      inline void
    seeding( TRecords1& seeds,
        TRecords2 const& reads,
        unsigned int k,
        unsigned int step )
    {
      clear( seeds );
      sdsl::bit_vector bv;
      seeding( seeds.str, reads.str, k, step, &bv );
      seeds.set_seedmap( std::move( bv ), step );
      seeds.set_record_offset( reads.get_record_offset() );
    }

  /**
   *  @brief  Seeding a set of sequence by reporting overlapping k-mers.
   *
   *  @param  seeds The resulting set of strings containing seeds.
   *  @param  string_set The string set from which seeds are extracted.
   *  @param  k The length of the seeds.
   *  @param  tag Tag for greedy fixed-length overlapping seeding strategy.
   *
   *  Extract a set of overlapping seeds of length k.
   */
  template< typename T1, typename T2,
    typename = std::enable_if_t< std::is_same< typename Ownership< T1 >::Type, seqan::Owner<> >::value, void > >
      inline void
    seeding( T1& seeds, const T2& string_set, unsigned int k, GreedyOverlapping )
    {
      seeding( seeds, string_set, k, 1 );
    }  /* -----  end of template function seeding  ----- */

  /**
   *  @brief  Seeding by partitioning each sequence into non-overlapping k-mers.
   *
   *  @param  seeds The resulting set of strings containing seeds.
   *  @param  string_set The string set from which seeds are extracted.
   *  @param  k The length of the seeds.
   *  @param  tag Tag for fixed-length non-overlapping seeding strategy.
   *
   *  Extract a set of non-overlapping seeds of length k.
   */
  template< typename T1, typename T2,
    typename = std::enable_if_t< std::is_same< typename Ownership< T1 >::Type, seqan::Owner<> >::value, void > >
      inline void
    seeding( T1& seeds, const T2& string_set, unsigned int k, NonOverlapping )
    {
      seeding( seeds, string_set, k, k );
    }  /* -----  end of function seeding  ----- */

  /**
   *  @brief  Seeding by greedy partitioning each sequence into non-overlapping k-mers.
   *
   *  @param  seeds The resulting set of strings containing seeds.
   *  @param  string_set The string set from which seeds are extracted.
   *  @param  k The length of the seeds.
   *  @param  tag Tag for greedy fixed-length non-overlapping seeding strategy.
   *
   *  Extract a set of non-overlapping seeds of length k.
   *
   *  NOTE: In case that the length of sequence is not dividable by k the last seed
   *        may overlap its previous (greedy).
   */
  template< typename TText, typename TStringSetSpec >
      inline void
    seeding( seqan::StringSet< TText, seqan::Owner<> >& seeds,
        const seqan::StringSet< TText, TStringSetSpec >& string_set,
        unsigned int k,
        GreedyNonOverlapping )
    {
      typedef typename seqan::Size< seqan::StringSet< TText, TStringSetSpec > >::Type size_type;
      typedef typename seqan::Position< TText >::Type pos_type;

      clear( seeds );
      reserve( seeds, static_cast<int>( lengthSum( string_set ) / k ) );

      for ( size_type idx = 0; idx < length( string_set ); ++idx ) {
        for ( pos_type i = 0; i < length( string_set[idx] ) - k; i += k ) {
          appendValue( seeds, infixWithLength( string_set[idx], i, k ) );
        }
        pos_type last = length( string_set[idx] ) - k;
        appendValue( seeds, infixWithLength( string_set[idx], last, k ) );
      }
    }  /* -----  end of function seeding  ----- */

  /* END OF Interface functions  ------------------------------------------------- */
}  /* -----  end of namespace grem  ----- */

namespace seqan {
  template< typename TStringSet, typename TSpec >
    struct Iterator< grem::Records< TStringSet >, TSpec > {
      typedef grem::RecordsIter< grem::Records< TStringSet >, TSpec > Type;
    };
}  /* -----  end of namespace seqan  ----- */

#endif  /* ----- #ifndef SEQUENCE_H__  ----- */

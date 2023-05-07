/**
 *    @file  ggsim.hpp
 *   @brief  Graph genome haplotype and reads simulator header file
 *
 *  This header file defines constants, traits, data structures, and utility functions
 *  for `ggsim`.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Mon Jan 18, 2021  22:17
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  PSI_TOOLS_GGSIM_HPP__
#define  PSI_TOOLS_GGSIM_HPP__

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>

#include <cxxopts.hpp>
#include <kseq++/kseq++.hpp>
#include <gum/graph.hpp>
#include <gum/io_utils.hpp>
#include <psi/graph.hpp>
#include <psi/graph_iter.hpp>
#include <psi/pathindex.hpp>
#include <psi/utils.hpp>

#include "vg/vg.pb.h"
#include "vg/stream.hpp"


using namespace psi;

/* ====== Constants ====== */

constexpr const char* LONG_DESC = "Simulate haplotypes or reads from a genome graph";
constexpr const char DEFAULT_QUAL_SCORE = 'I';
constexpr const int MAX_TRIES = 100;
constexpr const char READ_COMMENT_DELIMITER = ';';
constexpr const int READ_NAME_LENGTH = 16;
// Default values for command line arguments
constexpr const char* DEFAULT_RNDSEED = "0";
constexpr const char* DEFAULT_OUTPUT = "-";  // stdout
constexpr const char* DEFAULT_PLOIDY = "2";
constexpr const char* DEFAULT_READLEN = "0";
constexpr const char* DEFAULT_NUMREADS = "0";
constexpr const char* DEFAULT_ERRRATE = "0.0";
constexpr const char* DEFAULT_INDRATE = "0.0";
constexpr const char* DEFAULT_FORWARD = "false";
constexpr const char* DEFAULT_ALLOWNS = "false";
constexpr const char* DEFAULT_DISTANCE = "0";
constexpr const char* DEFAULT_DEVIATION = "50";

struct Parameters {
  std::string output;
  unsigned int ploidy;
  unsigned int readlen;
  unsigned long int numreads;
  double errorrate;
  double indelrate;
  unsigned int distance;
  unsigned int sd;
  unsigned int seed;
  bool fwd;
  bool allow_ns;
};

namespace fmt {
  struct Fasta {
    constexpr static const char* extension = ".fasta";
    constexpr static const char* short_extension = ".fa";
    constexpr static const char* extension_repr = ".fasta'/'.fa";
    constexpr static const char* type_string = "fasta";
    constexpr static unsigned char type_code = 1;
  };

  struct Fastq {
    constexpr static const char* extension = ".fastq";
    constexpr static const char* short_extension = ".fq";
    constexpr static const char* extension_repr = ".fastq'/'.fq";
    constexpr static const char* type_string = "fastq";
    constexpr static unsigned char type_code = 2;
  };

  struct Seq {
    constexpr static const char* extension = ".txt";
    constexpr static const char* short_extension = ".txt";
    constexpr static const char* extension_repr = ".txt";
    constexpr static const char* type_string = "plain";
    constexpr static unsigned char type_code = 3;
  };

  struct Gam {
    constexpr static const char* extension = ".gam";
    constexpr static const char* short_extension = ".gam";
    constexpr static const char* extension_repr = ".gam";
    constexpr static const char* type_string = "gam";
    constexpr static unsigned char type_code = 4;
  };

  class Type {
    public:
      Type( ) : code( 0 )
      { }

      Type( std::string const& type_str )
      {
        this->set( type_str );
      }

        void
      set( std::string const& type_str )
      {
        if ( type_str == Fasta::type_string ) {
          this->code = Fasta::type_code;
        }
        else if ( type_str == Fastq::type_string ) {
          this->code = Fastq::type_code;
        }
        else if ( type_str == Seq::type_string ) {
          this->code = Seq::type_code;
        }
        else if ( type_str == Gam::type_string ) {
          this->code = Gam::type_code;
        }
        else throw std::runtime_error( "Internal error: unknown type string" );
      }

        inline
      operator bool( ) const
      {
        return this->code != 0;
      }

      template< typename TFormat >
          inline bool
        operator==( TFormat ) const
        {
          return this->code == TFormat::type_code;
        }

        inline bool
      operator==( Type const& t ) const
      {
        return this->code == t.code;
      }
    private:
      unsigned char code;
  };

  template< typename TFormat >
      inline bool
    operator==( TFormat, Type const& t )
    {
      return t.operator==( TFormat() );
    }

    std::istream&
  operator>>( std::istream& is, Type& t )
  {
    std::string type_str;
    is >> type_str;
    t.set( type_str );
    return is;
  }

  template< typename TFormat >
      inline static bool
    check_extension( std::string const& filename, TFormat )
    {
      return ends_with( filename, TFormat::extension ) ||
        ends_with( filename, TFormat::short_extension );
    }

    inline fmt::Type
  get_type( std::string const& output )
  {
    if ( fmt::check_extension( output, fmt::Seq() ) ) {
      return fmt::Type( "plain" );
    }
    if ( fmt::check_extension( output, fmt::Fastq() ) ) {
      return fmt::Type( "fastq" );
    }
    if ( fmt::check_extension( output, fmt::Fasta() ) ) {
      return fmt::Type( "fasta" );
    }
    if ( fmt::check_extension( output, fmt::Gam() ) ) {
      return fmt::Type( "gam" );
    }

    std::string msg = std::string( "Output file extension must be either '" ) +
        fmt::Fasta::extension_repr + "', '" +
        fmt::Gam::extension_repr + "', '" +
        fmt::Fastq::extension_repr + "', or '" +
        fmt::Seq::extension_repr + "'";
    throw cxxopts::OptionParseException( msg );
  }

  template< typename TParseResult >
      inline fmt::Type
    get_type( TParseResult const& parsed )
    {
      fmt::Type type;
      if ( parsed.count( "type" ) ) type = parsed[ "type" ].template as< fmt::Type >();
      if ( type ) return type;
      return get_type( parsed[ "output" ].template as< std::string >() );
    }
}  /* -----  end of namespace format  ----- */

namespace event {
  enum Event { match=0, mismatch, insertion, deletion };
}

namespace rnd {
  thread_local static std::mt19937 lgen;
  thread_local static unsigned int lseed = std::mt19937::default_seed;
  std::atomic_uint iseed = 0;

  inline void
  init_gen( unsigned int seed=0 )
  {
    iseed.store( seed );
    if ( seed != 0 && seed != lseed ) {
      lseed = seed;
      lgen.seed( seed );
    }
  }

  inline std::mt19937&
  get_gen( )
  {
    if ( iseed.load() == 0 ) return psi::random::gen;
    else return lgen;
  }
}  /* -----  end of namespace rnd  ----- */

/* ====== Read type tags ====== */

struct SingleEnd {};
struct PairedEnd {};

/* ====== Writer ====== */

template< typename TType >
class Writer {
public:
  /* === TYPE MEMBERS === */
  typedef klibpp::KSeq value_type;
  typedef std::size_t size_type;
  /* === LIFECYCLE === */
  Writer( std::string const& output )
  {
    if ( output != "-" ) this->ost = std::make_unique< klibpp::SeqStreamOut >( output.c_str() );
    else this->ost = std::make_unique< klibpp::SeqStreamOut >( fileno( stdout ) );
    if ( std::is_same< TType, fmt::Fastq >::value ) *this->ost << klibpp::format::fastq;
    else if ( std::is_same< TType, fmt::Fasta >::value ) *this->ost << klibpp::format::fasta;
    else assert( false );
  }
  /* === METHODS === */
  Writer&
  operator<<( value_type const& record )
  {
    *this->ost << record;
    return *this;
  }
private:
  /* === DATA MEMBERS === */
  std::unique_ptr< klibpp::SeqStreamOut > ost;
};

template<>
class Writer< fmt::Seq > {
public:
  /* === TYPE MEMBERS === */
  typedef klibpp::KSeq value_type;
  typedef std::size_t size_type;
  /* === LIFECYCLE === */
  Writer( std::string const& output )
    : ost( nullptr )
  {
    if ( output == "-" ) this->ost.rdbuf( std::cout.rdbuf() );
    else {
      this->ofs.open( output, std::ofstream::out | std::ofstream::binary );
      this->ost.rdbuf( ofs.rdbuf() );
    }
  }
  /* === METHODS === */
  Writer&
  operator<<( value_type const& record )
  {
    this->ost << record.seq << std::endl;
    return *this;
  }
private:
  /* === DATA MEMBERS === */
  std::ofstream ofs;
  std::ostream ost;
};

template<>
class Writer< fmt::Gam > {
public:
  /* === TYPE MEMBERS === */
  typedef vg::Alignment value_type;
  typedef std::size_t size_type;
  /* === STATIC MEMBERS === */
  static constexpr const size_type DEFAULT_BUFFER_SIZE = 128;
  /* === LIFECYCLE === */
  Writer( std::string const& output, size_type buffer_size=DEFAULT_BUFFER_SIZE )
    : ost( nullptr )
  {
    if ( output == "-" ) this->ost.rdbuf( std::cout.rdbuf() );
    else {
      this->ofs.open( output, std::ofstream::out | std::ofstream::binary );
      this->ost.rdbuf( ofs.rdbuf() );
    }
  }

  ~Writer( ) noexcept
  {
    if ( !this->buffer.empty() ) write();
  }
  /* === METHODS === */
  Writer&
  operator<<( value_type alignment )
  {
    this->buffer.push_back( std::move( alignment ) );
    if ( this->buffer.size() >= this->buffer_size ) write();
    return *this;
  }
private:
  /* === DATA MEMBERS === */
  std::ofstream ofs;
  std::ostream ost;
  size_type buffer_size;
  std::vector< value_type > buffer;
  /* === METHODS === */
  inline void
  write()
  {
    std::function< value_type&( uint64_t ) > access = [this]( size_type i ) -> value_type& { return this->buffer[i]; };
    stream::write< value_type& >( this->ost, this->buffer.size(), access );
    this->buffer.clear();
  }
};

/* ====== Helper functions ====== */

template< typename TGraph >
inline void
_to( klibpp::KSeq& record, Path< TGraph > const& haplotype, std::string const& name )
{
  record.name = name;
  record.seq = sequence( haplotype );
  record.qual = std::string( record.seq.size(), DEFAULT_QUAL_SCORE );
}

template< typename TGraph >
inline void
_to( vg::Path& path, Path< TGraph > const& haplotype, std::string const& name )
{
  path.set_name( name );
  TGraph const* graph_ptr = haplotype.get_graph_ptr();
  int64_t rank = 1;
  for ( auto const& id : haplotype ) {
    auto label_len = graph_ptr->node_length( id );
    auto map_ptr = path.add_mapping();
    auto p = map_ptr->mutable_position();
    p->set_node_id( graph_ptr->coordinate_id( id ) );
    p->set_offset( 0 );
    auto edit_ptr = map_ptr->add_edit();
    edit_ptr->set_from_length( label_len );
    edit_ptr->set_to_length( label_len );
    map_ptr->set_rank( rank++ );
  }
}

template< typename TGraph >
inline void
_to( vg::Alignment& aln, Path< TGraph > const& haplotype, std::string const& name )
{
  aln.set_sequence( sequence( haplotype ) );
  aln.set_name( name );
  auto path_ptr = aln.mutable_path();
  _to( *path_ptr, haplotype, name );
}

inline void
_to( klibpp::KSeq& record, klibpp::KSeq segment, vg::Path )
{
  record = std::move( segment );
}

inline void
_to( vg::Alignment& aln, klibpp::KSeq segment, vg::Path path )
{
  aln.set_sequence( std::move( segment.seq ) );
  aln.set_name( std::move( segment.name ) );
  aln.mutable_path()->Swap( &path );
}

template< typename TValueType, typename ...TArgs >
inline TValueType
to( TArgs&&... args )
{
  TValueType retval;
  _to( retval, std::forward< TArgs >( args )... );
  return retval;
}

inline void
_as( std::string& out, std::string str )
{
  out = std::move( str );
}

inline void
_as( unsigned long long int& out, std::string str )
{
  out = std::stoull( str );
}

inline void
_as( unsigned long int& out, std::string str )
{
  out = std::stoul( str );
}

inline void
_as( long long int& out, std::string str )
{
  out = std::stoll( str );
}

template< typename T >
inline T
as( std::string const& str, std::string::size_type& n )
{
  auto m = str.find( READ_COMMENT_DELIMITER, n );
  std::string::size_type count = ( m == std::string::npos ? m : m - n );
  T retval;
  _as( retval, str.substr( n, count ) );
  n = m;
  return retval;
}

#endif  /* --- #ifndef PSI_TOOLS_GGSIM_HPP__ --- */

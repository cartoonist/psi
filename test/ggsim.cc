/**
 *    @file  ggsim.cc
 *   @brief  Graph genome haplotype and reads simulator
 *
 *  Simulate haplotypes or reads from a graph genome.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Jul 03, 2018  19:09
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cxxopts/cxxopts.h>
#include <kseq++/kseq++.h>

#include "vargraph.h"


using namespace grem;

constexpr const char* LONG_DESC = "Simulate haplotypes or reads from a graph genome";
constexpr const char CHAR_BP_DELETED = '-';
constexpr const char DEFAULT_QUAL_SCORE = 'I';
constexpr int MAX_TRIES = 100;
// Default values for command line arguments
constexpr const char* DEFAULT_RNDSEED = "0";
constexpr const char* DEFAULT_OUTPUT = "-";  // stdout
constexpr const char* DEFAULT_PLOIDY = "2";
constexpr const char* DEFAULT_SUBRATE = "0.0";
constexpr const char* DEFAULT_INDRATE = "0.0";
constexpr const char* DEFAULT_FORWARD = "false";
constexpr const char* DEFAULT_ALLOWNS = "false";


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
    constexpr static const char* extension = ".seq";
    constexpr static const char* short_extension = ".seq";
    constexpr static const char* extension_repr = ".seq";
    constexpr static const char* type_string = "sequence";
    constexpr static unsigned char type_code = 3;
  };

  class Type {
    public:
      Type( ) : code( 0 )
      { }

      Type( std::string const type_str )
      {
        this->set( type_str );
      }

        void
      set( std::string const type_str )
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
      return fmt::Type( "sequence" );
    }
    if ( fmt::check_extension( output, fmt::Fastq() ) ) {
      return fmt::Type( "fastq" );
    }
    if ( fmt::check_extension( output, fmt::Fasta() ) ) {
      return fmt::Type( "fasta" );
    }

    std::string msg = std::string( "Output file extension must be either '" ) +
      fmt::Fasta::extension_repr + "' for haplotypes, or '" +
      fmt::Fastq::extension_repr + "' or'" +
      fmt::Seq::extension_repr + "' for reads";
    throw cxxopts::OptionParseException( msg );
  }

  template< typename TResult >
      inline fmt::Type
    get_type( TResult const& result )
    {
      fmt::Type type;
      if ( result.count( "type" ) ) {
        type = result[ "type" ].template as< fmt::Type >();
      }
      if ( type ) return type;
      return get_type( result[ "output" ].template as< std::string >() );
    }

    inline bool
  is_reads( fmt::Type const& type )
  {
    return type == fmt::Seq() || type == fmt::Fastq();
  }
}  /* -----  end of namespace format  ----- */

  inline bool
has_n( std::string const& seq )
{
  return seq.find( 'N' ) != std::string::npos;
}

  inline void
simulate_haplotype( VarGraph const& vargraph, unsigned int ploidy,
    unsigned int seed, std::vector< klibpp::KSeq >& hapseqs )
{
  using TIterator = seqan::Iterator< VarGraph, Haplotyper< Random > >::Type;

  hapseqs.reserve( ploidy );
  TIterator hap_itr( &vargraph );
  Path< VarGraph > haplotype( &vargraph );
  for ( std::size_t rank = 1; rank <= vargraph.max_path_rank(); ++rank ) {
    const auto& pathname = vargraph.path_name( rank );
    VarGraph::nodeid_type start = vargraph.node_at_path_position( pathname, 0 );
    go_begin( hap_itr, start, seed );
    klibpp::KSeq record;
    for ( unsigned int i = 0; i < ploidy; ++i ) {
      get_rnd_full_haplotype( haplotype, hap_itr );
      record.name = pathname + "-" + std::to_string( i + 1 );
      record.seq = sequence( haplotype );
      hapseqs.push_back( std::move( record ) );
      haplotype.clear();
      record.clear();
    }
  }
}

  inline void
impose_errors( klibpp::KSeq& haps,
    double subrate, double indelrate )
{
}

  inline void
impose_errors( std::vector< klibpp::KSeq >& haps,
    double subrate, double indelrate )
{
  if ( subrate == 0 && indelrate == 0 ) return;
  throw std::runtime_error( "imposing errors is not implemented" );
}

  inline std::string
simulate_read( std::string const& haplotype, std::size_t pos,
    unsigned int readlen, bool fwd )
{
  /* FIXME: This is assumed that pos + readlen won't exceed the sequence. */
  std::string candidate;
  candidate.reserve( readlen );
  std::size_t cursor;
  if ( fwd ) cursor = pos;
  else if ( haplotype.size() > pos + readlen ) {
    cursor = haplotype.size() - pos - readlen;
  }
  else cursor = 0;
  unsigned int len = 0;
  for ( ; len < readlen && cursor < haplotype.size(); ++cursor ) {
    if ( haplotype[ cursor ] == CHAR_BP_DELETED ) continue;
    candidate += haplotype[ cursor ];
    len++;
  }
  if ( !fwd ) {
    candidate = complement( candidate );
    std::reverse( candidate.begin(), candidate.end() );
  }
  return candidate;
}

  inline void
simulate_all_reads( std::vector< klibpp::KSeq > const& haps, unsigned int seed,
    unsigned int readlen, unsigned int numreads, bool fwd, bool allow_ns,
    std::vector< klibpp::KSeq >& seqs )
{
  std::random_device rd;         // Will be used to obtain a seed for the random no. engine
  if ( seed == 0 ) seed = rd();  // use random_device to generate a seed if seed is not provided
  std::mt19937 gen( seed );      // Standard mersenne_twister_engine seeded with seed

  seqs.reserve( numreads );
  klibpp::KSeq candidate;
  std::size_t count = 0;
  bool dir = true;
  for ( const auto& h : haps ) {
    std::size_t ubound = ( h.seq.size() > readlen ? h.seq.size() - readlen : 0 );
    std::uniform_int_distribution<> dis( 0, ubound );
    for ( unsigned int i = 0; i < ceil( numreads / haps.size() ); ++i ) {
      int tries = MAX_TRIES;
      std::size_t pos;
      do {
        pos = dis( gen );
        candidate.seq = simulate_read( h.seq, pos, readlen, fwd || dir );
      } while ( !allow_ns && tries-- && has_n( candidate.seq ) );

      if ( tries == 0 ) {
        std::cerr << "Reads may contain 'N' since nothing found after "
                  << MAX_TRIES << " attemps!" << std::endl;
      }

      candidate.name = "read-" + std::to_string( count++ );
      candidate.comment = h.name + "@" + std::to_string( pos ) + " " +
        ( fwd || dir ? "F" : "R" );
      candidate.qual = std::string( candidate.seq.size(), DEFAULT_QUAL_SCORE );
      seqs.push_back( std::move( candidate ) );
      candidate.clear();
      dir = !dir;
    }
  }
  seqs.resize( numreads );
}

  inline void
simulate( VarGraph const& vargraph, unsigned int ploidy, unsigned int seed,
    unsigned int readlen, unsigned int numreads, double subrate,
    double indelrate, bool fwd, bool allow_ns, fmt::Type type,
    std::vector< klibpp::KSeq >& seqs )
{
  std::vector< klibpp::KSeq > haps_tmp;
  std::vector< klibpp::KSeq >* haps_ptr;
  if ( is_reads( type ) ) haps_ptr = &haps_tmp;
  else haps_ptr = &seqs;
  simulate_haplotype( vargraph, ploidy, seed, *haps_ptr );
  if ( !is_reads( type ) ) return;
  impose_errors( haps_tmp, subrate, indelrate );
  simulate_all_reads( haps_tmp, seed, readlen, numreads, fwd, allow_ns, seqs );
}

template< typename TType >
    inline void
  write_output( std::string const& output,
      std::vector< klibpp::KSeq > const& seqs, TType )
  {
    int fd;
    if ( output == "-" ) fd = STDOUT_FILENO;
    else fd = open( output.c_str(), O_CREAT | O_WRONLY, 0666 );
    auto ks = klibpp::make_okstream( fd, write );
    for ( const auto& rec : seqs ) ks << rec;
    ks << klibpp::kend;
    close( fd );
  }

  inline void
write_output( std::string const& output,
    std::vector< klibpp::KSeq > const& seqs, fmt::Seq )
{
  std::streambuf* buf;
  std::ofstream ofs;
  if ( output == "-" ) {
    buf = std::cout.rdbuf();
  }
  else {
    ofs.open( output, ofstream::out | ofstream::binary );
    buf = ofs.rdbuf();
  }
  std::ostream out( buf );
  for ( const auto& rec : seqs ) out << rec.seq << endl;
}

  inline void
write_output( std::string const& output,
    std::vector< klibpp::KSeq > const& seqs, fmt::Type const& t )
{
  if ( t == fmt::Seq() ) {
    write_output( output, seqs, fmt::Seq() );
  }
  else if ( t == fmt::Fasta() ) {
    write_output( output, seqs, fmt::Fasta() );
  }
  else if ( t == fmt::Fastq() ) {
    write_output( output, seqs, fmt::Fastq() );
  }
}

  void
config_parser( cxxopts::Options& options )
{
  options.positional_help( "GRAPH" );
  options.add_options()
    ( "o, output", "Write to this file instead of standard output",
      cxxopts::value< std::string >()->default_value( DEFAULT_OUTPUT ) )
    ( "t, type", "Output type indicating whether haplotypes or reads should be"
      " simulated (inferred from file extension if not provided). Values are: "
      "'sequence', 'fastq', or 'fasta'. The 'sequence' format simply puts each"
      " read in a line and the 'fasta' outputs haplotypes in FASTA format.",
      cxxopts::value< fmt::Type >() )
    ( "p, ploidy", "Set the ploidy",
      cxxopts::value< unsigned int >()->default_value( DEFAULT_PLOIDY ) )
    ( "l, read-length", "Read length", cxxopts::value< unsigned int >() )
    ( "n, num-reads", "Number of reads", cxxopts::value< unsigned int >() )
    ( "e, sub-rate", "Base substitution error rate",
      cxxopts::value< double >()->default_value( DEFAULT_SUBRATE ) )
    ( "i, indel-rate", "Indel error rate",
      cxxopts::value< double >()->default_value( DEFAULT_INDRATE ) )
    ( "s, random-seed", "Seed for random generator",
      cxxopts::value< unsigned int >()->default_value( DEFAULT_RNDSEED ) )
    ( "f, forward-only", "Simulate reads only from forward strand",
      cxxopts::value< bool >()->default_value( DEFAULT_FORWARD ) )
    ( "N, allow-Ns", "Allow reads to be sampled from the graph with Ns in them",
      cxxopts::value< bool >()->default_value( DEFAULT_ALLOWNS ) )
    ( "h, help", "Print this message and exit" )
    ;

  options.add_options( "positional" )
    ( "graph", "graph file (vg or xg)", cxxopts::value< std::string >() )
    ;
  options.parse_positional( { "graph" } );
}

  cxxopts::ParseResult
parse_opts( cxxopts::Options& options, int& argc, char**& argv )
{
  auto result = options.parse( argc, argv );

  if ( result.count( "help" ) ) {
    cout << options.help( { "" } ) << endl;
    throw EXIT_SUCCESS;
  }

  if ( ! result.count( "graph" ) ) {
    throw cxxopts::OptionParseException( "Graph file must be specified" );
  }
  if ( ! readable( result[ "graph" ].as< std::string >() ) ) {
    throw cxxopts::OptionParseException( "Graph file not found" );
  }

  if ( ! result.count( "type" ) &&
      result[ "output" ].as< std::string >() == DEFAULT_OUTPUT ) {
    throw cxxopts::OptionParseException( "File type must be specified" );
  }

  if ( is_reads( fmt::get_type( result ) ) ) {
    if ( ! result.count( "read-length" ) ) {
      throw cxxopts::OptionParseException( "Read length must be specified" );
    }
    if ( ! result.count( "num-reads" ) ) {
      throw cxxopts::OptionParseException( "Number of reads must be specified" );
    }
  }

  return result;
}

  int
main( int argc, char* argv[] )
{
  cxxopts::Options options( argv[0], LONG_DESC );
  config_parser( options );

  try {
    auto res = parse_opts( options, argc, argv );

    std::string graph_path = res[ "graph" ].as< std::string >();
    std::string output = res[ "output" ].as< std::string >();
    fmt::Type type = fmt::get_type( res );
    unsigned int ploidy = res[ "ploidy" ].as< unsigned int >();
    unsigned int readlen = 0;
    unsigned int numreads = 0;
    if ( is_reads( type ) ) {
      readlen = res[ "read-length" ].as< unsigned int >();
      numreads = res[ "num-reads" ].as< unsigned int >();
    }
    double subrate = res[ "sub-rate" ].as< double >();
    double indelrate = res[ "indel-rate" ].as< double >();
    unsigned int seed = res[ "random-seed" ].as< unsigned int >();
    bool forward = res[ "forward-only" ].as< bool >();
    bool allow_ns = res[ "allow-Ns" ].as< bool >();

    VarGraph vargraph;
    ifstream ifs( graph_path, ifstream::in | ifstream::binary );
    if ( ends_with( graph_path, ".vg" ) ) {
      vargraph.from_stream( ifs );
    }
    else {
      vargraph.load( ifs );
    }

    std::vector< klibpp::KSeq > seqs;
    simulate( vargraph, ploidy, seed, readlen, numreads, subrate, indelrate,
        forward, allow_ns, type, seqs );
    write_output( output, seqs, type );
  }
  catch ( const cxxopts::OptionException& e ) {
    cerr << "Error: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch ( const int& rv ) {
    return rv;
  }

  return EXIT_SUCCESS;
}

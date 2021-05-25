/**
 *    @file  alice.cpp
 *   @brief  ALignment InspeCtor and analysEr
 *
 *  A tool for analysis and inspection of resulting alignments by GraphAligner/PSI.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Mon Apr 19, 2021  17:37
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2020, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>

#include <cxxopts.hpp>
#include <gum/graph.hpp>
#include <gum/io_utils.hpp>
#include <psi/graph.hpp>
#include <psi/graph.hpp>
#include <psi/pathset.hpp>
#include <psi/utils.hpp>


using namespace psi;

constexpr const char* LONG_DESC = ( "ALICE\n"
                                    "-----\n"
                                    "ALignment InspeCtor and analysEr\n" );
// Default values for command line arguments
constexpr const char* DEFAULT_OUTPUT = "-";  // stdout

void
config_parser( cxxopts::Options& options )
{
  options.add_options( "general" )
      ( "o, output", "Write to this file instead of stdout",
        cxxopts::value< std::string >()->default_value( DEFAULT_OUTPUT ) )
      ( "g, graph", "Corresponding graph file (vg or gfa)",
        cxxopts::value< std::string >() )
      ( "h, help", "Print this message and exit" )
      ;

  options.add_options( "dstats" )
      ;

  options.add_options( "analyse" )
      ( "F, full-report", "Output full report" )
      ;

  options.add_options( "positional" )
      ( "command", "Operation type", cxxopts::value< std::string >() )
      ( "alignment", "Alignment file (GAF)", cxxopts::value< std::string >() )
      ;
  options.parse_positional( { "command", "alignment" } );
}

cxxopts::ParseResult
parse_opts( cxxopts::Options& options, int& argc, char**& argv )
{
  auto result = options.parse( argc, argv );

  if ( !result.count( "command" ) ) {  // no command specified
    options.positional_help( "COMMAND" );
    auto help_message = ( options.help( { "general" } )
                          + "\n COMMANDS:\n"
                          + "  dstats\tPrint statistics of inner-distances of read pairs\n"
                          + "  analyse\tAnalyse the given alignment file\n" );
    if ( result.count( "help" ) )
    {
      std::cout << help_message << std::endl;
      throw EXIT_SUCCESS;
    }
    std::cerr << help_message << "\n" /* extra vertical space */ << std::endl;
    throw cxxopts::OptionParseException( "No command specified" );
  }
  else if ( result[ "command" ].as< std::string >() == "dstats" ) {  // dstats
    options.custom_help( "dstats [OPTION...]" );
    options.positional_help( "ALIGNMENT" );
    if ( result.count( "help" ) ) {
      std::cout << options.help( { "general", "dstats" } ) << std::endl;
      throw EXIT_SUCCESS;
    }
  }
  else if ( result[ "command" ].as< std::string >() == "analyse" ) {  // analyse
    options.custom_help( "analyse [OPTION...]" );
    options.positional_help( "ALIGNMENT" );
    if ( result.count( "help" ) ) {
      std::cout << options.help( { "general", "analyse" } ) << std::endl;
      throw EXIT_SUCCESS;
    }
  }
  else {
    throw cxxopts::OptionParseException( "Unknown command '" +
                                         result[ "command" ].as< std::string >() + "'" );
  }

  /* Verifying general arguments */
  if ( ! result.count( "graph" ) ) {
    throw cxxopts::OptionParseException( "Graph file must be specified" );
  }
  if ( ! readable( result[ "graph" ].as< std::string >() ) ) {
    throw cxxopts::OptionParseException( "Graph file not found" );
  }

  /* Verifying positional arguments */
  if ( !result.count( "alignment" ) ) {
    throw cxxopts::OptionParseException( "Alignment file must be specified" );
  }
  if ( !readable( result[ "alignment" ].as< std::string >() ) ) {
    throw cxxopts::OptionParseException( "Alignment file not found" );
  }

  return result;
}

namespace gaf {
  struct GAFRecord {
    /* === CONSTANTS === */
    static constexpr const unsigned int GAF_MANDATORY_FIELDS_NUM = 12;

    static constexpr const unsigned int QNAME_IDX = 0;
    static constexpr const unsigned int QLEN_IDX = 1;
    static constexpr const unsigned int QSTART_IDX = 2;
    static constexpr const unsigned int QEND_IDX = 3;
    static constexpr const unsigned int QORIENT_IDX = 4;
    static constexpr const unsigned int PATH_IDX = 5;
    static constexpr const unsigned int PLEN_IDX = 6;
    static constexpr const unsigned int PSTART_IDX = 7;
    static constexpr const unsigned int PEND_IDX = 8;
    static constexpr const unsigned int MATCH_IDX = 9;
    static constexpr const unsigned int BLOCK_IDX = 10;
    static constexpr const unsigned int QUAL_IDX = 11;
    static constexpr const unsigned int AUX_START_IDX = 12;
    /* === DATA MEMBERS === */
    std::string q_name;
    std::size_t q_len;
    std::size_t q_start;
    std::size_t q_end;
    bool q_fwd;
    std::string path;
    std::size_t p_len;
    std::size_t p_start;
    std::size_t p_end;
    std::size_t match;
    std::size_t block;
    std::size_t qual;
    phmap::flat_hash_map< std::string, std::string > tag_az;
    phmap::flat_hash_map< std::string, int32_t > tag_i;
    phmap::flat_hash_map< std::string, float > tag_f;
    /* === STATIC METHODS === */
    static inline int64_t
    generate_path_id( )
    {
      static std::atomic< int64_t > path_id = 0;
      path_id.fetch_add( 1 );
      return path_id;
    }

    static inline bool
    is_orientation_char( int c )
    {
      if ( c == '>' || c == '<' ) return true;
      return false;
    }

    static inline bool
    parse_node_orientation( int c )
    {
      if ( c == '>' ) return false;
      else if ( c == '<' ) return true;
      else throw std::runtime_error( "expected '<' or '>', got '" +
                                     std::to_string( c ) + "' instead" );
    }

    template< typename TIter >
    static inline TIter
    parse_oriented_node( TIter begin, TIter end,
                         std::function< void( std::string, bool ) > callback )
    {
      thread_local static std::string segid;

      assert( is_orientation_char( *begin ) );
      bool reverse = GAFRecord::parse_node_orientation( *begin++ );
      segid.clear();
      for ( ; begin != end && !is_orientation_char( *begin ); ++begin ) segid += *begin;
      callback( segid, reverse );
      return begin;
    }
    /* === LIFECYCLE === */
    GAFRecord( )
      : q_name( "" ), q_len( 0 ), q_start( 0 ), q_end( 0 ), q_fwd( false ), path( "" ),
        p_len( 0 ), p_start( 0 ), p_end( 0 ), match( 0 ), block( 0 ), qual( 0 )
    { }

    GAFRecord( std::string line )
      : GAFRecord( )
    {
      psi::trim( line );

      const std::regex re(R"(\t)");
      std::vector< std::string > tokens(
          std::sregex_token_iterator( line.begin(), line.end(), re, -1 ),
          std::sregex_token_iterator()
        );

      if ( tokens.size() < GAF_MANDATORY_FIELDS_NUM ) {
        throw std::runtime_error( "missing mandatory field(s) in input GAF file" );
      }

      try {
        this->q_name = tokens[ QNAME_IDX ];
        this->q_len = std::stoull( tokens[ QLEN_IDX ] );
        this->q_start = ( tokens[ QSTART_IDX ] == "*" ?
                          this->q_start : std::stoull( tokens[ QSTART_IDX ] ) );
        this->q_end = ( tokens[ QEND_IDX ] == "*" ?
                        this->q_end : std::stoull( tokens[ QEND_IDX ] ) );

        if ( tokens[ QORIENT_IDX ] == "+" ) this->q_fwd = true;
        else if ( tokens[ QORIENT_IDX ] == "-" ) this->q_fwd = false;
        else if ( tokens[ QORIENT_IDX ] == "*" ) /* do nothing */;
        else throw std::invalid_argument( "invalid query orientation char" );

        this->path = ( tokens[ PATH_IDX ] == "*" ?
                       this->path : tokens[ PATH_IDX ] );
        this->p_len = ( tokens[ PLEN_IDX ] == "*" ?
                        this->p_len : std::stoull( tokens[ PLEN_IDX ] ) );
        this->p_start = ( tokens[ PSTART_IDX ] == "*" ?
                          this->p_start : std::stoull( tokens[ PSTART_IDX ] ) );
        this->p_end = ( tokens[ PEND_IDX ] == "*" ?
                        this->p_end : std::stoull( tokens[ PEND_IDX ] ) );
        this->match = ( tokens[ MATCH_IDX ] == "*" ?
                        this->match : std::stoull( tokens[ MATCH_IDX ] ) );
        this->block = ( tokens[ BLOCK_IDX ] == "*" ?
                        this->block : std::stoull( tokens[ BLOCK_IDX ] ) );
        this->qual = ( tokens[ QUAL_IDX ] == "*" ?
                       this->qual : std::stoull( tokens[ QUAL_IDX ] ) );
      }
      catch ( std::invalid_argument const& ) {
        std::cerr << "! Error in parsing input GAF:" << std::endl;
        std::cerr << "  === Record tokens ===\n"
                  << "  * QNAME: " << tokens[ QNAME_IDX ] << "\n"
                  << "  * QLEN: " << tokens[ QLEN_IDX ] << "\n"
                  << "  * QSTART: " << tokens[ QSTART_IDX ] << "\n"
                  << "  * QEND: " << tokens[ QEND_IDX ] << "\n"
                  << "  * QORIENT: " << tokens[ QORIENT_IDX ] << "\n"
                  << "  * PATH: " << tokens[ PATH_IDX ] << "\n"
                  << "  * PLEN: " << tokens[ PLEN_IDX ] << "\n"
                  << "  * PSTART: " << tokens[ PSTART_IDX ] << "\n"
                  << "  * PEND: " << tokens[ PEND_IDX ] << "\n"
                  << "  * MATCH: " << tokens[ MATCH_IDX ] << "\n"
                  << "  * BLOCK: " << tokens[ BLOCK_IDX ] << "\n"
                  << "  * QUAL: " << tokens[ QUAL_IDX ] << "\n"
                  << std::endl;
      }

      for ( std::size_t i = AUX_START_IDX; i < tokens.size(); ++i ) {
        this->parse_tag( tokens[ i ] );
      }
    }
    /* === OPERATORS === */
    inline operator bool() const
    {
      return !this->empty();
    }
    /* === METHODS === */
    inline bool
    empty( ) const
    {
      return this->q_name.empty();
    }

    inline bool
    is_valid( ) const
    {
      return this->block != 0;
    }

    template< typename TGraph >
    inline typename TGraph::dynamic_type::path_type
    parse_stable_path( TGraph const& graph ) const
    {
      typedef typename TGraph::dynamic_type::path_type path_type;

      assert( !this->empty() );
      assert( !GAFRecord::is_orientation_char( this->path[ 0 ] ) );
      path_type p( GAFRecord::generate_path_id(), this->q_name );
      _parse_stable_path( p, this->path.begin(), this->path.end(), graph );
      return p;
    }

    template< typename TGraph >
    inline typename TGraph::dynamic_type::path_type
    parse_oriented_path( TGraph const& graph ) const
    {
      typedef typename TGraph::dynamic_type::path_type path_type;

      assert( !this->empty() );
      assert( GAFRecord::is_orientation_char( this->path[ 0 ] ) );
      path_type p( GAFRecord::generate_path_id(), this->q_name );
      auto begin = this->path.begin();
      auto end = this->path.end();
      while ( begin != end ) {
        begin = GAFRecord::parse_oriented_node(
            begin,
            end,
            [&p, &graph, this]( std::string segid, bool reverse ) {
              auto id = graph.id_by_coordinate( std::stoll( segid ) );
              if ( !graph.has_node( id ) ) {
                this->_parse_stable_path( p, segid.begin(), segid.end(), graph, reverse );
              }
              p.add_node( id, reverse );
            } );
      }
      return p;
    }

    inline bool
    is_oriented_path( ) const
    {
      return GAFRecord::is_orientation_char( this->path[ 0 ] );
    }

    template< typename TGraph >
    inline typename TGraph::dynamic_type::path_type
    parse_path( TGraph const& graph ) const
    {
      assert( !this->empty() );
      if ( this->is_oriented_path() ) return parse_oriented_path( graph );
      return parse_stable_path( graph );
    }

  private:
    inline void
    parse_tag( std::string field ) {
      psi::trim( field );

      const std::regex re(R"(:)");
      std::vector< std::string > tokens(
          std::sregex_token_iterator( field.begin(), field.end(), re, -1 ),
          std::sregex_token_iterator()
        );

      if ( tokens.size() != 3 || tokens[ 0 ].size() != 2 || tokens[ 1 ].size() != 1 ) {
        std::cerr << "! Warning: ignoring tag '" << field << "' (wrong tokens)" << std::endl;
        return;
      }
      auto const& name = tokens[ 0 ];
      char type = tokens[ 1 ][ 0 ];
      auto&& value = tokens[ 2 ];

      try {
        switch( type ) {
        case 'i':
          this->tag_i[ name ] = ( value != "*" ? std::stoul( value ) : 0 );
          break;
        case 'f':
          this->tag_f[ name ] = ( value != "*" ? std::stof( value ) : 0.0 );
          break;
        case 'A':
          if ( value.size() != 1 ) throw std::length_error( "expected one character" );
        case 'Z':
          this->tag_az[ name ] = ( value != "*" ? std::move( value ) : "" );
          break;
        }
      }
      catch ( ... ) {
        std::cerr << "! Error in parsing tag value:" << std::endl;
        std::cerr << "  === Tag tokens ===\n"
                  << "  * NAME: " << tokens[ 0 ] << "\n"
                  << "  * TYPE: " << tokens[ 1 ] << "\n"
                  << "  * VALUE: " << tokens[ 2 ] << "\n"
                  << std::endl;
      }
    }

    template< typename TIter, typename TGraph >
    inline void
    _parse_stable_path( typename TGraph::dynamic_type::path_type& p,
                        TIter begin, TIter end, TGraph const& graph, bool reverse=false ) const
    {
      throw std::runtime_error( "parsing path with stable ID is not implemented" );
    }
  };

  inline std::ostream&
  operator<<( std::ostream& os, GAFRecord const& record )
  {
    os << "=== Record ===\n"
       << "* Query name: " << record.q_name << "\n"
       << "* Query length: " << record.q_len << "\n"
       << "* Query start: " << record.q_start << "\n"
       << "* Query end: " << record.q_end << "\n"
       << "* Query strand: " << ( record.q_fwd ? "Forward" : "Reverse" ) << "\n"
       << "* Path: " << record.path << "\n"
       << "* Path length: " << record.p_len << "\n"
       << "* Path start: " << record.p_start << "\n"
       << "* Path end: " << record.p_end << "\n"
       << "* No. of matches: " << record.match << "\n"
       << "* Alignment block length: " << record.block << "\n"
       << "* Mapping quality: " << record.qual;
    for ( auto const& elem : record.tag_az ) {
      os << "\n* " << elem.first << " (A/Z): " << elem.second;
    }
    for ( auto const& elem : record.tag_i ) {
      os << "\n* " << elem.first << " (i): " << elem.second;
    }
    for ( auto const& elem : record.tag_f ) {
      os << "\n* " << elem.first << " (f): " << elem.second;
    }
    return os;
  }

  inline GAFRecord
  next( std::istream& is )
  {
    std::string line;
    std::getline( is, line );
    if ( is ) return GAFRecord( std::move( line ) );
    else return GAFRecord();
  }
}

template< typename TPathSet, typename TGraph >
void
index_reference_paths( TPathSet& pathset, TGraph& graph )
{
  typedef typename TGraph::id_type id_type;

  std::vector< id_type > nodes;
  graph.for_each_path(
      [&nodes, &graph, &pathset]( auto rank, auto pid ) {
        typedef typename TPathSet::value_type path_type;
        for ( auto n : graph.path( pid ) ) nodes.push_back( n );
        typename path_type::nodes_type cnodes( nodes );
        pathset.push_back( path_type( &graph, std::move( cnodes ) ) );
        nodes.clear();
        return true;
      } );
  pathset.initialize();
}

template< typename TPathSet, typename TPath >
std::pair< std::size_t, std::size_t >
ref_pos( TPathSet& rpaths, typename TPathSet::graph_type graph, TPath& path )
{
  typedef typename TPathSet::graph_type graph_type;

  psi::Path< graph_type, psi::Dynamic > one_node_path( &graph );
  std::size_t offset = 0;
  for ( auto n : path ) {
    auto id = path.id_of( n );
    bool reverse = path.is_reverse( n );
    one_node_path.push_back( id );
    auto occs = rpaths.get_occurrences( one_node_path );
    one_node_path.pop_front();
    if ( occs.size() > 1 ) throw std::runtime_error( "node occurs on multiple paths" );
    if ( occs.size() == 0 || reverse ) {
      offset += ( 1 - 2*( int )( !reverse ) ) * graph.node_length( id );
    }
    if ( occs.size() == 0 ) continue;
    std::size_t idx = occs[0].first;
    std::size_t pos = rpaths[ idx ].select( occs[0].second ) + offset;
    return { idx, pos };
  }
  throw std::runtime_error( "position not found" );
}

template< typename TPathSet >
long long int
distance_estimate( TPathSet& rpaths, typename TPathSet::graph_type& graph,
                   gaf::GAFRecord& rec1, gaf::GAFRecord& rec2)
{
  std::size_t idx1;
  std::size_t pos1;
  std::size_t idx2;
  std::size_t pos2;
  auto path1 = rec1.parse_path( graph );
  auto path2 = rec2.parse_path( graph );
  std::tie( idx1, pos1 ) = ref_pos( rpaths, graph, path1 );
  std::tie( idx2, pos2 ) = ref_pos( rpaths, graph, path2 );
  if ( idx1 != idx2 ) throw std::runtime_error( "not in the same reference path" );

  long long int distance;
  if ( pos1 < pos2 ) {
    pos1 += rec1.p_end;
    pos2 -= rec2.p_end;
    distance = pos2 - pos1;
  }
  else if ( pos2 < pos1 ) {
    pos2 += rec2.p_end;
    pos1 -= rec1.p_end;
    distance = pos1 - pos2;
  }
  else throw std::runtime_error( "same reference positions" );

  if ( distance < 0 ) std::runtime_error( "negative distance" );

  return  distance;
}

void
dstats( cxxopts::ParseResult& res )
{
  typedef gum::SeqGraph< gum::Succinct > graph_type;

  // Fetching input parameters
  std::string output = res[ "output" ].as< std::string >();
  std::string graph_path = res[ "graph" ].as< std::string >();
  std::string aln_path = res[ "alignment" ].as< std::string >();

  // Opening output file for writing
  std::ostream ost( nullptr );
  std::ofstream ofs;
  if ( output == "-" ) ost.rdbuf( std::cout.rdbuf() );
  else{
    ofs.open( output, std::ofstream::out | std::ofstream::binary );
    ost.rdbuf( ofs.rdbuf() );
  }
  if ( !ost ) throw std::runtime_error( "output file cannot be opened" );

  // Loading input graph
  graph_type graph;
  std::cout << "Loading input graph..." << std::endl;
  gum::util::load( graph, graph_path, true );

  // Loading reference paths
  psi::PathSet< psi::Path< graph_type, psi::Compact > > rpaths( graph );
  std::cout << "Loading reference paths..." << std::endl;
  index_reference_paths( rpaths, graph );

  // Opening alignment file for reading
  std::ifstream ifs( aln_path, std::ifstream::in | std::ifstream::binary );
  std::cout << "Estimating inner-distances between aligned read pairs..." << std::endl;
  gaf::GAFRecord record1 = gaf::next( ifs );
  gaf::GAFRecord record2 = gaf::next( ifs );
  while ( record1 && record2 ) {
    auto distance = distance_estimate( rpaths, graph, record1, record2 );
    ost << distance << std::endl;
    record1 = gaf::next( ifs );
    record2 = gaf::next( ifs );
  }
}

void
analyse( cxxopts::ParseResult& res )
{
  typedef gum::SeqGraph< gum::Succinct > graph_type;

  // Fetching input parameters
  std::string output = res[ "output" ].as< std::string >();
  std::string graph_path = res[ "graph" ].as< std::string >();
  std::string aln_path = res[ "alignment" ].as< std::string >();
  bool full = res.count( "full-report" );

  // Opening output file for writing
  std::ostream ost( nullptr );
  std::ofstream ofs;
  if ( output == "-" ) ost.rdbuf( std::cout.rdbuf() );
  else {
    ofs.open( output, std::ofstream::out | std::ofstream::binary );
    ost.rdbuf( ofs.rdbuf() );
  }
  if ( !ost ) throw std::runtime_error( "output file cannot be opened" );

  // Loading input graph
  graph_type graph;
  std::cout << "Loading input graph..." << std::endl;
  gum::util::load( graph, graph_path, true );

  // Opening alignment file for reading
  std::ifstream ifs( aln_path, std::ifstream::in | std::ifstream::binary );
  std::cout << "Analysing..." << std::endl;

  std::string del = "\t";
  struct Tuple {
    std::size_t paired;
    std::size_t single;
    std::size_t invalid;
  };

  phmap::flat_hash_map< std::string, Tuple > counts;
  std::string name;
  gaf::GAFRecord record = gaf::next( ifs );
  std::size_t valids = 0;
  std::size_t invalids = 0;
  std::size_t nrecords = 0;
  std::size_t with_valid = 0;
  std::size_t paired = 0;
  std::size_t uniq_paired = 0;
  std::size_t uniq_single = 0;
  std::size_t multis = 0;
  std::size_t with_invalid = 0;

  while( record ) {
    ++nrecords;
    name = record.q_name;
    if ( !record.is_valid() ) {  // invalid
      ++counts[ name ].invalid;
      ++invalids;
    }
    else {
      ++valids;
      auto fnptr = record.tag_az.find( "fn" );
      if ( fnptr != record.tag_az.end() ) {  // paired
        record = gaf::next( ifs );
        ++nrecords;
        if ( !record.is_valid() ) {
          ++counts[ name ].single;
          ++counts[ name ].invalid;
          ++invalids;
        }
        else if ( record.q_name != name ) {
          ++counts[ name ].single;
          std::cerr << "! Warning: missing next fragment alignment of '" << name << "'"
                    << std::endl;
          continue;
        }
        else {
          ++valids;
          auto fpptr = record.tag_az.find( "fp" );
          if ( fpptr == record.tag_az.end() || fpptr->second != fnptr->second ) {
            std::cerr << "! Warning: missing proper 'fp' tag in next fragment alignment of '"
                      << name << "'" << std::endl;
            std::cerr << "  consider two alignments as unpaired" << std::endl;
            counts[ name ].single += 2;
          }
          else {
            ++counts[ name ].paired;
          }
        }
      }
      else {  // single
        ++counts[ name ].single;
      }
    }
    record = gaf::next( ifs );
  }

  if ( full ) {
    ost << "#RNAME: read name\n"
        << "#NP: number of paired alignments\n"
        << "#NS: number of single alignments\n"
        << "#NI: number of invalid alignments\n"
        << "RNAME" + del + "NP" + del + "NS" + del + "NI" << std::endl;
    for ( auto const& elem : counts ) {
      ost << elem.first << del << elem.second.paired << del << elem.second.single << del
          << elem.second.invalid << std::endl;
    }
  }
  else {
    for ( auto const& elem : counts ) {
      auto const& cnts = elem.second;
      if ( cnts.paired != 0 ) ++paired;
      if ( cnts.paired != 0 || cnts.single != 0 ) ++with_valid;
      if ( cnts.paired == 1 && cnts.single == 0 ) ++uniq_paired;
      else if ( cnts.paired == 0 && cnts.single == 1 ) ++uniq_single;
      else if ( cnts.paired >= 1 || cnts.single >= 1 ) ++multis;
      if ( cnts.invalid != 0 ) ++with_invalid;
    }

    ost << "#NREC: number of records\n"
        << "#NVAL: number of valid alignments\n"
        << "#NIVR: number of invalid alignments\n"
        << "#NALR: number of reads with at least one alignment\n"
        << "#NAPR: number of reads with at least one paired alignment\n"
        << "#NUQP: number of reads with a unique paired alignment\n"
        << "#NUQS: number of reads with a unique single alignment\n"
        << "#NMLT: number of reads with multiple alignments\n"
        << "#NRIN: number of reads with at least one invalid alignment\n"
        << "NREC" + del + "NVAL" + del + "NIVR" + del + "NALR" + del + "NAPR" + del
        << "NUQP" + del + "NUQS" + del + "NMLT" + del + "NRIN\n"
        << nrecords << ( valids + invalids != nrecords ? "*" : "") << del
        << valids << del << invalids << del << with_valid << del << paired << del
        << uniq_paired << del << uniq_single << del << multis << del << with_invalid
        << std::endl;
  }
}

int
main( int argc, char* argv[] )
{
  cxxopts::Options options( argv[0], LONG_DESC );
  config_parser( options );

  try {
    auto res = parse_opts( options, argc, argv );
    std::string command = res[ "command" ].as< std::string >();

    if ( command == "dstats" ) {
      dstats( res );
    }
    if ( command == "analyse" ) {
      analyse( res );
    }
    else {
      // should not reach here!
      assert( false );
    }
  }
  catch ( const cxxopts::OptionException& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch ( const int& rv ) {
    return rv;
  }

  return EXIT_SUCCESS;
}

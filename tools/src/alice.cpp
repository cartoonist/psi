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
constexpr const char* DEFAULT_ID_THRESHOLD = "0.9";

void
config_parser( cxxopts::Options& options )
{
  options.add_options( "general" )
      ( "o, output", "Write to this file instead of stdout",
        cxxopts::value< std::string >()->default_value( DEFAULT_OUTPUT ) )
      ( "g, graph", "Corresponding graph file (vg or gfa)",
        cxxopts::value< std::string >() )
      ( "P, progress", "Show progress" )
      ( "h, help", "Print this message and exit" )
      ;

  options.add_options( "dstats" )
      ;

  options.add_options( "analyse" )
      ( "F, full-report", "Output full report" )
      ( "I, identity-threshold", "Minimum identity score of a good alignment",
        cxxopts::value< float >()->default_value( DEFAULT_ID_THRESHOLD ) )
      ( "T, ground-truth", "Ground truth alignment (GAF)",
        cxxopts::value< std::string >()->default_value( "" ) )
      ( "m, trim-name",
        "Trim fragment numbers at the end of read names in the input ground truth set" )
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

    if ( result.count( "ground-truth" ) &&
         ! readable( result[ "ground-truth" ].as< std::string >() ) ) {
      throw cxxopts::OptionParseException( "Ground truth alignment file not found" );
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

      if ( field.size() < 6 || field[ 0 ] == ':' || field[ 1 ] == ':' ||
           field[ 2 ] != ':' || field[ 3 ] == ':' || field[4] != ':' ) {
        std::cerr << "! Warning: ignoring tag '" << field << "' (wrong tokens)" << std::endl;
        return;
      }

      std::string name = field.substr( 0, 2 );
      char type = field[ 3 ];
      std::string value = field.substr( 5 );

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
                  << "  * NAME: " << name << "\n"
                  << "  * TYPE: " << type << "\n"
                  << "  * VALUE: " << value << "\n"
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

  inline float
  get_identity( GAFRecord const& record )
  {
    auto id_itr = record.tag_f.find( "id" );
    if ( id_itr != record.tag_f.end() ) return id_itr->second;
    auto dv_itr = record.tag_f.find( "dv" );
    if ( dv_itr != record.tag_f.end() ) return 1 - dv_itr->second;
    return 0;
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
        std::cerr << "! INFO Fetching reference path " << rank << "..." << std::endl;
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
  long long int offset = 0;
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
  std::cerr << "Loading input graph..." << std::endl;
  gum::util::load( graph, graph_path, true );

  // Loading reference paths
  psi::PathSet< psi::Path< graph_type, psi::Compact > > rpaths( graph );
  std::cerr << "Loading reference paths..." << std::endl;
  index_reference_paths( rpaths, graph );

  // Opening alignment file for reading
  std::ifstream ifs( aln_path, std::ifstream::in | std::ifstream::binary );
  std::cerr << "Estimating inner-distances between aligned read pairs..." << std::endl;
  gaf::GAFRecord record1 = gaf::next( ifs );
  gaf::GAFRecord record2 = gaf::next( ifs );
  while ( record1 && record2 ) {
    if ( record1.is_valid() && record2.is_valid() ) {
      auto distance = distance_estimate( rpaths, graph, record1, record2 );
      ost << distance << std::endl;
    }
    record1 = gaf::next( ifs );
    record2 = gaf::next( ifs );
  }
}

template< typename TGraph >
auto
load_ground_truth( std::string const& truth_path, TGraph const& graph,
                   bool trim_name=false )
{
  struct OrientedPos {
    typename TGraph::dynamic_type::path_type::value_type oriented_id;
    std::size_t offset;
    OrientedPos( typename TGraph::dynamic_type::path_type::value_type i=0,
                 std::size_t o=0 )
      : oriented_id( i ), offset( o )
    { }
  };
  phmap::flat_hash_map< std::string, std::vector< OrientedPos > > truth;

  if ( truth_path.empty() ) return truth;

  // Opening ground truth alignment file for reading
  std::ifstream ifs( truth_path, std::ifstream::in | std::ifstream::binary );
  std::cerr << "Loading ground truth alignments..." << std::endl;

  gaf::GAFRecord record = gaf::next( ifs );
  while ( record ) {
    auto name = record.q_name;
    if ( trim_name ) name.resize( name.size() - 2 );
    auto&& loci = truth[ name ];
    if ( record.is_valid() ) {
      auto path = record.parse_path( graph );
      loci.push_back( { path.front(), record.p_start } );
      if ( loci.size() > 2 ) {
        std::cerr << "! Warning: '" << name << "' has more than two alignments"
                  << std::endl;
      }
    }
    record = gaf::next( ifs );
  }

  std::size_t multiple = 0;
  std::array< std::size_t, 3 > counters = { 0, 0, 0 };
  for ( auto const& element : truth ) {
    if ( element.second.size() > 2 ) ++multiple;
    else ++counters[ element.second.size() ];
  }

  std::cerr << "Loaded ground truth alignments with (0, 1, 2, 3+) fragments: ("
            << counters[ 0 ] << ", " << counters[ 1 ] << ", " << counters[ 2 ] << ", "
            << multiple << ")" << std::endl;

  return truth;
}

template< typename TRecord, typename TMap, typename TGraph >
unsigned char
get_truth_flag( TRecord const& record, TMap const& truth, TGraph const& graph,
                bool progress=false )
{
  unsigned char flag = 0;
  auto&& name = record.q_name;
  auto found = truth.find( name );
  if ( found == truth.end() ) {
    if ( progress ) std::cerr << std::endl;
    std::cerr << "! Warning: '" << name << "' has no ground truth alignment" << std::endl;
    return 0;
  }
  auto path = record.parse_path( graph );
  if ( found->second.size() > 4 ) throw std::runtime_error( "too many ground truth fragments" );
  for ( auto it = found->second.begin(); it != found->second.end(); ++it ) {
    if ( it->oriented_id == path.front() ) {
      ++flag;
      flag <<= 1;
      if ( it->offset == record.p_start )
      {
        ++flag;
      }
      flag <<= 2*( it - found->second.begin() );
      break;
    }
  }
  return flag;
}

void
analyse( cxxopts::ParseResult& res )
{
  typedef gum::SeqGraph< gum::Succinct > graph_type;

  // Fetching input parameters
  std::string output = res[ "output" ].as< std::string >();
  std::string graph_path = res[ "graph" ].as< std::string >();
  std::string truth_path = res[ "ground-truth" ].as< std::string >();
  std::string aln_path = res[ "alignment" ].as< std::string >();
  float id_threshold = res[ "identity-threshold" ].as< float >();
  bool full = res.count( "full-report" );
  bool trim = res.count( "trim-name" );
  bool progress = res.count( "progress" );

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
  std::cerr << "Loading input graph..." << std::endl;
  gum::util::load( graph, graph_path, true );

  // Loading ground truth if available
  auto truth = load_ground_truth( truth_path, graph, trim );

  // Opening alignment file for reading
  std::ifstream ifs( aln_path, std::ifstream::in | std::ifstream::binary );
  std::cerr << "Analysing..." << std::endl;

  std::string del = "\t";
  struct Tuple {
    std::size_t paired;
    std::size_t single;
    std::size_t hi_paired;
    std::size_t hi_single;
    std::size_t invalid;
    unsigned char truth_flag;
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
  std::size_t hi_paired = 0;
  std::size_t hi_single = 0;
  std::size_t hi_uniq_paired = 0;
  std::size_t hi_uniq_single = 0;
  std::size_t multis = 0;
  std::size_t with_invalid = 0;
  std::size_t nffm = 0;
  std::size_t nppm = 0;
  std::size_t nfpm = 0;
  std::size_t nfnm = 0;
  std::size_t npnm = 0;
  std::size_t nnnm = 0;

  while( record ) {
    ++nrecords;
    if ( progress ) std::cerr << "\rAnalysing record " << nrecords << "...";
    name = record.q_name;
    auto&& cnt = counts[ name ];
    if ( !record.is_valid() ) {  // invalid
      ++cnt.invalid;
      ++invalids;
    }
    else {
      ++valids;
      cnt.truth_flag |= get_truth_flag( record, truth, graph, progress );
      auto identity = gaf::get_identity( record );
      auto tagptr = record.tag_az.find( "fn" );
      if ( tagptr != record.tag_az.end() ) {  // paired
        auto fn = tagptr->second;
        record = gaf::next( ifs );
        if ( record.q_name != name ) {
          ++cnt.single;
          if ( identity >= id_threshold ) ++cnt.hi_single;
          if ( progress ) std::cerr << std::endl;
          std::cerr << "! Warning: missing next fragment alignment of '" << name << "'"
                    << std::endl;
          continue;
        }
        else if ( !record.is_valid() ) {
          ++nrecords;  // process fetched record
          ++cnt.single;
          ++cnt.invalid;
          ++invalids;
          if ( identity >= id_threshold ) ++cnt.hi_single;
        }
        else {
          ++nrecords;  // process fetched record
          ++valids;
          auto identity2 = gaf::get_identity( record );
          tagptr = record.tag_az.find( "fp" );
          if ( tagptr == record.tag_az.end() || tagptr->second != fn ) {
            if ( progress ) std::cerr << std::endl;
            std::cerr << "! Warning: missing proper 'fp' tag in next fragment alignment of '"
                      << name << "'" << std::endl;
            std::cerr << "  consider two alignments as unpaired" << std::endl;
            cnt.single += 2;
            if ( identity >= id_threshold ) ++cnt.hi_single;
            if ( identity2 >= id_threshold ) ++cnt.hi_single;
            cnt.truth_flag |= get_truth_flag( record, truth, graph, progress );  // pick the best map
          }
          else {
            ++cnt.paired;
            if ( identity >= id_threshold && identity2 >= id_threshold ) ++cnt.hi_paired;
            else if ( identity >= id_threshold || identity2 >= id_threshold ) ++cnt.hi_single;
            cnt.truth_flag |= get_truth_flag( record, truth, graph, progress );
          }
        }
      }
      else {  // single
        ++cnt.single;
        if ( identity >= id_threshold ) ++cnt.hi_single;
      }
    }
    record = gaf::next( ifs );
  }
  if ( progress ) std::cerr << "\rAnalysing " << nrecords << " record... Done.";

  if ( full ) {
    ost << "#RNAME: read name\n"
        << "#NP: number of paired alignments\n"
        << "#NS: number of single alignments\n"
        << "#NI: number of invalid alignments\n"
        << "#NHP: number of paired alignments with high identity score\n"
        << "#NHS: number of single alignments with high identity score\n"
        << "#MTF: alignment truth flag\n"
        << "RNAME" + del + "NP" + del + "NS" + del + "NI" + del + "NHP" + del + "NHS" + del + "MTF"
        << std::endl;
    for ( auto const& elem : counts ) {
      ost << elem.first << del << elem.second.paired << del << elem.second.single << del
          << elem.second.invalid << del << elem.second.hi_paired << del
          << elem.second.hi_single << del << std::bitset< 8 >( elem.second.truth_flag )
          << std::endl;
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
      if ( cnts.hi_paired != 0 ) ++hi_paired;
      if ( cnts.hi_paired == 1 ) ++hi_uniq_paired;
      if ( cnts.hi_single != 0 ) ++hi_single;
      if ( cnts.hi_single == 1 ) ++hi_uniq_single;
      if ( cnts.invalid != 0 ) ++with_invalid;
      if ( cnts.truth_flag == 0b1111 ) ++nffm;
      else if ( cnts.truth_flag == 0b1010 ) ++nppm;
      else if ( cnts.truth_flag == 0b1110 || cnts.truth_flag == 0b1011 ) ++nfpm;
      else if ( cnts.truth_flag == 0b1100 || cnts.truth_flag == 0b0011 ) ++nfnm;
      else if ( cnts.truth_flag == 0b1000 || cnts.truth_flag == 0b0010 ) ++npnm;
      else if ( cnts.truth_flag == 0b0000 ) ++nnnm;
    }

    ost << "#NREC: number of records\n"
        << "#NVAL: number of valid alignments\n"
        << "#NIVR: number of invalid alignments\n"
        << "#NALR: number of reads with at least one alignment\n"
        << "#NAPR: number of reads with at least one paired alignment\n"
        << "#NUQP: number of reads with a unique paired alignment\n"
        << "#NUQS: number of reads with a unique single alignment\n"
        << "#NHIP: number of reads with a at least one paired alignment with high identity score\n"
        << "#NHIS: number of reads with a at least one single alignment with high identity score\n"
        << "#NUHP: number of reads with a unique paired alignment with high identity score\n"
        << "#NUHS: number of reads with a unique single alignment with high identity score\n"
        << "#NMLT: number of reads with multiple alignments\n"
        << "#NRIN: number of reads with at least one invalid alignment\n"
        << "#NFFM: number of reads with fully-true alignments for both ends\n"
        << "#NPPM: number of reads with partially-true alignments for both ends\n"
        << "#NFPM: number of reads with fully-true alignments for one end and partial ones for the other\n"
        << "#NFNM: number of reads with fully-true alignments for one end and no true alignment for the other\n"
        << "#NPNM: number of reads with partially-true alignments for one end and no true alignment for the other\n"
        << "#NNNM: number of reads with no true alignments for both ends\n"
        << "NREC" + del + "NVAL" + del + "NIVR" + del + "NALR" + del + "NAPR" + del
        << "NUQP" + del + "NUQS" + del + "NHIP" + del + "NHIS" + del + "NUHP" + del
        << "NUHS" + del + "NMLT" + del + "NRIN" + del + "NFFM" + del + "NPPM" + del
        << "NFPM" + del + "NFNM" + del + "NPNM" + del + "NNNM" + "\n"
        << nrecords << ( valids + invalids != nrecords ? "*" : "") << del
        << valids << del << invalids << del << with_valid << del << paired << del
        << uniq_paired << del << uniq_single << del << hi_paired << del << hi_single
        << del << hi_uniq_paired << del << hi_uniq_single << del << multis << del
        << with_invalid << del << nffm << del << nppm << del << nfpm << del << nfnm
        << del << npnm << del << nnnm << std::endl;
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
    else if ( command == "analyse" ) {
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

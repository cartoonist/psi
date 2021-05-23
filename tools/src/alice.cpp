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
      catch ( std::invalid_argument const& e ) {
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

    template< typename TGraph >
    inline typename TGraph::dynamic_type::path_type
    parse_stable_path( TGraph const& graph ) const
    {
      typedef typename TGraph::dynamic_type::path_type path_type;

      assert( !this->empty() );
      assert( !GAFRecord::is_orientation_char( this->path[0] ) );
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
      assert( GAFRecord::is_orientation_char( this->path[0] ) );
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
      return GAFRecord::is_orientation_char( this->path[0] );
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
    template< typename TIter, typename TGraph >
    inline void
    _parse_stable_path( typename TGraph::dynamic_type::path_type& p,
                        TIter begin, TIter end, TGraph const& graph, bool reverse=false ) const
    {
      throw std::runtime_error( "parsing path with stable ID is not implemented" );
    }
  };

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
  std::ofstream ofs( output, std::ofstream::out | std::ofstream::binary );
  if ( !ofs ) throw std::runtime_error( "output file cannot be opened" );

  // Loading input graph
  graph_type graph;
  std::cout << "Loading input graph..." << std::endl;
  gum::util::load( graph, graph_path, true );

  // Loading reference paths
  psi::PathSet< psi::Path< graph_type, psi::Compact > > rpaths( graph );
  index_reference_paths( rpaths, graph );

  // Opening alignment file for reading
  std::ifstream ifs( aln_path, std::ifstream::in | std::ifstream::binary );
  gaf::GAFRecord record1 = gaf::next( ifs );
  gaf::GAFRecord record2 = gaf::next( ifs );
  while ( record1 && record2 ) {
    auto distance = distance_estimate( rpaths, graph, record1, record2 );
    ofs << distance << std::endl;
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

  // Opening output file for writing
  std::ofstream ofs( output, std::ofstream::out | std::ofstream::binary );
  if ( !ofs ) throw std::runtime_error( "output file cannot be opened" );

  // Loading input graph
  graph_type graph;
  std::cout << "Loading input graph..." << std::endl;
  gum::util::load( graph, graph_path, true );

  // Opening alignment file for reading
  std::ifstream ifs( aln_path, std::ifstream::in | std::ifstream::binary );
  phmap::flat_hash_map< std::string, std::size_t > counts;
  std::string name;
  gaf::GAFRecord record = gaf::next( ifs );
  std::size_t total_alns = 0;
  std::size_t total_multi_alns = 0;
  while( record ) {
    name = record.q_name;
    name.pop_back();  // remove 1 or 2 at the end
    name.pop_back();  // remove the separator
    auto found = counts.find( name );
    if ( found == counts.end() ) counts[ name ] = 0;
    found->second += 1;
    ++total_alns;
    record = gaf::next( ifs );
  }
  for ( auto entry : counts ) if ( entry.second > 2 ) ++total_multi_alns;
  std::cout << "Total number of alignments: " << total_alns << std::endl;
  std::cout << "Total pairs with valid alignments: " << counts.size() << std::endl;
  std::cout << "Total pairs with multiple alignments: " << total_multi_alns << std::endl;
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

/**
 *    @file  dindexctl.cpp
 *   @brief  Distance index hacking tool
 *
 *  A tool for hacking (e.g. compressing, merging) distance indices.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Thu Dec 17, 2020  01:22
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2020, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include <cxxopts.hpp>
#include <gum/graph.hpp>
#include <gum/io_utils.hpp>
#include <psi/graph.hpp>
#include <psi/crs_matrix.hpp>


using namespace psi;

constexpr const char* LONG_DESC = ( "dindexctl\n"
                                    "---------\n"
                                    "Hacking tool for distance indices\n" );
// Default values for command line arguments
constexpr const char* DEFAULT_OUTPUT = "-";  // stdout

void
config_parser( cxxopts::Options& options )
{
  options.add_options( "general" )
      ( "h, help", "Print this message and exit" )
      ;

  options.add_options( "compress" )
      ( "o, output", "Write to this file instead of standard output",
        cxxopts::value< std::string >()->default_value( DEFAULT_OUTPUT ) )
      ( "d, min-insert-size", "Distance index minimum read insert size",
        cxxopts::value< unsigned int >() )
      ( "D, max-insert-size", "Distance index maximum read insert size",
        cxxopts::value< unsigned int >() )
      ( "g, graph", "Corresponding graph file (vg or gfa)",
        cxxopts::value< std::string >() )
      ;

  options.add_options( "positional" )
      ( "command", "Operation type", cxxopts::value< std::string >() )
      ( "prefix", "Path index prefix", cxxopts::value< std::string >() )
      ;
  options.parse_positional( { "command", "prefix" } );
}

cxxopts::ParseResult
parse_opts( cxxopts::Options& options, int& argc, char**& argv )
{
  auto result = options.parse( argc, argv );

  if ( !result.count( "command" ) ) {  // no command specified
    options.positional_help( "COMMAND" );
    auto help_message = ( options.help( { "general" } )
                          + "\n COMMANDS:\n"
                          + "  compress\tCompress a distance index\n"
                          + "  merge\t\tMerge distance indices" );
    if ( result.count( "help" ) )
    {
      std::cout << help_message << std::endl;
      throw EXIT_SUCCESS;
    }
    std::cerr << help_message << "\n" /* extra vertical space */ << std::endl;
    throw cxxopts::OptionParseException( "No command specified" );
  }
  else if ( result[ "command" ].as< std::string >() == "compress" ) {  // compress
    options.custom_help( "compress [OPTION...]" );
    options.positional_help( "DINDEX" );
    if ( result.count( "help" ) ) {
      std::cout << options.help( { "general", "compress" } ) << std::endl;
      throw EXIT_SUCCESS;
    }

    if ( !result.count( "prefix" ) ) {
      throw cxxopts::OptionParseException( "Index prefix must be specified" );
    }
    if ( !readable( result[ "prefix" ].as< std::string >() ) ) {
      throw cxxopts::OptionParseException( "Index file not found" );
    }

    if ( ! result.count( "graph" ) ) {
      throw cxxopts::OptionParseException( "Graph file must be specified" );
    }
    if ( ! readable( result[ "graph" ].as< std::string >() ) ) {
      throw cxxopts::OptionParseException( "Graph file not found" );
    }

    if ( !result.count( "min-insert-size" ) ) {
      throw cxxopts::OptionParseException( "Minimum insert size must be specified" );
    }

    if ( !result.count( "max-insert-size" ) ) {
      throw cxxopts::OptionParseException( "Maximum insert size must be specified" );
    }
  }
  else if ( result[ "command" ].as< std::string >() == "merge" ) {
    throw std::runtime_error( "Operation merge has not been implemented" );
  }
  else {
    throw std::runtime_error( "Command not found" );
  }

  return result;
}

template< typename TCRSMatrix, typename TGraph >
bool
verify_distance_matrix( TCRSMatrix const& cdi, TCRSMatrix const& udi, TGraph const& g )
{
  typedef typename TCRSMatrix::ordinal_type ordinal_type;
  typedef typename TGraph::rank_type rank_type;

  rank_type cnode_rank = 0;  // current node rank
  ordinal_type start = 0;    // row start index
  ordinal_type cstart = 0;    // row start index
  ordinal_type end;          // row end index
  ordinal_type nloc = 0;     // next node loci index
  for ( ordinal_type nrow = 0; nrow < udi.numRows(); ++nrow ) {
    if ( nrow == nloc) {
      ++cnode_rank;
      if ( cnode_rank == g.get_node_count() ) nloc = udi.numRows();
      else nloc = gum::util::id_to_charorder( g, g.rank_to_id( cnode_rank + 1 ) );
    }
    assert( nrow < nloc );
    end = udi.rowMap( nrow + 1 );
    for ( ; start < end; ++start ) {
      if ( nrow <= udi.entry( start ) && udi.entry( start ) < nloc ) continue;
      else if ( udi.entry( start ) == cdi.entry( cstart ) ) ++cstart;
      else return false;
    }
  }
  return true;
}

void
compress( cxxopts::ParseResult& res )
{
  typedef CRSMatrix< crs_matrix::Compressed, bool, uint32_t, uint64_t > crsmat_type;
  typedef CRSMatrix< crs_matrix::Buffered, bool, uint32_t, uint64_t > crsmat_buffer_type;

  std::string graph_path = res[ "graph" ].as< std::string >();
  std::string pindex_prefix = res[ "prefix" ].as< std::string >();
  std::string output = res[ "output" ].as< std::string >();
  unsigned int min_size = res[ "min-insert-size" ].as< unsigned int >();
  unsigned int max_size = res[ "max-insert-size" ].as< unsigned int >();
  gum::SeqGraph< gum::Succinct > graph;
  crsmat_type dindex;
  auto index_path = pindex_prefix + "_dist_mat_" + "m" + std::to_string( min_size ) +
      "M" + std::to_string( max_size );

  std::cout << "Loading input graph..." << std::endl;
  gum::util::load( graph, graph_path );

  std::cout << "Loading distance index..." << std::endl;
  std::ifstream ifs( index_path, std::ifstream::in | std::ifstream::binary );
  if ( !ifs ) throw std::runtime_error( "distance matrix cannot be opened" );
  dindex.load( ifs );

  std::cout << "Compressing distance index..." << std::endl;
  auto cindex = util::compress_distance_index< crsmat_buffer_type >( dindex, graph );

  std::cout << "Verifying compressed distance index..." << std::endl;
  if ( !verify_distance_matrix( cindex, dindex, graph ) ) {
    std::cerr << "Verification failed!" << std::endl;
    throw EXIT_FAILURE;
  }

  std::cout << "Serialising compressed distance index..." << std::endl;
  std::ofstream ofs( output, std::ofstream::out | std::ofstream::binary );
  if ( !ofs ) throw std::runtime_error( "output file cannot be opened" );
  cindex.serialize( ofs );
}

void
merge( cxxopts::ParseResult& res )
{
}

int
main( int argc, char* argv[] )
{
  cxxopts::Options options( argv[0], LONG_DESC );
  config_parser( options );

  try {
    auto res = parse_opts( options, argc, argv );
    std::string command = res[ "command" ].as< std::string >();
    if ( command == "compress" ) {
      compress( res );
    }
    else if ( command == "merge" ) {
      merge( res );
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

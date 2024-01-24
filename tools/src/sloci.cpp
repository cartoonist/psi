/**
 *    @file  sloci.cpp
 *   @brief  Starting loci inspection tool
 *
 *  A tool for inspecting the index files storing starting loci (i.e. uncovered loci).
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Mar 14, 2023  16:23
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2023, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cxxopts.hpp>
#include <gum/io_utils.hpp>
#include <psi/graph.hpp>
#include <psi/utils.hpp>
#include <psi/seed_finder.hpp>

#include "vg/vg.pb.h"
#include "vg/stream.hpp"

using namespace psi;

constexpr const char* LONG_DESC = ( "Sloci\n"
                                    "-----\n"
                                    "Starting loci inspection tool\n" );
// Default values for command line arguments
constexpr const char* DEFAULT_OUTPUT = "-";  // stdout
constexpr const char* DEFAULT_STEP_SIZE = "1";

void
config_parser( cxxopts::Options& options )
{
  options.add_options( "general" )
      ( "o, output", "Write to this file instead of stdout",
        cxxopts::value< std::string >()->default_value( DEFAULT_OUTPUT ) )
      ( "g, graph", "Corresponding graph file (vg or gfa)",
        cxxopts::value< std::string >() )
      ( "l, seed-length", "Seed length",
        cxxopts::value< unsigned int >() )
      ( "e, step-size", "Step size",
        cxxopts::value< unsigned int >()->default_value( DEFAULT_STEP_SIZE ) )
      ( "h, help", "Print this message and exit" )
      ;

  options.add_options( "convert" )
      ( "P, from-proto", "Consider input to be a protobuf stream" )
      ( "p, to-proto", "Consider output to be a protobuf stream" )
      ( "N, from-native", "Consider input format to be PSI native seriliasation format" )
      ( "n, to-native", "Write output in PSI native seriliasation format" )
      ( "J, from-json", "Consider input format to be JSON" )
      ( "j, to-json", "Write output in JSON" )
      ;

  options.add_options( "positional" )
      ( "command", "Operation type", cxxopts::value< std::string >() )
      ( "index-prefix", "PSI index prefix", cxxopts::value< std::string >() )
      ;
  options.parse_positional( { "command", "index-prefix" } );
}

cxxopts::ParseResult
parse_opts( cxxopts::Options& options, int& argc, char**& argv )
{
  auto result = options.parse( argc, argv );

  if ( !result.count( "command" ) ) {  // no command specified
    options.positional_help( "COMMAND" );
    auto help_message = ( options.help( { "general" } )
                          + "\n COMMANDS:\n"
                          + "  convert\tConvert between different file formats\n" );
    if ( result.count( "help" ) )
    {
      std::cout << help_message << std::endl;
      throw EXIT_SUCCESS;
    }
    std::cerr << help_message << "\n" /* extra vertical space */ << std::endl;
    throw cxxopts::exceptions::parsing( "No command specified" );
  }
  else if ( result[ "command" ].as< std::string >() == "convert" ) {  // convert
    options.custom_help( "convert [OPTION...]" );
    options.positional_help( "INDEX-PREFIX" );
    if ( result.count( "help" ) ) {
      std::cout << options.help( { "general", "convert" } ) << std::endl;
      throw EXIT_SUCCESS;
    }

    short int sum_fargs = 0;
    sum_fargs += ( result.count( "from-proto" ) ? 1 : 0 );
    sum_fargs += ( result.count( "from-native" ) ? 1 : 0 );
    sum_fargs += ( result.count( "from-json" ) ? 1 : 0 );
    if ( sum_fargs > 1 ) throw cxxopts::exceptions::parsing( "Only one input format specifier can be used" );

    sum_fargs = 0;
    sum_fargs += ( result.count( "to-proto" ) ? 1 : 0 );
    sum_fargs += ( result.count( "to-native" ) ? 1 : 0 );
    sum_fargs += ( result.count( "to-json" ) ? 1 : 0 );
    if ( sum_fargs > 1 ) throw cxxopts::exceptions::parsing( "Only one output format specifier can be used" );

    if ( result.count( "from-proto" ) && result.count( "to-proto" ) ) {
      throw cxxopts::exceptions::parsing( "No conversion needed" );
    }

    if ( result.count( "from-native" ) && result.count( "to-native" ) ) {
      throw cxxopts::exceptions::parsing( "No conversion needed" );
    }

    if ( result.count( "from-json" ) && result.count( "to-json" ) ) {
      throw cxxopts::exceptions::parsing( "No conversion needed" );
    }
  }
  else {
    throw cxxopts::exceptions::parsing( "Unknown command '" +
                                         result[ "command" ].as< std::string >() + "'" );
  }

  /* Verifying general arguments */
  if ( ! result.count( "graph" ) ) {
    throw cxxopts::exceptions::parsing( "Graph file must be specified" );
  }
  if ( ! readable( result[ "graph" ].as< std::string >() ) ) {
    throw cxxopts::exceptions::parsing( "Graph file not found" );
  }

  /* Verifying positional arguments */
  if ( !result.count( "index-prefix" ) ) {
    throw cxxopts::exceptions::parsing( "No PSI index prefix has been specified" );
  }
  if ( !readable( result[ "index-prefix" ].as< std::string >() ) ) {
    throw cxxopts::exceptions::parsing( "PSI index file not found" );
  }

  if ( !result.count( "seed-length" ) ) {
    throw cxxopts::exceptions::parsing( "No seed length has been specified" );
  }

  return result;
}

  std::vector< psi::Position<> >
read_proto( const std::string& prefix, unsigned int seed_len, unsigned int step_size )
{
  std::string filepath = psi::SeedFinder<>::get_sloci_filepath( prefix, seed_len, step_size );
  std::ifstream ifs( filepath, std::ifstream::in | std::ifstream::binary );
  if ( !ifs ) throw std::runtime_error( "Index file not found" );

  std::vector< psi::Position<> > retval;

  std::function< void( vg::Position& ) > push_back =
    [&retval]( vg::Position& pos ) {
      typedef psi::Position<>::id_type id_type;
      typedef psi::Position<>::offset_type offset_type;
      retval.push_back( { ( id_type )( pos.node_id() ), ( offset_type )( pos.offset() ) } );
    };

  try {
    stream::for_each( ifs, push_back );
  }
  catch ( const std::runtime_error& ) {
    throw std::runtime_error( "Unexpected error while reading input protobuf stream" );
  }

  return retval;
}

  std::vector< psi::Position<> >
read_native( const std::string& prefix, unsigned int seed_len, unsigned int step_size )
{
  std::string filepath = psi::SeedFinder<>::get_sloci_filepath( prefix, seed_len, step_size );
  std::ifstream ifs( filepath, std::ifstream::in | std::ifstream::binary );
  if ( !ifs ) throw std::runtime_error( "Index file not found" );

  std::vector< psi::Position<> > retval;
  psi::SeedFinder<>::deserialize_starts( ifs, retval );
  return retval;
}

  std::vector< psi::Position<> >
read_json( const std::string& prefix, unsigned int seed_len, unsigned int step_size )
{
  throw std::runtime_error( "Parsing from JSON is not implemented" );
}

  void
write_proto( std::ostream& ofs, const std::vector< psi::Position<> >& loci )
{
  std::function< vg::Position( uint64_t ) > lambda =
    [&loci]( uint64_t i ) -> vg::Position {
      vg::Position pos;
      pos.set_node_id( loci.at( i ).node_id() );
      pos.set_offset( loci.at( i ).offset() );
      return pos;
    };

  try {
    stream::write< vg::Position >( ofs, loci.size(), lambda );
  }
  catch ( const std::runtime_error& ) {
    throw std::runtime_error( "Unexpected error while writing output protobuf stream" );
  }
}

  void
write_native( std::ostream& ofs, const std::vector< psi::Position<> >& loci )
{
  psi::SeedFinder<>::serialize_starts( ofs, loci );
}

  void
write_json( std::ostream& ofs, const std::vector< psi::Position<> >& loci )
{
  throw std::runtime_error( "Writing to JSON is not implemented" );
}

  void
convert( cxxopts::ParseResult& res )
{
  typedef gum::SeqGraph< gum::Succinct > graph_type;

  // Fetching input parameters
  std::string output = res[ "output" ].as< std::string >();
  std::string graph_path = res[ "graph" ].as< std::string >();
  std::string index_prefix = res[ "index-prefix" ].as< std::string >();
  unsigned int seed_len = res[ "seed-length" ].as< unsigned int >();
  unsigned int step_size = res[ "step-size" ].as< unsigned int >();
  bool from_proto = res[ "from-proto" ].as< bool >();
  bool to_proto = res[ "to-proto" ].as< bool >();
  bool from_native = res[ "from-native" ].as< bool >();
  bool to_native = res[ "to-native" ].as< bool >();
  bool from_json = res[ "from-json" ].as< bool >();
  bool to_json = res[ "to-json" ].as< bool >();

  int sum_fargs = 0;
  sum_fargs += ( from_proto ? 1 : 0 );
  sum_fargs += ( from_native ? 1 : 0 );
  sum_fargs += ( from_json ? 1 : 0 );
  if ( sum_fargs == 0 ) from_proto = true;

  sum_fargs = 0;
  sum_fargs += ( to_proto ? 1 : 0 );
  sum_fargs += ( to_native ? 1 : 0 );
  sum_fargs += ( to_json ? 1 : 0 );
  if ( sum_fargs == 0 ) to_native = true;

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

  auto parse_vg = []( std::istream& in ) -> vg::Graph {
    vg::Graph merged;
    std::function< void( vg::Graph& ) > handle_chunks =
      [&]( vg::Graph& other ) {
        gum::util::merge_vg( merged, static_cast< vg::Graph const& >( other ) );
      };
    stream::for_each( in, handle_chunks );
    return merged;
  };

  gum::ExternalLoader< vg::Graph > loader{ parse_vg };
  gum::util::load( graph, graph_path, loader, true );

  std::vector< psi::Position<> > sloci;
  if ( from_proto ) sloci = read_proto( index_prefix, seed_len, step_size );
  else if ( from_native ) sloci = read_native( index_prefix, seed_len, step_size );
  else if ( from_json ) sloci = read_json( index_prefix, seed_len, step_size );
  else assert( false );

  if ( to_proto ) write_proto( ost, sloci );
  else if ( to_native ) write_native( ost, sloci );
  else if ( to_json ) write_json( ost, sloci );
  else assert( false );
}

  int
main( int argc, char* argv[] )
{
  cxxopts::Options options( argv[0], LONG_DESC );
  config_parser( options );

  try {
    auto res = parse_opts( options, argc, argv );
    std::string command = res[ "command" ].as< std::string >();

    if ( command == "convert" ) {
      convert( res );
    }
    else {
      // should not reach here!
      assert( false );
    }
  }
  catch ( const cxxopts::exceptions::exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch ( const int& rv ) {
    return rv;
  }

  return EXIT_SUCCESS;
}

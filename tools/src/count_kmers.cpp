/**
 *    @file  count_kmers.cpp
 *   @brief  Count kmers in a graph.
 *
 *  This is a tool for counting k-mers in a variation graph.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Thu May 09, 2019  17:03
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2019, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "cxxopts/cxxopts.h"

#include "vargraph.h"
#include "logger.h"
#include "utils.h"

using namespace std;
using namespace grem;

const char * const LONG_DESC = "Count k-mers in a variation graph";

  void
config_parser( cxxopts::Options& options )
{
  options.positional_help( "GRAPH" );
  options.add_options()
    ( "k, length", "Value of k", cxxopts::value< unsigned int >() )
    ( "F, forward", "Only count k-mer in forward strand", cxxopts::value<bool>()->default_value( "false" ) )
    ( "h, help", "Print this message and exit" )
    ;

  options.add_options( "positional" )
    ( "graph", "Graph file in `xg` or `vg` format", cxxopts::value< string >() )
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

  if ( !result.count( "graph" ) ) {
    throw cxxopts::exceptions::parsing( "Graph file must be provided" );
  }

  if ( !result.count( "length" ) ) {
    throw cxxopts::exceptions::parsing( "k-mer length must be specified" );
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
    /* Configure loggers */
    config_logger( false, false, false, true, true, "" );

    string graph_path = res[ "graph" ].as< string >();
    unsigned int k = res[ "length" ].as< unsigned int >();
    bool forward = res[ "forward" ].as< bool >();

    VarGraph vargraph;
    ifstream ifs( graph_path, ifstream::in | ifstream::binary );
    if ( ends_with( graph_path, ".vg" ) ) {
      vargraph.from_stream( ifs );
    }
    else {
      vargraph.load( ifs );
    }

    cout << count_kmers( vargraph, k, forward ) << endl;
  }
  catch ( const cxxopts::exceptions::exception& e ) {
    cerr << "Error: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch ( const int& rv ) {
    return rv;
  }

  return EXIT_SUCCESS;
}

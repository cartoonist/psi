/**
 *    @file  loci_stats.cc
 *   @brief  Report statistics about starting loci.
 *
 *  This program reads the starting loci file as a part of path index and report
 *  some statistics.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Jun 12, 2018  13:58
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdlib>
#include <iostream>
#include <string>

#include "cxxopts/cxxopts.h"

#include "mapper.h"
#include "traverser.h"
#include "utils.h"

using namespace std;
using namespace grem;


const char * const LONG_DESC = "Report statistics about starting loci";

  void
config_parser( cxxopts::Options& options )
{
  options.positional_help( "INDEX_PREFIX" );
  options.add_options()
    ( "l, seed-length", "Seed length", cxxopts::value< unsigned int >() )
    ( "e, step-size", "Step size", cxxopts::value< unsigned int >() )
    ( "s, start-node", "Start node", cxxopts::value< int64_t >()->default_value( "1" ) )
    ( "t, end-node", "End node [0 means last node]", cxxopts::value< int64_t >()->default_value( "0" ) )
    ( "n, number", "Number of loci to be reported [0 means all]",
      cxxopts::value< unsigned int >()->default_value( "0" ) )
    ( "h, help", "Print this message and exit" )
    ;

  options.add_options( "positional" )
    ( "prefix", "Path index prefix", cxxopts::value< string >() )
    ;
  options.parse_positional( { "prefix" } );
}

  cxxopts::ParseResult
parse_opts( cxxopts::Options& options, int& argc, char**& argv )
{
  auto result = options.parse( argc, argv );

  if ( result.count( "help" ) ) {
    cout << options.help( { "" } ) << endl;
    throw EXIT_SUCCESS;
  }

  if ( ! result.count( "prefix" ) ) {
    throw cxxopts::OptionParseException( "Index prefix must be specified" );
  }
  if ( !readable( result[ "prefix" ].as< string >() ) ) {
    throw cxxopts::OptionParseException( "Index file not found" );
  }

  if ( ! result.count( "seed-length" ) ) {
    throw cxxopts::OptionParseException( "Seed length must be specified" );
  }

  if ( ! result.count( "step-size" ) ) {
    throw cxxopts::OptionParseException( "Step size must be specified" );
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

    typedef seqan::Index< Dna5QStringSet<>, seqan::IndexWotd<> > TIndex;
    typedef typename Traverser< TIndex, BFS, ExactMatching >::Type TTraverser;
    typedef Mapper< TTraverser > TMapper;
    unsigned int seedlen = res["seed-length"].as< unsigned int >();
    unsigned int stepsize = res["step-size"].as< unsigned int >();
    int64_t start = res["start-node"].as< int64_t >();
    int64_t end = res["end-node"].as< int64_t >();
    unsigned int num = res["number"].as< unsigned int >();
    TMapper mapper( NULL, seedlen );
    if ( ! mapper.open_starts( res["prefix"].as< string >(), seedlen, stepsize ) )
      throw cxxopts::OptionException( "Index file seems corrupted" );

    cout << "Number of loci: " << mapper.get_starting_loci().size() << endl;

    if ( ! mapper.get_starting_loci().empty() ) {
      unsigned int i = 1;
      bool found = false;
      cout << endl
        << "---------------" << endl
        << "num: id, offset" << endl
        << "---------------" << endl;
      for ( const auto& locus : mapper.get_starting_loci() ) {
        if ( !found && locus.node_id() >= start ) found = true;
        if ( found ){
          if ( ( num && i > num ) || ( end && locus.node_id() > end ) ) break;
          cout << i++ << ": " << locus.node_id() << ", " << locus.offset() << endl;
        }
      }
      cout << "---------------" << endl;
    }
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

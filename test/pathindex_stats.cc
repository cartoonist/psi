/**
 *    @file  pathindex_stats.cc
 *   @brief  Report statistics about path index.
 *
 *  This program reads the path index and report some statistics about the paths.
 *
 *  NOTE: In order to inspect the starting loci, see `loci_stats.cc`.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Wed Jun 13, 2018  10:54
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
#include <functional>

#include "cxxopts/cxxopts.h"
#include "stream.hpp"

#include "vargraph.h"
#include "pathindex.h"
#include "traverser.h"
#include "mapper.h"
#include "utils.h"

using namespace std;
using namespace grem;


const char * const LONG_DESC = "Report statistics about path index";

  void
config_parser( cxxopts::Options& options )
{
  options.positional_help( "INDEX_PREFIX" );
  options.add_options()
    ( "l, seed-length", "Seed length", cxxopts::value< unsigned int >() )
    ( "e, step-size", "Step size", cxxopts::value< unsigned int >() )
    ( "L, no-loci", "Do not include starting loci as SNP", cxxopts::value<bool>()->default_value( "false" ) )
    ( "o, output", "Output GAM file", cxxopts::value< string >()->default_value( "pathindex.gam" ) )
    ( "g, graph", "Corresponding graph (vg or xg)", cxxopts::value< string >() )
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
  if ( ! readable( result[ "prefix" ].as< string >() ) ) {
    throw cxxopts::OptionParseException( "Index file not found" );
  }

  if ( ! result.count( "graph" ) ) {
    throw cxxopts::OptionParseException( "Graph must be specified" );
  }
  if ( ! readable( result[ "graph" ].as< string >() ) ) {
    throw cxxopts::OptionParseException( "Graph file not found" );
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

    typedef seqan::Index< Dna5QStringSet< grem::Dependent >, seqan::IndexWotd<> > TIndex;
    typedef typename Traverser< TIndex, BFS, ExactMatching >::Type TTraverser;
    typedef Mapper< TTraverser > TMapper;

    string graph_path = res[ "graph" ].as< string >();
    string pindex_prefix = res[ "prefix" ].as< string >();
    string output = res[ "output" ].as< string >();
    unsigned int seedlen = res["seed-length"].as< unsigned int >();
    unsigned int stepsize = res["step-size"].as< unsigned int >();
    bool noloci = res["no-loci"].as< bool >();

    VarGraph vargraph;
    ifstream ifs( graph_path, ifstream::in | ifstream::binary );
    if ( ends_with( graph_path, ".vg" ) ) {
      vargraph.from_stream( ifs );
    }
    else {
      vargraph.load( ifs );
    }

    PathIndex< VarGraph, DiskString, grem::FMIndex<>, Forward > pindex;
    if ( ! pindex.load( pindex_prefix, &vargraph ) ) {
      throw cxxopts::OptionException( "Index file seems corrupted" );
    }

    TMapper mapper( &vargraph, seedlen );
    if ( ! mapper.open_starts( res["prefix"].as< string >(), seedlen, stepsize ) )
      throw cxxopts::OptionException( "Index file seems corrupted" );

    auto nofpaths = pindex.get_paths_set().size();
    auto totseqlen = getFibre( pindex.index, seqan::FibreText() ).raw_length();
    cout << "Number of paths: " << nofpaths << endl;
    cout << "Total sequence length: " << totseqlen << endl;
    cout << endl;

    vector< vg::Alignment > paths;
    vg::Alignment p;
    string pathname;
    for ( size_t i = 0; i < nofpaths; ++i ) {
      pathname = "path" + to_string( i + 1 );
      p.set_name( pathname );
      cout << "\rConverted " << to_string( i + 1 ) << "/" << to_string( nofpaths )
           << " paths to vg::Path.";
      if ( noloci ) {
        grem::convert( pindex.get_paths_set()[i], p.mutable_path() );
      }
      else {
        grem::convert( pindex.get_paths_set()[i], p.mutable_path(),
            mapper.get_starting_loci() );
      }
      p.mutable_path()->set_name( pathname );
      paths.push_back( move( p ) );
      p.Clear();
    }

    ofstream ofs( output, ofstream::out | ofstream::binary );
    function< vg::Alignment( uint64_t ) > lambda =
      [&paths]( uint64_t i ) { return paths.at( i ); };
    cout << "\nWriting all paths to a GAM file... ";
    stream::write( ofs, paths.size(), lambda );
    cout << "Done." << endl;
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

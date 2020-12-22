/**
 *    @file  pathindex_stats.cpp
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
#include "logger.h"
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
    ( "t, context", "Context size", cxxopts::value< unsigned int >()->default_value( "0" ) )
    ( "L, no-loci", "Do not include starting loci as SNP", cxxopts::value<bool>()->default_value( "false" ) )
    ( "F, forward", "Set if the path index is NOT from reversed sequence", cxxopts::value<bool>()->default_value( "false" ) )
    ( "m, max-nodes", "Maximum number of nodes allowed in a `vg::Graph` message", cxxopts::value< unsigned int >()->default_value( "1000" ) )
    ( "o, output", "Output GAM/vg file", cxxopts::value< string >()->default_value( "pathindex.gam" ) )
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

  if ( !result.count( "prefix" ) ) {
    throw cxxopts::OptionParseException( "Index prefix must be specified" );
  }
  if ( !readable( result[ "prefix" ].as< string >() ) ) {
    throw cxxopts::OptionParseException( "Index file not found" );
  }

  if ( !result.count( "graph" ) ) {
    throw cxxopts::OptionParseException( "Graph must be specified" );
  }
  if ( !readable( result[ "graph" ].as< string >() ) ) {
    throw cxxopts::OptionParseException( "Graph file not found" );
  }

  if ( !result.count( "seed-length" ) ) {
    throw cxxopts::OptionParseException( "Seed length must be specified" );
  }

  if ( !result.count( "step-size" ) ) {
    throw cxxopts::OptionParseException( "Step size must be specified" );
  }

  return result;
}

template< typename TGraph, typename TMapper >
    void
  to_gam( PathSet< Path< TGraph, Compact > > const& pathset, TMapper const& mapper,
      bool noloci, string const& output )
  {
    vector< vg::Alignment > paths;
    vg::Alignment p;
    string pathname;
    for ( size_t i = 0; i < pathset.size(); ++i ) {
      pathname = "path" + to_string( i + 1 );
      p.set_name( pathname );
      cout << "\rConverted " << to_string( i + 1 ) << "/" << to_string( pathset.size() )
           << " paths to vg::Path.";
      if ( noloci ) {
        grem::convert( pathset[i], p.mutable_path() );
      }
      else {
        grem::convert( pathset[i], p.mutable_path(), mapper.get_starting_loci() );
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

template< typename TGraph, typename TMapper >
    void
  to_vg( PathSet< Path< TGraph, Compact > > const& pathset, TGraph const* vargraph,
      TMapper const& /*mapper*/, bool /*noloci*/, unsigned int max_nodes,
      string const& output )
  {
    using nodeid_type = typename TGraph::nodeid_type;

    vector< nodeid_type > nodes;
    vector< tuple< nodeid_type, nodeid_type, char > > edges;
    vector< vg::Graph > graphset;
    ofstream ofs( output, ofstream::out | ofstream::binary );

    function< vg::Graph( uint64_t ) > graph_at =
      [&graphset]( uint64_t i ) { return graphset.at( i ); };

    function< void( vg::Graph& ) > accumulate =
      [&graphset]( vg::Graph& g ) {
        graphset.push_back( std::move( g ) );
      };

    cout << "Calculating the graph induced by paths set... " << endl;
    grem::induced_graph( pathset.begin(), pathset.end(), nodes, edges );
    cout << "Converting the induced graph to a set of `vg::Graph` messages... " << endl;
    vargraph->induced_graph( nodes.begin(), nodes.end(), edges.begin(), edges.end(),
        accumulate, max_nodes );
    cout << "Writing the induced graph to a vg file... " << endl;
    stream::write( ofs, graphset.size(), graph_at );
    cout << "Done." << endl;
  }

template< typename TSequenceDirection, typename TGraph, typename TMapper >
    void
  inspect_pathindex( const TGraph* vargraph, TMapper& mapper,
      const string& pindex_prefix, const string& output, unsigned int ctx,
      unsigned int seedlen, unsigned int stepsize, unsigned int max_nodes, bool noloci )
  {
    PathIndex< TGraph, DiskString, grem::FMIndex<>, TSequenceDirection > pindex( ctx, false );
    if ( !pindex.load( pindex_prefix, vargraph ) ) {
      throw cxxopts::OptionException( "Index file seems corrupted" );
    }

    if ( !mapper.open_starts( pindex_prefix, seedlen, stepsize ) ) {
      throw cxxopts::OptionException( "Starting loci file seems corrupted" );
    }

    auto nofpaths = pindex.get_paths_set().size();
    cout << "Number of paths: " << nofpaths << endl;
    auto totseqlen = getFibre( pindex.index, seqan::FibreText() ).raw_length();
    cout << "Total sequence length: " << totseqlen << endl;
    cout << "Context size: " << pindex.get_context() << endl;
    cout << "Number of uncovered loci: " << mapper.get_starting_loci().size() << endl;
    cout << "Number of total loci: " << vargraph->get_total_nof_loci() << endl;
//    auto nofuckmers = mapper.nof_uncovered_kmers( pindex.get_paths_set(), seedlen );
//    cout << "Number of uncovered k-mers: " << nofuckmers << endl;
    cout << endl;

    if ( ends_with( output, ".vg" ) ) {
      to_vg( pindex.get_paths_set(), vargraph, mapper, noloci, max_nodes, output );
    }
    else if ( ends_with( output, ".gam" ) ) {
      to_gam( pindex.get_paths_set(), mapper, noloci, output );
    }
    else {
      runtime_error( "Unsupported output format" );
    }
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

    typedef seqan::Index< Dna5QStringSet<>, seqan::IndexWotd<> > TIndex;
    typedef typename Traverser< TIndex, BFS, ExactMatching >::Type TTraverser;
    typedef Mapper< TTraverser > TMapper;

    string graph_path = res[ "graph" ].as< string >();
    string pindex_prefix = res[ "prefix" ].as< string >();
    string output = res[ "output" ].as< string >();
    unsigned int context = res["context"].as< unsigned int >();
    bool noloci = res["no-loci"].as< bool >();
    bool forward = res["forward"].as< bool >();
    unsigned int max_nodes = res[ "max-nodes" ].as< unsigned int >();
    unsigned int seedlen = res["seed-length"].as< unsigned int >();
    unsigned int stepsize = res["step-size"].as< unsigned int >();

    VarGraph vargraph;
    ifstream ifs( graph_path, ifstream::in | ifstream::binary );
    if ( ends_with( graph_path, ".vg" ) ) {
      vargraph.from_stream( ifs );
    }
    else {
      vargraph.load( ifs );
    }

    TMapper mapper( &vargraph, seedlen );

    if ( forward ) {
      inspect_pathindex< Forward >( &vargraph, mapper, pindex_prefix, output, context,
          seedlen, stepsize, max_nodes, noloci );
    }
    else {
      inspect_pathindex< Reversed >( &vargraph, mapper, pindex_prefix, output, context,
          seedlen, stepsize, max_nodes, noloci );
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

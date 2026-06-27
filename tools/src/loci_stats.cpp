/**
 *    @file  loci_stats.cpp
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

#include <cxxopts.hpp>
#include <gum/graph.hpp>
#include <gum/io_utils.hpp>
#include <psi/seed_finder.hpp>
#include <psi/utils.hpp>

#include "vg/vg.pb.h"
#include "vg/stream.hpp"


using namespace psi;

/* ====== Constants ====== */
constexpr const char* const LONG_DESC = "Report statistics about starting loci";

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
      ( "g, graph", "Corresponding graph (vg or gfa)", cxxopts::value< std::string >() )
      ( "h, help", "Print this message and exit" )
      ;

  options.add_options( "positional" )
      ( "prefix", "Path index prefix", cxxopts::value< std::string >() )
      ;
  options.parse_positional( { "prefix" } );
}

  cxxopts::ParseResult
parse_opts( cxxopts::Options& options, int& argc, char**& argv )
{
  auto result = options.parse( argc, argv );

  if ( result.count( "help" ) ) {
    std::cout << options.help( { "" } ) << std::endl;
    throw EXIT_SUCCESS;
  }

  if ( ! result.count( "prefix" ) ) {
    throw cxxopts::exceptions::parsing( "Index prefix must be specified" );
  }
  if ( !readable( result[ "prefix" ].as< std::string >() ) ) {
    throw cxxopts::exceptions::parsing( "Index file not found" );
  }

  if ( !result.count( "graph" ) ) {
    throw cxxopts::exceptions::parsing( "Graph must be specified" );
  }
  if ( !readable( result[ "graph" ].as< std::string >() ) ) {
    throw cxxopts::exceptions::parsing( "Graph file not found" );
  }

  if ( ! result.count( "seed-length" ) ) {
    throw cxxopts::exceptions::parsing( "Seed length must be specified" );
  }

  if ( ! result.count( "step-size" ) ) {
    throw cxxopts::exceptions::parsing( "Step size must be specified" );
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

    typedef SeedFinder<> TFinder;

    std::string graph_path = res[ "graph" ].as< std::string >();
    unsigned int seedlen = res["seed-length"].as< unsigned int >();
    unsigned int stepsize = res["step-size"].as< unsigned int >();
    int64_t start = res["start-node"].as< int64_t >();
    int64_t end = res["end-node"].as< int64_t >();
    unsigned int num = res["number"].as< unsigned int >();
    std::string prefix = res["prefix"].as< std::string >();

    auto parse_vg = []( std::istream& in ) -> vg::Graph {
      vg::Graph merged;
      std::function< void( vg::Graph& ) > handle_chunks =
        [&]( vg::Graph& other ) {
          gum::util::merge_vg( merged, static_cast< vg::Graph const& >( other ) );
        };
      stream::for_each( in, handle_chunks );
      return merged;
    };

    gum::SeqGraph< gum::Succinct > graph;
    gum::ExternalLoader< vg::Graph > loader{ parse_vg };
    gum::util::load( graph, graph_path, loader, true );
    std::string sort_status = gum::util::ids_in_topological_order( graph ) ? "" : "not ";
    std::cout << "Input graph node IDs are " << sort_status << "in topological sort order."
              << std::endl;
    TFinder finder( graph, seedlen );

    if ( ! finder.open_starts( prefix, seedlen, stepsize ) )
      throw cxxopts::exceptions::exception( "Index file seems corrupted" );

    std::cout << "Number of loci: " << finder.get_starting_loci().size() << std::endl;

    if ( ! finder.get_starting_loci().empty() ) {
      unsigned int i = 1;
      bool found = false;
      std::cout << std::endl
                << "---------------" << std::endl
                << "num: id, offset" << std::endl
                << "---------------" << std::endl;
      for ( const auto& locus : finder.get_starting_loci() ) {
        auto id = graph.coordinate_id( locus.node_id() );
        if ( !found && id >= start ) found = true;
        if ( found ){
          if ( ( num && i > num ) || ( end && id > end ) ) break;
          std::cout << i++ << ": " << id << ", " << locus.offset() << std::endl;
        }
      }
      std::cout << "---------------" << std::endl;
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

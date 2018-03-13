/**
 *    @file  tests_traverser.cc
 *   @brief  Test traverser module.
 *
 *  Test scenarios for traverser module.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Fri Mar 17, 2017  00:44
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

// :TODO:Fri Mar 17 00:46:\@cartoonist: traverser should be refactored, also this file.

#include <string>
#include <vector>

#include "tests_base.h"
#include "traverser.h"
#include "vargraph.h"
#include "logger.h"

INITIALIZE_EASYLOGGINGPP

using namespace grem;


// :TODO:Fri Mar 17 00:52:\@cartoonist: Mapper -> refactor
// :TODO:Fri Mar 17 00:52:\@cartoonist: traverser -> refactor
// :TODO:Sun Jun 11 15:14:\@cartoonist: Fix this test scenario.
//SCENARIO ( "Adding starting points for seed finding using Mapper", "[traverser]" )
//{
//  GIVEN ( "A small variation graph and a graph seed finder" )
//  {
//    typedef PathTraverser < seqan::IndexWotd<> > TPathTraverser;
//    typedef Mapper < TPathTraverser > TMapper;
//
//    std::string vgpath = _testdir + "/data/small/x.vg";
//    VarGraph vargraph ( vgpath.c_str() );
//    TMapper mapper ( vargraph );
//
//    unsigned int stepsize = 50;
//    WHEN ( "Step size is " + std::to_string ( stepsize ) + " and the exclude nodes list is non-empty." )
//    {
//      std::vector < VarGraph::NodeCoverage > excluded_nodes;
//      excluded_nodes.insert ( 1 );
//      mapper.add_all_loci ( excluded_nodes, stepsize );
//
//      THEN ( "all node should be picked up as starting points except the excluded ones." )
//      {
//        const std::vector < vg::Position > &st_points = mapper.get_starting_points ();
//
//        unsigned int counter = 2;
//        for ( auto const &pos : st_points ) {
//          REQUIRE ( pos.node_id() == counter );
//          ++counter;
//        }
//
//        REQUIRE ( --counter == vargraph.nodes_size() );
//      }
//    }
//  }
//}

SCENARIO ( "Serialize/deserialize paths nodes coverage into/from the file", "[traverser]" )
{
  GIVEN ( "Nodes coverage of two paths" )
  {
    std::vector< VarGraph::NodeCoverage > paths_node_coverage;
    VarGraph::NodeCoverage covered_nodes;
    unsigned int paths_num = 2;
    std::string file_path = "/tmp/test_path_coverage";

    for ( unsigned int i = 0; i < paths_num; ++i ) {
      for ( VarGraph::nodeid_type j = 100*(i+1); j < 100*(i+1) + 12; ++j ) {
        covered_nodes.insert ( j );
      }
      paths_node_coverage.push_back ( covered_nodes );
      covered_nodes.clear();
    }

    WHEN ( "Serialize it to the file" )
    {
      save_paths_coverage ( paths_node_coverage, file_path );

      THEN ( "Deserializing should yield the same coverage" )
      {
        std::vector< VarGraph::NodeCoverage > paths_coverage_reloaded;
        load_paths_coverage ( paths_coverage_reloaded, file_path, paths_num );

        std::vector< VarGraph::nodeid_type > sorted_node_ids;
        REQUIRE ( paths_coverage_reloaded.size() == paths_num );
        for ( unsigned int i = 0; i < paths_num; ++i ) {
          REQUIRE ( paths_coverage_reloaded[i].size() == 12 );
          for ( const auto &node_id : paths_coverage_reloaded[i] ) {
            sorted_node_ids.push_back ( node_id );
          }
          unsigned int j = 0;
          std::sort ( sorted_node_ids.begin(), sorted_node_ids.end() );
          for ( const auto &node_id : sorted_node_ids ) {
            REQUIRE ( node_id == 100*(i+1) + j );
            ++j;
          }
          sorted_node_ids.clear();
        }
      }
    }
  }
}

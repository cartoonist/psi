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


// :TODO:Fri Mar 17 00:52:\@cartoonist: GraphTraverser -> refactor
// :TODO:Fri Mar 17 00:52:\@cartoonist: traverser -> refactor
SCENARIO ( "Adding starting points for seed finding using GraphTraverser", "[traverser]" )
{
  GIVEN ( "A small variation graph and a graph seed finder" )
  {
    typedef PathTraverser < seqan::IndexWotd<> > TPathTraverser;
    typedef GraphTraverser < TPathTraverser > TGraphTraverser;

    std::string vgpath = _testdir + "/data/small/x.vg";
    VarGraph vargraph ( vgpath.c_str() );
    TGraphTraverser gtraverser ( vargraph );

    unsigned int stepsize = 50;
    WHEN ( "Step size is " + std::to_string ( stepsize ) + " and the exclude nodes list is non-empty." )
    {
      std::unordered_set < VarGraph::NodeID > excluded_nodes;
      excluded_nodes.insert ( 1 );
      gtraverser.add_all_loci ( stepsize, &excluded_nodes );

      THEN ( "all node should be picked up as starting points except the excluded ones." )
      {
        const std::vector < vg::Position > &st_points = gtraverser.get_starting_points ();

        unsigned int counter = 2;
        for ( auto const &pos : st_points ) {
          REQUIRE ( pos.node_id() == counter );
          ++counter;
        }

        REQUIRE ( --counter == vargraph.nodes_size() );
      }
    }
  }
}

/**
 *    @file  tests_bfs.cc
 *   @brief  Test scenarios for BFS graph iterator.
 *
 *  It contains all test scenarios for BFS graph iterator.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Wed Jan 18, 2017  16:07
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <string>

#include "tests_base.h"
#include "vargraph.h"
#include "logger.h"

INITIALIZE_EASYLOGGINGPP

using namespace grem;


SCENARIO ( "Variation graph breadth-first traverse (BFS)", "[graph][iterator]" )
{
  GIVEN ( "A small variation graph" )
  {
    std::string vgpath = _testdir + "/data/small/x.vg";
    VarGraph vargraph(vgpath.c_str());

    WHEN ( "traverse the graph using BFS graph iterator" )
    {
      seqan::Iterator < VarGraph, BFS<> >::Type bfs_itr (vargraph);

      THEN ( "nodes should be traversed in BFS order" )
      {
        VarGraph::NodeID truth = 1;
        while ( !at_end ( bfs_itr ) ) {
          REQUIRE ( *bfs_itr == truth );  // The graph is such that its BFS is in order.
          ++truth;
          ++bfs_itr;
        }
      }
    }
  }
}

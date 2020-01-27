/**
 *    @file  tests_pathset.cc
 *   @brief  Test PathSet class.
 *
 *  Test scenarios for PathSet class.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sun Sep 17, 2017  22:41
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <fstream>
#include <string>

#include "tests_base.h"
#include "vargraph.h"
#include "pathset.h"
#include "index.h"
#include "logger.h"


using namespace grem;

SCENARIO ( "Serialize/deserialize paths set into/from the file", "[pathset]" )
{
  GIVEN ( "Two paths from a small graph" )
  {
    typedef seqan::IndexEsa<> TIndexSpec;

    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream gifs( vgpath.c_str() );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );
    PathSet< TIndexSpec > paths_set;

    unsigned int paths_num = 2;
    std::string file_path = "/tmp/test_pathset";

    paths_set.reserve( paths_num );
    for ( unsigned int i = 0; i < paths_num; ++i ) {
      Path<> path( &vargraph );
      for ( VarGraph::nodeid_type j = 3+i; j <= 210; j+=(i+1)*4 ) {
        add_node( path, j );
      }
      paths_set.add_path( std::move( path ) );
    }

    WHEN ( "Serialize it to the file" )
    {
      paths_set.save( file_path );

      THEN ( "Deserializing should yield the same paths" )
      {
        PathSet< TIndexSpec > loaded_paths;
        loaded_paths.load( file_path, &vargraph );
        REQUIRE ( loaded_paths.size() == paths_num );
        REQUIRE( length( loaded_paths.paths_set.at(0) ) == 52 );
        REQUIRE( length( loaded_paths.paths_set.at(1) ) == 26 );

        for ( unsigned int i = 0; i < paths_num; ++i ) {
          unsigned int j = 0;
          for ( const auto &node_id : loaded_paths.paths_set.at(i).get_nodes() ) {
            REQUIRE ( node_id == 3+i+(i+1)*4*j );
            ++j;
          }
        }
      }
    }
  }
}

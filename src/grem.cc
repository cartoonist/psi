/*
 * =====================================================================================
 *
 * Filename: grem.cpp
 *
 * Created: Tue Nov 08, 2016  16:48
 * Last modified: Sun Nov 27, 2016  19:44
 *
 * Description: GREM main function.
 *
 * Copyright (c) 2016, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari, <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#include <cstdlib>
#include <iostream>
#include <ios>
#include <fstream>
#include <string>
#include <functional>

#include <seqan/seq_io.h>

#include "vargraph.h"
#include "traverser.h"
#include "types.h"
#include "release.h"

//#undef NDEBUG
#include <easyloggingpp/src/easylogging++.h>

using namespace seqan;
using namespace grem;

// TODO: Documentation.
// TODO: Memory footprint.
// TODO: Logging: replace all iostream calls with logging APIs.
// TODO: Logging messages.

INITIALIZE_EASYLOGGINGPP


int main(int argc, char *argv[])
{
  START_EASYLOGGINGPP(argc, argv);

  el::Configurations log_conf;
  log_conf.setToDefault();
  el::Loggers::reconfigureLogger("default", log_conf);

  // Verify that the version of the library that we linked against is
  // compatible with the version of the headers we compiled against.
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  // TODO: get command-line arguments. Now, assuming:
  //
  // Usage: GREM fastq vg length chunksize

  if (argc != 5) LOG(FATAL) << "too few or too many arguments.";

  CharString fqPath    = argv[1];
  CharString vgPath    = argv[2];
  unsigned int seedLen = std::stoul(argv[3]);
  unsigned int chkSize = std::stoul(argv[4]);

  LOG(INFO) << "Opening file '" << toCString(fqPath) << "'...";

  SeqFileIn readInFile;
  if (!open(readInFile, toCString(fqPath)))
  {
    LOG(FATAL) << "could not open the file '" << toCString(fqPath) << "'.";
  }

  try
  {
    LOG(INFO) << "Loading the vg graph from file '" << toCString(vgPath) << "'...";

    std::string graph_name = "graph-1";
    VarGraph vargraph(toCString(vgPath), graph_name);

    LOG(INFO) << "Loading the vg graph from file '" << toCString(vgPath) << "': Done.";

    long int found = 0;
    while (true)
    {
      ReadsChunk reads;

      {
#ifndef NDEBUG
        TIMED_SCOPE(loadChunkTimer, "load-chunk");
#endif
        readRecords(reads.ids, reads.seqs, reads.quals, readInFile, chkSize);
      }

      if (length(reads.ids) == 0) break;

      GraphTraverser< PathTraverser > gtraverser(vargraph);

      {
#ifndef NDEBUG
        TIMED_SCOPE(addStartsTimer, "add-starts");
#endif
        for (unsigned int i = 0; i < vargraph.nodes_size(); ++i)
        {
          const vg::Node &node = vargraph.node_at(i);
          for (unsigned int j = 0; j < node.sequence().length(); ++j)
          {
            vg::Position s_point;
            s_point.set_node_id(node.id());
            s_point.set_offset(j);
//          s_point.set_offset(0);

            gtraverser.add_start(s_point);
          }
        }
      }

      std::function< void(vg::Alignment &) > write = [&found](vg::Alignment &aln){
        ++found;
        if (found % 1000 == 0) LOG(DEBUG) << found << " seeds found so far.";
      };

      PathTraverser::Param params(reads, seedLen);
      gtraverser.traverse(params, write);
    }

    LOG(INFO) << "Total number of " << found << " seeds found.";
  }
  catch(std::ios::failure &e)
  {
    LOG(ERROR) << "failed to open the file '" << toCString(vgPath) << "'.";
    LOG(FATAL) << "Caught an ios_base::failure: " << e.what();
  }

  // Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return EXIT_SUCCESS;
}

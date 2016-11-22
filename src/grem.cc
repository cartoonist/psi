/*
 * =====================================================================================
 *
 * Filename: grem.cpp
 *
 * Created: Tue Nov 08, 2016  16:48
 * Last modified: Tue Nov 22, 2016  19:42
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
#include <easyloggingpp/src/easylogging++.h>

#include "vargraph.h"
#include "traverser.h"
#include "types.h"
#include "release.h"

using namespace seqan;
using namespace grem;

// TODO: Documentation.
// TODO: Memory footprint.
// TODO: Logging: replace all iostream calls with logging APIs.
// TODO: Logging messages.

INITIALIZE_EASYLOGGINGPP


int main(int argc, char *argv[])
{
  el::Configurations log_conf;
  log_conf.setToDefault();
  el::Loggers::reconfigureLogger("default", log_conf);

  // TODO: get command-line arguments. Now, assuming:
  //
  // Usage: GREM fastq vg length chunksize

  if (argc != 5)
  {
    std::cerr << "[" << PACKAGE << "] Error: "
              << "too few or too many arguments." << std::endl;
    exit(EXIT_FAILURE);
  }

  CharString fqPath    = argv[1];
  CharString vgPath    = argv[2];
  unsigned int seedLen = std::stoul(argv[3]);
  unsigned int chkSize = std::stoul(argv[4]);

#ifndef NDEBUG
  LOG(INFO) << "Opening file '" << toCString(fqPath) << "'...";
#endif

  SeqFileIn readInFile;
  if (!open(readInFile, toCString(fqPath)))
  {
    std::cerr << "[" << PACKAGE << "] Error: "
              << "could not open the file '" << toCString(fqPath)
              << "'." << std::endl;
    exit(EXIT_FAILURE);
  }

#ifndef NDEBUG
  LOG(INFO) << "Loading the vg graph from file '" << toCString(vgPath) << "'...";
#endif

  try
  {
    std::string graph_name = "graph-1";
    VarGraph vargraph(toCString(vgPath), graph_name);

    while (true)
    {
      ReadsChunk reads;

      readRecords(reads.ids, reads.seqs, reads.quals, readInFile, chkSize);
      if (length(reads.ids) == 0) break;

      GraphTraverser< PathTraverser > gtraverser(vargraph);
      for (unsigned int i = 0; i < vargraph.nodes_size(); ++i)
      {
        const vg::Node &node = vargraph.node_at(i);
        for (unsigned int j = 0; j < node.sequence().length(); ++j)
        {
          vg::Position s_point;
          s_point.set_node_id(node.id());
          s_point.set_offset(j);
//        s_point.set_offset(0);

          gtraverser.add_start(s_point);
        }
      }

      long int found = 0;
      std::function< void(vg::Alignment &) > write = [&found](vg::Alignment &aln){
#ifndef NDEBUG
        LOG(INFO) << ++found << " seeds found: "
                  << aln.name() << " @ " << aln.path().name();
#endif
      };

      PathTraverser::Param params(reads, seedLen);
      gtraverser.traverse(params, write);
    }
  }
  catch(std::ios::failure &e)
  {
    std::cerr << "[" << PACKAGE << "] Error: "
              << "failed to open the file '" << toCString(vgPath)
              << "'." << std::endl << "Caught an ios_base::failure: "
              << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}

/*
 * =====================================================================================
 *
 * Filename: grem.cpp
 *
 * Created: Tue Nov 08, 2016  16:48
 * Last modified: Mon Nov 21, 2016  01:29
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
// TODO: Logging.

INITIALIZE_EASYLOGGINGPP


int main(int argc, char *argv[])
{
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

  SeqFileIn readInFile;
  if (!open(readInFile, toCString(fqPath)))
  {
    std::cerr << "[" << PACKAGE << "] Error: "
              << "could not open the file '" << toCString(fqPath)
              << "'." << std::endl;
    exit(EXIT_FAILURE);
  }

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
        vg::Position s_point;
        s_point.set_node_id(vargraph.node_at(i).id());
        s_point.set_offset(0);
        gtraverser.add_start(s_point);
      }
      std::function< void(vg::Alignment &) > write = [](vg::Alignment &aln){ return; };
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

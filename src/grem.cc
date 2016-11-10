/*
 * =====================================================================================
 *
 * Filename: grem.cpp
 *
 * Created: Tue Nov 08, 2016  16:48
 * Last modified: Mon Nov 14, 2016  01:02
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

#include <seqan/seq_io.h>

#include "vargraph.h"
#include "types.h"
#include "release.h"

using namespace seqan;
using namespace grem;


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
  }
  catch(std::ios::failure &e)
  {
    std::cerr << "[" << PACKAGE << "] Error: "
              << "failed to open the file '" << toCString(vgPath)
              << "'." << std::endl << "Caught an ios_base::failure: "
              << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

  while (true)
  {
    ReadsChunk reads;

    readRecords(reads.ids, reads.seqs, reads.quals, readInFile, chkSize);
    if (length(reads.ids) == 0) break;

    for (unsigned int i = 0; i < length(reads.ids); ++i)
    {
      std::cout << reads.ids[i]   << '\t'
                << reads.seqs[i]  << '\t'
                << reads.quals[i] << std::endl;
    }
  }

  std::cout << seedLen << std::endl;

  return EXIT_SUCCESS;
}

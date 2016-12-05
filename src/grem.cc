/*
 * =====================================================================================
 *
 * Filename: grem.cpp
 *
 * Created: Tue Nov 08, 2016  16:48
 * Last modified: Mon Dec 05, 2016  01:54
 *
 * Description: grem main function.
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
#include <seqan/arg_parse.h>

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
// TODO: Logging messages.
// TODO: Source codes line limit: 88.
// TODO: performance logs' level should be DEBUG.

INITIALIZE_EASYLOGGINGPP

// Forwards
void                               setup_argparser(seqan::ArgumentParser & parser);
seqan::ArgumentParser::ParseResult parse_args(GremOptions & options, int argc, char *argv[]);


int main(int argc, char *argv[])
{
  START_EASYLOGGINGPP(argc, argv);

  el::Configurations log_conf;
  log_conf.setToDefault();
  el::Loggers::reconfigureLogger("default", log_conf);

  // Verify that the version of the library that we linked against is
  // compatible with the version of the headers we compiled against.
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  // Parse the command line.
  GremOptions options;
  auto res = parse_args(options, argc, argv);
  // If parsing was not successful then exit with code 1 if there were errors.
  // Otherwise, exit with code 0 (e.g. help was printed).
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

  CharString const & fqPath    = options.fq_path;
  CharString const & vgPath    = options.rf_path;
  unsigned int const & seedLen = options.seed_len;
  unsigned int const & chkSize = options.chunk_size;

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

    ReadsChunk reads;
    long int found = 0;
    GraphTraverser< PathTraverser > gtraverser(vargraph);
    gtraverser.add_all_loci();
    std::function< void(vg::Alignment &) > write = [&found](vg::Alignment &aln){
      ++found;
      if (found % 1000 == 0) LOG(DEBUG) << found << " seeds found so far.";
    };
    while (true)
    {
      {
#ifndef NDEBUG
        TIMED_SCOPE(loadChunkTimer, "load-chunk");
#endif
        readRecords(reads.ids, reads.seqs, reads.quals, readInFile, chkSize);
      }

      if (length(reads.ids) == 0) break;

      PathTraverser::Param params(reads, seedLen);
      gtraverser.traverse(params, write);

      clear(reads.ids);
      clear(reads.seqs);
      clear(reads.quals);
    }

    LOG(INFO) << "Total number of " << found << " seeds found.";
  }
  catch(std::ios::failure &e)
  {
    // FIXME: It may capture more general exceptions. This try/catch block is for
    //        exceptions thrown by VarGraph ctor. It should surround that line. In order
    //        to do so default ctor should be defined for VarGraph class.
    LOG(ERROR) << "failed to open the file '" << toCString(vgPath) << "'.";
    LOG(FATAL) << "Caught an ios_base::failure: " << e.what();
  }

  // Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return EXIT_SUCCESS;
}


  void
setup_argparser(seqan::ArgumentParser & parser)
{
  // positional arguments.
  std::string POSARG1 = "VG_FILE";

  // add usage line.
  addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fI" + POSARG1 + "\\fP\"");

  // vg file -- positional argument.
  seqan::ArgParseArgument vgfile_arg(seqan::ArgParseArgument::INPUT_FILE, POSARG1);
  setValidValues(vgfile_arg, "vg");
  addArgument(parser, vgfile_arg);

  // reads in FASTQ format -- **required** option.
  seqan::ArgParseOption fqfile_arg("f", "fastq", "Reads in FASTQ format.",
                                   seqan::ArgParseArgument::INPUT_FILE, "FASTQ_FILE");
  setValidValues(fqfile_arg, "fq fastq");
  addOption(parser, fqfile_arg);
  setRequired(parser, "f");

  // seed length -- **required** option.
  addOption(parser, seqan::ArgParseOption("l", "seed-length", "Seed length.",
                                          seqan::ArgParseArgument::INTEGER, "INT"));
  setRequired(parser, "l");

  // chunk size -- **required** option.
  addOption(parser, seqan::ArgParseOption("c", "chunk-size", "Reads chunk size.",
                                          seqan::ArgParseArgument::INTEGER, "INT"));
  setRequired(parser, "c");

  // verbosity options -- HANDLED BY EASYLOGGING++
  addOption(parser, seqan::ArgParseOption("v", "verbose",
                                          "Activates maximum verbosity."));
}


  seqan::ArgumentParser::ParseResult
parse_args(GremOptions & options, int argc, char *argv[])
{
  // setup ArgumentParser.
  seqan::ArgumentParser parser("grem");
  setup_argparser(parser);

  // Embedding program's meta data and build information.
  setShortDescription(parser, SHORT_DESC);
  setVersion(parser, VERSION);
  setDate(parser, __DATE__);
  addDescription(parser, LONG_DESC);

  // parse command line.
  auto res = seqan::parse(parser, argc, argv);

  // only extract options if the program will continue after parse_args()
  if (res != seqan::ArgumentParser::PARSE_OK) return res;

  getOptionValue(options.fq_path, parser, "fastq");
  getOptionValue(options.seed_len, parser, "seed-length");
  getOptionValue(options.chunk_size, parser, "chunk-size");
  getArgumentValue(options.rf_path, parser, 0);

  return seqan::ArgumentParser::PARSE_OK;
}

/*
 * =====================================================================================
 *
 * Filename: grem.cpp
 *
 * Created: Tue Nov 08, 2016  16:48
 * Last modified: Mon Dec 19, 2016  01:43
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

#include <easyloggingpp/src/easylogging++.h>

using namespace seqan;
using namespace grem;

// TODO: Documentation.
// TODO: Memory footprint.
// TODO: Logging messages.
// TODO: Source codes line limit: 88.
// TODO: performance logs' level should be DEBUG.
// TODO: handle 'N's.

INITIALIZE_EASYLOGGINGPP

// Forwards
void                               open_fastq(const CharString & fqpath, SeqFileIn & infile);
void                               load_graph(const CharString & vgpath, VarGraph & vargraph);
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

  CharString const & fqpath    = options.fq_path;
  CharString const & vgpath    = options.rf_path;
  unsigned int const & seedlen = options.seed_len;
  unsigned int const & chksize = options.chunk_size;

  SeqFileIn reads_infile;
  open_fastq(fqpath, reads_infile);
  VarGraph vargraph;
  load_graph(vgpath, vargraph);

  GraphTraverser< PathTraverser > gtraverser(vargraph);
  gtraverser.add_all_loci();

  long int found = 0;
  std::function< void(PathTraverser::Output &) > write = [&found]
    (PathTraverser::Output & seed_hit){
    ++found;
    if (found % SEEDHITS_REPORT_BUF == 0)
    {
      LOG(DEBUG) << found << " seeds found so far.";
    }
  };

  ReadsChunk reads;
  while (true)
  {
    {
      TIMED_SCOPE(loadChunkTimer, "load-chunk");
      readRecords(reads.ids, reads.seqs, reads.quals, reads_infile, chksize);
    }

    if (length(reads.ids) == 0) break;

    PathTraverser::Param params(reads, seedlen);
    gtraverser.traverse(params, write);

    clear(reads.ids);
    clear(reads.seqs);
    clear(reads.quals);
  }

  LOG(INFO) << "Total number of " << found << " seeds found.";

  // Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return EXIT_SUCCESS;
}


  inline void
open_fastq(const CharString & fqpath, SeqFileIn & infile)
{
  LOG(INFO) << "Opening file '" << toCString(fqpath) << "'...";

  if (!open(infile, toCString(fqpath)))
  {
    LOG(FATAL) << "could not open the file '" << toCString(fqpath) << "'.";
  }
}


  inline void
load_graph(const CharString & vgpath, VarGraph & vargraph)
{
  try
  {
    LOG(INFO) << "Loading the vg graph from file '" << toCString(vgpath) << "'...";

    vargraph.extend_from_file(toCString(vgpath));

    LOG(INFO) << "Loading the vg graph from file '" << toCString(vgpath) << "': Done.";
  }
  catch(std::ios::failure &e)
  {
    LOG(ERROR) << "failed to open the file '" << toCString(vgpath) << "'.";
    LOG(FATAL) << "Caught an ios_base::failure: " << e.what();
  }
}


  inline void
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


  inline seqan::ArgumentParser::ParseResult
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

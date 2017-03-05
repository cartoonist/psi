/*
 * =====================================================================================
 *
 * Filename: bfs.cc
 *
 * Created: Wed Jan 18, 2017  16:07
 * Last modified: Thu Mar 02, 2017  13:16
 *
 * Description: Test variation graph iterator class.
 *
 * Copyright (c) 2017, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari (cartoonist), ali.ghaffaari@mpi-inf.mpg.de
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#include <string>

#include <seqan/arg_parse.h>

#include "vargraph.h"
#include "vargraph_iter.h"
#include "base.h"
#include "release.h"

#include <easyloggingpp/src/easylogging++.h>

INITIALIZE_EASYLOGGINGPP

using namespace grem;

typedef struct
{
  std::string vgpath;
  unsigned long int step;
  grem::VarGraph::NodeID start;
} Options;


void
setup_argparser(seqan::ArgumentParser& parser)
{
  // positional arguments.
  std::string POSARG1 = "VG_FILE";

  // add usage line.
  addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fI" + POSARG1 + "\\fP\"");

  // graph file -- positional argument.
  seqan::ArgParseArgument vg_file_arg(seqan::ArgParseArgument::INPUT_FILE, POSARG1);
  setValidValues(vg_file_arg, "vg");
  addArgument(parser, vg_file_arg);

  // start node
  addOption(parser, seqan::ArgParseOption("s", "start-node", "Starting node id.",
                                          seqan::ArgParseArgument::INTEGER, "INT"));
  setMinValue(parser, "s", "1");
  setDefaultValue(parser, "s", 1);

  // start every
  addOption(parser, seqan::ArgParseOption("e", "start-every",
                                          "Add starting position at every this number.",
                                          seqan::ArgParseArgument::INTEGER, "INT"));
  setMinValue(parser, "e", "1");
  setRequired(parser, "e");

  // verbosity options -- HANDLED BY EASYLOGGING++
  addOption(parser, seqan::ArgParseOption("v", "verbose",
                                          "Activates maximum verbosity."));
}


seqan::ArgumentParser::ParseResult
parse_args(Options & options, int argc, char *argv[])
{
  // setup ArgumentParser.
  seqan::ArgumentParser parser("test_bfs");
  setup_argparser(parser);

  // Embedding program's meta data and build information.
  setShortDescription(parser, "Breadth-first Graph Traverse");
  setVersion(parser, VERSION);
  // TODO: Set date to the date of the commit at HEAD.
  //       Use `-D` option of gcc to define a macro equal to the date of the commit at
  //       HEAD which is captured by git command.
  setDate(parser, __DATE__);
  addDescription(parser,
                 "Testing variation graph iterator by implementing "
                 "BFS graph traversal algorithm.");

  // parse command line.
  auto res = seqan::parse(parser, argc, argv);

  // only extract options if the program will continue after parse_args()
  if (res != seqan::ArgumentParser::PARSE_OK) return res;

  getOptionValue(options.step, parser, "start-every");
  getOptionValue(options.start, parser, "start-node");
  getArgumentValue(options.vgpath, parser, 0);

  return seqan::ArgumentParser::PARSE_OK;
}


int main(int argc, char *argv[])
{
  START_EASYLOGGINGPP(argc, argv);

  Options options;
  auto res = parse_args(options, argc, argv);
  // If parsing was not successful then exit with code 1 if there were errors.
  // Otherwise, exit with code 0 (e.g. help was printed).
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

  VarGraph vargraph(options.vgpath);
  Iterator<VarGraph, BFS<>> itr(vargraph, options.start);

  unsigned long int prenode_remain = 0;
  unsigned long int remain_estimate = 0;
  grem::VarGraph::NodeID prenode_level = 0;
  std::string seq;
  while (!at_end(itr))
  {
    if (prenode_level != level(itr))
    {
      prenode_remain = remain_estimate;
      remain_estimate = 0;
      prenode_level = level(itr);
    }

    seq = vargraph.node_by(*itr).sequence();

    unsigned long int cursor = (options.step - prenode_remain) % options.step;
    while (cursor < seq.length())
    {
      std::cout << *itr << ", " << cursor << std::endl;
      cursor += options.step;
    }

    unsigned long int new_remain;
    if (options.step - prenode_remain > seq.length())
    {
      new_remain = prenode_remain + seq.length();
    }
    else
    {
      new_remain = (seq.length() - options.step + prenode_remain) % options.step;
    }

    if (remain_estimate < new_remain)
    {
      remain_estimate = new_remain;
    }

    ++itr;
  }

  return EXIT_SUCCESS;
}

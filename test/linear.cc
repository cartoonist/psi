/**
 *    @file  linear.cc
 *   @brief  Finding seed hits in a linear sequence.
 *
 *  This program does exactly same algorithm on linear sequence for comparison.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Nov 29, 2016  15:04
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdlib>

#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

#include "sequence.h"
#include "index.h"
#include "utils.h"
#include "options.h"
#include "stat.h"
#include "logger.h"
#include "release.h"

#define SEEDHITS_REPORT_BUF 1000
#define TRAVERSE_CHECKPOINT_LOCI_NO 1000000

using namespace std;
using namespace seqan;
using namespace grem;

// TODO: Localize Options (it is written for the main program).


typedef struct
{
  typedef seqan::Index< Dna5QStringSet<>, IndexWotd<> > TIndex;
  typedef typename Iterator < TIndex, TopDown<> >::Type TIndexIter;

  TIndexIter index_iter;
  unsigned int   boffset;
  unsigned int   ref_len;
} IterState;


void
setup_argparser(seqan::ArgumentParser& parser)
{
  // positional arguments.
  std::string POSARG1 = "REF_FILE";

  // add usage line.
  addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fI" + POSARG1 + "\\fP\"");

  // reference file -- positional argument.
  seqan::ArgParseArgument ref_file_arg(seqan::ArgParseArgument::INPUT_FILE, POSARG1);
  setValidValues(ref_file_arg, "fa fasta");
  addArgument(parser, ref_file_arg);

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

  // quiet -- no output to console
  addOption(parser, seqan::ArgParseOption("q", "quiet", "Quiet mode. No output will be printed to console."));

  // no colored output
  addOption(parser, seqan::ArgParseOption("C", "no-color", "Do not use a colored output."));

  // verbosity options
  addOption(parser, seqan::ArgParseOption("v", "verbose",
                                          "Activates maximum verbosity."));
}


seqan::ArgumentParser::ParseResult
parse_args(Options & options, int argc, char *argv[])
{
  // setup ArgumentParser.
  seqan::ArgumentParser parser("test_linear");
  setup_argparser(parser);

  // Embedding program's meta data and build information.
  setShortDescription(parser, "Find seed hits for a linear sequence.");
  setVersion(parser, VERSION);
  // TODO: Set date to the date of the commit at HEAD.
  //       Use `-D` option of gcc to define a macro equal to the date of the commit at
  //       HEAD which is captured by git command.
  setDate(parser, __DATE__);
  addDescription(parser,
                 "Instead of using graph, this simple program uses a linear reference "
                 "in order to find seed hits.");

  // parse command line.
  auto res = seqan::parse(parser, argc, argv);

  // only extract options if the program will continue after parse_args()
  if (res != seqan::ArgumentParser::PARSE_OK) return res;

  getOptionValue(options.fq_path, parser, "fastq");
  getOptionValue(options.seed_len, parser, "seed-length");
  getOptionValue(options.chunk_size, parser, "chunk-size");
  options.quiet = isSet( parser, "quiet" );
  options.nocolor = isSet( parser, "no-color" );
  options.verbose = isSet( parser, "verbose" );
  getArgumentValue(options.rf_path, parser, 0);

  return seqan::ArgumentParser::PARSE_OK;
}


bool go_down(IterState &its, seqan::Value< Dna5QString >::Type c)
{
  if ( c == 'N' ) return false;

  if (its.boffset == 0) {
    if (!seqan::goDown(its.index_iter, c)) return false;

    its.boffset = parentEdgeLength(its.index_iter) - 1;
  } else if (c == parentEdgeLabel(its.index_iter)[ parentEdgeLength(its.index_iter) - its.boffset ]) {
    --its.boffset;
  } else {
    return false;
  }

  ++its.ref_len;

  return true;
}


int main(int argc, char *argv[])
{
  // Parse the command line.
  Options options;
  auto res = parse_args(options, argc, argv);
  // If parsing was not successful then exit with code 1 if there were errors.
  // Otherwise, exit with code 0 (e.g. help was printed).
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

  options.nolog = false;
  options.nologfile = true;
  config_logger(options);
  auto log = get_logger("main");

  SeqFileIn refInFile;
  if (!open(refInFile, options.rf_path.c_str()))
  {
    std::string msg = "could not open the file '" + options.rf_path + "'.";
    log->error(msg);
    throw std::runtime_error(msg);
  }

  CharString ref_id;
  Dna5QString     ref_seq;

  {
    auto timer = Timer("load-ref");
    readRecord(ref_id, ref_seq, refInFile);
  }
  log->info("Reference loaded in {} us.", Timer::get_duration("load-ref").count());
  log->info("Reference sequence length: {}.", length(ref_seq));

  SeqFileIn readsInFile;
  if (!open(readsInFile, options.fq_path.c_str()))
  {
    std::string msg = "could not open the file '" + options.fq_path + "'.";
    log->error(msg);
    throw std::runtime_error(msg);
  }

  Records< Dna5QStringSet<> > reads;
  bool found;
  unsigned int nof_found = 0;
  log->info("Seed finding...");
  {
    auto timer = Timer("seed-finding");
    while (true)
    {
      {
        auto timer = Timer("load-reads");
        readRecords(reads, readsInFile, options.chunk_size);
      }
      log->info("Reads loaded in {} us.", Timer::get_duration("load-reads").count());
      if (length(reads.name) == 0)
      {
        log->info("All reads are processed.");
        break;
      }
      log->info( "Reading {} reads...", options.chunk_size );

      {
        auto timer = Timer("traverse");

        unsigned int i;
        IterState::TIndex reads_index( reads.str );
        for (unsigned int pos = 0; pos < length(ref_seq); ++pos)
        {
          IterState iter_state = {IterState::TIndexIter(reads_index), 0, 0};
          // Explicit is better than implicit.
          if (pos > length(ref_seq) - options.seed_len /*+ allowed_diffs*/) break;

          found = true;
          for (i = pos; iter_state.ref_len < options.seed_len; ++i)
          {
            if(!go_down(iter_state, ref_seq[i]))
            {
              found = false;
              break;
            }
          }
          if (found)
          {
            ++nof_found;
            if (nof_found % SEEDHITS_REPORT_BUF == 0)
              log->info("{} seed hits so far.", nof_found);
          }

          if (pos % TRAVERSE_CHECKPOINT_LOCI_NO == 0)
          {
            log->info("Traversing lap: {} us.", Timer::get_lap("traverse").count());
          }
        }
      }
      log->info("Traversed in {} us.", Timer::get_duration("traverse").count());
      clear(reads.name);
      clear(reads.str);
    }
  }
  log->info("Seed finding was done in {} us.", Timer::get_duration("seed-finding").count());

  return EXIT_SUCCESS;
}

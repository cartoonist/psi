/**
 *    @file  grem.cc
 *   @brief  grem main program.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.de>
 *
 *  @internal
 *       Created:  Tue Nov 08, 2016  16:48
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdlib>
#include <csignal>
#include <iostream>
#include <ios>
#include <fstream>
#include <sstream>
#include <string>
#include <functional>
#include <unordered_set>

#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

#include "vargraph.h"
#include "mapper.h"
#include "traverser.h"
#include "pathset.h"
#include "sequence.h"
#include "utils.h"
#include "options.h"
#include "stat.h"
#include "logger.h"
#include "release.h"

INITIALIZE_EASYLOGGINGPP

using namespace seqan;
using namespace grem;

// TODO: Documentation.
// TODO: Memory footprint.
// TODO: Logging messages.
// TODO: Source codes line limit: 88.
// TODO: performance logs' level should be DEBUG.
// TODO: handle 'N's.
// TODO: comments' first letter?
// TODO: fix code style: spacing.
// TODO: inconsistency: some public methods are interface functions, some are members.

// Forwards
  void
startup ( const Options & options );

template < typename TIndexSpec >
  void
find_seeds ( VarGraph & vargraph, SeqFileIn & reads_infile, unsigned int seed_len,
    unsigned int chunk_size, unsigned int start_every, unsigned int path_num,
    std::string paths_index_file, bool nomapping, TIndexSpec const /* Tag */ );

template< typename TMapper >
    void
  signal_handler( int signal );

  void
default_signal_handler( int signal ) { /* Do nothing */ }

  seqan::ArgumentParser::ParseResult
parse_args ( Options & options, int argc, char *argv[] );

  inline void
get_option_values ( Options & options, seqan::ArgumentParser & parser );

  void
setup_argparser ( seqan::ArgumentParser & parser );

  void
config_logger ( const Options & options );

  void
config_logger ( const char * logger_name, const Options & options );


int main(int argc, char *argv[])
{
  START_EASYLOGGINGPP(argc, argv);

  // Parse the command line.
  Options options;
  auto res = parse_args(options, argc, argv);
  // If parsing was not successful then exit with code 1 if there were errors.
  // Otherwise, exit with code 0 (e.g. help was printed).
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

  // Verify that the version of the library that we linked against is
  // compatible with the version of the headers we compiled against.
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  std::signal( SIGUSR1, default_signal_handler );

  /* Configure loggers */
  config_logger(options);
  config_logger("default", options);
  config_logger("performance", options);

  startup ( options );

  // Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return EXIT_SUCCESS;
}

  void
startup ( const Options & options )
{
  LOG(INFO) << "Parameters:";
  LOG(INFO) << "- Seed length: " << options.seed_len;
  LOG(INFO) << "- Reads index type: " << index_to_str(options.index);
  LOG(INFO) << "- Starting points interval: " << options.start_every;
  LOG(INFO) << "- Reads chunk size: " << options.chunk_size;
  LOG(INFO) << "- Number of paths: " << options.path_num;
  LOG(INFO) << "- Paths index file: " << options.paths_index_file;

  SeqFileIn reads_infile;
  LOG(INFO) << "Opening file '" << options.fq_path << "'...";
  if (!open(reads_infile, options.fq_path.c_str()))
  {
    LOG(FATAL) << "could not open the file '" << options.fq_path << "'.";
  }

  std::ifstream ifs( options.rf_path, std::ifstream::in | std::ifstream::binary );
  LOG(INFO) << "Loading the graph from file '" << options.rf_path << "'...";
  if( !ifs )
  {
    LOG(FATAL) << "could not open the file '" << options.rf_path << "'.";
  }

  bool xg_format = true;
  if ( ends_with( options.rf_path, ".vg" ) ) {
    xg_format = false;
  }

  VarGraph vargraph( ifs, xg_format );

  switch ( options.index ) {
    case IndexType::Wotd: find_seeds ( vargraph,
                              reads_infile,
                              options.seed_len,
                              options.chunk_size,
                              options.start_every,
                              options.path_num,
                              options.paths_index_file,
                              options.nomapping,
                              UsingIndexWotd() );
                          break;
    case IndexType::Esa: find_seeds ( vargraph,
                             reads_infile,
                             options.seed_len,
                             options.chunk_size,
                             options.start_every,
                             options.path_num,
                             options.paths_index_file,
                             options.nomapping,
                             UsingIndexEsa() );
                         break;
    default: throw std::runtime_error("Index not implemented.");
             break;
  }

}


template<typename TIndexSpec >
  void
find_seeds ( VarGraph & vargraph, SeqFileIn & reads_infile, unsigned int seed_len,
    unsigned int chunk_size, unsigned int start_every, unsigned int path_num,
    std::string paths_index_file, bool nomapping, TIndexSpec const /* Tag */ )
{
  typedef typename Traverser< TIndexSpec, BFS, ExactMatching >::Type TTraverser;
  Mapper< TTraverser > mapper( &vargraph, seed_len );

  typedef decltype( mapper ) TMapper;
  std::signal( SIGUSR1, signal_handler< TMapper > );

  // :TODO:Mon Mar 06 13:00:\@cartoonist: IndexEsa<> -> IndexFM<>
  PathSet< seqan::IndexEsa<> > paths;
  if ( path_num == 0 ) {
    LOG(INFO) << "Specified number of path is 0. Skipping paths indexing...";
  }
  else if ( !paths.load( paths_index_file, &vargraph, path_num ) ) {
    LOG(INFO) << "No valid paths index found. Picking paths...";
    mapper.pick_paths( paths, path_num );
    {
      auto timer = Timer( "path-indexing" );
      LOG(INFO) << "Indexing the paths...";
      paths.create_index();
      LOG(INFO) << "Saving paths index...";
      if ( !paths.save( paths_index_file.c_str() ) ) {
        LOG(WARNING) << "Paths index file is not specified or is not writable.";
        LOG(INFO) << "Skipping...";
      }
    }
  }
  else {
    LOG(INFO) << "Paths index found. Loaded.";
  }

  if ( nomapping ) {
    LOG(INFO) << "Skipping mapping as requested...";
    return;
  }

  mapper.add_all_loci( paths, seed_len, start_every );

   // :TODO:Mon May 08 12:02:\@cartoonist: read id reported in the seed is relative to the chunk.
  long int found = 0;
  std::unordered_set< Dna5QStringSetPosition > covered_reads;
  std::function< void(typename TTraverser::output_type const &) > write =
    [&found, &covered_reads] (typename TTraverser::output_type const & seed_hit){
    ++found;
    covered_reads.insert(seqan::beginPositionV(seed_hit));
  };

  Dna5QRecords reads_chunk;

  {
    auto timer = Timer( "seed-finding" );
    while (true)
    {
      {
        auto timer = Timer( "load-chunk" );
        readRecords(reads_chunk, reads_infile, chunk_size);
      }

      if (length(reads_chunk.id) == 0) break;

      mapper.set_reads( std::move( reads_chunk ) );
      mapper.seeds_on_paths( paths, write );
      LOG(INFO) << "Total number of seeds found on paths: " << found;
      mapper.traverse ( write );
      clear( reads_chunk.id );
      clear( reads_chunk.str );
    }
  }

  LOG(INFO) << "Total number of seeds found: " << found;
  LOG(INFO) << "Total number of reads covered: " << covered_reads.size();
  LOG(INFO) << "Total number of starting points: " << mapper.get_starting_loci().size();
#ifndef NDEBUG
  LOG(INFO) << "Total number of 'godown' operations: "
            << TTraverser::stats_type::get_total_nof_godowns();
#endif
}


template< typename TMapper >
  void
signal_handler( int signal )
{
  std::cout << std::endl << "Report requested by SIGUSR1" << std::endl
            << "---------------------------" << std::endl;
  std::cout << "Elapsed time in traversal phase: "
            << Stat< TMapper >::Type::get_lap( "traverse" ).count() << " us"
            << std::endl;
  std::cout << "Current node: ("
            << Stat< TMapper >::Type::get_current_locus().node_id() << ", "
            << Stat< TMapper >::Type::get_current_locus().offset() << ")" << std::endl;
  auto idx = Stat< TMapper >::Type::get_current_locus_idx();
  auto total = Stat< TMapper >::Type::get_total_nof_loci();
  unsigned int wlen = std::to_string( total ).length();
  std::cout << "Progress: " << std::setw(wlen) << idx << " / " << std::setw(wlen)
            << total << " [%" << std::setw(3) << idx * 100 / total << "]" << std::endl;
}


  inline void
setup_argparser(seqan::ArgumentParser & parser)
{
  // positional arguments.
  std::string POSARG1 = "VG_FILE";

  // add usage line.
  addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fI" + POSARG1 + "\\fP\"");

  // graph file -- positional argument.
  seqan::ArgParseArgument vgfile_arg(seqan::ArgParseArgument::INPUT_FILE, POSARG1);
  setValidValues(vgfile_arg, "vg xg");
  addArgument(parser, vgfile_arg);

  // reads in FASTQ format -- **required** option.
  seqan::ArgParseOption fqfile_arg("f", "fastq", "Reads in FASTQ format.",
                                   seqan::ArgParseArgument::INPUT_FILE, "FASTQ_FILE");
  setValidValues(fqfile_arg, "fq fastq");
  addOption(parser, fqfile_arg);
  setRequired(parser, "f");

  // paths index file
  seqan::ArgParseOption paths_index_arg ( "I", "paths-index", "Paths index file.",
      seqan::ArgParseArgument::STRING, "PATHS_INDEX_FILE" );
  addOption ( parser, paths_index_arg );
  setDefaultValue ( parser, "I", "/tmp/GREM.XXXXXX/paths_index" );

  // seed length -- **required** option.
  addOption(parser, seqan::ArgParseOption("l", "seed-length", "Seed length.",
                                          seqan::ArgParseArgument::INTEGER, "INT"));
  setRequired(parser, "l");

  // chunk size -- **required** option.
  addOption(parser, seqan::ArgParseOption("c", "chunk-size", "Reads chunk size.",
                                          seqan::ArgParseArgument::INTEGER, "INT"));
  setRequired(parser, "c");

  // starting points
  addOption(parser, seqan::ArgParseOption("e", "start-every", "Start from every given "
                                          "number of loci in all nodes. If it is set to"
                                          " 1, it means start from all positions in a "
                                          "node.", seqan::ArgParseArgument::INTEGER,
                                          "INT"));
  setDefaultValue(parser, "e", 1);

  // number of paths
  addOption ( parser, seqan::ArgParseOption ( "n", "path-num", "Number of paths from "
        "the variation graph in hybrid approach.", seqan::ArgParseArgument::INTEGER,
        "INT" ) );
  setDefaultValue(parser, "n", 0);

  // index
  addOption(parser, seqan::ArgParseOption("i", "index", "Index type for indexing reads.",
                                          seqan::ArgParseArgument::STRING, "INDEX"));
  setValidValues(parser, "i", "SA ESA WOTD DFI QGRAM FM");
  setDefaultValue(parser, "i", "WOTD");

  // no map -- just do the paths indexing phase
  addOption(parser, seqan::ArgParseOption("x", "only-index", "Only build paths index and skip mapping."));

  // log file
  seqan::ArgParseOption logfile_arg("L", "log-file",
                                    "Sets default log file for existing and future loggers.",
                                    seqan::ArgParseArgument::OUTPUT_FILE, "LOG_FILE");
  addOption(parser, logfile_arg);
  setDefaultValue(parser, "L", "logs/grem.log");

  // no log to file
  addOption(parser, seqan::ArgParseOption("Q", "no-log-file", "Disable writing logs to file (overrides \\fB-L\\fP)."));

  // quiet -- no output to console
  addOption(parser, seqan::ArgParseOption("q", "quiet", "Quiet mode. No output will be printed to console."));

  // no colored output
  addOption(parser, seqan::ArgParseOption("C", "no-color", "Do not use a colored output."));

  // disable logging
  addOption(parser, seqan::ArgParseOption("D", "disable-log", "Disable logging completely."));

  // verbosity options -- HANDLED BY EASYLOGGING++
  addOption(parser, seqan::ArgParseOption("v", "verbose",
                                          "Activates maximum verbosity."));
}


  inline void
get_option_values ( Options & options, seqan::ArgumentParser & parser )
{
  std::string indexname;

  getOptionValue(options.fq_path, parser, "fastq");
  getOptionValue(options.seed_len, parser, "seed-length");
  getOptionValue(options.chunk_size, parser, "chunk-size");
  getOptionValue(options.start_every, parser, "start-every");
  getOptionValue(options.path_num, parser, "path-num");
  getOptionValue(options.paths_index_file, parser, "paths-index");
  getOptionValue(indexname, parser, "index");
  options.nomapping = isSet(parser, "only-index");
  getOptionValue(options.log_path, parser, "log-file");
  options.nologfile = isSet(parser, "no-log-file");
  options.quiet = isSet(parser, "quiet");
  options.nocolor = isSet(parser, "no-color");
  options.nolog = isSet(parser, "disable-log");
  getArgumentValue(options.rf_path, parser, 0);

  options.index = index_from_str(indexname);
}


  inline seqan::ArgumentParser::ParseResult
parse_args(Options & options, int argc, char *argv[])
{
  // setup ArgumentParser.
  seqan::ArgumentParser parser( PACKAGE );
  setup_argparser(parser);

  // Embedding program's meta data and build information.
  setShortDescription(parser, SHORT_DESC);
  setVersion(parser, VERSION);
  // :TODO:Thu Apr 06 02:06:\@cartoonist: date should be captured from git.
  setDate(parser, __DATE__);
  addDescription(parser, LONG_DESC);

  std::ostringstream hold_buf_stdout;
  std::ostringstream hold_buf_stderr;
  // parse command line.
  auto res = seqan::parse(parser, argc, argv, hold_buf_stdout, hold_buf_stderr);
  // print the banner in help or version messages.
  if (res == seqan::ArgumentParser::PARSE_HELP ||
      res == seqan::ArgumentParser::PARSE_VERSION)
  {
    std::cout << BANNER << std::endl;
  }
  // print the buffer.
  std::cout << hold_buf_stderr.str();
  std::cout << hold_buf_stdout.str();

  // only extract options if the program will continue after parse_args()
  if (res != seqan::ArgumentParser::PARSE_OK) return res;

  get_option_values ( options, parser );

  return seqan::ArgumentParser::PARSE_OK;
}

  inline void
config_logger ( const Options & options )
{
  // Configure color output.
  if (!options.nocolor) el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
  // Configure time format.
  el::Loggers::addFlag(el::LoggingFlag::FixedTimeFormat);
}

  inline void
config_logger(const char * logger_name, const Options & options)
{
  el::Configurations conf;
  conf.setToDefault();
  // Configure log file.
  if (!options.nologfile)
  {
    conf.setGlobally(el::ConfigurationType::Filename, options.log_path);
  }
  // Enabling quiet mode.
  if (options.quiet)
  {
    conf.setGlobally(el::ConfigurationType::ToStandardOutput , std::string("false"));
  }
  // Enable or disable the performance logger.
  if (options.nolog)
  {
    conf.setGlobally(el::ConfigurationType::Enabled, std::string("false"));
  }
  // Reconfigure the logger.
  el::Loggers::reconfigureLogger(logger_name, conf);
}

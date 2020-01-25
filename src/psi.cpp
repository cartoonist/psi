/**
 *    @file  psi.cpp
 *   @brief  PSI command-line interface.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
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

#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <gum/io_utils.hpp>

#include "graph.hpp"
#include "seed_finder.hpp"
#include "sequence.hpp"
#include "seed.hpp"
#include "utils.hpp"
#include "options.hpp"
#include "stat.hpp"
#include "logger.hpp"
#include "release.hpp"

using namespace klibpp;
using namespace seqan;
using namespace psi;

// TODO: Documentation.
// TODO: Memory footprint.
// TODO: Logging messages.
// TODO: Source codes line limit: 88.
// TODO: performance logs' level should be DEBUG.
// TODO: handle 'N's.
// TODO: comments' first letter?
// TODO: fix code style: spacing.
// TODO: inconsistency: some public methods are interface functions, some are members.
// TODO: Add value_t< T > typedef as typename seqan::Value< T >::Type


template< typename TSeedFinder >
  void
signal_handler( int signal )
{
  std::cout << std::endl << "Report requested by SIGUSR1" << std::endl
            << "---------------------------" << std::endl;
  std::cout << "Elapsed time in traversal phase: "
            << Stat< TSeedFinder >::Type::get_lap_str( "seeds-off-paths" ) << std::endl;
  auto pos = Stat< TSeedFinder >::Type::get_lastproc_locus().load();
  std::cout << "Current node: (" << pos.node_id << ", " << pos.offset << ")"
            << std::endl;
  auto idx = Stat< TSeedFinder >::Type::get_lastdone_locus_idx().load();
  auto total = Stat< TSeedFinder >::Type::get_total_nof_loci();
  unsigned int wlen = std::to_string( total ).length();
  std::cout << "Progress: " << std::setw(wlen) << idx << " / " << std::setw(wlen)
            << total << " [%" << std::setw(3) << idx * 100 / total << "]" << std::endl;
}


  void
default_signal_handler( int signal ) { /* Do nothing */ }

template< typename TSeedFinder, typename TSet >
    void
  report( TSeedFinder& finder, TSet& covered_reads, unsigned long long int found )
  {
    /* Get the main logger. */
    auto log = get_logger( "main" );
    log->info( "Total number of starting loci: {}", finder.get_starting_loci().size() );
    log->info( "Total number of seeds found: {}", found );
    log->info( "Total number of reads covered: {}", covered_reads.size() );
    log->info( "Total number of 'godown' operations: {}",
      TSeedFinder::traverser_type::stats_type::get_total_nof_godowns() );

    log->info( "All Timers" );
    log->info( "----------" );
    for ( const auto& timer : Timer<>::get_timers() ) {
      log->info( "{}: {}", timer.first, Timer<>::get_duration_str( timer.first ) );
    }
  }


template< class TGraph, typename TReadsIndexSpec >
    void
  find_seeds( TGraph& graph, SeqStreamIn& reads_iss, seqan::File<>& output_file,
              Options const& params, TReadsIndexSpec const )
  {
    /* typedefs */
    typedef Dna5QStringSet<> TReadsStringSet;
    typedef seqan::Index< TReadsStringSet, TReadsIndexSpec > TReadsIndex;
    typedef typename Traverser< TGraph, TReadsIndex, BFS, ExactMatching >::Type TTraverser;
    typedef SeedFinder< TTraverser > TSeedFinder;

    /* Get the main logger. */
    auto log = get_logger( "main" );
    /* Install seed finder singal handler for getting progress report. */
    std::signal( SIGUSR1, signal_handler< TSeedFinder > );

    /* The seed finder for the input graph. */
    TSeedFinder finder( graph, params.seed_len );
    /* Prepare (load or create) genome-wide paths. */
    log->info( "Looking for an existing path index..." );
    /* Load the genome-wide path index for the graph if available. */
    if ( finder.load_path_index( params.pindex_path,
                                 params.context,
                                 params.step_size ) ) {
      log->info( "The path index has been found and loaded." );
      return;
    }
    /* No genome-wide path index requested. */
    if ( params.path_num == 0 ) {
      log->info( "No path has been specified. Skipping path indexing..." );
    }
    else {
      log->info( "No valid path index found. Creating the path index..." );
      log->info( "Selecting {} different path(s) in the graph...", params.path_num );
      finder.create_path_index( params.path_num, params.context,
                                params.patched, params.step_size,
                                [&log]( std::string const& msg ) {
                                  log->info( msg );
                                },
                                [&log]( std::string const& msg ) {
                                  log->warn( msg );
                                } );
      log->info( "Picked paths in {}.", Timer<>::get_duration_str( "pick-paths" ) );
      log->info( "Indexed paths in {}.", Timer<>::get_duration_str( "index-paths" ) );
      log->info( "Selected starting loci in {}.",
                 Timer<>::get_duration_str( "add-starts" ) );
      log->info( "Saving path index..." );
      /* Serialize the indexed paths. */
      if ( params.pindex_path.empty() ) {
        log->warn( "No path index file is specified. Skipping..." );
      } else if ( !finder.serialize_path_index( params.pindex_path, params.step_size ) ) {
        log->warn( "Specified path index file is not writable. Skipping..." );
      } else {
        log->info( "Saved path index in {}.", Timer<>::get_duration_str( "save-paths" ) );
      }
    }
    log->info( "Number of uncovered loci (in {} nodes of total {}): {}",
        finder.get_nof_uniq_nodes(), finder.get_graph_ptr()->get_node_count(),
        finder.get_starting_loci().size() );

    if ( params.indexonly ) {
      log->info( "Skipping seed finding as requested..." );
      return;
    }

    unsigned long long int found = 0;
    unsigned long long int total_found = 0;
    std::unordered_set< Records< TReadsStringSet >::TPosition > covered_reads;
    std::function< void(typename TTraverser::output_type const &) > write_callback =
      [&found, &output_file, &covered_reads]
      (typename TTraverser::output_type const & seed_hit) {
      ++found;
      write( output_file, &seed_hit.node_id, 1 );
      write( output_file, &seed_hit.node_offset, 1 );
      write( output_file, &seed_hit.read_id, 1 );
      write( output_file, &seed_hit.read_offset, 1 );
      covered_reads.insert(seed_hit.read_id);
    };

    /* Found seeds in chunks. */
    {
      auto chunk = finder.create_readrecord();
      log->info( "Finding seeds..." );
      auto timer = Timer<>( "seed-finding" );
      while ( true ) {
        log->info( "Loading a read chunk..." );
        {
          auto timer = Timer<>( "load-chunk" );
          /* Load a chunk from reads set. */
          if ( !readRecords( chunk, reads_iss, params.chunk_size ) ) break;
        }
        log->info( "Fetched {} reads in {}.", length( chunk ),
            Timer<>::get_duration_str( "load-chunk" ) );
        /* Give the current chunk to the finder. */
        finder.set_reads( chunk, params.distance );
        log->info( "Seeding done in {}.", Timer<>::get_duration_str( "seeding" ) );
        log->info( "Finding seeds on paths..." );
        /* Find seeds on genome-wide paths. */
        finder.seeds_on_paths( write_callback );
        log->info( "Found seeds on paths in {}.",
            Timer<>::get_duration_str( "seeds-on-paths" ) );
        log->info( "Total number of seeds found on paths: {}", found );
        total_found += found;
        found = 0;
        log->info( "Finding seeds off paths..." );
        /* Find seeds on the graph by traversing starting loci. */
        finder.seeds_off_paths( write_callback );
        log->info( "Found seeds off paths in {}.", Timer<>::get_duration_str( "seeds-off-paths" ) );
        log->info( "Total number of seeds found off paths: {}", found );
        total_found += found;
        found = 0;
      }
    }
    log->info( "Found seed in {}.", Timer<>::get_duration_str( "seed-finding" ) );
    report( finder, covered_reads, total_found );
  }


  void
startup( const Options & options )
{
  auto log = get_logger( "main" );
  log->info( "Parameters:" );
  log->info( "- Seed length: {}", options.seed_len );
  log->info( "- Seed distance: {}", options.distance );
  log->info( "- Number of paths: {}", options.path_num );
  log->info( "- Context size (used in patching): {}", options.context );
  log->info( "- Patched: {}", ( options.patched ? "yes" : "no" ) );
  log->info( "- Path index file: '{}'", options.pindex_path );
  log->info( "- Reads chunk size: {}", options.chunk_size );
  log->info( "- Reads index type: {}", index_to_str(options.index) );
  log->info( "- Step size: {}", options.step_size );
  log->info( "- Temporary directory: '{}'", get_tmpdir() );
  log->info( "- Output file: '{}'", options.output_path );

  log->info( "Loading input graph from file '{}'...", options.rf_path );
  gum::SeqGraph< gum::Succinct > graph;
  gum::util::load( graph, options.rf_path );

  log->info( "Opening reads file '{}'...", options.fq_path );
  SeqStreamIn reads_iss( options.fq_path.c_str() );
  if ( !reads_iss ) {
    std::string msg = "could not open file '" + options.fq_path + "'!";
    log->error( msg );
    throw std::runtime_error( msg );
  }

  seqan::File<> output_file;
  auto mode = seqan::OPEN_CREATE | seqan::OPEN_WRONLY;
  if ( !open( output_file, options.output_path.c_str(), mode ) ) {
    std::string msg = "could not open file '" + options.output_path + "'!";
    log->error( msg );
    throw std::runtime_error( msg );
  }

  switch ( options.index ) {
    case IndexType::Wotd: find_seeds( graph,
                              reads_iss,
                              output_file,
                              options,
                              UsingIndexWotd() );
                          break;
    case IndexType::Esa: find_seeds( graph,
                             reads_iss,
                             output_file,
                             options,
                             UsingIndexEsa() );
                         break;
    default: throw std::runtime_error("Index not implemented.");
             break;
  }
}


  inline void
setup_argparser( seqan::ArgumentParser& parser )
{
  // positional arguments.
  std::string POSARG1 = "VG_FILE";
  // add usage line.
  addUsageLine( parser, "[\\fIOPTIONS\\fP] \"\\fI" + POSARG1 + "\\fP\"" );

  // graph file -- positional argument.
  seqan::ArgParseArgument vgfile_arg( seqan::ArgParseArgument::INPUT_FILE, POSARG1 );
  setValidValues( vgfile_arg, "vg xg" );
  addArgument( parser, vgfile_arg );

  // Options
  // reads in FASTQ format -- **required** option.
  addOption( parser,
      seqan::ArgParseOption( "f", "fastq",
        "Reads in FASTQ format.",
        seqan::ArgParseArgument::INPUT_FILE, "FASTQ_FILE" ) );
  setValidValues( parser, "f", "fq fastq" );
  setRequired( parser, "f" );
  // output file
  // :TODO:Sat Oct 21 00:06:\@cartoonist: output should be alignment in the GAM format.
  addOption( parser,
      seqan::ArgParseOption( "o", "output",
        "Output file.",
        seqan::ArgParseArgument::OUTPUT_FILE, "OUTPUT_FILE" ) );
  setDefaultValue( parser, "o", "out.gam" );
  // path index file
  addOption( parser,
      seqan::ArgParseOption( "I", "path-index",
        "Path index file.",
        seqan::ArgParseArgument::STRING, "PATH_INDEX_FILE" ) );
  // seed length -- **required** option.
  addOption( parser,
      seqan::ArgParseOption( "l", "seed-length",
        "Seed length.",
        seqan::ArgParseArgument::INTEGER, "INT" ) );
  setRequired( parser, "l" );
  // chunk size -- **required** option.
  addOption( parser,
      seqan::ArgParseOption( "c", "chunk-size",
        "Reads chunk size. Set it to 0 to consider all reads as one chunk (default).",
        seqan::ArgParseArgument::INTEGER, "INT" ) );
  setDefaultValue( parser, "c", 0 );
  // step size
  addOption( parser,
      seqan::ArgParseOption( "e", "step-size",
        "Minimum approximate distance allowed between two consecutive loci.",
        seqan::ArgParseArgument::INTEGER, "INT" ) );
  setDefaultValue( parser, "e", 1 );
  // seed distance
  addOption( parser,
      seqan::ArgParseOption( "d", "distance",
        "Distance between seeds",
        seqan::ArgParseArgument::INTEGER, "INT" ) );
  setDefaultValue( parser, "d", 0 );  /* Default value is seed length. */
  // number of paths
  addOption( parser,
      seqan::ArgParseOption( "n", "path-num",
        "Number of paths from the graph included in the path index.",
        seqan::ArgParseArgument::INTEGER, "INT" ) );
  setDefaultValue( parser, "n", 0 );
  // whether use patched paths or full genome-wide paths
  addOption( parser,
      seqan::ArgParseOption( "P", "no-patched",
        "Use full genome-wide paths." ) );
  // context in patching the paths
  addOption( parser,
      seqan::ArgParseOption( "t", "context",
        "Context length in patching.",
        seqan::ArgParseArgument::INTEGER, "INT" ) );
  setDefaultValue( parser, "t", 0 );
  // index
  addOption( parser,
      seqan::ArgParseOption( "i", "index",
        "Index type for indexing reads.",
        seqan::ArgParseArgument::STRING, "INDEX" ) );
  setValidValues( parser, "i", "SA ESA WOTD DFI QGRAM FM" );
  setDefaultValue( parser, "i", "WOTD" );
  // index only
  addOption( parser,
      seqan::ArgParseOption( "x", "index-only",
        "Only build path index and skip seed finding." ) );
  // log file
  addOption( parser,
      seqan::ArgParseOption( "L", "log-file",
        "Sets default log file for existing and future loggers.",
        seqan::ArgParseArgument::OUTPUT_FILE, "LOG_FILE" ) );
  setDefaultValue( parser, "L", "psi.log" );
  // no log to file
  addOption( parser,
      seqan::ArgParseOption( "Q", "no-log-file",
        "Disable writing logs to file (overrides \\fB-L\\fP)." ) );
  // quiet -- no output to console
  addOption( parser,
      seqan::ArgParseOption( "q", "quiet",
        "Quiet mode. No output will be printed to console." ) );
  // no colored output
  addOption( parser,
      seqan::ArgParseOption( "C", "no-color",
        "Do not use a colored output." ) );
  // disable logging
  addOption( parser,
      seqan::ArgParseOption( "D", "disable-log",
        "Disable logging completely." ) );
  // verbosity option
  addOption( parser,
      seqan::ArgParseOption( "v", "verbose",
        "Activates maximum verbosity." ) );
}


  inline void
get_option_values( Options & options, seqan::ArgumentParser & parser )
{
  std::string indexname;

  getOptionValue( options.fq_path, parser, "fastq" );
  getOptionValue( options.output_path, parser, "output" );
  getOptionValue( options.seed_len, parser, "seed-length" );
  getOptionValue( options.chunk_size, parser, "chunk-size" );
  getOptionValue( options.step_size, parser, "step-size" );
  getOptionValue( options.distance, parser, "distance" );
  getOptionValue( options.path_num, parser, "path-num" );
  getOptionValue( options.context, parser, "context" );
  options.patched = !isSet( parser, "no-patched" );
  getOptionValue( options.pindex_path, parser, "path-index" );
  getOptionValue( indexname, parser, "index" );
  options.indexonly = isSet( parser, "index-only" );
  getOptionValue( options.log_path, parser, "log-file" );
  options.nologfile = isSet( parser, "no-log-file" );
  options.quiet = isSet( parser, "quiet" );
  options.nocolor = isSet( parser, "no-color" );
  options.nolog = isSet( parser, "disable-log" );
  options.verbose = isSet( parser, "verbose" );
  getArgumentValue( options.rf_path, parser, 0 );

  options.index = index_from_str( indexname );
  if ( options.distance == 0 ) options.distance = options.seed_len;
}


  inline seqan::ArgumentParser::ParseResult
parse_args( Options& options, int argc, char* argv[] )
{
  // setup ArgumentParser.
  seqan::ArgumentParser parser( PACKAGE );
  setup_argparser( parser );

  // Embedding program's meta data and build information.
  setShortDescription( parser, SHORT_DESC );
  setVersion( parser, GIT_VERSION );
  setDate( parser, UPDATE_DATE );
  addDescription( parser, LONG_DESC );

  std::ostringstream hold_buf_stdout;
  std::ostringstream hold_buf_stderr;
  // parse command line.
  auto res = seqan::parse( parser, argc, argv, hold_buf_stdout, hold_buf_stderr );
  // print the banner in help or version messages.
  if ( res == seqan::ArgumentParser::PARSE_HELP ||
      res == seqan::ArgumentParser::PARSE_VERSION ) {
    std::cout << BANNER << std::endl;
  }
  // print the buffer.
  std::cout << hold_buf_stderr.str();
  std::cout << hold_buf_stdout.str();

  // only extract options if the program will continue after parse_args()
  if ( res != seqan::ArgumentParser::PARSE_OK ) return res;

  get_option_values( options, parser );

  return seqan::ArgumentParser::PARSE_OK;
}


  int
main( int argc, char *argv[] )
{
  /* Parse the command line. */
  Options options;
  auto res = parse_args( options, argc, argv );
  /* If parsing was not successful then exit with code 1 if there were errors.
   * Otherwise, exit with code 0 (e.g. help was printed).
   */
  if ( res != seqan::ArgumentParser::PARSE_OK ) {
    return res == seqan::ArgumentParser::PARSE_ERROR;
  }

  /* Verify that the version of the library that we linked against is
   * compatible with the version of the headers we compiled against.
   */
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  /* Install signal handler -- default handler */
  std::signal( SIGUSR1, default_signal_handler );
  /* Configure loggers */
  config_logger( options );
  /* Start seed finding... */
  startup( options );

  /* Close all loggers. */
  drop_all_loggers();
  /* Delete all global objects allocated by libprotobuf. */
  google::protobuf::ShutdownProtobufLibrary();

  return EXIT_SUCCESS;
}

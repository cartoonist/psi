/**
 *    @file  psikt.cpp
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
#include <psi/graph.hpp>
#include <psi/seed_finder.hpp>
#include <psi/sequence.hpp>
#include <psi/seed.hpp>
#include <psi/utils.hpp>
#include <psi/stats.hpp>
#include <psi/release.hpp>

#include "options.hpp"
#include "logger.hpp"

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

template< typename TSeedFinder, typename TSet >
    void
  report( TSeedFinder& finder, TSet& covered_reads, unsigned long long int found )
  {
    typedef typename TSeedFinder::stats_type::timer_type timer_type;

    /* Get the main logger. */
    auto log = get_logger( "main" );
    log->info( "Total number of starting loci: {}", finder.get_starting_loci().size() );
    log->info( "Total number of seeds found: {}", found );
    log->info( "-> of which found off paths: {}",
      TSeedFinder::traverser_type::stats_type::get_total_seeds_off_paths() );
    log->info( "Total number of reads covered: {}", covered_reads.size() );
    log->info( "Total number of 'godown' operations: {}",
      TSeedFinder::traverser_type::stats_type::get_total_nof_godowns() );

    log->info( "All Timers" );
    log->info( "----------" );
    for ( const auto& timer : timer_type::get_timers() ) {
      log->info( "{}: {}", timer.first, timer.second.str() );
    }
  }


template< class TGraph, typename TReadsIndexSpec >
    void
  find_seeds( TGraph& graph, SeqStreamIn& reads_iss, seqan::File<>& output_file,
              Options const& params, TReadsIndexSpec const )
  {
    /* typedefs */
    typedef Dna5QStringSet<> readsstringset_type;
    typedef SeedFinderTraits< typename TGraph::spec_type,
                              readsstringset_type, TReadsIndexSpec > finder_traits_type;
#ifdef PSI_STATS
    typedef SeedFinder< WithStats, finder_traits_type > finder_type;
#else
    typedef SeedFinder< NoStats, finder_traits_type > finder_type;
#endif
    typedef typename finder_type::traverser_type traverser_type;
    typedef typename finder_type::stats_type stats_type;
    typedef typename stats_type::timer_type timer_type;

    /* Get the main logger. */
    auto log = get_logger( "main" );
    auto tid = get_thread_id();
    /* Install seed finder singal handler for getting progress report. */
    std::signal( SIGUSR1, finder_type::stats_type::signal_handler );

    /* The seed finder for the input graph. */
    finder_type finder( graph, params.seed_len, params.gocc_threshold, params.max_mem );
    auto const& stats = finder.get_stats();
    /* Prepare (load or create) genome-wide paths. */
    log->info( "Looking for an existing path index..." );
    /* Load the genome-wide path index for the graph if available. */
    if ( finder.load_path_index( params.pindex_path,
                                 params.context,
                                 params.step_size,
                                 params.dindex_min_ris,
                                 params.dindex_max_ris ) ) {
      log->info( "The path index has been found and loaded." );
    }
    /* No genome-wide path index requested. */
    else if ( params.path_num == 0 ) {
      log->info( "No path has been specified. Skipping path indexing..." );
    }
    else {
      log->info( "No valid path index found. Creating the path index..." );
      log->info( "Selecting {} different path(s) in the graph...", params.path_num );
      finder.create_path_index( params.path_num, params.patched,
                                params.context, params.step_size,
                                params.dindex_min_ris, params.dindex_max_ris,
                                [&log]( std::string const& msg ) {
                                  log->info( msg );
                                },
                                [&log]( std::string const& msg ) {
                                  log->warn( msg );
                                } );
      log->info( "Picked paths in {}.", stats.get_timer( "pick-paths", tid ).str() );
      log->info( "Indexed paths in {}.", stats.get_timer( "index-paths", tid ).str() );
      log->info( "Found uncovered loci in {}.", stats.get_timer( "find-uncovered", tid ).str() );
      log->info( "Created distance index in {}.", stats.get_timer( "index-distances", tid ).str() );
      log->info( "Saving path index..." );
      /* Serialize the indexed paths. */
      if ( params.pindex_path.empty() ) {
        log->warn( "No path index file is specified. Skipping..." );
      } else if ( !finder.serialize_path_index( params.pindex_path, params.step_size ) ) {
        log->warn( "Specified path index file is not writable. Skipping..." );
      } else {
        log->info( "Saved path index in {}.", stats.get_timer( "save-pindex", tid ).str() );
        log->info( "Saved distance index in {}.", stats.get_timer( "save-dindex", tid ).str() );
      }
    }
    log->info( "Number of starting loci (in {} nodes of total {}): {}",
        finder.get_nof_uniq_nodes(), finder.get_graph_ptr()->get_node_count(),
        finder.get_starting_loci().size() );

    if ( params.indexonly ) {
      log->info( "Skipping seed finding as requested..." );
      return;
    }

    unsigned long long int found = 0;
    std::unordered_set< Records< readsstringset_type >::TPosition > covered_reads;
    std::function< void(typename traverser_type::output_type const &) > write_callback =
      [&found, &output_file, &covered_reads]
      (typename traverser_type::output_type const & seed_hit) {
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
      auto seeds = finder.create_readrecord();
      auto traverser = finder.create_traverser();
      log->info( "Finding seeds..." );
      [[maybe_unused]] auto timer = timer_type( "seed-finding" );
      while ( true ) {
        log->info( "Loading a read chunk..." );
        {
          [[maybe_unused]] auto timer = timer_type( "load-chunk" );
          /* Load a chunk from reads set. */
          if ( !readRecords( chunk, reads_iss, params.chunk_size ) ) break;
        }
        log->info( "Fetched {} reads in {}.", length( chunk ),
                   timer_type::get_duration_str( "load-chunk" ) );
        /* Give the current chunk to the finder. */
        finder.get_seeds( seeds, chunk, params.distance );
        auto seeds_index = finder.index_reads( seeds );
        log->info( "Seeding done in {}.", stats.get_timer( "seeding", tid ).str() );
        log->info( "Finding all seeds..." );
        finder.seeds_all( seeds, seeds_index, traverser, write_callback );
        log->info( "Found seeds on paths in {}.", stats.get_timer( "seeds-on-paths", tid ).str() );
        log->info( "Found seeds off paths in {}.", stats.get_timer( "seeds-off-paths", tid ).str() );
        log->info( "Verified distance constraints in {}.", stats.get_timer( "query-dindex", tid ).str() );
      }
    }
    log->info( "Found seed in {}.", timer_type::get_duration_str( "seed-finding" ) );
    report( finder, covered_reads, found );
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
  log->info( "- Seed genome occurrence count threshold: {}", options.gocc_threshold );
  log->info( "- Maximum number of MEMs on paths: {}", options.max_mem );
  log->info( "- Distance index minimum read insert size: {}", options.dindex_min_ris );
  log->info( "- Distance index maximum read insert size: {}", options.dindex_max_ris );
  log->info( "- Temporary directory: '{}'", get_tmpdir() );
  log->info( "- Output file: '{}'", options.output_path );

  log->info( "Loading input graph from file '{}'...", options.rf_path );
  gum::SeqGraph< gum::Succinct > graph;
  gum::util::load( graph, options.rf_path, true );
  if ( gum::util::ids_in_topological_order( graph ) ) {
    log->info( "Input graph node IDs are in topological sort order." );
  }
  else log->warn( "Input graph node IDs are NOT in topological sort order." );

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
  setValidValues( vgfile_arg, "vg gfa" );
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
  // seed genome occurrence count threshold
  addOption( parser,
             seqan::ArgParseOption( "r", "gocc-threshold",
                                    "Seed genome occurrence count threshold (no threshold by default).",
                                    seqan::ArgParseArgument::INTEGER, "INT" ) );
  setDefaultValue( parser, "r", 0 );
  // maximum number of MEMs on paths
  addOption( parser,
             seqan::ArgParseOption( "E", "max-mem",
                                    "Maximum number of MEMs on paths (default: find all).",
                                    seqan::ArgParseArgument::INTEGER, "INT" ) );
  setDefaultValue( parser, "E", 0 );
  // distance index minimum read insert size
  addOption( parser,
             seqan::ArgParseOption( "m", "min-insert-size",
                                    "Distance index minimum read insert size (no distance indexing by default).",
                                    seqan::ArgParseArgument::INTEGER, "INT" ) );
  setDefaultValue( parser, "m", 0 );
  // distance index maximum read insert size
  addOption( parser,
             seqan::ArgParseOption( "M", "max-insert-size",
                                    "Distance index maximum read insert size (minimum insert size by default).",
                                    seqan::ArgParseArgument::INTEGER, "INT" ) );
  setDefaultValue( parser, "M", 0 );
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
  getOptionValue( options.gocc_threshold, parser, "gocc-threshold" );
  getOptionValue( options.max_mem, parser, "max-mem" );
  getOptionValue( options.dindex_min_ris, parser, "min-insert-size" );
  getOptionValue( options.dindex_max_ris, parser, "max-insert-size" );
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
  if ( options.dindex_max_ris == 0 ) options.dindex_max_ris = options.dindex_min_ris;
}


  inline seqan::ArgumentParser::ParseResult
parse_args( Options& options, int argc, char* argv[] )
{
  // setup ArgumentParser.
  seqan::ArgumentParser parser( "psikt" );
  setup_argparser( parser );

  // Embedding program's meta data and build information.
  setShortDescription( parser, SHORT_DESC );
#ifdef PSI_GIT_REVISION
    setVersion( parser, REVISION+1 /* skip the leading 'v' before version number */ );
#else
    setVersion( parser, VERSION );
#endif
#ifdef PSI_GIT_COMMIT_DATE
  setDate( parser, PSI_GIT_COMMIT_DATE );
#endif  // PSI_GIT_COMMIT_DATE
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

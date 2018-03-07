/**
 *    @file  grem.cc
 *   @brief  grem main program.
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

#include "vargraph.h"
#include "mapper.h"
#include "traverser.h"
#include "pathset.h"
#include "sequence.h"
#include "seed.h"
#include "utils.h"
#include "options.h"
#include "stat.h"
#include "logger.h"
#include "release.h"

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
// TODO: Add value_t< T > typedef as typename seqan::Value< T >::Type


template< typename TMapper >
  void
signal_handler( int signal )
{
  std::cout << std::endl << "Report requested by SIGUSR1" << std::endl
            << "---------------------------" << std::endl;
  std::cout << "Elapsed time in traversal phase: "
            << Stat< TMapper >::Type::get_lap( "traverse" ).count() << " us"
            << std::endl;
  auto pos = Stat< TMapper >::Type::get_lastproc_locus().load();
  std::cout << "Current node: (" << pos.node_id << ", " << pos.offset << ")"
            << std::endl;
  auto idx = Stat< TMapper >::Type::get_lastdone_locus_idx().load();
  auto total = Stat< TMapper >::Type::get_total_nof_loci();
  unsigned int wlen = std::to_string( total ).length();
  std::cout << "Progress: " << std::setw(wlen) << idx << " / " << std::setw(wlen)
            << total << " [%" << std::setw(3) << idx * 100 / total << "]" << std::endl;
}


  void
default_signal_handler( int signal ) { /* Do nothing */ }

template< typename TMapper, typename TSet >
    void
  report( TMapper& mapper, TSet& covered_reads, unsigned long long int found )
  {
    auto log = get_logger( "main" );
    log->info( "Total number of starting loci: {}", mapper.get_starting_loci().size() );
    log->info( "Total number of seeds found: {}", found );
    log->info( "Total number of reads covered: {}", covered_reads.size() );
    log->info( "Total number of 'godown' operations: {}",
      TMapper::traverser_type::stats_type::get_total_nof_godowns() );

    log->info( "All Timers" );
    log->info( "----------" );
    for ( const auto& timer : Timer::get_timers() ) {
      log->info( "{}: {} us", timer.first, Timer::get_duration( timer.first ).count() );
    }
  }

template< typename TPathSet, typename TMapper >
    void
  prepare_paths_index( TPathSet& paths, TMapper& mapper, bool paths_index, bool patched,
      const std::string& paths_index_file, unsigned int path_num )
  {
    /* Get the main logger. */
    auto log = get_logger( "main" );
    log->info( "Loading paths index..." );
    /* Load the genome-wide paths for the variation graph if available. */
    if ( paths_index && paths.load( paths_index_file, mapper.get_vargraph() ) ) {
      log->info( "Paths index found. Loaded." );
      return;
    }
    /* No genome-wide paths requested. */
    if ( path_num == 0 ) {
      log->info( "Specified number of path is 0. Skipping paths indexing..." );
    }
    else {
      log->info( "No valid paths index found. Picking paths..." );
      log->info( "Picking {} different path(s) on the graph...", path_num );
      /* Generate the requested number of genome-wide paths. */
      mapper.pick_paths( paths, path_num, patched );
      log->info( "Picked paths in {} us.", Timer::get_duration( "pick-paths" ).count() );
      {
        auto timer = Timer( "sort-paths" );
        log->info( "Sorting the paths..." );
        /* Sort the paths. */
        paths.sort();
      }
      log->info( "Sorted paths in {} us.", Timer::get_duration( "sort-paths" ).count() );
      {
        auto timer = Timer( "index-paths" );
        log->info( "Indexing the paths..." );
        /* Index the paths. */
        paths.create_index();
      }
      log->info( "Indexed paths in {} us.", Timer::get_duration( "index-paths" ).count() );
      {
        auto timer = Timer( "save-paths" );
        log->info( "Saving paths index..." );
        /* Serialize the indexed paths. */
        if ( !paths_index ) {
          log->warn( "No paths index file is specified. Skipping..." );
        } else if ( !paths.serialize( paths_index_file ) ) {
          log->warn( "Specified paths index file is not writable. Skipping..." );
        }
      }
      log->info( "Saved paths in {} us.", Timer::get_duration( "save-paths" ).count() );
    }
  }


template< typename TIndexSpec  >
    void
  find_seeds( VarGraph& vargraph, SeqFileIn& reads_infile, seqan::File<>& output_file,
      unsigned int seed_len, unsigned int chunk_size, unsigned int start_every,
      unsigned int path_num, unsigned int context, bool paths_index, bool patched,
      const std::string& paths_index_file, bool nomapping, bool dumpstarts,
      const std::string& starts_path, TIndexSpec const /* Tag */ )
  {
    /* Get the main logger. */
    auto log = get_logger( "main" );
    /* Method typedefs. */
    typedef seqan::Index< Dna5QStringSet< grem::Dependent >, TIndexSpec > TIndex;
    typedef typename Traverser< TIndex, BFS, ExactMatching >::Type TTraverser;
    typedef Mapper< TTraverser > TMapper;
    /* The mapper for input variation graph. */
    TMapper mapper( &vargraph, seed_len );
    /* Install mapper singal handler for getting progress report. */
    std::signal( SIGUSR1, signal_handler< TMapper > );
    /* Genome-wide paths set in lazy mode. */
    PathSet< VarGraph, DiskString, grem::FMIndex<>, Forward > paths( context, true );
    /* Prepare (load or create) genome-wide paths. */
    prepare_paths_index( paths, mapper, paths_index, patched, paths_index_file, path_num );

    log->info( "Selecting starting loci..." );
    /* Locate starting loci. */
    mapper.add_all_loci( paths, seed_len, start_every );
    log->info( "Selected starting loci in {} us.",
        Timer::get_duration( "add-starts" ).count() );
    log->info( "Number of starting loci selected (in {} nodes): {}",
        mapper.get_vargraph()->node_count, mapper.get_starting_loci().size() );

    // :TODO:Tue Oct 17 22:37:\@cartoonist: save starting points as a part of paths index?
    if ( dumpstarts ) {
      log->info( "Dumping starting loci..." );
      std::ofstream starts_ofs( starts_path, std::ofstream::out );
      for ( const auto& p : mapper.get_starting_loci() ) {
        starts_ofs << p.node_id() << "\t" << p.offset() << "\t"
                   << vargraph.node_length( p.node_id() ) << std::endl;
      }
    }

    if ( nomapping ) {
      log->info( "Skipping mapping as requested..." );
      return;
    }

    /* All DNA reads to be mapped. */
    Records< Dna5QStringSet<> > reads;
    log->info( "Reading all reads from file..." );
    {
      auto timer = Timer( "read-reads" );
      /* Load all reads from file. */
      readRecords( reads, reads_infile );
    }
    log->info( "Read {} reads in {} us.", length( reads ),
        Timer::get_duration( "read-reads" ).count() );

    /* If chunksize is zero (unspecified), set it to the no of reads; i.e. one chunk. */
    if ( chunk_size == 0 ) {
      chunk_size = length( reads );
    }

    unsigned long long int found = 0;
    std::unordered_set< Records< Dna5QStringSet<> >::TPosition > covered_reads;
    std::function< void(typename TTraverser::output_type const &) > write_callback =
      [&found, &output_file, &covered_reads]
      (typename TTraverser::output_type const & seed_hit){
      ++found;
      write( output_file, &seed_hit.node_id, 1 );
      write( output_file, &seed_hit.node_offset, 1 );
      write( output_file, &seed_hit.read_id, 1 );
      write( output_file, &seed_hit.read_offset, 1 );
      covered_reads.insert(seed_hit.read_id);
    };

    /* Reads are mapped in chunks. */
    Records< Dna5QStringSet< grem::Dependent > > reads_chunk;
    log->info( "Finding seeds..." );
    {
      auto timer = Timer( "seed-finding" );
      while (true)
      {
        log->info( "Loading the next reads chunk..." );
        {
          auto timer = Timer( "load-chunk" );
          /* Load a chunk from reads set. */
          if ( ! load_chunk( reads_chunk, reads, chunk_size ) ) break;
        }
        log->info( "Loaded reads chunk in {} us.",
            Timer::get_duration( "load-chunk" ).count() );
        /* Give the current chunk to the mapper. */
        mapper.set_reads( std::move( reads_chunk ) );
        log->info( "Finding seeds on paths..." );
        auto pre_found = found;
        /* Find seeds on genome-wide paths. */
        mapper.seeds_on_paths( paths, write_callback );
        log->info( "Found seed on paths in {} us.",
            Timer::get_duration( "paths-seed-find" ).count() );
        log->info( "Total number of seeds found on paths: {}", found - pre_found );
        log->info( "Traversing..." );
        /* Find seeds on variation graph by traversing starting loci. */
        mapper.traverse ( write_callback );
        log->info( "Traversed in {} us.", Timer::get_duration( "traverse" ).count() );
      }
    }
    log->info( "found seed in {} us.", Timer::get_duration( "seed-finding" ).count() );
    report( mapper, covered_reads, found );
  }


  void
startup( const Options & options )
{
  auto log = get_logger( "main" );
  log->info( "Parameters:" );
  log->info( "- Seed length: {}", options.seed_len );
  log->info( "- Number of paths: {}", options.path_num );
  log->info( "- Context length (used in patching): {}", options.context );
  log->info( "- Patched: {}", ( options.patched ? "yes" : "no" ) );
  log->info( "- Paths index file: '{}'", options.paths_index_file );
  log->info( "- Reads chunk size: {}", options.chunk_size );
  log->info( "- Reads index type: {}", index_to_str(options.index) );
  log->info( "- Starting loci interval: {}", options.start_every );

  log->info( "Opening file '{}'...", options.fq_path );
  SeqFileIn reads_infile;
  if ( !open( reads_infile, options.fq_path.c_str() ) )
  {
    std::string msg = "could not open file '" + options.fq_path + "'!";
    log->error( msg );
    throw std::runtime_error( msg );
  }

  log->info( "Loading the graph from file '{}'...", options.rf_path );
  std::ifstream ifs( options.rf_path, std::ifstream::in | std::ifstream::binary );
  if( !ifs )
  {
    std::string msg = "could not open file '" + options.rf_path + "'!";
    log->error( msg );
    throw std::runtime_error( msg );
  }

  bool xg_format = true;
  if ( ends_with( options.rf_path, ".vg" ) ) {
    xg_format = false;
  }
  VarGraph vargraph( ifs, xg_format );

  seqan::File<> output_file;
  auto mode = seqan::OPEN_CREATE | seqan::OPEN_WRONLY;
  if ( !open( output_file, options.output_path.c_str(), mode ) ) {
    std::string msg = "could not open file '" + options.output_path + "'!";
    log->error( msg );
    throw std::runtime_error( msg );
  }

  switch ( options.index ) {
    case IndexType::Wotd: find_seeds ( vargraph,
                              reads_infile,
                              output_file,
                              options.seed_len,
                              options.chunk_size,
                              options.start_every,
                              options.path_num,
                              options.context,
                              options.paths_index,
                              options.patched,
                              options.paths_index_file,
                              options.nomapping,
                              options.dumpstarts,
                              options.starts_path,
                              UsingIndexWotd() );
                          break;
    case IndexType::Esa: find_seeds ( vargraph,
                             reads_infile,
                             output_file,
                             options.seed_len,
                             options.chunk_size,
                             options.start_every,
                             options.path_num,
                             options.context,
                             options.paths_index,
                             options.patched,
                             options.paths_index_file,
                             options.nomapping,
                             options.dumpstarts,
                             options.starts_path,
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

  // output file
  // :TODO:Sat Oct 21 00:06:\@cartoonist: output should be alignment in the GAM format.
  seqan::ArgParseOption outfile_arg( "o", "output", "Output file.",
      seqan::ArgParseArgument::OUTPUT_FILE, "OUTPUT_FILE" );
  addOption( parser, outfile_arg );
  setDefaultValue( parser, "o", "out.gam" );

  // paths index file
  seqan::ArgParseOption paths_index_arg ( "I", "paths-index", "Paths index file.",
      seqan::ArgParseArgument::STRING, "PATHS_INDEX_FILE" );
  addOption ( parser, paths_index_arg );

  // seed length -- **required** option.
  addOption(parser, seqan::ArgParseOption("l", "seed-length", "Seed length.",
                                          seqan::ArgParseArgument::INTEGER, "INT"));
  setRequired(parser, "l");

  // chunk size -- **required** option.
  addOption(parser, seqan::ArgParseOption("c", "chunk-size",
        "Reads chunk size. Set it to zero to consider all reads as one chunk (default).",
        seqan::ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "c", 0);

  // starting loci interval
  addOption(parser, seqan::ArgParseOption("e", "start-every", "Start from every given "
                                          "number of loci in all nodes. If it is set to"
                                          " 1, it means start from all positions in a "
                                          "node.", seqan::ArgParseArgument::INTEGER,
                                          "INT"));
  setDefaultValue(parser, "e", 1);

  // number of paths
  addOption( parser, seqan::ArgParseOption( "n", "path-num", "Number of paths from "
        "the variation graph in hybrid approach.", seqan::ArgParseArgument::INTEGER,
        "INT" ) );
  setDefaultValue(parser, "n", 0);

  // whether use patched paths or full genome-wide paths
  addOption(parser, seqan::ArgParseOption("P", "no-patched", "Use full genome-wide paths."));

  // context in patching the paths
  addOption( parser, seqan::ArgParseOption( "t", "context", "Context length in patching.",
        seqan::ArgParseArgument::INTEGER, "INT" ) );
  setDefaultValue( parser, "t", 0 );

  // index
  addOption(parser, seqan::ArgParseOption("i", "index", "Index type for indexing reads.",
                                          seqan::ArgParseArgument::STRING, "INDEX"));
  setValidValues(parser, "i", "SA ESA WOTD DFI QGRAM FM");
  setDefaultValue(parser, "i", "WOTD");

  // no map -- just do the paths indexing phase
  addOption(parser, seqan::ArgParseOption("x", "only-index", "Only build paths index and skip mapping."));

  // dump starting loci
  seqan::ArgParseOption startsfile_arg( "s", "dump-starts", "Dump starting loci.",
      seqan::ArgParseArgument::OUTPUT_FILE, "DUMP_FILE" );
  addOption( parser, startsfile_arg );

  // log file
  seqan::ArgParseOption logfile_arg("L", "log-file",
                                    "Sets default log file for existing and future loggers.",
                                    seqan::ArgParseArgument::OUTPUT_FILE, "LOG_FILE");
  addOption(parser, logfile_arg);
  setDefaultValue(parser, "L", "grem.log");

  // no log to file
  addOption(parser, seqan::ArgParseOption("Q", "no-log-file", "Disable writing logs to file (overrides \\fB-L\\fP)."));

  // quiet -- no output to console
  addOption(parser, seqan::ArgParseOption("q", "quiet", "Quiet mode. No output will be printed to console."));

  // no colored output
  addOption(parser, seqan::ArgParseOption("C", "no-color", "Do not use a colored output."));

  // disable logging
  addOption(parser, seqan::ArgParseOption("D", "disable-log", "Disable logging completely."));

  // verbosity options
  addOption(parser, seqan::ArgParseOption("v", "verbose",
                                          "Activates maximum verbosity."));
}


  inline void
get_option_values ( Options & options, seqan::ArgumentParser & parser )
{
  std::string indexname;

  getOptionValue( options.fq_path, parser, "fastq" );
  getOptionValue( options.output_path, parser, "output" );
  getOptionValue( options.seed_len, parser, "seed-length" );
  getOptionValue( options.chunk_size, parser, "chunk-size" );
  getOptionValue( options.start_every, parser, "start-every" );
  getOptionValue( options.path_num, parser, "path-num" );
  getOptionValue( options.context, parser, "context" );
  options.paths_index = isSet( parser, "paths-index" );
  options.patched = !isSet( parser, "no-patched" );
  getOptionValue( options.paths_index_file, parser, "paths-index" );
  getOptionValue( indexname, parser, "index" );
  options.nomapping = isSet( parser, "only-index" );
  options.dumpstarts = isSet( parser, "dump-starts");
  getOptionValue( options.starts_path, parser, "dump-starts" );
  getOptionValue( options.log_path, parser, "log-file" );
  options.nologfile = isSet( parser, "no-log-file" );
  options.quiet = isSet( parser, "quiet" );
  options.nocolor = isSet( parser, "no-color" );
  options.nolog = isSet( parser, "disable-log" );
  options.verbose = isSet( parser, "verbose" );
  getArgumentValue( options.rf_path, parser, 0 );

  options.index = index_from_str( indexname );
}


  inline seqan::ArgumentParser::ParseResult
parse_args( Options& options, int argc, char* argv[] )
{
  // setup ArgumentParser.
  seqan::ArgumentParser parser( PACKAGE );
  setup_argparser( parser );

  // Embedding program's meta data and build information.
  setShortDescription( parser, SHORT_DESC );
  setVersion( parser, VERSION );
  // :TODO:Thu Apr 06 02:06:\@cartoonist: date should be specified by autoconf.
  setDate( parser, __DATE__ );
  addDescription( parser, LONG_DESC );

  std::ostringstream hold_buf_stdout;
  std::ostringstream hold_buf_stderr;
  // parse command line.
  auto res = seqan::parse( parser, argc, argv, hold_buf_stdout, hold_buf_stderr );
  // print the banner in help or version messages.
  if ( res == seqan::ArgumentParser::PARSE_HELP ||
      res == seqan::ArgumentParser::PARSE_VERSION )
  {
    std::cout << BANNER << std::endl;
  }
  // print the buffer.
  std::cout << hold_buf_stderr.str();
  std::cout << hold_buf_stdout.str();

  // only extract options if the program will continue after parse_args()
  if ( res != seqan::ArgumentParser::PARSE_OK ) return res;

  get_option_values ( options, parser );

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
  if ( res != seqan::ArgumentParser::PARSE_OK )
    return res == seqan::ArgumentParser::PARSE_ERROR;

  /* Verify that the version of the library that we linked against is
   * compatible with the version of the headers we compiled against.
   */
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  /* Install signal handler -- default handler */
  std::signal( SIGUSR1, default_signal_handler );
  /* Configure loggers */
  config_logger( options );
  /* Start mapping... */
  startup( options );

  /* Close all loggers. */
  drop_all_loggers();
  /* Delete all global objects allocated by libprotobuf. */
  google::protobuf::ShutdownProtobufLibrary();

  return EXIT_SUCCESS;
}

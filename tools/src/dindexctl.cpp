/**
 *    @file  dindexctl.cpp
 *   @brief  Distance index hacking tool
 *
 *  A tool for hacking (e.g. compressing, merging) distance indices.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Thu Dec 17, 2020  01:22
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2020, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include <cxxopts.hpp>
#include <gum/graph.hpp>
#include <gum/io_utils.hpp>
#include <psi/graph.hpp>
#include <psi/crs_matrix.hpp>
#include <psi/seed_finder.hpp>


using namespace psi;

constexpr const char* LONG_DESC = ( "dindexctl\n"
                                    "---------\n"
                                    "Hacking tool for distance indices\n" );
// Default values for command line arguments
constexpr const char* DEFAULT_OUTPUT = "-";  // stdout
constexpr const char* DEFAULT_SAMPLING_RATE = "0.001";
constexpr const char* DEFAULT_RNDSEED = "0";


namespace rnd {
  thread_local static std::mt19937 lgen;
  thread_local static unsigned int iseed = 0;

  inline void
  init_gen( unsigned int seed=0 )
  {
    if ( seed != 0 ) {
      iseed = seed;
      lgen.seed( seed );
    }
  }

  inline std::mt19937&
  get_gen( )
  {
    if ( iseed == 0 ) return psi::random::gen;
    else return lgen;
  }
}  /* ---  end of namespace rnd  --- */

void
config_parser( cxxopts::Options& options )
{
  options.add_options( "general" )
      ( "o, output", "Write to this file instead of stdout",
        cxxopts::value< std::string >()->default_value( DEFAULT_OUTPUT ) )
      ( "b, basic-mode", "Consider indices as Basic CRS matrices" )
      ( "h, help", "Print this message and exit" )
      ;

  options.add_options( "compress" )
      ( "d, min-insert-size", "Distance index minimum read insert size",
        cxxopts::value< unsigned int >() )
      ( "D, max-insert-size", "Distance index maximum read insert size",
        cxxopts::value< unsigned int >() )
      ( "g, graph", "Corresponding graph file (vg or gfa)",
        cxxopts::value< std::string >() )
      ( "V, verify", "Verify if the distance index is compressed" )
      ( "r, sample-rate", "Node sampling rate for verification",
        cxxopts::value< float >()->default_value( DEFAULT_SAMPLING_RATE ) )
      ( "S, random-seed", "Seed for random generator",
        cxxopts::value< unsigned int >()->default_value( DEFAULT_RNDSEED ) )
      ;

  options.add_options( "merge" )
      ( "1, first-range", "Distance constraint range of the first index (comma-separated: min,max)",
        cxxopts::value< std::vector< unsigned int > >() )
      ( "2, second-range", "Distance constraint range of the second index (comma-separated: min,max)",
        cxxopts::value< std::vector< unsigned int > >() )
      ;

  options.add_options( "positional" )
      ( "command", "Operation type", cxxopts::value< std::string >() )
      ( "prefix", "Path index prefix", cxxopts::value< std::string >() )
      ;
  options.parse_positional( { "command", "prefix" } );
}

cxxopts::ParseResult
parse_opts( cxxopts::Options& options, int& argc, char**& argv )
{
  auto result = options.parse( argc, argv );

  if ( !result.count( "command" ) ) {  // no command specified
    options.positional_help( "COMMAND" );
    auto help_message = ( options.help( { "general" } )
                          + "\n COMMANDS:\n"
                          + "  compress\tCompress a distance index\n"
                          + "  merge\t\tMerge two distance indices" );
    if ( result.count( "help" ) )
    {
      std::cout << help_message << std::endl;
      throw EXIT_SUCCESS;
    }
    std::cerr << help_message << "\n" /* extra vertical space */ << std::endl;
    throw cxxopts::OptionParseException( "No command specified" );
  }
  else if ( result[ "command" ].as< std::string >() == "compress" ) {  // compress
    options.custom_help( "compress [OPTION...]" );
    options.positional_help( "PREFIX" );
    if ( result.count( "help" ) ) {
      std::cout << options.help( { "general", "compress" } ) << std::endl;
      throw EXIT_SUCCESS;
    }

    if ( result.count( "basic-mode" ) && !result.count( "graph" ) ) {
      throw cxxopts::OptionParseException( "Graph file must be specified" );
    }
    if ( result.count( "basic-mode" ) && !readable( result[ "graph" ].as< std::string >() ) ) {
      throw cxxopts::OptionParseException( "Graph file not found" );
    }

    if ( !result.count( "min-insert-size" ) ) {
      throw cxxopts::OptionParseException( "Minimum insert size must be specified" );
    }

    if ( !result.count( "max-insert-size" ) ) {
      throw cxxopts::OptionParseException( "Maximum insert size must be specified" );
    }
  }
  else if ( result[ "command" ].as< std::string >() == "merge" ) {  // merge
    options.custom_help( "merge [OPTION...]" );
    options.positional_help( "PREFIX" );
    if ( result.count( "help" ) ) {
      std::cout << options.help( { "general", "merge" } ) << std::endl;
      throw EXIT_SUCCESS;
    }

    if ( !result.count( "first-range" ) ) {
      throw cxxopts::OptionParseException( "First distance constraint range must be specified" );
    }
    else {
      auto range = result[ "first-range" ].as< std::vector< unsigned int > >();
      if ( range.size() != 2 ) {
        throw cxxopts::OptionParseException( "Invalid range for the first constraint" );
      }
    }

    if ( !result.count( "second-range" ) ) {
      throw cxxopts::OptionParseException( "Second distance constraint range must be specified" );
    }
    else {
      auto range = result[ "second-range" ].as< std::vector< unsigned int > >();
      if ( range.size() != 2 ) {
        throw cxxopts::OptionParseException( "Invalid range for the second constraint" );
      }
    }
  }
  else {
    throw cxxopts::OptionParseException( "Unknown command '" +
                                         result[ "command" ].as< std::string >() + "'" );
  }

  /* Verifying positional arguments */
  if ( !result.count( "prefix" ) ) {
    throw cxxopts::OptionParseException( "Index prefix must be specified" );
  }
  if ( !readable( result[ "prefix" ].as< std::string >() ) ) {
    throw cxxopts::OptionParseException( "Index file not found" );
  }

  return result;
}

template< typename TCRSMatrix, typename TGraph >
bool
verify_compressed_distance_matrix( TCRSMatrix const& cdi, TCRSMatrix const& udi,
                                   TGraph const& g )
{
  typedef typename TCRSMatrix::ordinal_type ordinal_type;
  typedef typename TCRSMatrix::size_type size_type;
  typedef typename TGraph::rank_type rank_type;

  rank_type cnode_rank = 0;  // current node rank
  size_type start = 0;    // row start index
  size_type cstart = 0;    // row start index
  size_type end;          // row end index
  ordinal_type nloc = 0;     // next node loci index
  for ( ordinal_type nrow = 0; nrow < udi.numRows(); ++nrow ) {
    if ( nrow == nloc ) {
      ++cnode_rank;
      if ( cnode_rank == g.get_node_count() ) nloc = udi.numRows();
      else nloc = gum::util::id_to_charorder( g, g.rank_to_id( cnode_rank + 1 ) );
    }
    assert( nrow < nloc );
    end = udi.rowMap( nrow + 1 );
    for ( ; start < end; ++start ) {
      if ( nrow <= udi.entry( start ) && udi.entry( start ) < nloc ) continue;
      else if ( udi.entry( start ) == cdi.entry( cstart ) ) ++cstart;
      else return false;
    }
  }
  assert( start == udi.nnz() );
  assert( cstart == cdi.nnz() );
  std::cout << "Reduced the distance matrix by " << start - cstart << " elements."
            << std::endl;
  return true;
}

template< typename TCRSMatrixRange, typename TCRSMatrixBasic,
          /* restrict it to range group CRS matrices */
          typename=std::enable_if_t< std::is_same< typename crs_matrix::Group< typename TCRSMatrixRange::spec_type >::type, crs_matrix::RangeGroup >::value >,
          typename=std::enable_if_t< std::is_same< typename crs_matrix::Group< typename TCRSMatrixBasic::spec_type >::type, crs_matrix::BasicGroup >::value > >
void
compress( cxxopts::ParseResult& res, crs_matrix::RangeGroup /* tag */ )
{
  typedef TCRSMatrixRange crsmat_range_type;
  typedef TCRSMatrixBasic crsmat_basic_type;

  std::string pindex_prefix = res[ "prefix" ].as< std::string >();
  std::string output = res[ "output" ].as< std::string >();
  unsigned int min_size = res[ "min-insert-size" ].as< unsigned int >();
  unsigned int max_size = res[ "max-insert-size" ].as< unsigned int >();
  crsmat_basic_type basic_dindex;
  crsmat_range_type range_dindex;
  auto index_path = psi::SeedFinder<>::get_distance_index_path( pindex_prefix, min_size, max_size );

  std::cout << "Loading distance index..." << std::endl;
  std::ifstream ifs( index_path, std::ifstream::in | std::ifstream::binary );
  if ( !ifs ) throw std::runtime_error( "distance matrix cannot be opened" );

  if ( res.count( "verify" ) ) {
    std::cout << "[WARNING] There is no verification procedure for "
              << "Range-based distance indices" << std::endl;
    range_dindex.load( ifs );
    std::cout << "Loaded distance index ("
              << range_dindex.numRows() << "x" << range_dindex.numCols() << ") has "
              << range_dindex.nnz() << " non-zero elements." << std::endl;
  }
  else {
    basic_dindex.load( ifs );
    std::cout << "Loaded distance index ("
              << basic_dindex.numRows() << "x" << basic_dindex.numCols() << ") has "
              << basic_dindex.nnz() << " non-zero elements." << std::endl;

    std::cout << "Compressing distance index..." << std::endl;
    range_dindex.assign( basic_dindex );
    std::cout << "Compressed distance index ("
              << range_dindex.numRows() << "x" << range_dindex.numCols() << ") has "
              << range_dindex.nnz() << " non-zero elements." << std::endl;

    std::cout << "Verifying compressed distance index..." << std::endl;
    if ( range_dindex.numRows() != basic_dindex.numRows() ||
         range_dindex.numCols() != basic_dindex.numCols() ||
         range_dindex.nnz() != basic_dindex.nnz() )  {
      std::cerr << "Verification failed!" << std::endl;
      throw EXIT_FAILURE;
    }

    std::cout << "Serialising compressed distance index..." << std::endl;
    std::ofstream ofs( output, std::ofstream::out | std::ofstream::binary );
    if ( !ofs ) throw std::runtime_error( "output file cannot be opened" );
    range_dindex.serialize( ofs );
  }
}

template< typename TCRSMatrix, typename TMutableCRSMatrix,
          /* restrict it to basic group CRS matrices */
          typename=std::enable_if_t< std::is_same< typename crs_matrix::Group< typename TCRSMatrix::spec_type >::type, crs_matrix::BasicGroup >::value >,
          typename=std::enable_if_t< std::is_same< typename crs_matrix::Group< typename TMutableCRSMatrix::spec_type >::type, crs_matrix::BasicGroup >::value > >
void
compress( cxxopts::ParseResult& res, crs_matrix::BasicGroup /* tag */ )
{
  typedef TCRSMatrix crsmat_type;
  typedef TMutableCRSMatrix mutable_crsmat_type;

  typedef gum::SeqGraph< gum::Succinct > graph_type;
  typedef typename graph_type::id_type id_type;
  typedef typename graph_type::rank_type rank_type;
  typedef typename graph_type::offset_type offset_type;

  std::string graph_path = res[ "graph" ].as< std::string >();
  std::string pindex_prefix = res[ "prefix" ].as< std::string >();
  std::string output = res[ "output" ].as< std::string >();
  unsigned int min_size = res[ "min-insert-size" ].as< unsigned int >();
  unsigned int max_size = res[ "max-insert-size" ].as< unsigned int >();
  graph_type graph;
  crsmat_type dindex;
  auto index_path = psi::SeedFinder<>::get_distance_index_path( pindex_prefix, min_size, max_size );

  std::cout << "Loading input graph..." << std::endl;
  gum::util::load( graph, graph_path, true );
  std::string sort_status = gum::util::ids_in_topological_order( graph ) ? "" : "not ";
  std::cout << "Input graph node IDs are " << sort_status << "in topological sort order."
            << std::endl;

  std::cout << "Loading distance index..." << std::endl;
  std::ifstream ifs( index_path, std::ifstream::in | std::ifstream::binary );
  if ( !ifs ) throw std::runtime_error( "distance matrix cannot be opened" );
  dindex.load( ifs );
  std::cout << "Loaded distance index ("
            << dindex.numRows() << "x" << dindex.numCols() << ") has " << dindex.nnz()
            << " non-zero elements." << std::endl;

  if ( res.count( "verify" ) ) {
    std::cout << "Verifying distance index for compression..." << std::endl;
    std::uniform_real_distribution< float > dis( 0, 1 );
    float srate = res[ "sample-rate" ].as< float >();
    unsigned int seed = res[ "random-seed" ].as< unsigned int >();
    ::rnd::init_gen( seed );
    bool success = graph.for_each_node(
      [&graph, &dindex, &dis, srate]( rank_type rank, id_type id ) {
        if ( dis( ::rnd::get_gen() ) >= srate ) return true;
        auto label_len = graph.node_length( id );
        auto charid = gum::util::id_to_charorder( graph, id );
        for ( offset_type i = 0; i < label_len; ++i ) {
          for ( offset_type j = i+1; j < label_len; ++j ) {
            if ( dindex( charid + i, charid + j ) ) return false;
          }
        }
        return true;
      });
    if ( success ) std::cout << "[PASS] Input distance index is compressed." << std::endl;
    else std::cerr << "[FAIL] Input distance index is not compressed!" << std::endl;
  }
  else {
    std::cout << "Compressing distance index..." << std::endl;
    crsmat_type cindex;
    cindex.assign( util::compress_distance_index< mutable_crsmat_type >( dindex, graph ) );
    std::cout << "Compressed distance index ("
              << cindex.numRows() << "x" << cindex.numCols() << ") has " << cindex.nnz()
              << " non-zero elements." << std::endl;

    std::cout << "Verifying compressed distance index..." << std::endl;
    if ( !verify_compressed_distance_matrix( cindex, dindex, graph ) ) {
      std::cerr << "Verification failed!" << std::endl;
      throw EXIT_FAILURE;
    }

    std::cout << "Serialising compressed distance index..." << std::endl;
    std::ofstream ofs( output, std::ofstream::out | std::ofstream::binary );
    if ( !ofs ) throw std::runtime_error( "output file cannot be opened" );
    cindex.serialize( ofs );
  }
}

template< typename TCRSMatrix1, typename TCRSMatrix2 >
void
compress( cxxopts::ParseResult& res )
{
  compress< TCRSMatrix1, TCRSMatrix2 >
    ( res, typename crs_matrix::Group< typename TCRSMatrix1::spec_type >::type{} );
}

template< typename TCRSMatrix >
bool
verify_merged_distance_matrix( TCRSMatrix& mdi, TCRSMatrix& di1, TCRSMatrix& di2,
                               crs_matrix::BasicGroup /* tag */ )
{
  typedef typename TCRSMatrix::ordinal_type ordinal_type;
  typedef typename TCRSMatrix::size_type size_type;

  size_type start = 0;    // row start index
  size_type end;          // row end index
  size_type sum_of_nnz = 0;
  for ( ordinal_type nrow = 0; nrow < mdi.numRows(); ++nrow ) {
    end = mdi.rowMap( nrow + 1 );
    for ( ; start < end; ++start ) {
      auto col = mdi.entry( start );
      unsigned int sum = di1( nrow, col ) + di2( nrow, col );
      sum_of_nnz += sum;
      if ( sum == 0 ) {
        std::cerr << "Merged matrix contains an invalid non-zero element at ("
                  << nrow << ", " << col << ")!" << std::endl;
        return false;
      }
    }
  }
  assert( start == mdi.nnz() );
  if ( sum_of_nnz != di1.nnz() + di2.nnz() ) return false;
  return true;
}

template< typename TCRSMatrix >
bool
verify_merged_distance_matrix( TCRSMatrix& mdi, TCRSMatrix& di1, TCRSMatrix& di2,
                               crs_matrix::RangeGroup /* tag */ )
{
  std::cout << "Skipping verification for Range CRS indices..." << std::endl;
  return true;
}

template< typename TCRSMatrix, typename TMutableCRSMatrix >
void
merge( cxxopts::ParseResult& res )
{
  typedef TCRSMatrix crsmat_type;
  typedef TMutableCRSMatrix crsmat_mutable_type;

  std::string pindex_prefix = res[ "prefix" ].as< std::string >();
  std::string output = res[ "output" ].as< std::string >();
  auto range1 = res[ "first-range" ].as< std::vector< unsigned int > >();
  auto range2 = res[ "second-range" ].as< std::vector< unsigned int > >();

  crsmat_type mindex;
  crsmat_mutable_type bmerged;

  {
    crsmat_type dindex1;
    crsmat_type dindex2;
    auto index_path1 = psi::SeedFinder<>::get_distance_index_path( pindex_prefix, range1[0], range1[1] );
    auto index_path2 = psi::SeedFinder<>::get_distance_index_path( pindex_prefix, range2[0], range2[1] );

    {
      std::cout << "Loading the first distance index '" << index_path1 << "'..." << std::endl;
      std::ifstream ifs( index_path1, std::ifstream::in | std::ifstream::binary );
      if ( !ifs ) throw std::runtime_error( "The first distance matrix cannot be opened" );
      dindex1.load( ifs );
      std::cout << "Loaded the first distance index ("
                << dindex1.numRows() << "x" << dindex1.numCols() << ") with "
                << dindex1.nnz() << " non-zero elements." << std::endl;
    }

    {
      std::cout << "Loading the second distance index '" << index_path2 << "'..." << std::endl;
      std::ifstream ifs( index_path2, std::ifstream::in | std::ifstream::binary );
      if ( !ifs ) throw std::runtime_error( "The second distance matrix cannot be opened" );
      dindex2.load( ifs );
      std::cout << "Loaded the second distance index ("
                << dindex2.numRows() << "x" << dindex2.numCols() << ") with "
                << dindex2.nnz() << " non-zero elements." << std::endl;
    }

    std::cout << "Merging distance indices..." << std::endl;
    bmerged = merge_distance_index< crsmat_mutable_type >( dindex1, dindex2 );
  }

  mindex.assign( bmerged );
  std::cout << "Merged distance index ("
            << mindex.numRows() << "x" << mindex.numCols() << ") has " << mindex.nnz()
            << " non-zero elements." << std::endl;

  std::cout << "Serialising merged distance index..." << std::endl;
  std::ofstream ofs( output, std::ofstream::out | std::ofstream::binary );
  if ( !ofs ) throw std::runtime_error( "output file cannot be opened" );
  mindex.serialize( ofs );
}

int
main( int argc, char* argv[] )
{
  cxxopts::Options options( argv[0], LONG_DESC );
  config_parser( options );

  try {
    auto res = parse_opts( options, argc, argv );
    std::string command = res[ "command" ].as< std::string >();
    bool basic_mode = res[ "basic-mode" ].as< bool >();

    typedef typename SeedFinder<>::crsmat_type crsmat_type;
    typedef make_basic_t< crsmat_type > crsmat_basic_type;

    typedef typename SeedFinder<>::mutable_crsmat_type crsmat_mut_type;
    typedef make_buffered_t< crsmat_basic_type > crsmat_basic_mut_type;

    if ( command == "compress" ) {
      if ( basic_mode ) compress< crsmat_basic_type, crsmat_basic_mut_type >( res );
      else compress< crsmat_type, crsmat_basic_type >( res );
    }
    else if ( command == "merge" ) {
      if ( basic_mode ) merge< crsmat_basic_type, crsmat_basic_mut_type >( res );
      else merge< crsmat_type, crsmat_mut_type >( res );
    }
    else {
      // should not reach here!
      assert( false );
    }
  }
  catch ( const cxxopts::OptionException& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch ( const int& rv ) {
    return rv;
  }

  return EXIT_SUCCESS;
}

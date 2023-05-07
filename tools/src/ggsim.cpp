/**
 *    @file  ggsim.cpp
 *   @brief  Graph genome haplotype and reads simulator
 *
 *  Simulate haplotypes or reads from a graph genome.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Jul 03, 2018  19:09
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include "ggsim.hpp"


template< typename TGraph >
inline void
simulate_haplotypes( TGraph const& graph, unsigned int ploidy, unsigned int seed,
                     std::function< void( std::string, Path< TGraph > ) > callback )
{
  // NOTE: This function assumes that there is a path corresponding to each region.
  auto hap_itr = begin( graph, Haplotyper< Random >() );
  auto hap_end = end( graph, Haplotyper< Random >() );
  Path< TGraph > hap_path( &graph );
  graph.for_each_path(
      [&]( auto path_rank, auto path_id ) {
        auto path_name = graph.path_name( path_id );
        typename TGraph::id_type start = *graph.path( path_id ).begin();
        hap_itr.reset( start, seed );
        for ( unsigned int i = 0; i < ploidy; ++i ) {
          get_rnd_full_haplotype( hap_path, hap_itr, hap_end );
          psi::initialize( hap_path );
          callback( path_name + "-" + std::to_string( i + 1 ), std::move( hap_path ) );
          hap_path.clear();
        }
        return true;
      } );
}

template< typename TReadType >
inline std::size_t
read_span( unsigned int readlen )
{
  std::size_t is_paired = std::is_same< TReadType, PairedEnd >::value;
  return ( 1 + is_paired ) * readlen + is_paired /* min one-base inner distance if paired */;
}

template< typename TReadType, typename TGraph >
inline std::vector< unsigned long int >
reads_dist( std::unordered_map< std::string, Path< TGraph > > const& haplotypes,
                    unsigned int readlen, unsigned long int numreads )
{
  std::vector< unsigned long int > dist;
  unsigned long int genome_size = 0;
  for ( const auto& hs : haplotypes ) {
    auto length = hs.second.get_sequence_len();
    if ( length >= read_span< TReadType >( readlen ) ) {
      genome_size += length;
      dist.push_back( length );
    }
    else {
      std::cerr << "Skipped haplotype '" << hs.first << "' due to its length."
                << std::endl;
      dist.push_back( 0 );
    }
  }

  unsigned long int sum = 0;
  for ( auto& d : dist ) {
    d = numreads * d / genome_size;
    sum += d;
  }
  assert( numreads >= sum );
  auto extra = numreads - sum;
  assert( extra < haplotypes.size() );
  for ( unsigned long int i = 0; i < extra; ++i ) ++dist[i];
  assert( static_cast< unsigned long int >( std::accumulate( dist.begin(), dist.end(), 0 ) ) == numreads );
  return dist;
}

template< typename TInputIter, typename TOutputIter >
inline std::pair< bool, TOutputIter >
impose_error( TInputIter first, TInputIter last, TOutputIter output, vg::Mapping* map_ptr,
              double errorrate, double indelrate, bool fwd, bool allow_ns, bool forced )
{
  std::uniform_real_distribution<> dis( 0.0, 1.0 );
  char bases[4] = { 'A', 'C', 'G', 'T' };
  event::Event e;
  int alt;

  auto count = last - first;
  auto edit_ptr = map_ptr->add_edit();
  while ( count --> 0 ) {
    assert( first != last );
    if ( dis( ::rnd::get_gen() ) < errorrate ) {
      if ( dis( ::rnd::get_gen() ) < indelrate ) {  // indel
        if ( dis( ::rnd::get_gen() ) < 0.5 || first == last - 1 ) {  // insertion
          alt = static_cast< int >( dis( ::rnd::get_gen() ) * 4 );
          --first;
          *output++ = bases[ alt ];
          e = event::insertion;
        }
        else {
          *output++ = *( ++first );  // deletion
          e = event::deletion;
        }
      }
      else {  // substitution
        alt = static_cast< int >( dis( ::rnd::get_gen() ) * 4 );
        if ( bases[ alt ] == *first ) alt = ( alt + 1 ) % 4;
        *output++ = bases[ alt ];
        e = event::mismatch;
      }
    }
    else {
      *output++ = *first++;  // no error
      e = event::match;
    }

    auto last_char = *( output - 1 );
    if ( !allow_ns && !forced && last_char == 'N' ) return { false, output };

    if ( e == event::match ) {
      if ( edit_ptr == nullptr ) edit_ptr = map_ptr->add_edit();
      edit_ptr->set_from_length( edit_ptr->from_length() + 1 );
      edit_ptr->set_to_length( edit_ptr->to_length() + 1 );
    }
    else if ( e == event::mismatch ) {
      edit_ptr = map_ptr->add_edit();
      edit_ptr->set_from_length( 1 );
      edit_ptr->set_to_length( 1 );
      edit_ptr->set_sequence( std::string( 1, last_char ) );
      edit_ptr = nullptr;
    }
    else if ( e == event::deletion ) {
      edit_ptr = map_ptr->add_edit();
      edit_ptr->set_to_length( 1 );
      edit_ptr = nullptr;
    }
    else if ( e == event::insertion ) {
      edit_ptr = map_ptr->add_edit();
      edit_ptr->set_from_length( 1 );
      edit_ptr->set_sequence( std::string( 1, last_char ) );
      edit_ptr = nullptr;
    }
  }

  if ( !fwd ) {
    std::reverse( map_ptr->mutable_edit()->begin(), map_ptr->mutable_edit()->end() );
  }
  return { true, output };
}

template< typename TGraph >
  inline bool
_simulate_read( klibpp::KSeq& read, vg::Path& read_path, klibpp::KSeq const& hapseq,
                Path< TGraph > const& haplotype, std::size_t pos, unsigned int readlen,
                double errorrate, double indelrate, bool fwd, bool allow_ns, bool forced,
                SingleEnd )
{
  /* NOTE: This is assumed that `pos + readlen` won't exceed the sequence size. */
  assert( pos + readlen <= hapseq.seq.size() );

  auto graph_ptr = haplotype.get_graph_ptr();
  if ( !fwd ) pos = hapseq.seq.size() - pos - readlen;

  auto cursor = pos;
  vg::Mapping* map_ptr = nullptr;
  read.seq.resize( readlen );
  auto output = read.seq.begin();
  auto remaining = read.seq.size();
  auto start = hapseq.seq.begin() + pos;
  auto next = start;
  while ( remaining > 0 ) {
    auto id = position_to_id( haplotype, cursor );
    auto offset = position_to_offset( haplotype, cursor );
    auto label_len = graph_ptr->node_length( id );
    auto range = std::min< std::size_t >( label_len - offset, remaining );
    assert( hapseq.seq.end() >= range + start );
    next = start + range;
    map_ptr = read_path.add_mapping();
    vg::Position* p = map_ptr->mutable_position();
    p->set_node_id( graph_ptr->coordinate_id( id ) );
    if ( fwd ) p->set_offset( offset );
    else if ( cursor == pos /* first node when not forward */ ) p->set_offset( 0 );
    else p->set_offset( label_len - range );
    p->set_is_reverse( !fwd );
    bool success;
    std::tie( success, output ) = impose_error( start, next, output, map_ptr, errorrate,
                                                indelrate, fwd, allow_ns, forced );
    if ( !success ) return false;
    start = next;
    cursor += range;
    remaining = read.seq.end() - output;
  }
  if ( !fwd ) {
    read.seq = complement( read.seq );
    std::reverse( read.seq.begin(), read.seq.end() );
    std::reverse( read_path.mutable_mapping()->begin(), read_path.mutable_mapping()->end() );
  }

  auto mappings_ptr = read_path.mutable_mapping();
  std::size_t rank = 0;
  for ( auto it = mappings_ptr->begin(); it != mappings_ptr->end(); ++it ) {
    it->set_rank( ++rank );
  }

  read.name = random::random_string( READ_NAME_LENGTH );
  read.comment = ( hapseq.name + READ_COMMENT_DELIMITER +
                   std::to_string( pos ) + READ_COMMENT_DELIMITER +
                   ( fwd ? "F" : "R" ) );
  read.qual = std::string( read.seq.size(), DEFAULT_QUAL_SCORE );
  return true;
}

template< typename TGraph >
inline bool
_simulate_read( klibpp::KSeq& read, vg::Path& read_path, klibpp::KSeq const& hapseq,
                Path< TGraph > const& haplotype, std::size_t pos, unsigned int readlen,
                double errorrate, double indelrate, bool fwd, bool allow_ns, bool forced,
                PairedEnd )
{
  return false;
}

template< typename TReadType, typename TGraph, typename TCallback >
inline void
simulate_reads( std::string const& name, Path< TGraph > const& haplotype,
                unsigned int readlen, unsigned long int n_reads, double errorrate,
                double indelrate, bool fwd, bool allow_ns, TCallback callback )
{
  if ( n_reads == 0 ) return;

  auto hapseq = to< klibpp::KSeq >( haplotype, name );
  assert( hapseq.seq.size() >= read_span< TReadType >( readlen ) );

  std::size_t ubound = hapseq.seq.size() - read_span< TReadType >( readlen );
  std::uniform_int_distribution<> dis( 0, ubound );
  klibpp::KSeq read;
  vg::Path read_path;
  bool dir = true;
  for ( unsigned long int i = 0; i < n_reads; ++i ) {
    int tries = MAX_TRIES;
    std::size_t pos;
    bool success = false;
    do {
      pos = dis( ::rnd::get_gen() );
      success = _simulate_read( read, read_path, hapseq, haplotype, pos, readlen,
                                errorrate, indelrate, fwd || dir, allow_ns,
                                !static_cast< bool >( tries ), TReadType() );
    } while ( tries-- > 0 && !success );
    if ( tries < 0 ) std::cerr << "Reads may contain 'N' since nothing found after "
                               << MAX_TRIES << " attemps!" << std::endl;
    callback( std::move( read ), std::move( read_path ) );
    read.clear();
    read_path.Clear();
    dir = !dir;
  }
}

template< typename TReadType, typename TType, typename TGraph >
inline void
_simulate( TGraph const& graph, Parameters const& params )
{
  typedef Writer< TType > writer_type;
  typedef typename writer_type::value_type value_type;

  ::rnd::init_gen( params.seed );

  writer_type writer( params.output );
  std::unordered_map< std::string, Path< TGraph > > haplotypes;

  // Simulate haplotypes
  haplotypes.reserve( params.ploidy * graph.get_path_count() /* == no. of regions */ );
  simulate_haplotypes< TGraph >(
      graph, params.ploidy, params.seed,
      [&haplotypes]( std::string name, Path< TGraph > path ) -> void {
        [[maybe_unused]] auto res = haplotypes.insert( { name, std::move( path ) } );
        assert( res.second );
      } );

  if ( params.numreads == 0 ) {  // Output haplotypes
    for ( auto const& h : haplotypes ) writer << to< value_type >( h.second, h.first );
  }
  else {  // Output reads
    // Compute reads distribution across haplotypes
    auto dist = reads_dist< TReadType >( haplotypes, params.readlen, params.numreads );
    assert( dist.size() == haplotypes.size() );
    // Simulate the reads
    auto h_it = haplotypes.begin();
    auto d_it = dist.begin();
    for ( ; h_it != haplotypes.end() /* && d_it != dist.end() */; ++h_it, ++d_it ) {
      simulate_reads< TReadType >(
          h_it->first, h_it->second, params.readlen, *d_it, params.errorrate,
          params.indelrate, params.fwd, params.allow_ns,
          [&writer]( klibpp::KSeq record, vg::Path path ) {
            writer << to< value_type >( std::move( record ), std::move( path ) );
          }
        );
    }
  }
}

template< typename TReadType, typename ...TArgs >
inline void
simulate( fmt::Type const& type, TArgs&&... args )
{
  if ( type == fmt::Seq() ) {
    _simulate< TReadType, fmt::Seq >( std::forward< TArgs >( args )... );
  }
  else if ( type == fmt::Fasta() ) {
    _simulate< TReadType, fmt::Fasta >( std::forward< TArgs >( args )... );
  }
  else if ( type == fmt::Fastq() ) {
    _simulate< TReadType, fmt::Fastq >( std::forward< TArgs >( args )... );
  }
  else if ( type == fmt::Gam() ) {
    _simulate< TReadType, fmt::Gam >( std::forward< TArgs >( args )... );
  }
  else assert( false );
}

  void
config_parser( cxxopts::Options& options )
{
  options.positional_help( "GRAPH" );
  options.add_options()
      ( "o, output", "Write to this file instead of standard output",
        cxxopts::value< std::string >()->default_value( DEFAULT_OUTPUT ) )
      ( "t, type", "Output type (inferred from file extension if not provided). "
        "Values are: 'plain', 'gam', 'fastq', or 'fasta'.",
        cxxopts::value< fmt::Type >() )
      ( "p, ploidy", "Set the ploidy",
        cxxopts::value< unsigned int >()->default_value( DEFAULT_PLOIDY ) )
      ( "l, read-length", "Read length",
        cxxopts::value< unsigned int >()->default_value( DEFAULT_READLEN ) )
      ( "n, num-reads", "Number of reads",
        cxxopts::value< unsigned long int >()->default_value( DEFAULT_NUMREADS ) )
      ( "e, error-rate", "Base error rate",
        cxxopts::value< double >()->default_value( DEFAULT_ERRRATE ) )
      ( "i, indel-rate", "Fraction of indels",
        cxxopts::value< double >()->default_value( DEFAULT_INDRATE ) )
      ( "d, distance", "Outer distance between the two ends (implies paired-end reads)",
        cxxopts::value< unsigned int >()->default_value( DEFAULT_DISTANCE ) )
      ( "s, std-deviation", "Standard deviation (in paired-end reads)",
        cxxopts::value< unsigned int >()->default_value( DEFAULT_DEVIATION ) )
      ( "S, random-seed", "Seed for random generator",
        cxxopts::value< unsigned int >()->default_value( DEFAULT_RNDSEED ) )
      ( "f, forward-only", "Simulate reads only from forward strand",
        cxxopts::value< bool >()->default_value( DEFAULT_FORWARD ) )
      ( "N, allow-Ns", "Allow reads to be sampled from the graph with Ns in them",
        cxxopts::value< bool >()->default_value( DEFAULT_ALLOWNS ) )
      ( "h, help", "Print this message and exit" )
      ;
  options.add_options( "positional" )
      ( "graph", "graph file (vg or gfa)", cxxopts::value< std::string >() )
      ;
  options.parse_positional( { "graph" } );
}

  cxxopts::ParseResult
parse_opts( cxxopts::Options& options, int& argc, char**& argv )
{
  auto result = options.parse( argc, argv );
  // help message
  if ( result.count( "help" ) ) {
    std::cout << options.help( { "" } ) << std::endl;
    throw EXIT_SUCCESS;
  }
  // `graph` option (positional)
  if ( ! result.count( "graph" ) ) {
    throw cxxopts::OptionParseException( "Graph file must be specified" );
  }
  if ( ! readable( result[ "graph" ].as< std::string >() ) ) {
    throw cxxopts::OptionParseException( "Graph file not found" );
  }
  // `type` option
  if ( ! result.count( "type" ) &&
      result[ "output" ].as< std::string >() == DEFAULT_OUTPUT ) {
    throw cxxopts::OptionParseException( "File type must be specified" );
  }
  // `read-length` and `num-reads` options
  bool has_readlen = result[ "read-length" ].as< unsigned int >();
  bool has_numreads = result[ "num-reads" ].as< unsigned long int >();
  if ( has_readlen != has_numreads ) {
    throw cxxopts::OptionParseException(
        "Options `read-length` and `num-reads` should be either both defined "
        "indicating to output simulated reads or not defined at all, in which case it "
        "outputs simulated haplotypes." );
  }

  return result;
}

  int
main( int argc, char* argv[] )
{
  cxxopts::Options options( argv[0], LONG_DESC );
  config_parser( options );

  try {
    auto res = parse_opts( options, argc, argv );

    Parameters params;
    fmt::Type type = fmt::get_type( res );
    std::string graph_path = res[ "graph" ].as< std::string >();

    params.output = res[ "output" ].as< std::string >();
    params.ploidy = res[ "ploidy" ].as< unsigned int >();
    params.readlen = res[ "read-length" ].as< unsigned int >();
    params.numreads = res[ "num-reads" ].as< unsigned long int >();
    params.errorrate = res[ "error-rate" ].as< double >();
    params.indelrate = res[ "indel-rate" ].as< double >();
    params.distance = res[ "distance" ].as< unsigned int >();
    params.sd = res[ "std-deviation" ].as< unsigned int >();
    params.seed = res[ "random-seed" ].as< unsigned int >();
    params.fwd = res[ "forward-only" ].as< bool >();
    params.allow_ns = res[ "allow-Ns" ].as< bool >();

    auto parse_vg = []( std::istream& in ) -> vg::Graph {
      vg::Graph merged;
      std::function< void( vg::Graph& ) > handle_chunks =
        [&]( vg::Graph& other ) {
          gum::util::merge_vg( merged, static_cast< vg::Graph const& >( other ) );
        };
      stream::for_each( in, handle_chunks );
      return merged;
    };

    gum::SeqGraph< gum::Succinct > graph;
    gum::ExternalLoader< vg::Graph > loader{ parse_vg };
    gum::util::load( graph, graph_path, loader, true );
    std::string sort_status = gum::util::ids_in_topological_order( graph ) ? "" : "not ";
    std::cout << "Input graph node IDs are " << sort_status << "in topological sort order."
              << std::endl;

    if ( params.distance != 0 ) simulate< PairedEnd >( type, graph, params );
    else simulate< SingleEnd >( type, graph, params );
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

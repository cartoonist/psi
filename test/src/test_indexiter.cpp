/**
 *    @file  test_indexiter.cpp
 *   @brief  Test index iterators.
 *
 *  It includes all tests related to modules in `index_iter.h`.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Fri Mar 17, 2017  11:37
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <string>
#include <functional>

#include <psi/sequence.hpp>
#include <psi/seed.hpp>
#include <psi/index.hpp>
#include <psi/index_iter.hpp>

#include "test_base.hpp"


using namespace psi;

TEMPLATE_SCENARIO_SIG( "Test if index iterator is movable", "[index][iterator]",
                      ( ( typename T, typename U, int V /* dummy type unless compile error */ ), T, U, V ),
                      ( seqan::IndexWotd<>, seqan::Preorder, 0 ),
                      ( seqan::IndexEsa<>, seqan::Preorder, 0 ),
                      ( seqan::IndexWotd<>, seqan::ParentLinks<>, 0 ),
                      ( seqan::IndexEsa<>, seqan::ParentLinks<>, 0 ) )
{
  typedef T indexspec_type;
  typedef U iterspec_type;
  typedef seqan::String< seqan::Dna, seqan::External<> > text_type;
  typedef seqan::Index< text_type, indexspec_type > index_type;
  typedef TFineIndexIter< index_type, iterspec_type > iterator_type;

  std::string textpath = test_data_dir + "/text/sample_long_sequence.txt";
  text_type  text( textpath.c_str(), seqan::OPEN_RDONLY );
  index_type index( text );
  iterator_type iter( index );

  REQUIRE( go_down( iter, 'A' ) );
  REQUIRE( go_down( iter, 'A' ) );
  REQUIRE( representative( iter.get_iter_() ) == "AA" );

  iterator_type iter2 = std::move( iter );
  iter = iterator_type( index );
  REQUIRE( go_down( iter2, 'A' ) );
  REQUIRE( go_down( iter2, 'G' ) );
  REQUIRE( representative( iter2.get_iter_() ) == "AAAG" );

  iterator_type iter3 = std::move( iter2 );
  iter2 = iterator_type( index );
  REQUIRE( go_down( iter3, 'G' ) );
  REQUIRE( go_down( iter3, 'G' ) );
  REQUIRE( representative( iter3.get_iter_() ) == "AAAGGG" );

  iter = std::move( iter3 );
  iter3 = iterator_type( index );
  REQUIRE( go_down( iter, 'G' ) );
  REQUIRE( representative( iter.get_iter_() ) == "AAAGGGG" );
}

TEMPLATE_SCENARIO( "Fine top-down index iterator basic functionalities", "[index][iterator]",
                   ( seqan::IndexWotd<> ),
                   ( seqan::IndexEsa<> ) )
{
  GIVEN( "A sample small path and ESA index" )
  {
    typedef TestType TIndexSpec;
    typedef seqan::Dna5QString TString;
    typedef seqan::Index< TString, TIndexSpec > TIndex;
    typedef TFineIndexIter< TIndex, seqan::ParentLinks<> > TIter;

    TString str = "GATAGACTAGCCA";
    TIndex index( str );
    TIter itr( index );

    REQUIRE( go_down( itr, 'A' ) );
    REQUIRE( go_down( itr, 'G' ) );
    REQUIRE( representative( itr.get_iter_() ) == "AG" );
    REQUIRE( go_right( itr ) );
    REQUIRE( representative( itr.get_iter_() ) == "ATAGACTAGCCA" );
    REQUIRE( go_down( itr, 'A' ) );
    REQUIRE( go_up( itr ) );
    REQUIRE( go_up( itr ) );
    REQUIRE( go_down( itr, 'C' ) );
    REQUIRE( go_down( itr, 'T' ) );
    REQUIRE( !go_right( itr ) );
    REQUIRE( go_up( itr ) );
    REQUIRE( go_right( itr ) );
    REQUIRE( representative( itr.get_iter_() ) == "AG" );
  }

  GIVEN( "A sample small path and FM-index" )
  {
    typedef psi::MemString TString;
    typedef psi::FMIndex<> TIndexSpec;
    typedef seqan::Index< TString, TIndexSpec > TIndex;
    typedef TFineIndexIter< TIndex, seqan::ParentLinks<> > TIter;

    TString str = "ACCGATCAGATAG";
    TIndex index( str );
    indexRequire( index, seqan::FibreSALF() );
    TIter itr( index );

    REQUIRE( go_down( itr, 'A' ) );
    REQUIRE( go_down( itr, 'G' ) );
    REQUIRE( representative( itr.get_iter_() ) == "GA" );
    REQUIRE( go_right( itr ) );
    REQUIRE( representative( itr.get_iter_() ) == "TA" );
    REQUIRE( go_down( itr, 'A' ) );
    REQUIRE( go_up( itr ) );
    REQUIRE( go_up( itr ) );
    REQUIRE( go_down( itr, 'C' ) );
    REQUIRE( go_down( itr, 'T' ) );
    REQUIRE( !go_right( itr ) );
    REQUIRE( go_up( itr ) );
    REQUIRE( go_right( itr ) );
    REQUIRE( representative( itr.get_iter_() ) == "GA" );
  }
}

TEMPLATE_SCENARIO( "Find k-mer exact matches between two texts using top-down index iterators", "[index][iterator]",
                   ( seqan::IndexWotd<> ),
                   ( seqan::IndexEsa<> ) )
{
  typedef TestType TIndexSpec;
  std::vector< Seed<> > seeds;
  std::function< void(const Seed<>&) > callback = [&seeds]( const Seed<>& hit )
  {
    seeds.push_back( hit );
  };

  GIVEN( "Given two small string set -- set 1" )
  {
    Records< Dna5QStringSet<> > rec1;
    appendValue(rec1.str, "GATAGACTAGCCA");
    appendValue(rec1.str, "GGGCGTAGCCA");
    Records< Dna5QStringSet<> > rec2;
    appendValue(rec2.str, "GGGCGTAGCCA");

    using TIndex1 = seqan::Index< Dna5QStringSet<>, TIndexSpec >;
    using TIndex2 = seqan::Index< Dna5QStringSet<>, TIndexSpec >;

    TIndex1 index1(rec1.str);
    TIndex2 index2(rec2.str);
    TFineIndexIter< TIndex1, seqan::ParentLinks<> > itr1( index1 );
    TFineIndexIter< TIndex2, seqan::ParentLinks<> > itr2( index2 );

    seeds.clear();

    WHEN( "Query all k-mers in one to index of the other" )
    {
      typedef typename seqan::Iterator< decltype( rec1 ), GreedyOverlapping >::Type TSeedsIterator1;
      typedef typename seqan::Iterator< decltype( rec2 ), GreedyOverlapping >::Type TSeedsIterator2;

      TSeedsIterator1 seeds_itr1( &rec1, 4 );
      kmer_exact_matches( index2, &rec2.str, seeds_itr1, callback );
      auto seeds1 = seeds;
      seeds.clear();
      TSeedsIterator2 seeds_itr2( &rec2, 4 );
      kmer_exact_matches( index1, &rec1.str, seeds_itr2, callback );
      auto seeds2 = seeds;
      seeds.clear();

      AND_WHEN( "Traverse both index suffix tree to find all k-mers" )
      {
        kmer_exact_matches( index1, index2, &rec1.str, &rec2, 4, callback );
        auto seeds3 = seeds;
        seeds.clear();
        kmer_exact_matches( itr1, itr2, &rec1.str, &rec2, 4, callback );
        THEN( "Both found all 4-mer exact matches" )
        {
          REQUIRE( seeds.size() == 11 );
          REQUIRE( seeds1.size() == 11 );
          REQUIRE( seeds2.size() == 11 );
          REQUIRE( seeds3.size() == 11 );
        }
      }
    }
  }

  GIVEN( "Given two small string set -- set 2" )
  {
    Records< Dna5QStringSet<> > rec1;
    appendValue(rec1.str, "CATATA");
    Records< Dna5QStringSet<> > rec2;
    appendValue(rec2.str, "ATATAC");

    using TIndex1 = seqan::Index< Dna5QStringSet<>, TIndexSpec >;
    using TIndex2 = seqan::Index< Dna5QStringSet<>, TIndexSpec >;

    TIndex1 index1(rec1.str);
    TIndex2 index2(rec2.str);
    TFineIndexIter< TIndex1, seqan::ParentLinks<> > itr1( index1 );
    TFineIndexIter< TIndex2, seqan::ParentLinks<> > itr2( index2 );

    seeds.clear();

    WHEN( "Query all k-mers in one to index of the other" )
    {
      typedef typename seqan::Iterator< decltype( rec1 ), GreedyOverlapping >::Type TSeedsIterator1;
      typedef typename seqan::Iterator< decltype( rec2 ), GreedyOverlapping >::Type TSeedsIterator2;

      TSeedsIterator1 seeds_itr1( &rec1, 3 );
      kmer_exact_matches( index2, &rec2.str, seeds_itr1, callback );
      auto seeds1 = seeds;
      seeds.clear();
      TSeedsIterator2 seeds_itr2( &rec2, 3 );
      kmer_exact_matches( index1, &rec1.str, seeds_itr2, callback );
      auto seeds2 = seeds;
      seeds.clear();

      AND_WHEN( "Traverse both index suffix tree to find all k-mers" )
      {
        THEN( "Both found all 3-mer exact matches" )
        {
          kmer_exact_matches( index1, index2, &rec1.str, &rec2, 3, callback );
          auto seeds3 = seeds;
          seeds.clear();
          kmer_exact_matches( itr1, itr2, &rec1.str, &rec2, 3, callback );
          REQUIRE( seeds.size() == 5 );
          REQUIRE( seeds1.size() == 5 );
          REQUIRE( seeds2.size() == 5 );
          REQUIRE( seeds3.size() == 5 );
        }
      }
    }
  }

  GIVEN( "Given two small string set -- set 3" )
  {
    Records< Dna5QStringSet<> > rec1;
    appendValue(rec1.str, "TAGGCTACCGATTTAAATAGGCACAC");
    appendValue(rec1.str, "TAGGCTACGGATTTAAATCGGCACAC");
    Records< Dna5QStringSet<> > rec2;
    appendValue(rec2.str, "GGATTTAAATA");
    appendValue(rec2.str, "CGATTTAAATC");
    appendValue(rec2.str, "GGATTTAAATC");
    appendValue(rec2.str, "CGATTTAAATA");

    using TIndex1 = seqan::Index< Dna5QStringSet<>, TIndexSpec >;
    using TIndex2 = seqan::Index< Dna5QStringSet<>, TIndexSpec >;

    TIndex1 index1(rec1.str);
    TIndex2 index2(rec2.str);
    TFineIndexIter< TIndex1, seqan::ParentLinks<> > itr1( index1 );
    TFineIndexIter< TIndex2, seqan::ParentLinks<> > itr2( index2 );

    seeds.clear();

    WHEN( "Query all k-mers in one to index of the other" )
    {
      typedef typename seqan::Iterator< decltype( rec1 ), GreedyOverlapping >::Type TSeedsIterator1;
      typedef typename seqan::Iterator< decltype( rec2 ), GreedyOverlapping >::Type TSeedsIterator2;

      TSeedsIterator1 seeds_itr1( &rec1, 10 );
      kmer_exact_matches( index2, &rec2.str, seeds_itr1, callback );
      auto seeds1 = seeds;
      seeds.clear();
      TSeedsIterator2 seeds_itr2( &rec2, 10 );
      kmer_exact_matches( index1, &rec1.str, seeds_itr2, callback );
      auto seeds2 = seeds;
      seeds.clear();

      AND_WHEN( "Traverse both index suffix tree to find all k-mers" )
      {
        THEN( "Both found all 10-mer exact matches" )
        {
          kmer_exact_matches( index1, index2, &rec1.str, &rec2, 10, callback );
          auto seeds3 = seeds;
          seeds.clear();
          kmer_exact_matches( itr1, itr2, &rec1.str, &rec2, 10, callback );
          REQUIRE( seeds.size() == 8 );
          REQUIRE( seeds1.size() == 8 );
          REQUIRE( seeds2.size() == 8 );
          REQUIRE( seeds3.size() == 8 );
        }
      }
    }
  }

  GIVEN( "Given two small string set with Ns -- set 4" )
  {
    Records< Dna5QStringSet<> > rec1;
    appendValue(rec1.str, "TAGGCTACCGATTNAAATAGGCACAC");
    appendValue(rec1.str, "TAGGCTACGGATTNAAATCGGCACAC");
    Records< Dna5QStringSet<> > rec2;
    appendValue(rec2.str, "GGATTNAAATA");
    appendValue(rec2.str, "CGATTNAAATC");
    appendValue(rec2.str, "GGATTNAAATC");
    appendValue(rec2.str, "CGATTNAAATA");

    using TIndex1 = seqan::Index< Dna5QStringSet<>, TIndexSpec >;
    using TIndex2 = seqan::Index< Dna5QStringSet<>, TIndexSpec >;

    TIndex1 index1(rec1.str);
    TIndex2 index2(rec2.str);
    TFineIndexIter< TIndex1, seqan::ParentLinks<> > itr1( index1 );
    TFineIndexIter< TIndex2, seqan::ParentLinks<> > itr2( index2 );

    seeds.clear();

    WHEN( "Query all k-mers in one to index of the other" )
    {
      typedef typename seqan::Iterator< decltype( rec1 ), GreedyOverlapping >::Type TSeedsIterator1;
      typedef typename seqan::Iterator< decltype( rec2 ), GreedyOverlapping >::Type TSeedsIterator2;

      TSeedsIterator1 seeds_itr1( &rec1, 10 );
      kmer_exact_matches( index2, &rec2.str, seeds_itr1, callback );
      auto seeds1 = seeds;
      seeds.clear();
      TSeedsIterator2 seeds_itr2( &rec2, 10 );
      kmer_exact_matches( index1, &rec1.str, seeds_itr2, callback );
      auto seeds2 = seeds;
      seeds.clear();

      AND_WHEN( "Traverse both index suffix tree to find all k-mers" )
      {
        THEN( "Both found all 10-mer exact matches" )
        {
          kmer_exact_matches( index1, index2, &rec1.str, &rec2, 10, callback );
          auto seeds3 = seeds;
          seeds.clear();
          kmer_exact_matches( itr1, itr2, &rec1.str, &rec2, 10, callback );

          REQUIRE( seeds.size() == 0 );
          REQUIRE( seeds1.size() == 8 );
          REQUIRE( seeds2.size() == 8 );
          REQUIRE( seeds3.size() == 8 );
        }
      }
    }
  }

  GIVEN( "Given two small string set -- set 5" )
  {
    Records< Dna5QStringSet<> > rec1;
    appendValue(rec1.str, "TGCAGTATAGTCGTCGCACGCCTTCTGGCCGCTGGCGGCAGTACAGGATCCTCTTGCTCACAGT"
        "GTAGGGCCCTCTTGCTCCCGGTGTGACGGCTGGCGTGCAGCTGGCTCCCCCGCTGGCAGCTGGGGACACTGACGGGCCC"
        "TCTTGCTCCCCTACTGGCCGCCTCCTGCACCAATTAAAGTCGGAGCACCGGTTACGC");
    appendValue(rec1.str, "TGCAGTATAGTCGTCGCACGCCTTCTGGCCGCTGGCGGCAGTACAGGATCCTCTTGCTCACAGT"
        "GTAGGGCCCTCTTGCTCCCGGTGTGACGGCTGGCGTGCAGCTGGCTCCCCCGCTCGCAGGTGGCGACACAAACGGGCCC"
        "TCTTGCTCCCCTACTGGCCGCCTCCTGCACCAATTAAAGTCGGAGCACCGGTTACGC");

    Records< Dna5QStringSet<> > rec2;
    appendValue(rec2.str, "CATTGCAGAGCCCTCTTGCTCACAGTGTAGTGGCAGCACGCCCGCCTCCTGGCAGCTAGGGACA"
        "GTGCCAGGCCCTCTTGCTCCAAGTGTAGTGGCAGCTGGCTCCCCCGCTGGCAGCTGGGGACACTGACGGGCCCTCTTGC"
        "TTGCAGT");
    appendValue(rec2.str, "TAGGGCAACTGCAGGGCTATCTTGCTTACAGTGGTGTCCAGCGCCCTCTGCTGGCGTCGGAGCA"
        "TTGCAGGGCTCTCTTGCTCGCAGTGTAGTGGCGGCACGCCGCCTGCTGGCAGCTAGGGACATTGCAGAGCCCTCTTGCT"
        "CACAGTG");

    using TIndex1 = seqan::Index< Dna5QStringSet<>, TIndexSpec >;
    using TIndex2 = seqan::Index< Dna5QStringSet<>, TIndexSpec >;

    TIndex1 index1(rec1.str);
    TIndex2 index2(rec2.str);
    TFineIndexIter< TIndex1, seqan::ParentLinks<> > itr1( index1 );
    TFineIndexIter< TIndex2, seqan::ParentLinks<> > itr2( index2 );

    seeds.clear();

    WHEN( "Query all k-mers in one to index of the other" )
    {
      typedef typename seqan::Iterator< decltype( rec1 ), GreedyOverlapping >::Type TSeedsIterator1;
      typedef typename seqan::Iterator< decltype( rec2 ), GreedyOverlapping >::Type TSeedsIterator2;

      TSeedsIterator1 seeds_itr1( &rec1, 30 );
      kmer_exact_matches( index2, &rec2.str, seeds_itr1, callback );
      auto seeds1 = seeds;
      seeds.clear();
      TSeedsIterator2 seeds_itr2( &rec2, 30 );
      kmer_exact_matches( index1, &rec1.str, seeds_itr2, callback );
      auto seeds2 = seeds;
      seeds.clear();

      AND_WHEN( "Traverse both index suffix tree to find all k-mers" )
      {
        THEN( "Both found all 30-mer exact matches" )
        {
          kmer_exact_matches( index1, index2, &rec1.str, &rec2, 30, callback );
          auto seeds3 = seeds;
          seeds.clear();
          kmer_exact_matches( itr1, itr2, &rec1.str, &rec2, 30, callback );
          REQUIRE( seeds.size() == 21 );
          REQUIRE( seeds1.size() == 21 );
          REQUIRE( seeds2.size() == 21 );
          REQUIRE( seeds3.size() == 21 );
        }
      }
    }
  }
}

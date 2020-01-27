/**
 *    @file  tests_indexiter.cc
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

#include "tests_base.h"
#include "sequence.h"
#include "seed.h"
#include "index.h"
#include "index_iter.h"
#include "logger.h"


using namespace grem;

SCENARIO( "Fine top-down index iterator basic functionalities", "[index][iterator]" )
{
  GIVEN( "A sample small path" )
  {
    seqan::Dna5QString str = "GATAGACTAGCCA";
    seqan::Index< seqan::Dna5QString, seqan::IndexEsa<> > index(str);
    TFineIndexIter< seqan::Index< seqan::Dna5QString, seqan::IndexEsa<> >, seqan::ParentLinks<> > itr(index);

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
}

SCENARIO( "Find k-mer exact matches between two texts using top-down index iterators", "[index][iterator]" )
{
  std::vector< Seed<> > seeds;
  std::function< void(const Seed<>&) > callback = [&seeds]( const Seed<>& hit )
  {
    seeds.push_back( hit );
  };

  GIVEN( "Given two small string set -- set 1" )
  {
    Dna5QStringSet<> str1;
    appendValue(str1, "GATAGACTAGCCA");
    appendValue(str1, "GGGCGTAGCCA");
    Dna5QStringSet<> str2;
    appendValue(str2, "GGGCGTAGCCA");

    using TIndex1 = seqan::Index< Dna5QStringSet<>, seqan::IndexEsa<> >;
    using TIndex2 = seqan::Index< Dna5QStringSet<>, seqan::IndexEsa<> >;

    TIndex1 index1(str1);
    TIndex2 index2(str2);
    TFineIndexIter< TIndex1, seqan::ParentLinks<> > itr1( index1 );
    TFineIndexIter< TIndex2, seqan::ParentLinks<> > itr2( index2 );

    seeds.clear();

    WHEN( "Query all k-mers in one to index of the other" )
    {
      kmer_exact_matches( index2, str1, 4, 1, callback );
      auto seeds1 = seeds;
      seeds.clear();
      kmer_exact_matches( index1, str2, 4, 1, callback );
      auto seeds2 = seeds;
      seeds.clear();

      AND_WHEN( "Traverse both index suffix tree to find all k-mers" )
      {
        kmer_exact_matches( index1, index2, &str2, 4, callback );
        auto seeds3 = seeds;
        seeds.clear();
        kmer_exact_matches( itr1, itr2, &str2, 4, callback );
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
    Dna5QStringSet<> str1;
    appendValue(str1, "CATATA");
    Dna5QStringSet<> str2;
    appendValue(str2, "ATATAC");

    using TIndex1 = seqan::Index< Dna5QStringSet<>, seqan::IndexEsa<> >;
    using TIndex2 = seqan::Index< Dna5QStringSet<>, seqan::IndexEsa<> >;

    TIndex1 index1(str1);
    TIndex2 index2(str2);
    TFineIndexIter< TIndex1, seqan::ParentLinks<> > itr1( index1 );
    TFineIndexIter< TIndex2, seqan::ParentLinks<> > itr2( index2 );

    seeds.clear();

    WHEN( "Query all k-mers in one to index of the other" )
    {
      kmer_exact_matches( index2, str1, 3, 1, callback );
      auto seeds1 = seeds;
      seeds.clear();
      kmer_exact_matches( index1, str2, 3, 1, callback );
      auto seeds2 = seeds;
      seeds.clear();

      AND_WHEN( "Traverse both index suffix tree to find all k-mers" )
      {
        THEN( "Both found all 3-mer exact matches" )
        {
          kmer_exact_matches( index1, index2, &str2, 3, callback );
          auto seeds3 = seeds;
          seeds.clear();
          kmer_exact_matches( itr1, itr2, &str2, 3, callback );
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
    Dna5QStringSet<> str1;
    appendValue(str1, "TAGGCTACCGATTTAAATAGGCACAC");
    appendValue(str1, "TAGGCTACGGATTTAAATCGGCACAC");
    Dna5QStringSet<> str2;
    appendValue(str2, "GGATTTAAATA");
    appendValue(str2, "CGATTTAAATC");
    appendValue(str2, "GGATTTAAATC");
    appendValue(str2, "CGATTTAAATA");

    using TIndex1 = seqan::Index< Dna5QStringSet<>, seqan::IndexEsa<> >;
    using TIndex2 = seqan::Index< Dna5QStringSet<>, seqan::IndexEsa<> >;

    TIndex1 index1(str1);
    TIndex2 index2(str2);
    TFineIndexIter< TIndex1, seqan::ParentLinks<> > itr1( index1 );
    TFineIndexIter< TIndex2, seqan::ParentLinks<> > itr2( index2 );

    seeds.clear();

    WHEN( "Query all k-mers in one to index of the other" )
    {
      kmer_exact_matches( index2, str1, 10, 1, callback );
      auto seeds1 = seeds;
      seeds.clear();
      kmer_exact_matches( index1, str2, 10, 1, callback );
      auto seeds2 = seeds;
      seeds.clear();

      AND_WHEN( "Traverse both index suffix tree to find all k-mers" )
      {
        THEN( "Both found all 10-mer exact matches" )
        {
          kmer_exact_matches( index1, index2, &str2, 10, callback );
          auto seeds3 = seeds;
          seeds.clear();
          kmer_exact_matches( itr1, itr2, &str2, 10, callback );
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
    Dna5QStringSet<> str1;
    appendValue(str1, "TAGGCTACCGATTNAAATAGGCACAC");
    appendValue(str1, "TAGGCTACGGATTNAAATCGGCACAC");
    Dna5QStringSet<> str2;
    appendValue(str2, "GGATTNAAATA");
    appendValue(str2, "CGATTNAAATC");
    appendValue(str2, "GGATTNAAATC");
    appendValue(str2, "CGATTNAAATA");

    using TIndex1 = seqan::Index< Dna5QStringSet<>, seqan::IndexEsa<> >;
    using TIndex2 = seqan::Index< Dna5QStringSet<>, seqan::IndexEsa<> >;

    TIndex1 index1(str1);
    TIndex2 index2(str2);
    TFineIndexIter< TIndex1, seqan::ParentLinks<> > itr1( index1 );
    TFineIndexIter< TIndex2, seqan::ParentLinks<> > itr2( index2 );

    seeds.clear();

    WHEN( "Query all k-mers in one to index of the other" )
    {
      kmer_exact_matches( index2, str1, 10, 1, callback );
      auto seeds1 = seeds;
      seeds.clear();
      kmer_exact_matches( index1, str2, 10, 1, callback );
      auto seeds2 = seeds;
      seeds.clear();

      AND_WHEN( "Traverse both index suffix tree to find all k-mers" )
      {
        THEN( "Both found all 10-mer exact matches" )
        {
          kmer_exact_matches( index1, index2, &str2, 10, callback );
          auto seeds3 = seeds;
          seeds.clear();
          kmer_exact_matches( itr1, itr2, &str2, 10, callback );

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
    Dna5QStringSet<> str1;
    appendValue(str1, "TGCAGTATAGTCGTCGCACGCCTTCTGGCCGCTGGCGGCAGTACAGGATCCTCTTGCTCACAGT"
        "GTAGGGCCCTCTTGCTCCCGGTGTGACGGCTGGCGTGCAGCTGGCTCCCCCGCTGGCAGCTGGGGACACTGACGGGCCC"
        "TCTTGCTCCCCTACTGGCCGCCTCCTGCACCAATTAAAGTCGGAGCACCGGTTACGC");
    appendValue(str1, "TGCAGTATAGTCGTCGCACGCCTTCTGGCCGCTGGCGGCAGTACAGGATCCTCTTGCTCACAGT"
        "GTAGGGCCCTCTTGCTCCCGGTGTGACGGCTGGCGTGCAGCTGGCTCCCCCGCTCGCAGGTGGCGACACAAACGGGCCC"
        "TCTTGCTCCCCTACTGGCCGCCTCCTGCACCAATTAAAGTCGGAGCACCGGTTACGC");

    Dna5QStringSet<> str2;
    appendValue(str2, "CATTGCAGAGCCCTCTTGCTCACAGTGTAGTGGCAGCACGCCCGCCTCCTGGCAGCTAGGGACA"
        "GTGCCAGGCCCTCTTGCTCCAAGTGTAGTGGCAGCTGGCTCCCCCGCTGGCAGCTGGGGACACTGACGGGCCCTCTTGC"
        "TTGCAGT");
    appendValue(str2, "TAGGGCAACTGCAGGGCTATCTTGCTTACAGTGGTGTCCAGCGCCCTCTGCTGGCGTCGGAGCA"
        "TTGCAGGGCTCTCTTGCTCGCAGTGTAGTGGCGGCACGCCGCCTGCTGGCAGCTAGGGACATTGCAGAGCCCTCTTGCT"
        "CACAGTG");

    using TIndex1 = seqan::Index< Dna5QStringSet<>, seqan::IndexEsa<> >;
    using TIndex2 = seqan::Index< Dna5QStringSet<>, seqan::IndexEsa<> >;

    TIndex1 index1(str1);
    TIndex2 index2(str2);
    TFineIndexIter< TIndex1, seqan::ParentLinks<> > itr1( index1 );
    TFineIndexIter< TIndex2, seqan::ParentLinks<> > itr2( index2 );

    seeds.clear();

    WHEN( "Query all k-mers in one to index of the other" )
    {
      kmer_exact_matches( index2, str1, 30, 1, callback );
      auto seeds1 = seeds;
      seeds.clear();
      kmer_exact_matches( index1, str2, 30, 1, callback );
      auto seeds2 = seeds;
      seeds.clear();

      AND_WHEN( "Traverse both index suffix tree to find all k-mers" )
      {
        THEN( "Both found all 30-mer exact matches" )
        {
          kmer_exact_matches( index1, index2, &str2, 30, callback );
          auto seeds3 = seeds;
          seeds.clear();
          kmer_exact_matches( itr1, itr2, &str2, 30, callback );
          REQUIRE( seeds.size() == 21 );
          REQUIRE( seeds1.size() == 21 );
          REQUIRE( seeds2.size() == 21 );
          REQUIRE( seeds3.size() == 21 );
        }
      }
    }
  }
}

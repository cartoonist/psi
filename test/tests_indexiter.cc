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

#include "tests_base.h"
#include "sequence.h"
#include "index.h"
#include "index_iter.h"
#include "logger.h"

INITIALIZE_EASYLOGGINGPP

using namespace grem;


SCENARIO ( "Fine top-down index iterator basic functionalities", "[index][iterator]" )
{
  GIVEN ( "A sample small path" )
  {
    seqan::Dna5QString str = "GATAGACTAGCCA";
    seqan::Index < seqan::Dna5QString, seqan::IndexEsa<> > index (str);
    TFineIndexIter < seqan::Index < seqan::Dna5QString, seqan::IndexEsa<> >, seqan::ParentLinks<> > itr(index);

    REQUIRE ( go_down ( itr, 'A' ) );
    REQUIRE ( go_down ( itr, 'G' ) );
    REQUIRE ( representative ( itr.get_iter_() ) == "AG" );
    REQUIRE ( go_right ( itr ) );
    REQUIRE ( representative ( itr.get_iter_() ) == "ATAGACTAGCCA" );
    REQUIRE ( go_down ( itr, 'A' ) );
    REQUIRE ( go_up ( itr ) );
    REQUIRE ( go_up ( itr ) );
    REQUIRE ( go_down ( itr, 'C' ) );
    REQUIRE ( go_down ( itr, 'T' ) );
    REQUIRE ( !go_right ( itr ) );
    REQUIRE ( go_up ( itr ) );
    REQUIRE ( go_right ( itr ) );
    REQUIRE ( representative ( itr.get_iter_() ) == "AG" );
  }
}

SCENARIO ( "Find k-mer exact matches between two fine top-down index iterators", "[index][iterator]" )
{
  GIVEN ( "Given two fine top-down index iterators of two small string set 1" )
  {
    Dna5QStringSet str1;
    appendValue (str1, "GATAGACTAGCCA");
    appendValue (str1, "GGGCGTAGCCA");
    seqan::Index < Dna5QStringSet, seqan::IndexEsa<> > index1 (str1);
    TFineIndexIter < seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >, seqan::ParentLinks<> > itr1(index1);
    Dna5QStringSet str2;
    appendValue (str2, "GGGCGTAGCCA");
    seqan::Index < Dna5QStringSet, seqan::IndexEsa<> > index2 (str2);
    TFineIndexIter < seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >, seqan::ParentLinks<> > itr2(index2);

    seqan::SeedSet < seqan::Seed < seqan::Simple > > seeds;

    THEN ( "Find 4-mer exact matches" )
    {
      kmer_exact_matches <seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >, seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >> ( seeds, itr1, itr2, 4 );
      REQUIRE ( length ( seeds ) == 11 );
    }
  }

  GIVEN ( "Given two fine top-down index iterators of two small string set 2" )
  {
    Dna5QStringSet str1;
    appendValue (str1, "CATATA");
    seqan::Index < Dna5QStringSet, seqan::IndexEsa<> > index1 (str1);
    TFineIndexIter < seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >, seqan::ParentLinks<> > itr1(index1);
    Dna5QStringSet str2;
    appendValue (str2, "ATATAC");
    seqan::Index < Dna5QStringSet, seqan::IndexEsa<> > index2 (str2);
    TFineIndexIter < seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >, seqan::ParentLinks<> > itr2(index2);

    seqan::SeedSet < seqan::Seed < seqan::Simple > > seeds;

    THEN ( "Find 3-mer exact matches" )
    {
      kmer_exact_matches <seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >, seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >> ( seeds, itr1, itr2, 3 );
      REQUIRE ( length ( seeds ) == 5 );
    }
  }

  GIVEN ( "Given two fine top-down index iterators of two small string set 3" )
  {
    Dna5QStringSet str1;
    appendValue (str1, "TAGGCTACCGATTTAAATAGGCACAC");
    appendValue (str1, "TAGGCTACGGATTTAAATCGGCACAC");
    seqan::Index < Dna5QStringSet, seqan::IndexEsa<> > index1 (str1);
    TFineIndexIter < seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >, seqan::ParentLinks<> > itr1(index1);
    Dna5QStringSet str2;
    appendValue (str2, "GGATTTAAATA");
    appendValue (str2, "CGATTTAAATC");
    appendValue (str2, "GGATTTAAATC");
    appendValue (str2, "CGATTTAAATA");
    seqan::Index < Dna5QStringSet, seqan::IndexEsa<> > index2 (str2);
    TFineIndexIter < seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >, seqan::ParentLinks<> > itr2(index2);

    seqan::SeedSet < seqan::Seed < seqan::Simple > > seeds;

    THEN ( "Find 10-mer exact matches" )
    {
      kmer_exact_matches <seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >, seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >> ( seeds, itr1, itr2, 10 );
      REQUIRE ( length ( seeds ) == 8 );
    }
  }

  GIVEN ( "Given two fine top-down index iterators of two small string set 4 with Ns" )
  {
    Dna5QStringSet str1;
    appendValue (str1, "TAGGCTACCGATTNAAATAGGCACAC");
    appendValue (str1, "TAGGCTACGGATTNAAATCGGCACAC");
    seqan::Index < Dna5QStringSet, seqan::IndexEsa<> > index1 (str1);
    TFineIndexIter < seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >, seqan::ParentLinks<> > itr1(index1);
    Dna5QStringSet str2;
    appendValue (str2, "GGATTNAAATA");
    appendValue (str2, "CGATTNAAATC");
    appendValue (str2, "GGATTNAAATC");
    appendValue (str2, "CGATTNAAATA");
    seqan::Index < Dna5QStringSet, seqan::IndexEsa<> > index2 (str2);
    TFineIndexIter < seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >, seqan::ParentLinks<> > itr2(index2);

    seqan::SeedSet < seqan::Seed < seqan::Simple > > seeds;

    THEN ( "Find 10-mer exact matches" )
    {
      kmer_exact_matches <seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >, seqan::Index < Dna5QStringSet, seqan::IndexEsa<> >> ( seeds, itr1, itr2, 10 );
      REQUIRE ( length ( seeds ) == 0 );
    }
  }

  GIVEN ( "Given two fine top-down index iterators of two small string set 5" )
  {
    Dna5QStringSet str1;
    appendValue (str1, "TGCAGTATAGTCGTCGCACGCCTTCTGGCCGCTGGCGGCAGTACAGGATCCTCTTGCTCACAGT"
        "GTAGGGCCCTCTTGCTCCCGGTGTGACGGCTGGCGTGCAGCTGGCTCCCCCGCTGGCAGCTGGGGACACTGACGGGCCC"
        "TCTTGCTCCCCTACTGGCCGCCTCCTGCACCAATTAAAGTCGGAGCACCGGTTACGC");
    appendValue (str1, "TGCAGTATAGTCGTCGCACGCCTTCTGGCCGCTGGCGGCAGTACAGGATCCTCTTGCTCACAGT"
        "GTAGGGCCCTCTTGCTCCCGGTGTGACGGCTGGCGTGCAGCTGGCTCCCCCGCTCGCAGGTGGCGACACAAACGGGCCC"
        "TCTTGCTCCCCTACTGGCCGCCTCCTGCACCAATTAAAGTCGGAGCACCGGTTACGC");

    Dna5QStringSet str2;
    appendValue (str2, "CATTGCAGAGCCCTCTTGCTCACAGTGTAGTGGCAGCACGCCCGCCTCCTGGCAGCTAGGGACA"
        "GTGCCAGGCCCTCTTGCTCCAAGTGTAGTGGCAGCTGGCTCCCCCGCTGGCAGCTGGGGACACTGACGGGCCCTCTTGC"
        "TTGCAGT");
    appendValue (str2, "TAGGGCAACTGCAGGGCTATCTTGCTTACAGTGGTGTCCAGCGCCCTCTGCTGGCGTCGGAGCA"
        "TTGCAGGGCTCTCTTGCTCGCAGTGTAGTGGCGGCACGCCGCCTGCTGGCAGCTAGGGACATTGCAGAGCCCTCTTGCT"
        "CACAGTG");

    typedef Dna5QStringSetIndex < seqan::IndexEsa<> > TIndexEsa;

    TIndexEsa index1 (str1);
    TIndexEsa index2 (str2);

    TFineIndexIter < TIndexEsa, seqan::ParentLinks<> > itr1(index1);
    TFineIndexIter < TIndexEsa, seqan::ParentLinks<> > itr2(index2);

    seqan::SeedSet < seqan::Seed < seqan::Simple > > seeds;

    THEN ( "Find 10-mer exact matches" )
    {
      kmer_exact_matches < TIndexEsa, TIndexEsa > ( seeds, itr1, itr2, 30 );
      REQUIRE ( length ( seeds ) == 21 );
    }
  }

  GIVEN ( "Given two top-down index iterators of two small string set" )
  {
    Dna5QStringSet str1;
    appendValue (str1, "TAGGCTACCGATTTAAATAGGCACAC");
    appendValue (str1, "TAGGCTACGGATTTAAATCGGCACAC");

    Dna5QStringSet str2;
    appendValue (str2, "GGATTTAAATA");
    appendValue (str2, "CGATTTAAATC");
    appendValue (str2, "GGATTTAAATC");
    appendValue (str2, "CGATTTAAATA");

    typedef Dna5QStringSetIndex < seqan::IndexEsa<> > TIndexEsa;

    TIndexEsa index1 (str1);
    TIndexEsa index2 (str2);

    TIndexIter < TIndexEsa, seqan::TopDown< seqan::ParentLinks<> > > itr1(index1);
    TIndexIter < TIndexEsa, seqan::TopDown< seqan::ParentLinks<> > > itr2(index2);

    seqan::SeedSet < seqan::Seed < seqan::Simple > > seeds;

    THEN ( "Find 10-mer exact matches" )
    {
      kmer_exact_matches < TIndexEsa, TIndexEsa > ( seeds, itr1, itr2, 10 );
      REQUIRE ( length ( seeds ) == 8 );
    }
  }
}

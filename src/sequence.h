/**
 *    @file  sequence.h
 *   @brief  Sequence abstract data types.
 *
 *  Sequence header file contains type definitions, abstract data types, and helper
 *  interface functions (mostly implemented at the top of SeqAn sequence library; i.e.
 *  `seqan/sequence.h`) to work with sequences and sequence sets.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Mar 07, 2017  12:50
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  SEQUENCE_H__
#define  SEQUENCE_H__

#include <seqan/seq_io.h>
#include <seqan/sequence.h>


namespace grem {
  // :TODO:Tue Sep 05 09:36:\@cartoonist: clear the code up from direct usage of seqan::Dna5QString.
  /* Typedefs  ------------------------------------------------------------------- */
  typedef seqan::StringSet< seqan::CharString > CharStringSet;
  typedef seqan::StringSet< seqan::Dna5QString > Dna5QStringSet;
  typedef seqan::Position< Dna5QStringSet >::Type Dna5QStringSetPosition;
  /* END OF Typedefs  ------------------------------------------------------------ */

  /* Data structures  ------------------------------------------------------------ */

  /**
   *  @brief  String set with an ID associated with each string.
   *
   *  It is a wrapper class on StringSet associating an ID to each string.
   */
  template< typename TStringSet >
    class NamedStringSet {
      public:
        /* ====================  DATA MEMBERS  ======================================= */
        TStringSet str;
        CharStringSet id;
    };  /* ----------  end of template class NamedStringSet  ---------- */

  typedef NamedStringSet< Dna5QStringSet > Dna5QRecords;

  /* END OF Data structures  ----------------------------------------------------- */

  /* Interface functions  -------------------------------------------------------- */
  /**
   *  @brief  Read records from the input file into a sequence record set.
   *
   *  @param[out]  records Sequence record set to store records in the input file.
   *  @param[in,out]  infile The input file.
   *  @param[in]  num_record Read this number of record from the input file.
   *
   *  A wrapper function for `seqan::readRecords` method to read the records into
   *  sequence record set. If `num_record` is equal to zero, it reads all recrods.
   */
    inline void
  readRecords( Dna5QRecords& records, seqan::SeqFileIn& infile,
      unsigned int num_record = 0 )
  {
    CharStringSet quals;
    if ( num_record != 0 ) {
      seqan::readRecords( records.id, records.str, quals, infile, num_record );
    }
    else {
      seqan::readRecords( records.id, records.str, quals, infile );
    }
    assignQualities( records.str, quals );
    return;
  }  /* -----  end of function readRecords  ----- */

  /* Seeding strategies */
  struct FixedLengthNonOverlappingSeeding;

  /* Seeding strategy tags */
  typedef seqan::Tag< FixedLengthNonOverlappingSeeding > FixedLengthNonOverlapping;

  /**
   *  @brief  Seeding by partitioning each sequence into non-overlapping k-mers.
   *
   *  @param  string_set The string set from which seeds are extracted.
   *  @param  k The length of the seeds.
   *  @param  tag Tag for fixed-length non-overlapping seeding strategy.
   *  @return A set of strings containing seeds.
   *
   *  Extract a set of non-overlapping seeds of length k.
   *
   *  NOTE: In case that the length of sequence is not dividable by k the last seed
   *        may overlap its previous.
   */
  inline Dna5QStringSet
    seeding( const Dna5QStringSet& string_set, unsigned int k, FixedLengthNonOverlapping )
    {
      Dna5QStringSet seeds;
      reserve ( seeds, static_cast<int>( lengthSum ( string_set ) / k ) );

      for ( unsigned int idx = 0; idx < length ( string_set ); ++idx ) {
        for ( unsigned int i = 0; i < length ( string_set[idx] ) - k; i += k ) {
          appendValue ( seeds, infixWithLength ( string_set[idx], i, k ) );
        }
        unsigned int last = length ( string_set[idx] ) - k;
        appendValue ( seeds, infixWithLength ( string_set[idx], last, k ) );
      }

      return seeds;
    }  /* -----  end of function seeding  ----- */
  /* END OF Interface functions  ------------------------------------------------- */

}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef SEQUENCE_H__  ----- */

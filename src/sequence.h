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
  template< typename TSpec = seqan::Owner<> >
    using CharStringSet = seqan::StringSet< seqan::CharString, TSpec >;
  template< typename TSpec = seqan::Owner<> >
    using Dna5QStringSet = seqan::StringSet< seqan::Dna5QString, TSpec >;
  typedef seqan::Dependent< seqan::Generous > Dependent;
  /* END OF Typedefs  ------------------------------------------------------------ */

  /* Meta-functions  ------------------------------------------------------------- */

  template< typename TStringSet >
    class MakeOwner;

  template< typename TStringSet >
    class MakeDependent;

  template< typename TText >
    class MakeOwner< seqan::StringSet< TText, grem::Dependent > > {
      public:
        typedef seqan::StringSet< TText, seqan::Owner<> > Type;
    };

  template< typename TText >
    class MakeOwner< seqan::StringSet< TText, seqan::Owner<> > > {
      public:
        typedef seqan::StringSet< TText, seqan::Owner<> > Type;
    };

  template< typename TText >
    class MakeDependent< seqan::StringSet< TText, seqan::Owner<> > > {
      public:
        typedef seqan::StringSet< TText, grem::Dependent > Type;
    };

  template< typename TText >
    class MakeDependent< seqan::StringSet< TText, grem::Dependent > > {
      public:
        typedef seqan::StringSet< TText, grem::Dependent > Type;
    };

  /* END OF Meta-functions  ------------------------------------------------------ */

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
        CharStringSet<> name;
    };  /* ----------  end of template class NamedStringSet  ---------- */

  /* Forwards */
  template< typename TStringSet >
    class Records;

  /* Records interface functions */
  template< typename TText >
      inline typename Records< seqan::StringSet< TText, seqan::Owner<> > >::TId
    position_to_id( const Records< seqan::StringSet< TText, seqan::Owner<> > >& records,
        typename Records< seqan::StringSet< TText, seqan::Owner<> > >::TPosition pos )
    {
      if ( pos >= length( records.str ) || pos < 0 ) {
        throw std::runtime_error( "position out of range" );
      }
      return pos;          /**< @brief In Owner records ID and position are identical. */
    }

  template< typename TText >
      inline typename Records< seqan::StringSet< TText, grem::Dependent > >::TId
    position_to_id( const Records< seqan::StringSet< TText, grem::Dependent > >& records,
        typename Records< seqan::StringSet< TText, grem::Dependent > >::TPosition pos )
    {
      if ( pos >= length( records.str ) || pos < 0 ) {
        throw std::runtime_error( "position out of range" );
      }
      assert( records.o_str != nullptr );
      return records.offset + pos;
    }


  template< typename TText >
      inline void
    clear( Records< seqan::StringSet< TText, seqan::Owner<> > >&records )
    {
      clear( records.str );
      clear( records.name );
    }

  template< typename TText >
      inline void
    clear( Records< seqan::StringSet< TText, grem::Dependent > >&records )
    {
      clear( records.str );
      records.offset = 0;
      records.o_str = nullptr;
    }

  template< typename TText, typename TStringSetSpec >
      inline typename Records< seqan::StringSet< TText, TStringSetSpec > >::TSize
    length( const Records< seqan::StringSet< TText, TStringSetSpec > >& records )
    {
      return length( records.str );
    }

  template< typename TText >
      inline bool
    load_chunk( Records< seqan::StringSet< TText, grem::Dependent > >& records,
        const Records< seqan::StringSet< TText, seqan::Owner<> > >& ref,
        typename Records< seqan::StringSet< TText, seqan::Owner<> > >::TPosition n,
        typename Records< seqan::StringSet< TText, seqan::Owner<> > >::TPosition start_pos )
    {
      if ( start_pos >= length( ref.str ) || start_pos < 0 )
      {
        return false;
      }

      clear( records.str );
      records.offset = start_pos;
      records.o_str = &ref.str;
      auto i = start_pos;
      for ( ; i < n + start_pos && i < length( ref.str ); ++i ) {
        appendValue( records.str, ref.str[ i ] );
      }
      return true;
    }

  template< typename TText >
      inline bool
    load_chunk( Records< seqan::StringSet< TText, grem::Dependent > >& records,
        const Records< seqan::StringSet< TText, seqan::Owner<> > >& ref,
        typename Records< seqan::StringSet< TText, seqan::Owner<> > >::TPosition n )
    {
      if ( length( records.str ) == 0 ) {
        return load_chunk( records, ref, n, 0 );
      }
      return load_chunk( records, ref, n, records.offset + n );
    }

  template< typename TText >
    class Records< seqan::StringSet< TText, seqan::Owner<> > >
    : public NamedStringSet< seqan::StringSet< TText, seqan::Owner<> > > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef seqan::StringSet< TText, seqan::Owner<> > TStringSet;
        typedef typename MakeOwner< TStringSet >::Type TRefStringSet;
        typedef typename seqan::Position< TStringSet >::Type TPosition;
        typedef typename seqan::Id< TStringSet >::Type TId;
        typedef typename seqan::Size< TStringSet >::Type TSize;
        typedef typename seqan::Value< TStringSet >::Type TValue;
    };

  template< typename TText >
    class Records< seqan::StringSet< TText, grem::Dependent > > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef seqan::StringSet< TText, grem::Dependent > TStringSet;
        typedef typename MakeOwner< TStringSet >::Type TRefStringSet;
        typedef typename seqan::Position< TStringSet >::Type TPosition;
        typedef typename seqan::Id< TStringSet >::Type TId;
        typedef typename seqan::Size< TStringSet >::Type TSize;
        typedef typename seqan::Value< TStringSet >::Type TValue;
        /* ====================  DATA MEMBERS  ======================================= */
        TStringSet str;
        /* ====================  LIFECYCLE     ======================================= */
        Records( ) : offset( 0 ), o_str( nullptr ) { }
        /* ====================  INTERFACE FUNCTIONS  ================================ */
          friend TId
        position_to_id< TText >( const Records& records, TPosition pos );

          friend void
        clear< TText >( Records& records );

          friend bool
        load_chunk< TText >( Records& records,
            const Records< TRefStringSet >& ref,
            typename Records< TRefStringSet >::TPosition n,
            typename Records< TRefStringSet >::TPosition start_pos );

          friend bool
        load_chunk< TText >( Records& records,
            const Records< TRefStringSet >& ref,
            typename Records< TRefStringSet >::TPosition n );

      protected:
        /* ====================  DATA MEMBERS  ======================================= */
        std::size_t offset;  /**< @brief First string ID: id(i) = offset + pos(i). */
        const TRefStringSet* o_str;  /**< @brief original string set. */
    };

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
  template< typename TText >
      inline void
    readRecords( Records< seqan::StringSet< TText, seqan::Owner<> > >& records,
        seqan::SeqFileIn& infile,
        unsigned int num_record = 0 )
    {
      CharStringSet<> quals;
      if ( num_record != 0 ) {
        seqan::readRecords( records.name, records.str, quals, infile, num_record );
      }
      else {
        seqan::readRecords( records.name, records.str, quals, infile );
      }
      assignQualities( records.str, quals );
      return;
    }  /* -----  end of template function readRecords  ----- */

  /* Seeding strategies */
  struct GreedyOverlapStrategy;
  struct NonOverlapStrategy;
  struct GreedyNonOverlapStrategy;

  /* Seeding strategy tags */
  typedef seqan::Tag< GreedyOverlapStrategy > GreedyOverlapping;
  typedef seqan::Tag< NonOverlapStrategy > NonOverlapping;
  typedef seqan::Tag< GreedyNonOverlapStrategy > GreedyNonOverlapping;

  /**
   *  @brief  Add any k-mers from the given string set with `step` distance to seed set.
   *
   *  @param  seeds The seed set.
   *  @param  string_set The string set from which seeds are extracted.
   *  @param  k The length of the seeds.
   *  @param  step The step size.
   *
   *  For each string in string set, it add all substring of length `k` starting from 0
   *  to end of string with `step` distance with each other. If `step` is equal to `k`,
   *  it gets non-overlapping substrings of length k.
   */
  template< typename TText, typename TStringSetSpec >
      inline void
    _seeding( seqan::StringSet< TText, seqan::Owner<> >& seeds,
        const seqan::StringSet< TText, TStringSetSpec >& string_set,
        unsigned int k,
        unsigned int step )
    {
      for ( unsigned int idx = 0; idx < length( string_set ); ++idx ) {
        for ( unsigned int i = 0; i < length( string_set[idx] ) - k + 1; i += step ) {
          appendValue( seeds, infixWithLength( string_set[idx], i, k ) );
        }
      }
    }  /* -----  end of template function _seeding  ----- */

  /**
   *  @brief  Seeding a set of sequence by reporting overlapping k-mers.
   *
   *  @param  seeds The resulting set of strings containing seeds.
   *  @param  string_set The string set from which seeds are extracted.
   *  @param  k The length of the seeds.
   *  @param  tag Tag for greedy fixed-length overlapping seeding strategy.
   *
   *  Extract a set of overlapping seeds of length k.
   */
  template< typename TText, typename TStringSetSpec >
      inline void
    seeding( seqan::StringSet< TText, seqan::Owner<> >& seeds,
        const seqan::StringSet< TText, TStringSetSpec >& string_set,
        unsigned int k,
        GreedyOverlapping )
    {
      clear( seeds );
      unsigned int avg_read_len = lengthSum( string_set ) / length( string_set );
      reserve( seeds, static_cast<int>( length( string_set ) * ( avg_read_len - k ) ) );
      _seeding( seeds, string_set, k, 1 );
    }  /* -----  end of template function seeding  ----- */

  /**
   *  @brief  Seeding by partitioning each sequence into non-overlapping k-mers.
   *
   *  @param  seeds The resulting set of strings containing seeds.
   *  @param  string_set The string set from which seeds are extracted.
   *  @param  k The length of the seeds.
   *  @param  tag Tag for fixed-length non-overlapping seeding strategy.
   *
   *  Extract a set of non-overlapping seeds of length k.
   */
  template< typename TText, typename TStringSetSpec >
      inline void
    seeding( seqan::StringSet< TText, seqan::Owner<> >& seeds,
        const seqan::StringSet< TText, TStringSetSpec >& string_set,
        unsigned int k,
        NonOverlapping )
    {
      clear( seeds );
      reserve( seeds, static_cast<int>( lengthSum( string_set ) / k ) );
      _seeding( seeds, string_set, k, k );
    }  /* -----  end of function seeding  ----- */

  /**
   *  @brief  Seeding by greedy partitioning each sequence into non-overlapping k-mers.
   *
   *  @param  seeds The resulting set of strings containing seeds.
   *  @param  string_set The string set from which seeds are extracted.
   *  @param  k The length of the seeds.
   *  @param  tag Tag for greedy fixed-length non-overlapping seeding strategy.
   *
   *  Extract a set of non-overlapping seeds of length k.
   *
   *  NOTE: In case that the length of sequence is not dividable by k the last seed
   *        may overlap its previous (greedy).
   */
  template< typename TText, typename TStringSetSpec >
      inline void
    seeding( seqan::StringSet< TText, seqan::Owner<> >& seeds,
        const seqan::StringSet< TText, TStringSetSpec >& string_set,
        unsigned int k,
        GreedyNonOverlapping )
    {
      clear( seeds );
      reserve( seeds, static_cast<int>( lengthSum( string_set ) / k ) );

      for ( unsigned int idx = 0; idx < length( string_set ); ++idx ) {
        for ( unsigned int i = 0; i < length( string_set[idx] ) - k; i += k ) {
          appendValue( seeds, infixWithLength( string_set[idx], i, k ) );
        }
        unsigned int last = length( string_set[idx] ) - k;
        appendValue( seeds, infixWithLength( string_set[idx], last, k ) );
      }
    }  /* -----  end of function seeding  ----- */

  /* END OF Interface functions  ------------------------------------------------- */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef SEQUENCE_H__  ----- */

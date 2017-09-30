/**
 *    @file  utils.h
 *   @brief  Utility and helper functions.
 *
 *  This header file contains general utility and helper functions.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Mar 07, 2017  20:11
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  UTILS_H__
#define  UTILS_H__

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include <seqan/seq_io.h>

#include "sequence.h"

namespace grem {
  /**
   *  @brief  Read records from the input file into named string set.
   *
   *  @param[out]  records Named string set to store records in the input file.
   *  @param[in,out]  infile The input file.
   *  @param[in]  num_record Read this number of record from the input file.
   *
   *  A wrapper function for `seqan::readRecords` method to read the records into named
   *  string set.
   */
    inline void
  readRecords( Dna5QRecords& records, seqan::SeqFileIn& infile, unsigned int num_record )
  {
    CharStringSet quals;
    seqan::readRecords( records.id, records.str, quals, infile, num_record );
    assignQualities( records.str, quals );
    return;
  }  /* -----  end of function readRecords  ----- */


  /**
   *  @brief  Check whether a string ends with another string.
   *
   *  @param  str The first string.
   *  @param  suf The second string.
   *  @return `true` if `suf` is a suffix of `str`; otherwise `false`.
   *
   *  It checks the first string whether the second one is one of its suffixes or not.
   */
    inline bool
  ends_with( std::string str, std::string suf )
  {
    if ( suf.length() <= str.length() ) {
      if ( std::equal( suf.rbegin(), suf.rend(), str.rbegin() ) ) {
        return true;
      }
    }
    return false;
  }  /* -----  end of function ends_with  ----- */


  /**
   *  @brief  Check if the given file exists and is readable.
   *
   *  @param  file_name The name of the file to be checked.
   *  @return `true` if exists and is readable; otherwise `false`.
   *
   *  Use `std::ifstream` to check for readability of the given file because of
   *  portability.
   */
    inline bool
  readable( const std::string& file_name )
  {
    std::ifstream ifs( file_name );
    return ifs.good();
  }  /* -----  end of function exists  ----- */


  /**
   *  @brief  Check if the given file exists and is writable.
   *
   *  @param  file_name The name of the file to be checked.
   *  @return `true` if exists and is writable; otherwise `false`.
   *
   *  Use `std::ofstream` to check for writability of the given file because of
   *  portability.
   */
    inline bool
  writable( const std::string& file_name )
  {
    std::ofstream ofs( file_name );
    if ( ofs.good() ) {
      std::remove( file_name.c_str() );
      return true;
    }
    return false;
  }  /* -----  end of function writable  ----- */


  typedef uint64_t TContainerSize;

  /**
   *  @brief  Container serialization implementation.
   *
   *  @param[out]  out Output stream.
   *  @param[in]  size Size of the container.
   *  @param[in]  begin The begin input iterator.
   *  @param[in]  end The end input iterator.
   *
   *  First, the size of the container is written represented by `TSize` type followed
   *  by all items between two input iterators are written.
   */
  template< typename TIter, typename TSize >
    void
  _serialize( std::ostream& out, TSize size, TIter begin, TIter end )
  {
    out.write( reinterpret_cast< char* >( &size ), sizeof( TSize ) );

    for ( ; begin != end; ++begin ) {
      out.write( reinterpret_cast< const char* >( &*begin ),
          sizeof( typename TIter::value_type ));
    }

    out.flush();
  }  /* -----  end of template function _serialize  ----- */


  /**
   *  @brief  Serialize an ordered container to an output stream.
   *
   *  @param[out]  out Output stream.
   *  @param[in]  begin The begin input iterator.
   *  @param[in]  end The end input iterator.
   *
   *  It gets two iterators of a container and serialize it to the given output stream
   *  by writing each item to the stream following the size of the container. The
   *  `TSize` represents the size type.
   */
  template< typename TIter, typename TSize = TContainerSize >
    void
  serialize( std::ostream& out, TIter begin, TIter end )
  {
    TSize size = end - begin;
    _serialize( out, size, begin, end );
  }  /* -----  end of template function serialize  ----- */


  /**
   *  @brief  Serialize an unordered container to an output stream.
   *
   *  @param[out]  out Output stream.
   *  @param[in]  container The container.
   *  @param[in]  begin The begin input iterator.
   *  @param[in]  end The end input iterator.
   *
   *  It gets two iterators of a container and serialize it to the given output stream
   *  by writing each item to the stream following the size of the container. The
   *  `TSize` represents the size type. Since the container is unordered, it needs
   *  container itself to infer the size.
   */
  template< typename TContainer, typename TIter, typename TSize = TContainerSize >
    void
  serialize( std::ostream& out, TContainer container, TIter begin, TIter end )
  {
    TSize size = container.size();
    _serialize( out, size, begin, end );
  }  /* -----  end of template function serialize  ----- */


  /**
   *  @brief  Deserialize a container from an input stream.
   *
   *  @param[in,out]  in The input stream.
   *  @param[out]  container The output container.
   *  @param[out]  itr The output iterator.
   *
   *  It gets an output container, and an output iterator. Then reads from input stream
   *  and populate the container from serialized data.
   */
  template< typename TContainer, typename TOutIter, typename TSize = TContainerSize >
    void
  deserialize ( std::istream& in, TContainer& container, TOutIter itr )
  {
    TSize size;
    in.read( reinterpret_cast< char* >( &size ), sizeof( TSize ) );
    container.reserve( size );
    for ( unsigned int i = 0; i < size; ++i ) {
      typename TContainer::value_type item;
      in.read( reinterpret_cast< char* >( &item ),
          sizeof( typename TContainer::value_type ) );
      *itr++ = item;
    }
  }  /* -----  end of template function deserialize  ----- */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef UTILS_H__  ----- */

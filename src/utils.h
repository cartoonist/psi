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

#include "sequence.h"

namespace grem {
  /**
   *  @brief  Check whether a string ends with another string.
   *
   *  @param  str The first string.
   *  @param  suf The second string.
   *  @return `true` if `suf` is a suffix of `str`; otherwise `false`.
   *
   *  It checks the first string whether the second one is one of its suffixes or not.
   */
  template< typename TText >
      inline bool
    ends_with( const TText& str, const TText& suf )
    {
      typedef typename seqan::Size< TText >::Type TSize;
      bool retval = false;
      if( length( suf ) <= length( str ) ) {
        retval = true;
        for ( TSize i = 1; i <= length( suf ); ++i ) {
          if ( suf[ length(suf) - i ] != str[ length(str) - i ] ) {
            retval = false;
            break;
          }
        }
      }
      return retval;
    }  /* -----  end of template function ends_with  ----- */

  /**
   *  @brief  Check whether a string ends with another string.
   *
   *  @param  str The first string.
   *  @param  suf The second string.
   *  @return `true` if `suf` is a suffix of `str`; otherwise `false`.
   *
   *  It checks the first string whether the second one is one of its suffixes or not.
   *
   *  Specialized for `std::string` and `std::string`.
   *
   *  NOTE: `suf` is passed by VALUE!
   */
    inline bool
  ends_with( const std::string& str, std::string suf )
  {
    if ( suf.length() <= str.length() &&
        std::equal( suf.rbegin(), suf.rend(), str.rbegin() ) ) {
      return true;
    }
    return false;
  }  /* -----  end of function ends_with  ----- */

  /**
   *  @brief  Check whether a string starts with another string.
   *
   *  @param  str The first string.
   *  @param  pre The second string.
   *  @return `true` if `pre` is a prefix of `str`; otherwise `false`.
   *
   *  It checks the first string whether the second one is one of its prefixes or not.
   */
  template< typename TText >
      inline bool
    starts_with( const TText& str, const TText& pre )
    {
      typedef typename seqan::Size< TText >::Type TSize;
      bool retval = false;
      if( length( pre ) <= length( str ) ) {
        retval = true;
        for ( TSize i = 0; i < length( pre ); ++i ) {
          if ( pre[ i ] != str[ i ] ) {
            retval = false;
            break;
          }
        }
      }
      return retval;
    }  /* -----  end of template function starts_with  ----- */

  /**
   *  @brief  Check whether a string starts with another string.
   *
   *  @param  str The first string.
   *  @param  pre The second string.
   *  @return `true` if `pre` is a prefix of `str`; otherwise `false`.
   *
   *  It checks the first string whether the second one is one of its prefixes or not.
   *
   *  Specialized for `std::string` and `std::string`.
   *
   *  NOTE: `pre` is passed by VALUE!
   */
    inline bool
  starts_with( const std::string& str, std::string pre )
  {
    if ( pre.length() <= str.length() &&
        std::equal( pre.begin(), pre.end(), str.begin() ) ) {
      return true;
    }
    return false;
  }  /* -----  end of function starts_with  ----- */


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
   *  @brief  Simple object serialization implementation.
   *
   *  @param[out]  out The output stream.
   *  @param[in]  obj The object to be written into the output stream.
   *
   *  It simply writes the object bitwise into the stream.
   */
  template< typename TObject >
      inline void
    serialize( std::ostream& out, const TObject& obj )
    {
      out.write( reinterpret_cast< const char* >( &obj ), sizeof( TObject ) );
    }


  /**
   *  @brief  Container serialization implementation.
   *
   *  @param[out]  out The output stream.
   *  @param[in]  size The size of the container.
   *  @param[in]  begin The begin input iterator.
   *  @param[in]  end The end input iterator.
   *
   *  First, the size of the container is written represented by `TSize` type followed
   *  by all items between two input iterators are written.
   */
  template< typename TIter, typename TSize >
    inline void
  _serialize( std::ostream& out, TSize size, TIter begin, TIter end )
  {
    serialize( out, size );
    for ( ; begin != end; ++begin ) serialize( out, *begin );
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
    inline void
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
    inline void
  serialize( std::ostream& out, const TContainer& container, TIter begin, TIter end )
  {
    TSize size = container.size();
    _serialize( out, size, begin, end );
  }  /* -----  end of template function serialize  ----- */


  /**
   *  @brief  Deserialize a simple object from an input stream.
   *
   *  @param[in,out]  in The input stream.
   *  @param[out]  obj The deserialized value will be written to this variable.
   *
   *  It simply reads the object from an input stream.
   */
  template< typename TObject >
      inline void
    deserialize( std::istream& in, TObject& obj )
    {
      in.read( reinterpret_cast< char* >( &obj ), sizeof( TObject ) );
    }


  /**
   *  @brief  Reserve the required memory for the container of size `size`.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   *  It calls reserve member function of the container.
   */
  template< typename TContainer, typename TSize >
      inline void
    _reserve( TContainer& container, TSize size )
    {
      container.reserve( size );
    }


  /**
   *  @overload The `std::deque` cannot be reserved.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   *  Do nothing.
   */
  template< typename TObject, typename TSize >
      inline void
    _reserve( std::deque< TObject >& container, TSize size )
    {
      /* NOOP */
    }


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
    inline void
  deserialize( std::istream& in, TContainer& container, TOutIter itr )
  {
    TSize size;
    deserialize( in, size );
    _reserve( container, size );
    for ( unsigned int i = 0; i < size; ++i ) {
      typename TContainer::value_type item;
      deserialize( in, item );
      *itr++ = item;
    }
  }  /* -----  end of template function deserialize  ----- */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef UTILS_H__  ----- */

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
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <memory>

#include <seqan/basic.h>
#include <sdsl/enc_vector.hpp>

#include "base.h"


#define GREM_DEFAULT_TMPDIR "/tmp"
#define GREM_TMPFILE_TEMPLATE "/grem-XXXXXX"

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
   *  @brief  Copy all bits in `src[idx...idx+len)` to `dst[idx...idx+len)`.
   *
   *  @param  src The source bit vector.
   *  @param  dst The destination bit vector.
   *  @param  idx The start index in `src` to copy to the "identical" index in `dst`.
   *
   *  Bitvector identical copy. The `src[idx...idx+len)` is copied to the same index of
   *  `dst`.
   */
  template< typename TBitVector >
      inline void
    bv_icopy( const TBitVector& src, TBitVector& dst, typename TBitVector::size_type idx=0,
        typename TBitVector::size_type len=0 )
    {
      static const short int WLEN = 64;

      assert( idx < src.size() );
      assert( dst.size() >= src.size() );

      if ( len == 0 ) len = src.size();
      if ( len + idx > src.size() ) len = src.size() - idx;

      auto i = idx + WLEN;
      for ( ; i < idx + len /* && i < src.size() */; i += WLEN ) {
        dst.set_int( i - WLEN, src.get_int( i - WLEN, WLEN ), WLEN );
      }
      i -= WLEN;
      for ( ; i < idx + len; ++i ) dst[ i ] = src[ i ];
    }  /* -----  end of function bv_icopy  ----- */


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


  /**
   *  @brief  Check if the given file exists and is appendable.
   *
   *  @param  file_name The name of the file to be checked.
   *  @return `true` if exists and is appendable; otherwise `false`.
   *
   *  Use `std::ofstream` to check for writability of the given file because of
   *  portability.
   */
    inline bool
  appendable( const std::string& file_name )
  {
    if ( !readable( file_name ) ) return false;
    std::ofstream ofs( file_name, std::ofstream::app );
    return ofs.good();
  }  /* -----  end of function appendable  ----- */


    inline std::string
  get_env( const std::string& var )
  {
    const char* val = ::getenv( var.c_str() );
    if ( val == 0 ) {
      return "";
    }
    else {
      return val;
    }
  }


    inline std::string
  get_tmpdir_env( )
  {
    return get_env( "TMPDIR" );
  }


    inline std::string
  get_tmpdir( )
  {
    std::string tmpdir = get_tmpdir_env();
    if ( tmpdir.size() == 0 ) tmpdir = GREM_DEFAULT_TMPDIR;
    return tmpdir;
  }

    inline std::string
  get_tmpfile( )
  {
    std::string tmpfile_templ = get_tmpdir() + GREM_TMPFILE_TEMPLATE;
    char* tmpl = new char [ tmpfile_templ.size() + 1 ];
    std::strcpy( tmpl, tmpfile_templ.c_str() );
    int fd = mkstemp( tmpl );
    tmpfile_templ = tmpl;

    ::close( fd );
    delete[] tmpl;
    return tmpfile_templ;
  }


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
   *  @brief  Serialize an `deque` to an output stream.
   *
   *  @param[out]  out Output stream.
   *  @param[in]  v The `deque`.
   *  @param[in]  begin The begin input iterator [unused].
   *  @param[in]  end The end input iterator [unused].
   *
   *  A wraper to `serialize` member function of `deque`.
   */
  template< typename T >
    inline void
  serialize( std::ostream& out, const std::deque< T >& v )
  {
    serialize( out, v.begin(), v.end() );
  }  /* -----  end of template function serialize  ----- */


  /**
   *  @brief  Serialize an `vector` to an output stream.
   *
   *  @param[out]  out Output stream.
   *  @param[in]  v The `vector`.
   *  @param[in]  begin The begin input iterator [unused].
   *  @param[in]  end The end input iterator [unused].
   *
   *  A wraper to `serialize` member function of `vector`.
   */
  template< typename T >
    inline void
  serialize( std::ostream& out, const std::vector< T >& v )
  {
    serialize( out, v.begin(), v.end() );
  }  /* -----  end of template function serialize  ----- */


  /**
   *  @brief  Serialize an `enc_vector` to an output stream.
   *
   *  @param[out]  out Output stream.
   *  @param[in]  ev The `enc_vector`.
   *  @param[in]  begin The begin input iterator [unused].
   *  @param[in]  end The end input iterator [unused].
   *
   *  A wraper to `serialize` member function of `enc_vector`.
   */
  template< typename TCoder >
    inline void
  serialize( std::ostream& out, const sdsl::enc_vector< TCoder >& ev )
  {
    ev.serialize( out );
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


  template< typename TObject >
      inline void
    open( TObject& obj, std::istream& in )
    {
      obj.load( in );
    }


  template< typename TObject >
      inline void
    save( TObject& obj, std::ostream& out )
    {
      obj.serialize( out );
    }


  /**
   *  @brief  Assign a container to another one.
   *
   *  @param  a The container to be assigned.
   *  @param  b The container to assign.
   *
   *  Provide an interface function for assignment.
   */
  template< typename TContainer1, typename TContainer2 >
      inline void
    _assign( TContainer1& a, const TContainer2& b )
    {
      a.resize( b.size() );
      std::copy( b.begin(), b.end(), a.begin() );
    }


  /**
   *  @brief  Assign a container to another one.
   *
   *  @param  a The container to be assigned.
   *  @param  b The container to assign.
   *
   *  Provide an interface function for assignment.
   */
  template< typename T, typename TContainer >
      inline void
    assign( std::vector< T >& a, const TContainer& b )
    {
      _assign( a, b );
    }


  /**
   *  @brief  Assign a container to another one.
   *
   *  @param  a The container to be assigned.
   *  @param  b The container to assign.
   *
   *  Provide an interface function for assignment.
   */
  template< typename T, typename TContainer >
      inline void
    assign( std::deque< T >& a, const TContainer& b )
    {
      _assign( a, b );
    }


  /**
   *  @brief  Assign a container to another one.
   *
   *  @param  a The container to be assigned.
   *  @param  b The container to assign.
   *
   *  Provide an interface function for assignment.
   *
   *  @overload for `sdsl::enc_vector`.
   */
  template< typename TCoder, typename TContainer >
      inline void
    assign( sdsl::enc_vector< TCoder >& a, const TContainer& b )
    {
      sdsl::util::assign( a, b );
    }


  /**
   *  @brief  Clear the given container.
   *
   *  @param  c The container to be cleaned.
   *
   *  A wrapper function to provide an interface to `clean` member function of
   *  `std::vector`.
   */
  template< typename T >
      inline void
    clear( std::vector< T >& c )
    {
      c.clear();
    }


  /**
   *  @brief  Clear the given container.
   *
   *  @param  c The container to be cleaned.
   *
   *  A wrapper function to provide an interface to `clean` member function of
   *  `std::vector`.
   */
  template< typename T >
      inline void
    clear( std::deque< T >& c )
    {
      c.clear();
    }


  /**
   *  @brief  Clear the given container.
   *
   *  @param  c The container to be cleaned.
   *
   *  A wrapper function to provide an interface to `clean` member function of
   *  `sdsl::enc_vector`.
   */
  template< typename TCoder >
      inline void
    clear( sdsl::enc_vector< TCoder >& ev )
    {
      sdsl::util::clear( ev );
    }

  /**
   *  @brief  Reserve the required memory for the vector of size `size`.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   *  It calls reserve member function of the container.
   */
  template< typename TObject, typename TSize >
      inline void
    reserve( std::vector< TObject >& container, const TSize size )
    {
      container.reserve( size );
    }

  /**
   *  @brief  Reserve the required memory for the string of size `size`.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   *  It calls reserve member function of the container.
   */
  template< typename TSize >
      inline void
    reserve( std::string& container, const TSize size )
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
    reserve( std::deque< TObject >& container, const TSize size )
    {
      /* NOOP */
    }


  /**
   *  @overload The `std::set` cannot be reserved.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   *  Do nothing.
   */
  template< typename TObject, typename TSize >
      inline void
    reserve( std::set< TObject >& container, const TSize size )
    {
      /* NOOP */
    }


  /**
   *  @overload The `sdsl::enc_vector` cannot be reserved.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   *  Do nothing.
   */
  template< typename TCoder, typename TSize >
      inline void
    reserve( sdsl::enc_vector< TCoder >& container, const TSize size )
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
    reserve( container, size );
    for ( unsigned int i = 0; i < size; ++i ) {
      typename TContainer::value_type item;
      deserialize( in, item );
      *itr++ = item;
    }
  }  /* -----  end of template function deserialize  ----- */

  template< typename T >
    inline void
  deserialize( std::istream& in, std::vector< T >& container )
  {
    deserialize( in, container, std::back_inserter( container ) );
  }  /* -----  end of template function deserialize  ----- */

  template< typename T >
    inline void
  deserialize( std::istream& in, std::deque< T >& container )
  {
    deserialize( in, container, std::back_inserter( container ) );
  }  /* -----  end of template function deserialize  ----- */

  template< typename TCoder >
    inline void
  deserialize( std::istream& in, sdsl::enc_vector< TCoder >& ev )
  {
    ev.load( in );
  }  /* -----  end of template function deserialize  ----- */

  /* Meta-functions */
  template< typename T1, typename T2 >
    using enable_if_equal = std::enable_if< std::is_same< T1, T2 >::value, T2 >;

  template< typename T1, typename T2 >
    using enable_if_not_equal = std::enable_if< !std::is_same< T1, T2 >::value, T2 >;

  template< typename T1, typename T2 >
    using enable_if_equal_t = typename enable_if_equal< T1, T2 >::type;

  template< typename T1, typename T2 >
    using enable_if_not_equal_t = typename enable_if_not_equal< T1, T2 >::type;
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef UTILS_H__  ----- */

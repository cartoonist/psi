/**
 *    @file  tests_fmindex.h
 *   @brief  Test grem::FM-Index module.
 *
 *  Test scenarios for FM-Index class in grem (not SeqAn).
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Fri Feb 16, 2018  20:08
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <string>

#include "tests_base.h"
#include "fmindex.h"
#include "logger.h"


using namespace grem;

SCENARIO( "Find occurrences of a pattern in a text using FM-index", "[fmindex]" )
{
  GIVEN( "An index based on a disk-based string" )
  {
    typedef YaString< DiskBased > string_type;
    typedef seqan::Index< string_type, grem::FMIndex<> > index_type;

    string_type text( "a-mississippian-lazy-fox-sits-on-a-pie" );
    index_type index( text );

    GIVEN( "A Finder based on that index" )
    {
      seqan::Finder< index_type > finder( index );

      WHEN( "An existing pattern is searched using that Finder" )
      {
        std::set< index_type::savalue_type > true_occs = { 5, 8, 25 };
        std::vector< index_type::savalue_type > occs;
        while ( find( finder, "si" ) ) occs.push_back( beginPosition( finder ) );

        THEN( "It should be find all occurrences" )
        {
          for ( auto loc : occs ) {
            REQUIRE( true_occs.find( loc ) != true_occs.end() );
          }
        }
      }

      WHEN( "Another existing pattern is searched using that Finder" )
      {
        std::string pattern( "pi" );
        std::set< index_type::savalue_type > true_occs = { 11, 35 };
        std::vector< index_type::savalue_type > occs;
        while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

        THEN( "It should be find all occurrences" )
        {
          for ( auto loc : occs ) {
            REQUIRE( true_occs.find( loc ) != true_occs.end() );
          }
        }
      }

      WHEN( "Another existing pattern is searched using that Finder" )
      {
        std::string pattern( "misissipian" );
        std::set< index_type::savalue_type > true_occs;
        std::vector< index_type::savalue_type > occs;
        while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

        THEN( "It should be find all occurrences" )
        {
          REQUIRE( occs.size() == 0 );
        }
      }
    }
  }

  GIVEN( "An index based on a in-memory string" )
  {
    typedef YaString< InMemory > string_type;
    typedef seqan::Index< string_type, grem::FMIndex<> > index_type;

    string_type text( "a-mississippian-lazy-fox-sits-on-a-pie" );
    index_type index( text );

    GIVEN( "A Finder based on that index" )
    {
      seqan::Finder< index_type > finder( index );

      WHEN( "An existing pattern is searched using that Finder" )
      {
        std::set< index_type::savalue_type > true_occs = { 5, 8, 25 };
        std::vector< index_type::savalue_type > occs;
        while ( find( finder, "si" ) ) occs.push_back( beginPosition( finder ) );

        THEN( "It should be find all occurrences" )
        {
          for ( auto loc : occs ) {
            REQUIRE( true_occs.find( loc ) != true_occs.end() );
          }
        }
      }

      WHEN( "Another existing pattern is searched using that Finder" )
      {
        std::string pattern( "pi" );
        std::set< index_type::savalue_type > true_occs = { 11, 35 };
        std::vector< index_type::savalue_type > occs;
        while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

        THEN( "It should be find all occurrences" )
        {
          for ( auto loc : occs ) {
            REQUIRE( true_occs.find( loc ) != true_occs.end() );
          }
        }
      }

      WHEN( "Another existing pattern is searched using that Finder" )
      {
        std::string pattern( "misissipian" );
        std::set< index_type::savalue_type > true_occs;
        std::vector< index_type::savalue_type > occs;
        while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

        THEN( "It should be find all occurrences" )
        {
          REQUIRE( occs.size() == 0 );
        }
      }
    }
  }
}

SCENARIO( "Find occurrences of a pattern in a string set using FM-index", "[fmindex]" )
{
  GIVEN( "An index based on a disk-based string set" )
  {
    typedef seqan::StringSet< DiskString > stringset_type;
    typedef seqan::Index< stringset_type, grem::FMIndex<> > index_type;

    std::string str1 = "a-mississippian-lazy-fox-sits-on-a-pie";
    std::string str2 = "another-brazilian-cute-beaver-builds-a-dam";
    std::string str3 = "some-african-stupid-chimps-eat-banana";
    stringset_type text;
    text.push_back( str1 );
    text.push_back( str2 );
    text.push_back( str3 );
    index_type index( text );

    GIVEN( "A Finder based on that index" )
    {
      seqan::Finder< index_type > finder( index );

      WHEN( "An existing pattern is searched using that Finder" )
      {
        std::set< index_type::pos_type > true_occs = { { 2, 32 }, { 2, 34 } };
        std::vector< index_type::pos_type > occs;
        while ( find( finder, "ana" ) ) occs.push_back( beginPosition( finder ) );

        THEN( "It should be find all occurrences" )
        {
          for ( auto loc : occs ) {
            REQUIRE( true_occs.find( loc ) != true_occs.end() );
          }
        }
      }

      WHEN( "Another existing pattern is searched using that Finder" )
      {
        std::string pattern( "pi" );
        std::set< index_type::pos_type > true_occs = { { 0, 11 }, { 0, 35 }, { 2, 16 } };
        std::vector< index_type::pos_type > occs;
        while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

        THEN( "It should be find all occurrences" )
        {
          for ( auto loc : occs ) {
            REQUIRE( true_occs.find( loc ) != true_occs.end() );
          }
        }
      }

      WHEN( "Another existing pattern is searched using that Finder" )
      {
        std::string pattern( "pieano" );
        std::set< index_type::pos_type > true_occs;
        std::vector< index_type::pos_type > occs;
        while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

        THEN( "It should be find all occurrences" )
        {
          REQUIRE( occs.size() == 0 );
        }
      }
    }
  }

  GIVEN( "An index based on a in-memory string set" )
  {
    typedef seqan::StringSet< MemString > stringset_type;
    typedef seqan::Index< stringset_type, grem::FMIndex<> > index_type;

    std::string str1 = "a-mississippian-lazy-fox-sits-on-a-pie";
    std::string str2 = "another-brazilian-cute-beaver-builds-a-dam";
    std::string str3 = "some-african-stupid-chimps-eat-banana";
    stringset_type text;
    text.push_back( str1 );
    text.push_back( str2 );
    text.push_back( str3 );
    index_type index( text );

    GIVEN( "A Finder based on that index" )
    {
      seqan::Finder< index_type > finder( index );

      WHEN( "An existing pattern is searched using that Finder" )
      {
        std::set< index_type::pos_type > true_occs = { { 2, 32 }, { 2, 34 } };
        std::vector< index_type::pos_type > occs;
        while ( find( finder, "ana" ) ) occs.push_back( beginPosition( finder ) );

        THEN( "It should be find all occurrences" )
        {
          for ( auto loc : occs ) {
            REQUIRE( true_occs.find( loc ) != true_occs.end() );
          }
        }
      }

      WHEN( "Another existing pattern is searched using that Finder" )
      {
        std::string pattern( "pi" );
        std::set< index_type::pos_type > true_occs = { { 0, 11 }, { 0, 35 }, { 2, 16 } };
        std::vector< index_type::pos_type > occs;
        while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

        THEN( "It should be find all occurrences" )
        {
          for ( auto loc : occs ) {
            REQUIRE( true_occs.find( loc ) != true_occs.end() );
          }
        }
      }

      WHEN( "Another existing pattern is searched using that Finder" )
      {
        std::string pattern( "pieano" );
        std::set< index_type::pos_type > true_occs;
        std::vector< index_type::pos_type > occs;
        while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

        THEN( "It should be find all occurrences" )
        {
          REQUIRE( occs.size() == 0 );
        }
      }
    }
  }
}

SCENARIO( "Save and load FM-index on string", "[fmindex]" )
{
  GIVEN( "An index based on a disk-based string serialized to the disk" )
  {
    typedef YaString< DiskBased > string_type;
    typedef seqan::Index< string_type, grem::FMIndex<> > index_type;

    string_type text( "a-mississippian-lazy-fox-sits-on-a-pie" );
    index_type index1( text );
    indexRequire( index1, seqan::FibreSALF() );
    std::string fpath = SEQAN_TEMP_FILENAME();
    save( index1, fpath );

    GIVEN( "Another index constructed by default constructor" )
    {
      index_type index2;

      WHEN( "It is opened" )
      {
        open( index2, fpath );

        THEN( "Requesting for construction should be ignored without throwing exception" )
        {
          REQUIRE_NOTHROW( indexRequire( index2, seqan::FibreSALF() ) );
        }
      }
    }

    GIVEN( "An index loaded from the saved file" )
    {
      index_type index2;
      open( index2, fpath );

      GIVEN( "A Finder based on that index" )
      {
        seqan::Finder< index_type > finder( index2 );

        WHEN( "An existing pattern is searched using that Finder" )
        {
          std::set< index_type::savalue_type > true_occs = { 5, 8, 25 };
          std::vector< index_type::savalue_type > occs;
          while ( find( finder, "si" ) ) occs.push_back( beginPosition( finder ) );

          THEN( "It should be find all occurrences" )
          {
            for ( auto loc : occs ) {
              REQUIRE( true_occs.find( loc ) != true_occs.end() );
            }
          }
        }

        WHEN( "Another existing pattern is searched using that Finder" )
        {
          std::string pattern( "pi" );
          std::set< index_type::savalue_type > true_occs = { 11, 35 };
          std::vector< index_type::savalue_type > occs;
          while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

          THEN( "It should be find all occurrences" )
          {
            for ( auto loc : occs ) {
              REQUIRE( true_occs.find( loc ) != true_occs.end() );
            }
          }
        }

        WHEN( "Another existing pattern is searched using that Finder" )
        {
          std::string pattern( "misissipian" );
          std::set< index_type::savalue_type > true_occs;
          std::vector< index_type::savalue_type > occs;
          while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

          THEN( "It should be find all occurrences" )
          {
            REQUIRE( occs.size() == 0 );
          }
        }
      }
    }
  }

  GIVEN( "An index based on a in-memory string serialized to the disk" )
  {
    typedef YaString< InMemory > string_type;
    typedef seqan::Index< string_type, grem::FMIndex<> > index_type;

    string_type text( "a-mississippian-lazy-fox-sits-on-a-pie" );
    index_type index1( text );
    indexRequire( index1, seqan::FibreSALF() );
    std::string fpath = SEQAN_TEMP_FILENAME();
    save( index1, fpath );

    GIVEN( "Another index constructed by default constructor" )
    {
      index_type index2;

      WHEN( "It is opened" )
      {
        open( index2, fpath );

        THEN( "Requesting for construction should be ignored without throwing exception" )
        {
          REQUIRE_NOTHROW( indexRequire( index2, seqan::FibreSALF() ) );
        }
      }
    }

    GIVEN( "An index loaded from the saved file" )
    {
      index_type index2;
      open( index2, fpath );

      GIVEN( "A Finder based on that index" )
      {
        seqan::Finder< index_type > finder( index2 );

        WHEN( "An existing pattern is searched using that Finder" )
        {
          std::set< index_type::savalue_type > true_occs = { 5, 8, 25 };
          std::vector< index_type::savalue_type > occs;
          while ( find( finder, "si" ) ) occs.push_back( beginPosition( finder ) );

          THEN( "It should be find all occurrences" )
          {
            for ( auto loc : occs ) {
              REQUIRE( true_occs.find( loc ) != true_occs.end() );
            }
          }
        }

        WHEN( "Another existing pattern is searched using that Finder" )
        {
          std::string pattern( "pi" );
          std::set< index_type::savalue_type > true_occs = { 11, 35 };
          std::vector< index_type::savalue_type > occs;
          while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

          THEN( "It should be find all occurrences" )
          {
            for ( auto loc : occs ) {
              REQUIRE( true_occs.find( loc ) != true_occs.end() );
            }
          }
        }

        WHEN( "Another existing pattern is searched using that Finder" )
        {
          std::string pattern( "misissipian" );
          std::set< index_type::savalue_type > true_occs;
          std::vector< index_type::savalue_type > occs;
          while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

          THEN( "It should be find all occurrences" )
          {
            REQUIRE( occs.size() == 0 );
          }
        }
      }
    }
  }
}

SCENARIO( "Save and load FM-index on stringset", "[fmindex]" )
{
  GIVEN( "An index based on a disk-based string set serialized to the disk" )
  {
    typedef seqan::StringSet< DiskString > stringset_type;
    typedef seqan::Index< stringset_type, grem::FMIndex<> > index_type;

    std::string str1 = "a-mississippian-lazy-fox-sits-on-a-pie";
    std::string str2 = "another-brazilian-cute-beaver-builds-a-dam";
    std::string str3 = "some-african-stupid-chimps-eat-banana";
    stringset_type text;
    text.push_back( str1 );
    text.push_back( str2 );
    text.push_back( str3 );
    index_type index1( text );
    indexRequire( index1, seqan::FibreSALF() );
    std::string fpath = SEQAN_TEMP_FILENAME();
    save( index1, fpath );

    GIVEN( "Another index constructed by default constructor" )
    {
      index_type index2;

      WHEN( "It is opened" )
      {
        open( index2, fpath );

        THEN( "Requesting for construction should be ignored without throwing exception" )
        {
          REQUIRE_NOTHROW( indexRequire( index2, seqan::FibreSALF() ) );
        }
      }
    }

    GIVEN( "An index loaded from the saved file" )
    {
      index_type index2;
      open( index2, fpath );

      GIVEN( "A Finder based on that index" )
      {
        seqan::Finder< index_type > finder( index2 );

        WHEN( "An existing pattern is searched using that Finder" )
        {
          std::set< index_type::pos_type > true_occs = { { 2, 32 }, { 2, 34 } };
          std::vector< index_type::pos_type > occs;
          while ( find( finder, "ana" ) ) occs.push_back( beginPosition( finder ) );

          THEN( "It should be find all occurrences" )
          {
            for ( auto loc : occs ) {
              REQUIRE( true_occs.find( loc ) != true_occs.end() );
            }
          }
        }

        WHEN( "Another existing pattern is searched using that Finder" )
        {
          std::string pattern( "pi" );
          std::set< index_type::pos_type > true_occs = { { 0, 11 }, { 0, 35 }, { 2, 16 } };
          std::vector< index_type::pos_type > occs;
          while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

          THEN( "It should be find all occurrences" )
          {
            for ( auto loc : occs ) {
              REQUIRE( true_occs.find( loc ) != true_occs.end() );
            }
          }
        }

        WHEN( "Another existing pattern is searched using that Finder" )
        {
          std::string pattern( "pieano" );
          std::set< index_type::pos_type > true_occs;
          std::vector< index_type::pos_type > occs;
          while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

          THEN( "It should be find all occurrences" )
          {
            REQUIRE( occs.size() == 0 );
          }
        }
      }
    }
  }

  GIVEN( "An index based on a in-memory string set serialized to the disk" )
  {
    typedef seqan::StringSet< MemString > stringset_type;
    typedef seqan::Index< stringset_type, grem::FMIndex<> > index_type;

    std::string str1 = "a-mississippian-lazy-fox-sits-on-a-pie";
    std::string str2 = "another-brazilian-cute-beaver-builds-a-dam";
    std::string str3 = "some-african-stupid-chimps-eat-banana";
    stringset_type text;
    text.push_back( str1 );
    text.push_back( str2 );
    text.push_back( str3 );
    index_type index1( text );
    indexRequire( index1, seqan::FibreSALF() );
    std::string fpath = SEQAN_TEMP_FILENAME();
    save( index1, fpath );

    GIVEN( "Another index constructed by default constructor" )
    {
      index_type index2;

      WHEN( "It is opened" )
      {
        open( index2, fpath );

        THEN( "Requesting for construction should be ignored without throwing exception" )
        {
          REQUIRE_NOTHROW( indexRequire( index2, seqan::FibreSALF() ) );
        }
      }
    }

    GIVEN( "An index loaded from the saved file" )
    {
      index_type index2;
      open( index2, fpath );

      GIVEN( "A Finder based on that index" )
      {
        seqan::Finder< index_type > finder( index2 );

        WHEN( "An existing pattern is searched using that Finder" )
        {
          std::set< index_type::pos_type > true_occs = { { 2, 32 }, { 2, 34 } };
          std::vector< index_type::pos_type > occs;
          while ( find( finder, "ana" ) ) occs.push_back( beginPosition( finder ) );

          THEN( "It should be find all occurrences" )
          {
            for ( auto loc : occs ) {
              REQUIRE( true_occs.find( loc ) != true_occs.end() );
            }
          }
        }

        WHEN( "Another existing pattern is searched using that Finder" )
        {
          std::string pattern( "pi" );
          std::set< index_type::pos_type > true_occs = { { 0, 11 }, { 0, 35 }, { 2, 16 } };
          std::vector< index_type::pos_type > occs;
          while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

          THEN( "It should be find all occurrences" )
          {
            for ( auto loc : occs ) {
              REQUIRE( true_occs.find( loc ) != true_occs.end() );
            }
          }
        }

        WHEN( "Another existing pattern is searched using that Finder" )
        {
          std::string pattern( "pieano" );
          std::set< index_type::pos_type > true_occs;
          std::vector< index_type::pos_type > occs;
          while ( find( finder, pattern ) ) occs.push_back( beginPosition( finder ) );

          THEN( "It should be find all occurrences" )
          {
            REQUIRE( occs.size() == 0 );
          }
        }
      }
    }
  }
}

SCENARIO( "Traverse prefix tree of a string using FM-index index iterator", "[fmindex][iterator]" )
{
  GIVEN( "An index based on a disk-based string" )
  {
    typedef YaString< DiskBased > string_type;
    typedef seqan::Index< string_type, grem::FMIndex<> > index_type;

    string_type text( "a-mississippian-lazy-fox-sits-on-a-pie" );
    index_type index( text );

    GIVEN( "A top-down iterator for virtual prefix tree of the string" )
    {
      typedef typename seqan::Iterator< index_type, seqan::TopDown<> >::Type iterator_type;

      iterator_type it( index );

      REQUIRE( isRoot( it ) );
      REQUIRE( !goUp( it ) );
      REQUIRE( !goRight( it ) );
      REQUIRE( parentEdgeLength( it ) == 0 );
      REQUIRE( parentEdgeLabel( it ) == "" );

      WHEN( "Get occurrences for root node" )
      {
        auto occs = getOccurrences( it );

        THEN( "It should report nothing" ){
          REQUIRE( occs.size() == 0 );
        }
      }

      WHEN( "Calling goDown by a character" )
      {
        REQUIRE( goDown( it, 'i' ) );

        THEN( "It should go down prefix tree along the given character" )
        {
          REQUIRE( !isRoot( it ) );
          REQUIRE( repLength( it ) == 1 );
          REQUIRE( representative( it ) == "i" );
          REQUIRE( parentEdgeLength( it ) == 1 );
          REQUIRE( parentEdgeLabel( it ) == "i" );
        }

        WHEN( "Getting occurrences of the iterator" )
        {
          std::set< index_type::savalue_type > truth = { 3, 6, 9, 12, 26, 36 };
          auto occs = getOccurrences( it );

          THEN( "It should report all occurrences indicated by the iterator" )
          {
            REQUIRE( occs.size() == 6 );
            for ( auto it = occs.begin(); it != occs.end(); ++it ) {
              REQUIRE( truth.find( *it ) != truth.end() );
            }
          }
        }

        WHEN( "Going down more" )
        {
          REQUIRE( goDown( it, 's' ) );
          REQUIRE( goDown( it, 's' ) );

          THEN( "It should point to corresponding node in the virtual prefix tree" )
          {
            REQUIRE( !isRoot( it ) );
            REQUIRE( repLength( it ) == 3 );
            REQUIRE( representative( it ) == "ssi" );
            REQUIRE( parentEdgeLength( it ) == 1 );
            REQUIRE( parentEdgeLabel( it ) == "s" );
          }

          WHEN( "Getting occurrences of the iterator" )
          {
            std::set< index_type::savalue_type > truth = { 4, 7 };
            auto occs = getOccurrences( it );

            THEN( "It should report all occurrences indicated by the iterator" )
            {
              REQUIRE( occs.size() == 2 );
              for ( auto it = occs.begin(); it != occs.end(); ++it ) {
                REQUIRE( truth.find( *it ) != truth.end() );
              }
            }
          }

          WHEN( "Going right" )
          {
            THEN( "It cannot go right" )
            {
              REQUIRE( !goRight( it ) );
              REQUIRE( representative( it ) == "ssi" );
            }
          }

          WHEN( "Going up then go down" )
          {
            REQUIRE( goUp( it ) );
            REQUIRE( goDown( it, '-' ) );

            THEN( "It should point to corresponding node in the virtual prefix tree" )
            {
              REQUIRE( repLength( it ) == 3 );
              REQUIRE( representative( it ) == "-si" );
              REQUIRE( parentEdgeLength( it ) == 1 );
              REQUIRE( parentEdgeLabel( it ) == "-" );
            }

            WHEN( "Getting occurrences of the iterator" )
            {
              auto occs = getOccurrences( it );

              THEN( "It should report all occurrences indicated by the iterator" )
              {
                REQUIRE( occs.size() == 1 );
                REQUIRE( *occs.begin() == 24 );
              }
            }

            WHEN( "Going right" )
            {
              REQUIRE( goRight( it ) );

              THEN( "It should point to corresponding node in the virtual prefix tree" )
              {
                REQUIRE( repLength( it ) == 3 );
                REQUIRE( representative( it ) == "ssi" );
                REQUIRE( parentEdgeLength( it ) == 1 );
                REQUIRE( parentEdgeLabel( it ) == "s" );
              }

              WHEN( "Getting occurrences of the iterator" )
              {
                std::set< index_type::savalue_type > truth = { 4, 7 };
                auto occs = getOccurrences( it );

                THEN( "It should report all occurrences indicated by the iterator" )
                {
                  REQUIRE( occs.size() == 2 );
                  for ( auto it = occs.begin(); it != occs.end(); ++it ) {
                    REQUIRE( truth.find( *it ) != truth.end() );
                  }
                }
              }
            }
          }

          WHEN( "Go root by calling goRoot interface function" )
          {
            goRoot( it );

            THEN( "It should point root node" )
            {
              REQUIRE( isRoot( it ) );
              REQUIRE( !goUp( it ) );
              REQUIRE( !goRight( it ) );
              REQUIRE( parentEdgeLength( it ) == 0 );
              REQUIRE( parentEdgeLabel( it ) == "" );
            }
          }

          WHEN( "Go root by going up" )
          {
            goUp( it );
            goUp( it );
            goUp( it );

            THEN( "It should point root node" )
            {
              REQUIRE( isRoot( it ) );
              REQUIRE( !goUp( it ) );
              REQUIRE( !goRight( it ) );
              REQUIRE( parentEdgeLength( it ) == 0 );
              REQUIRE( parentEdgeLabel( it ) == "" );
            }
          }
        }
      }
    }
  }

  GIVEN( "An index based on a in-memory string" )
  {
    typedef YaString< InMemory > string_type;
    typedef seqan::Index< string_type, grem::FMIndex<> > index_type;

    string_type text( "a-mississippian-lazy-fox-sits-on-a-pie" );
    index_type index( text );

    GIVEN( "A top-down iterator for virtual prefix tree of the string" )
    {
      typedef typename seqan::Iterator< index_type, seqan::TopDown< seqan::ParentLinks<> > >::Type iterator_type;

      iterator_type it( index );

      REQUIRE( isRoot( it ) );
      REQUIRE( !goUp( it ) );
      REQUIRE( !goRight( it ) );
      REQUIRE( parentEdgeLength( it ) == 0 );
      REQUIRE( parentEdgeLabel( it ) == "" );

      WHEN( "Get occurrences for root node" )
      {
        auto occs = getOccurrences( it );

        THEN( "It should report nothing" ){
          REQUIRE( occs.size() == 0 );
        }
      }

      WHEN( "Calling goDown by a character" )
      {
        REQUIRE( goDown( it, 'i' ) );

        THEN( "It should go down prefix tree along the given character" )
        {
          REQUIRE( !isRoot( it ) );
          REQUIRE( repLength( it ) == 1 );
          REQUIRE( representative( it ) == "i" );
          REQUIRE( parentEdgeLength( it ) == 1 );
          REQUIRE( parentEdgeLabel( it ) == "i" );
        }

        WHEN( "Getting occurrences of the iterator" )
        {
          std::set< index_type::savalue_type > truth = { 3, 6, 9, 12, 26, 36 };
          auto occs = getOccurrences( it );

          THEN( "It should report all occurrences indicated by the iterator" )
          {
            REQUIRE( occs.size() == 6 );
            for ( auto it = occs.begin(); it != occs.end(); ++it ) {
              REQUIRE( truth.find( *it ) != truth.end() );
            }
          }
        }

        WHEN( "Going down more" )
        {
          REQUIRE( goDown( it, 's' ) );
          REQUIRE( goDown( it, 's' ) );

          THEN( "It should point to corresponding node in the virtual prefix tree" )
          {
            REQUIRE( !isRoot( it ) );
            REQUIRE( repLength( it ) == 3 );
            REQUIRE( representative( it ) == "ssi" );
            REQUIRE( parentEdgeLength( it ) == 1 );
            REQUIRE( parentEdgeLabel( it ) == "s" );
          }

          WHEN( "Getting occurrences of the iterator" )
          {
            std::set< index_type::savalue_type > truth = { 4, 7 };
            auto occs = getOccurrences( it );

            THEN( "It should report all occurrences indicated by the iterator" )
            {
              REQUIRE( occs.size() == 2 );
              for ( auto it = occs.begin(); it != occs.end(); ++it ) {
                REQUIRE( truth.find( *it ) != truth.end() );
              }
            }
          }

          WHEN( "Going right" )
          {
            THEN( "It cannot go right" )
            {
              REQUIRE( !goRight( it ) );
              REQUIRE( representative( it ) == "ssi" );
            }
          }

          WHEN( "Going up then go down" )
          {
            REQUIRE( goUp( it ) );
            REQUIRE( goDown( it, '-' ) );

            THEN( "It should point to corresponding node in the virtual prefix tree" )
            {
              REQUIRE( repLength( it ) == 3 );
              REQUIRE( representative( it ) == "-si" );
              REQUIRE( parentEdgeLength( it ) == 1 );
              REQUIRE( parentEdgeLabel( it ) == "-" );
            }

            WHEN( "Getting occurrences of the iterator" )
            {
              auto occs = getOccurrences( it );

              THEN( "It should report all occurrences indicated by the iterator" )
              {
                REQUIRE( occs.size() == 1 );
                REQUIRE( *occs.begin() == 24 );
              }
            }

            WHEN( "Going right" )
            {
              REQUIRE( goRight( it ) );

              THEN( "It should point to corresponding node in the virtual prefix tree" )
              {
                REQUIRE( repLength( it ) == 3 );
                REQUIRE( representative( it ) == "ssi" );
                REQUIRE( parentEdgeLength( it ) == 1 );
                REQUIRE( parentEdgeLabel( it ) == "s" );
              }

              WHEN( "Getting occurrences of the iterator" )
              {
                std::set< index_type::savalue_type > truth = { 4, 7 };
                auto occs = getOccurrences( it );

                THEN( "It should report all occurrences indicated by the iterator" )
                {
                  REQUIRE( occs.size() == 2 );
                  for ( auto it = occs.begin(); it != occs.end(); ++it ) {
                    REQUIRE( truth.find( *it ) != truth.end() );
                  }
                }
              }
            }
          }

          WHEN( "Go root by calling goRoot interface function" )
          {
            goRoot( it );

            THEN( "It should point root node" )
            {
              REQUIRE( isRoot( it ) );
              REQUIRE( !goUp( it ) );
              REQUIRE( !goRight( it ) );
              REQUIRE( parentEdgeLength( it ) == 0 );
              REQUIRE( parentEdgeLabel( it ) == "" );
            }
          }

          WHEN( "Go root by going up" )
          {
            goUp( it );
            goUp( it );
            goUp( it );

            THEN( "It should point root node" )
            {
              REQUIRE( isRoot( it ) );
              REQUIRE( !goUp( it ) );
              REQUIRE( !goRight( it ) );
              REQUIRE( parentEdgeLength( it ) == 0 );
              REQUIRE( parentEdgeLabel( it ) == "" );
            }
          }
        }
      }
    }
  }
}

SCENARIO( "Traverse prefix tree of a string set using FM-index index iterator", "[fmindex][iterator]" )
{
  GIVEN( "An index based on a disk-based string set serialized to the disk" )
  {
    typedef seqan::StringSet< DiskString > stringset_type;
    typedef seqan::Index< stringset_type, grem::FMIndex<> > index_type;

    std::string str1 = "a-mississippian-lazy-fox-sits-on-a-pie";
    std::string str2 = "another-brazilian-cute-beaver-builds-a-dam";
    std::string str3 = "some-african-stupid-chimps-eat-banana";
    stringset_type text;
    text.push_back( str1 );
    text.push_back( str2 );
    text.push_back( str3 );
    index_type index1( text );
    indexRequire( index1, seqan::FibreSALF() );
    std::string fpath = SEQAN_TEMP_FILENAME();
    save( index1, fpath );

    GIVEN( "An index loaded from the saved file" )
    {
      index_type index2;
      open( index2, fpath );

      GIVEN( "A top-down iterator for virtual prefix tree of the string" )
      {
        typedef typename seqan::Iterator< index_type, seqan::TopDown<> >::Type iterator_type;

        iterator_type it( index2 );

        REQUIRE( isRoot( it ) );
        REQUIRE( !goUp( it ) );
        REQUIRE( !goRight( it ) );
        REQUIRE( parentEdgeLength( it ) == 0 );
        REQUIRE( parentEdgeLabel( it ) == "" );

        WHEN( "Get occurrences for root node" )
        {
          auto occs = getOccurrences( it );

          THEN( "It should report nothing" ){
            REQUIRE( occs.size() == 0 );
          }
        }

        WHEN( "Calling goDown by a character" )
        {
          REQUIRE( goDown( it, 'n' ) );

          THEN( "It should go down prefix tree along the given character" )
          {
            REQUIRE( !isRoot( it ) );
            REQUIRE( repLength( it ) == 1 );
            REQUIRE( representative( it ) == "n" );
            REQUIRE( parentEdgeLength( it ) == 1 );
            REQUIRE( parentEdgeLabel( it ) == "n" );
          }

          WHEN( "Getting occurrences of the iterator" )
          {
            std::set< index_type::pos_type > truth = { { 0, 14 }, { 0, 31 }, { 1, 1 }, { 1, 16 }, { 2, 11 }, { 2, 33 }, { 2, 35 } };
            auto occs = getOccurrences( it );

            THEN( "It should report all occurrences indicated by the iterator" )
            {
              REQUIRE( occs.size() == 7 );
              for ( auto it = occs.begin(); it != occs.end(); ++it ) {
                  REQUIRE( truth.find( *it ) != truth.end() );
              }
            }
          }

          WHEN( "Going down more" )
          {
            REQUIRE( goDown( it, 'a' ) );
            REQUIRE( goDown( it, 'n' ) );

            THEN( "It should point to corresponding node in the virtual prefix tree" )
            {
              REQUIRE( !isRoot( it ) );
              REQUIRE( repLength( it ) == 3 );
              REQUIRE( representative( it ) == "nan" );
              REQUIRE( parentEdgeLength( it ) == 1 );
              REQUIRE( parentEdgeLabel( it ) == "n" );
            }

            WHEN( "Getting occurrences of the iterator" )
            {
              auto occs = getOccurrences( it );

              THEN( "It should report all occurrences indicated by the iterator" )
              {
                REQUIRE( occs.size() == 1 );
                REQUIRE( (*occs.begin()).i1 == 2 );
                REQUIRE( (*occs.begin()).i2 == 33 );
              }
            }

            WHEN( "Going right" )
            {
              THEN( "It cannot go right" )
              {
                REQUIRE( !goRight( it ) );
                REQUIRE( representative( it ) == "nan" );
              }
            }

            WHEN( "Going up then go down" )
            {
              REQUIRE( goUp( it ) );
              REQUIRE( goDown( it, 'c' ) );

              THEN( "It should point to corresponding node in the virtual prefix tree" )
              {
                REQUIRE( repLength( it ) == 3 );
                REQUIRE( representative( it ) == "can" );
                REQUIRE( parentEdgeLength( it ) == 1 );
                REQUIRE( parentEdgeLabel( it ) == "c" );
              }

              WHEN( "Getting occurrences of the iterator" )
              {
                auto occs = getOccurrences( it );

                THEN( "It should report all occurrences indicated by the iterator" )
                {
                  REQUIRE( occs.size() == 1 );
                  REQUIRE( (*occs.begin()).i1 == 2 );
                  REQUIRE( (*occs.begin()).i2 == 9 );
                }
              }

              WHEN( "Going right" )
              {
                REQUIRE( goRight( it ) );

                THEN( "It should point to corresponding node in the virtual prefix tree" )
                {
                  REQUIRE( repLength( it ) == 3 );
                  REQUIRE( representative( it ) == "ian" );
                  REQUIRE( parentEdgeLength( it ) == 1 );
                  REQUIRE( parentEdgeLabel( it ) == "i" );
                }

                WHEN( "Getting occurrences of the iterator" )
                {
                  std::set< index_type::pos_type > truth = { { 0, 12 }, { 1, 14 } };
                  auto occs = getOccurrences( it );

                  THEN( "It should report all occurrences indicated by the iterator" )
                  {
                    REQUIRE( occs.size() == 2 );
                    for ( auto it = occs.begin(); it != occs.end(); ++it ) {
                      REQUIRE( truth.find( *it ) != truth.end() );
                    }
                  }
                }
              }
            }

            WHEN( "Go root by calling goRoot interface function" )
            {
              goRoot( it );

              THEN( "It should point root node" )
              {
                REQUIRE( isRoot( it ) );
                REQUIRE( !goUp( it ) );
                REQUIRE( !goRight( it ) );
                REQUIRE( parentEdgeLength( it ) == 0 );
                REQUIRE( parentEdgeLabel( it ) == "" );
              }
            }

            WHEN( "Go root by going up" )
            {
              goUp( it );
              goUp( it );
              goUp( it );

              THEN( "It should point root node" )
              {
                REQUIRE( isRoot( it ) );
                REQUIRE( !goUp( it ) );
                REQUIRE( !goRight( it ) );
                REQUIRE( parentEdgeLength( it ) == 0 );
                REQUIRE( parentEdgeLabel( it ) == "" );
              }
            }
          }
        }
      }
    }
  }

  GIVEN( "An index based on a in-memory string set serialized to the disk" )
  {
    typedef seqan::StringSet< MemString > stringset_type;
    typedef seqan::Index< stringset_type, grem::FMIndex<> > index_type;

    std::string str1 = "a-mississippian-lazy-fox-sits-on-a-pie";
    std::string str2 = "another-brazilian-cute-beaver-builds-a-dam";
    std::string str3 = "some-african-stupid-chimps-eat-banana";
    stringset_type text;
    text.push_back( str1 );
    text.push_back( str2 );
    text.push_back( str3 );
    index_type index1( text );
    indexRequire( index1, seqan::FibreSALF() );
    std::string fpath = SEQAN_TEMP_FILENAME();
    save( index1, fpath );

    GIVEN( "An index loaded from the saved file" )
    {
      index_type index2;
      open( index2, fpath );

      GIVEN( "A top-down iterator for virtual prefix tree of the string" )
      {
        typedef typename seqan::Iterator< index_type, seqan::TopDown<> >::Type iterator_type;

        iterator_type it( index2 );

        REQUIRE( isRoot( it ) );
        REQUIRE( !goUp( it ) );
        REQUIRE( !goRight( it ) );
        REQUIRE( parentEdgeLength( it ) == 0 );
        REQUIRE( parentEdgeLabel( it ) == "" );

        WHEN( "Get occurrences for root node" )
        {
          auto occs = getOccurrences( it );

          THEN( "It should report nothing" ){
            REQUIRE( occs.size() == 0 );
          }
        }

        WHEN( "Calling goDown by a character" )
        {
          REQUIRE( goDown( it, 'n' ) );

          THEN( "It should go down prefix tree along the given character" )
          {
            REQUIRE( !isRoot( it ) );
            REQUIRE( repLength( it ) == 1 );
            REQUIRE( representative( it ) == "n" );
            REQUIRE( parentEdgeLength( it ) == 1 );
            REQUIRE( parentEdgeLabel( it ) == "n" );
          }

          WHEN( "Getting occurrences of the iterator" )
          {
            std::set< index_type::pos_type > truth = { { 0, 14 }, { 0, 31 }, { 1, 1 }, { 1, 16 }, { 2, 11 }, { 2, 33 }, { 2, 35 } };
            auto occs = getOccurrences( it );

            THEN( "It should report all occurrences indicated by the iterator" )
            {
              REQUIRE( occs.size() == 7 );
              for ( auto it = occs.begin(); it != occs.end(); ++it ) {
                  REQUIRE( truth.find( *it ) != truth.end() );
              }
            }
          }

          WHEN( "Going down more" )
          {
            REQUIRE( goDown( it, 'a' ) );
            REQUIRE( goDown( it, 'n' ) );

            THEN( "It should point to corresponding node in the virtual prefix tree" )
            {
              REQUIRE( !isRoot( it ) );
              REQUIRE( repLength( it ) == 3 );
              REQUIRE( representative( it ) == "nan" );
              REQUIRE( parentEdgeLength( it ) == 1 );
              REQUIRE( parentEdgeLabel( it ) == "n" );
            }

            WHEN( "Getting occurrences of the iterator" )
            {
              auto occs = getOccurrences( it );

              THEN( "It should report all occurrences indicated by the iterator" )
              {
                REQUIRE( occs.size() == 1 );
                REQUIRE( (*occs.begin()).i1 == 2 );
                REQUIRE( (*occs.begin()).i2 == 33 );
              }
            }

            WHEN( "Going right" )
            {
              THEN( "It cannot go right" )
              {
                REQUIRE( !goRight( it ) );
                REQUIRE( representative( it ) == "nan" );
              }
            }

            WHEN( "Going up then go down" )
            {
              REQUIRE( goUp( it ) );
              REQUIRE( goDown( it, 'c' ) );

              THEN( "It should point to corresponding node in the virtual prefix tree" )
              {
                REQUIRE( repLength( it ) == 3 );
                REQUIRE( representative( it ) == "can" );
                REQUIRE( parentEdgeLength( it ) == 1 );
                REQUIRE( parentEdgeLabel( it ) == "c" );
              }

              WHEN( "Getting occurrences of the iterator" )
              {
                auto occs = getOccurrences( it );

                THEN( "It should report all occurrences indicated by the iterator" )
                {
                  REQUIRE( occs.size() == 1 );
                  REQUIRE( (*occs.begin()).i1 == 2 );
                  REQUIRE( (*occs.begin()).i2 == 9 );
                }
              }

              WHEN( "Going right" )
              {
                REQUIRE( goRight( it ) );

                THEN( "It should point to corresponding node in the virtual prefix tree" )
                {
                  REQUIRE( repLength( it ) == 3 );
                  REQUIRE( representative( it ) == "ian" );
                  REQUIRE( parentEdgeLength( it ) == 1 );
                  REQUIRE( parentEdgeLabel( it ) == "i" );
                }

                WHEN( "Getting occurrences of the iterator" )
                {
                  std::set< index_type::pos_type > truth = { { 0, 12 }, { 1, 14 } };
                  auto occs = getOccurrences( it );

                  THEN( "It should report all occurrences indicated by the iterator" )
                  {
                    REQUIRE( occs.size() == 2 );
                    for ( auto it = occs.begin(); it != occs.end(); ++it ) {
                      REQUIRE( truth.find( *it ) != truth.end() );
                    }
                  }
                }
              }
            }

            WHEN( "Go root by calling goRoot interface function" )
            {
              goRoot( it );

              THEN( "It should point root node" )
              {
                REQUIRE( isRoot( it ) );
                REQUIRE( !goUp( it ) );
                REQUIRE( !goRight( it ) );
                REQUIRE( parentEdgeLength( it ) == 0 );
                REQUIRE( parentEdgeLabel( it ) == "" );
              }
            }

            WHEN( "Go root by going up" )
            {
              goUp( it );
              goUp( it );
              goUp( it );

              THEN( "It should point root node" )
              {
                REQUIRE( isRoot( it ) );
                REQUIRE( !goUp( it ) );
                REQUIRE( !goRight( it ) );
                REQUIRE( parentEdgeLength( it ) == 0 );
                REQUIRE( parentEdgeLabel( it ) == "" );
              }
            }
          }
        }
      }
    }
  }
}

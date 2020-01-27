/**
 *    @file  pathset.h
 *   @brief  PathSet template class definition.
 *
 *  This header files contains PathSet template class definition and its interface
 *  functions.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Wed Sep 13, 2017  01:07
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  PATHSET_H__
#define  PATHSET_H__

#include <vector>

#include "vargraph.h"
#include "sequence.h"
#include "index.h"
#include "utils.h"


namespace grem {
  /**
   *  @brief  Represent a set of paths.
   *
   *  This class encapsulate the data structures to represent a set of path, its strings
   *  set, and its index. It also provides load/save functionalities.
   */
  template< typename TGraph, typename TText, typename TIndexSpec >
    class PathSet {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef seqan::StringSet< TText, seqan::Owner<> > TStringSet;
        typedef seqan::Index< TStringSet, TIndexSpec > TIndex;
        typedef uint64_t size_type;    /* The max size type can be (de)serialized now. */
        /* ====================  DATA MEMBERS  ======================================= */
        TStringSet string_set;
        TIndex index;
        std::vector< Path< TGraph > > paths_set;
        /* ====================  LIFECYCLE     ======================================= */
        PathSet( ) = default;
        /* ====================  METHODS       ======================================= */
        /**
         *  @brief  Load the set of paths and its index from file.
         *
         *  @param[in]   filepath_prefix The file path prefix of the saved set of paths.
         *  @param[in]   vargraph A pointer to the graph whose paths is loading.
         *  @return `true` if the set of paths are successfully loaded from file;
         *          otherwise `false`.
         *
         *  It first loads the index from the file. Then, the nodes and string sets will
         *  be loaded.
         */
          inline bool
        load( const std::string& filepath_prefix, const TGraph* vargraph )
        {
          clear( this->string_set );
          clear( this->index );

          if ( open( this->index, filepath_prefix.c_str() ) &&
              this->load_paths_set( filepath_prefix + "_paths", vargraph ) ) {
            this->string_set = indexText( this->index );
            return true;
          }
          return false;
        }  /* -----  end of method load  ----- */

        /**
         *  @brief  Save the set of paths and its index into file.
         *
         *  @param[in]   filepath_prefix The file path prefix to save the set of paths.
         *  @return `true` if the set of paths are successfully saved into file;
         *          otherwise `false`.
         *
         *  It first saves the index into the file. Then, the nodes sets will be written.
         */
          inline bool
        save( const std::string& filepath_prefix )
        {
          if ( seqan::save( this->index, filepath_prefix.c_str() ) &&
              this->save_paths_set( filepath_prefix + "_paths" ) ) {
            return true;
          }
          return false;
        }  /* -----  end of method save  ----- */

        /**
         *  @brief  Add a path to the set.
         *
         *  @param  new_path The new path to be added.
         *
         *  It appends the new path to the paths set, adds its string representation
         *  to the string set, and creates string set index.
         */
          inline void
        add_path( Path< TGraph >&& new_path )
        {
          TText path_str( sequence( new_path ) );
          // :TODO:Mon Mar 06 13:00:\@cartoonist: faked quality score.
          char fake_qual = 'I';
          assignQualities( path_str, std::string( length( path_str ), fake_qual ) );
          appendValue( this->string_set, path_str );
          this->index = TIndex( this->string_set );
          initialize( new_path );
          paths_set.push_back( std::move( new_path ) );
        }  /* -----  end of method add_path  ----- */

        /**
         *  @brief  Add a path to the set.
         *
         *  @param  new_path The new path to be added.
         *
         *  Overloaded.
         */
          inline void
        add_path( const Path< TGraph >& new_path )
        {
          this->add_path( Path< TGraph >( new_path ) );
        }

        /**
         *  @brief  Get the size of the paths set.
         *
         *  @return The size of the paths set.
         *
         *  It returns the size of the vector of `Path` objects.
         */
          inline size_type
        size( ) const
        {
          return this->paths_set.size();
        }  /* -----  end of method size  ----- */

        /**
         *  @brief  Reserve memory for the set of paths.
         *
         *  @param  size The target capacity.
         *
         *  It reserves memory for paths set and strings set.
         */
          void
        reserve( size_type size )
        {
          seqan::reserve( this->string_set, size, seqan::Exact() );
          this->paths_set.reserve( size );
        }  /* -----  end of method reserve  ----- */

        /**
         *  @brief  Create index fibres.
         *
         *  Initializing the index does not create index fibres. Index fibres are
         *  generated on-demand. This function forces to create these fibres in
         *  advance.
         */
          inline void
        create_index( )
        {
          grem::create_index( this->index );
        }  /* -----  end of method create_index  ----- */
      private:
        /**
         *  @brief  Load nodes sets of paths set from file.
         *
         *  @param[in]   filepath The file path of the saved paths.
         *  @param[in]   vargraph A pointer to the graph whose paths is loading.
         *  @return `true` if the nodes sets are successfully loaded from file; otherwise
         *          `false`.
         *
         *  It deserializes the paths set. The paths are prefixed by the number of paths.
         */
          inline bool
        load_paths_set( const std::string& filepath, const TGraph* vargraph )
        {
          std::ifstream ifs( filepath, std::ifstream::in | std::ifstream::binary );
          if ( !ifs ) return false;
          size_type path_num = 0;

          try {
            deserialize( ifs, path_num );
            this->paths_set.clear();
            this->paths_set.reserve( path_num );
            for ( unsigned int i = 0; i < path_num; ++i ) {
              Path< TGraph > path( vargraph );
              grem::load( path, ifs );
              this->paths_set.push_back( std::move( path ) );
            }
          }
          catch ( const std::runtime_error& ) {
            this->paths_set.clear();
            return false;
          }

          return true;
        }  /* -----  end of function load_paths_set  ----- */

        /**
         *  @brief  Save nodes sets of paths set into file.
         *
         *  @param[in]   filepath The file path of the file to be written.
         *  @return `true` if the nodes sets are successfully written to file; otherwise
         *          `false`.
         *
         *  It serializes the paths to the files followed by the number of paths.
         */
          inline bool
        save_paths_set( const std::string& filepath )
        {
          std::ofstream ofs( filepath, std::ofstream::out | std::ofstream::binary );
          if( !ofs ) return false;

          try {
            serialize( ofs, static_cast< size_type >( this->paths_set.size() ) );
            for ( const auto& path : this->paths_set ) {
              grem::save( path, ofs );
            }
          }
          catch ( const std::runtime_error& ) {
            return false;
          }

          return true;
        }  /* -----  end of function save_paths_set  ----- */
    };  /* -----  end of template class PathSet  ----- */

  /* Typedefs  ----------------------------------------------------------------- */

  template< typename TGraph, typename TIndexSpec >
    using Dna5QPathSet = PathSet< TGraph, seqan::Dna5QString, TIndexSpec >;

}  /* -----  end of namespace grem  ----- */

#endif  // end of PATHSET_H__

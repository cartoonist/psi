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
  template< typename TIndexSpec >
    class PathSet {
      public:
        /* ====================  DATA MEMBERS  ======================================= */
        Dna5QStringSet string_set;
        Dna5QStringSetIndex< TIndexSpec > index;
        std::vector< Path<> > paths_set;
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename decltype( paths_set )::size_type size_type;
        /* ====================  LIFECYCLE     ======================================= */
        PathSet( ) = default;
        /* ====================  METHODS       ======================================= */
        /**
         *  @brief  Load the set of paths and its index from file.
         *
         *  @param[in]   filepath_prefix The file path prefix of the saved set of paths.
         *  @param[in]   path_num The number of paths.
         *  @return `true` if the set of paths are successfully loaded from file;
         *          otherwise `false`.
         *
         *  It first loads the index from the file. Then, the nodes and string sets will
         *  be loaded. The number of saved paths and requested number of paths should be
         *  the same.
         */
          inline bool
        load( const std::string& filepath_prefix, const VarGraph* vargraph,
            unsigned int path_num )
        {
          if ( open( this->index, filepath_prefix.c_str() ) &&
              this->load_paths_set( filepath_prefix + "_path_", vargraph, path_num ) ) {
            this->string_set = indexText( this->index );
            const auto& nof_loaded_paths = length( this->string_set );
            if ( nof_loaded_paths == path_num ) {
              return true;
            }
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
          if ( !writable( filepath_prefix ) ) return false;
          seqan::save( this->index, filepath_prefix.c_str() );
          this->save_paths_set( filepath_prefix + "_path_" );
          return true;
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
        add_path( Path<>&& new_path )
        {
          typedef typename seqan::Value< Dna5QStringSet >::Type TText;
          TText path_str( sequence( new_path ) );
          // :TODO:Mon Mar 06 13:00:\@cartoonist: faked quality score.
          char fake_qual = 'I';
          assignQualities( path_str, std::string( length( path_str ), fake_qual ) );
          appendValue( this->string_set, path_str );
          this->index = Dna5QStringSetIndex< TIndexSpec >( this->string_set );
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
        add_path( const Path<>& new_path )
        {
          this->add_path( Path<>( new_path ) );
        }

        /**
         *  @brief  Get the size of the paths set.
         *
         *  @return The size of the paths set.
         *
         *  It returns the size of the vector of `Path` objects.
         */
          inline size_type
        size ( )
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
         *  @param[in]   filepath_prefix The file path prefix of the saved paths.
         *  @param[in]   path_num The number of paths.
         *  @return `true` if the nodes sets are successfully loaded from file; otherwise
         *          `false`.
         *
         *  It deserializes the node IDs saved in the files whose names can be obtained by
         *  appending the path number to the given path prefix.
         */
          inline bool
        load_paths_set( const std::string& filepath_prefix, const VarGraph* vargraph,
            unsigned int path_num )
        {
          this->paths_set.clear();
          this->paths_set.reserve( path_num );
          for ( unsigned int i = 0; i < path_num; ++i ) {
            std::string file_path = filepath_prefix + std::to_string( i );
            Path<> path( vargraph );
            try {
              grem::load( path, file_path );
            } catch ( const std::runtime_error& ) {
              this->paths_set.clear();
              return false;
            }
            this->paths_set.push_back( std::move( path ) );
          }

          return true;
        }  /* -----  end of function load_paths_set  ----- */

        /**
         *  @brief  Save nodes sets of paths set into file.
         *
         *  @param[in]   filepath_prefix The file path prefix of the file to be written.
         *
         *  It serializes the node IDs of the paths to the files - a file for each path -
         *  whose names can be obtained by appending the path number to the given path prefix.
         */
          inline void
        save_paths_set( const std::string& filepath_prefix )
        {
          unsigned int i = 0;
          for ( const auto& path : this->paths_set ) {
            std::string file_path = filepath_prefix + std::to_string ( i );
            grem::save( path, file_path );
            ++i;
          }
        }  /* -----  end of function save_paths_set  ----- */
    };  /* -----  end of template class PathSet  ----- */
}  /* -----  end of namespace grem  ----- */

#endif  // end of PATHSET_H__

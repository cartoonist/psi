/**
 *    @file  pathindex.h
 *   @brief  PathIndex template class definition.
 *
 *  This header files contains PathIndex template class definition and its interface
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

#ifndef  PATHINDEX_H__
#define  PATHINDEX_H__

#include <vector>
#include <algorithm>

#include "vargraph.h"
#include "sequence.h"
#include "pathset.h"
#include "index.h"
#include "utils.h"


namespace grem {
  /**
   *  @brief  Represent a path index.
   *
   *  This class encapsulate the data structures to represent an index of a set of
   *  paths, its strings set. It also provides load/save functionalities.
   */
  template< typename TGraph, typename TText, typename TIndexSpec, typename TSequenceDirection = Forward >
    class PathIndex {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef seqan::StringSet< TText, seqan::Owner<> > TStringSet;
        typedef seqan::Index< TStringSet, TIndexSpec > TIndex;
        typedef Path< TGraph, Compact > TPath;
        typedef uint64_t size_type;    /* The max size type can be (de)serialized now. */
        typedef uint64_t context_type;
        typedef TPath value_type;
        /* ====================  DATA MEMBERS  ======================================= */
        TIndex index;
        std::vector< TPath > paths_set;
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        TStringSet string_set;
        /**< @brief The extreme nodes' sequence will be trimmed to the length of `context-1`. */
        context_type context;
        bool lazy_mode;
      public:
        /* ====================  LIFECYCLE     ======================================= */
        PathIndex( context_type ct=0, bool l=false )
          : context( ct ), lazy_mode( l ) { };
        PathIndex( bool lazy ) : context( 0 ), lazy_mode( lazy ) { };
        /* ====================  METHODS       ======================================= */
        /**
         *  @brief  Get the shift required to add to trimmed position to get original one.
         *
         *  @param  pos Position in the path set.
         *  @return the shift required to add to the position in the trimmed path in
         *          order to get the position in the original path.
         *
         *  The parameter `pos` indicates one path in the path set. If `context` is
         *  non-zero the path sequence in the string set is left-trimmed by the length
         *  of `length( sequence( path.nodes[0] ) ) - context + 1`. So position found in
         *  the trimmed path sequence needs to be converted to the original position in
         *  the path by adding this amount.
         */
          inline unsigned int
        get_context_shift( typename seqan::SAValue< TIndex >::Type const& pos ) const
        {
          assert( pos.i1 < this->paths_set.size() );
          auto&& the_path = this->paths_set[ pos.i1 ];
          auto flen = the_path.get_vargraph()->node_length( the_path.get_nodes()[0] );
          if ( this->context != 0 && flen + 1 > this->context ) {
            return flen - this->context + 1;
          }
          return 0;
        }  /* -----  end of method get_context_shift  ----- */

        /**
         *  @brief  Load the path index from file.
         *
         *  @param[in]   filepath_prefix The file path prefix of the saved path index.
         *  @param[in]   vargraph A pointer to the graph whose paths is loading.
         *  @return `true` if the path index is successfully loaded from file; otherwise
         *          `false`.
         *
         *  It first loads the index from the file, then the paths set and other
         *  attributes.
         */
          inline bool
        load( const std::string& filepath_prefix, const TGraph* vargraph )
        {
          clear( this->string_set );
          clear( this->index );

          if ( open( this->index, filepath_prefix ) &&
              this->load_paths_set( filepath_prefix + "_paths", vargraph ) ) {
            //this->string_set = indexText( this->index );  /* XXX: Unnecessary copy. Use `FibreText` of `this->index` directly. */
            return true;
          }
          return false;
        }  /* -----  end of method load  ----- */

        /**
         *  @brief  Save the path index into file.
         *
         *  @param[in]   filepath_prefix The file path prefix to save the path index.
         *  @return `true` if the path index is successfully saved into file; otherwise
         *          `false`.
         *
         *  It first saves the index into the file, then the paths set and other
         *  attributes.
         */
          inline bool
        serialize( const std::string& filepath_prefix )
        {
          if ( save( this->index, filepath_prefix ) &&
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
         *  to the string set, and creates string set index if it is not in lazy mode.
         */
          inline void
        add_path( TPath new_path )
        {
          this->paths_set.push_back( std::move( new_path ) );
          initialize( this->paths_set.back() );
          if ( !this->lazy_mode ) {
            this->add_path_sequence( std::prev(this->paths_set.end()), this->paths_set.end() );
          }
        }  /* -----  end of method add_path  ----- */

        /**
         *  @brief  Add a path of different type to the set.
         *
         *  @param  new_path The new path to be added of different type.
         *
         *  It appends the new path to the paths set, adds its string representation
         *  to the string set, and creates string set index if it is not in lazy mode.
         */
        template< typename TPathSpec >
            inline void
          add_path( const Path< VarGraph, TPathSpec >& new_path,
              enable_if_not_equal_t< typename TPath::spec_type, TPathSpec > tag = TPathSpec() )
          {
            TPath native_path( new_path.get_vargraph() );
            native_path = new_path;
            this->add_path( std::move( native_path ) );
          }

        /**
         *  @brief  Add a path of different type to the set.
         *
         *  @param  new_path The new path to be added of different type.
         *
         *  It appends the new path to the paths set, adds its string representation
         *  to the string set, and creates string set index if it is not in lazy mode.
         */
        template< typename TPathSpec >
            inline void
          add_path( Path< VarGraph, TPathSpec >&& new_path,
              enable_if_not_equal_t< typename TPath::spec_type, TPathSpec > tag = TPathSpec() )
          {
            TPath native_path( new_path.get_vargraph() );
            native_path = std::move( new_path );
            this->add_path( std::move( native_path ) );
          }

        /**
         *  @brief  alias: See `add_path`.
         */
          inline void
        push_back( TPath new_path )
        {
          this->add_path( std::move( new_path ) );
        }

        /**
         *  @brief  alias: See `add_path`.
         */
        template< typename TPathSpec >
          inline void
        push_back( const Path< VarGraph, TPathSpec >& new_path )
        {
          this->add_path( new_path );
        }

        /**
         *  @brief  alias: See `add_path`.
         */
        template< typename TPathSpec >
          inline void
        push_back( Path< VarGraph, TPathSpec >&& new_path )
        {
          this->add_path( std::move( new_path ) );
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
         *  @brief  Reserve memory for the set of paths and corresponding string set.
         *
         *  @param  size The target capacity.
         *
         *  It reserves memory for paths set and strings set.
         */
          void
        reserve( size_type size )
        {
          seqan::reserve( this->string_set, size );
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
          if ( lazy_mode ) {
            this->add_path_sequence( this->paths_set.begin(), this->paths_set.end() );
          }
          grem::create_index( this->index );
        }  /* -----  end of method create_index  ----- */
      private:
        /**
         *  @brief  Add sequence of the paths in the set from `begin` to `end`.
         *
         *  @param  begin The begin iterator of the paths set to add.
         *  @param  end The end iterator of the paths set to add.
         *
         *  It appends the sequence of the paths in the range `[begin, end)` to the
         *  string set, and update string set index.
         *
         *  @note The quality score of the paths are considered `I`.
         */
        template< typename TIter >
            inline void
          add_path_sequence( TIter begin, TIter end )
          {
            for ( ; begin != end; ++begin ) {
              //TText path_str( sequence( *begin, TSequenceDirection(), this->context ) );
              // :TODO:Mon Mar 06 13:00:\@cartoonist: quality score?
              //char fake_qual = 'I';
              //assignQualities( path_str, std::string( length( path_str ), fake_qual ) );
              appendValue( this->string_set,
                  sequence( *begin, TSequenceDirection(), this->context ) );
            }
            this->index = TIndex( this->string_set );
          }

        /**
         *  @brief  Load the set of paths and other attributes from file.
         *
         *  @param[in]   filepath The file path of the saved paths.
         *  @param[in]   vargraph A pointer to the graph whose paths is loading.
         *  @return `true` if the paths set is successfully loaded from file; otherwise
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
            deserialize( ifs, this->context );
            deserialize( ifs, path_num );
            this->paths_set.clear();
            this->paths_set.reserve( path_num );
            for ( unsigned int i = 0; i < path_num; ++i ) {
              TPath path( vargraph );
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
         *  @brief  Save set of paths and other attributes into file.
         *
         *  @param[in]   filepath The file path of the file to be written.
         *  @return `true` if the paths set is successfully written to file; otherwise
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
            grem::serialize( ofs, this->context );
            grem::serialize( ofs, static_cast< size_type >( this->paths_set.size() ) );
            for ( const auto& path : this->paths_set ) {
              grem::save( path, ofs );
            }
          }
          catch ( const std::runtime_error& ) {
            return false;
          }

          return true;
        }  /* -----  end of function save_paths_set  ----- */
    };  /* -----  end of template class PathIndex  ----- */

  /* Typedefs  ----------------------------------------------------------------- */

  template< typename TGraph, typename TIndexSpec, typename TSequenceDirection = Forward >
    using Dna5QPathIndex = PathIndex< TGraph, seqan::Dna5QString, TIndexSpec, TSequenceDirection >;

  /**
   *  @brief  Meta-function getting direction of paths in a `PathIndex`.
   */
  template< typename TGraph, typename TText, typename TIndexSpec, typename TSequenceDirection >
    struct Direction< PathIndex< TGraph, TText, TIndexSpec, TSequenceDirection > > {
      typedef TSequenceDirection Type;
    };

  /* PathIndex interface functions  ------------------------------------------------ */

  template< typename TGraph, typename TText, typename TIndexSpec >
      inline typename PathIndex< TGraph, TText, TIndexSpec >::size_type
    length( PathIndex< TGraph, TText, TIndexSpec >& pindex )
    {
      return pindex.size();
    }

  template< typename TIndex >
    using TSAValue = typename seqan::SAValue< TIndex >::Type;

  template< typename TGraph, typename TText, typename TIndexSpec >
      inline typename TGraph::offset_type
    position_to_offset( PathIndex< TGraph, TText, TIndexSpec, Forward > const& pindex,
        TSAValue< typename PathIndex< TGraph, TText, TIndexSpec, Forward >::TIndex > const& pos )
    {
      auto context_shift = pindex.get_context_shift( pos );
      assert( pos.i1 < pindex.paths_set.size() );
      return position_to_offset( pindex.paths_set[ pos.i1 ], context_shift + pos.i2 );
    }

  /**
   *  @note The input `pos` parameter should be the "end position" of the occurrence; e.g.:
   *
   *  The pattern 'ttc' is found in the reversed string:
   *        0123 456 7890123
   *        acga ctt taggtcc
   *  the input `pos` parameter should be 6 (not 4) where the `real_pos` computed inside
   *  the function would be the real position in forward sequence.
   */
  template< typename TGraph, typename TText, typename TIndexSpec >
      inline typename TGraph::offset_type
    position_to_offset( PathIndex< TGraph, TText, TIndexSpec, Reversed > const& pindex,
        TSAValue< typename PathIndex< TGraph, TText, TIndexSpec, Reversed >::TIndex > const& pos )
    {
      auto real_pos = pos;
      real_pos.i2 = length( indexText( pindex.index )[ pos.i1 ] ) - pos.i2 - 1;
      auto context_shift = pindex.get_context_shift( real_pos );
      assert( real_pos.i1 < pindex.paths_set.size() );
      return position_to_offset( pindex.paths_set[ real_pos.i1 ], context_shift + real_pos.i2 );
    }

  template< typename TGraph, typename TText, typename TIndexSpec >
      inline typename TGraph::nodeid_type
    position_to_id( PathIndex< TGraph, TText, TIndexSpec, Forward > const& pindex,
        TSAValue< typename PathIndex< TGraph, TText, TIndexSpec, Forward >::TIndex > const& pos )
    {
      auto context_shift = pindex.get_context_shift( pos );
      assert( pos.i1 < pindex.paths_set.size() );
      return position_to_id( pindex.paths_set[ pos.i1 ], context_shift + pos.i2 );
    }

  /**
   *  @note The input `pos` parameter should be the "end position" of the occurrence; e.g.:
   *
   *  The pattern 'ttc' is found in the reversed string:
   *        0123 456 7890123
   *        acga ctt taggtcc
   *  the input `pos` parameter should be 6 (not 4) where the `real_pos` computed inside
   *  the function would be the real position in forward sequence.
   */
  template< typename TGraph, typename TText, typename TIndexSpec >
      inline typename TGraph::nodeid_type
    position_to_id( PathIndex< TGraph, TText, TIndexSpec, Reversed > const& pindex,
        TSAValue< typename PathIndex< TGraph, TText, TIndexSpec, Reversed >::TIndex > const& pos )
    {
      auto real_pos = pos;
      real_pos.i2 = length( indexText( pindex.index )[ pos.i1 ] ) - pos.i2 - 1;
      auto context_shift = pindex.get_context_shift( real_pos );
      assert( real_pos.i1 < pindex.paths_set.size() );
      return position_to_id( pindex.paths_set[ real_pos.i1 ], context_shift + real_pos.i2 );
    }

  /* END OF PathIndex interface functions  ----------------------------------------- */

}  /* -----  end of namespace grem  ----- */

#endif  // end of PATHINDEX_H__

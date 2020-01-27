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
#include <algorithm>

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
  template< typename TGraph, typename TText, typename TIndexSpec, typename TSequenceDirection = Forward >
    class PathSet {
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
        bool sorted;
      public:
        /* ====================  LIFECYCLE     ======================================= */
        PathSet( context_type ct=0, bool l=false )
          : context( ct ), lazy_mode( l ), sorted( false ) { };
        PathSet( bool lazy ) : context( 0 ), lazy_mode( lazy ), sorted( false ) { };
        /* ====================  ACCESSORS     ======================================= */
        /**
         *  @brief  return true if the paths in the set is sorted by their minimum node ID.
         */
          inline bool
        is_sorted() const
        {
          return this->sorted;
        }  /* -----  end of method is_sorted  ----- */
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
          auto&& the_path = this->paths_set.at( pos.i1 );
          auto flen = the_path.get_vargraph()->node_length( the_path.get_nodes()[0] );
          if ( this->context != 0 && flen + 1 > this->context ) {
            return flen - this->context + 1;
          }
          return 0;
        }  /* -----  end of method get_context_shift  ----- */

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

          if ( open( this->index, filepath_prefix ) &&
              this->load_paths_set( filepath_prefix + "_paths", vargraph ) ) {
            //this->string_set = indexText( this->index );  /* XXX: Unnecessary copy. Use `FibreText` of `this->index` directly. */
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
          this->sorted = false;
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
         *  @brief  Reserve memory for the set of paths.
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

        /**
         *  @brief  Sort the paths set by the minimum node ID in each path.
         *
         *  Sorting the paths in ascending order of minimum node ID in the paths. This
         *  helps the `covered_by` query to search in the paths local to the query path.
         */
          inline void
        sort()
        {
          auto compare = []( TPath& a, TPath& b ) {
            return *a.get_nodes_set().begin() < *b.get_nodes_set().begin();
          };

          std::sort( this->paths_set.begin(), this->paths_set.end(), compare );
          this->sorted = true;
        }
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
          size_type sort_int = 0;

          try {
            deserialize( ifs, this->context );
            deserialize( ifs, sort_int );
            this->sorted = static_cast< bool >( sort_int );
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
            grem::serialize( ofs, this->context );
            grem::serialize( ofs, static_cast< size_type >( this->sorted ) );
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
    };  /* -----  end of template class PathSet  ----- */

  /* Typedefs  ----------------------------------------------------------------- */

  template< typename TGraph, typename TIndexSpec, typename TSequenceDirection = Forward >
    using Dna5QPathSet = PathSet< TGraph, seqan::Dna5QString, TIndexSpec, TSequenceDirection >;

  /**
   *  @brief  Meta-function getting direction of paths in a `PathSet`.
   */
  template< typename TGraph, typename TText, typename TIndexSpec, typename TSequenceDirection >
    struct Direction< PathSet< TGraph, TText, TIndexSpec, TSequenceDirection > > {
      typedef TSequenceDirection Type;
    };

  /* PathSet interface functions  ------------------------------------------------ */

  template< typename TGraph, typename TText, typename TIndexSpec >
      inline typename PathSet< TGraph, TText, TIndexSpec >::size_type
    length( PathSet< TGraph, TText, TIndexSpec >& set )
    {
      return set.size();
    }

  template< typename TIndex >
    using TSAValue = typename seqan::SAValue< TIndex >::Type;

  template< typename TGraph, typename TText, typename TIndexSpec >
      inline typename TGraph::offset_type
    position_to_offset( PathSet< TGraph, TText, TIndexSpec, Forward > const& set,
        TSAValue< typename PathSet< TGraph, TText, TIndexSpec, Forward >::TIndex > const& pos )
    {
      auto context_shift = set.get_context_shift( pos );
      return position_to_offset( set.paths_set.at( pos.i1 ), context_shift + pos.i2 );
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
    position_to_offset( PathSet< TGraph, TText, TIndexSpec, Reversed > const& set,
        TSAValue< typename PathSet< TGraph, TText, TIndexSpec, Reversed >::TIndex > const& pos )
    {
      auto real_pos = pos;
      real_pos.i2 = length( indexText( set.index )[ pos.i1 ] ) - pos.i2 - 1;
      auto context_shift = set.get_context_shift( real_pos );
      return position_to_offset( set.paths_set.at( real_pos.i1 ), context_shift + real_pos.i2 );
    }

  template< typename TGraph, typename TText, typename TIndexSpec >
      inline typename TGraph::nodeid_type
    position_to_id( PathSet< TGraph, TText, TIndexSpec, Forward > const& set,
        TSAValue< typename PathSet< TGraph, TText, TIndexSpec, Forward >::TIndex > const& pos )
    {
      auto context_shift = set.get_context_shift( pos );
      return position_to_id( set.paths_set.at( pos.i1 ), context_shift + pos.i2 );
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
    position_to_id( PathSet< TGraph, TText, TIndexSpec, Reversed > const& set,
        TSAValue< typename PathSet< TGraph, TText, TIndexSpec, Reversed >::TIndex > const& pos )
    {
      auto real_pos = pos;
      real_pos.i2 = length( indexText( set.index )[ pos.i1 ] ) - pos.i2 - 1;
      auto context_shift = set.get_context_shift( real_pos );
      return position_to_id( set.paths_set.at( real_pos.i1 ), context_shift + real_pos.i2 );
    }

  /**
   *  @brief  Check whether a path is covered by a PathSet.
   *
   *  @param  path The given path as Path class instance.
   *  @param  paths_set A set of paths as PathSet class instance.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  Overloaded. See `covered_by( TIter1, TIter1, TIter2, TIter2 )`.
   */
  template< typename TGraph, typename TPathSpec, typename TText, typename TIndexSpec, typename TSequenceDirection, typename TStrategy >
      inline bool
    covered_by( const Path< TGraph, TPathSpec >& path,
        const PathSet< TGraph, TText, TIndexSpec, TSequenceDirection >& set, TStrategy )
    {
      typedef typename PathSet< TGraph, TText, TIndexSpec, TSequenceDirection >::TPath TPath;

      /* Return `true` if the min node ID of the query path is less than the min node ID
       * of the element path */
      auto mn = []( Path< TGraph, TPathSpec > const& query, TPath const& element ) {
        return *query.get_nodes_set().begin() < *element.get_nodes_set().begin();
      };
      /* Return `true` if the max node ID of the element path is less than the max node
       * ID of the query path */
      auto mx = []( TPath const& element, Path< TGraph, TPathSpec > const& query ) {
        return *element.get_nodes_set().rbegin() < *query.get_nodes_set().rbegin();
      };

      if ( set.is_sorted() ) {
        /* Find the lower bound and upper bound of the search. */
        auto lbound = set.paths_set.begin();
        while ( lbound != set.paths_set.end() && mx( *lbound, path ) ) ++lbound;
        auto ubound = std::upper_bound( set.paths_set.begin(), set.paths_set.end(), path, mn );

        return lbound < ubound &&
          covered_by( path.get_nodes().begin(), path.get_nodes().end(), lbound, ubound, TStrategy() );
      }
      else {
        return covered_by( path.get_nodes().begin(), path.get_nodes().end(),
              set.paths_set.begin(), set.paths_set.end(), TStrategy() );
      }
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Check whether a path is covered by a PathSet.
   *
   *  @param  path The given path as Path class instance.
   *  @param  paths_set A set of paths as PathSet class instance.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  Overloaded. with no strategy specified it calls with Default strategy.
   */
  template< typename TGraph, typename TPathSpec, typename TText, typename TIndexSpec, typename TSequenceDirection >
      inline bool
    covered_by( const Path< TGraph, TPathSpec >& path,
        const PathSet< TGraph, TText, TIndexSpec, TSequenceDirection >& set )
    {
      return covered_by( path, set, Default() );
    }

  template< typename TGraph, typename TText, typename TIndexSpec, typename TSequenceDirection, typename TPathSet >
      inline void
    compress( PathSet< TGraph, TText, TIndexSpec, TSequenceDirection > const& set, TPathSet& out )
    {
      typedef typename TPathSet::value_type TPath;

      auto pre_itr = set.paths_set.begin();
      auto cur_itr = pre_itr + 1;
      auto vargraph = (*pre_itr).get_vargraph();

      if ( pre_itr == set.paths_set.end() ) return;

      Path< VarGraph > path( vargraph );
      path = *pre_itr;
      while ( cur_itr != set.paths_set.end() ) {
        if ( *( (*pre_itr).get_nodes().end() - 1 ) > *(*cur_itr).get_nodes().begin() ) {
          TPath final_path( vargraph );
          final_path = std::move( path );
          out.push_back( std::move( final_path ) );
          clear( path );
        }
        path += *cur_itr;
        pre_itr = cur_itr;
        ++cur_itr;
      }

      TPath final_path( vargraph );
      final_path = std::move( path );
      out.push_back( std::move( final_path ) );
    }

  /* END OF PathSet interface functions  ----------------------------------------- */

}  /* -----  end of namespace grem  ----- */

#endif  // end of PATHSET_H__

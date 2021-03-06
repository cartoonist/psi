/**
 *    @file  pathindex.hpp
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

#ifndef  PSI_PATHINDEX_HPP__
#define  PSI_PATHINDEX_HPP__

#include <vector>
#include <algorithm>

#include "sequence.hpp"
#include "index.hpp"
#include "pathset.hpp"
#include "graph_iter.hpp"
#include "utils.hpp"


namespace psi {
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
        typedef TGraph graph_type;
        typedef Path< graph_type, Compact > value_type;
        typedef PathSet< value_type > container_type;
        typedef seqan::StringSet< TText, seqan::Owner<> > stringset_type;
        typedef seqan::Index< stringset_type, TIndexSpec > index_type;
        typedef uint64_t size_type;    /* The max size type can be (de)serialized now. */
        typedef uint64_t context_type;
        /* ====================  DATA MEMBERS  ======================================= */
        index_type index;
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        container_type paths_set;
        stringset_type string_set;
        context_type context;
        bool lazy_mode;
      public:
        /* ====================  LIFECYCLE     ======================================= */
        PathIndex( graph_type const& graph,  context_type ct=0, bool l=false )
          : paths_set( graph ), context( ct ), lazy_mode( l ) { };

        PathIndex( graph_type const& graph, bool lazy )
          : paths_set( graph ), context( 0 ), lazy_mode( lazy ) { };
        /* ====================  ACCESSORS     ======================================= */
          inline container_type&
        get_paths_set( )
        {
          return this->paths_set;
        }

          inline container_type const&
        get_paths_set( ) const
        {
          return this->paths_set;
        }

          inline context_type
        get_context( ) const
        {
          return this->context;
        }
        /* ====================  MUTATORS      ======================================= */
          inline void
        set_context( context_type value )
        {
          this->context = value;
        }
        /* ====================  METHODS       ======================================= */
          inline void
        clear( )
        {
          seqan::clear( this->index );
          this->paths_set.clear();
          seqan::clear( this->string_set );
        }

        /**
         *  @brief  Load the path index from file.
         *
         *  @param[in]   filepath_prefix The file path prefix of the saved path index.
         *  @return `true` if the path index is successfully loaded from file; otherwise
         *          `false`.
         *
         *  It first loads the index from the file, then the paths set and other
         *  attributes.
         */
          inline bool
        load( const std::string& filepath_prefix )
        {
          this->clear();

          if ( open( this->index, filepath_prefix ) &&
              this->load_paths_set( filepath_prefix + "_paths" ) ) {
            /* XXX: Unnecessary copy. Use `FibreText` of `this->index` directly. */
            //this->string_set = indexText( this->index );
            return true;
          }

          this->clear();
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
        add_path( value_type new_path )
        {
          this->paths_set.push_back( std::move( new_path ) );
          if ( !this->lazy_mode ) {
            this->add_path_sequence( std::prev( this->paths_set.end() ),
                this->paths_set.end() );
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
          add_path( Path< graph_type, TPathSpec > new_path,
              enable_if_not_equal_t< typename value_type::spec_type, TPathSpec > tag = TPathSpec() )
          {
            value_type native_path( new_path.get_graph_ptr() );
            native_path = std::move( new_path );
            this->add_path( std::move( native_path ) );
          }

        /**
         *  @brief  alias: See `add_path`.
         */
          inline void
        push_back( value_type new_path )
        {
          this->add_path( std::move( new_path ) );
        }

        /**
         *  @brief  alias: See `add_path`.
         */
        template< typename TPathSpec >
            inline void
          push_back( Path< graph_type, TPathSpec > new_path,
              enable_if_not_equal_t< typename value_type::spec_type, TPathSpec > tag = TPathSpec() )
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
          inline void
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
          psi::create_index( this->index );
          this->paths_set.initialize();
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
              //TText path_str( sequence( *begin, TSequenceDirection() ) );
              // :TODO:Mon Mar 06 13:00:\@cartoonist: quality score?
              //char fake_qual = 'I';
              //assignQualities( path_str, std::string( length( path_str ), fake_qual ) );
              appendValue( this->string_set, sequence( *begin, TSequenceDirection() ) );
            }
            this->index = index_type( this->string_set );
          }

        /**
         *  @brief  Load the set of paths and other attributes from file.
         *
         *  @param[in]   filepath The file path of the saved paths.
         *  @return `true` if the paths set is successfully loaded from file; otherwise
         *          `false`.
         *
         *  It deserializes the paths set. The paths are prefixed by the number of paths.
         */
          inline bool
        load_paths_set( const std::string& filepath )
        {
          std::ifstream ifs( filepath, std::ifstream::in | std::ifstream::binary );
          if ( !ifs ) return false;

          try {
            context_type c;
            size_type dir;
            /* Context */
            deserialize( ifs, c );
            if ( this->context == 0 ) this->context = c;
            else if ( this->context != c ) return false;
            /* Direction */
            deserialize( ifs, dir );
            if ( dir != std::is_same< TSequenceDirection, Forward >::value ) {
              return false;
            }
            this->paths_set.load( ifs );
          }
          catch ( const std::runtime_error& ) {
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
            psi::serialize( ofs, this->context );
            size_type dir = std::is_same< TSequenceDirection, Forward >::value;
            psi::serialize( ofs, dir );
            this->paths_set.serialize( ofs );
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
    length( PathIndex< TGraph, TText, TIndexSpec > const& pindex )
    {
      return pindex.size();
    }

  template< typename TIndex >
    using TSAValue = typename seqan::SAValue< TIndex >::Type;

  template< typename TGraph, typename TText, typename TIndexSpec >
      inline typename TGraph::offset_type
    position_to_offset( PathIndex< TGraph, TText, TIndexSpec, Forward > const& pindex,
        TSAValue< typename PathIndex< TGraph, TText, TIndexSpec, Forward >::index_type > const& pos )
    {
      assert( pos.i1 < pindex.get_paths_set().size() );
      return position_to_offset( pindex.get_paths_set()[ pos.i1 ], pos.i2 );
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
        TSAValue< typename PathIndex< TGraph, TText, TIndexSpec, Reversed >::index_type > const& pos )
    {
      auto real_pos = pos;
      real_pos.i2 = pindex.get_paths_set()[ pos.i1 ].get_sequence_len() - pos.i2 - 1;
      assert( real_pos.i1 < pindex.get_paths_set().size() );
      return position_to_offset( pindex.get_paths_set()[ real_pos.i1 ], real_pos.i2 );
    }

  template< typename TGraph, typename TText, typename TIndexSpec >
      inline typename TGraph::id_type
    position_to_id( PathIndex< TGraph, TText, TIndexSpec, Forward > const& pindex,
        TSAValue< typename PathIndex< TGraph, TText, TIndexSpec, Forward >::index_type > const& pos )
    {
      assert( pos.i1 < pindex.get_paths_set().size() );
      return position_to_id( pindex.get_paths_set()[ pos.i1 ], pos.i2 );
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
      inline typename TGraph::id_type
    position_to_id( PathIndex< TGraph, TText, TIndexSpec, Reversed > const& pindex,
        TSAValue< typename PathIndex< TGraph, TText, TIndexSpec, Reversed >::index_type > const& pos )
    {
      auto real_pos = pos;
      real_pos.i2 = pindex.get_paths_set()[ pos.i1 ].get_sequence_len() - pos.i2 - 1;
      assert( real_pos.i1 < pindex.get_paths_set().size() );
      return position_to_id( pindex.get_paths_set()[ real_pos.i1 ], real_pos.i2 );
    }

  /**
   *  @brief  Check whether a path is covered by a PathIndex.
   *
   *  @param  path The given path as Path class instance.
   *  @param  pindex The path index as PathIndex class instance.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  It calls the `covered_by` interface function of the paths set of the path index.
   */
  template< typename TGraph, typename TPathSpec, typename TText, typename TIndexSpec, typename TSequenceDirection >
      inline bool
    covered_by( const Path< TGraph, TPathSpec >& path,
        PathIndex< TGraph, TText, TIndexSpec, TSequenceDirection >& pindex )
    {
      return covered_by( path.get_nodes().begin(), path.get_nodes().end(),
          pindex.get_paths_set() );
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Simulate a unique haplotype.
   *
   *  @param[out]  haplotype The simulated haplotype as a Path.
   *  @param[in,out]  iter Haplotyper graph iterator.
   *  @param[in]  tries Number of tries if the generated haplotype is not unique.
   *
   *  This function gets a Haplotyper graph iterator and generate a unique haplotype
   *  if available. The input Haplotyper iterator stores required information of the
   *  previous simulated haplotypes for which the iterator is used. So, in order to
   *  simulate multiple unique haplotypes use the same iterator as the input. It tries
   *  `tries` times to generated a unique haplotype.
   */
  template< class TGraph, typename TSpec >
      inline void
    get_uniq_full_haplotype( Path< TGraph >& haplotype,
        GraphIter< TGraph, Haplotyper< TSpec > >& hp_itr,
        typename GraphIter< TGraph, Haplotyper< TSpec > >::end_type& hp_end,
        int tries=0 )
    {
      do {
        while ( hp_itr != hp_end ) {
          add_node( haplotype, *hp_itr );
          ++hp_itr;
        }
        if ( tries-- && hp_itr[ haplotype.get_nodes() ] ) {
          hp_itr--;  // discard the traversed path and reset the Haplotyper iterator.
        }
        else {
          --hp_itr;  // save the traversed path and reset the Haplotyper iterator.
          break;
        }
        /* trying again */
      } while ( true );
    }

  /**
   *  @brief  Add a simulated haplotype to the given pathset.
   *
   *  @param  paths The pathset.
   *  @param[in,out]  iter Haplotyper graph iterator.
   *  @param[in]  tries Number of tries if the generated haplotype is not unique.
   *
   *  @overload for `PathSet`.
   */
  template< typename TPathSet, class TGraph, typename TSpec >
      inline void
    get_uniq_full_haplotype( TPathSet& paths,
        GraphIter< TGraph, Haplotyper< TSpec > >& hp_itr,
        typename GraphIter< TGraph, Haplotyper< TSpec > >::end_type& hp_end,
        int tries=0 )
    {
      Path< TGraph > haplotype( hp_itr.get_graph_ptr() );
      get_uniq_full_haplotype( haplotype, hp_itr, hp_end, tries );
      if ( length( haplotype ) != 0 ){
        paths.push_back( std::move( haplotype ) );
      }
    }

  template< typename TGraph, typename TPathSet, typename TSpec >
      inline void
    get_uniq_patches( TPathSet& paths,
        GraphIter< TGraph, Haplotyper< TSpec > >& hp_itr,
        typename GraphIter< TGraph, Haplotyper< TSpec > >::end_type& hp_end,
        unsigned int k )
    { // :TODO:Fri Dec 01 21:26:\@cartoonist: the length of pre-context sequence won't be k all the time.
      hp_itr.set_raise_on_end( true );
      Path< TGraph > patch( hp_itr.get_graph_ptr() );
      Path< TGraph, Dynamic > frontier( hp_itr.get_graph_ptr() );
      typename Path< TGraph, Dynamic >::nodes_type::value_type marked;
      try {
        while ( true ) {
          marked = 0;
          if ( !frontier.empty() ) marked = frontier.get_nodes().back();
          // Bootstrap.
          if ( !marked ) util::extend_to_k( frontier, hp_itr, hp_end, k );
          else util::extend_to_k( frontier, hp_itr, hp_end,
              2*k + frontier.get_sequence_len() - frontier.get_seqlen_tail() );
          // Check the next patch distance to merge with previous patch if is less than k.
          if ( !patch.empty() && hp_itr[ frontier.get_nodes() ] ) {
            patch.set_right_by_len( k - 1 );
            paths.push_back( std::move( patch ) );
            clear( patch );
            rtrim_front_by_len( frontier, k, true );
          }
          else if ( !patch.empty() ) {
            // Nodes from first to the `marked` are already added.
            trim_front( frontier, marked );
            marked = 0;
            util::extend_to_k( frontier, hp_itr, hp_end, k );
          }
          if ( patch.empty() ) {
            // Search for a patch of length k that is not covered by visited paths of `hp_itr`.
            while ( hp_itr[ frontier.get_nodes() ] ) {
              add_node( frontier, *hp_itr );
              ltrim_front_by_len( frontier, k, true );
              ++hp_itr;
            }
          }
          // Extend the patch.
          patch += frontier;
          rtrim_front_by_len( frontier, k );
          while ( !hp_itr[ frontier.get_nodes() ] ) {
            add_node( frontier, *hp_itr );
            add_node( patch, *hp_itr );
            rtrim_front_by_len( frontier, k );
            ++hp_itr;
          }
        }
      }
      catch( const std::range_error& ) {
        if ( length( patch ) > 0 ) {
          if ( !hp_itr[ frontier.get_nodes() ] &&
              !rcontains( patch, frontier.get_nodes().rbegin(), frontier.get_nodes().rend() ) )
          {
            if ( marked != 0 ) trim_front( frontier, marked );
            patch += frontier;
          }
          paths.push_back( std::move( patch ) );
        }
        --hp_itr;  // save the traversed path and reset the Haplotyper iterator.
      }
      hp_itr.set_raise_on_end( false );
    }

  template< typename TPathSet, class TGraph, typename TSpec >
      inline bool
    get_uniq_patched_haplotype( TPathSet& paths,
        GraphIter< TGraph, Haplotyper< TSpec > >& hp_itr,
        typename GraphIter< TGraph, Haplotyper< TSpec > >::end_type& hp_end,
        unsigned int context_len )
    {
      assert( context_len != 0 );
      if ( hp_itr.level() == 0 ) {
        get_uniq_full_haplotype( paths, hp_itr, hp_end );
        return true;
      }
      unsigned int paths_no = paths.size();
      get_uniq_patches( paths, hp_itr, hp_end, context_len );
      if ( paths_no == paths.size() ) return false;
      return true;
    }

  template< typename TPathSet, class TGraph, typename TSpec >
      inline bool
    get_uniq_haplotype( TPathSet& paths,
        GraphIter< TGraph, Haplotyper< TSpec > >& hp_itr,
        typename GraphIter< TGraph, Haplotyper< TSpec > >::end_type& hp_end,
        unsigned int context_len,
        bool patched )
    {
      if ( patched )
        return get_uniq_patched_haplotype( paths, hp_itr, hp_end, context_len );
      get_uniq_full_haplotype( paths, hp_itr, hp_end );
      return true;
    }

  template< typename TPath, class TGraph >
      inline void
    get_rnd_full_haplotype( TPath& haplotype,
        GraphIter< TGraph, Haplotyper< Random > >& hp_itr,
        typename GraphIter< TGraph, Haplotyper< Random > >::end_type& hp_end )
    {
      while ( hp_itr != hp_end ) {
        haplotype.push_back( *hp_itr );
        ++hp_itr;
      }
      --hp_itr;
    }

  /* END OF PathIndex interface functions  ----------------------------------------- */

}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_PATHINDEX_HPP__ --- */

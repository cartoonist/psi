/**
 *    @file  pathset.hpp
 *   @brief  PathSet template class definition.
 *
 *  This header files contains PathSet template class definition and its interface
 *  functions.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Mar 27, 2018  17:50
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  PSI_PATHSET_HPP__
#define  PSI_PATHSET_HPP__

#include <vector>

#include <sdsl/bit_vectors.hpp>

#include "sequence.hpp"
#include "index.hpp"
#include "path.hpp"
#include "utils.hpp"


#define PATHSET_ID_SEPARATOR_CHAR ','


namespace psi {
  /**
   *  @brief  Represent a set of path with node ID query functionalities.
   *
   *  This class keeps a set of paths with some auxiliary data structure that allows
   *  more efficient node ID querie.
   */
  template< typename TPath, typename TSpec = DiskBased >
    class PathSet {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TPath value_type;
        typedef uint64_t size_type;     /* The max size type to be (de)serialized now. */
        typedef typename value_type::graph_type graph_type;
        typedef std::vector< value_type > container_type;
        typedef typename container_type::iterator iterator;
        typedef typename container_type::const_iterator const_iterator;
        typedef seqan::StringSet< YaString< TSpec > > stringset_type;
        /* XXX: `stringset_type::string_type` is underlying internal (in-memory) string
         *      ... type which might not be equivalent to `YaString<TSpec>`. */
        typedef typename stringset_type::string_type string_type;
        typedef seqan::Index< stringset_type, psi::FMIndex<> > index_type;
        typedef std::pair< size_type, typename value_type::size_type > pos_type;
        /* ====================  CONST MEMBERS  ====================================== */
        const char ID_SEPARATOR = PATHSET_ID_SEPARATOR_CHAR;
        /* ====================  LIFECYCLE      ====================================== */
        PathSet( ) = default;

        PathSet( PathSet const& ) = delete;
        PathSet& operator=( PathSet const& ) = delete;

        PathSet( PathSet&& other )
        {
          this->set = std::move( other.set );
          this->encids_set = std::move( other.encids_set );
          if ( other.encids_index.owns_text() ) {
            this->encids_index = std::move( other.encids_index );
          }
          else if ( other.encids_index.empty() ) {
            this->encids_index = index_type( this->encids_set );
          }
          else {
            this->encids_index = std::move( other.encids_index );
            this->encids_index.set_text_fibre( &this->encids_set, false );
          }
          this->bv_ids_set = std::move( other.bv_ids_set );

          for ( const auto& bv_id_breaks : this->bv_ids_set ) {
            this->rs_ids_set.push_back( sdsl::bit_vector::rank_1_type() );
            sdsl::util::init_support( this->rs_ids_set.back(), &bv_id_breaks );
          }
        }

        PathSet& operator=( PathSet&& other )
        {
          this->set = std::move( other.set );
          this->encids_set = std::move( other.encids_set );
          if ( other.encids_index.owns_text() ) {
            this->encids_index = std::move( other.encids_index );
          }
          else if ( other.encids_index.empty() ) {
            this->encids_index = index_type( this->encids_set );
          }
          else {
            this->encids_index = std::move( other.encids_index );
            this->encids_index.set_text_fibre( &this->encids_set, false );
          }
          this->bv_ids_set = std::move( other.bv_ids_set );

          this->rs_ids_set.clear();
          for ( const auto& bv_id_breaks : this->bv_ids_set ) {
            this->rs_ids_set.push_back( sdsl::bit_vector::rank_1_type() );
            sdsl::util::init_support( this->rs_ids_set.back(), &bv_id_breaks );
          }

          return *this;
        }

        ~PathSet( ) = default;
        /* ====================  OPERATORS     ======================================= */
          inline const value_type&
        operator[]( size_type idx ) const
        {
          return this->set[ idx ];
        }

          inline value_type&
        operator[]( size_type idx )
        {
          return this->set[ idx ];
        }
        /* ====================  METHODS       ======================================= */
          inline const_iterator
        begin( ) const
        {
          return this->set.begin();
        }

          inline const_iterator
        end( ) const
        {
          return this->set.end();
        }

          inline iterator
        begin( )
        {
          return this->set.begin();
        }

          inline iterator
        end( )
        {
          return this->set.end();
        }

          inline void
        push_back( value_type path )
        {
          if ( length( path ) == 0 ) {
            throw std::runtime_error( "attempting to add an empty path" );
          }

          this->set.push_back( std::move( path ) );
          psi::initialize( this->set.back() );

          std::string encids_str = this->get_encids_str( this->set.back() );
          appendValue( this->encids_set, encids_str );
          this->encids_index = index_type( this->encids_set );
          this->set_id_breaks( encids_str );
        }

        template< typename TPathSpec >
            inline void
          push_back( Path< graph_type, TPathSpec > path,
              enable_if_not_equal_t< typename value_type::spec_type, TPathSpec > tag = TPathSpec() )
          {
            value_type native_path( path.get_graph_ptr() );
            native_path = std::move( path );
            this->push_back( std::move( native_path ) );
          }

          inline std::vector< pos_type >
        get_occurrences( const string_type& idstr )
        {
          std::vector< pos_type > retval;
          seqan::Finder< index_type > finder( this->encids_index );
          while( find( finder, idstr ) ) {
            auto str_pos = beginPosition( finder );
            pos_type pos = { 0, 0 };
            pos.first = str_pos.first;
            pos.second = this->rank( str_pos );
            retval.push_back( std::move( pos ) );
          }

          return retval;
        }

        template< typename TPath2 >
            inline std::vector< pos_type >
          get_occurrences( const TPath2& path )
          {
            std::string encids_str = this->get_encids_str( path );
            return this->get_occurrences( encids_str );
          }

          inline bool
        found( const string_type& idstr )
        {
          seqan::Finder< index_type > finder( this->encids_index );
          find( finder, idstr );
          return finder.count() != 0;
        }

        template< typename TPath2 >
            inline bool
          found( const TPath2& path )
          {
            std::string encids_str = this->get_encids_str( path );
            return this->found( encids_str );
          }

          inline typename string_type::size_type
        rank( typename stringset_type::pos_type pos ) const
        {
          return this->rs_ids_set[ pos.first ]( pos.second );
        }

          inline void
        reserve( size_type size )
        {
          this->set.reserve( size );
          this->encids_set.reserve( size );
          this->bv_ids_set.reserve( size );
          this->rs_ids_set.reserve( size );
        }

          inline size_type
        size( ) const
        {
          return this->set.size();
        }

          inline void
        clear( )
        {
          this->set.clear();
          this->encids_set.clear();
          this->encids_index.clear();
          sdsl::util::clear( this->bv_ids_set );
          sdsl::util::clear( this->rs_ids_set );
        }

          inline void
        initialize( )
        {
          indexRequire( this->encids_index, seqan::FibreSALF() );
        }

          inline void
        serialize( std::ostream& out )
        {
          psi::serialize( out, static_cast< size_type >( this->size() ) );
          for ( auto&& path : this->set ) {
            psi::save( path, out );
          }

          indexRequire( this->encids_index, seqan::FibreSALF() );
          save( this->encids_index, out );

          for ( auto&& bv_id_breaks : this->bv_ids_set ) {
            bv_id_breaks.serialize( out );
          }
        }

          inline void
        load( std::istream& in, const graph_type* graph_ptr )
        {
          this->clear();
          size_type paths_num = 0;
          deserialize( in, paths_num );
          this->set.reserve( paths_num );
          this->encids_set.reserve( paths_num );
          this->bv_ids_set.reserve( paths_num );
          this->rs_ids_set.reserve( paths_num );
          for ( size_type i = 0; i < paths_num; ++i ) {
            value_type path( graph_ptr );
            psi::open( path, in );
            this->set.push_back( std::move( path ) );
          }

          open( this->encids_index, in );

          for ( size_type i = 0; i < paths_num; ++i ) {
            sdsl::bit_vector bv;
            bv.load( in );
            this->bv_ids_set.push_back( std::move( bv ) );
            this->rs_ids_set.push_back( sdsl::bit_vector::rank_1_type() );
            sdsl::util::init_support( this->rs_ids_set.back(), &this->bv_ids_set.back() );
          }
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        container_type set;
        stringset_type encids_set;
        index_type encids_index;
        std::vector< sdsl::bit_vector > bv_ids_set;
        std::vector< sdsl::bit_vector::rank_1_type > rs_ids_set;
        /* ====================  METHODS       ======================================= */
        template< typename TPath2 >
            inline string_type
          get_encids_str( const TPath2& path ) const
          {
            string_type encids = std::string( 1, ID_SEPARATOR );
            for ( const auto& node_id : path.get_nodes() ) {
              encids += std::to_string( node_id ) + std::string( 1, ID_SEPARATOR );
            }

            return encids;
          }

          inline void
        set_id_breaks( const string_type& encids )
        {
          sdsl::bit_vector bv( encids.size(), 0 );
#ifdef PSI_DEBUG_ENABLED
          bool flag = false;
          assert( encids[ 0 ] == ID_SEPARATOR );
          assert( encids[ 1 ] != ID_SEPARATOR );
          assert( encids[ encids.size() - 1 ] == ID_SEPARATOR );
#endif
          for ( sdsl::bit_vector::size_type i = 2; i < bv.size(); ++i ) {
            if ( encids[ i ] == ID_SEPARATOR ) {
              bv[ i - 1 ] = 1;
#ifdef PSI_DEBUG_ENABLED
              /* Sanity check: there should not be two 1s subsequently. */
              assert( !flag );
              flag = true;
            }
            else {
              flag = false;
#endif
            }
          }
          this->bv_ids_set.push_back( std::move( bv ) );
          this->rs_ids_set.push_back( sdsl::bit_vector::rank_1_type() );
          sdsl::util::init_support( this->rs_ids_set.back(), &this->bv_ids_set.back() );
        }
    };

  /* PathSet interface functions  -------------------------------------------------- */

  template< typename TPath, typename TSpec >
      inline bool
    save( PathSet< TPath, TSpec >& set, const std::string& file_path )
    {
      std::ofstream ofs( file_path, std::ofstream::out | std::ofstream::binary );
      if( !ofs ) return false;
      save( set, ofs );
      return true;
    }

  template< typename TPath, typename TSpec >
      inline bool
    open( PathSet< TPath, TSpec >& set, typename TPath::graph_type const* graph_ptr,
        const std::string& file_path )
    {
      std::ifstream ifs( file_path, std::ifstream::in | std::ifstream::binary );
      if( !ifs ) return false;
      set.load( ifs, graph_ptr );
      return true;
    }

  /**
   *  @brief  Check whether a path is covered by a PathSet.
   *
   *  @param  path The given path as Path class instance.
   *  @param  pset The path index as PathSet class instance.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  Using created index on node IDs string set, it finds the node IDs string of the
   *  query path in the paths set. If the number of occurrences is not zero, the query
   *  path is covered by the paths set.
   */
  template< typename TGraph, typename TPathSpec1, typename TPathSpec2, typename TSpec >
      inline bool
    covered_by( const Path< TGraph, TPathSpec1 >& path,
        PathSet< Path< TGraph, TPathSpec2 >, TSpec >& pset )
    {
      return pset.found( path );
    }

  /* END OF PathSet interface functions  ------------------------------------------- */
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_PATHSET_HPP__ --- */

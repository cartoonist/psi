/*
 * =====================================================================================
 *
 * Filename: vargraph.h
 *
 * Created: Fri Nov 11, 2016  01:08
 *
 * Description: VarGraph class definition.
 *
 * Copyright (c) 2016, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari, <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#ifndef VARGRAPH_H__
#define VARGRAPH_H__

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <deque>
#include <utility>
#include <cstdint>
#include <random>

#ifdef HAVE_MEMCPY
#  define MACRO_STACK
#  pragma push_macro("HAVE_MEMCPY")
#  undef HAVE_MEMCPY
#endif

#include <xg/xg.hpp>

#ifdef MACRO_STACK
#  pragma pop_macro("HAVE_MEMCPY")
#  undef MACRO_STACK
#endif

#include <seqan/basic.h>

#include "graph_iter.h"
#include "utils.h"


namespace grem
{
  class VarGraph : public xg::XG
  {
    public:
      // typedefs
      typedef vg::Node node_type;                               /**< @brief Node type. */
      typedef decltype( vg::Node().id() ) nodeid_type;          /**< @brief Node ID type. */
      typedef std::size_t rank_type;                            /**< @brief Node ID type. */
      typedef decltype( vg::Position().offset() ) offset_type;  /**< @brief Node offset type. */

      // Constructors
      VarGraph( void ) : XG( ) { }

      VarGraph( std::ifstream &ifs, bool fmt_xg = true )
      {
        if ( fmt_xg ) {
          load( ifs );
        }
        else {
          from_stream( ifs );
        }
      }

      VarGraph( vg::Graph &vg_graph ) : XG( vg_graph )
      { }

      // Public methods
        inline bool
      is_branch ( nodeid_type node_id ) const
      {
        if ( this->edges_from( node_id ).size() > 1 ) {
          return true;
        }
        return false;
      }  /* -----  end of method is_branch  ----- */

        inline bool
      is_merge ( nodeid_type node_id ) const
      {
        if ( this->edges_to( node_id ).size() > 1 ) {
          return true;
        }
        return false;
      }  /* -----  end of method is_merge  ----- */

        inline bool
      has_edges_from( nodeid_type node_id ) const
      {
        return this->edges_from( node_id ).size() != 0;
      }

        inline bool
      has_edges_to( nodeid_type node_id ) const
      {
        return this->edges_to( node_id ).size() != 0;
      }
  };

  /* Path interface functions forwards  ---------------------------------------- */

  template< typename TSpec >
    class Path;

  template< typename TSpec >
      void
    save( const Path< TSpec >& path, const std::string& file_name );
  template< typename TSpec >
      void
    _add_node( Path< TSpec >& path, const VarGraph::nodeid_type& node_id );
  template< typename TSpec >
      void
    add_node( Path< TSpec >& path, const VarGraph::nodeid_type& node_id );
  template< typename TSpec >
      std::string
    sequence( const Path< TSpec >& path );
  template< typename TSpec >
      void
    _clear( Path< TSpec >& path );
  template< typename TSpec >
      void
    clear( Path< TSpec >& path );
  template< typename TSpec >
      void
    reserve( Path< TSpec >& path, typename Path< TSpec >::size_type size );
  template< typename TSpec >
      typename Path< TSpec >::size_type
    length( const Path< TSpec >& path );
  template< typename TSpec >
      bool
    contains( const Path< TSpec >& path, VarGraph::nodeid_type node_id );
  template< typename TSpec >
      void
    trim( Path< TSpec >& path, VarGraph::nodeid_type node_id=0 );

  /* END OF path interface functions forwards  --------------------------------- */

  /* Path specialization tags. */
  struct CompactStrategy;
  struct FullStrategy;
  typedef seqan::Tag< CompactStrategy > Compact;
  typedef seqan::Tag< FullStrategy > Full;

  /**
   *  @brief  Path template class.
   *
   *  Represent a path in the variation graph with efficient node ID query.
   */
  template< typename TSpec = void >
    class Path
    {
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        const VarGraph* vargraph;
        std::vector< VarGraph::nodeid_type > nodes;
        std::unordered_set< VarGraph::nodeid_type > nodes_set;
        /* ====================  METHODS       ======================================= */
        /**
         * @brief  pop the last node from the path.
         */
          inline void
        pop_back( )
        {
          this->nodes_set.erase( this->nodes_set.find( this->nodes.back() ) );
          this->nodes.pop_back();
        }  /* -----  end of method pop_back  ----- */
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename decltype( nodes )::size_type size_type;
        /* ====================  LIFECYCLE     ======================================= */
        Path( const VarGraph* g )
          : vargraph( g )
        { }

        Path( const VarGraph* g, std::vector< VarGraph::nodeid_type >&& p )
          : Path( g )
        {
          this->set_nodes( std::move( p ) );
        }

        Path( const VarGraph* g, const std::vector< VarGraph::nodeid_type >& p )
          : Path( g, std::vector< VarGraph::nodeid_type >( p ) )
        { }
        /* ====================  ACCESSORS     ======================================= */

        /**
         *  @brief  getter function for vargraph.
         */
          inline const VarGraph*
        get_vargraph( ) const
        {
          return this->vargraph;
        }  /* -----  end of method get_vargraph  ----- */

        /**
         *  @brief  getter function for nodes.
         */
          inline const std::vector< VarGraph::nodeid_type >&
        get_nodes( ) const
        {
          return this->nodes;
        }  /* -----  end of method get_nodes  ----- */
        /* ====================  MUTATORS      ======================================= */
        /**
         *  @brief  setter function for vargraph.
         */
          inline void
        set_vargraph( const VarGraph* value )
        {
          this->vargraph = value;
        }  /* -----  end of method set_vargraph  ----- */

        /**
         *  @brief  setter function for nodes.
         */
          inline void
        set_nodes( std::vector< VarGraph::nodeid_type >&& value )
        {
          this->nodes = std::move( value );
          this->nodes_set.clear();
          this->nodes_set.reserve( this->nodes.size() );
          std::copy( this->nodes.begin(), this->nodes.end(),
              std::inserter( this->nodes_set, this->nodes_set.end() ) );
        }  /* -----  end of method set_nodes  ----- */

        /**
         *  @brief  setter function for nodes.
         */
          inline void
        set_nodes( const std::vector< VarGraph::nodeid_type >& value )
        {
          this->set_nodes( std::vector< VarGraph::nodeid_type >( value ) );
        }  /* -----  end of method set_nodes  ----- */
        /* ====================  INTERFACE FUNCTIONS  ================================ */
          friend void
        save< TSpec >( const Path< TSpec >& path, const std::string& file_name );
          friend void
        _add_node< TSpec >( Path< TSpec >& path, const VarGraph::nodeid_type& node_id );
          friend std::string
        sequence< TSpec >( const Path< TSpec >& path );
          friend void
        _clear< TSpec >( Path< TSpec >& path );
          friend void
        reserve< TSpec >( Path< TSpec >& path, size_type size );
          friend bool
        contains< TSpec >( const Path< TSpec >& path, VarGraph::nodeid_type node_id );
          friend void
        trim< TSpec >( Path< TSpec >& path, VarGraph::nodeid_type node_id );
    };  /* -----  end of template class Path  ----- */

  /**
   *  @brief  Path template class specialized for Full path.
   *
   *  Represent a path in the variation graph with efficient node ID query and stored
   *  path sequence.
   */
  template< >
    class Path< Full >
    {
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        const VarGraph* vargraph;
        std::vector< VarGraph::nodeid_type > nodes;
        std::unordered_set< VarGraph::nodeid_type > nodes_set;
        std::string path_sequence;
        /* ====================  METHODS       ======================================= */
        /**
         * @brief  pop the last node from the path.
         */
          inline void
        pop_back( )
        {
          auto path_new_size = this->path_sequence.size()
            - this->vargraph->node_length( this->nodes.back() );
          this->path_sequence.resize( path_new_size );
          this->nodes_set.erase( this->nodes_set.find( this->nodes.back() ) );
          this->nodes.pop_back();
        }  /* -----  end of method pop_back  ----- */
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename decltype( nodes )::size_type size_type;
        /* ====================  LIFECYCLE     ======================================= */
        Path( const VarGraph* g )
          : vargraph( g )
        { }

        Path( const VarGraph* g, std::vector< VarGraph::nodeid_type >&& p )
          : Path( g )
        {
          this->set_nodes( std::move( p ) );
        }

        Path( const VarGraph* g, const std::vector< VarGraph::nodeid_type >& p )
          : Path( g, std::vector< VarGraph::nodeid_type >( p ) )
        { }
        /* ====================  ACCESSORS     ======================================= */

        /**
         *  @brief  getter function for vargraph.
         */
          inline const VarGraph*
        get_vargraph( ) const
        {
          return this->vargraph;
        }  /* -----  end of method get_vargraph  ----- */

        /**
         *  @brief  getter function for nodes.
         */
          inline const std::vector< VarGraph::nodeid_type >&
        get_nodes( ) const
        {
          return this->nodes;
        }  /* -----  end of method get_nodes  ----- */

        /**
         *  @brief  getter function for path_sequence.
         */
          inline const std::string&
        get_sequence( ) const
        {
          return this->path_sequence;
        }  /* -----  end of method get_sequence  ----- */
        /* ====================  MUTATORS      ======================================= */
        /**
         *  @brief  setter function for vargraph.
         */
          inline void
        set_vargraph( const VarGraph* value )
        {
          this->vargraph = value;
        }  /* -----  end of method set_vargraph  ----- */

        /**
         *  @brief  setter function for nodes.
         */
          inline void
        set_nodes( std::vector< VarGraph::nodeid_type >&& value )
        {
          this->nodes = std::move( value );
          this->nodes_set.clear();
          this->nodes_set.reserve( this->nodes.size() );
          std::copy( this->nodes.begin(), this->nodes.end(),
              std::inserter( this->nodes_set, this->nodes_set.end() ) );
          this->path_sequence = sequence( *this );
        }  /* -----  end of method set_nodes  ----- */

        /**
         *  @brief  setter function for nodes.
         */
          inline void
        set_nodes( const std::vector< VarGraph::nodeid_type >& value )
        {
          this->set_nodes( std::vector< VarGraph::nodeid_type >( value ) );
        }  /* -----  end of method set_nodes  ----- */
        /* ====================  INTERFACE FUNCTIONS  ================================ */
          friend void
        save< Full >( const Path< Full >& path, const std::string& file_name );
          friend void
        _add_node< Full >( Path< Full >& path, const VarGraph::nodeid_type& node_id );
          friend void
        add_node< Full >( Path< Full >& path, const VarGraph::nodeid_type& node_id );
          friend std::string
        sequence< Full >( const Path< Full >& path );
          friend void
        _clear< Full >( Path< Full >& path );
          friend void
        clear< Full >( Path< Full >& path );
          friend void
        reserve< Full >( Path< Full >& path, size_type size );
          friend bool
        contains< Full >( const Path< Full >& path, VarGraph::nodeid_type node_id );
          friend void
        trim< Full >( Path< Full >& path, VarGraph::nodeid_type node_id );
    };  /* -----  end of specialized template class Path  ----- */

  /* Normal and Full Path interface functions  --------------------------------- */

  /**
   *  @brief  Save the path to file.
   *
   *  @param  path The path.
   *  @param  file_name The name of the file to be written.
   *
   *  It saves the sequence of the node IDs into the file.
   */
  template< typename TSpec >
      inline void
    save( const Path< TSpec >& path, const std::string& file_name )
    {
      std::ofstream ofs( file_name, std::ofstream::out | std::ofstream::binary );
      if( !ofs ) {
        throw std::runtime_error( "cannot open file '" + file_name + "'" );
      }
      serialize( ofs, path.nodes, path.nodes.begin(), path.nodes.end() );
    }  /* -----  end of function save  ----- */

  /**
   *  @brief  Extend the path forward.
   *
   *  @param  path The path.
   *  @param  node_id The new node ID.
   *
   *  Internal function to add the `node_id` to the end of the current path.
   */
  template< typename TSpec >
      inline void
    _add_node( Path< TSpec >& path, const VarGraph::nodeid_type& node_id )
    {
      path.nodes.push_back( node_id );
      path.nodes_set.insert( node_id );
    }  /* -----  end of template function add_node  ----- */

  /**
   *  @brief  Extend the path forward.
   *
   *  @param  path The path.
   *  @param  node_id The new node ID.
   *
   *  Add the `node_id` to the end of the current path. Specialized for Normal Path.
   */
  template< >
      inline void
    add_node( Path<>& path, const VarGraph::nodeid_type& node_id )
    {
      _add_node( path, node_id );
    }

  /**
   *  @brief  Extend the path forward.
   *
   *  @param  path The path.
   *  @param  node_id The new node ID.
   *
   *  Add the `node_id` to the end of the current path. Specialized for Full Path.
   */
  template< >
      inline void
    add_node( Path< Full >& path, const VarGraph::nodeid_type& node_id )
    {
      _add_node( path, node_id );

      assert( path.vargraph != nullptr );

      path.path_sequence += path.vargraph->node_sequence( node_id );
    }

  /**
   *  @brief  Compute sequence from the nodes vector.
   *
   *  @param  path The path.
   *  @return sequence represented by the path.
   *
   *  Compute the sequence of the path from nodes vector.
   */
  template< typename TSpec >
      inline std::string
    sequence( const Path< TSpec >& path )
    {
      assert( path.vargraph != nullptr );

      std::string repr_str;
      for ( const auto& node_id : path.nodes ) {
        repr_str += path.vargraph->node_sequence( node_id );
      }
      return repr_str;
    }  /* -----  end of template function sequence  ----- */

  /**
   *  @brief  Clear the path.
   *
   *  @param  path The path.
   *
   *  Internal function to clear nodes vector and nodes set.
   */
  template< typename TSpec >
      inline void
    _clear( Path< TSpec >& path )
    {
      path.nodes.clear();
      path.nodes_set.clear();
    }  /* -----  end of template function _clear  ----- */

  /**
   *  @brief  Clear the path.
   *
   *  @param  path The path.
   *
   *  Specialized for Normal Path.
   */
  template< >
      inline void
    clear( Path< >& path )
    {
      _clear( path );
    }  /* -----  end of template function clear  ----- */

  /**
   *  @brief  Clear the path.
   *
   *  @param  path The path.
   *
   *  Specialized for Full Path.
   */
  template< >
      inline void
    clear( Path< Full >& path )
    {
      _clear( path );
      path.path_sequence.clear();
    }  /* -----  end of template function clear  ----- */

  /**
   *  @brief  Reserve memory for the path.
   *
   *  @param  path The path.
   *  @param  size The target capacity.
   *
   *  Internal function to reserve memory for nodes and nodes set.
   */
  template< typename TSpec >
      inline void
    reserve( Path< TSpec >& path, typename Path< TSpec >::size_type size )
    {
      path.nodes.reserve( size );
      path.nodes_set.reserve( size );
    }  /* -----  end of template function reserve  ----- */

  /**
   *  @brief  Get path length (number of nodes).
   *
   *  @param  path The path.
   *  @return the number of nodes on the path.
   *
   *  It gets the length of the path which is the number of nodes in the path.
   */
  template< typename TSpec >
      inline typename Path< TSpec >::size_type
    length( const Path< TSpec >& path )
    {
      return path.get_nodes().size();
    }  /* -----  end of template function length  ----- */

  /**
   *  @brief  Trim the path from the node whose ID matches the given ID.
   *
   *  @param  path The path to be trimmed.
   *  @param  node_id The node ID from where the path should be trimmed.
   *
   *  It pops the last nodes in the nodes list until it reaches to the given node ID
   *  inclusive. It will yield an empty path if `node_id` is not among the path nodes.
   *  If `node_id` is equal to zero, an invalid node ID, it pops only the last node from
   *  the path.
   */
  template< typename TSpec >
      inline void
    trim( Path< TSpec >& path, VarGraph::nodeid_type node_id )
    {
      VarGraph::nodeid_type last_node;
      bool found = false;
      while ( !found && length( path ) != 0 ) {
        last_node = path.get_nodes().back();
        if ( node_id == 0 || last_node == node_id ) found = true;
        path.pop_back();
      }
    }  /* -----  end of template function trim  ----- */

  /* END OF Normal and Full Path interface functions  -------------------------- */

  /**
   *  @brief  Path template class Compact specialization.
   *
   *  Represent a path in the variation graph with efficient node ID query.
   */
  template< >
    class Path< Compact >
    {
      public:
        /* ====================  DATA MEMBERS  ======================================= */
        std::unordered_set< VarGraph::nodeid_type > nodes_set;
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename decltype( nodes_set )::size_type size_type;
        /* ====================  LIFECYCLE     ======================================= */
        Path( ) = default;
        Path( const std::vector< VarGraph::nodeid_type >& p )
        {
          this->set_nodes( p );
        }
        /* ====================  MUTATORS      ======================================= */
        /**
         *  @brief  Set the nodes in the path (clear the previous state).
         */
          inline void
        set_nodes( const std::vector< VarGraph::nodeid_type >& value )
        {
          this->nodes_set.clear();
          this->nodes_set.reserve( value.size() );
          std::copy( value.begin(), value.end(),
              std::inserter( this->nodes_set, this->nodes_set.end() ) );
        }  /* -----  end of method set_nodes  ----- */
        /* ====================  INTERFACE FUNCTIONS  ================================ */
          friend void
        save< Compact >( const Path< Compact >& path, const std::string& file_name );
          friend void
        add_node< Compact >( Path< Compact >& path,
            const VarGraph::nodeid_type& node_id );
          friend void
        clear< Compact >( Path< Compact >& path );
          friend void
        reserve< Compact >( Path< Compact >& path, size_type size );
          friend size_type
        length< Compact >( const Path< Compact >& path );
          friend bool
        contains< Compact >( const Path< Compact >& path,
            VarGraph::nodeid_type node_id );
    };  /* -----  end of specialized template class Path  ----- */

  /* Compact Path interface functions definitions  ----------------------------- */

  /**
   *  @brief  Save the path to file -- Compact specialization.
   *
   *  @param  path The path.
   *  @param  file_name The name of the file to be written.
   *
   *  It saves the sequence of the node IDs into the file.
   */
  template< >
      inline void
    save( const Path< Compact >& path, const std::string& file_name )
    {
      std::ofstream ofs( file_name, std::ofstream::out | std::ofstream::binary );
      if( !ofs ) {
        throw std::runtime_error( "cannot open file '" + file_name + "'" );
      }
      serialize( ofs, path.nodes_set, path.nodes_set.begin(), path.nodes_set.end() );
    }  /* -----  end of function save  ----- */

  /**
   *  @brief  Extend the path.
   *
   *  @param  path The path.
   *  @param  node_id The new node ID to be added into the path node set.
   */
  template< >
      inline void
    add_node( Path< Compact >& path, const VarGraph::nodeid_type& node_id )
    {
      path.nodes_set.insert( node_id );
    }

  /**
   *  @brief  Clear the path.
   *
   *  @param  path The path.
   */
  template< >
      inline void
    clear( Path< Compact >& path )
    {
      path.nodes_set.clear();
    }

  /**
   *  @brief  Reserve memory for the path.
   *
   *  @param  path The path.
   *  @param  size The target capacity.
   *
   *  It reserves memory for nodes set.
   */
  template< >
      inline void
    reserve( Path< Compact >& path, Path< Compact >::size_type size )
    {
      path.nodes_set.reserve( size );
    }

  /**
   *  @brief  Get path length.
   *
   *  @param  path The path.
   *  @return The path length.
   *
   *  It gets the length of the path.
   */
  template< >
      inline Path< Compact >::size_type
    length( const Path< Compact >& path )
    {
      return path.nodes_set.size();
    }  /* -----  end of function length  ----- */

  /* END OF Compact Path interface functions  ---------------------------------- */

  /* Path interface functions  ------------------------------------------------- */

  /**
   *  @brief  Load the path from file.
   *
   *  @param  path The path.
   *  @param  file_name The name of the file to be read.
   *
   *  It loads the sequence of the node IDs from the file.
   */
  template< typename TSpec >
      void
    load( Path< TSpec >& path, const std::string& file_name )
    {
      std::ifstream ifs( file_name, std::ifstream::in | std::ifstream::binary );
      if( !ifs ) {
        throw std::runtime_error( "cannot open file '" + file_name + "'" );
      }
      std::vector< VarGraph::nodeid_type > nodes;
      deserialize( ifs, nodes, std::back_inserter( nodes ) );
      path.set_nodes( std::move( nodes ) );
    }  /* -----  end of template function load  ----- */

  /**
   *  @brief  Check whether the given node ID is on the path or not.
   *
   *  @param  path The path.
   *  @param  node_id The ID of the node.
   *  @return `true` if the node ID is on the path.
   *
   *  Check if the given node ID presents in the path.
   */
  template< typename TSpec >
      inline bool
    contains( const Path< TSpec >& path, VarGraph::nodeid_type node_id )
    {
      return path.nodes_set.find( node_id ) != path.nodes_set.end();
    }  /* -----  end of template function contains  ----- */

  /**
   *  @brief  Check whether this path contains another path.
   *
   *  @param  path The path.
   *  @param  begin The begin iterator of node IDs set of the path to be checked.
   *  @param  end The end iterator of node IDs set of the path to be checked.
   *  @return `true` if this path is a superset of the given path; otherwise `false`
   *          -- including the case that the given path is empty; i.e.
   *          `end == begin`.
   *
   *  The input path should be smaller than the path. It checks if the nodes of the
   *  given path is present in the node set. It does not check the order of the
   *  nodes.
   */
  template< typename TSpec, typename TIter >
      inline bool
    contains( const Path< TSpec >& path, TIter begin, TIter end )
    {
      if ( end - begin > 0 ) {
        auto&& on_path = [&path]( const VarGraph::nodeid_type& i ) {
          return contains( path, i );
        };

        if ( std::all_of( begin, end, on_path ) ) return true;
      }

      return false;
    }  /* -----  end of template function contains  ----- */

  /**
   *  @brief  Check whether this path contains another path.
   *
   *  @param  path The path.
   *  @param  path_nodes The node IDs of the given path to be checked.
   *  @return `true` if this path is a superset of the given path.
   *
   *  Overloaded. See `contains( const Path< TSpec >&, TIter, TIter )`.
   */
  template< typename TSpec, typename TContainer >
      inline bool
    contains( const Path< TSpec >& path, const TContainer& path_nodes )
    {
      return contains( path, path_nodes.begin(), path_nodes.end() );
    }  /* -----  end of template function contains  ----- */

  /* END OF Path interface functions  ------------------------------------------ */

  /* Graph interface functions  ------------------------------------------------ */

  /**
   *  @brief  Check whether a path is covered by a set of paths.
   *
   *  @param  begin The begin iterator of the path as set of node IDs.
   *  @param  end The end iterator of the path as set of node IDs.
   *  @param  path_set_begin The begin iterator of the paths set.
   *  @param  path_set_end The end iterator of the paths set.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  This function checks if all nodes of the given path is covered by at least one
   *  path in the given set of paths. The path set should be a container of the Path
   *  class.
   */
  template< class Iter1, class Iter2 >
      inline bool
    covered_by( Iter1 begin, Iter1 end, Iter2 paths_set_begin, Iter2 paths_set_end )
    {
      for ( auto itr = paths_set_begin; itr != paths_set_end; ++itr ) {
        if ( contains( (*itr), begin, end ) )
        {
          return true;
        }
      }
      return false;
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Check whether a path is covered by a set of paths.
   *
   *  @param  path_nodes The given path to be checked as a container of node IDs
   *  @param  paths_set A set of paths as a container of `Path`.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  Overloaded. See `covered_by( TIter1, TIter1, TIter2, TIter2 )`.
   */
  template< class TContainer1, class TContainer2 >
      inline bool
    covered_by( const TContainer1& path_nodes, const TContainer2& paths_set )
    {
      return covered_by( path_nodes.begin(), path_nodes.end(),
          paths_set.begin(), paths_set.end() );
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Check whether a node is covered by a set of paths.
   *
   *  @param  node_id The ID of the given node to be checked.
   *  @param  paths_set A set of paths as a container of `Path`.
   *  @return `true` if the node is on one of the path in the paths set, otherwise
   *          `false`.
   *
   *  This function simply checks if a node is on any of the path in the given set of
   *  paths.
   */
  template< typename TContainer >
      inline bool
    covered_by( VarGraph::nodeid_type node_id, const TContainer& paths_set )
    {
      for ( const auto& path : paths_set ) {
        if ( contains( path, node_id ) ) return true;
      }

      return false;
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Return the path coverage of a given node by its ID.
   *
   *  @param  node_id The node ID of the node to be checked.
   *  @param  paths_set A set of paths as a container of `Path`.
   *  @return the number of paths that cover the given node.
   *
   *  This function get a node ID and return the number of paths that cover the node.
   */
  template< typename TContainer >
      inline std::size_t
    get_path_coverage( const VarGraph::nodeid_type& node_id, const TContainer& paths_set )
    {
      std::size_t coverage = 0;
      for ( const auto & path : paths_set ) {
        if ( contains( path, node_id ) ) {
          ++coverage;
        }
      }
      return coverage;
    }  /* -----  end of template function get_path_coverage  ----- */

  /**
   *  @brief  Get the node ID of an ajacent node randomly.
   *
   *  @param  vargraph The variation graph.
   *  @param  node_id The ID of the node whose an adjacent node should be returned.
   *  @return an adjacent node ID if available any; otherwise 0 -- an invalid node ID.
   *
   *  Picking one of the adjacent nodes by generating a pseudo-random number with
   *  uniform distribution in [0, out-degree(v)].
   */
    inline VarGraph::nodeid_type
  get_random_adjacent ( const VarGraph &vargraph, VarGraph::nodeid_type node_id )
  {
    if ( !vargraph.has_edges_from( node_id ) ) {
      return 0;
    }

    auto fwd_edges = vargraph.edges_from(node_id);

    std::random_device rd;  // Will be used to obtain a seed for the random no. engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, fwd_edges.size() - 1);

    return fwd_edges.at(dis(gen)).to();
  }  /* -----  end of function get_random_adjacent  ----- */

  /**
   *  @brief  Get the ID of an adjacent node with least coverage.
   *
   *  @param  vargraph The variation graph.
   *  @param  node_id The ID ot the node whose an adjacent node should be returned.
   *  @param  paths_set A set of paths as a container of `Path`.
   *  @return the ID of the one of adjacent nodes with least coverage. If there are
   *          multiple nodes with minimum coverage, it returns one of them. If all nodes
   *          are covered equally or no forward edge exists, it returns zero.
   *
   *  Calculate coverage for all adjacent nodes and find the smallest one.
   */
  template< typename TContainer >
      inline VarGraph::nodeid_type
    least_covered_adjacent( const VarGraph& vargraph, VarGraph::nodeid_type node_id,
        const TContainer& paths_set )
    {
      if ( vargraph.has_edges_from( node_id ) ) {
        auto fwd_edges = vargraph.edges_from( node_id );

        VarGraph::nodeid_type first_adj_id = ( *fwd_edges.begin() ).to();

        VarGraph::nodeid_type lc_node_id = first_adj_id;
        unsigned int lc = get_path_coverage( first_adj_id, paths_set );
        bool equally_covered = true;

        for ( auto e_itr = ++fwd_edges.begin(); e_itr != fwd_edges.end(); ++e_itr ) {
          const VarGraph::nodeid_type& next_node = ( *e_itr ).to();
          unsigned int next_node_cov = get_path_coverage( next_node, paths_set );

          if ( equally_covered && lc != next_node_cov ) {
            equally_covered = false;
          }

          if ( next_node_cov < lc ) {
            lc = next_node_cov;
            lc_node_id = next_node;
          }
        }

        if ( !equally_covered ) {
          return lc_node_id;
        }
      }

      return 0;
    }  /* -----  end of template function least_covered_adjacent  ----- */

  /* END OF graph interface functions  ----------------------------------------- */

  /* GRAPH ITERATORS  ============================================================ */

  /* Traits template specialization  ------------------------------------------- */

  struct BFSStrategy;
  struct DFSStrategy;
  struct BacktrackStrategy;
  struct HaplotypeStrategy;

  typedef seqan::Tag< BFSStrategy > BFS;
  typedef seqan::Tag< DFSStrategy > DFS;
  typedef seqan::Tag< BacktrackStrategy > Backtracker;
  typedef seqan::Tag< HaplotypeStrategy > Haplotyper;

  /**
   *  @brief  Breadth-first search graph iterator tag.
   *
   *  Specialization of generic graph iterator tag BFS for VarGraph.
   */
  template< >
    struct GraphIterTraits < VarGraph, BFS >
    {
      typedef VarGraph::nodeid_type Value;
      typedef Value Level;
      typedef std::deque< std::pair< Value, Level > > TContainer;

      struct pair_hash
      {
        inline std::size_t operator()(const std::pair< Value, Level > & v) const
        {
          return std::hash< Value >()(v.first);
        }
      };

      struct pair_pred
      {
        inline bool operator()
          (const std::pair< Value, Level > & v,
           const std::pair< Value, Level > & u)
          const noexcept
        {
          if (v.first == u.first) return true;
          else return false;
        }
      };

      typedef std::unordered_set< TContainer::value_type, pair_hash, pair_pred > TSet;
      typedef void* TState;
    };

  /**
   *  @brief  Backtracker graph iterator tag.
   *
   *  Specialization of generic graph iterator tag BacktrackerIter for VarGraph.
   */
  template < >
    struct GraphIterTraits < VarGraph, Backtracker > {
      typedef VarGraph::nodeid_type Value;
      typedef Value Level;
      typedef std::deque< std::pair< Value, Value > > TContainer;
      typedef void* TSet;
      typedef struct {
        Value start;                            /**< @brief Start node ID. */
        Value buffer;                           /**< @brief Buffer node ID. 0=nothing */
        bool end;                               /**< @brief End flag. */
      } TState;
    };  /* ----------  end of struct Backtracker  ---------- */

  /**
   *  @brief  Haplotyper graph iterator tag.
   *
   *  Specialization of generic graph iterator tag HaplotyperIter for VarGraph.
   */
  template < >
    struct GraphIterTraits < VarGraph, Haplotyper > {
      typedef VarGraph::nodeid_type Value;
      typedef Value Level;
      typedef std::deque< Value > TContainer;
      /**< @brief Set of visited paths. */
      typedef std::vector< Path< Compact > > TSet;
      typedef struct {
        Value start;                            /**< @brief Start node ID. */
        bool end;                               /**< @brief End flag. */
        Path< Compact > current_path;
        unsigned int setback;
      } TState;
    };  /* ----------  end of struct HaplotyperIter  ---------- */

  /* END OF tags template specialization  -------------------------------------- */

}  /* -----  end of namespace grem  ----- */

namespace seqan {
  /**
   *  @brief  Iterator template class specialization for VarGraph.
   *
   *  This class extends existing Iterator class in seqan namespace.
   */
  template < typename TSpec >
    class Iterator < grem::VarGraph, TSpec >
    {
      public:
        typedef grem::GraphIter < grem::VarGraph, TSpec > Type;
        /* ====================  TYPEDEFS      ======================================= */
    };  /* ----------  end of template class Iterator  ---------- */
}  /* -----  end of namespace seqan  ----- */

namespace grem {
  /* Haplotyper iterator interface function declaration  ------------------------- */

    void
  get_uniq_haplotype( Path<>& haplotype,
      typename seqan::Iterator< VarGraph, Haplotyper >::Type& iter,
      int tries=0 );

  /* END OF Haplotyper iterator interface function declaration  ------------------ */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef VARGRAPH_H__  ----- */

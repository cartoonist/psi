/*
 * =====================================================================================
 *
 * Filename: vargraph_iter.h
 *
 * Created: Wed Jan 18, 2017  14:40
 * Last modified: Mon Jan 30, 2017  23:27
 *
 * Description: Variation graph iterator for traversing the graph.
 *
 * Copyright (c) 2017, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari (cartoonist), ali.ghaffaari@mpi-inf.mpg.de
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#ifndef VARGRAPH_ITER_H__
#define VARGRAPH_ITER_H__

#include <unordered_set>
#include <deque>
#include <utility>

#include "vargraph.h"

namespace grem
{
  // forwards
  template<typename TGraph, typename TSpec>
    class Iterator;
  template<typename TGraph, typename TSpec>
    bool at_end(Iterator<TGraph, TSpec> it);
  template<typename TGraph, typename TSpec>
    Iterator<TGraph, TSpec> begin(const TGraph & g,
        typename TSpec::Value start=0);
  template<typename TGraph, typename TSpec>
    typename TSpec::Level level(Iterator<TGraph, TSpec> & it);

  template<typename TGraph, typename TSpec>
    class Iterator
    {
      friend bool at_end<TGraph, TSpec>(Iterator<TGraph, TSpec> it);
      friend Iterator<TGraph, TSpec> begin<TGraph, TSpec>(const TGraph & g,
                                                          typename TSpec::Value start);

      public:
        Iterator(const TGraph & vargraph, typename TSpec::Value start=0) :
          Iterator(&vargraph, start) {}
        Iterator(const TGraph * _vargraph_ptr, typename TSpec::Value start=0)
        {
          *this = begin<TGraph, TSpec>(*_vargraph_ptr, start);  // use move assign.
        }

        typename TSpec::Value operator*() { return this->itr_value; }
        typename TSpec::TContainer::value_type buffer_front();
        Iterator & operator++();
      private:
        const TGraph * vargraph_ptr;
        typename TSpec::Value itr_value;  // e.g. current node id pointed by iterator.
        typename TSpec::TContainer visiting_buffer;
        typename TSpec::TSet visited;

        Iterator() : vargraph_ptr(nullptr) {}
    };

  template<typename TSpec = void>
    struct BfsIterator
    {
      typedef id_t Value;
      typedef Value Level;
      typedef std::deque< std::pair< Value, Level >> TContainer;

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
    };

  template<>
    BfsIterator<>::TContainer::value_type
    Iterator<VarGraph, BfsIterator<>>::buffer_front()
    {
      return this->visiting_buffer.front();
    }

  template<>
    BfsIterator<>::Level level(Iterator<VarGraph, BfsIterator<>> & it)
    {
      return it.buffer_front().second;
    }

  template<>
    Iterator<VarGraph, BfsIterator<>> &
    Iterator<VarGraph, BfsIterator<>>::operator++()
    {
      BfsIterator<>::Level plevel = level(*this);
      if (this->vargraph_ptr->has_fwd_edge(this->itr_value))
      {
        auto edges = this->vargraph_ptr->fwd_edges(this->itr_value);
        for (auto it = edges.begin(); it != edges.end(); ++it)
        {
          BfsIterator<>::Value adj_node = (*it)->to();
          if (visited.find(std::make_pair(adj_node, 0)) == // level doesn't matter (see
              visited.end())                               //   method `pair_pred`).
          {
            this->visiting_buffer.push_back(
                std::make_pair(adj_node, plevel + 1));
            this->visited.insert(std::make_pair(adj_node, plevel + 1));
          }
        }
      }
      this->visiting_buffer.pop_front();
      this->itr_value = this->visiting_buffer.front().first;
      if (this->visiting_buffer.front().second != plevel)
      {
        for (auto vn_it = this->visited.begin();
             vn_it != this->visited.end();
             /* noop */)
        {
          if ((*vn_it).second == plevel)
          {
            this->visited.erase(vn_it++);
          }
          else
          {
            ++vn_it;
          }
        }
      }

      return *this;
    }

  template<>
  bool at_end(Iterator<VarGraph, BfsIterator<>> it)
  {
    return it.visiting_buffer.empty();
  }

  template<>
    Iterator<VarGraph, BfsIterator<>>
    begin(const VarGraph & g, BfsIterator<>::Value start)
    {
      Iterator<VarGraph, BfsIterator<>> begin_it;
      BfsIterator<>::Value start_node_id;
      if (start != 0) start_node_id = start;
      else start_node_id = g.node_at(0).id();

      begin_it.vargraph_ptr = &g;
      begin_it.visiting_buffer.push_back(std::make_pair(start_node_id, 0));
      begin_it.visited.insert(std::make_pair(start_node_id, 0));
      begin_it.itr_value = begin_it.visiting_buffer.front().first;

      return begin_it;
    }
}  // namespace grem

#endif  // VARGRAPH_ITER_H__

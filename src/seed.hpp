/**
 *    @file  seed.hpp
 *   @brief  Define the Seed class.
 *
 *  Seed class definition.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Aug 29, 2017  00:40
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_SEED_HPP__
#define PSI_SEED_HPP__

#include <cstdlib>

#include "base.hpp"


namespace grem {
  /**
   *  @brief  Seed class.
   *
   *  Represent a seed hit on the graph.
   */
  template< typename TId = std::size_t, typename TOffset = std::size_t >
  class Seed
  {
    public:
      /* ====================  TYPEDEFS      ======================================= */
      typedef TId id_type;
      typedef TOffset offset_type;

      /* ====================  DATA MEMBERS  ======================================= */
      TId node_id;                              /**< @brief Graph node ID. */
      TOffset node_offset;                      /**< @brief Graph node offset. */
      TId read_id;                              /**< @brief Read id. */
      TOffset read_offset;                      /**< @brief Read offset. */
  };  /* -----  end of class Seed  ----- */

}  /* -----  end of namespace grem  ----- */

#endif  /* --- #ifndef PSI_SEED_HPP__ --- */

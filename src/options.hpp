/**
 *    @file  options.hpp
 *   @brief  Options class definition.
 *
 *  It contains data structures storing program options.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Fri Nov 11, 2016  09:40
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_OPTIONS_HPP__
#define PSI_OPTIONS_HPP__

#include <stdexcept>
#include <string>

#include <seqan/index.h>


namespace psi {
  // :TODO:Tue Mar 07 20:34:\@cartoonist: Option class.
  enum class IndexType {
    Sa = 1,               /**< @brief Suffix array index. */
    Esa,                  /**< @brief Enhanced suffix array index. */
    Wotd,                 /**< @brief Lazy suffix tree (write only, top down) index. */
    Dfi,                  /**< @brief Deferred frequency index. */
    QGram,                /**< @brief An index based on an array of sorted q-grams. */
    FM                    /**< @brief FM index. */
  };

  typedef seqan::IndexWotd<> UsingIndexWotd;
  typedef seqan::IndexEsa<> UsingIndexEsa;

    inline IndexType
  index_from_str(std::string str)
  {
    if (str == "SA") return IndexType::Sa;
    if (str == "ESA") return IndexType::Esa;
    if (str == "WOTD") return IndexType::Wotd;
    if (str == "DFI") return IndexType::Dfi;
    if (str == "QGRAM") return IndexType::QGram;
    if (str == "FM") return IndexType::FM;

    throw std::runtime_error("Undefined index type.");
  }

    inline std::string
  index_to_str(IndexType index)
  {
    if (index == IndexType::Sa) return std::string("SA");
    if (index == IndexType::Esa) return std::string("ESA");
    if (index == IndexType::Wotd) return std::string("WOTD");
    if (index == IndexType::Dfi) return std::string("DFI");
    if (index == IndexType::QGram) return std::string("QGRAM");
    if (index == IndexType::FM) return std::string("FM");

    throw std::runtime_error("Undefined index type.");
  }

  typedef struct
  {
    unsigned int seed_len;
    unsigned int chunk_size;
    unsigned int step_size;
    unsigned int distance;
    unsigned int path_num;
    unsigned int context;
    unsigned int gocc_threshold;
    unsigned int dindex_min_ris;
    unsigned int dindex_max_ris;
    IndexType index;
    std::string rf_path;
    std::string fq_path;
    std::string output_path;
    std::string log_path;
    std::string pindex_path;
    bool patched;
    bool indexonly;
    bool nologfile;
    bool nolog;
    bool quiet;
    bool nocolor;
    bool verbose;
  } Options;
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_OPTIONS_HPP__ --- */

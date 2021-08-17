/**
 *    @file  release.hpp
 *   @brief  Release information.
 *
 *  This header file contains release information, such as version, short and long
 *  description, and other constant values.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sun Nov 13, 2016  00:08
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_RELEASE_HPP__
#define PSI_RELEASE_HPP__

#include "base.hpp"


namespace psi {
  const char * const VERSION = PSI_PROJECT_VERSION;         /**< @brief Version number. */
#ifdef PSI_GIT_REVISION
  const char * const REVISION = PSI_GIT_REVISION;           /**< @brief Git revision. */
#else
  const char * const REVISION = "";
#endif
  const char * const PACKAGE = "PSI";                       /**< @brief Package name. */
  const char * const SHORT_DESC = "Pan-genome Seed Index";  /**< @brief Short description. */
  /** @brief Long description. */
  const char * const LONG_DESC =
      "Fully-sensitive seed finder in sequence graphs. "
      "This is an implementation of the method introduced in:\n\n"
      "  Ghaffaari, A. & Marschall, T. "
      "Fully-sensitive seed finding in sequence graphs using a hybrid index. "
      "Bioinformatics 35, i81-i89 (2019).\n\n"
      "When using PSI library for a publication please cite above article.";
  /** @brief Banner to be printed on help and version messages. */
  const char * const BANNER =
      "\n"
      ".---..---..-.\n"
      "| |-' \\ \\ | |\n"
      "`-'  `---'`-'\n";
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_RELEASE_HPP__ --- */

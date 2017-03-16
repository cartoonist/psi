/**
 *    @file  utils.h
 *   @brief  Utility and helper functions.
 *
 *  This header file contains general utility and helper functions.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Mar 07, 2017  20:11
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  UTILS_H__
#define  UTILS_H__

#include <seqan/seq_io.h>

#include "sequence.h"

namespace grem {
  /**
   *  @brief  Read records from the input file into named string set.
   *
   *  @param[out]  records Named string set to store records in the input file.
   *  @param[in,out]  infile The input file.
   *  @param[in]  num_record Read this number of record from the input file.
   *
   *  A wrapper function for `seqan::readRecords` method to read the records into named
   *  string set.
   */
  inline void
    readRecords ( Dna5QRecords &records,
        seqan::SeqFileIn &infile, unsigned int num_record )
    {
      CharStringSet quals;
      seqan::readRecords ( records.str, records.id, quals, infile, num_record );
      assignQualities ( records.str, quals );
      return;
    }  /* -----  end of function readRecords  ----- */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef UTILS_H__  ----- */

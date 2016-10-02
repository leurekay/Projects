/*
 * Copyright (c) 2013, Robert Rueger <rueger@itp.uni-frankfurt.de>
 *
 * This file is part of hVMC.
 *
 * hVMC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hVMC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hVMC.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FPTYPE_H_INCLUDED
#define FPTYPE_H_INCLUDED

#ifdef USE_FP_DBLPREC
typedef double fptype;
#else
typedef float fptype;
#endif


#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

namespace Eigen {
  typedef Matrix<fptype, Dynamic, Dynamic> MatrixXfp;
  typedef Array<fptype, Dynamic, Dynamic> ArrayXfp;
  typedef Matrix<fptype, Dynamic, 1> VectorXfp;
}

#endif // FPTYPE_H_INCLUDED

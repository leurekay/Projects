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

#include "analysis.hpp"

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <cmath>

#include <boost/archive/text_iarchive.hpp>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "serialization_eigen.hpp"
#include "lattice.hpp"
#include "mccrun_prepare.hpp"

using namespace std;
namespace fs = boost::filesystem;
namespace ar = boost::archive;


void analysis_static_structure_factor( const Options& opts, const fs::path& dir )
{
  cout << ":: Calculating static structure factor" << endl;

  // check if we have the required data
  if ( !fs::exists( dir / "sim_res_nncorr.dat" ) ) {
    cout << "ERROR: Calculation of the static structure factor requires "
            "measurements of the density density correlation, but no such "
            "result file was found in the working directory" << endl;
    return;
  }

  // make lattice and all the q vectors that we will use in the Fourier transform
  const shared_ptr<Lattice>& lat = prepare_lattice( opts );
  const vector<Eigen::VectorXd>& allq = lat->get_qvectors();

  // read the density-density correlation data from file
  ifstream nn_file( ( dir / "sim_res_nncorr.dat" ).string() );
  ar::text_iarchive nn_archive( nn_file );
  Eigen::MatrixXd nn( lat->L, lat->L );
  nn_archive >> nn;

  // open a file to output our results to
  ofstream ssfac_file( ( dir / "ana_res_ssfac.txt" ).string() );

  // calculate the onsite part of the sum
  const double Nq_onsite = nn.diagonal().sum() / static_cast<double>( lat->L );

  for ( auto q = allq.begin(); q != allq.end(); ++q ) {
    double sum = 0.0;
    for ( Lattice::index l = 0; l < lat->L; ++l ) {
      for ( Lattice::index k = l + 1; k < lat->L; ++k ) {
        if ( lat->include_r_in_ssfac( l, k ) ) {
          sum += cos( q->dot( lat->r( 0, k ) - lat->r( 0, l ) ) ) * nn( l, k);
        }
      }
    }
    const double Nq = 2.0 / static_cast<double>( lat->L ) * sum + Nq_onsite;
    ssfac_file << q->transpose() << " " << Nq << endl;
    if ( opts.count("verbose") ) {
      cout << q->transpose() << " " << Nq << endl;
    }
  }
}

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

#include "obs_dkdkp.hpp"

#if VERBOSE >= 1
# include <iostream>
#endif

#include <boost/mpi/collectives.hpp>

#include "macros.h"
#include "serialization_eigen.hpp"

using namespace std;
namespace mpi = boost::mpi;


ObservableDeltaKDeltaKPrime::ObservableDeltaKDeltaKPrime(
  unsigned int num_vpar, unsigned int optimizers_init )
  : optimizers( optimizers_init ),
    thisbin_DkDkp_sum( Eigen::MatrixXd::Zero( num_vpar, num_vpar ) ),
    thisbin_count( 0 ),
    binmean_DkDkp_sum( Eigen::MatrixXd::Zero( num_vpar, num_vpar ) ),
    binmean_count( 0 ) { }


void ObservableDeltaKDeltaKPrime::measure(
  const ModelManager& model, ObservableCache& cache )
{
  if ( !cache.DeltaK ) {
    cache.DeltaK = model.Delta_k( optimizers );
  }
  const Eigen::VectorXd& Dk_current = cache.DeltaK.get();

  thisbin_DkDkp_sum += Dk_current * Dk_current.transpose();
  ++thisbin_count;

#if VERBOSE >= 1
  cout << "ObservableDeltaKDeltaKPrime::measure() : thisbin_DkDkp_sum = " << endl
       << thisbin_DkDkp_sum << endl;
#endif
}


void ObservableDeltaKDeltaKPrime::completebin()
{
  binmean_DkDkp_sum += thisbin_DkDkp_sum / static_cast<double>( thisbin_count );
  ++binmean_count;

  thisbin_DkDkp_sum.setZero();
  thisbin_count = 0;
}


void ObservableDeltaKDeltaKPrime::collect_and_write_results(
  const mpi::communicator& mpicomm,
  MCCResults& results ) const
{
  assert( mpicomm.rank() == 0 );

  vector<Eigen::MatrixXd> binmeans_collector(
      mpicomm.size(),
      Eigen::MatrixXd( binmean_DkDkp_sum.rows(), binmean_DkDkp_sum.cols() )
  );
  mpi::gather( mpicomm, binmean_DkDkp_sum, binmeans_collector, 0 );
  vector<unsigned int> binmeans_collector_count;
  mpi::gather( mpicomm, binmean_count, binmeans_collector_count, 0 );

  results.Deltak_Deltakprime
    = accumulate(
        binmeans_collector.begin() + 1,
        binmeans_collector.end(),
        binmeans_collector.front()
      ) / static_cast<double>(
        accumulate(
          binmeans_collector_count.begin(),
          binmeans_collector_count.end(),
          0
        )
      );
}


void ObservableDeltaKDeltaKPrime::send_results_to_master(
  const mpi::communicator& mpicomm ) const
{
  assert( mpicomm.rank() != 0 );

  mpi::gather( mpicomm, binmean_DkDkp_sum, 0 );
  mpi::gather( mpicomm, binmean_count, 0 );
}

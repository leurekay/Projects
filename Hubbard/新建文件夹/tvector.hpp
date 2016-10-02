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

#ifndef TVECTOR_H_INCLUDED
#define TVECTOR_H_INCLUDED

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "lattice.hpp"
#include "pconf.hpp"
#include "jastrow.hpp"
#include "fpctrl.hpp"


class TVector final
{
  protected:

    const Lattice* const lat;
    const Jastrow& v;
    const ParticleConfiguration& pconf;

    Eigen::VectorXd T;

    const unsigned int updates_until_recalc;
    unsigned int updates_since_recalc;
    FPDevStat devstat;

    Eigen::VectorXd calc_new() const;
    Eigen::VectorXd calc_qupdated( const ParticleHop& hop ) const;

  public:

    TVector(
      const Lattice* lat_init,
      const Jastrow& v_init,
      const ParticleConfiguration& pconf_init,
      double deviation_target,
      unsigned int updates_until_recalc_init
    );

    void init();
    void update( const ParticleHop& hop );

    const Eigen::VectorXd& get() const;

    FPDevStat get_devstat() const;

};

#endif // TVECTOR_H_INCLUDED

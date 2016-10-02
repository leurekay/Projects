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

#include "options.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>

#include "analysis.hpp"
#include "obs.hpp"
#include "lattice.hpp"
#include "utils.hpp"

using namespace std;
namespace po  = boost::program_options;
namespace fs  = boost::filesystem;


Options read_options( int argc, char* argv[], bool is_master )
{
  Options vm;

  // define command line only options

  po::options_description clionly( "command line options" );
  clionly.add_options()
  ( "help,h", "print this help message and exit" )
  ( "version,V", "print hVMC's version and exit" )
  ( "verbose,v", "makes hVMC write additional information to stdout and disk" )
  ( "job-file,J", po::value<fs::path>(), "job file to execute" );
  po::positional_options_description p;
  p.add( "job-file", -1 );


  // define command line and jobfile options

  po::options_description physparam( "physical parameters" );
  physparam.add_options()

  ( "phys.nn-hopping,1",
    po::value<double>()->required(),
    "nearest neighbor hopping matrix element t" )

  ( "phys.2nd-nn-hopping,2",
    po::value<double>()->default_value( 0.0 ),
    "2nd nearest neighbor hopping matrix element t'" )

  ( "phys.3rd-nn-hopping,3",
    po::value<double>()->default_value( 0.0 ),
    "3rd nearest neighbor hopping matrix element t''" )

  ( "phys.onsite-energy,U",
    po::value<double>()->required(),
    "on-site energy U" )

  ( "phys.lattice,l",
    po::value<Lattice::type>()->required(),
    "lattice type (1dchain, 2dsquare, 2dsquare2layer)" )

  ( "phys.num-lattice-sites,L",
    po::value<unsigned int>()->required(),
    "number of lattice sites" )

  ( "phys.num-electrons,N",
    po::value<unsigned int>()->required(),
    "total number of electrons" )

  ( "phys.pairing-symmetry,y",
    po::value<optpairsym_t>()->default_value( OPTION_PAIRING_SYMMETRY_SWAVE ),
    "symmetry of the pairing term (swave, dwave, dwave_twisted)" )

  ( "phys.vpar-file,P",
    po::value<fs::path>(),
    "file to load the variational parameters from" )

  ( "phys.vpar-ovwrt-t2",
    po::value<double>(),
    "variational parameter overwrite: t'" )

  ( "phys.vpar-ovwrt-t3",
    po::value<double>(),
    "variational parameter overwrite: t''" )

  ( "phys.vpar-ovwrt-D0",
    po::value<double>(),
    "variational parameter overwrite: Delta_onsite" )

  ( "phys.vpar-ovwrt-D1",
    po::value<double>(),
    "variational parameter overwrite: Delta" )

  ( "phys.vpar-ovwrt-D2",
    po::value<double>(),
    "variational parameter overwrite: Delta'" )

  ( "phys.vpar-ovwrt-D3",
    po::value<double>(),
    "variational parameter overwrite: Delta""'" )

  ( "phys.vpar-ovwrt-mu",
    po::value<double>(),
    "variational parameter overwrite: mu" )

  ( "phys.vpar-ovwrt-mu-m",
    po::value<double>(),
    "variational parameter overwrite: mu_m" )

  ( "phys.vpar-ovwrt-J0",
    po::value<double>(),
    "variational parameter overwrite: J_00" );

  po::options_description calcset( "calculation settings" );
  calcset.add_options()

  ( "calc.mode,m",
    po::value<optmode_t>()->required(),
    "mode (optimization, simulation, analysis)" )

  ( "calc.working-dir,D",
    po::value<fs::path>()->default_value( "." ),
    "[opt+sim+ana]: input/output directory" )

  ( "calc.observable,O",
    po::value< std::vector<observables_t> >(),
    "[sim]: measured observables "
    "(E, Dk, DkDkp, DkE, dblocc, nncorr, sscorr, pconfs)" )

  ( "calc.analysis,A",
    po::value< std::vector<analysis_t> >(),
    "[ana]: selected analysis modules (ssfac)" )

  ( "calc.update-hop-maxdistance,H",
    po::value<unsigned int>()->default_value( 1 ),
    "[opt+sim]: maximum hopping distance for electronic configuration updates" )

  ( "calc.num-mcs-equil,E",
    po::value<unsigned int>()->default_value( 100 ),
    "[opt+sim]: number of Monte Carlo steps for equilibration" )

  ( "calc.num-bins,B",
    po::value<unsigned int>()->default_value( 50 ),
    "[opt+sim]: number of measurement bins" )

  ( "calc.num-binmcs,M",
    po::value<unsigned int>()->default_value( 50 ),
    "[opt+sim]: number of Monte Carlo steps per bin" )

  ( "calc.rng-seed,S",
    po::value<unsigned int>(),
    "[opt+sim]: random number generator seed" )

  ( "calc.optimizers,Z",
    po::value<unsigned int>()->default_value( 511 ),
    "[opt]: which variational parameters to optimize" )

  ( "calc.sr-dt,d",
    po::value<double>()->default_value( 1.0 ),
    "[opt]: controls the SR convergence: vpar += dt * dvpar" )

  ( "calc.sr-dt-Jboost,J",
    po::value<double>()->default_value( 1.0 ),
    "[opt]: factor to boost the convergence speed of the Jastrow" )

  ( "calc.sr-mkthreshold,T",
    po::value<double>()->default_value( 0.5 ),
    "[opt]: Mann-Kendall threshold for convergence detection" )

  ( "calc.sr-max-refinements,R",
    po::value<unsigned int>()->default_value( 4 ),
    "[opt]: number of refinements during the SR cycle" )

  ( "calc.sr-averaging-cycles,A",
    po::value<unsigned int>()->default_value( 10 ),
    "[opt]: number of SR cycles to average the converged variational parameters" )

  ( "calc.vpar-minabs-select",
    po::value<unsigned int>()->default_value( 0 ),
    "[opt]: variational parameters that have a fixed minimum absolute value" )

  ( "calc.vpar-minabs-value",
    po::value<double>()->default_value( 0.001 ),
    "[opt]: fixed minimum absolute value for the selected parameters" );

  po::options_description fpctrl( "[opt+sim]: floating point precision control" );
  fpctrl.add_options()

  ( "fpctrl.W-deviation-target",
    po::value<double>()->default_value( 0.001, "0.001" ),
    "deviation target for the matrix W" )

  ( "fpctrl.W-updates-until-recalc",
    po::value<unsigned int>()
#ifdef USE_FP_DBLPREC
      ->default_value( 5000 ),
#else
      ->default_value( 500 ),
#endif
    "number of quick updates until recalculation of the matrix W" )

  ( "fpctrl.T-deviation-target",
    po::value<double>()->default_value( 0.001, "0.001" ),
    "deviation target for the vector T" )

  ( "fpctrl.T-updates-until-recalc",
    po::value<unsigned int>()->default_value( 500000 ),
    "number of quick updates until recalculation of the vector T" );

  // define option groups for cli and jobfile
  po::options_description cmdline_options;
  cmdline_options.add( clionly ).add( physparam ).add( calcset ).add( fpctrl );
  po::options_description jobfile_options;
  jobfile_options.add( physparam ).add( calcset ).add( fpctrl );



  try {
    // parse options from the command line
    if ( is_master ) {
      cout << ":: Parsing command line ..." << endl;
    }
    po::store( po::command_line_parser( argc, argv ).
               options( cmdline_options ).positional( p ).run(), vm );
  } catch ( const po::error& e ) {
    if ( is_master ) {
      cerr << "Error while parsing the command line: " << e.what() << endl;
    }
    exit( 1 );
  }


  // display help or hVMC version information

  if ( vm.count( "help" ) ) {
    if ( is_master ) {
      cout << endl;
      cout << "usage: hVMC [OPTIONS] JOBFILE -o OUTDIR" << endl;
      cout << cmdline_options << endl;
    }
    exit( 0 );
  }

  if ( vm.count( "version" ) ) {
    if ( is_master ) {
      cout << endl;

      cout << "hVMC - built from git commit " << GIT_HASH << endl;

#ifndef NDEBUG
      cout << "=========== DEBUG BUILD ============" << endl;
#else
      cout << "========== RELEASE BUILD ===========" << endl;
#endif

      cout << "compiled " << __DATE__ " with ";
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#  if defined(__GNUC_PATCHLEVEL__)
      cout << "GCC " << __GNUC__ << "." << __GNUC_MINOR__
           << "." << __GNUC_PATCHLEVEL__;
#  else
      cout << "GCC " << __GNUC__ << "." << __GNUC_MINOR__;
#  endif
#elif defined(__INTEL_COMPILER)
      cout << "ICC " << __INTEL_COMPILER;
#elif defined(_MSC_VER)
      cout << "MSC " << _MSC_VER;
#else
      cout << "unknown compiler";
#endif
      cout << endl;

#ifdef USE_FP_DBLPREC
      cout << "WMatrix implementation: double prec." << endl;
#else
      cout << "WMatrix implementation: single prec." << endl;
#endif
#ifdef USE_CBLAS
      cout << "external CBLAS: enabled" << endl;
#else
      cout << "external CBLAS: disabled" << endl;
#endif
#ifdef EIGEN_DEFAULT_TO_ROW_MAJOR
      cout << "matrix storage order: row major" << endl;
#else
      cout << "matrix storage order: column major" << endl;
#endif
      cout << endl;

      cout
          << "Copyright (C) 2013, Robert Rueger <rueger@itp.uni-frankfurt.de>"
          << endl
          << "License GPLv3+: GNU GPL version 3 or later"
          " <http://gnu.org/licenses/gpl.html>" << endl
          << "This is free software: you are free to change and redistribute it."
          << endl
          << "There is NO WARRANTY, to the extent permitted by law."
          << endl << endl;
    }
    exit( 0 );
  }


  if ( vm.count( "job-file" ) ) {
    // parse the jobfile
    if ( is_master ) {
      cout << ":: Parsing jobfile ..." << endl;
    }
    ifstream jobifs( vm["job-file"].as<fs::path>().string() );
    if ( jobifs.is_open() ) {
      try {
        po::store( po::parse_config_file( jobifs, jobfile_options ), vm );
      } catch ( const po::error& e ) {
        if ( is_master ) {
          cerr << "Error while parsing the job file: " << e.what() << endl;
        }
        exit( 1 );
      }
    } else {
      if ( is_master ) {
        cerr << "Error: unable to open jobfile " << vm["job-file"].as<string>()
             << endl;
      }
      exit( 1 );
    }
  }

  try {
    // finalize the variable map
    // (will throw exception on missing options)
    po::notify( vm );
  } catch ( const po::error& e ) {
    if ( is_master ) {
      cerr << "Error in program options: " << e.what() << endl;
    }
    exit( 1 );
  }


  // manual modifications of the options

  // set phys.vpar-file to working dir / opt_vpar_final.dat
  // (can't be a default value because the working dir is not known
  // before the command line is parsed)
  // (insert will do nothing if vpar-file is already set ...)
  vm.insert(
    make_pair(
      "phys.vpar-file",
      po::variable_value(
        vm["calc.working-dir"].as<fs::path>() / "opt_vpar_final.dat",
        true
      )
    )
  );


  // check for logical errors in the physical parameters
  try {

    if ( vm["phys.num-electrons"].as<unsigned int>() % 2 != 0 ) {
      throw logic_error( "electron number must be even" );
    }

    if ( vm["phys.num-electrons"].as<unsigned int>() >
         2 * vm["phys.num-lattice-sites"].as<unsigned int>() ) {
      throw logic_error(
        "too many electrons (num-electrons > 2 * num-lattice-sites) "
      );
    }

    if ( vm["phys.lattice"].as<Lattice::type>() == Lattice::type::square2d &&
         !is_perfect_square( vm["phys.num-lattice-sites"].as<unsigned int>() ) ) {
      throw logic_error( "the number of lattice sites must be a perfect square" );
    }

    if ( vm["phys.lattice"].as<Lattice::type>() == Lattice::type::square2d2layer &&
         !is_perfect_square( vm["phys.num-lattice-sites"].as<unsigned int>() / 2 ) ) {
      throw logic_error(
              "the number of lattice sites must be two times a perfect square"
            );
    }

    if ( vm["phys.lattice"].as<Lattice::type>() == Lattice::type::square2d2layer &&
         vm["phys.3rd-nn-hopping"].as<double>() == 0.0 ) {
      throw logic_error(
              "running the square lattice bilayer model with no interplane hopping does not make sense, use the monolayer instead"
            );
    }

    // TODO: minimum lattice size checks (Robert Rueger, 2012-11-17 22:57)

  } catch ( const logic_error& e ) {
    if ( is_master ) {
      cerr << "Logical error in physical parameters: " << e.what() << endl;
    }
    exit( 1 );
  }


  // check for logical errors in the calculation settings
  try {

    if ( vm["calc.update-hop-maxdistance"].as<unsigned int>() == 0 ) {
      throw logic_error(
        "electronic configuration updates need at least nearest neighbor hopping"
      );
    }

    if ( vm["calc.update-hop-maxdistance"].as<unsigned int>() > 3 ) {
      throw logic_error(
        "electronic configuration updates with hopping > 3rd nearest neighbors"
        "are not supported"
      );
    }

    if ( vm["calc.mode"].as<optmode_t>() == OPTION_MODE_OPTIMIZATION &&
         vm["calc.optimizers"].as<unsigned int>() == 0 ) {
      throw logic_error( "you need to optimize at least one variational parameter" );
    }

    if ( vm["calc.mode"].as<optmode_t>() == OPTION_MODE_OPTIMIZATION &&
         vm["calc.optimizers"].as<unsigned int>() > 511 ) {
      throw logic_error(
        "hVMC only has 8+Jastrow variational parameters "
        "-> calc.optimizers can not be larger than 511"
      );
    }

    if ( vm["calc.mode"].as<optmode_t>() == OPTION_MODE_OPTIMIZATION &&
         vm["calc.vpar-minabs-select"].as<unsigned int>() > 255 ) {
      throw logic_error(
        "the minimum vpar threshold is only intended for the det. parameters "
        "-> calc.vpar-minabs-select can not be larger than 255"
      );
    }

    if ( vm["calc.mode"].as<optmode_t>() == OPTION_MODE_SIMULATION &&
         vm.count( "calc.observable" ) == 0 ) {
      throw logic_error( "you need to measure at least one observable" );
    }

    if ( vm["calc.mode"].as<optmode_t>() == OPTION_MODE_ANALYSIS &&
         vm.count( "calc.analysis" ) == 0 ) {
      throw logic_error( "you need to select at least one analysis module" );
    }

  } catch ( const logic_error& e ) {
    if ( is_master ) {
      cerr << "Logical error in calculation settings: " << e.what() << endl;
    }
    exit( 1 );
  }

  if ( is_master ) {
    cout << endl;
  }

  return vm;
}


istream& operator>>( std::istream& in, Lattice::type& lat )
{
  string token;
  in >> token;
  if ( token == "1dchain" ) {
    lat = Lattice::type::chain1d;
  } else if ( token == "2dsquare" ) {
    lat = Lattice::type::square2d;
  } else if ( token == "2dsquare2layer" ) {
    lat = Lattice::type::square2d2layer;
  } else {
    throw po::validation_error( po::validation_error::invalid_option_value );
  }
  return in;
}


istream& operator>>( std::istream& in, optmode_t& m )
{
  string token;
  in >> token;
  if ( token == "opt" || token == "optimization" ) {
    m = OPTION_MODE_OPTIMIZATION;
  } else if ( token == "sim" || token == "simulation" ) {
    m = OPTION_MODE_SIMULATION;
  } else if ( token == "ana" || token == "analysis" ) {
    m = OPTION_MODE_ANALYSIS;
  } else {
    throw po::validation_error( po::validation_error::invalid_option_value );
  }
  return in;
}

istream& operator>>( std::istream& in, optpairsym_t& s )
{
  string token;
  in >> token;
  if ( token == "s" || token == "swave" ) {
    s = OPTION_PAIRING_SYMMETRY_SWAVE;
  } else if ( token == "d" || token == "dwave" ) {
    s = OPTION_PAIRING_SYMMETRY_DWAVE;
  } else if ( token == "dt" || token == "dwave_twisted" ) {
    s = OPTION_PAIRING_SYMMETRY_DWAVE_TWISTED;
  } else {
    throw po::validation_error( po::validation_error::invalid_option_value );
  }
  return in;
}

istream& operator>>( std::istream& in, analysis_t& a )
{
  string token;
  in >> token;
  if ( token == "ssfac" ) {
    a = ANALYSIS_STATIC_STRUCTURE_FACTOR;
  } else {
    throw po::validation_error( po::validation_error::invalid_option_value );
  }
  return in;
}

istream& operator>>( std::istream& in, observables_t& obs )
{
  string token;
  in >> token;
  if ( token == "E" ) {
    obs = OBSERVABLE_E;
  } else if ( token == "Dk" ) {
    obs = OBSERVABLE_DELTAK;
  } else if ( token == "DkDkp" ) {
    obs = OBSERVABLE_DELTAK_DELTAKPRIME;
  } else if ( token == "DkE" ) {
    obs = OBSERVABLE_DELTAK_E;
  } else if ( token == "dblocc" ) {
    obs = OBSERVABLE_DOUBLE_OCCUPANCY_DENSITY;
  } else if ( token == "nncorr" ) {
    obs = OBSERVABLE_DENSITY_DENSITY_CORRELATION;
  } else if ( token == "sscorr" ) {
    obs = OBSERVABLE_SPIN_SPIN_CORRELATION;
  } else if ( token == "pconfs" ) {
    obs = OBSERVABLE_PARTICLE_CONFIGURATIONS;
  } else {
    throw po::validation_error( po::validation_error::invalid_option_value );
  }
  return in;
}

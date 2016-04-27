/**
 * \file ECCQECalculator.h
 *
 * \ingroup Utilities
 *
 * \brief A static function to calculate CCQE energy.
 *
 * @author davidkaleko
 */

/** \addtogroup Utilities

    @{*/
#ifndef LEE_ECCQECALCULATOR_H
#define LEE_ECCQECALCULATOR_H

#include <vector>
#include <iostream>
#include <math.h> //pow
#include "TMath.h"
#include "ERTool/Base/AnaBase.h"

/**
   \class ECCQECalculator
   A static function to calculate CCQE energy. Returns in units of GEV.
 */
namespace lee {
  namespace util {

    class ECCQECalculator {

    public:

      /// Default constructor
      ECCQECalculator() {};

      /// Default destructor
      virtual ~ECCQECalculator() {};

      /// Method using only truth momentum 4-vector (wrapper)
      /// Momentum vector should be in MeV
      static double ComputeECCQE(const std::vector<double> &lepton_4momentum);

      /// Method using manually-input energy (IE if you smear energy first)
      /// Energy should be in MeV, direction can be (doesn't have to be) unit-normalized
      static double ComputeECCQE(double energy, const std::vector<double> &lepton_dir, bool is_electron = true);

      /// Method using ERTool Shower
      static double ComputeECCQE(const ::ertool::Shower &ertshower);

    };
  }// end namespace util
}// end namespace ubsens
#endif
/** @} */ // end of doxygen group


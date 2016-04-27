#ifndef LEE_ECCQECALCULATOR_CXX
#define LEE_ECCQECALCULATOR_CXX

#include "ECCQECalculator.h"
namespace lee {
  namespace util {

    double ECCQECalculator::ComputeECCQE(const std::vector<double> &lepton_4momentum) {

      double l_energy = -999.;
      try {l_energy = lepton_4momentum.at(3);}
      catch (...) {
        std::cout << "Exception within ECCQECalculator. "
                  << "Did you feed in the right size 4-momentum? "
                  << "Returning -999. CCQE energy to punish your mistakes."
                  << std::endl;
      }

      std::vector<double> lepton_dir;
      for (size_t i = 0; i < 3; i++)
        lepton_dir.push_back(lepton_4momentum.at(i));

      return ComputeECCQE(l_energy, lepton_dir);

    }

    double ECCQECalculator::ComputeECCQE(double totalenergy, const std::vector<double> &lepton_dir, bool is_electron) {

      if ( lepton_dir.size() != 3 ) {
        std::cerr << "From ComputeECCQE: input direction vector doesn't have size 3! Quitting..." << std::endl;
        return -99.;
      }

      ///.at(0) is x momentum in MeV/C

      double M_n = 939.565;    // MeV/c2
      double M_p = 938.272;    // MeV/c2
      double leptonmass = 0.511;      // MeV/c2
      if(!is_electron) leptonmass = 105.6583;
      double bindingE = 30.0;  // MeV

      double l_energy = totalenergy;
      double l_mom = pow(pow(l_energy, 2) - pow(leptonmass, 2), 0.5);

      // Only truth info goes into theta calculation
      double l_theta =
        TMath::ACos(lepton_dir.at(2) /
                    pow(
                      (
                        pow(lepton_dir.at(0), 2) +
                        pow(lepton_dir.at(1), 2) +
                        pow(lepton_dir.at(2), 2)
                      ), 0.5
                    )
                   );

      double nu_energy_num = pow(M_p, 2) - pow(M_n - bindingE, 2)
                             - pow(leptonmass, 2) + 2.0 * (M_n - bindingE) * l_energy;
      double nu_energy_den = 2.0 * (M_n - bindingE - l_energy + l_mom * TMath::Cos(l_theta));

      // For a result in GEV, divide by 1000.
      return (nu_energy_num / nu_energy_den) / 1000.;

    }

    double ECCQECalculator::ComputeECCQE(const ::ertool::Shower &ertshower){

      return ComputeECCQE( ertshower._energy, ertshower.Dir() );
    
    }
  }// end namespace util
}// end namespace ubsens
#endif

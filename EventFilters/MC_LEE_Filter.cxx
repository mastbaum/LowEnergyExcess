#ifndef LARLITE_MC_LEE_FILTER_CXX
#define LARLITE_MC_LEE_FILTER_CXX

#include "MC_LEE_Filter.h"

namespace larlite {

  bool MC_LEE_Filter::initialize() {

    total_events = 0;
    kept_events = 0;

    return true;

  }

  bool MC_LEE_Filter::analyze(storage_manager* storage) {

    //Grab the MCTruth
    auto ev_mctruth = storage->get_data<event_mctruth>("generator");
    if (!ev_mctruth) {
      print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
      return false;
    }

    total_events++;


    //Enforce that there is exactly 1 electron, above 20MeV kinetic energy
    //Don't care about neutrons, protons.
    //Don't care about other particles if they are below 20MeV KE
    //I don't care about neutrons, weird quarks, the neutrino, etc.
    size_t n_electrons = 0;

    for (auto const& particle : ev_mctruth->at(0).GetParticles()) {

      // Only particles with status code 1 are relevant
      if ( particle.StatusCode() != 1 ) continue;

      //Note: this KE is in units of GEV!
      double KE = particle.Trajectory().at(0).E() - particle.Mass();

      //Don't care about any particles with less than 20 MeV KE
      if ( KE < 0.02 ) continue;

      //Count up the number of electrons
      if ( abs(particle.PdgCode()) == 11 )
        n_electrons++;

      //Skip this event if any other particles exist... this filter
      //is to pick out 1eNpNn events!
      // 22  => gammas
      // 211 => charged pions
      // 111 => pi0s
      // 13  => muons
      // 321 => kaons
      if ( abs(particle.PdgCode()) == 22  ||
           abs(particle.PdgCode()) == 211 ||
           abs(particle.PdgCode()) == 111 ||
           abs(particle.PdgCode()) == 13  ||
           abs(particle.PdgCode()) == 321 )
        return false;

      /*
      //temp cout to check for unexpected particles
      if (abs(particle.PdgCode()) != 2212 && //proton
      abs(particle.PdgCode()) != 2112 && //neutron
      abs(particle.PdgCode()) != 11  &&  //electron
      abs(particle.PdgCode()) != 12  &&  //nue
      abs(particle.PdgCode()) != 2000000101 && //bindino
      abs(particle.PdgCode()) != 311  &&  //K0
      abs(particle.PdgCode()) != 3222  &&  //sigma+
      abs(particle.PdgCode()) != 3122 )  //lamda+
      std::cout<<"Found particle with pdg = "<<particle.PdgCode()<<std::endl;
      */
    }

    //Skip this event if the number of electrons is not 1
    if ( n_electrons != 1 )
      return false;

    kept_events++;
    return true;
  }

  bool MC_LEE_Filter::finalize() {

    std::cout << "MC_LEE_Filter: Total events = " << total_events << std::endl;
    std::cout << "MC_LEE_Filter: Final kept events = " << kept_events << std::endl;

    return true;
  }

}
#endif

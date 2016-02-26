#ifndef LARLITE_TEST_LEERW_CXX
#define LARLITE_TEST_LEERW_CXX

#include "test_LEERW.h"
#include "DataFormat/mctruth.h"

namespace larlite {

bool test_LEERW::initialize() {

    _rw.set_debug(_debug);

    _rw.set_source_filename("$LARLITE_USERDEVDIR/LowEnergyExcess/LEEReweight/source/whatttt/LEE_Reweight_plots.root");
    // _rw.set_source_filename("$LARLITE_USERDEVDIR/LowEnergyExcess/LEEReweight/source/whatttt/LEE_Reweight_plots_mcc7temp.root");
    // _rw.set_source_filename("$LARLITE_USERDEVDIR/LowEnergyExcess/LEEReweight/source/fuckfuck.root");

    // _rw.set_n_generated_events(6637);
    _rw.set_n_generated_events(369);

    _rw.set_generated_evis_uz_corr_name("initial_evis_uz_corr");
    // _rw.set_generated_evis_uz_corr_name("temp_mcc7_lowstat");

    _rw.initialize();

    initial_nu_spectrum = new TH1F("initial_nu_spectrum", "Raw Electron-Neutrino Energy Spectrum", 100, 0, 3);
    reweighted_nu_spectrum = new TH1F("reweighted_nu_spectrum", "Reweighted Electron-Neutrino Energy Spectrum", 100, 0, 3);
    initial_electron_spectrum = new TH1F("initial_electron_spectrum", "Raw Electron Energy Spectrum", 30, 0, 3);
    reweighted_electron_spectrum = new TH1F("reweighted_electron_spectrum", "Reweighted Electron Energy Spectrum", 30, 0, 3);
    initial_evis_uz_corr = new TH2D("initial_evis_uz_corr", "Raw Electron Uz vs. Evis Plot", 19, 100, 2000, 10, -1, 1);
    sculpted_evis_uz_corr = new TH2D("sculpted_evis_uz_corr", "Raw Electron Uz vs. Evis Plot", 19, 100, 2000, 10, -1, 1);
  
    summed_sculpting_weight = 0.;

    return true;
}

bool test_LEERW::analyze(storage_manager* storage) {

    auto ev_mctruth = storage->get_data<event_mctruth>("generator");
    if (!ev_mctruth) {
        print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, MCtruth!"));
        return false;
    }
    if (!ev_mctruth->size())
        return false;

    ::lee::EventInfo_t hacky_evt_info = _rw.extract_event_info(&(ev_mctruth->at(0)));

    double true_nu_energy = hacky_evt_info.nue_energy_GEV;
    initial_nu_spectrum->Fill(true_nu_energy);
    initial_evis_uz_corr->Fill(hacky_evt_info.electron_energy_MEV, hacky_evt_info.electron_uz);
    initial_electron_spectrum->Fill(hacky_evt_info.electron_energy_MEV / 1000.);
    double sculpting_weight = _rw.get_sculpting_weight(&(ev_mctruth->at(0)));
    summed_sculpting_weight += sculpting_weight;
    // std::cout << "(non-normalized) sculpting weight is computed to be: " << sculpting_weight << std::endl;

    double normalized_weight = _rw.get_normalized_weight(&(ev_mctruth->at(0)));
    // std::cout << "(normalized) scaling weight is computed to be: " << normalized_weight << std::endl;

    sculpted_evis_uz_corr->Fill(hacky_evt_info.electron_energy_MEV, hacky_evt_info.electron_uz, sculpting_weight);
    reweighted_nu_spectrum->Fill(true_nu_energy, normalized_weight * sculpting_weight);
    reweighted_electron_spectrum->Fill(hacky_evt_info.electron_energy_MEV / 1000., normalized_weight * sculpting_weight);
    return true;
}

bool test_LEERW::finalize() {
    std::cout<<"test_LEERW::finalize: summed sculpt weight is "<<summed_sculpting_weight<<std::endl;
    std::cout<<"test_LEERW::reweighted_nu_spectrum integral is "<<reweighted_nu_spectrum->Integral()<<std::endl;
    if (_fout) {
        _fout->cd();
        initial_nu_spectrum->Write();
        // reweighted_nu_spectrum->Scale( 1. / summed_sculpting_weight);
        reweighted_nu_spectrum->Write();
        initial_electron_spectrum->Write();
        // reweighted_electron_spectrum->Scale( 1. / summed_sculpting_weight);
        reweighted_electron_spectrum->Write();
        initial_evis_uz_corr->Write();
        // sculpted_evis_uz_corr->Scale( 1. / summed_sculpting_weight);
        sculpted_evis_uz_corr->Write();
    }
    return true;
}

}
#endif

#ifndef LARLITE_GENLEERWSOURCEPLOT_CXX
#define LARLITE_GENLEERWSOURCEPLOT_CXX

#include "GenLEERWSourcePlot.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctruth.h"

namespace larlite {

    bool GenLEERWSourcePlot::initialize() {
        if (_plot_name.empty()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not set plot name in GenLEERWSourcePlot!"));
            return false;
        }
        _out_hist = new TH2D(_plot_name.c_str(), "Raw Electron Uz vs. Evis Plot", 19, 100, 2000, 10, -1, 1);
        return true;
    }

    bool GenLEERWSourcePlot::analyze(storage_manager* storage) {

        //Grab the MCShowers
        auto ev_mcshower = storage->get_data<event_mcshower>("mcreco");
        if (!ev_mcshower) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mcshower!"));
            return false;
        }

        //Grab the MCTruth
        auto ev_mctruth = storage->get_data<event_mctruth>("generator");
        if (!ev_mctruth) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
            return false;
        }

        if (ev_mcshower->size() != 1) return false;
        if (abs(ev_mctruth->at(0).GetNeutrino().Nu().PdgCode()) != 12) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("You handed a non-nue file to GenLEERWSourcePlot!! What the heck?!"));
            return false;
        }

        auto const& mcs = ev_mcshower->at(0);

        double uz = mcs.Start().Momentum().Vect().Unit().Z();
        double evis = mcs.Start().E();
        _out_hist->Fill(evis, uz);

        return true;
    }

    bool GenLEERWSourcePlot::finalize() {

        if (_fout) {
            _fout->cd();
            _out_hist->Write();
        }
        return true;
    }

}
#endif

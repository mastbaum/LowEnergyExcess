#ifndef LEERW_CXX
#define LEERW_CXX

#include "LEERW.h"

namespace lee {

	bool LEERW::initialize() {

		if (_source_filename.empty())
			throw std::runtime_error("LEERW needs to be told the input root file (containing scaling graphs and histos).");

		if (!_n_generated_evts)
			throw std::runtime_error("LEERW needs to be told how many nue events were generated.");

		if(_generated_evis_uz_corr_name.empty())
			throw std::runtime_error("LEERW needs to be told the name of an input histogram!");
		// "_n_generated_evts" is the number of events generated in the entire cryostat.
		// If you generated events only inside of the TPC, you need to take into account the ratio
		// of entire cryostat volume to TPC volume
		// For now, I have the ratio of volumes as 171 tons / 86 tons = 1.99 ... this should probably be
		// replaced with a more precise value.
		if ( _events_generated_only_in_TPC ) _n_generated_evts /= 1.99;

		util::PlotReader::GetME()->SetFileName(_source_filename.c_str());
		util::PlotReader::GetME()->SetObjectName(_flux_ratio_name.c_str());
		util::PlotReader::GetME()->GetObject(_flux_ratio);
		util::PlotReader::GetME()->SetObjectName(_xsec_ratio_name.c_str());
		util::PlotReader::GetME()->GetObject(_xsec_ratio);
		util::PlotReader::GetME()->SetObjectName(_MB_evis_uz_corr_name.c_str());
		util::PlotReader::GetME()->GetObject(_MB_evis_uz_corr);
		util::PlotReader::GetME()->SetObjectName(_generated_evis_uz_corr_name.c_str());
		util::PlotReader::GetME()->GetObject(_generated_evis_uz_corr);
		//normalize evis correlation plot to unit area
		//so its integral is total # of expected LEE events in MINIboone
		_MB_evis_uz_corr.Scale(1. / 1000.);
		// _generated_evis_uz_corr->Scale(1. / _generated_evis_uz_corr->Integral());

		std::cout<<"LEERW:Initialize: _MB_evis_uz_corr integral is "<<_MB_evis_uz_corr.Integral()<<std::endl;
		std::cout<<"LEERW:Initialize: _generated_evis_uz_corr integral is "<<_generated_evis_uz_corr.Integral()<<std::endl;
		return true;
	}

	double LEERW::get_sculpting_weight(const ::larlite::mctruth* mytruth) {

		//Check to make sure LEERW instance is fully initialized
		check_is_initialized();

		//Grab relevant stuff for reweighting from mctruth object
		EventInfo_t evt_info = extract_event_info(mytruth);

		// std::cout << "is there nue? " << is_there_a_neutrino << ", number of electrons? " << n_electrons << std::endl;
		if (!evt_info.is_there_a_neutrino || evt_info.n_electrons != 1) {
			print(::larlite::msg::kWARNING, __FUNCTION__, "LEERW Package was handed an event with either no nue or no electron in the final state!");
			print_evt_info(evt_info);
			return 0;
		}

		// if (_debug) print_evt_info(evt_info);

		return get_sculpting_weight(evt_info.electron_energy_MEV, evt_info.electron_uz);

	}

	double LEERW::get_sculpting_weight(double electron_energy_MEV, double electron_uz) {

		//Check to make sure LEERW instance is fully initialized
		check_is_initialized();

		//Poll the 2D histogram to get an initial weight to sculpt energy and angle

		double weight_numerator = _MB_evis_uz_corr.GetBinContent(
		                              _MB_evis_uz_corr.GetXaxis()->FindBin(electron_energy_MEV),
		                              _MB_evis_uz_corr.GetYaxis()->FindBin(electron_uz)
		                          );
		double weight_denominator = _generated_evis_uz_corr.GetBinContent(
		                                _generated_evis_uz_corr.GetXaxis()->FindBin(electron_energy_MEV),
		                                _generated_evis_uz_corr.GetYaxis()->FindBin(electron_uz)
		                            );

		if (!weight_denominator) {
			if (_debug)
				std::cout << "WARNING: Can't find any entry in generated 2d hist with electron energy = "
				          << electron_energy_MEV
				          << " and uz = "
				          << electron_uz
				          << ". Returning 0 weight!" << std::endl;
			return 0.;
		}

		double weight = weight_numerator / weight_denominator;

		if (_debug) {
			std::cout << "(non-normalized) Sculpting (evis = " << electron_energy_MEV
			          << ", uz = " << electron_uz << ") weight is " << weight << "." << std::endl;
			std::cout << "This came from a numerator (miniboone) of " << weight_numerator
			          << " and a denominator (generated) of " << weight_denominator << std::endl;
		}

		return weight;
	}

	double LEERW::get_normalized_weight(const ::larlite::mctruth* mytruth) {

		//Check to make sure LEERW instance is fully initialized
		check_is_initialized();

		//Grab relevant stuff for reweighting from mctruth object
		EventInfo_t evt_info = extract_event_info(mytruth);

		// std::cout << "is there nue? " << is_there_a_neutrino << ", number of electrons? " << n_electrons << std::endl;
		if (!evt_info.is_there_a_neutrino || evt_info.n_electrons != 1) {
			print(::larlite::msg::kWARNING, __FUNCTION__, "LEERW Package was handed an event with either no nue or no electron in the final state!");
			return 0;
		}

		// if (_debug) print_evt_info(evt_info);

		return get_normalized_weight(evt_info.nue_energy_GEV);


	}

	double LEERW::get_normalized_weight(double nue_energy_GEV) {

		//Check to make sure LEERW instance is fully initialized
		check_is_initialized();

		double weight = 1.;

		//Poll the 1D histograms and use POT scaling/etc to get a scaling weight
		// POT weight is just the ratio of micro to miniboone POT
		weight *= _pot_weight;
		if (_debug)
			std::cout << "POT ratio weight is " << _pot_weight << "." << std::endl;

		// Tonnage weight is the ratio of micro to miniboone tonnage
		weight *= _tonnage_weight;
		if (_debug)
			std::cout << "Tonnage ratio weight is " << _tonnage_weight << "." << std::endl;

		//XSec weight uses neutrino energy in GEV. It also takes into account the molecular density of different materials.
		double dummy = _xsec_ratio.Eval(nue_energy_GEV);
		if (_debug)
			std::cout << "XSec ratio weight is " << dummy << "." << std::endl;
		weight *= dummy;

		//Flux weight uses neutrino energy in GEV. It comes from total neutrino flux ratio microboone to miniboone.
		dummy = _flux_ratio.Eval(nue_energy_GEV);
		if (_debug)
			std::cout << "Flux ratio weight is " << dummy << "." << std::endl;
		weight *= dummy;

		// Efficiency is not included in the weight calculation, this will come from whatever analysis is using this reweighter.

		//Overall normalization comes from the fact MiniBooNE saw 1212 excess events (MB efficiency unfolded)
		//Use the number of generated events to compute an overall normalization weight
		//weight *= _true_MB_excess_evts / _n_generated_evts;
		double temptemp = _true_MB_excess_evts / _n_generated_evts;
		// double temptemp = (_true_MB_excess_evts / _MB_evis_uz_corr.Integral() ) / ( _n_generated_evts / _generated_evis_uz_corr.Integral());
		

		// In the case of mcc7 files + mcc7 histo + mcc7 n events generated
		// If I don't multiply temptemp in, the reweighted neutrino spectrum integral is equal to
		// the number of generated events.
		// In the case of mcc6 files + mcc6 histo + mcc6 n events generated,
		// If I don't multiply temptemp in, the reweighted neutrino spectrum integral is equal to 442 ??
		weight *= temptemp;
		if (_debug)
			std::cout << "Overall normalization weight is " << temptemp << "." << std::endl;

		return weight;
	}

	const EventInfo_t LEERW::extract_event_info(const ::larlite::mctruth* mytruth) {

		EventInfo_t my_event_info;
		my_event_info.electron_energy_MEV = -1.;
		my_event_info.electron_uz = -2.;
		my_event_info.nue_energy_GEV = -1.;
		my_event_info.is_there_a_neutrino = false;
		my_event_info.n_electrons = 0;

		//First inspect the MCTruth object make sure it's ok (nue interaction, one electron in final state with status 1)
		if (!mytruth)
			throw std::invalid_argument("LEERW was handed a nonexistant mctruth!");

		auto &particles = mytruth->GetParticles();
		//Loop over the particles in the mctruth, make sure there is one nue, one electron.
		//Save the electron's energy [note: Should I be using deposited energy maybe?]
		for (auto const& particle : particles) {

			size_t PDG = abs(particle.PdgCode());
			if (PDG == 14)
				throw std::runtime_error("LEERW Package was handed a numu interaction?!");
			if (PDG == 12) {
				my_event_info.is_there_a_neutrino = true;
				my_event_info.nue_energy_GEV = particle.Trajectory().at(0).E();
			}
			//Note, when counting final state particles, ignore any particles that don't have status code == 1
			if (particle.StatusCode() != 1) continue;
			if (PDG == 11) {
				my_event_info.n_electrons++;
				my_event_info.electron_energy_MEV = particle.Trajectory().at(0).E() * 1000.;
				my_event_info.electron_uz = particle.Trajectory().at(0).Momentum().CosTheta();
			}
		}//end loop over mctruth particles

		return my_event_info;
	}

	void LEERW::check_is_initialized() {

		//Check to make sure LEERW instance is fully initialized
		if (!_flux_ratio.GetN() || !_xsec_ratio.GetN() || !_MB_evis_uz_corr.GetEntries() || !_generated_evis_uz_corr.GetEntries())
			throw std::runtime_error("LEERW Package not fully initialized: Missing input graphs/histos!");

		if (!_n_generated_evts)
			throw std::runtime_error("LEERW Package not fully initialized: Missing number of generated nue events!");
	}

	void LEERW::print_evt_info(const EventInfo_t evt_info) {

		std::cout << " -- DUMPING EVENT INFO -- " << std::endl;
		std::cout << "Is there a nue? " << evt_info.is_there_a_neutrino << std::endl;
		std::cout << "Number of electrons? " << evt_info.n_electrons << std::endl;
		std::cout << "Neutrino info: Energy = " << evt_info.nue_energy_GEV * 1000. << std::endl;
		std::cout << "Electron info: Energy = " << evt_info.electron_energy_MEV << ", Uz = " << evt_info.electron_uz << std::endl;
		std::cout << " -- END EVENT INFO DUMP -- " << std::endl;
	}
}
#endif

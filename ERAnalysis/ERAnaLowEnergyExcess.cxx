#ifndef ERTOOL_ERANALOWENERGYEXCESS_CXX
#define ERTOOL_ERANALOWENERGYEXCESS_CXX

#include "ERAnaLowEnergyExcess.h"

namespace ertool {

	ERAnaLowEnergyExcess::ERAnaLowEnergyExcess(const std::string& name)
		: AnaBase(name)
		, _result_tree(nullptr)
	{
		PrepareTreeVariables();
	}

	void ERAnaLowEnergyExcess::ProcessBegin() {

		/// Initialize the LEE reweighting package, if in LEE sample mode...
		if (_LEESample_mode) {
			_rw.set_debug(false);
			_rw.set_source_filename("$LARLITE_USERDEVDIR/LowEnergyExcess/LEEReweight/source/LEE_Reweight_plots.root");
			_rw.set_n_generated_events(6637);
			/// The input neutrinos were generated only in the TPC, not the entire cryostat
			_rw.set_events_generated_only_in_TPC(true);
			_rw.initialize();
		}

		// Build Box for TPC active volume
		_vactive  = ::geoalgo::AABox(0,
		                             -larutil::Geometry::GetME()->DetHalfHeight(),
		                             0,
		                             2 * larutil::Geometry::GetME()->DetHalfWidth(),
		                             larutil::Geometry::GetME()->DetHalfHeight(),
		                             larutil::Geometry::GetME()->DetLength());

	}


	bool ERAnaLowEnergyExcess::Analyze(const EventData &data, const ParticleGraph &graph)
	{
		_result_tree->SetName(Form("%s", _treename.c_str()));

		// Reset tree variables
		ResetTreeVariables();

		/// This variable seems to indicate whether a neutrino was reconstructed
		/// (IE a "ccsingleE" was found)
		bool reco = false;

		// size of ParticleSet should be the number of neutrinos found, each associated with a single electron
		auto const& particles = graph.GetParticleArray();

		// First off, if no nue was reconstructed, skip this event entirely.
		// Also if more than one nues were reconstructed, keep track of that too
		for ( auto const & p : particles )
			if ( abs(p.PdgCode()) == 12 ) {
				_n_nues_in_evt++;
				reco = true;
			}
		if (!reco) {
			// std::cout<<"No reconstructed nue in this event."<<std::endl;
			return false;
		}

		// Store the # of ertool showers in the event
		_n_ertool_showers = graph.GetParticleNodes(RecoType_t::kShower).size();

		// Reset the ertool::Shower copy of the ccsingleE-identified ertool::Shower
		singleE_shower.Reset();

		// Loop over particles and find the nue
		for ( auto const & p : particles ) {

			if ( abs(p.PdgCode()) == 12 ) {

				if (p.ProcessType() == kPiZeroMID) _maybe_pi0_MID = true;

				// Get the event timing from the most ancestor particle
				try {
					// Careful: p.Ancestor() returns a NODEID, but the data.Flash() function wants either a flash ID
					// or an actual particle. Instead, use data.Flash(graph.GetParticle(p.Ancestor()))
					_flash_time = data.Flash(graph.GetParticle(p.Ancestor()))._t;
					_summed_flash_PE = data.Flash(graph.GetParticle(p.Ancestor())).TotalPE();
				}
				catch ( ERException &e ) {}

				// Save the neutrino vertex to the ana tree
				_x_vtx = p.Vertex().at(0);
				_y_vtx = p.Vertex().at(1);
				_z_vtx = p.Vertex().at(2);

				// Save the reconstructed neutrino direction and momentum information to the ana tree
				_nu_theta = p.Momentum().Theta();
				_nu_p = p.Momentum().Length();
				_nu_pt = _nu_p * std::sin(_nu_theta);

				/// There are various ways to compute the neutrino energy.
				/// This function fills all the different reconstructed nue energy variables in the ttree
				FillRecoNuEnergies(p, graph, data);

				// get all reconstructed descendants of the neutrino and fill some relevant variables
				// "descendants" mean immediate children, their children, their children... all the way down
				auto const descendants = graph.GetAllDescendantNodes(p.ID());
				_n_children = descendants.size();
				for ( auto const & desc : descendants) {
					auto const & part = graph.GetParticle(desc);
					if (part.PdgCode() == 22) std::cout << "WTF gamma is daughter of neutrino?" << std::endl;
					if (abs(part.PdgCode()) == 13) _has_muon_child = true;
				}// for all neutrino descendants


				// Loop over the neutrinos immediate children and fill some relevant variables
				for (auto const& d : p.Children()) {

					auto daught = graph.GetParticle(d);

					// This is the "ccsinglee" electron.
					if (daught.PdgCode() == 11) {

						// Make a copy of the shower that is the ccsingleE
						singleE_shower = data.Shower(daught.RecoID());

						// Some info about the shower to store in the analysis ttree
						_e_theta = singleE_shower.Dir().Theta();
						_e_phi = singleE_shower.Dir().Phi();
						_e_Edep = singleE_shower._energy;
						_is_simple = isInteractionSimple(daught, graph);
						_dedx = data.Shower(daught.RecoID())._dedx;

						/// Fills _dist_2wall_shr and _dist_2wall_vtx
						FillBITEVariables(singleE_shower, p);
					}

					/// Compute longest track length associated with the immediate neutrino intxn
					if (daught.HasRecoObject() == true) {
						if (daught.RecoType() == kTrack) {
							auto mytrack = data.Track(daught.RecoID());
							double current_tracklen = ( mytrack.back() - mytrack.front() ).Length();
							if (current_tracklen > _longestTrackLen) _longestTrackLen = current_tracklen;
						}
					}// if the particle has a reco object
				} // End loop over neutrino children
			}// if we found the neutrino
		}// End loop over particles



		//// Now we loop over the MC particle graph and extract some MC information
		//// Note, the MC particle graph is made from mctracks and mcshowers.
		//// If your reco particle graph is also made from mctracks and showers,
		//// you can find a one-to-one match between the reco particle graph nodes
		//// and the MC particlegraph nodes in order to do reco-to-MC comparisons

		// Get MC particle set
		auto const& mc_graph = MCParticleGraph();
		// Get the MC data
		auto const& mc_data = MCEventData();

		if (!mc_graph.GetParticleArray().size())
			std::cout << "WARNING: Size of mc particle graph is zero! Perhaps you forgot to include mctruth/mctrack/mcshower?" << std::endl;

		// Compute a reweight:
		// in the case of BNB files, this is flux reweighting
		// in case of LEE sample, this is the LEERW package to make scaled excess
		// (note this also fills the _ptype variable)
		_weight = GetWeight(mc_graph);


		for ( auto const & mc : mc_graph.GetParticleArray() ) {

			// Find the shower particle in the mcparticlegraph that matches the object CCSingleE identified
			// as the single electron (note, the mcparticlegraph object could be a gamma, for example)
			// To do this, grab the ertool Shower in mc_data associated with each mcparticlegraph
			// shower particle and compare
			// (note this works for perfect-reco, but a more sophisticated method is needed for reco-reco)
			if (mc.RecoType() == kShower) {

				ertool::Shower mc_ertoolshower = mc_data.Shower(mc.RecoID());
				// We match ertool showers from mc particle graph to ertool showers from reco particle graphs
				// By comparing the energy to double precision... can consider also comparing _dedx and _time as well
				if (mc_ertoolshower._energy == singleE_shower._energy) {
					if (mc_ertoolshower._energy < 0 || mc_ertoolshower._energy > 5000)
						std::cout << "wtf? shower energy in mcgraph is stupid. it's " << mc_ertoolshower._energy << std::endl;
					auto parent = mc_graph.GetParticle(mc.Parent());
					_parentPDG = parent.PdgCode();
					_mcPDG = mc.PdgCode();
				}
			} // end loop if the reco type was mcshower

		} // end loop over mc particle graph


		/// Actually fill the analysis tree
		_result_tree->Fill();

		return true;
	}

	void ERAnaLowEnergyExcess::ProcessEnd(TFile * fout)
	{
		if (fout) {
			fout->cd();
			_result_tree->Write();
		}

		return;

	}

	double ERAnaLowEnergyExcess::GetWeight(const ParticleGraph mc_graph) {

		double nu_E_GEV = 1.;
		double e_E_MEV = -1.;
		double e_uz = -2.;

		for ( auto const & mc : mc_graph.GetParticleArray() ) {

			if (abs(mc.PdgCode()) == 12)
				nu_E_GEV = mc.Energy() / 1000.;
			if (abs(mc.PdgCode()) == 11) {
				e_E_MEV = mc.Energy();
				e_uz = std::cos(mc.Momentum().Theta());
			}

			if (!_LEESample_mode) {
				/// This stuff takes the truth neutrino information and computes a flux RW
				/// Make sure the neutrino is its own ancestor (it wasn't from something decaying IN the event)
				if ( (abs(mc.PdgCode()) == 12 || abs(mc.PdgCode()) == 14) && mc.Ancestor() == mc.ID() ) {

					int ntype = 0;
					int ptype = 0;
					double E = mc.Energy() / 1e3;

					if (mc.PdgCode() == 12)       ntype = 1;
					else if (mc.PdgCode() == -12) ntype = 2;
					else if (mc.PdgCode() ==  14) ntype = 3;
					else if (mc.PdgCode() == -14) ntype = 4;

					if (mc.ProcessType() == ::ertool::kK0L) ptype = 3;
					else if (mc.ProcessType() == ::ertool::kKCharged) ptype = 4;
					else if (mc.ProcessType() == ::ertool::kMuDecay) ptype = 1;
					else if (mc.ProcessType() == ::ertool::kPionDecay) ptype = 2;

					if (mc.ProcessType() != ::ertool::kK0L &&
					        mc.ProcessType() != ::ertool::kKCharged &&
					        mc.ProcessType() != ::ertool::kMuDecay &&
					        mc.ProcessType() != ::ertool::kPionDecay) {

						std::cout << " PDG : " << mc.PdgCode() << " Process Type : " << mc.ProcessType() << " from " <<
						          ::ertool::kK0L <<  " or " <<
						          ::ertool::kKCharged << " or " <<
						          ::ertool::kMuDecay << " or " <<
						          ::ertool::kPionDecay << std::endl;
					}

					_ptype = ptype;
					return _fluxRW.get_weight(E, ntype, ptype);
				}
			} // end if you aren't using LEESample mode
		} // end loop over mc particle graph

		if (_LEESample_mode) {
			if (e_E_MEV < 0 || e_uz < -1 || nu_E_GEV < 0)
				std::cout << "wtf i don't understand" << std::endl;
			return _rw.get_sculpting_weight(e_E_MEV, e_uz) * _rw.get_normalized_weight(nu_E_GEV);
		}

		/// You get here if you are running on cosmics (no truth neutrino in the event)
		return 1;
	}

	void ERAnaLowEnergyExcess::PrepareTreeVariables() {

		if (_result_tree) { delete _result_tree; }

		_result_tree = new TTree(Form("%s", _treename.c_str()), "Result Tree");
		_result_tree->Branch("_e_nuReco", &_e_nuReco, "e_nuReco/D");
		_result_tree->Branch("_e_nuReco_better", &_e_nuReco_better, "e_nuReco_better/D");
		_result_tree->Branch("_e_dep", &_e_dep, "e_dep/D");
		_result_tree->Branch("_weight", &_weight, "weight/D");
		_result_tree->Branch("_ptype", &_ptype, "ptype/I");
		_result_tree->Branch("_parentPDG", &_parentPDG, "parent_PDG/I");
		_result_tree->Branch("_mcPDG", &_mcPDG, "mc_PDG/I");
		_result_tree->Branch("_mcGeneration", &_mcGeneration, "mc_Generation/I");
		_result_tree->Branch("_longestTrackLen", &_longestTrackLen, "longest_tracklen/D");
		_result_tree->Branch("_x_vtx", &_x_vtx, "x_vtx/D");
		_result_tree->Branch("_y_vtx", &_y_vtx, "y_vtx/D");
		_result_tree->Branch("_z_vtx", &_z_vtx, "z_vtx/D");
		_result_tree->Branch("_closest_perpendicular_dist2wall_shr", &_closest_perpendicular_dist2wall_shr, "closest_perpendicular_dist2wall_shr/D");
		_result_tree->Branch("_closest_perpendicular_dist2wall_vtx", &_closest_perpendicular_dist2wall_vtx, "closest_perpendicular_dist2wall_vtx/D");
		_result_tree->Branch("_e_theta", &_e_theta, "_e_theta/D");
		_result_tree->Branch("_e_phi", &_e_phi, "_e_phi/D");
		_result_tree->Branch("_e_Edep", &_e_Edep, "_e_Edep/D");
		_result_tree->Branch("_e_CCQE", &_e_CCQE, "_e_CCQE/D");
		_result_tree->Branch("_nu_theta", &_nu_theta, "_nu_theta/D");
		_result_tree->Branch("_nu_pt", &_nu_pt, "_nu_pt/D");
		_result_tree->Branch("_nu_p", &_nu_p, "_nu_p/D");
		_result_tree->Branch("_n_children", &_n_children, "_n_children/I");
		_result_tree->Branch("_is_simple", &_is_simple, "_is_simple/O");
		_result_tree->Branch("_dedx", &_dedx, "dedx/D");
		_result_tree->Branch("_flash_time", &_flash_time, "flash_time/D");
		_result_tree->Branch("_summed_flash_PE", &_summed_flash_PE, "summed_flash_PE/D");
		_result_tree->Branch("_dist_2wall_shr", &_dist_2wall_shr, "dist_2wall_shr/D");
		_result_tree->Branch("_dist_2wall_vtx", &_dist_2wall_vtx, "dist_2wall_vtx/D");
		_result_tree->Branch("_maybe_pi0_MID", &_maybe_pi0_MID, "_maybe_pi0_MID/O");
		_result_tree->Branch("_n_ertool_showers", &_n_ertool_showers, "_n_ertool_showers/I");
		_result_tree->Branch("_n_nues_in_evt", &_n_nues_in_evt, "n_nues_in_evt/I");
		_result_tree->Branch("_has_muon_child", &_has_muon_child, "_has_muon_child/O");

		return;
	}

	void ERAnaLowEnergyExcess::ResetTreeVariables() {

		_e_nuReco = 0.;
		_e_nuReco_better = 0.;
		_e_dep = 0;
		_parentPDG = -99999;
		_ptype = -1;
		_mcPDG = -99999;
		_mcGeneration = -99999;
		_longestTrackLen = 0.;
		_x_vtx = -999.;
		_y_vtx = -999.;
		_z_vtx = -999.;
		_e_theta = -999.;
		_e_phi = -999.;
		_e_Edep = -999.;
		_e_CCQE = -999.;
		_nu_p = -999.;
		_nu_pt = -999.;
		_nu_theta = -999.;
		_n_children = -999;
		_is_simple = false;
		_dedx = -999.;
		_flash_time = -999999999.;
		_summed_flash_PE = -999999999.;
		_dist_2wall_vtx = -999.;
		_dist_2wall_shr = -999.;
		_maybe_pi0_MID = false;
		_n_ertool_showers = -1;
		_n_nues_in_evt = 0;
		_has_muon_child = false;

		return;

	}

	double ERAnaLowEnergyExcess::EnuCaloMissingPt(const std::vector< ::ertool::NodeID_t >& Children, const ParticleGraph &graph) {

		double Enu = 0;          //MeV
		double Elep = 0;         //MeV
		double Ehad = 0;         //MeV
		double pT = 0;           //MeV
		double Es = 30.5;        //MeV
		double mAr = 37211.3261; //MeV
		double mp =  938.28;     //MeV
		double mn = 939.57;      //MeV
		double Emdefect = 8.5;   //MeV //why 8.5? Because en.wikipedia.org/wiki/Nuclear_binding_energy, find something better.
		int nP = 0, nN = 0;
		auto XY = ::geoalgo::Vector(1, 1, 0);


		for (auto const& d : Children) {

			auto daught = graph.GetParticle(d);

			if (daught.PdgCode() == 11 || daught.PdgCode()) {
				Elep += daught.KineticEnergy();
			}
			else if (daught.RecoType() == kTrack || daught.RecoType() == kShower) {
				Ehad += daught.KineticEnergy();
			}

			pT += daught.Momentum().Dot(XY);

			if (daught.PdgCode() == 2212) { nP++;}
			if (daught.PdgCode() == 2112) { nN++;}
		}

		mAr -= nP * mp + nN * mn + (nN + nP) * Emdefect;

		// Enu = Elep + Ehad + Es +
		//       sqrt(pow(pT, 2) + pow(mAr, 2)) - mAr;
		Enu = Elep + Ehad + Es +
		      pow(pT, 2) / (2 * mAr);

		return Enu;
	}

	bool ERAnaLowEnergyExcess::isInteractionSimple(const Particle &singleE, const ParticleGraph &ps) {

		auto const &kids = ps.GetAllDescendantNodes(singleE.ID());
		auto const &bros = ps.GetSiblingNodes(singleE.ID());

		// Number of particles associated with this electron that are not protons, or the single e itself
		size_t _n_else = 0;
		for ( auto const& kid : kids )
			if (ps.GetParticle(kid).PdgCode() != 2212) _n_else++;
		for ( auto const& bro : bros )
			if (ps.GetParticle(bro).PdgCode() != 2212) _n_else++;

		return _n_else ? false : true;

	}

	void ERAnaLowEnergyExcess::FillBITEVariables(const Shower &singleE_shower, const Particle &p) {

	  ///###### B.I.T.E Analysis Start #####
	  // Build backward halflines
	  ::geoalgo::Vector inverse_shr_dir(-singleE_shower.Dir()[0],
					    -singleE_shower.Dir()[1],
					    -singleE_shower.Dir()[2]);
	  ::geoalgo::Vector inverse_vtx_dir(-p.Momentum().Dir()[0],
					    -p.Momentum().Dir()[1],
					    -p.Momentum().Dir()[2]);
	  ::geoalgo::HalfLine ext9_shr(singleE_shower.Start(), inverse_shr_dir);
	  ::geoalgo::HalfLine ext9_vtx(p.Vertex(), inverse_vtx_dir);
	  
	  
	  auto crs_tpc_ext9_shr  = _geoalg.Intersection(ext9_shr, _vactive);
	  auto crs_tpc_ext9_vtx  = _geoalg.Intersection(ext9_vtx, _vactive);
	  
	  double dist9_shr = 999.;
	  double dist9_vtx = 999.;
	  if (crs_tpc_ext9_shr.size())     dist9_shr  = crs_tpc_ext9_shr[0].Dist(singleE_shower.Start());
	  if (crs_tpc_ext9_vtx.size())     dist9_vtx  = crs_tpc_ext9_vtx[0].Dist(p.Vertex());
	  
	  _dist_2wall_shr = dist9_shr;
	  _dist_2wall_vtx = dist9_vtx;

	  //Calculate cloest perpendiculat distance to TPC wall
	  double PD2W_tmp_x;
	  double PD2W_tmp_y;
	  double PD2W_tmp_z;

	  //Shower
	  if (2 * larutil::Geometry::GetME()->DetHalfWidth()-singleE_shower.Start()[0]>singleE_shower.Start()[0]-0.0) PD2W_tmp_x = singleE_shower.Start()[0]-0.0;
	  else PD2W_tmp_x = 2 * larutil::Geometry::GetME()->DetHalfWidth() - singleE_shower.Start()[0];

	  if (singleE_shower.Start()[1]+larutil::Geometry::GetME()->DetHalfHeight()>-singleE_shower.Start()[1]+larutil::Geometry::GetME()->DetHalfHeight()) {
	    PD2W_tmp_y = -singleE_shower.Start()[1]+larutil::Geometry::GetME()->DetHalfHeight();
	  }
	  else PD2W_tmp_y =singleE_shower.Start()[1]+larutil::Geometry::GetME()->DetHalfHeight();

	  if (larutil::Geometry::GetME()->DetLength() - singleE_shower.Start()[2]>singleE_shower.Start()[2]-0.0) PD2W_tmp_z = singleE_shower.Start()[2] - 0.0;
	  else PD2W_tmp_z = larutil::Geometry::GetME()->DetLength() - singleE_shower.Start()[2];
	  
	  if (PD2W_tmp_x<=PD2W_tmp_y && PD2W_tmp_x<=PD2W_tmp_z)  _closest_perpendicular_dist2wall_shr  = PD2W_tmp_x;
	  if (PD2W_tmp_y<=PD2W_tmp_x && PD2W_tmp_x<=PD2W_tmp_z)  _closest_perpendicular_dist2wall_shr  = PD2W_tmp_y;
	  if (PD2W_tmp_z<=PD2W_tmp_x && PD2W_tmp_x<=PD2W_tmp_y)  _closest_perpendicular_dist2wall_shr  = PD2W_tmp_z;

	  //Vertex
	  if (2 * larutil::Geometry::GetME()->DetHalfWidth() - p.Vertex()[0]>p.Vertex()[0]-0.0) PD2W_tmp_x = p.Vertex()[0]-0.0;
	  else PD2W_tmp_x = 2 * larutil::Geometry::GetME()->DetHalfWidth()-p.Vertex()[0];

	  if (p.Vertex()[1]+larutil::Geometry::GetME()->DetHalfHeight()>-p.Vertex()[1]+larutil::Geometry::GetME()->DetHalfHeight()) {
	    PD2W_tmp_y = -p.Vertex()[1]+larutil::Geometry::GetME()->DetHalfHeight();
	  }
	  else PD2W_tmp_y =p.Vertex()[1]+larutil::Geometry::GetME()->DetHalfHeight();

	  if (larutil::Geometry::GetME()->DetLength() - p.Vertex()[2]>p.Vertex()[2]-0.0) PD2W_tmp_z = p.Vertex()[2] - 0.0;
	  else PD2W_tmp_z = larutil::Geometry::GetME()->DetLength() - p.Vertex()[2];

	  if (PD2W_tmp_x<=PD2W_tmp_y && PD2W_tmp_x<=PD2W_tmp_z)  _closest_perpendicular_dist2wall_vtx  = PD2W_tmp_x;
	  if (PD2W_tmp_y<=PD2W_tmp_x && PD2W_tmp_x<=PD2W_tmp_z)  _closest_perpendicular_dist2wall_vtx  = PD2W_tmp_y;
	  if (PD2W_tmp_z<=PD2W_tmp_x && PD2W_tmp_x<=PD2W_tmp_y)  _closest_perpendicular_dist2wall_vtx  = PD2W_tmp_z;
	  ///###### B.I.T.E Analysis END #####
	  
	}// End FillBITEVariables
  

	void ERAnaLowEnergyExcess::FillRecoNuEnergies(const Particle &nue, const ParticleGraph &graph, const EventData &data) {

		// get all reconstructed descendants of the neutrino in order to calculate total energy deposited
		// "descendants" mean immediate children, their children, their children... all the way down
		_e_dep = 0;
		_e_nuReco_better = 0;

		auto const descendants = graph.GetAllDescendantNodes(nue.ID());
		_n_children = descendants.size();

		for ( auto const & desc : descendants) {

			auto const & part = graph.GetParticle(desc);

			///haven't yet figured out how to use kINVALID_INT or whatever
			/// _e_nuReco_better adds just KE of protons but total energy (w/ mass) of pions
			if (!part.Children().size() && part.PdgCode() != 2212 && part.PdgCode() < 999999) {
				_e_nuReco_better += part.Mass();
			}

			// get the reco object's dep. energy
			// is this particle associated with a reconstructed track or shower?
			if (part.HasRecoObject() == true) {
				// if shower
				if (part.RecoType() == kShower) {
					_e_dep += data.Shower(part.RecoID())._energy;
					_e_nuReco_better += data.Shower(part.RecoID())._energy;
				}
				if (part.RecoType() == kTrack) {
					_e_dep += data.Track(part.RecoID())._energy;
					_e_nuReco_better += data.Track(part.RecoID())._energy;
				}
			}// if the particle has a reco object
		}// for all neutrino descendants

		// Compute "_e_nuReco" which is the neutrino energy from just the immediate children
		_e_nuReco = 0.;

		// Loop over the neutrinos immediate children
		for (auto const& d : nue.Children()) {

			auto daught = graph.GetParticle(d);

			// This is the "ccsinglee" electron.
			if (daught.PdgCode() == 11)
				_e_CCQE = _eccqecalc.ComputeECCQE(data.Shower(daught.RecoID())) * 1000.;


			//Note sometimes particle.KineticEnergy() is infinite!
			//however Track._energy is fine, so we'll use that for _e_nuReco
			if (daught.HasRecoObject() == true) {
				// get the reco object's dep. energy
				if (daught.RecoType() == kTrack)
					_e_nuReco += data.Track(daught.RecoID())._energy;
				if (daught.RecoType() == kShower)
					_e_nuReco += data.Shower(daught.RecoID())._energy;

			}// if the particle has a reco object
		} // End loop over neutrino children
	}// End FillRecoNuEnergies

} // end namespace ertool

#endif

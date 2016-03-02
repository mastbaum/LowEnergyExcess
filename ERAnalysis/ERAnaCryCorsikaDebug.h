/**
 * \file ERAnaCryCorsikaDebug.h
 *
 * \ingroup ERAnalysis
 * 
 * \brief Class def header for a class ERAnaCryCorsikaDebug
 *
 * @author davidkaleko
 */

/** \addtogroup ERAnalysis

    @{*/

#ifndef ERTOOL_ERANACRYCORSIKADEBUG_H
#define ERTOOL_ERANACRYCORSIKADEBUG_H

#include "ERTool/Base/AnaBase.h"
 #include "GeoAlgo/GeoAABox.h"
 #include "LArUtil/Geometry.h"
 #include "GeoAlgo/GeoSphere.h"

namespace ertool {

  /**
     \class ERAnaCryCorsikaDebug
     User custom Analysis class made by kaleko
   */
  class ERAnaCryCorsikaDebug : public AnaBase {
  
  public:

    /// Default constructor
    ERAnaCryCorsikaDebug(const std::string& name="ERAnaCryCorsikaDebug");

    /// Default destructor
    virtual ~ERAnaCryCorsikaDebug(){}

    /// Reset function
    virtual void Reset();

    /// Function to accept fclite::PSet
    void AcceptPSet(const ::fcllite::PSet& cfg);

    /// Called @ before processing the first event sample
    void ProcessBegin();

    /// Function to evaluate input showers and determine a score
    bool Analyze(const EventData &data, const ParticleGraph &ps);

    /// Called after processing the last event sample
    void ProcessEnd(TFile* fout=nullptr);

private:

  bool IsVertexActivity(const Particle &shower, const ParticleGraph &ps, const EventData &data);
  
  ::geoalgo::AABox _vactive;
  int n_with_parent_track;
  int n_with_parent_shower;
  int n_with_parent_uncontained_track;
  int n_with_parent_contained_track;
  int n_showers_with_start_not_contained;
  int n_showers_with_no_vtx_activity;
  int n_showers_from_pizeros;

  int n_evts;

  ::geoalgo::Sphere _vtx_sphere;
  };
}
#endif

/** @} */ // end of doxygen group 

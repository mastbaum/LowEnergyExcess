/**
 * \file ERSelSaveSingleEEvents.h
 *
 * \ingroup LowEnergyExcess
 * 
 * \brief Class def header for a class ERSelSaveSingleEEvents
 *
 * @author kaleko
 */

/** \addtogroup LowEnergyExcess

    @{*/

#ifndef LARLITE_ERSELSAVESINGLEEEVENTS_H
#define LARLITE_ERSELSAVESINGLEEEVENTS_H


#include "ERToolBackend/ERToolAnaBase.h"
namespace larlite {
  /**
     \class ERSelSaveSingleEEvents
     Example analysis unit for running singleE algos and saving events with reco nues in them
   */
  class ERSelSaveSingleEEvents : public ERToolAnaBase {
  
  public:

    /// Default constructor
    ERSelSaveSingleEEvents( const ::ertool::io::StreamType_t in_strm = ::ertool::io::kEmptyStream,
			const ::ertool::io::StreamType_t out_strm = ::ertool::io::kEmptyStream);

    /// Default destructor
    virtual ~ERSelSaveSingleEEvents(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 

/**
 * \file CosmicTriggerHacker.h
 *
 * \ingroup EventFilters
 * 
 * \brief Class def header for a class CosmicTriggerHacker
 *
 * @author davidkaleko
 */

/** \addtogroup EventFilters

    @{*/

#ifndef LARLITE_COSMICTRIGGERHACKER_H
#define LARLITE_COSMICTRIGGERHACKER_H

#include "Analysis/ana_base.h"
 #include "DataFormat/mctrack.h"

namespace larlite {
  /**
     \class CosmicTriggerHacker
     User custom analysis class made by SHELL_USER_NAME
   */
  class CosmicTriggerHacker : public ana_base{
  
  public:

    /// Default constructor
    CosmicTriggerHacker(){ _name="CosmicTriggerHacker"; _fout=0; _trigger_hack_tree=0;}

    /// Default destructor
    virtual ~CosmicTriggerHacker(){}

    /** IMPLEMENT in CosmicTriggerHacker.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CosmicTriggerHacker.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CosmicTriggerHacker.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:


    TTree* _trigger_hack_tree;
    double hacked_trig_time;
    
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

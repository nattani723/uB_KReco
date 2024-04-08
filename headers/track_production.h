#include "../CCKaonAnalyzerRebuild_module.h"


#include "art/Framework/Core/EDProducer.h"
//#include "art/Framework/Core/ModuleMacros.h"
//#include "art/Framework/Principal/Event.h"

//#include "fhiclcpp/ParameterSet.h"

//#include "lardataobj/RecoBase/Hit.h"
//#include "lardataobj/RecoBase/SpacePoint.h"
//#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"


//#include "art/Framework/Principal/Handle.h"
//#include "art/Framework/Principal/Run.h"
//#include "art/Framework/Principal/SubRun.h"

#include "art/Persistency/Common/PtrMaker.h"

//#include "canvas/Utilities/InputTag.h"

//#include "larcore/Geometry/Geometry.h"

//#include "lardata/Utilities/AssociationUtil.h"

//#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"

//#include "larpandora/LArPandoraInterface/Detectors/GetDetectorType.h"
//#include "larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h"
//#include "larpandora/LArPandoraEventBuilding/LArPandoraTrackCreation_module.cc"

#include <iostream>

using namespace pandora;
namespace Kaon_Analyzer
{
  /*
  void trackRebuid(std::vector<art::Ptr<recob::Hit>>& hit_list, 
				   art::FindManyP<recob::SpacePoint>& spacepoint_per_hit, 
				   const recob::Track& track);
  */

  double trackRebuid(std::vector<art::Ptr<recob::Hit>>& hit_list, 
				   art::FindManyP<recob::SpacePoint>& spacepoint_per_hit, 
				   const recob::Track& track);

  recob::Track trackRebuid2(std::vector<art::Ptr<recob::Hit>>& hit_list, 
				   art::FindManyP<recob::SpacePoint>& spacepoint_per_hit, 
				   const recob::Track& track);

  recob::Track buildTrack(int id,
			  lar_content::LArTrackStateVector& trackStateVector);

  //------------------------------------------------

  /*
  void trackRebuid(std::vector<art::Ptr<recob::Hit>>& hit_list,
		   art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
		   const recob::Track& track)
  */



  recob::Track trackRebuid2(std::vector<art::Ptr<recob::Hit>>& hit_list,
			   art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
			   const recob::Track& track)
  {
    
    //if(hit_list.size()!=0) continue;

    art::ServiceHandle<geo::Geometry> theGeometry;
    const unsigned int nWirePlanes(theGeometry->MaxPlanes());

    if (nWirePlanes > 3)
        throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- More than three wire planes present ";

    if ((0 == theGeometry->Ncryostats()) || (0 == theGeometry->NTPC(0)))
        throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- unable to access first tpc in first cryostat ";

    std::unordered_set<geo::_plane_proj> planeSet;
    for (unsigned int iPlane = 0; iPlane < nWirePlanes; ++iPlane)
        (void) planeSet.insert(theGeometry->TPC(0, 0).Plane(iPlane).View());

    if ((nWirePlanes != planeSet.size()) || !planeSet.count(geo::kU) || !planeSet.count(geo::kV) || (planeSet.count(geo::kW) && planeSet.count(geo::kY)))
        throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- expect to find u and v views; if there is one further view, it must be w or y ";

    const bool useYPlane((nWirePlanes > 2) && planeSet.count(geo::kY));  

    const float wirePitchU(theGeometry->WirePitch(geo::kU));
    const float wirePitchV(theGeometry->WirePitch(geo::kV));
    const float wirePitchW((nWirePlanes < 3) ? 0.5f * (wirePitchU + wirePitchV) : (useYPlane) ? theGeometry->WirePitch(geo::kY) : theGeometry->WirePitch(geo::kW));

    //int trackCounter=1000;
    int trackCounter(0);

    const int m_slidingFitHalfWindow = 20;

    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec; 
    pandora::CartesianPointVector pandora_hit_positions;
    
    for(size_t i_h=0; i_h<hit_list.size(); i_h++){
      
      spacepoint_vec.clear();
      spacepoint_vec = spacepoint_per_hit.at(hit_list[i_h].key()); 
      if(spacepoint_vec.size()!=1) continue;
      
      const TVector3 hit_position = spacepoint_vec[0]->XYZ();
      const pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
      pandora_hit_positions.push_back(pandora_hit_position);
      
    }
    
    std::unique_ptr<std::vector<recob::Track>> outputTracks(new std::vector<recob::Track>);
    //const float wirePitchW = 0.30;
    
    /*
    recob::tracking::SMatrixSym33 mtrx;// = track.SMatrixSym55();
    recob::Vertex vertex(track.End(), mtrx, 0, 0, 0);
    
    double vertexXYZ[3] = {0., 0., 0.};
    vertex.XYZ(vertexXYZ);
    const pandora::CartesianVector vertexPosition(vertexXYZ[0], vertexXYZ[1], vertexXYZ[2]);
    */
    const pandora::CartesianVector vertexPosition(track.End().X(), track.End().Y(), track.End().Z());
    
    
    lar_content::LArTrackStateVector trackStateVector;
    pandora::IntVector indexVector;
    
    lar_content::LArPfoHelper::GetSlidingFitTrajectory(pandora_hit_positions,
						       vertexPosition,
						       m_slidingFitHalfWindow,
						       wirePitchW,
						       trackStateVector,
						       &indexVector);

    //cout << "trackStateVector.size(): " << trackStateVector.size() << endl;
    //cout << "hit_list.size(): " << hit_list.size() << endl;
    
    //outputTracks->emplace_back(buildTrack(0, trackStateVector));
    outputTracks->emplace_back(buildTrack(trackCounter++, trackStateVector));

    //cout << "Rebuilt track length: " <<  outputTracks->at(0).Length() << endl;
    //double trkln = outputTracks->at(0).Length();

    recob::Track reco_track = outputTracks->at(0);
    return reco_track;

    /*
    art::Handle< std::vector<recob::Track> > trackListHandle_new;  
    std::vector<art::Ptr<recob::Track> > tracklist_new;
    if(evt.getByLabel(fTrackModuleLabel_new,trackListHandle_new)) { 
      art::fill_ptr_vector(tracklist_new, trackListHandle_new);
    }

    std::type_info const& typeInfoOfWrapper{typeid(art::Wrapper< std::vector<recob::Track> >)};
    auto res = getByLabel(typeInfoOfWrapper, inputTag);

    std::vector<recob::Track> tracklistvec = trackListHandle_new.product();
    tracklistvec.push_back(reco_track);
    trackListHandle_new = art::Handle< std::vector<recob::Track> >{tracklistvec, res.second};
    */

  }



  double trackRebuid(std::vector<art::Ptr<recob::Hit>>& hit_list,
		     art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
		     const recob::Track& track)
  {
    
    if(hit_list.size()==0) return -999;

    art::ServiceHandle<geo::Geometry> theGeometry;
    const unsigned int nWirePlanes(theGeometry->MaxPlanes());

    if (nWirePlanes > 3)
        throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- More than three wire planes present ";

    if ((0 == theGeometry->Ncryostats()) || (0 == theGeometry->NTPC(0)))
        throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- unable to access first tpc in first cryostat ";

    std::unordered_set<geo::_plane_proj> planeSet;
    for (unsigned int iPlane = 0; iPlane < nWirePlanes; ++iPlane)
        (void) planeSet.insert(theGeometry->TPC(0, 0).Plane(iPlane).View());

    if ((nWirePlanes != planeSet.size()) || !planeSet.count(geo::kU) || !planeSet.count(geo::kV) || (planeSet.count(geo::kW) && planeSet.count(geo::kY)))
        throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- expect to find u and v views; if there is one further view, it must be w or y ";

    const bool useYPlane((nWirePlanes > 2) && planeSet.count(geo::kY));  

    const float wirePitchU(theGeometry->WirePitch(geo::kU));
    const float wirePitchV(theGeometry->WirePitch(geo::kV));
    const float wirePitchW((nWirePlanes < 3) ? 0.5f * (wirePitchU + wirePitchV) : (useYPlane) ? theGeometry->WirePitch(geo::kY) : theGeometry->WirePitch(geo::kW));

    const int m_slidingFitHalfWindow = 20;

    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec; 
    pandora::CartesianPointVector pandora_hit_positions;
    
    for(size_t i_h=0; i_h<hit_list.size(); i_h++){
      
      spacepoint_vec.clear();
      spacepoint_vec = spacepoint_per_hit.at(hit_list[i_h].key()); 
      if(spacepoint_vec.size()!=1) continue;
      
      const TVector3 hit_position = spacepoint_vec[0]->XYZ();
      const pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
      pandora_hit_positions.push_back(pandora_hit_position);
      
    }
    
    std::unique_ptr<std::vector<recob::Track>> outputTracks(new std::vector<recob::Track>);
    //std::unique_ptr<recob::Track> outputTrack(new recob::Track);
    //recob::Track outputTrack(new recob::Track);
    //recob::Track outputTrack;

    //LArPandoraDetectorType* detType(detector_functions::GetDetectorType());
    //const float wirePitchW(detType->WirePitchW());
    //const float wirePitchW = 0.30;
    
    //int track_counter = 0;
    //const art::PtrMaker<recob::Track> makeTrackPtr(evt);
    
    //organise inputs
    
    //lar_pandora::PFParticleVector pfParticleVector, extraPfParticleVector;
    //lar_pandora::PFParticlesToSpacePoints pfParticlesToSpacePoints;
    //lar_pandora::PFParticlesToClusters pfParticlesToClusters;
    
    //collect pfparticle
    
    //Consider the end of K track candidate as pi+ 
    //recob::Vertex vertex = track.End();

    //Point_t pos = track.End();
    //recob::tracking::SMatrixSym33 = track.SMatrixSym55();
    //recob::Vertex vertex(pos, track.SMatrixSym55(), 0, 0, 0);
    //recob::Vertex vertex(track.End(), track.SMatrixSym55(), 0, 0, 0);

    recob::tracking::SMatrixSym33 mtrx;// = track.SMatrixSym55();
    recob::Vertex vertex(track.End(), mtrx, 0, 0, 0);
    //lar_pandora::LArPandoraOutput::BuildVertex(vertex, 0);
    
    double vertexXYZ[3] = {0., 0., 0.};
    vertex.XYZ(vertexXYZ);
    const pandora::CartesianVector vertexPosition(vertexXYZ[0], vertexXYZ[1], vertexXYZ[2]);
    
    
    lar_content::LArTrackStateVector trackStateVector;
    pandora::IntVector indexVector;
    
    lar_content::LArPfoHelper::GetSlidingFitTrajectory(pandora_hit_positions,
						       vertexPosition,
						       m_slidingFitHalfWindow,
						       wirePitchW,
						       trackStateVector,
						       &indexVector);
    
    //outputTracks->emplace_back(lar_pandora::LArPandoraTrackCreation().BuildTrack(0, trackStateVector));

    if(trackStateVector.empty()) return -999;

    outputTracks->emplace_back(buildTrack(0, trackStateVector));
    //outputTrack().lar_pandora::LArPandoraTrackCreation::BuildTrack(0, trackStateVector);

    //cout << "Rebuilt track length: " <<  outputTracks->at(0).Length() << endl;
    double trkln = outputTracks->at(0).Length();

    recob::Track reco_track = outputTracks->at(0);
    //art::Ptr<recob::Track>& reco_ptrack = *reco_track;
    //const recob::Track& track2 = *ptrack2;  
    //cout << "ptrack.key(): " << reco_track.key() << endl;

    return trkln;

  }




  recob::Track buildTrack(int id,
			  lar_content::LArTrackStateVector& trackStateVector)
    {
      
      if (trackStateVector.empty()) cout << "BuildTrack - No input trajectory points provided" << endl;

      recob::tracking::Positions_t xyz;
      recob::tracking::Momenta_t pxpypz;
      recob::TrackTrajectory::Flags_t flags;

      for (const lar_content::LArTrackState& trackState : trackStateVector) {

	xyz.emplace_back(recob::tracking::Point_t(trackState.GetPosition().GetX(),
						  trackState.GetPosition().GetY(),
						  trackState.GetPosition().GetZ()));
	pxpypz.emplace_back(recob::tracking::Vector_t(trackState.GetDirection().GetX(),
						      trackState.GetDirection().GetY(),
						      trackState.GetDirection().GetZ()));

	// Set flag NoPoint if point has bogus coordinates, otherwise use clean flag set
	if (std::fabs(trackState.GetPosition().GetX() - util::kBogusF) <
            std::numeric_limits<float>::epsilon() &&
	    std::fabs(trackState.GetPosition().GetY() - util::kBogusF) <
            std::numeric_limits<float>::epsilon() &&
	    std::fabs(trackState.GetPosition().GetZ() - util::kBogusF) <
            std::numeric_limits<float>::epsilon()) {
	  flags.emplace_back(recob::TrajectoryPointFlags(recob::TrajectoryPointFlags::InvalidHitIndex,
							 recob::TrajectoryPointFlagTraits::NoPoint));
	}
	else {
	  flags.emplace_back(recob::TrajectoryPointFlags());
	}
      }

      // note from gc: eventually we should produce a TrackTrajectory, not a Track with empty covariance matrix and bogus chi2, etc.
      return recob::Track(
			  recob::TrackTrajectory(std::move(xyz), std::move(pxpypz), std::move(flags), false),
			  util::kBogusI,
			  util::kBogusF,
			  util::kBogusI,
			  recob::tracking::SMatrixSym55(),
			  recob::tracking::SMatrixSym55(),
			  id);

    }
}

#include "CCKaonProducer_module.h"
#include "headers/particle_split_basetool_producer.h"
#include "headers/track_production_producer.h"

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/EDProducer.h"
//#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Principal/View.h"
#include "art/Utilities/make_tool.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h" 

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h" 
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "/uboone/app/users/taniuchi/51_pandora/srcs/larpandoracontent/larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h" 

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/Deprecated/BezierTrack.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "larcoreobj/SummaryData/POTSummary.h"


#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"


#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "ubana/AnalysisTree/MCTruth/IMCTruthMatching.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "ubobj/Optical/UbooneOpticalFilter.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "PID/LLR_PID.h"
#include "PID/LLRPID_proton_muon_lookup.h"
#include "PID_K/LLR_PID_K.h"
#include "PID_K/LLRPID_kaon_proton_lookup.h"

#include "Pandora/PdgTable.h" 

#include "TCanvas.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphDelaunay.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TTimeStamp.h"

#include <iostream>
#include <TH2D.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TVectorT.h>
#include <THStack.h>
#include <TPDF.h>
#include <TLegend.h>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include <array>
#include <map>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <string>
#include <numeric>
#include <algorithm>
#include <functional>
#include <typeinfo>
#include <memory>
#include <chrono>
#include <sys/stat.h>
#include "TDirectory.h"

//#ifdef __MAKECINT__
#ifdef __CLING__
#pragma link C++ class std::vector < std::vector<Float_t> >+; 
#pragma link C++ class std::vector < std::vector< std::vector<Float_t> > >+; 
#endif


using namespace std;
namespace Kaon_Analyzer{

  CCKaonProducer::CCKaonProducer(fhicl::ParameterSet const& pset) :
    EDProducer(pset),
    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","gaushit")),
    fLArG4ModuleLabel         (pset.get< std::string >("LArG4ModuleLabel","largeant")),
    fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel","generator")),  
    fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","pandora")),
    fShowerModuleLabel        (pset.get< std::string >("ShowerModuleLabel","pandora")),
    fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","pandoracaliSCE")),
    fPIDLabel                 (pset.get< std::string >("PIDLabel","pandoracalipidSCE")),
    fHitTruthAssns            (pset.get< std::string >("HitTruthAssn","gaushitTruthMatch")), 
    fHitTrackAssns            (pset.get< std::string >("HitTrackAssn","pandora")), 
    fHitShowerAssns            (pset.get< std::string >("HitShowerAssn","pandora")), 
    m_pfp_producer            (pset.get< std::string >("pfp_producer","pandora")),
    fPFParticleLabel          (pset.get< std::string >("PFParticleLabel", "pandora")),
    fSpacePointproducer       (pset.get< std::string >("SpacePointproducer", "pandora")),
    fPandoraLabel             (pset.get< std::string >("PandoraLabel", "pandora")),
    isMC                      (pset.get< bool >("IsMC",true))
    
  {
    
    theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
    SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    geom = lar::providerFrom<geo::Geometry>();
    
    produces< std::vector<recob::Track> >(); 
    produces< art::Assns<recob::Track, recob::Hit> >();

  }
  
  CCKaonProducer::~CCKaonProducer()
  {
    //destructor
  }


  void CCKaonProducer::beginJob()
  {
    //art::ServiceHandle<art::TFileService> tfs;
    //fEventTree->Branch();
  }


  void CCKaonProducer::produce(art::Event& evt)
  {
    art::Handle< std::vector<simb::MCTruth> > mctruths;
    std::vector< art::Ptr<simb::MCParticle> > ptList;
   
    std::unique_ptr< std::vector<recob::Track> > anaTrackCollection(new std::vector<recob::Track> );
    std::unique_ptr<art::Assns<recob::Track, recob::Hit>>  anaTrackHitAssociations(new art::Assns<recob::Track, recob::Hit>);
   
    // get MCTruth
    evt.getByLabel("generator", mctruths);
    if (mctruths->size()!=1) {
      //return;
      //continue;
    }
    simb::MCTruth mctruth = mctruths->at(0);
    

    // get MCParticles
    art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
    if (evt.getByLabel(fLArG4ModuleLabel, mcParticleHandle)){
      art::fill_ptr_vector(ptList, mcParticleHandle); 
    }   

        
    // Collect all recontructed particles
    art::Handle<std::vector<recob::PFParticle>> pfparticles;
    std::vector< art::Ptr<recob::PFParticle> > pfpVector;
    if(evt.getByLabel(m_pfp_producer, pfparticles)){
      art::fill_ptr_vector(pfpVector, pfparticles);
    }
    
    //event has no particle
    if (pfparticles->size()==0) {
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }

    // Get PFParticle associations
    art::FindManyP<anab::T0> pfp_muon_assn(pfparticles, evt, "NuCCproducer");
    if(!pfp_muon_assn.isValid()){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }
    
    art::FindManyP<recob::Track> pfparticleTrackAssn(pfparticles, evt, "pandora");
    if(!pfparticleTrackAssn.isValid()){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }
    
    art::FindManyP<recob::Vertex> pfparticleVertexAssn(pfparticles, evt, "pandora");
    if(!pfparticleVertexAssn.isValid()){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }

    

    // Find recontructed neutrino (there should be one)
    lar_pandora::PFParticleVector pfneutrinos(0);
    for (unsigned int i=0; i<pfparticles->size(); ++i) {
      
      art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);      
      if (pfparticle->IsPrimary() && pfparticle->PdgCode()==14) {
	pfneutrinos.push_back(pfparticle);
      }
      
    }
    
    if (pfneutrinos.size() != 1) {
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }
    
    art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();

    
    // Find CC muon and daughters
    lar_pandora::PFParticleVector pfmuons(0);
    vector<int> reco_nu_daughters_id(0);

    for (unsigned int i=0; i<pfparticles->size(); ++i) {
      
      art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);
      
      // look at particles with neutrino parent and one associated track
      if (pfparticle->Parent()==pfnu->Self() && pfparticleTrackAssn.at(i).size()==1) {
	
	art::Ptr<recob::Track> track = pfparticleTrackAssn.at(i).front();
	reco_nu_daughters_id.push_back(track.key());
	
	// CC muon has a T0 associated
	if (pfp_muon_assn.at(i).size()==1) {
	  pfmuons.push_back(pfparticle);
	}

      }      
    }
    
    reco_nu_ndaughters = reco_nu_daughters_id.size();
    
    if (pfmuons.size()!=1) {
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }
    
    art::Ptr<recob::PFParticle> pfmuon = pfmuons.front();
    art::Ptr<recob::Track> trkmuon = pfparticleTrackAssn.at(pfmuon.key()).front();
    

    // get hits
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if(evt.getByLabel(fHitsModuleLabel,hitListHandle)){
      art::fill_ptr_vector(hitlist, hitListHandle);
    }

    // get tracks
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) {
      art::fill_ptr_vector(tracklist, trackListHandle);
    }

    // get showers
    art::Handle< std::vector<recob::Shower> > showerListHandle;
    std::vector<art::Ptr<recob::Shower> > showerlist;
    if(evt.getByLabel(fShowerModuleLabel,showerListHandle)) {
      art::fill_ptr_vector(showerlist, showerListHandle);
    }

    /*
    art::Handle< std::vector<recob::Cluster> > clusterHandle;
    std::vector< art::Ptr<recob::Cluster> > clusterVector;
    if(evt.getByLabel(fPandoraLabel,clusterHandle)){
      art::fill_ptr_vector(clusterVector,clusterHandle);
    }
    */

    art::Handle< std::vector<recob::SpacePoint> > spacepointHandle;
    std::vector< art::Ptr<recob::SpacePoint> > spacepointVector;
    if(evt.getByLabel(fPandoraLabel,spacepointHandle)){
      art::fill_ptr_vector(spacepointVector,spacepointHandle);
    }


    art::FindManyP<recob::Hit> hits_from_tracks(trackListHandle, evt, fHitTrackAssns);
    art::FindManyP<recob::Hit> hits_from_showers(showerListHandle, evt, fHitShowerAssns);
    if(!hits_from_tracks.isValid()){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }


    //art::FindManyP<recob::Cluster> clusters_per_pfparticle(pfparticles, evt, fPandoraLabel);
    //art::FindManyP<recob::Hit> hits_per_cluster(clusterHandle, evt, fPandoraLabel);
    //std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > pfParticleToClustersMap;
    //std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> > clusterToHitsMap;

    ///

    art::FindOneP<recob::Hit> findSPToHit(spacepointVector, evt, fSpacePointproducer); 
    art::FindManyP<recob::Hit> findTrackToHit(trackVector, evt, fTrackModuleLabel);
    art::FindManyP<recob::Hit> findShowerToHit(showerVector, evt, fShowerModuleLabel);
    //art::FindManyP<recob::Hit> findSPToHit(spacepointHandle, evt, fSpacePointproducer);

    std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>> hitToSpacePointMap;
    std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>> spacepointToHitMap;
    //std::map<art::Ptr<recob::SpacePoint>,  std::vector<art::Ptr<recob::Hit>> > spacepointToHitsMap;


    art::FindManyP<recob::SpacePoint> findPFParticlesToSPs(pfparticles, evt, fSpacePointproducer);
    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>> > pfParticleToSpacePointsMap;
    //art::FindManyP<recob::SpacePoint> spacepoints_per_pfparticle(pfparticles, evt, fSpacePointproducer);
    //until here




    art::FindManyP<recob::SpacePoint> spacepoints_per_pfparticle(pfparticles, evt, fSpacePointproducer);
    //art::FindManyP<recob::Hit> hits_per_spacepoint(spacepointHandle, evt, fSpacePointproducer);
    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>> > pfParticleToSpacePointsMap;
    //std::map<art::Ptr<recob::SpacePoint>,  std::vector<art::Ptr<recob::Hit>> > spacepointToHitsMap;


    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMap;
    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsFromSpacePointsMap;
    //std::map<art::Ptr<recob::PFParticle>, int> pfParticleToNHitsFromSpacePoints;

    if( !hits_per_cluster.isValid() || !spacepoints_per_pfparticle.isValid() || !hits_per_spacepoint.isValid()  ){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }


    const art::PtrMaker<recob::Track> makeTrackPtr(evt);

    for (unsigned int ipfp=0; ipfp < pfpVector.size(); ++ipfp) {
      const art::Ptr<recob::PFParticle> pfp = pfpVector.at(ipfp);
      pfParticleToSpacePointsMap[pfp] = spacepoints_per_pfparticle.at(pfp.key());
    }


    for (unsigned int iSP = 0; iSP < spacepointVector.size(); ++iSP) { 
      const art::Ptr<recob::SpacePoint> spacepoint = spacepointVector.at(iSP);
      const art::Ptr<recob::Hit> hit = findSPToHit.at(iSP);
      spacepointToHitMap[spacepoint] = hit;
      hitToSpacePointMap[hit] = spacepoint;
      
      //const std::vector<art::Ptr<recob::Hit>> hits = findSPToHit[spacepoint.key()];
      //spacepointToHitsMap[spacepoint] = hits;
      /*
      for(unsigned int ihit=0; ihit < hits.size(); ++ihit){
	hitToSpacePointMap[hits.at(ihit)] = spacepoint;
      }
      */
    }

    for (unsigned int ipfp=0; ipfp < pfpVector.size(); ++ipfp) {

      auto pfp = pfpVector[i];
      std::vector<art::Ptr<recob::SpacePoint>> spacepoints_vec  = pfParticleToSpacePointsMap[pfp];
      std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};

      //lar_pandora::HitSet hitsInParticleSet;

      for (art::Ptr<recob::SpacePoint> spacepoint: spacepoints_vec){
	std::vector<art::Ptr<recob::Hit>> hits_vec =  spacepointToHitsMap[spacepoint];
	hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );

	/*
	for(auto hit : hits_vec){
	  (void) hitsInParticleSet.insert(hit);
	}
	*/
	
      }

      //pfParticleToNHitsFromSpacePoints[pfp] = hits_for_pfp.size();
      pfParticleToHitsFromSpacePointsMap[pfp] = hits_for_pfp;

    }


    //old
    for(size_t i=0; i< spacepointVector.size(); ++i){
      auto spacepoint = spacepointVector[i];
      spacepointToHitsMap[spacepoint] = hits_per_spacepoint.at(spacepoint.key());

      //auto const& spkey = spacepointVector.at(i).key();
    }

    for (unsigned int i=0; i<pfpVector.size(); ++i) {

      auto pfp = pfpVector[i];
      std::vector<art::Ptr<recob::Cluster>> clusters_vec  = pfParticleToClustersMap[pfp];
      std::vector<art::Ptr<recob::SpacePoint>> spacepoints_vec  = pfParticleToSpacePointsMap[pfp];
      std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};

      lar_pandora::HitSet hitsInParticleSet;

      for (art::Ptr<recob::SpacePoint> spacepoint: spacepoints_vec){
	std::vector<art::Ptr<recob::Hit>> hits_vec =  spacepointToHitsMap[spacepoint];
	hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );

	for(auto hit : hits_vec){
	  (void) hitsInParticleSet.insert(hit);
	}
	
      }
      pfParticleToNHitsFromSpacePoints[pfp] = hits_for_pfp.size();
      pfParticleToHitsFromSpacePointsMap[pfp] = hits_for_pfp;

      for (art::Ptr<recob::Cluster> cluster: clusters_vec){
	std::vector<art::Ptr<recob::Hit>> hits_vec =  clusterToHitsMap[cluster];

	for(auto hit : hits_vec){
	  if(hitsInParticleSet.count(hit) == 0) hits_for_pfp.push_back(hit);
	}
      }
      cout << "hits_for_pfp size after cluster " << hits_for_pfp.size() <<endl;
      pfParticleToHitsMap[pfp] = hits_for_pfp;
    }
  

    
    // find track multiplicity around vertex
    //int NTracks=tracklist.size();
    //int NShowers=showerlist.size();
    //int ntracks = 0;
    
    //loop over primary nu track
    for (int i=0; i<reco_nu_ndaughters; i++) {

  //SPList spacepointVector;

      /*
    vector<bool> v_trk_flg_peak;
    std::map<double, TVector2, std::greater<>> view_peak_map;
    vector<TVector2> best_peak_bins;

    std::map<int, std::map<int, double>> angular_distribution_map_3D_track;
    std::map<int, std::map<int, double>> angular_distribution_map_3D_shower;
    std::map<int, std::map<int, double>> angular_distribution_map_3D_pfparticle;

    std::vector<art::Ptr<recob::Hit>> unavailable_hit_list;
    std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list;
    std::vector<std::vector<art::Ptr<recob::Hit>>> shower_spine_hit_list_vector;

    std::vector<art::Ptr<recob::Hit>> hits_from_reco;
      */

      //get list of sps
      //get list of unavailable list

    art::Ptr<recob::Track> ptrack(trackListHandle,reco_nu_daughters_id[i]);
    const recob::Track& track = *ptrack;

    /*
    lar_pandora::HitVector hitsInParticle;
    unsigned int hitsFromSpacePointsSize = 0;
    
    for (unsigned int i=0; i<pfparticles->size(); ++i) {  
      //for (unsigned int i=0; i<pfpVector.size(); ++i) {
      //auto pfp = pfpVector[i];
      art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);
      if (pfparticle->Parent()==pfnu->Self() && pfparticleTrackAssn.at(i).size()==1){ 
	art::Ptr<recob::Track> ppfptrack = pfparticleTrackAssn.at(i).front();
	//const recob::PFParticle& pfp = *pfparticle;
	//const recob::Track& pfptrack = *ppfptrack;
	
	if(ppfptrack.key() == ptrack.key()){
	  //hitsInParticle.clear();
	  hitsInParticle = pfParticleToHitsMap[pfparticle];
	  //hitsInParticle = pfParticleToHitsFromSpacePointsMap[pfparticle];
	  hitsFromSpacePointsSize = pfParticleToNHitsFromSpacePoints[pfparticle];
	  //util::CreateAssn(*this, evt, ptrack, pfparticle, *(anaParticleTrackAssociations.get()));
	  cout << "hitsInParticle.size()" << hitsInParticle.size() << endl;
	  break;
	}
      }
    }
    */



    /*
      ParticleDirectionFinder a;
      vector<TVector3> peak_direction_vector;
      a.Run(spacepointVector, ptrack, hits_from_track, peak_direction_vector);
      TrackHitCollector b;
      b.Run();
    */

    //art::FindManyP<recob::SpacePoint> spacepoint_per_hit_test(hitListHandle, evt, fSpacePointproducer);

    //recob::Track primary_track = trackRebuid(trackcount, hits_from_track, spacepoint_per_hit_test, track);

    //std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec_test;

    //for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
    //spacepoint_vec_test.clear();
    //spacepoint_vec_test = spacepoint_per_hit_test.at(hits_from_track[i_h].key());

    //}
    //cout << "hits_from_track.size(): " << hits_from_track.size() << endl;
    //cout << "hitsFromSpacePointsSize: " << hitsFromSpacePointsSize << endl;


    //push back the hit
    //size_t prev_size = anaHitCollection->size(); 
    //std::unique_ptr< std::vector<recob::Hit> > anaHitCollection_tmp( new std::vector<recob::Hit> );
    //lar_pandora::HitVector anaHitCollection_tmp;
    //anaHitCollection_tmp.clear();
    /*
    for(auto hitptr : hits_from_track){
      //anaHitCollection->push_back(*hitptr);
      //cout << "hit.key(): " << hitptr->key() << endl;
      anaHitCollection_tmp.push_back(hitptr);
    }
    */

    //art::Ptr<recob::Track> pTrack(makeTrackPtr(anaTrackCollection->size() - 1));

    //ATTN: metadata added with index from space points if available, null for others
    /*                                                                        
    for (unsigned int hitIndex = 0; hitIndex < hitsInParticle.size(); hitIndex++){
	const art::Ptr<recob::Hit> phit(hitsInParticle.at(hitIndex));
	//std::vector<art::Ptr<recob::SpacePoint>> hitsFromSpacePoints = spacepoint_per_hit_test.at(hits_from_track[hitIndex].key());
	const int index((hitIndex < hitsFromSpacePointsSize) ? hitIndex : std::numeric_limits<int>::max());
	////recob::TrackHitMeta metadata(hitIndex, -std::numeric_limits<double>::max());
	recob::TrackHitMeta metadata(index, -std::numeric_limits<double>::max());
	////anaTrackHitAssociationsMeta->addSingle(ptrack, phit, metadata);
	//auto vmeta = fmthm.data(ptrack.key());
	//std::vector<art::Ptr<recob::TrackHitMeta>> vmeta = fmthm.data(ptrack.key());
	//anaTrackHitAssociationsMeta->addSingle(pTrack, phit, metadata);
      }
    */


    // skip cc muon track
    if (ptrack.key()==trkmuon.key()) continue;
    
    std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());

    ReconstructionOrchestrator orchestrator;
    orchestrator.runReconstruction(sp_list, k_track);
    std::vector<Reco::Track> rebuildTrackList = orchestrator.getRebuildTrackList();

    //for(Reco::Track reco_track : rebuildTrackList) {
    for(int i=0; i < rebuildTrackList.size(); i++) {

	 anaTrackCollection->push_back(rebuildTrackList[i]);

	 std::vector<art::Ptr<recob::Hit>> hits_from_track_rebuild = trackHitLists[i];
	 
	 lar_pandora::HitVector anaHitCollection_rebuild_tmp;
	 for(auto hitptr : hits_from_track_rebuild){
	   anaHitCollection_rebuild_tmp.push_back(hitptr);
	 }
	 util::CreateAssn(*this, evt, *(anaTrackCollection.get()), anaHitCollection_rebuild_tmp, *(anaTrackHitAssociations.get()));

    }

   
    }
    evt.put(std::move(anaTrackCollection));
    evt.put(std::move(anaTrackHitAssociations));
  }
}



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
    
    //produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::Track> >(); 
    //produces< std::vector<recob::Hit> >();
    //produces< art::Assns<recob::PFParticle, recob::Track> >();
    produces< art::Assns<recob::Track, recob::Hit> >();
    //produces< art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();
    //produces< std::vector<recob::SpacePoint> >();
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
   
    //std::unique_ptr< std::vector<recob::Vertex> >          anaVertexCollection(new std::vector<recob::Vertex> ); 
    std::unique_ptr< std::vector<recob::Track> > anaTrackCollection(new std::vector<recob::Track> );
    //std::unique_ptr< std::vector<recob::Hit> >             anaHitCollection(new std::vector<recob::Hit> );
    //std::unique_ptr< art::Assns<recob::PFParticle, recob::Track> > anaParticleTrackAssociations( new art::Assns<recob::PFParticle, recob::Track> );
    std::unique_ptr<art::Assns<recob::Track, recob::Hit>>  anaTrackHitAssociations(new art::Assns<recob::Track, recob::Hit>);
    //std::unique_ptr< art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> > anaTrackHitAssociationsMeta( new art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> );

    //std::unique_ptr< std::vector<recob::SpacePoint> > anaSpacePointCollection(new std::vector<recob::SpacePoint> );
   
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

    
    // Check if event passed the NuCC inclusive filter
    /*
    reco_nu_cc_filter = false;
    string process(isMC ? "OverlayFiltersPostStage2" : "DataFiltersPostStage2");
    art::InputTag trigResInputTag("TriggerResults","",process.data()); // the last is the name of process where the filters were run
    art::ValidHandle<art::TriggerResults> trigRes = evt.getValidHandle<art::TriggerResults>(trigResInputTag);
    fhicl::ParameterSet pset;
    if (!fhicl::ParameterSetRegistry::get(trigRes->parameterSetID(), pset)) { throw cet::exception("PSet Not Found???"); }
    std::vector<std::string> trigger_path_names = pset.get<std::vector<std::string> >("trigger_paths", {});
    if (trigger_path_names.size()!=trigRes->size()) { throw cet::exception("Size mismatch???"); }
    for (size_t itp=0;itp<trigRes->size();itp++) {
      //cout << "Filter name " << trigger_path_names.at(itp) << " decision=" << trigRes->at(itp).accept() << endl;
      if (trigger_path_names.at(itp)=="NuCC") {
	reco_nu_cc_filter = trigRes->at(itp).accept();
      }
    }

    
    if (!reco_nu_cc_filter) {
      //evt.put(std::move(anaTrackCollection));
      //return;
    }
    */
    
    // Collect all recontructed particles
    art::Handle<std::vector<recob::PFParticle>> pfparticles;
    std::vector< art::Ptr<recob::PFParticle> > pfpVector;
    if(evt.getByLabel(m_pfp_producer, pfparticles)){
      art::fill_ptr_vector(pfpVector, pfparticles);
    }
    
    if (pfparticles->size()==0) {
      std::cout << "No PFParticles found" << std::endl;
      //evt.put(std::move(anaVertexCollection));
      evt.put(std::move(anaTrackCollection));
      //evt.put(std::move(anaHitCollection));
      //evt.put(std::move(anaParticleTrackAssociations)); 
      evt.put(std::move(anaTrackHitAssociations));
      //evt.put(std::move(anaTrackHitAssociationsMeta));
      //evt.put(std::move(anaSpacePointCollection));
      return;
      //continue;
    }

    // Get PFParticle associations
    art::FindManyP<anab::T0> pfp_muon_assn(pfparticles, evt, "NuCCproducer");
    if(!pfp_muon_assn.isValid()){
      //evt.put(std::move(anaVertexCollection));
      evt.put(std::move(anaTrackCollection));
      //evt.put(std::move(anaHitCollection));
      //evt.put(std::move(anaParticleTrackAssociations));
      evt.put(std::move(anaTrackHitAssociations));
      //evt.put(std::move(anaTrackHitAssociationsMeta));
      //evt.put(std::move(anaSpacePointCollection));
      return;
      //continue;
    }
    
    art::FindManyP<recob::Track> pfparticleTrackAssn(pfparticles, evt, "pandora");
    if(!pfparticleTrackAssn.isValid()){
      //evt.put(std::move(anaVertexCollection));
      evt.put(std::move(anaTrackCollection));
      //evt.put(std::move(anaHitCollection));
      //evt.put(std::move(anaParticleTrackAssociations));
      evt.put(std::move(anaTrackHitAssociations));
      //evt.put(std::move(anaTrackHitAssociationsMeta));
      //evt.put(std::move(anaSpacePointCollection));
      return;
    }
    
    art::FindManyP<recob::Vertex> pfparticleVertexAssn(pfparticles, evt, "pandora");
    if(!pfparticleVertexAssn.isValid()){
      //evt.put(std::move(anaVertexCollection));
      evt.put(std::move(anaTrackCollection));
      //evt.put(std::move(anaHitCollection));
      //evt.put(std::move(anaParticleTrackAssociations));
      evt.put(std::move(anaTrackHitAssociations));
      //evt.put(std::move(anaTrackHitAssociationsMeta));
      //evt.put(std::move(anaSpacePointCollection)); 
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
      //evt.put(std::move(anaVertexCollection));
      evt.put(std::move(anaTrackCollection));
      //evt.put(std::move(anaHitCollection));
      //evt.put(std::move(anaParticleTrackAssociations));
      evt.put(std::move(anaTrackHitAssociations));
      //evt.put(std::move(anaTrackHitAssociationsMeta));
      //evt.put(std::move(anaSpacePointCollection)); 
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
    //reco_nu_cc_nmue = pfmuons.size();
    
    if (pfmuons.size()!=1) {
      //evt.put(std::move(anaVertexCollection));
      evt.put(std::move(anaTrackCollection));
      //evt.put(std::move(anaHitCollection));
      //evt.put(std::move(anaParticleTrackAssociations));
      evt.put(std::move(anaTrackHitAssociations));
      //evt.put(std::move(anaTrackHitAssociationsMeta));
      //evt.put(std::move(anaSpacePointCollection));
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

    art::Handle< std::vector<recob::Cluster> > clusterHandle;
    std::vector< art::Ptr<recob::Cluster> > clusterVector;
    if(evt.getByLabel(fPandoraLabel,clusterHandle)){
      art::fill_ptr_vector(clusterVector,clusterHandle);
    }

    art::Handle< std::vector<recob::SpacePoint> > spacepointHandle;
    std::vector< art::Ptr<recob::SpacePoint> > spacepointVector;
    if(evt.getByLabel(fPandoraLabel,spacepointHandle)){
      art::fill_ptr_vector(spacepointVector,spacepointHandle);
    }

    //fSpacePointproducer

    //art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel);
    
    // get track associations
    /*
    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
    if(!fmcal.isValid()){
      //return;
    }
    */

    art::FindManyP<recob::Hit> hits_from_tracks(trackListHandle, evt, fHitTrackAssns);
    art::FindManyP<recob::Hit> hits_from_showers(showerListHandle, evt, fHitShowerAssns);
    if(!hits_from_tracks.isValid()){
      evt.put(std::move(anaTrackCollection));
      //evt.put(std::move(anaHitCollection));
      //evt.put(std::move(anaParticleTrackAssociations));
      evt.put(std::move(anaTrackHitAssociations));
      //evt.put(std::move(anaTrackHitAssociationsMeta));
      //evt.put(std::move(anaSpacePointCollection));
      return;
      //return;
    }

    /*
    art::FindManyP<anab::ParticleID> trackPIDAssn(trackListHandle, evt, fPIDLabel);
    if(!trackPIDAssn.isValid()){
      //return;
    }
    */

    art::FindManyP<recob::Cluster> clusters_per_pfparticle(pfparticles, evt, fPandoraLabel);
    art::FindManyP<recob::Hit> hits_per_cluster(clusterHandle, evt, fPandoraLabel);
    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > pfParticleToClustersMap;
    std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> > clusterToHitsMap;

    art::FindManyP<recob::SpacePoint> spacepoints_per_pfparticle(pfparticles, evt, fSpacePointproducer);
    art::FindManyP<recob::Hit> hits_per_spacepoint(spacepointHandle, evt, fSpacePointproducer);
    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>> > pfParticleToSpacePointsMap;
    std::map<art::Ptr<recob::SpacePoint>,  std::vector<art::Ptr<recob::Hit>> > spacepointToHitsMap;

    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMap;
    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsFromSpacePointsMap;
    std::map<art::Ptr<recob::PFParticle>, int> pfParticleToNHitsFromSpacePoints;

    if( !clusters_per_pfparticle.isValid() || !hits_per_cluster.isValid() || !spacepoints_per_pfparticle.isValid() || !hits_per_spacepoint.isValid()  ){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      //evt.put(std::move(anaTrackHitAssociationsMeta));
      return;
    }

    //std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMapSPC;

    //art::ValidHandle<std::vector<recob::Cluster>> const & clusterHandle = evt.getValidHandle<std::vector<recob::Cluster>>(m_pandoraLabel);
    //std::vector< art::Ptr<recob::Cluster> > clusterVector;
    //art::fill_ptr_vector(clusterVector,clusterHandle);

    const art::PtrMaker<recob::Track> makeTrackPtr(evt);

    for (unsigned int i=0; i<pfpVector.size(); ++i) {
      auto pfp = pfpVector[i];
      //art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);
      //pfParticleToClustersMap[*pfparticle] = clusters_per_pfparticle.at(pfparticle->key());
      pfParticleToClustersMap[pfp] = clusters_per_pfparticle.at(pfp.key());
      pfParticleToSpacePointsMap[pfp] = spacepoints_per_pfparticle.at(pfp.key());
    }

    for(size_t i=0; i< clusterVector.size(); ++i){
      auto cluster = clusterVector[i];
      clusterToHitsMap[cluster] = hits_per_cluster.at(cluster.key());
    }

    for(size_t i=0; i< spacepointVector.size(); ++i){
      auto spacepoint = spacepointVector[i];
      spacepointToHitsMap[spacepoint] = hits_per_spacepoint.at(spacepoint.key());

      //auto const& spkey = spacepointVector.at(i).key();
    }

    for (unsigned int i=0; i<pfpVector.size(); ++i) {
      auto pfp = pfpVector[i];
      //art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);
      //art::Ptr<recob::Track> ppfptrack = pfparticleTrackAssn.at(i).front();
      //std::vector<art::Ptr<recob::Cluster>> clusters_vec  = pfParticleToClustersMap[*pfparticle];
      std::vector<art::Ptr<recob::Cluster>> clusters_vec  = pfParticleToClustersMap[pfp];
      std::vector<art::Ptr<recob::SpacePoint>> spacepoints_vec  = pfParticleToSpacePointsMap[pfp];
      std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};

      lar_pandora::HitSet hitsInParticleSet;

      for (art::Ptr<recob::SpacePoint> spacepoint: spacepoints_vec){
	std::vector<art::Ptr<recob::Hit>> hits_vec =  spacepointToHitsMap[spacepoint];
	hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );
	//pfParticleToNHitsFromSpacePoints[pfp] += hits_vec.size();
	
	//cout << hits_vec.size() << endl;
	for(auto hit : hits_vec){
	  (void) hitsInParticleSet.insert(hit);
	}
	
      }
      pfParticleToNHitsFromSpacePoints[pfp] = hits_for_pfp.size();
      pfParticleToHitsFromSpacePointsMap[pfp] = hits_for_pfp;
      //pfParticleToHitsMapSPC[pfp] = hits_for_pfp;
      cout << "hits_for_pfp size before cluster " << hits_for_pfp.size() << endl;

      for (art::Ptr<recob::Cluster> cluster: clusters_vec){
	std::vector<art::Ptr<recob::Hit>> hits_vec =  clusterToHitsMap[cluster];
	//hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );

	for(auto hit : hits_vec){
	  //cout << "hitsInParticleSet.count(hit): " << hitsInParticleSet.count(hit) << endl;
	  if(hitsInParticleSet.count(hit) == 0) hits_for_pfp.push_back(hit);
	}
      }
      cout << "hits_for_pfp size after cluster " << hits_for_pfp.size() <<endl;
      //pfParticleToHitsMap[*pfparticle] = hits_for_pfp;
      pfParticleToHitsMap[pfp] = hits_for_pfp;
    }


    
    // find track multiplicity around vertex
    int NTracks=tracklist.size();
    int NShowers=showerlist.size();
    //int ntracks = 0;
    
    for (int i=0; i<reco_nu_ndaughters; i++) {

      //std::map<int,int> mrgidpdg;
      //std::map<int,int> mrgidpdg_3D;
      //std::map<int,TVector3> mrgidmom;
      //vector<TH1D*> h_angular_distribution_pfparticle_cheat;
      //vector<TH2D*> h_angular_distribution_pfparticle_cheat_3D;
      //vector<TH1D*> h_angular_distribution_pfparticle;
      //vector<TH2D*> h_angular_distribution_pfparticle_3D;
    //vector<TH2D*> h_angular_distribution_pfparticle_surface;
    //vector<TH3D*> h_angular_distribution_pfparticle_sphere;
    //vector<bool> v_trk_flg;
    //vector<bool> v_trk_flg_3D;
    //vector<bool> v_trk_flg_surface;
    vector<bool> v_trk_flg_peak;
    //vector<int> v_pdg;
    //vector<int> v_pdg_3D;
    //vector<int> v_pdg_peak;
    //vector<TVector2> view_peak_vector_cheat;
    std::map<double, TVector2, std::greater<>> view_peak_map;
    //std::map<int, std::map<double, TVector2, std::greater<>>> view_peak_map_cheat;
    //vector<TVector2> view_peak_vector;
    vector<TVector2> best_peak_bins;
    //std::map<int, vector<TVector2>> best_peak_bins_cheat;

    //std::map<int, std::map<int, double>> angular_distribution_mrgid_map;
    //std::map<int, std::map<int, std::map<int, double>>> angular_distribution_mrgid_map_3D;
    //std::map<int, double> angular_distribution_map_track;
    //std::map<int, double> angular_distribution_map_shower;
    std::map<int, std::map<int, double>> angular_distribution_map_3D_track;
    std::map<int, std::map<int, double>> angular_distribution_map_3D_shower;
    std::map<int, std::map<int, double>> angular_distribution_map_3D_pfparticle;

    std::vector<art::Ptr<recob::Hit>> unavailable_hit_list;
    std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list;
    //std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list_cheat;
    //std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list_cheat_pip;
    //std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list_cheat_mup;
    std::vector<std::vector<art::Ptr<recob::Hit>>> shower_spine_hit_list_vector;
    //std::vector<art::Ptr<recob::Hit>> true_pi0_hit_list;

    std::vector<art::Ptr<recob::Hit>> hits_from_reco;

    //int trackcount = 0;

    //int currentMergedId = 1;
    //int currentMergedId_3D = 1;

    //lar_pandora::PFParticleVector pfParticleVector;
    //lar_pandora::PFParticlesToSpacePoints pfParticlesToSpacePoints;
    //lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPFParticleLabel, pfParticleVector, pfParticlesToSpacePoints);


    art::Ptr<recob::Track> ptrack(trackListHandle,reco_nu_daughters_id[i]);
    const recob::Track& track = *ptrack;
    lar_pandora::HitVector hitsInParticle;

    //anaTrackCollection->push_back(track);
    //anaTrackCollection->emplace_back(track);
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
    

    /*
    for (unsigned int i=0; i<pfparticles->size(); ++i) {
      art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);
      if (pfparticle->Parent()==pfnu->Self() && pfparticleTrackAssn.at(i).size()==1){
	art::Ptr<recob::Track> ppfptrack = pfparticleTrackAssn.at(i).front();
	const recob::Track& pfptrack = *ppfptrack; 

	//int id = 0;
	if(pfptrack.ID() == track.ID()){
	  
	  //pandora::IntVector indexVector;
	  //indexVector.pushback(id);
	  //lar_content::LArPfoHelper::GetSlidingFitTrajectory(cartesianPointVector, track.Vertex(), 20, wirePitchW, trackStateVector, &indexVector);

	  lar_pandora::PFParticlesToSpacePoints::const_iterator particleToSpacePointIter(pfParticlesToSpacePoints.find(pfparticle));
	  // PFParticlesToSpacePoints::const_iterator
	  //auto particleToSpacePointIter(pfParticlesToSpacePoints.find(pfparticle));
	  
	  lar_pandora::HitVector hitsFromSpacePoints;// hitsInParticle;
	  hitsInParticle.clear();
	  lar_pandora::LArPandoraHelper::GetAssociatedHits(evt, fPFParticleLabel, particleToSpacePointIter->second, hitsFromSpacePoints, &indexVector);
	  //lar_pandora::LArPandoraHelper::GetAssociatedHits(evt, fPFParticleLabel, particleToSpacePointIter, hitsFromSpacePoints, &indexVector);
	  
	  for (unsigned int hitIndex = 0; hitIndex < hitsFromSpacePoints.size(); hitIndex++){
	    hitsInParticle.push_back(hitsFromSpacePoints.at(hitIndex));
	  }
	  
	  hitsInParticle.clear();
	  hitsInParticle = pfParticleToHitsMap[*pfparticle];
	  continue;
	}

      }
    }
    */

    //art::FindManyP<recob::Hit> hits_from_tracks(trackListHandle, evt, fHitTrackAssns); 
    std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());



    //int hitsFromSpacePointsSize = 0;

    //unsigned int hitsFromSpacePointsSize = 0;
    art::FindManyP<recob::SpacePoint> spacepoint_per_hit_test(hitListHandle, evt, fSpacePointproducer);

    //recob::Track primary_track = trackRebuid(trackcount, hits_from_track, spacepoint_per_hit_test, track);

    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec_test;

    for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
      spacepoint_vec_test.clear();
      spacepoint_vec_test = spacepoint_per_hit_test.at(hits_from_track[i_h].key());

      //hitsFromSpacePointsSize += spacepoint_vec_test.size();
      //for(auto spaceptr : spacepoint_vec_test){
	//hitsFromSpacePointsSize++;
  //anaSpacePointCollection->push_back(*spaceptr);
      //}
    }
    cout << "hits_from_track.size(): " << hits_from_track.size() << endl;
    cout << "hitsFromSpacePointsSize: " << hitsFromSpacePointsSize << endl;

    //cout << "hitPtrVec[0].key(): " << hits_from_track[0].key() << endl;

    //get associated hits to track
    //std::vector<art::Ptr<recob::Hit>> hitPtrVec = trackToHitAssns.at(ptrack.key());
    //cout << "hitPtrVec[0].key(): " << hitPtrVec[0].key() << endl;
    
    //art::FindMany<recob::Hit> trackToHitAssns(trackListHandle,  evt, fHitTrackAssns);
    //auto hitPtrVec = trackToHitAssns.at(ptrack.key());
    //size_t prev_size = anaHitCollection->size();

    /*
    for(auto hitptr : hitPtrVec){
      anaHitCollection->push_back(*hitptr);
    }
    */

    //push back the hit
    //size_t prev_size = anaHitCollection->size(); 
    //std::unique_ptr< std::vector<recob::Hit> > anaHitCollection_tmp( new std::vector<recob::Hit> );
    lar_pandora::HitVector anaHitCollection_tmp;
    anaHitCollection_tmp.clear();

    for(auto hitptr : hits_from_track){
      //anaHitCollection->push_back(*hitptr);
      //cout << "hit.key(): " << hitptr->key() << endl;
      anaHitCollection_tmp.push_back(hitptr);
    }

    /*
    PFParticleVector pfParticleVector;
    PFParticlesToSpacePoints pfParticlesToSpacePoints;
    LArPandoraHelper::CollectPFParticles(evt, fPFParticleLabel, pfParticleVector, pfParticlesToSpacePoints); 

    PFParticlesToSpacePoints::const_iterator particleToSpacePointIter(pfParticlesToSpacePoints.find(pPFParticle));

    for (const art::Ptr<recob::PFParticle> pPFParticle : pfParticleVector){
      // Obtain associated spacepoints
      PFParticlesToSpacePoints::const_iterator particleToSpacePointIter(pfParticlesToSpacePoints.find(pPFParticle));

      HitVector hitsFromSpacePoints;
      LArPandoraHelper::GetAssociatedHits(evt, m_pfParticleLabel, particleToSpacePointIter->second, hitsFromSpacePoints, track.key());
    }
    */


    //util::CreateAssn(*this, evt, *(anaTrackCollection.get()), hitsInParticle, *(anaTrackHitAssociations.get()));

    art::Ptr<recob::Track> pTrack(makeTrackPtr(anaTrackCollection->size() - 1));

    //ATTN: metadata added with index from space points if available, null for others                                                                        
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

    //fmthm.at(ptrack.key())
    //std::vector<art::Ptr<recob::TrackHitMeta>> vmeta = fmthm.at(ptrack.key();
    //								  vmeta[hitIndex];
    
    /*

    util::CreateAssn(*this, evt, pTrack, pPFParticle, *(outputParticlesToTracks.get()));
    util::CreateAssn(*this, evt, *(outputTracks.get()), hitsInParticle, *(outputTracksToHits.get()));
    */

    //util::CreateAssn(*this,evt,*anaTrackCollection,*anaHitCollection,*anaTrackHitAssociations,prev_size,anaHitCollection->size()); 
    //util::CreateAssn(*this, evt, *(anaTrackCollection.get()), anaHitCollection_tmp, *(anaTrackHitAssociations.get()));
    //cout <<"create assns" << endl;
    //util::CreateAssn(*this, evt, *(anaTrackCollection.get()), hitsInParticle, *(anaTrackHitAssociations.get()));
    //cout <<"end assns" << endl;

    // skip cc muon track
    if (ptrack.key()==trkmuon.key()) continue;


    

    // check track start and end
    TVector3 pos(track.Vertex().X(),track.Vertex().Y(),track.Vertex().Z());
    TVector3 end(track.End().X(),track.End().Y(),track.End().Z());

    /*
    double st_vtx=TMath::Sqrt((reco_nu_vtx_x-pos.X())*(reco_nu_vtx_x-pos.X()) +
                              (reco_nu_vtx_y-pos.Y())*(reco_nu_vtx_y-pos.Y()) +
                              (reco_nu_vtx_z-pos.Z())*(reco_nu_vtx_z-pos.Z()));
    */

    // track length and angles
    //double trklen=track.Length();
    /*
    TVector3 endmom = track.EndMomentumVector<TVector3>(); 
    TVector3 strmom = track.StartMomentumVector<TVector3>(); 
    TVector3 vtxmom = track.VertexMomentumVector<TVector3>(); 
    */

    // check kaon track start and end distance from vertex
    /*
    double end_dis=TMath::Sqrt((reco_nu_vtx_x-end.X())*(reco_nu_vtx_x-end.X()) +
                               (reco_nu_vtx_y-end.Y())*(reco_nu_vtx_y-end.Y()) +
                               (reco_nu_vtx_z-end.Z())*(reco_nu_vtx_z-end.Z()));
    */

    //fillCalorimetry(fmcal.at(ptrack.key()),ntracks);




    int ndaughters = 0;
    int ndaughters_sh = 0;

    cout << "before NTracks loop" << endl;
    cout << "Ntracks: " << NTracks << endl;

    for (int j=0; j<NTracks; j++) {
    
      art::Ptr<recob::Track> ptrack_dau(trackListHandle,j);
      const recob::Track& track_dau = *ptrack_dau;
    
      // skip all primary tracks
      bool skip = false;
      for (int k=0; k<reco_nu_ndaughters; k++) {
        if (int(ptrack_dau.key())==reco_nu_daughters_id[k]) {
          skip=true;
          break;
        }
      }
      if (skip) continue;

      TVector3 pos2(track_dau.Vertex().X(),track_dau.Vertex().Y(),track_dau.Vertex().Z());

      double track_dau_distance=TMath::Sqrt((end.X()-pos2.X())*(end.X()-pos2.X()) +
                                         (end.Y()-pos2.Y())*(end.Y()-pos2.Y()) +
                                         (end.Z()-pos2.Z())*(end.Z()-pos2.Z()));

      cout << "distance to vertex: " << track_dau_distance << endl;

      // check distance to vertex
      //if (track_dau_distance<10) { //7cm, 
      if (track_dau_distance<40) { //7cm, 
      //if (track_dau_distance<100) { //7cm, 

	//erase-remove this daughter track to avoid double counting
	/*
	auto customComparator = [&track_dau](const recob::Track& element) {
	  return track_dau.ID() == element.ID();
	};
	auto removeEnd = std::remove_if(anaTrackCollection->begin(), anaTrackCollection->end(), customComparator);
	anaTrackCollection->erase(removeEnd, anaTrackCollection->end());
	*/

        TVector3 end2(track_dau.End().X(),track_dau.End().Y(),track_dau.End().Z());

	//if(true_kaon_end_process == 1 || true_kaon_end_process == 1){
	//if(reco_track_true_pdg[ntracks] == 321){
	std::vector<art::Ptr<recob::Hit>> hits_from_dau_track = hits_from_tracks.at(ptrack_dau.key());
	art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);
	hits_from_reco.insert(hits_from_reco.end(), hits_from_dau_track.begin(), hits_from_dau_track.end());
	
	//fillAngularDistributionMap(hits_from_track, end, spacepoint_per_hit, angular_distribution_map_track);
	fillAngularDistributionMap3D(hits_from_dau_track, end, spacepoint_per_hit, angular_distribution_map_3D_track);
	//fillAngularDistributionMapCheat(hits_from_track, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map, mrgidpdg, currentMergedId);
	//fillAngularDistributionMapCheat3D(hits_from_track, true_pi0_hit_list, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map_3D, mrgidpdg_3D, mrgidmom, currentMergedId_3D);
	//}
	//}

	ndaughters++;

      }//gap between K+ end and daughter

    }//NTrack loop
  


    for (int j_s=0; j_s<NShowers; j_s++) {
      art::Ptr<recob::Shower> pshower_dau(showerListHandle,j_s);
      const recob::Shower& shower_dau = *pshower_dau;
      
      // skip all primary showers                                                                                                                           
      bool skip2 = false;
      for (int k_sh=0; k_sh<reco_nu_ndaughters; k_sh++) {
	if (int(pshower_dau.key())==reco_nu_daughters_id[k_sh]) {
	  skip2=true;
	  break;
	}
      }
      if (skip2) continue;
      TVector3 pos2_sh(shower_dau.ShowerStart().X(),shower_dau.ShowerStart().Y(),shower_dau.ShowerStart().Z());

      double shower_dau_distance=TMath::Sqrt((end.X()-pos2_sh.X())*(end.X()-pos2_sh.X()) +
                                            (end.Y()-pos2_sh.Y())*(end.Y()-pos2_sh.Y()) +
                                            (end.Z()-pos2_sh.Z())*(end.Z()-pos2_sh.Z()));

      //if (ndaughters_sh<=50 && shower_dau_distance<100) {  
      if (shower_dau_distance<40) {  

	if(ndaughters_sh>20) break;

	std::vector<art::Ptr<recob::Hit>> hits_from_shower = hits_from_showers.at(pshower_dau.key());
	art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);
	hits_from_reco.insert(hits_from_reco.end(), hits_from_shower.begin(), hits_from_shower.end()); 
	
	//fillTrueMatching_sh(hits_from_shower, particles_per_hit, ntracks, ndaughters_sh);
	//fillShowerMatching(hits_from_shower, particles_per_hit, ntracks, ndaughters_sh);
	
	//if(true_kaon_end_process == 1 || true_kaon_end_process == 1){
	//if(reco_track_true_pdg[ntracks] == 321){	
	//fillAngularDistributionMap(hits_from_shower, end, spacepoint_per_hit, angular_distribution_map_shower);
	fillAngularDistributionMap3D(hits_from_shower, end, spacepoint_per_hit, angular_distribution_map_3D_shower);
	//fillAngularDistributionMapCheat(hits_from_shower, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map, mrgidpdg, currentMergedId);
	//fillAngularDistributionMapCheat3D(hits_from_shower, true_pi0_hit_list, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map_3D, mrgidpdg_3D, mrgidmom, currentMergedId_3D);
	
	//}
	//}

	ndaughters_sh++;
      }
    }

    //smoothAngularDistributionMapCheat3D(angular_distribution_mrgid_map_3D);
    smoothAngularDistributionMap3D(angular_distribution_map_3D_track);
    smoothAngularDistributionMap3D(angular_distribution_map_3D_shower);

    accumulateAngularDistributionMap3D(angular_distribution_map_3D_track, angular_distribution_map_3D_shower, angular_distribution_map_3D_pfparticle);

    //obtainPeakVectorCheat3D(angular_distribution_mrgid_map_3D, mrgidpdg_3D, view_peak_map_cheat, v_pdg_peak);
  
    obtainPeakVector3D(angular_distribution_map_3D_track, v_trk_flg_peak, view_peak_map, true);
    obtainPeakVector3D(angular_distribution_map_3D_shower, v_trk_flg_peak, view_peak_map, false);
    //obtainPeakVector3D(angular_distribution_map_3D_pfparticle, v_trk_flg_peak, view_peak_map, false);

    //findBestAngularPeakCheat3D(view_peak_map_cheat, best_peak_bins_cheat);
    findBestAngularPeak3D(angular_distribution_map_3D_pfparticle, view_peak_map, best_peak_bins);

    art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);
    findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_vector, end, view_peak_map, best_peak_bins);


    /*
   if(best_peak_bins_cheat.size()!=0){
      for(auto& entry : best_peak_bins_cheat){
	cout << "PEAK OF PDG IS " << entry.first << endl;
	if(reco_track_true_pdg[ntracks] == 321){
	  if(entry.first == 211) findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat_pip, end, view_peak_map, entry.second);
	  if(entry.first == -13) findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat_mup, end, view_peak_map, entry.second);
	}
	//findShowerSpine3D(true_pi0_hit_list, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat, end, view_peak_map, entry.second);
	//findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat, end, view_peak_map, entry.second);

	//findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat, end, view_peak_map, entry.second, true_dau_pip_length);
      }
    }
    */
   

   //double trkln;

   cout << "best_peak_bins.size(): " << best_peak_bins.size() << endl;
   cout << "shower_spine_hit_list_vector.size(): " << shower_spine_hit_list_vector.size() << endl;

   if(best_peak_bins.size()!=0){
     int i=0;
     for(auto const& entry : best_peak_bins){	
       
       cout << "entry.X(): " << entry.X() << endl;
       
       if(shower_spine_hit_list_vector.size()==0) continue;
       cout << "shower_spine_hit_list_vector[i]: " << shower_spine_hit_list_vector[i].size() << endl;
	 if(shower_spine_hit_list_vector[i].size()<10 || shower_spine_hit_list_vector[i].size() > 10000) continue;
	 cout << "shower_spine_hit_list_vector[i]: " << shower_spine_hit_list_vector[i].size() << endl;
	 //trkln = trackRebuid(shower_spine_hit_list_vector[i], spacepoint_per_hit, track);
	 
	 recob::Track reco_track = trackRebuid(shower_spine_hit_list_vector[i], spacepoint_per_hit, track);
	 
	 cout << "adding reco_track" << endl;
	 anaTrackCollection->push_back(reco_track);

	 //get associated hits to track
	 //art::FindMany<recob::Hit> trackToHitAssns(trackListHandle,  evt, fTrackModuleLabel);
	 //auto hitPtrVec_reco = shower_spine_hit_list_vector[i];
	 std::vector<art::Ptr<recob::Hit>> hits_from_track_rebuild = shower_spine_hit_list_vector[i];
	 // = trackToHitAssns.at(ptrack.key());
	 
	 //push back the hit
	 /*
	 size_t prev_size_reco = anaHitCollection->size(); 
	 for(auto hitptr : hitPtrVec_reco)
	   anaHitCollection->push_back(*hitptr);
	 util::CreateAssn(*this,evt,*anaTrackCollection,*anaHitCollection,*anaTrackHitAssociations,prev_size_reco,anaHitCollection->size()); 
	 */

	 
	 art::FindManyP<recob::SpacePoint> spacepoint_per_hit_reco(hitListHandle, evt, fSpacePointproducer);
	 std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec_reco;
	 unsigned int hitsFromSpacePointsRecoSize = 0;
	 for(size_t i_h=0; i_h<hits_from_track_rebuild.size(); i_h++){
	   spacepoint_vec_reco.clear();
	   spacepoint_vec_reco = spacepoint_per_hit_reco.at(hits_from_track_rebuild[i_h].key());
	   //for(spacepoint_vec_reco_it : spacepoint_vec_reco){
	   hitsFromSpacePointsRecoSize += spacepoint_vec_reco.size();
	     //}
	 }
	 

	 /*
	 for (unsigned int hitIndex = 0; hitIndex < hits_from_track_rebuild.size(); hitIndex++){
	   const art::Ptr<recob::Hit> phit(hits_from_track_rebuild.at(hitIndex));
	   recob::TrackHitMeta metadata(hitIndex, -std::numeric_limits<double>::max());
	   //recob::Track reco_track
	   art::Ptr<recob::Track> preco_track = &reco_track;
	   anaTrackHitAssociationsMeta->addSingle(preco_track, phit, metadata);
	 }
	 */
	 
	 art::Ptr<recob::Track> pTrackdau(makeTrackPtr(anaTrackCollection->size() - 1));

	 lar_pandora::HitVector anaHitCollection_rebuild_tmp;
	 for(auto hitptr : hits_from_track_rebuild){
	   //anaHitCollection->push_back(*hitptr);
	   anaHitCollection_rebuild_tmp.push_back(hitptr);
	 }
	 util::CreateAssn(*this, evt, *(anaTrackCollection.get()), anaHitCollection_rebuild_tmp, *(anaTrackHitAssociations.get()));

	 for (unsigned int hitIndex = 0; hitIndex < hits_from_track_rebuild.size(); hitIndex++){
	   const art::Ptr<recob::Hit> phit(hits_from_track_rebuild.at(hitIndex));
	   const int index((hitIndex < hitsFromSpacePointsRecoSize) ? hitIndex : std::numeric_limits<int>::max());
	   recob::TrackHitMeta metadata(index, -std::numeric_limits<double>::max());
	   //anaTrackHitAssociationsMeta->addSingle(pTrackdau, phit, metadata);
	 }
	 
	 //best_peak_trkln[ntracks][i] = trkln;
	 
	 i++;
     }
   }
   
    }
    evt.put(std::move(anaTrackCollection));
    //evt.put(std::move(anaHitCollection));
    //evt.put(std::move(anaParticleTrackAssociations));
    evt.put(std::move(anaTrackHitAssociations));
    //evt.put(std::move(anaTrackHitAssociationsMeta));
    //evt.put(std::move(anaSpacePointCollection));
    //return;
  }
}


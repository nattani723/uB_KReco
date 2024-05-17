#include "CCKaonProducer_module.h"
#include "TrackRebuilder/ReconstructionOrchestrator.cc"

//#ifdef __MAKECINT__
#ifdef __CLING__
#pragma link C++ class std::vector < std::vector<Float_t> >+; 
#pragma link C++ class std::vector < std::vector< std::vector<Float_t> > >+; 
#endif


namespace kaon_reconstruction {

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
    
    produces< std::vector<recob::Track> >(); 
    produces< art::Assns<recob::Track, recob::Hit> >();

    theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
    //SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    //geom = lar::providerFrom<geo::Geometry>();    

  }
  
  CCKaonProducer::~CCKaonProducer()
  {
    //destructor
  }


  void CCKaonProducer::beginJob()
  {
    //set output tree
  }


  void CCKaonProducer::produce(art::Event& evt)
  {
    std::vector< art::Ptr<simb::MCParticle> > ptList;
   
    std::unique_ptr< std::vector<recob::Track> > anaTrackCollection(new std::vector<recob::Track>);
    std::unique_ptr<art::Assns<recob::Track, recob::Hit>>  anaTrackHitAssociations(new art::Assns<recob::Track, recob::Hit>);
   
    // get MCTruth information
    /*
    art::Handle< std::vector<simb::MCTruth> > mctruths;
    evt.getByLabel(fGenieGenModuleLabel, mctruths);
    if (!mctruths.isValid() || mctruths->empty()) return;
    const simb::MCTruth& mctruth = mctruths->front();
    */
    /*
    if (mctruths->size()!=1) {
      //return;
      //continue;
    }
    simb::MCTruth mctruth = mctruths->at(0);
    */


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
    art::FindManyP<recob::Track> pfparticleTrackAssn(pfparticles, evt, fTrackModuleLabel);
    art::FindManyP<recob::Vertex> pfparticleVertexAssn(pfparticles, evt, fTrackModuleLabel);
    if(!pfp_muon_assn.isValid() || !pfparticleTrackAssn.isValid() || !pfparticleVertexAssn.isValid()){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }


    // Find the primary neutrino particle and associated muon
    lar_pandora::PFParticleVector pfneutrinos(0);
    lar_pandora::PFParticleVector pfmuons(0);
    std::vector<int> reco_nu_daughters_id(0);

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


    if (pfmuons.size()!=1) {
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }

    
    reco_nu_ndaughters = reco_nu_daughters_id.size();

    art::Ptr<recob::PFParticle> pfmuon = pfmuons.front();
    art::Ptr<recob::Track> trkmuon = pfparticleTrackAssn.at(pfmuon.key()).front();
    

    // get hits
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitList;
    if(evt.getByLabel(fHitsModuleLabel,hitListHandle)){
      art::fill_ptr_vector(hitList, hitListHandle);
    }

    // get tracks
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > trackList;
    if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) {
      art::fill_ptr_vector(trackList, trackListHandle);
    }

    // get showers
    art::Handle< std::vector<recob::Shower> > showerListHandle;
    std::vector<art::Ptr<recob::Shower> > showerList;
    if(evt.getByLabel(fShowerModuleLabel,showerListHandle)) {
      art::fill_ptr_vector(showerList, showerListHandle);
    }

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

    art::FindOneP<recob::Hit> findSPToHit(spacepointVector, evt, fSpacePointproducer); 
    art::FindManyP<recob::Hit> findTrackToHit(trackList, evt, fTrackModuleLabel);
    art::FindManyP<recob::Hit> findShowerToHit(showerList, evt, fShowerModuleLabel);

    std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>> hitToSpacePointMap;
    std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>> spacepointToHitMap;

    if( !findSPToHit.isValid()  ){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }

    
    for (unsigned int iSP = 0; iSP < spacepointVector.size(); ++iSP) { 
      const art::Ptr<recob::SpacePoint> spacepoint = spacepointVector.at(iSP);
      const art::Ptr<recob::Hit> hit = findSPToHit.at(iSP);
      spacepointToHitMap[spacepoint] = hit;
      hitToSpacePointMap[hit] = spacepoint;
    }
    

    //loop over primary nu track
    for (int i=0; i<reco_nu_ndaughters; i++) {

      art::Ptr<recob::Track> ptrack(trackListHandle,reco_nu_daughters_id[i]);
      //const recob::Track& track = *ptrack;
      
      
      // skip cc muon track
      if (ptrack.key()==trkmuon.key()) continue;
      
      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());
            
      ReconstructionOrchestrator orchestrator;
      orchestrator.runReconstruction(spacepointList, spacepointToHitMap, hitToSpacePointMap, ptrack, hits_from_track);
      //orchestrator.runReconstruction(spacepointVector, spacepointToHitMap, hitToSpacePointMap, ptrack, hits_from_track);

      std::vector<recob::Track> rebuildTrackList = orchestrator.getRebuildTrackList();

      std::vector<std::vector<art::Ptr<recob::Hit>>> trackHitLists = orchestrator.getHitLists();

      //for(Reco::Track reco_track : rebuildTrackList) {
      for(unsigned int i=0; i < rebuildTrackList.size(); i++) {
	
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
}// namespace kaon_reconstruction


//DEFINE_ART_MODULE(CCKaonProducer)

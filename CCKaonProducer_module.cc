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

    /*
    for (unsigned int iSP = 0; iSP < spacepointVector.size(); ++iSP) { 
      const art::Ptr<recob::SpacePoint> spacepoint = spacepointVector.at(iSP);
      const art::Ptr<recob::Hit> hit = findSPToHit.at(iSP);
      spacepointToHitMap[spacepoint] = hit;
      hitToSpacePointMap[hit] = spacepoint;
    }
    */
    
    art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, evt, fHitTruthAssns);
    std::vector< art::Ptr<recob::SpacePoint> > spacepointFromRecoObject;
    std::vector< art::Ptr<recob::SpacePoint> > spacepointFromMu;
    std::vector< art::Ptr<recob::SpacePoint> > spacepointFromPi;
    std::vector< art::Ptr<recob::Hit> > hitFromTrack;
    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;


    for (unsigned int itrk=0; itrk < trackList.size(); ++itrk) { 
      
      const art::Ptr<recob::Track> track = trackList.at(itrk);

      //skip primary track
      bool primary = false; 
      for (int i_nutrk=0; i_nutrk<reco_nu_ndaughters; i_nutrk++) {
	if (int(track.key())==reco_nu_daughters_id[i_nutrk]) {
	  primary=true;
	  break;
	}
      }
      
      if (track.key()==trkmuon.key() || primary){
	std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(track.key());
	hitFromTrack.insert(hitFromTrack.end(), hits_from_track.begin(), hits_from_track.end());
      }
    }
    
    for (unsigned int iSP = 0; iSP < spacepointVector.size(); ++iSP) { 
      art::Ptr<recob::SpacePoint> spacepoint = spacepointVector.at(iSP);
      art::Ptr<recob::Hit> hit = findSPToHit.at(iSP);
      spacepointToHitMap[spacepoint] = hit;
      hitToSpacePointMap[hit] = spacepoint;
      
      //skip ccmu primary hits
      if( std::find(hitFromTrack.begin(), hitFromTrack.end(), hit) != hitFromTrack.end()) continue;
      
      spacepointFromRecoObject.push_back(spacepoint);
      
      simb::MCParticle const* mcparticle = truthMatchHit(hit, particles_per_hit);
      if(!mcparticle) continue;
      if(mcparticle->PdgCode()==-13) spacepointFromMu.push_back(spacepoint);
      if(mcparticle->PdgCode()==211) spacepointFromPi.push_back(spacepoint);
      
    }

    /*
    //fill the SP-Hit map 
    for (unsigned int itrk=0; itrk < trackList.size(); ++itrk) {

      const art::Ptr<recob::Track> track = trackList.at(itrk);

      if (track.key()==trkmuon.key()) continue;

      //skip primary track
      bool skip = false;
      for (int i_nutrk=0; i_nutrk<reco_nu_ndaughters; i_nutrk++) {
	if (int(track.key())==reco_nu_daughters_id[i_nutrk]) { 
	  skip=true;
	  break;
	}
      }
      if (skip) continue;

      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(track.key());

      art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);

      for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){

	spacepoint_vec.clear();
	spacepoint_vec = spacepoint_per_hit.at(hits_from_track[i_h].key());
	
	if(spacepoint_vec.size()!=1) continue;
	art::Ptr<recob::SpacePoint> spacepoint = spacepoint_vec.at(0);
	art::Ptr<recob::Hit> hit = hits_from_track.at(i_h);
			
	spacepointFromRecoObject.push_back(spacepoint);
	spacepointToHitMap[spacepoint] = hit;
	hitToSpacePointMap[hit] = spacepoint;

	simb::MCParticle const* mcparticle = truthMatchHit(hit, particles_per_hit);
	if(!mcparticle) continue;
	if(mcparticle->PdgCode()==-13) spacepointFromMu.push_back(spacepoint);
	if(mcparticle->PdgCode()==211) spacepointFromPi.push_back(spacepoint);

      }

    }

    
    for (unsigned int ishw=0; ishw < showerList.size(); ++ishw) {

      const art::Ptr<recob::Shower> shower = showerList.at(ishw);

      std::vector<art::Ptr<recob::Hit>> hits_from_shower = hits_from_showers.at(shower.key());

      art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);

      for(size_t i_h=0; i_h<hits_from_shower.size(); i_h++){

	spacepoint_vec.clear();
	spacepoint_vec = spacepoint_per_hit.at(hits_from_shower[i_h].key());
	
	if(spacepoint_vec.size()!=1) continue;
	art::Ptr<recob::SpacePoint> spacepoint = spacepoint_vec.at(0);
	art::Ptr<recob::Hit> hit = hits_from_shower.at(i_h);

	spacepointFromRecoObject.push_back(spacepoint);
	spacepointToHitMap[spacepoint] = hit;
	hitToSpacePointMap[hit] = spacepoint;

	simb::MCParticle const* mcparticle = truthMatchHit(hit, particles_per_hit);
	if(!mcparticle) continue;
	if(mcparticle->PdgCode()==-13) spacepointFromMu.push_back(spacepoint);
	if(mcparticle->PdgCode()==211) spacepointFromPi.push_back(spacepoint);

      }

    }
    */


    //loop over primary nu track
    for (int i=0; i<reco_nu_ndaughters; i++) {

      art::Ptr<recob::Track> ptrack(trackListHandle,reco_nu_daughters_id[i]);
      //const recob::Track& track = *ptrack;
      
      
      // skip cc muon track
      if (ptrack.key()==trkmuon.key()) continue;
      
      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());

      //simb::MCParticle const* mcparticle = truthMatchTrack(hits_from_track, particles_per_hit);
      //if(mcparticle && mcparticle->PdgCode()==321) std::cout << "This is Primary Kaon" << std::endl;
      //if(mcparticle->PdgCode()!=321) continue;
      //if(true_kaon_end_process!=0) continue; 
      /*
      simb::MCParticle const* mcparticle = truthMatchTrack(hits_from_track, particles_per_hit);
      if(mcparticle){
	if(mcparticle->PdgCode()!=321) continue;
      } 
      */
 
      ReconstructionOrchestrator orchestrator;
      orchestrator.runReconstruction(spacepointFromRecoObject, spacepointToHitMap, hitToSpacePointMap, ptrack, hits_from_track);
      //orchestrator.runReconstruction(spacepointFromMu, spacepointToHitMap, hitToSpacePointMap, ptrack, hits_from_track);

      std::vector<recob::Track> rebuildTrackList = orchestrator.getRebuildTrackList();

      std::vector<std::vector<art::Ptr<recob::Hit>>> trackHitLists = orchestrator.getHitLists();

      //for(Reco::Track reco_track : rebuildTrackList) {
      for(unsigned int j=0; j < rebuildTrackList.size(); j++) {
	
	if(trackHitLists[j].empty()) continue;

	anaTrackCollection->push_back(rebuildTrackList[j]);

	std::vector<art::Ptr<recob::Hit>> hits_from_track_rebuild = trackHitLists[j];
	
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

  simb::MCParticle const* CCKaonProducer::truthMatchTrack(std::vector<art::Ptr<recob::Hit>>& hits_from_track, art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit){

    simb::MCParticle const* matched_mcparticle = NULL;

    std::unordered_map<int,double> trkide;
    double maxe=-1, tote=0;
    std::vector<simb::MCParticle const*> particle_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
  
    for(size_t i_h=0; i_h<hits_from_track.size(); i_h++) {
      particle_vec.clear(); match_vec.clear();
      particles_per_hit.get(hits_from_track[i_h].key(), particle_vec, match_vec);

      //cout << "track's particle_vec.size(): " << particle_vec.size() << endl;

      for(size_t i_p=0; i_p<particle_vec.size(); ++i_p) {
	trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
	//trkide[ particle_vec[i_p]->TrackId() ] ++;
	tote += match_vec[i_p]->energy;
	if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
	  maxe = trkide[ particle_vec[i_p]->TrackId() ];
	  matched_mcparticle = particle_vec[i_p];
	}
      }
    }
    if(matched_mcparticle) return matched_mcparticle;
    else return NULL;
  }

  simb::MCParticle const* CCKaonProducer::truthMatchHit(art::Ptr<recob::Hit>& hit, art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit){

    simb::MCParticle const* matched_mcparticle = NULL;

    double maxe=-1;
    std::vector<simb::MCParticle const*> particle_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

    particle_vec.clear();
    particles_per_hit.get(hit.key(), particle_vec, match_vec);

    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p) {
      if(match_vec[i_p]->energy>maxe){
	maxe = match_vec[i_p]->energy;
	matched_mcparticle = particle_vec[i_p];
      }
    }
    if(matched_mcparticle) return matched_mcparticle;
    else return NULL;
  }

}// namespace kaon_reconstruction


//DEFINE_ART_MODULE(CCKaonProducer)

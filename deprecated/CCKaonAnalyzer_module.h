#ifndef CCKAON_ANALYSIS
#define CCKAON_ANALYSIS

#include "TH1.h"
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
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <string.h>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h" 

#include "larcoreobj/SummaryData/POTSummary.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include  "nusimdata/SimulationBase/GTruth.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "/exp/uboone/app/users/taniuchi/51_pandora/srcs/larpandoracontent/larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

//#include "larpandora/LArPandoraEventBuilding/LArPandoraTrackCreation_module.cc"

#include "larcore/Geometry/Geometry.h"

#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Helper function for PID stuff
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphDelaunay.h"
#include "TRandom3.h"

#include "Pandora/PdgTable.h"
#include <chrono>

#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <sys/stat.h>

//------------------------------------------------------------------------------------------------------------------------------------------


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/TriggerResults.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "larcore/Geometry/Geometry.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "headers/LLR_PID.h"
#include "headers/LLRPID_proton_muon_lookup.h"
#include "headers/LLR_PID_K.h"
#include "headers/LLRPID_kaon_proton_lookup.h"
#include "objects/RecoParticle.h"
#include "objects/Helpers.h"
#include "headers/ParticleTypes.h"
#include "headers/FV.h"
#include "algorithms/MeandEdX.h"
#include "algorithms/ThreePlaneMeandEdX.h"

//#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "TTree.h"
#include "TMath.h"

#include <array>
#include <vector>
#include <map>
//#include "LinkDef.h"


#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/make_tool.h"


#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/Deprecated/BezierTrack.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "ubobj/Optical/UbooneOpticalFilter.h"
#include "ubana/AnalysisTree/MCTruth/IMCTruthMatching.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
    

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "canvas/Persistency/Common/TriggerResults.h" 
#include "fhiclcpp/ParameterSetRegistry.h" 

#include <cstddef> // std::ptrdiff_t
#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <typeinfo>
#include <memory> // std::unique_ptr<>

#include "TTree.h"
#include "TTimeStamp.h"

//#include "ubana/SinglePhotonAnalysis/SinglePhoton_module.h"

//#ifdef __MAKECINT__
#ifdef __CLING__
#pragma link C++ class std::vector < std::vector<Float_t> >+; 
#pragma link C++ class std::vector < std::vector< std::vector<Float_t> > >+; 
#endif

const int kMaxTracks=20;
const int kMaxShowers=20;
const int kMaxParticles=10;
const int kMaxMerge=10;

namespace Kaon_Analyzer
{
  //namespace microboone{

  using namespace std;
  
   class CCKaonAnalyzer : public art::EDAnalyzer { 
   public: 

     explicit CCKaonAnalyzer(fhicl::ParameterSet const& pset); 
     virtual ~CCKaonAnalyzer(); 

     void endSubRun(const art::SubRun &subrun); 
     void beginJob(); 
     void filter(const art::Event &evt); 
     void analyze(const art::Event& evt); 
     void reset(); 

     //void fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, int track_i=-1, int daughter_i=-1, bool wTrackRebuilder=false); 
     void fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, art::Ptr<recob::Track>& ptrack, int track_i=-1, int daughter_i=-1, bool wTrackRebuilder=false); 
     void fillPID(const std::vector<art::Ptr<anab::ParticleID>> &trackPID, double angle_y, int track_i=-1, int daughter_i=-1, bool wTrackRebuilder=false);   
     void fillTrueMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_recoobj,
                           art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                           int track_i=-1,
                           int daughter_i=-1,
			   bool isTrack=false,
			   bool wTrackRebuilder=false);

     void mergeChecker(std::vector<art::Ptr<recob::Hit>>& hits_from_recoobj,
		       art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
		       int track_i,
		       int daughter_i,
		       bool isTrack=false,
		       bool wTrackRebuilder=false);
       
     double length(const simb::MCParticle& part, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi); 

     bool isInsideVolume(string volume, double x, double y, double z); 
     bool isInsideVolume(string volume, const TVector3& v) { 
       return isInsideVolume(volume, v.X(), v.Y(), v.Z()); 
      } 

     /*
    typedef art::ValidHandle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
    typedef std::vector< art::Ptr<recob::Track> > TrackVector;
    typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;
    typedef std::map< size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;
     */

  //void 

  private:

     //basic event info
     unsigned int fEventID;
     int run,subrun,event;
     
     //output trees
     TTree* fEventTree;
     //TTree * fMetaTree;


     ///////////////////////////
     //   Truth level info    //
     ///////////////////////////

     //generator truth info

     int fNMCTruths=0;

     std::string fMode; //interaction mode
     std::string fCCNC; //charged current/neutral current

     bool fIsKaon = false;
     bool fIsKaonPlus = false;
     bool fIsKaonMinus = false;
     bool fIsHyperon = false;

     //flags for K+ decay modes
     bool fIsKaonPlus_NuMuP = false;
     bool fIsKaonPlus_PiPPiN = false;
     bool fIsKaonPlus_2PiPPiM = false;
     bool fIsKaonPlus_ENuE = false;
     bool fIsKaonPlus_2PiNPiP = false;
     bool fIsKaonPlus_Others = false;
     bool fIsInelastic_KaonPlus = false;
     bool fIsInelastic_Others = false;

     bool fIsSignal_NuMuP = false;
     bool fIsSignal_PiPPiN = false;

     //neutrino
     std::vector<SimParticle> fNeutrino;

     //K+ from primary vertex
     std::vector<SimParticle> fKaonPlus;

     //K- from primary vertex
     std::vector<SimParticle> fKaonMinus;

     //kaons produced at primary vtx
     std::vector<SimParticle> fKaonOthers;

     //lepton produced at primary vtx
     std::vector<SimParticle> fLepton;

     //hyperon produced at primary vtx
     std::vector<SimParticle> fPrimaryHyperon;

     //nucleons produced at primary vtx
     std::vector<SimParticle> fPrimaryNucleon;

     //pions produced at primary vtx
     std::vector<SimParticle> fPrimaryPion;

     //vertex information
     TVector3 fTruePrimaryVertex;

     //g4 truth info
     TVector3 fDecayVertex;

     std::vector<SimParticle> fDecay; //kaon decay products
     std::vector<SimParticle> fHyperonDecay; //hyperon decay products

     //data storage (should not be written to trees)

     //used by G4 to track particles
     std::vector<int>daughter_IDs; //ids of semistable K+ decay products
     std::vector<int>primary_IDs; //ids of particles produced at primary vertex
     std::vector<int>KaonPlus_daughter_IDs; //ids of K+ decay products
     std::vector<int>KaonMinus_daughter_IDs; //ids of K- decay products
     std::vector<int>KaonOthers_daughter_IDs;
     std::vector<int>Hyperon_daughter_IDs;
     std::vector<int>KaonPlus_Inelastic_daughter_IDs;

     //create map between particles and their ID's
     std::map<int,art::Ptr<simb::MCParticle>> partByID;
     std::pair<int,art::Ptr<simb::MCParticle>>  part_and_id;

     int mode; //interaction mode for event generator
     int ccnc;



     ///////////////////////////
     //   Reco level info    //
     ///////////////////////////


     bool fSelectedEvent = false; //true if event passes some selection criteria

     TVector3 fRecoPrimaryVertex;

     int fNPrimaryDaughters; //num of primary daughters
     int fNPrimaryTrackDaughters_NuMuP; //num of track like primary daughters
     int fNPrimaryTrackDaughters_PiPPiN; //num of track like primary daughters
     int fNPrimaryShowerDaughters; //num of shower like primary daughters


     //Primary daughters
     std::vector<RecoParticle> fTrackPrimaryDaughters_NuMuP;
     std::vector<RecoParticle> fTrackPrimaryDaughters_PiPPiN;
     std::vector<RecoParticle> fShowerPrimaryDaughters;
     int fMuonIndex=-1; //index of muon candidate in fTrackPrimaryDaughters , -1 if no muon candidate found



     ///////////////////////////
     //  Truth Matching info  //
     ///////////////////////////

     //indices in the reco track vector of the true muon
     //and proton and pion from hyperon decay if the exist

     //will have values of -1 if they do not exist

     int fTrueMuonIndex=-1;
     int fTrueKaonIndex=-1;
     int fTrueDecayMuonIndex=-1;
     int fTrueDecayPionIndex=-1;



     /////////////////////////
     // Metadata for sample //
     /////////////////////////

     //truth level metadata

     int fNEvents; //total events in sample
     int fNEventsInActiveVol; //total events in active volume (at truth level)

     int fNChargedCurrent; //number of cc events
     int fNNeutralCurrent; //number of nc events
     int fNnuMu; //number of numu events
     int fNnue; //number of nue events
     int fNnuMuBar; //number of numubar events
     int fNnueBar; //number of nuebar events

     int fNHyperons;

     int fNSignal; //number of signal events

     int fNGoodReco; //number of signal events with both pion and proton reco'd

     //reco level metadata

     int fNSelectedEvents; //total events passing selection

     int fNSelectedHyperons; //total true hyperon events passing selection

     int fNSelectedSignal; //number of signal events passing selection

     int fNSelectedGoodReco; //number of signal events passing selection

     //selection performance metrics

     double fHyperonEfficiency; //hyperon selection efficiency
     double fHyperonPurity;  //hyperon selection purity
     double fHyperonTruePurity; //hyperon selection efficiency x purity
     doublefHyperonEfficiencyTimesPurity; //hyperon selection efficiency x purity
     doublefHyperonEfficiencyTimesTruePurity; //hyperon selection efficiency x true purity

     doublefSignalEfficiency=0; //hyperon selection efficiency
     doublefSignalPurity=0;  //hyperon selection purity
     doublefSignalTruePurity=0; //hyperon purity after converting from enriched sample to real sample
     doublefSignalEfficiencyTimesPurity=0; //hyperon selection efficiency x purity
     doublefSignalEfficiencyTimesTruePurity=0; //hyperon selection efficiency x true purity

     doublefGoodRecoEfficiency=0; //hyperon selection efficiency
     doublefGoodRecoPurity=0;  //hyperon selection purity
     doublefGoodRecoTruePurity=0; //hyperon purity after converting from enriched sample to real sample
     doublefGoodRecoEfficiencyTimesPurity=0; //hyperon selection efficiency x purity
     doublefGoodRecoEfficiencyTimesTruePurity=0; //hyperon selection efficiency x true purity

     double fPOT=0; //total POT of the sample


     // PID
    searchingfornues::LLRPID llr_pid_calculator;
    searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;
    searchingfornuesk::LLRPIDK llr_pid_calculator_k;
    searchingfornuesk::KaonProtonLookUpParameters protonmuon_parameters_k;



    //////////////////////////
    //   FHICL PARAMETERS   //
    //////////////////////////

    bool isMC;
    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fTrackModuleLabel;
    std::string fShowerModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fPIDLabel;
    std::string fHitTruthAssns;
    std::string fHitTrackAssns;
    std::string fHitShowerAssns;
    std::string m_pfp_producer;
    std::string fPFParticleLabel;
    std::string fSpacePointproducer;
    std::string fRebuiltTrackModuleLabel;
    std::string fRebuiltCalorimetryModuleLabel;
    std::string fRebuiltPIDLabel;
    std::string fRebuiltHitTrackAssns;


    //// Tree for the POT subrun info
    TTree *fSubrunTree;
    uint m_run, m_subrun;
    float m_pot;
    //float event_weight;

    //TString output_pdf;
    TCanvas * c = new TCanvas("c", "c", 800, 600); 

    std::map<std::string, std::vector<float>> event_weight;

    std::vector<std::string> evtwgt_funcname;          // the name of the functions used
    std::vector<std::vector<float>> evtwgt_weight;    // the weights (a vector for each function used)
    std::vector<int> evtwgt_nweight;                   // number of weights for each function
    //Int_t evtwgt_nfunc;                                // number of functions used

    bool m_is_verbose;
    std::string m_pandoraLabel;

  };






    DEFINE_ART_MODULE(CCKaonAnalyzer)

      } // namespace lar_pandora

#endif
//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

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

#include "PID/LLR_PID.h"
#include "PID/LLRPID_proton_muon_lookup.h"
#include "PID/LLR_PID_K.h"
#include "PID/LLRPID_kaon_proton_lookup.h"

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

     void fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, int track_i=-1, int daughter_i=-1, bool wTrackRebuilder=false); 
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


    typedef art::ValidHandle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
    typedef std::vector< art::Ptr<recob::Track> > TrackVector;
    typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;
    typedef std::map< size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;

  //void 

  private:

    searchingfornues::LLRPID llr_pid_calculator;
    searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;

    searchingfornuesk::LLRPIDK llr_pid_calculator_k;
    searchingfornuesk::KaonProtonLookUpParameters protonmuon_parameters_k;

    TTree* fEventTree;
    
    Int_t   run;                  
    Int_t   subrun;               
    Int_t   event;

     Int_t n_prip;
     std::vector<int> prip{vector<int>(20,-999)}; 
     Int_t n_prip_k_dau;
     std::vector<int> prip_k_dau{vector<int>(20,-999)};
     Bool_t IsKaon;
     Bool_t IsSingleKaon;
     Bool_t IsAssociatedKaon;
     Bool_t IsMuBR; 
     Bool_t IsPiBR; 
     Bool_t flag_prim_k;

     Bool_t flg_reco_k  = false;
     Bool_t flg_reco_mu  = false;
     Bool_t flg_reco_pi  = false;

     std::string Process;
     std::string EndProcess;

    Float_t true_nu_energy;
    Int_t   true_nu_pdg;
    Int_t   true_nu_mode;
    Int_t   true_nu_ccnc;
    Float_t true_nu_vtx_x;
    Float_t true_nu_vtx_y;
    Float_t true_nu_vtx_z;
    Bool_t  true_nu_vtx_inTPC;
    Bool_t  true_nu_vtx_in5cmTPC;
    Bool_t  true_nu_vtx_inCCInclusiveTPC;

    Int_t   true_lepton_pdg;
    Float_t true_lepton_p;
    Float_t true_lepton_ke;
    Float_t true_lepton_theta;
    Float_t true_lepton_costheta;
    Float_t true_lepton_phi;

    Int_t   true_nkaons;
    Int_t   true_nprimary;
    Int_t   true_nhyperons;

    Float_t true_kaon_length;
    Float_t true_kaon_p;
    Float_t true_kaon_ke;
    Float_t true_kaon_theta;
    Float_t true_kaon_costheta;
    Float_t true_kaon_phi;
    Float_t true_kaon_ccmuon_angle;
    Float_t true_kaon_ccmuon_cosangle;
    Int_t   true_kaon_end_process;
    Float_t true_kaon_end_ke;
    Float_t true_kaon_end_x;
    Float_t true_kaon_end_y;
    Float_t true_kaon_end_z;
    Bool_t  true_kaon_end_inTPC;
    Bool_t  true_kaon_end_in5cmTPC;
    Bool_t  true_kaon_end_inCCInclusiveTPC;

    Float_t true_dau_muon_length;
    Float_t true_dau_muon_p;
    Float_t true_dau_muon_ke;
    Float_t true_dau_muon_theta;
    Float_t true_dau_muon_costheta;
    Float_t true_dau_muon_phi;
    Float_t true_dau_muon_ccmuon_angle;
    Float_t true_dau_muon_ccmuon_cosangle;
    Int_t   true_dau_muon_end_process;
    Float_t true_dau_muon_end_ke;
    Float_t true_dau_muon_start_x;
    Float_t true_dau_muon_start_y;
    Float_t true_dau_muon_start_z;
    Float_t true_dau_muon_end_x;
    Float_t true_dau_muon_end_y;
    Float_t true_dau_muon_end_z;
    Bool_t  true_dau_muon_end_inTPC;
    Bool_t  true_dau_muon_end_in5cmTPC;
    Bool_t  true_dau_muon_end_inCCInclusiveTPC;

    Float_t true_dau_pip_length;
    Float_t true_dau_pip_p;
    Float_t true_dau_pip_ke;
    Float_t true_dau_pip_theta;
    Float_t true_dau_pip_costheta;
    Float_t true_dau_pip_phi;
    Float_t true_dau_pip_ccmuon_angle;
    Float_t true_dau_pip_ccmuon_cosangle;
    Int_t   true_dau_pip_end_process;
    Float_t true_dau_pip_end_ke;
    Float_t true_dau_pip_start_x;
    Float_t true_dau_pip_start_y;
    Float_t true_dau_pip_start_z;
    Float_t true_dau_pip_end_x;
    Float_t true_dau_pip_end_y;
    Float_t true_dau_pip_end_z;
    Bool_t  true_dau_pip_end_inTPC;
    Bool_t  true_dau_pip_end_in5cmTPC;
    Bool_t  true_dau_pip_end_inCCInclusiveTPC;

    Float_t true_dau_pin_length;
    Float_t true_dau_pin_p;
    Float_t true_dau_pin_ke;
    Float_t true_dau_pin_theta;
    Float_t true_dau_pin_costheta;
    Float_t true_dau_pin_phi;
    Float_t true_dau_pin_ccmuon_angle;
    Float_t true_dau_pin_ccmuon_cosangle;
    Int_t   true_dau_pin_end_process;
    Float_t true_dau_pin_end_ke;
    Float_t true_dau_pin_start_x;
    Float_t true_dau_pin_start_y;
    Float_t true_dau_pin_start_z;
    Float_t true_dau_pin_end_x;
    Float_t true_dau_pin_end_y;
    Float_t true_dau_pin_end_z;
    Bool_t  true_dau_pin_end_inTPC;
    Bool_t  true_dau_pin_end_in5cmTPC;
    Bool_t  true_dau_pin_end_inCCInclusiveTPC;

    //should be [kMaxParticles] -> [kMaxTracks][kMaxParticles]
    Float_t cheat_peak_pdg[kMaxTracks][kMaxParticles];
    Float_t cheat_peak_theta[kMaxTracks][kMaxParticles];
    Float_t cheat_peak_phi[kMaxTracks][kMaxParticles];
    Float_t best_peak_theta[kMaxTracks][kMaxParticles];
    Float_t best_peak_phi[kMaxTracks][kMaxParticles];
    Float_t best_peak_trkln[kMaxTracks][kMaxParticles];
    Float_t cheat_pip_trkln;
    Float_t cheat_mup_trkln;

    Float_t true_length[kMaxParticles];
    Float_t true_p[kMaxParticles];
    Float_t true_ke[kMaxParticles];
    Float_t true_theta[kMaxParticles];
    Float_t true_costheta[kMaxParticles];
    Float_t true_phi[kMaxParticles];
    Float_t true_ccmuon_angle[kMaxParticles];
    Float_t true_ccmuon_cosangle[kMaxParticles];
    Int_t   true_end_process[kMaxParticles];
    Float_t true_end_ke[kMaxParticles];
    Float_t true_end_x[kMaxParticles];
    Float_t true_end_y[kMaxParticles];
    Float_t true_end_z[kMaxParticles];
    Bool_t  true_end_inTPC[kMaxParticles];
    Bool_t  true_end_in5cmTPC[kMaxParticles];
    Bool_t  true_end_inCCInclusiveTPC[kMaxParticles];

    Int_t   true_kaon_ndaughters;
    Int_t   true_kaon_ndaughters_decay;
    Int_t   true_kaon_ndaughters_inelastic;
    Int_t   true_kaon_ndecmup;
    Int_t   true_kaon_ndecpip;
    Int_t   true_kaon_ninekap;
    Int_t   true_kaon_ninepip;
    Int_t   true_kaon_ninepro;
    Float_t true_kaon_daughter_length;
    Float_t true_kaon_daughter_p;
    Float_t true_kaon_daughter_ke;
    Float_t true_kaon_daughter_theta;
    Float_t true_kaon_daughter_costheta;
    Float_t true_kaon_daughter_angle;
    Float_t true_kaon_daughter_cosangle;
    Int_t   true_kaon_daughter_pdg;
    Float_t true_kaon_daughter_end_x;
    Float_t true_kaon_daughter_end_y;
    Float_t true_kaon_daughter_end_z;
    Bool_t  true_kaon_daughter_end_inTPC;
    Bool_t  true_kaon_daughter_end_in5cmTPC;
    Bool_t  true_kaon_daughter_end_inCCInclusiveTPC;


    Int_t   true_nkaons_anti;
    Float_t true_kaon_length_anti;
    Float_t true_kaon_p_anti;
    Float_t true_kaon_ke_anti;
    Float_t true_kaon_theta_anti;
    Float_t true_kaon_costheta_anti;
    Float_t true_kaon_phi_anti;
    Float_t true_kaon_ccmuon_angle_anti;
    Float_t true_kaon_ccmuon_cosangle_anti;
    Int_t   true_kaon_end_process_anti;
    Float_t true_kaon_end_ke_anti;
    Float_t true_kaon_end_x_anti;
    Float_t true_kaon_end_y_anti;
    Float_t true_kaon_end_z_anti;
    Bool_t  true_kaon_end_inTPC_anti;
    Bool_t  true_kaon_end_in5cmTPC_anti;
    Bool_t  true_kaon_end_inCCInclusiveTPC_anti;

    Int_t   true_kaon_ndaughters_anti;
    Int_t   true_kaon_ndaughters_decay_anti;
    Int_t   true_kaon_ndaughters_inelastic_anti;
    Int_t   true_kaon_ndecmup_anti;
    Int_t   true_kaon_ndecpip_anti;
    Int_t   true_kaon_ninekap_anti;
    Int_t   true_kaon_ninepip_anti;
    Int_t   true_kaon_ninepro_anti;
    Float_t true_kaon_daughter_length_anti;
    Float_t true_kaon_daughter_p_anti;
    Float_t true_kaon_daughter_ke_anti;
    Float_t true_kaon_daughter_theta_anti;
    Float_t true_kaon_daughter_costheta_anti;
    Float_t true_kaon_daughter_angle_anti;
    Float_t true_kaon_daughter_cosangle_anti;
    Int_t   true_kaon_daughter_pdg_anti;
    Float_t true_kaon_daughter_end_x_anti;
    Float_t true_kaon_daughter_end_y_anti;
    Float_t true_kaon_daughter_end_z_anti;
    Bool_t  true_kaon_daughter_end_inTPC_anti;
    Bool_t  true_kaon_daughter_end_in5cmTPC_anti;
    Bool_t  true_kaon_daughter_end_inCCInclusiveTPC_anti;

    Bool_t  reco_nu_cc_filter;

    Float_t reco_nu_vtx_x;
    Float_t reco_nu_vtx_y;
    Float_t reco_nu_vtx_z;
    Bool_t  reco_nu_vtx_inTPC;
    Bool_t  reco_nu_vtx_in5cmTPC;
    Bool_t  reco_nu_vtx_inCCInclusiveTPC;
    Int_t   reco_nu_ndaughters;
    Int_t   reco_nu_cc_nmue;

    Float_t reco_ccmu_vtx_x;
    Float_t reco_ccmu_vtx_y;
    Float_t reco_ccmu_vtx_z;
    Bool_t  reco_ccmu_vtx_inTPC;
    Bool_t  reco_ccmu_vtx_in5cmTPC;
    Bool_t  reco_ccmu_vtx_inCCInclusiveTPC;
    Int_t   reco_ccmu_true_pdg;
    Int_t   reco_ccmu_true_origin;
    Bool_t  reco_ccmu_true_primary;
    Bool_t  reco_ccmu_true_end_inTPC;
    Bool_t  reco_ccmu_true_end_in5cmTPC;
    Bool_t  reco_ccmu_true_end_inCCInclusiveTPC;
    Float_t reco_ccmu_true_length;

     Int_t   reco_ntracks; 
     Int_t   reco_nshowers;
     Int_t   nTracks;
     Int_t   nShowers; 

     Float_t true_kaon_start_x;
     Float_t true_kaon_start_y;
     Float_t true_kaon_start_z;

     Float_t true_kaon_daughter_start_x;
     Float_t true_kaon_daughter_start_y;
     Float_t true_kaon_daughter_start_z;
     
     Float_t reco_track_start_x[kMaxTracks];//
     Float_t reco_track_start_y[kMaxTracks];
    Float_t reco_track_start_z[kMaxTracks];
    Float_t reco_track_end_x[kMaxTracks];
    Float_t reco_track_end_y[kMaxTracks];
    Float_t reco_track_end_z[kMaxTracks];


    Float_t reco_track_daughter_start_x[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_start_y[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_start_z[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_end_x[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_end_y[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_end_z[kMaxTracks][kMaxTracks];

    Float_t reco_track_daughter_start_x_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_start_y_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_start_z_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_end_x_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_end_y_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_end_z_rebuild[kMaxTracks][kMaxTracks];


    Float_t reco_track_distance[kMaxTracks];

    Float_t   reco_track_kin0[kMaxTracks];
    Float_t   reco_track_kin1[kMaxTracks];
    Float_t   reco_track_kin2[kMaxTracks];

    Float_t reco_track_length[kMaxTracks];
    Float_t reco_track_theta[kMaxTracks];
    Float_t reco_track_phi[kMaxTracks];
    Bool_t  reco_track_dir[kMaxTracks];

    Float_t reco_track_P_vtx[kMaxTracks];
    Float_t reco_track_P_str[kMaxTracks];
    Float_t reco_track_P_end[kMaxTracks];

    Float_t reco_track_chi2ka_pl0[kMaxTracks];
    Float_t reco_track_chi2pr_pl0[kMaxTracks];
    Float_t reco_track_chi2pi_pl0[kMaxTracks];
    Float_t reco_track_chi2mu_pl0[kMaxTracks];
    Float_t reco_track_chi2ka_pl1[kMaxTracks];
    Float_t reco_track_chi2pr_pl1[kMaxTracks];
    Float_t reco_track_chi2pi_pl1[kMaxTracks];
    Float_t reco_track_chi2mu_pl1[kMaxTracks];
    Float_t reco_track_chi2ka_pl2[kMaxTracks];
    Float_t reco_track_chi2pr_pl2[kMaxTracks];
    Float_t reco_track_chi2pi_pl2[kMaxTracks];
    Float_t reco_track_chi2mu_pl2[kMaxTracks];
    Float_t reco_track_chi2ka_3pl[kMaxTracks];
    Float_t reco_track_chi2pr_3pl[kMaxTracks];
    Float_t reco_track_chi2pi_3pl[kMaxTracks];
    Float_t reco_track_chi2mu_3pl[kMaxTracks];
    Float_t reco_track_likepr_3pl[kMaxTracks];


    Float_t reco_track_Bragg_fwd_ka_pl0[kMaxTracks];
    Float_t reco_track_Bragg_fwd_pr_pl0[kMaxTracks];
    Float_t reco_track_Bragg_fwd_pi_pl0[kMaxTracks];
    Float_t reco_track_Bragg_fwd_mu_pl0[kMaxTracks];
    Float_t reco_track_Bragg_fwd_ka_pl1[kMaxTracks];
    Float_t reco_track_Bragg_fwd_pr_pl1[kMaxTracks];
    Float_t reco_track_Bragg_fwd_pi_pl1[kMaxTracks];
    Float_t reco_track_Bragg_fwd_mu_pl1[kMaxTracks];
    Float_t reco_track_Bragg_fwd_ka_pl2[kMaxTracks];
    Float_t reco_track_Bragg_fwd_pr_pl2[kMaxTracks];
    Float_t reco_track_Bragg_fwd_pi_pl2[kMaxTracks];
    Float_t reco_track_Bragg_fwd_mu_pl2[kMaxTracks];

    Float_t reco_track_MIP_pl0[kMaxTracks];
    Float_t reco_track_MIP_pl1[kMaxTracks];
    Float_t reco_track_MIP_pl2[kMaxTracks];


    Float_t reco_track_daughter_Bragg_fwd_ka_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pr_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pi_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_mu_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_ka_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pr_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pi_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_mu_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_ka_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pr_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pi_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_mu_pl2[kMaxTracks][kMaxTracks];

    Float_t reco_track_daughter_MIP_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_MIP_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_MIP_pl2[kMaxTracks][kMaxTracks];


    Float_t reco_track_llrpid_3pl[kMaxTracks];
    Float_t reco_track_total_llrpid_3pl[kMaxTracks];
    Float_t reco_track_llrpid_k_3pl[kMaxTracks];
    Bool_t  reco_track_vtx_inTPC[kMaxTracks];
    Bool_t  reco_track_vtx_in5cmTPC[kMaxTracks];
    Bool_t  reco_track_vtx_inCCInclusiveTPC[kMaxTracks];
    Bool_t  reco_track_end_inTPC[kMaxTracks];
    Bool_t  reco_track_end_in5cmTPC[kMaxTracks];
    Bool_t  reco_track_end_inCCInclusiveTPC[kMaxTracks];
    Int_t   reco_track_true_pdg[kMaxTracks];
    Int_t   reco_track_true_origin[kMaxTracks];
    Bool_t  reco_track_true_primary[kMaxTracks];
    Bool_t  reco_track_true_end_inTPC[kMaxTracks];
    Bool_t  reco_track_true_end_in5cmTPC[kMaxTracks];
    Bool_t  reco_track_true_end_inCCInclusiveTPC[kMaxTracks];
    Float_t reco_track_true_length[kMaxTracks];

    Float_t reco_track_match_e[kMaxTracks][kMaxMerge];
    Float_t reco_track_match_hit[kMaxTracks][kMaxMerge];
    Int_t reco_track_match_epdg[kMaxTracks][kMaxMerge];
    Int_t reco_track_match_hitpdg[kMaxTracks][kMaxMerge];

    Float_t reco_track_daughter_match_e[kMaxTracks][kMaxTracks][kMaxMerge];
    Float_t reco_track_daughter_match_hit[kMaxTracks][kMaxTracks][kMaxMerge];
    Int_t reco_track_daughter_match_epdg[kMaxTracks][kMaxTracks][kMaxMerge];
    Int_t reco_track_daughter_match_hitpdg[kMaxTracks][kMaxTracks][kMaxMerge];

    Float_t reco_track_daughter_match_e_rebuild[kMaxTracks][kMaxTracks][kMaxMerge];
    Float_t reco_track_daughter_match_hit_rebuild[kMaxTracks][kMaxTracks][kMaxMerge];
    Int_t reco_track_daughter_match_epdg_rebuild[kMaxTracks][kMaxTracks][kMaxMerge];
    Int_t reco_track_daughter_match_hitpdg_rebuild[kMaxTracks][kMaxTracks][kMaxMerge];

    Float_t reco_track_daughter_shower_match_e[kMaxTracks][kMaxShowers][kMaxMerge];
    Float_t reco_track_daughter_shower_match_hit[kMaxTracks][kMaxShowers][kMaxMerge];
    Int_t reco_track_daughter_shower_match_epdg[kMaxTracks][kMaxShowers][kMaxMerge];
    Int_t reco_track_daughter_shower_match_hitpdg[kMaxTracks][kMaxShowers][kMaxMerge];


    Int_t   reco_track_ndaughters[kMaxTracks];
    Float_t reco_track_daughter_distance[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_vtx_distance[kMaxTracks][kMaxTracks];
    Float_t reco_angle_track_daughter[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_nhits0[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_nhits1[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_nhits2[kMaxTracks][kMaxTracks];


    Int_t   reco_track_ndaughters_rebuild[kMaxTracks];
    Float_t reco_track_daughter_distance_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_vtx_distance_rebuild[kMaxTracks][kMaxShowers]; 
    Int_t   reco_track_daughter_true_pdg_rebuild[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_true_origin_rebuild[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_true_primary_rebuild[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_true_end_inTPC_rebuild[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_true_end_in5cmTPC_rebuild[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_true_end_inCCInclusiveTPC_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_true_length_rebuild[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_true_mother_rebuild[kMaxTracks][kMaxTracks];

    Float_t reco_track_daughter_length_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_theta_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_phi_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2ka_pl0_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pr_pl0_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi_pl0_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu_pl0_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2ka_pl1_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pr_pl1_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi_pl1_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu_pl1_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2ka_pl2_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pr_pl2_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi_pl2_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu_pl2_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2ka_3pl_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pr_3pl_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi_3pl_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu_3pl_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_likepr_3pl_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_llrpid_3pl_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_llrpid_k_3pl_rebuild[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_vtx_inTPC_rebuild[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_vtx_in5cmTPC_rebuild[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_vtx_inCCInclusiveTPC_rebuild[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_end_inTPC_rebuild[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_end_in5cmTPC_rebuild[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_end_inCCInclusiveTPC_rebuild[kMaxTracks][kMaxTracks];

    Float_t reco_track_daughter_Bragg_fwd_ka_pl0_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pr_pl0_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pi_pl0_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_mu_pl0_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_ka_pl1_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pr_pl1_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pi_pl1_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_mu_pl1_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_ka_pl2_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pr_pl2_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pi_pl2_rebuild[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_mu_pl2_rebuild[kMaxTracks][kMaxTracks];



     Int_t   reco_shower_ndaughters[kMaxTracks]; 
     Float_t reco_track_daughter_distance_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_angle_track_daughter_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_angle_daughter_track_daughter_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_length_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_theta_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_phi_sh[kMaxTracks][kMaxShowers]; 
     Int_t reco_track_daughter_true_pdg_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_open_angle_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_dedx_pl0_sh[kMaxTracks][kMaxShowers];
     Float_t reco_track_daughter_dedx_pl1_sh[kMaxTracks][kMaxShowers];
     Float_t reco_track_daughter_dedx_pl2_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_energy_pl0_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_energy_pl1_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_energy_pl2_sh[kMaxTracks][kMaxShowers]; 
 
     Bool_t  reco_track_daughter_vtx_inTPC_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_vtx_in5cmTPC_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_vtx_inCCInclusiveTPC_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_end_inTPC_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_end_in5cmTPC_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_end_inCCInclusiveTPC_sh[kMaxTracks][kMaxShowers]; 

     Int_t   reco_track_true_pdg_sh[kMaxShowers]; 
     Int_t   reco_track_true_origin_sh[kMaxShowers];
     Bool_t  reco_track_true_primary_sh[kMaxShowers]; 
     Bool_t  reco_track_true_end_inTPC_sh[kMaxShowers];
     Bool_t  reco_track_true_end_in5cmTPC_sh[kMaxShowers];
     Bool_t  reco_track_true_end_inCCInclusiveTPC_sh[kMaxShowers]; 
     Float_t reco_track_true_length_sh[kMaxShowers]; 

     Int_t   reco_track_daughter_true_origin_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_true_primary_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_true_end_inTPC_sh[kMaxTracks][kMaxShowers];
     Bool_t  reco_track_daughter_true_end_in5cmTPC_sh[kMaxTracks][kMaxShowers];
     Bool_t  reco_track_daughter_true_end_inCCInclusiveTPC_sh[kMaxTracks][kMaxShowers];
     Float_t reco_track_daughter_true_length_sh[kMaxTracks][kMaxShowers];
     Int_t   reco_track_daughter_true_mother_sh[kMaxTracks][kMaxShowers];


    Float_t reco_track_daughter_length[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_theta[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_phi[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2ka_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pr_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2ka_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pr_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2ka_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pr_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2ka_3pl[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pr_3pl[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi_3pl[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu_3pl[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_likepr_3pl[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_llrpid_3pl[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_llrpid_k_3pl[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_vtx_inTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_vtx_in5cmTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_vtx_inCCInclusiveTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_end_inTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_end_in5cmTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_end_inCCInclusiveTPC[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_true_pdg[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_true_origin[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_true_primary[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_true_end_inTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_true_end_in5cmTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_true_end_inCCInclusiveTPC[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_true_length[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_true_mother[kMaxTracks][kMaxTracks];


    Int_t   k_can_trkid; //track id of reco kaon
    Int_t   mu_can_trkid; //track id of reco muon
    Float_t k_mu_can_dis; //distance between kaon and muon end
    Float_t k_mu_open_angle; //angle between kaon and muon end
    Float_t k_vtx_dis; //distance between kaon and muon start


    Int_t cut_1;
    Int_t cut_2;
    Int_t cut_3;
    Int_t cut_4;
    Int_t cut_5;
    Int_t cut_6;
    Int_t cut_7;
    Int_t cut_8;
    Int_t cut_9;
    Int_t cut_10;
    Int_t cut_11;
    Int_t cut_12;
    Int_t cut_13;

    //const trkf::TrackMomentumCalculator trkmom;

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
    bool isMC;

    std::string fRebuiltTrackModuleLabel;
    std::string fRebuiltCalorimetryModuleLabel;
    std::string fRebuiltPIDLabel;
    std::string fRebuiltHitTrackAssns;

  //    std::string m_pandoraLabel;
  //   std::string m_is_verbose;  

    //std::vector<Float_t> adEdx;
  std::vector<Float_t> dv0;
  std::vector<Float_t> rv0;
  std::vector<Float_t> dv1;
  std::vector<Float_t> rv1;
  std::vector<Float_t> dv2;
  std::vector<Float_t> rv2;

    //// Tree for the POT subrun info
    TTree *fSubrunTree;
    uint m_run, m_subrun;
    float m_pot;
    //float event_weight;

    //TString output_pdf;
    TCanvas * c = new TCanvas("c", "c", 800, 600); 
    //bool isFirstpdf = true;
    //bool isFirstpdf_cheat = true;
    //int i_pdf=0;
    //int i_pdf_cheat=0;
    //const int num_pdf=3;
    //const int num_pdf_cheat=3;

    //std::map<int, int> angular_distribution_map_track;
    //std::map<int, int> angular_distribution_map_shower;
    //std::map<int, std::map<int, int>> angular_distribution_mrgid_map;
    //vector<TH1D*> h_angular_distribution_pfparticle;
    //vector<bool> v_trk_flg;

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

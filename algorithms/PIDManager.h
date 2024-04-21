#ifndef _PIDManager_h_
#define _PIDManager_h_

#include <string>
#include <vector>
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "cetlib_except/exception.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "../headers/LLR_PID.h"
#include "../headers/LLRPID_proton_muon_lookup.h"
#include "../headers/LLR_PID_K.h"
#include "../headers/LLRPID_kaon_proton_lookup.h"
#include "../headers/LLRPID_correction_lookup.h"
#include "TVector3.h"
#include <TH3.h>
#include <TFile.h>

namespace hyperon {

  //Define wire plane angles
  const double plane0_wireangle = 30*6.28/360.0;
  const double plane1_wireangle = -30*6.28/360.0;
  const double plane2_wireangle = 90*6.28/360.0;

  struct PIDStore {

    double Weight_Plane0 = 0;
    double Weight_Plane1 = 0;
    double Weight_Plane2 = 0;

    double MeandEdX_Plane0 = -1;
    double MeandEdX_Plane1 = -1;
    double MeandEdX_Plane2 = -1;
    double MeandEdX_3Plane = -1;

    double LLR;
    double LLR_Kaon;
    double LLR_Kaon_Partial; 

    double Bragg_Kaon_Plane0;
    double Bragg_Kaon_Plane1;
    double Bragg_Kaon_Plane2;    
    double Bragg_Kaon_3Plane;

    std::vector<float> dEdX_Corrected_Plane0;
    std::vector<float> dEdX_Corrected_Plane1;
    std::vector<float> dEdX_Corrected_Plane2;

    double LLR_SigmaKaon;
    double LLR_SigmaProton;
    double LLR_SigmaMuon;

  };

  class PIDManager {

  public: 

    PIDManager(const fhicl::ParameterSet& p);

    double GetMeandEdX(art::Ptr<anab::Calorimetry> calo) const;
    void ThreePlaneMeandEdX(art::Ptr<recob::Track> track,std::vector<art::Ptr<anab::Calorimetry>> calo_v,PIDStore& store) const;
    void LLRPID(std::vector<art::Ptr<anab::Calorimetry>> calo_v,PIDStore& store);
    void BraggPID(art::Ptr<recob::Track> track,std::vector<anab::sParticleIDAlgScores> algscores_v,PIDStore& store) const;
    double GetGenericLLRPID(std::vector<art::Ptr<anab::Calorimetry>> calo_v,std::pair<int,int> hypotheses) const;
    PIDStore GetPIDs(art::Ptr<recob::Track> track,std::vector<art::Ptr<anab::Calorimetry>> calo_v,std::vector<anab::sParticleIDAlgScores> algscores_v);   
               
    double PlaneWeight(TVector3 dir,int i_pl) const ;
    double PlaneWeight(art::Ptr<recob::Track> track,int i_pl) const;

    void SetTophatThresh(double thresh){ TophatThresh = thresh; }

  private:

    enum Planes {kPlane0,kPlane1,kPlane2,kInvalid};

    searchingfornues::LLRPID llr_pid_calculator;
    searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;
    searchingfornuesk::LLRPIDK llr_pid_calculator_kaon;
    searchingfornuesk::KaonProtonLookUpParameters kaonproton_parameters;
    searchingfornues::CorrectionLookUpParameters correction_parameters;

    // Miniumum value of sin2(angle between track and wires)
    double TophatThresh = 0.175;
    double ResRangeCutoff=5; 

    const std::string PIDReferenceHists;
    const std::vector<int> pdg_v = {3222,3112,321,2212,13,211};
    std::vector<TH3D*> h_dEdx_Reference_Plane0;
    std::vector<TH3D*> h_dEdx_Reference_Plane1;
    std::vector<TH3D*> h_dEdx_Reference_Plane2;

    void LoadGenericLLRPID();

  };
}

#endif

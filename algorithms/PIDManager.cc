#ifndef _PIDManager_cxx_
#define _PIDManager_cxx_

#include "PIDManager.h"

using namespace hyperon;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

PIDManager::PIDManager(const fhicl::ParameterSet& p) :
  PIDReferenceHists(p.get<std::string>("PIDReferenceHists",""))
{

  llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
  llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
  llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);
  llr_pid_calculator.set_corr_par_binning(0,correction_parameters.parameter_correction_edges_pl_0);
  llr_pid_calculator.set_correction_tables(0,correction_parameters.correction_table_pl_0);
      
  llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
  llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
  llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);
  llr_pid_calculator.set_corr_par_binning(1,correction_parameters.parameter_correction_edges_pl_1);
  llr_pid_calculator.set_correction_tables(1,correction_parameters.correction_table_pl_1);

  llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
  llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
  llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);
  llr_pid_calculator.set_corr_par_binning(2,correction_parameters.parameter_correction_edges_pl_2);
  llr_pid_calculator.set_correction_tables(2,correction_parameters.correction_table_pl_2);

  llr_pid_calculator_kaon.set_dedx_binning(0, kaonproton_parameters.dedx_edges_pl_0);
  llr_pid_calculator_kaon.set_par_binning(0, kaonproton_parameters.parameters_edges_pl_0);
  llr_pid_calculator_kaon.set_lookup_tables(0, kaonproton_parameters.dedx_pdf_pl_0);

  llr_pid_calculator_kaon.set_dedx_binning(1, kaonproton_parameters.dedx_edges_pl_1);
  llr_pid_calculator_kaon.set_par_binning(1, kaonproton_parameters.parameters_edges_pl_1);
  llr_pid_calculator_kaon.set_lookup_tables(1, kaonproton_parameters.dedx_pdf_pl_1);

  llr_pid_calculator_kaon.set_dedx_binning(2, kaonproton_parameters.dedx_edges_pl_2);
  llr_pid_calculator_kaon.set_par_binning(2, kaonproton_parameters.parameters_edges_pl_2);
  llr_pid_calculator_kaon.set_lookup_tables(2, kaonproton_parameters.dedx_pdf_pl_2);
     
  if(PIDReferenceHists != "") LoadGenericLLRPID();   

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double PIDManager::GetMeandEdX(art::Ptr<anab::Calorimetry> calo) const {

  double totalE=0;
  double totalX=0;

  // Sometimes this vector is empty, causes crash below, skip plane if it is
  if(calo->XYZ().size() < 2) return -1;

  for(size_t i_point = 0;i_point < calo->XYZ().size()-1;i_point++){

    anab::Point_t thisPos = calo->XYZ().at(i_point);
    anab::Point_t nextPos = calo->XYZ().at(i_point+1);

    // Step vector
    TVector3 D(thisPos.X()-nextPos.X(),thisPos.X()-nextPos.X(),thisPos.X()-nextPos.X());

    totalX += D.Mag();
    totalE += calo->dEdx().at(i_point)*D.Mag();
  }

  return totalE/totalX;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PIDManager::ThreePlaneMeandEdX(art::Ptr<recob::Track> track,std::vector<art::Ptr<anab::Calorimetry>> calo_v,PIDStore& store) const {

  double TotaldEdX=0;
  double TotalWeight=0;

  for(size_t i_pl=0;i_pl<calo_v.size();i_pl++){

    int plane = calo_v.at(i_pl)->PlaneID().Plane;

    if(plane != 0 && plane != 1 && plane != 2) continue;        

    double dEdX = GetMeandEdX(calo_v.at(i_pl));

    // Catch default fills
    if(dEdX < 0) continue;

    double thisPlaneWeight = PlaneWeight(track,plane);

    /*
      if(plane == 0){
         store.Weight_Plane0 = thisPlaneWeight;       
         store.MeandEdX_Plane0 = dEdX;
         store.dEdX_Plane0 = calo_v.at(i_pl)->dEdx();
         store.ResidualRange_Plane0 = calo_v.at(i_pl)->ResidualRange();
         store.Pitch_Plane0 = calo_v.at(i_pl)->TrkPitchVec();
         store.dEdX_Corrected_Plane0 = llr_pid_calculator.correct_many_hits_one_plane(calo_v.at(i_pl),*track,true,true);
      }
      if(plane == 1){
         store.Weight_Plane1 = thisPlaneWeight;       
         store.MeandEdX_Plane1 = dEdX;
         store.dEdX_Plane1 = calo_v.at(i_pl)->dEdx();
         store.ResidualRange_Plane1 = calo_v.at(i_pl)->ResidualRange();
         store.Pitch_Plane1 = calo_v.at(i_pl)->TrkPitchVec();
         store.dEdX_Corrected_Plane1 = llr_pid_calculator.correct_many_hits_one_plane(calo_v.at(i_pl),*track,true,true);
      }
      if(plane == 2){
         store.Weight_Plane2 = thisPlaneWeight;       
         store.MeandEdX_Plane2 = dEdX;
         store.dEdX_Plane2 = calo_v.at(i_pl)->dEdx();
         store.ResidualRange_Plane2 = calo_v.at(i_pl)->ResidualRange();
         store.Pitch_Plane2 = calo_v.at(i_pl)->TrkPitchVec();
         store.dEdX_Corrected_Plane2 = llr_pid_calculator.correct_many_hits_one_plane(calo_v.at(i_pl),*track,true,true);
      }
    */

    TotaldEdX += dEdX*thisPlaneWeight;
    TotalWeight += thisPlaneWeight;
  }

  if(TotalWeight > 0) store.MeandEdX_3Plane = TotaldEdX/TotalWeight;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PIDManager::LLRPID(std::vector<art::Ptr<anab::Calorimetry>> calo_v,PIDStore& store){

  double this_llr_pid=0;
  double this_llr_pid_score=0;
  double this_llr_pid_kaon=0;
  double this_llr_pid_score_kaon=0;

  double this_llr_pid_kaon_partial = 0;
  double this_llr_pid_score_kaon_partial = 0;

  for(auto const &calo : calo_v){

    auto const &plane = calo->PlaneID().Plane;
    auto const &dedx_values = calo->dEdx();
    auto const &rr = calo->ResidualRange();
    auto const &pitch = calo->TrkPitchVec();
    std::vector<std::vector<float>> par_values;
    par_values.push_back(rr);
    par_values.push_back(pitch);

    // Get parital length PIDs
    std::vector<std::vector<float>> par_values_partial;
    std::vector<float> dedx_values_partial,rr_partial,pitch_partial;      
    if(calo->dEdx().size() != calo->ResidualRange().size() || calo->ResidualRange().size() != calo->TrkPitchVec().size())
      throw cet::exception("SubModuleReco") << "Track calo point list size mismatch" << std::endl;
    for(size_t i_p=0;i_p<calo->dEdx().size();i_p++){
      if(rr.at(i_p) > ResRangeCutoff) continue;
      dedx_values_partial.push_back(calo->dEdx().at(i_p));
      rr_partial.push_back(calo->ResidualRange().at(i_p));
      pitch_partial.push_back(calo->TrkPitchVec().at(i_p));        
    }
    par_values_partial.push_back(rr_partial);
    par_values_partial.push_back(pitch_partial);

    if(calo->ResidualRange().size() == 0) continue;

    float calo_energy = 0;
    for(size_t i=0;i<dedx_values.size();i++)
      calo_energy += dedx_values[i] * pitch[i];

    float llr_pid = llr_pid_calculator.LLR_many_hits_one_plane(dedx_values,par_values,plane);
    float llr_pid_kaon = llr_pid_calculator_kaon.LLR_many_hits_one_plane(dedx_values,par_values,plane);
    this_llr_pid += llr_pid;
    this_llr_pid_kaon += llr_pid_kaon;

    // Partial length calculation
    float calo_energy_partial = 0;
    for(size_t i=0;i<dedx_values_partial.size();i++)
      calo_energy_partial += dedx_values_partial[i] * pitch_partial[i];

    float llr_pid_kaon_partial = llr_pid_calculator_kaon.LLR_many_hits_one_plane(dedx_values_partial,par_values_partial,plane);
    this_llr_pid_kaon_partial += llr_pid_kaon_partial;     
  }

  this_llr_pid_score = atan(this_llr_pid/100.)*2/3.14159266;
  this_llr_pid_score_kaon = atan(this_llr_pid_kaon/100.)*2/3.14159266;
  this_llr_pid_score_kaon_partial = atan(this_llr_pid_kaon_partial/100.)*2/3.14159266;

  store.LLR = this_llr_pid_score;
  store.LLR_Kaon = this_llr_pid_score_kaon;
  store.LLR_Kaon_Partial = this_llr_pid_score_kaon_partial;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PIDManager::BraggPID(art::Ptr<recob::Track> track,std::vector<anab::sParticleIDAlgScores> algscores_v,PIDStore& store) const {

  for(size_t i_algscore=0;i_algscore<algscores_v.size();i_algscore++){
    anab::sParticleIDAlgScores algscore = algscores_v.at(i_algscore);
    if(algscore.fAssumedPdg == 321 && algscore.fAlgName=="BraggPeakLLH" && anab::kTrackDir(algscore.fTrackDir) == anab::kForward){
      if(UBPID::uB_getSinglePlane(algscore.fPlaneMask) == 0) store.Bragg_Kaon_Plane0 = algscore.fValue;
      if(UBPID::uB_getSinglePlane(algscore.fPlaneMask) == 1) store.Bragg_Kaon_Plane1 = algscore.fValue;
      if(UBPID::uB_getSinglePlane(algscore.fPlaneMask) == 2) store.Bragg_Kaon_Plane2 = algscore.fValue;
    }
  }

  store.Bragg_Kaon_3Plane = store.Bragg_Kaon_Plane0*PlaneWeight(track,0) + store.Bragg_Kaon_Plane1*PlaneWeight(track,1) + store.Bragg_Kaon_Plane2*PlaneWeight(track,2);
  store.Bragg_Kaon_3Plane /= (PlaneWeight(track,0) + PlaneWeight(track,1) + PlaneWeight(track,2));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

PIDStore PIDManager::GetPIDs(art::Ptr<recob::Track> track,std::vector<art::Ptr<anab::Calorimetry>> calo_v,std::vector<anab::sParticleIDAlgScores> algscores_v){

  PIDStore theStore;
  ThreePlaneMeandEdX(track,calo_v,theStore);
  LLRPID(calo_v,theStore);
  BraggPID(track,algscores_v,theStore);
  theStore.LLR_SigmaKaon = GetGenericLLRPID(calo_v,std::make_pair(3112,321));
  theStore.LLR_SigmaProton = GetGenericLLRPID(calo_v,std::make_pair(3112,2212));
  theStore.LLR_SigmaMuon = GetGenericLLRPID(calo_v,std::make_pair(3112,13));

  return theStore;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double PIDManager::PlaneWeight(TVector3 dir,int i_pl) const {

  TVector3 trackvec(0, dir.Y(), dir.Z());
  trackvec = trackvec.Unit();
  TVector3 zaxis(0, 0, 1);
  double costhetayz = trackvec.Dot(zaxis);
  double thetayz = TMath::ACos(costhetayz);
  if ((dir.Y() < 0) && (thetayz > 0)) thetayz *= -1;

  double theta_towires = 0;
  if (i_pl == 0) theta_towires = std::min(std::abs(plane0_wireangle - thetayz), std::abs((-1*(6.28-plane0_wireangle) - thetayz)));
  if (i_pl == 1) theta_towires = std::min(std::abs(plane1_wireangle - thetayz), std::abs((-1*(6.28-plane1_wireangle) - thetayz)));
  if (i_pl == 2) theta_towires = std::min(std::abs(plane2_wireangle - thetayz), std::abs((-1*(6.28-plane2_wireangle) - thetayz)));

  double angle_planeweight = sin(theta_towires)*sin(theta_towires);
  if (angle_planeweight < TophatThresh) angle_planeweight = 0;
  if (angle_planeweight != 0) angle_planeweight = 1;

  return angle_planeweight;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double PIDManager::PlaneWeight(art::Ptr<recob::Track> track,int i_pl) const {

  TVector3 trackdir(track->End().x()-track->Start().x(),track->End().y()-track->Start().y(),track->End().z()-track->Start().z());
  return PlaneWeight(trackdir,i_pl);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//void PIDManager::LoadGenericLLRPID(const fhicl::ParameterSet& p){
void PIDManager::LoadGenericLLRPID(){

  TFile* f = TFile::Open(PIDReferenceHists.c_str());

  for(size_t i_pdg=0;i_pdg<pdg_v.size();i_pdg++){
    const std::string pdg = std::to_string(pdg_v.at(i_pdg));
    h_dEdx_Reference_Plane0.push_back(static_cast<TH3D*>(f->Get((pdg+"_Plane0").c_str()))); 
    h_dEdx_Reference_Plane1.push_back(static_cast<TH3D*>(f->Get((pdg+"_Plane1").c_str()))); 
    h_dEdx_Reference_Plane2.push_back(static_cast<TH3D*>(f->Get((pdg+"_Plane2").c_str()))); 
  }
 
  f->Close();

  std::cout << "PIDManager: Loaded PID reference hists from " << PIDReferenceHists << std::endl;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double PIDManager::GetGenericLLRPID(std::vector<art::Ptr<anab::Calorimetry>> calo_v,std::pair<int,int> hypotheses) const {

  if(PIDReferenceHists == "")
    throw cet::exception("PIDManager:") << "Generic LLR PID reference hists not loaded" << std::endl;

  int pdg_index_first = -1;
  int pdg_index_second = -1;
  for(size_t i_pdg=0;i_pdg<pdg_v.size();i_pdg++){
    if(pdg_v.at(i_pdg) == hypotheses.first) pdg_index_first = i_pdg;
    if(pdg_v.at(i_pdg) == hypotheses.second) pdg_index_second = i_pdg;
  }
 
  // TODO: Add a check to catch if hypothesis isn't in reference

  double llr = 0;

  for(auto const &calo : calo_v){

    auto const &plane = calo->PlaneID().Plane;

    TH3D* h_hyp_first = nullptr;
    TH3D* h_hyp_second = nullptr;

    if(plane == kPlane0){
      h_hyp_first = h_dEdx_Reference_Plane0.at(pdg_index_first);
      h_hyp_second = h_dEdx_Reference_Plane0.at(pdg_index_second);
    }
    if(plane == kPlane1){
      h_hyp_first = h_dEdx_Reference_Plane1.at(pdg_index_first);
      h_hyp_second = h_dEdx_Reference_Plane1.at(pdg_index_second);
    }
    if(plane == kPlane2){
      h_hyp_first = h_dEdx_Reference_Plane2.at(pdg_index_first);
      h_hyp_second = h_dEdx_Reference_Plane2.at(pdg_index_second);
    }

    auto const &dedx_values = calo->dEdx();
    auto const &rr = calo->ResidualRange();
    auto const &pitch = calo->TrkPitchVec();
    for(size_t i_p=0;i_p<dedx_values.size();i_p++){ 
      int bin_x = h_hyp_first->GetXaxis()->FindBin(rr.at(i_p));
      int bin_y = h_hyp_first->GetYaxis()->FindBin(dedx_values.at(i_p));
      int bin_z = h_hyp_first->GetZaxis()->FindBin(180/3.142*acos(0.3/pitch.at(i_p)));
      double l_first = h_hyp_first->GetBinContent(bin_x,bin_y,bin_z);
      double l_second = h_hyp_second->GetBinContent(bin_x,bin_y,bin_z);
      //std::cout << l_first <<  "  " << l_second << std::endl;
      if(l_first > 0 && l_second > 0 && !std::isnan(l_first) && !std::isnan(l_second)) llr += log(l_first) - log(l_second);
    } 
  }
    
  //std::cout << "Score: " << 2/3.1415*atan(llr) << std::endl;

  return 2.0/3.1415*atan(llr);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

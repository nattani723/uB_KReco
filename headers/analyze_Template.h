#include "../CCKaonAnalyzerRebuild_module.h"

namespace Kaon_Analzyer
{
    void CCKaonAnalyzerRebuild::ClearTemplates(){
        m_reco_num_templates = 0;
        m_reco_template.clear();
    }

    void CCKaonAnalyzerRebuild::ResizeTemplates(size_t size){
        m_reco_template.resize(size);
    }


    void CCKaonAnalyzerRebuild::CreateTemplateBranches(){
        vertex_tree->Branch("reco_template",&m_reco_template);
    }

    void CCKaonAnalyzerRebuild::AnalyzeTemplates(){
        m_reco_num_templates = 1;
        this->ResizeTemplates(m_reco_num_templates);

        m_reco_template[0]=99;

    }
}

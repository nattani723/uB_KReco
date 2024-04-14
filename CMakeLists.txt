link_libraries( -L$ENV{BOOST_LIB} -lboost_system -lboost_filesystem)
art_make(
          LIBRARY_NAME uboonecode_uboone_Filters

          LIB_LIBRARIES 

		ubraw_RawData
		ubobj_Trigger
            	larcorealg_Geometry
            	larcore_Geometry_Geometry_service
		ubana_BasicTool_GeoAlgo
            	lardata_Utilities
            	larevt_Filters
            	lardataobj_RawData
		lardataobj_RecoBase
            	lardataobj_MCBase
		larcoreobj_SummaryData
		larsim_MCCheater_BackTrackerService_service
		larsim_MCCheater_ParticleInventoryService_service
		larpandora_LArPandoraInterface
            	${PANDORASDK}
            	${PANDORAMONITORING}
            	LArPandoraContent 
		nusimdata_SimulationBase
		${UB_DATA_TYPES}
		${Boost_SERIALIZATION_LIBRARY}
		nutools_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                art_Persistency_Common canvas
                art_Persistency_Provenance canvas
                art_Utilities canvas
                ${MF_MESSAGELOGGER}
                        
                ${FHICLCPP}
		cetlib cetlib_except
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
                ${ROOT_GEOM} 

                        lardata_Utilities_LArFFT_service
                        ${Boost_SYSTEM_LIBRARY}
          MODULE_LIBRARIES 
                        larpandora_LArPandoraInterface
                        hep_concurrency
                        ubobj_UBXSec
			lardataobj_RawData
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larsim_MCCheater_ParticleInventoryService_service
                        larsim_MCCheater_BackTrackerService_service
			larsim_Simulation nutools_ParticleNavigation
			lardataobj_Simulation
			lardataobj_RecoBase
			lardataobj_AnalysisBase
			lardata_RecoObjects larreco_Deprecated  
			larreco_RecoAlg
                        larreco_Calorimetry
			lardata_Utilities
			larevt_Filters
                        lardata_Utilities_LArFFT_service
			SignalShapingServiceMicroBooNE_service
                        larevt_CalibrationDBI_Providers
                        nusimdata_SimulationBase
			${LARRECO_LIB}
			${LARDATA_LIB}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
			${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
			${ROOT_MINUIT}
			${ROOT_MINUIT2}	

			larcoreobj_SummaryData
			larsim_EventWeight_Base
			${ART_FRAMEWORK_SERVICES_BASIC}
			art_Framework_Services_Registry
			${Boost_SYSTEM_LIBRARY}
			CCKaonHeaders
        )
add_subdirectory(headers)
install_headers()
install_fhicl()
install_source()

#include "flashgg/Validation/interface/SimpleTreeMaker.h"


// ******************************************************************************************
// constructors and destructor
//
/*SimpleTreeMaker::SimpleTreeMaker( const edm::ParameterSet &iConfig, TFileDirectory& fs ):
  edm::BasicAnalyzer::BasicAnalyzer(iConfig, fs),
  genParticleToken_( iConfig.getParameter<InputTag>( "genParticleTag" ) ),
  genInfoToken_( iConfig.getParameter<InputTag> ( "generatorInfo" ) ),
  PileUpToken_( iConfig.getParameter<InputTag> ( "PileUpTag" ) ), 
  vertexToken_( iConfig.getParameter<InputTag> ( "VertexTag" ) ),
  diphotonToken_( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ),
  mvaResultToken_( iConfig.getParameter<InputTag> ( "MVAResultTag" ) ),
  inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
  genJetToken_( iConfig.getParameter<InputTag> ( "GenJetTag" ) ),
  electronToken_( iConfig.getParameter<InputTag>( "ElectronTag" ) ),
  muonToken_( iConfig.getParameter<InputTag>( "MuonTag" ) ),
  METToken_( iConfig.getParameter<InputTag> ( "METTag" ) ),
  triggerBitsToken_(  iConfig.getParameter<InputTag>( "triggerBits" ) ),
  rhoToken_( iConfig.getParameter <edm::InputTag>("rhoFixedGridCollection") ) 
{
  jetPtThreshold_ = iConfig.getUntrackedParameter<double>( "jetPtThreshold", 20. );
  bTag_ = iConfig.getUntrackedParameter<string> ( "bTag", "pfCombinedInclusiveSecondaryVertexV2BJetTags" );
  electronPtThreshold_ = iConfig.getUntrackedParameter<double>( "electronPtThreshold", 20. );
  muonPtThreshold_ = iConfig.getUntrackedParameter<double>( "muonPtThreshold", 20. );
  isControlSample_ = iConfig.getUntrackedParameter<bool>( "isControlSample", false );
  lumiWeight_ = iConfig.getUntrackedParameter<double>( "lumiWeight", 1000. ); //pb
  globalVarsDumper_ = new GlobalVariablesDumper( iConfig.getParameter<edm::ParameterSet>( "globalVariables" ) );
  
  eventTree = fs.make<TTree>( "event", "event" );
  
}*/

SimpleTreeMaker::SimpleTreeMaker( const edm::ParameterSet &iConfig, TFileDirectory& fs, edm::ConsumesCollector && cc ):
  edm::BasicAnalyzer::BasicAnalyzer(iConfig, fs),
  genParticleToken_( cc.consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag>( "genParticleTag" ) ) ),
  genInfoToken_(cc.consumes<GenEventInfoProduct>( iConfig.getParameter<InputTag> ( "generatorInfo" ) ) ),
  PileUpToken_(cc.consumes<View<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ),
  vertexToken_( cc.consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
  diphotonToken_( cc.consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
  mvaResultToken_( cc.consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag> ( "MVAResultTag" ) ) ),
  inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
  genJetToken_( cc.consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "GenJetTag" ) ) ),
  electronToken_( cc.consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
  muonToken_( cc.consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
  METToken_( cc.consumes<View<pat::MET> >( iConfig.getParameter<InputTag> ( "METTag" ) ) ),
  triggerBitsToken_( cc.consumes<edm::TriggerResults>( iConfig.getParameter<InputTag>( "triggerBits" ) ) ),
  rhoToken_(cc.consumes<double>(iConfig.getParameter <edm::InputTag>("rhoFixedGridCollection" ) ) )
{
  jetPtThreshold_ = iConfig.getUntrackedParameter<double>( "jetPtThreshold", 20. );
  bTag_ = iConfig.getUntrackedParameter<string> ( "bTag", "pfCombinedInclusiveSecondaryVertexV2BJetTags" );
  electronPtThreshold_ = iConfig.getUntrackedParameter<double>( "electronPtThreshold", 20. );
  muonPtThreshold_ = iConfig.getUntrackedParameter<double>( "muonPtThreshold", 20. );
  isControlSample_ = iConfig.getUntrackedParameter<bool>( "isControlSample", false );
  lumiWeight_ = iConfig.getUntrackedParameter<double>( "lumiWeight", 1000. ); //pb                                                                                                                              
  
  globalVarsDumper_ = new GlobalVariablesDumper( iConfig.getParameter<edm::ParameterSet>( "globalVariables" ), std::forward<edm::ConsumesCollector>(cc) );



  for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
    auto token = cc.consumes<View<flashgg::Jet> >(inputTagJets_[i]);
    tokenJets_.push_back(token);
  }

  
  eventTree = fs.make<TTree>( "event", "event" );
  
}



SimpleTreeMaker::~SimpleTreeMaker()
{
}
// ******************************************************************************************


// ******************************************************************************************
// analyzer
//
void SimpleTreeMaker::analyze(const edm::EventBase& evt)
{
  const edm::Event *fullEvent = dynamic_cast<const edm::Event *>(&evt);
  const edm::Event &iEvent = (*fullEvent);  
  

  // access edm objects
  Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken( triggerBitsToken_, triggerBits );
  
  Handle<double> rhoHandle;
  iEvent.getByToken( rhoToken_, rhoHandle );
  double rho = *( rhoHandle.product() );
  
  Handle<View<reco::Vertex> > vertices;
  iEvent.getByToken( vertexToken_, vertices );
  
  Handle<View<flashgg::DiPhotonCandidate> > diphotons;
  iEvent.getByToken( diphotonToken_, diphotons );
    
  Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
  iEvent.getByToken( mvaResultToken_, mvaResults );
  
  //JetCollectionVector Jets( inputTagJets_.size() );
  //for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
  //  iEvent.getByLabel( inputTagJets_[j], Jets[j] );
  //}

  JetCollectionVector Jets( inputTagJets_.size() );
  for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
    iEvent.getByToken( tokenJets_[j], Jets[j] );
  }

  
  Handle<View<flashgg::Electron> > electrons;
  iEvent.getByToken( electronToken_, electrons );
  
  Handle<View<flashgg::Muon> > muons;
  iEvent.getByToken( muonToken_, muons );
    
    Handle<View<pat::MET> > METs;
    iEvent.getByToken( METToken_, METs );
    if( METs->size() != 1 )
        { std::cout << "WARNING number of MET is not equal to 1" << std::endl; }
    Ptr<pat::MET> theMET = METs->ptrAt( 0 );

    // only if MC
    Handle<GenEventInfoProduct> genInfo;
    Handle<View<reco::GenJet> > genJets;
    Handle<View< PileupSummaryInfo> > PileupInfos;
    Handle<View<reco::GenParticle> > genParticles;
    if ( !iEvent.isRealData() ) {
        iEvent.getByToken( genInfoToken_, genInfo );
        iEvent.getByToken( genJetToken_, genJets );
        iEvent.getByToken( PileUpToken_, PileupInfos );
        iEvent.getByToken( genParticleToken_, genParticles );
    }
    
    //for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
    //    Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
    //    cout << " pdgId = "<< gen->pdgId()<< " prompt final state = "<< gen->isPromptFinalState() << "  status = " << gen->status() << "   isPrompt = " << gen->statusFlags().isPrompt() <<endl;
    //}    
    
    // -- initialize tree
    initEventStructure();
    
    // -- check if event passes HLT: "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*"  
    const edm::TriggerNames &triggerNames = iEvent.triggerNames( *triggerBits );
    //vector<std::string> const &names = triggerNames.triggerNames();  
    for( unsigned index = 0; index < triggerNames.size(); ++index ) {
      //if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Diphoton") ) 
      //cout << (triggerNames.triggerName( index )).c_str() <<endl;
      if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v") ) {
	//cout << TString::Format((triggerNames.triggerName( index )).c_str()) << " " << triggerBits->accept( index ) << endl;
            evInfo.passHLT =  triggerBits->accept( index );
        }
    }
       
    // -- pre-select best di-photon pair
    //    * pt cut, diphoton preselection, id mva cut on leading and subleading photons 
    //    * if more then one di-photon candidate, take the one with highest sumpt = pt_lead+pt_sublead (DiPhotonCandidates are ordered by decreasing sumpt) 
    //    * di-pho mva cut not applied, needs optimization
    int bestIndex = -1;
    for ( unsigned int idipho = 0; idipho < diphotons->size(); idipho++){
        edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( idipho );        
        //if (! iEvent.isRealData() && (dipho->leadingPhoton()->genMatchType()!=1  || dipho->subLeadingPhoton()->genMatchType()!=1 )) continue;
        //ngen++;
        //if (! passDiphotonPreselection(dipho, rho)) continue; 
        //npre++;
        if (! iEvent.isRealData() && (dipho->leadingPhoton()->genMatchType()==1 && dipho->subLeadingPhoton()->genMatchType()==1 )){
            //cout << dipho->leadingPhoton()->pt() << "  " << dipho->subLeadingPhoton()-> pt() <<endl;
            ngen++;
            npre++;
        }
                
        // if not control sample, apply loose photon id mva cut 
        //- photon id mva cut (~99% efficient on signal photons after preselection)
        if (!isControlSample_) {
            if (dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() ) < -0.9 ) continue;
            if (dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() ) < -0.9 ) continue;
        }

        // - pt threshold
        if (dipho->leadingPhoton()->pt() < dipho->mass()/3. ) continue;
        if (dipho->subLeadingPhoton()->pt() < dipho->mass()/4. ) continue;

        bestIndex = idipho;
        if (! iEvent.isRealData() && (dipho->leadingPhoton()->genMatchType()==1 && dipho->subLeadingPhoton()->genMatchType()==1 )){
            nfullpre++;
        }
        break;
    }

    
    
    // -- analyze event if there is at least one good di-photon candidate
    
    if ( bestIndex > -1) {
        
        // -- event weight (Lumi x cross section x gen weight)
        float w = 1.;
        if( ! iEvent.isRealData() ) {
            w = lumiWeight_;
            if( genInfo.isValid() ) {
                const auto &weights = genInfo->weights();
                if( ! weights.empty() ) {
                    w *= weights[0];
                }
            }
        }
        evInfo.weight = w;

        // -- pileup weights
        globalVarsDumper_->fill( iEvent );
        evInfo.run = globalVarsDumper_->cache().run;
        evInfo.lumi = globalVarsDumper_->cache().lumi;
        evInfo.event = globalVarsDumper_->cache().event;

        if( globalVarsDumper_->puReWeight() ) {
            evInfo.puweight = globalVarsDumper_->cache().puweight;
        }


        // -- number of pileup events
        float pu = 0.; 
        if( ! iEvent.isRealData() ) {
            for( unsigned int PVI = 0; PVI < PileupInfos->size(); ++PVI ) {
                Int_t pu_bunchcrossing = PileupInfos->ptrAt( PVI )->getBunchCrossing();
                if( pu_bunchcrossing == 0 ) {
                    pu = PileupInfos->ptrAt( PVI )->getPU_NumInteractions();
                }
            }
        }
        evInfo.npu = pu;
        
        // -- number of reco vertices
        evInfo.nvtx = vertices->size() ;

        // -- photons
        edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( bestIndex );
        edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( bestIndex );
        
        evInfo.pho1_e  = dipho->leadingPhoton()->energy();
        evInfo.pho1_pt  = dipho->leadingPhoton()->pt();
        evInfo.pho1_eta = dipho->leadingPhoton()->eta();
        evInfo.pho1_phi = dipho->leadingPhoton()->phi();
        evInfo.pho1_idmva = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        evInfo.pho1_genMatchType = dipho->leadingPhoton()->genMatchType();
        
        evInfo.pho2_e  = dipho->subLeadingPhoton()->energy();
        evInfo.pho2_pt  = dipho->subLeadingPhoton()->pt();
        evInfo.pho2_eta = dipho->subLeadingPhoton()->eta();
        evInfo.pho2_phi = dipho->subLeadingPhoton()->phi();
        evInfo.pho2_idmva = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        evInfo.pho2_genMatchType = dipho->subLeadingPhoton()->genMatchType();
                
        evInfo.dipho_pt  = dipho->pt();
        evInfo.dipho_m   = dipho->mass();
        evInfo.dipho_mva = mvares->result ;
        
        // -- jets
        // take the jets corresponding to the diphoton candidate
        unsigned int jetCollectionIndex = diphotons->ptrAt( bestIndex )->jetCollectionIndex();
        
        int njets = 0;
        
        for( UInt_t ijet = 0; ijet < Jets[jetCollectionIndex]->size() ; ijet++ ) {
            
            Ptr<flashgg::Jet> jet  = Jets[jetCollectionIndex]->ptrAt( ijet );
            
            float dRJetPhoLead = deltaR(jet->eta(), jet->phi(), dipho->leadingPhoton()->eta(), dipho->leadingPhoton()->phi());
            float dRJetPhoSubLead = deltaR(jet->eta(), jet->phi(), dipho->subLeadingPhoton()->eta(), dipho->subLeadingPhoton()->phi());

            if( dRJetPhoLead < 0.4 || dRJetPhoSubLead < 0.4 ) { continue; }
            if( !jet->passesJetID( flashgg::Loose ) ) continue;// pass jet id (reject surios detector noise)
            if( !jet->passesPuJetId(diphotons->ptrAt( bestIndex ))){ continue;} // pass PU jet id (always = 1, not implemented yet)
            if( jet->pt() < jetPtThreshold_ ) { continue; }
            
            njets++;
            
            // matching to gen jets
            int isMatchedToGen = 0; 
            if( ! iEvent.isRealData() ) {
                for( unsigned int jg = 0 ; jg < genJets->size() ; jg++ ) {
                    float dr = deltaR(jet->eta(), jet->phi(), genJets->ptrAt( jg )->eta() , genJets->ptrAt( jg )->phi() );
                    if (dr > 0.4) continue;
                    isMatchedToGen = 1;
                }
            }
            
            evInfo.jet_e.push_back(jet->energy());
            evInfo.jet_pt.push_back(jet->pt());
            evInfo.jet_eta.push_back(jet->eta());
            evInfo.jet_phi.push_back(jet->phi());
            evInfo.jet_bdiscriminant.push_back(jet->bDiscriminator( bTag_ ));
            evInfo.jet_partonFlavour.push_back(jet->partonFlavour());
            evInfo.jet_hadronFlavour.push_back(jet->hadronFlavour());
            evInfo.jet_isMatchedToGen.push_back(isMatchedToGen);
        }
        
        // -- leptons (e, mu)

        // -- electrons
        for (UInt_t iele = 0 ; iele < electrons->size(); iele++){
            edm::Ptr<flashgg::Electron> electron = electrons->ptrAt( iele );
            if (fabs(electron->eta()) > 2.4) { continue; }
            if (electron->pt()  < electronPtThreshold_) { continue; }
            if( electron->hasMatchedConversion() ) { continue; } // remove conversions
            // missing hits: from cut-based selection: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
            //if( electron->isEB() && electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS ) > 2 ) { continue; } 
            //if( electron->isEE() && electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS ) > 1 ) { continue; }
	    
            /*
            // dr ele, photon
            TLorentzVector elec_p4;
            elec_p4.SetXYZT( electron->px(), electron->py(), electron->pz(), electron->energy() );
            
            float phi = electron->superCluster()->phi();
            float theta = ( 2 * atan( exp( -electron->superCluster()->eta() ) ) );
            float energy = electron->ecalEnergy();
            float px = energy * sin( theta ) * cos( phi );
            float py = energy * sin( theta ) * sin( phi );
            float pz = energy * cos( theta );
            
            TLorentzVector elec_superClusterVect;
            elec_superClusterVect.SetXYZT( px, py, pz, energy );
            TLorentzVector p1, p2;
            p1.SetXYZT( dipho->leadingPhoton()->px(), dipho->leadingPhoton()->py(), dipho->leadingPhoton()->pz(), dipho->leadingPhoton()->energy() );
            p2.SetXYZT( dipho->subLeadingPhoton()->px(), dipho->subLeadingPhoton()->py(), dipho->subLeadingPhoton()->pz(), dipho->subLeadingPhoton()->energy() );
            cout << "drSCPho1 = " << p1.DeltaR( elec_superClusterVect) <<endl;;
            cout << "drSCPho2 = " << p2.DeltaR( elec_superClusterVect) <<endl;;
            */

            //save dr(gamma, gsf) for checking electron/photon overlap
            float TrkElecSCDeltaR1 = 999.;
            float TrkElecSCDeltaR2 = 999.;
            if( &( *(dipho->leadingPhoton())->superCluster() ) == &( *electron->superCluster() ) ) {
                TrkElecSCDeltaR1 = sqrt( electron->deltaEtaSuperClusterTrackAtVtx() * electron->deltaEtaSuperClusterTrackAtVtx() +
                                              electron->deltaPhiSuperClusterTrackAtVtx() * electron->deltaPhiSuperClusterTrackAtVtx() );
            }
            if( &( *(dipho->subLeadingPhoton())->superCluster() ) == &( *electron->superCluster() ) ) {
                TrkElecSCDeltaR2 = sqrt( electron->deltaEtaSuperClusterTrackAtVtx() * electron->deltaEtaSuperClusterTrackAtVtx() +
                                              electron->deltaPhiSuperClusterTrackAtVtx() * electron->deltaPhiSuperClusterTrackAtVtx() );
            }
            evInfo.ele_drGsfToPho1.push_back(TrkElecSCDeltaR1);
            evInfo.ele_drGsfToPho2.push_back(TrkElecSCDeltaR2);

            //cout << "drGsfPho1 = " << TrkElecSCDeltaR1 <<endl;;
            //cout << "drGsfPho2 = " << TrkElecSCDeltaR2 <<endl;;
            
            Ptr<reco::Vertex> ele_vtx = chooseElectronVertex( electron,  vertices->ptrs() );
            float d0 = electron->gsfTrack()->dxy( ele_vtx->position() );
            float dz = electron->gsfTrack()->dz( ele_vtx->position() );
            float isol = electronIsolation(electron, rho); 
	    int nhits = electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS ); 
	    
	    int passCutBasedIdLoose = passCutBasedElectronIdLoose( electron, vertices->ptrs() );
	    
            int mcMatch = -1;
            if( ! iEvent.isRealData() )
                mcMatch = electronMatchingToGen(electron, genParticles); 

            evInfo.ele_e.push_back(electron->energy());
            evInfo.ele_pt.push_back(electron->pt());
            evInfo.ele_eta.push_back(electron->eta());
            evInfo.ele_phi.push_back(electron->phi());
            evInfo.ele_idmva.push_back(electron->nonTrigMVA());
            evInfo.ele_iso.push_back(isol);
            evInfo.ele_d0.push_back(d0);
            evInfo.ele_dz.push_back(dz);
            evInfo.ele_nMissingHits.push_back(nhits);
	    evInfo.ele_passCutBasedIdLoose.push_back(passCutBasedIdLoose);
            evInfo.ele_isMatchedToGen.push_back(mcMatch);
        }       


        for (UInt_t imu = 0 ; imu < muons->size(); imu++){
            edm::Ptr<flashgg::Muon> muon = muons->ptrAt( imu );
            if (fabs(muon->eta()) > 2.4) { continue; }
            if (muon->pt()  < muonPtThreshold_) { continue; }
            // muon ID and isolation: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
            float muPFCombRelIso = ( muon->pfIsolationR04().sumChargedHadronPt + max( 0.,muon->pfIsolationR04().sumNeutralHadronEt + muon->pfIsolationR04().sumPhotonEt - 0.5 * muon->pfIsolationR04().sumPUPt ) ) / ( muon->pt() );

            int vtxInd = 0;
            double dzmin = 9999;
            for( size_t ivtx = 0 ; ivtx < vertices->size(); ivtx++ ) {
                Ptr<reco::Vertex> vtx = vertices->ptrAt(ivtx);
                if( !muon->innerTrack() ) { continue; }
                if( fabs( muon->innerTrack()->vz() - vtx->position().z() ) < dzmin ) {
                    dzmin = fabs( muon->innerTrack()->vz() - vtx->position().z() );
                    vtxInd = ivtx;
                }
            }
            Ptr<reco::Vertex> muonVtx = vertices->ptrAt(vtxInd);

            int mcMatch =  -1;
            if( ! iEvent.isRealData() ) mcMatch = muonMatchingToGen(muon, genParticles); 

            evInfo.mu_pt.push_back(muon->pt());
            evInfo.mu_eta.push_back(muon->eta());
            evInfo.mu_phi.push_back(muon->phi());
            evInfo.mu_iso.push_back(muPFCombRelIso);
            evInfo.mu_isTight.push_back(muon::isTightMuon( *muon, *muonVtx ));
            evInfo.mu_isMedium.push_back(muon::isMediumMuon( *muon ));
            evInfo.mu_isLoose.push_back(muon::isLooseMuon( *muon ));
            evInfo.mu_isMatchedToGen.push_back(mcMatch); 
        }


        // -- MET 
        evInfo.met = theMET->corPt();
        evInfo.metx = theMET->corPx();
        evInfo.mety = theMET->corPy();
        evInfo.metphi = theMET->corPhi();
        evInfo.metSumEt = theMET->corSumEt();

        
        // --- fill the tree
        //if ( njets > 0. ) // fill only if min number of jets?
        eventTree->Fill();
    }
}
// ******************************************************************************************


// ******************************************************************************************
void
SimpleTreeMaker::beginJob()
{
  ngen = 0;
  npre = 0;
  nfullpre =0 ;

  // per-event tree

  eventTree->Branch( "run", &evInfo.run, "run/I" );
  eventTree->Branch( "lumi", &evInfo.lumi, "lumi/I" );
  eventTree->Branch( "event", &evInfo.event, "event/I" );

  eventTree->Branch( "weight", &evInfo.weight, "weight/F" );
  eventTree->Branch( "puweight", &evInfo.puweight,"puweight/F");

  eventTree->Branch( "npu", &evInfo.npu, "npu/I" );
  eventTree->Branch( "nvtx", &evInfo.nvtx, "nvtx/I" );

  eventTree->Branch( "passHLT", &evInfo.passHLT, "passHLT/I" );
  
  eventTree->Branch( "pho1_e", &evInfo.pho1_e, "pho1_e/F" );
  eventTree->Branch( "pho1_pt", &evInfo.pho1_pt, "pho1_pt/F" );
  eventTree->Branch( "pho1_eta", &evInfo.pho1_eta, "pho1_eta/F" );
  eventTree->Branch( "pho1_phi", &evInfo.pho1_phi, "pho1_phi/F" );
  eventTree->Branch( "pho1_idmva", &evInfo.pho1_idmva, "pho1_idmva/F" );
  eventTree->Branch( "pho1_genMatchType", &evInfo.pho1_genMatchType, "pho1_genMatchType/I" );

  eventTree->Branch( "pho2_e", &evInfo.pho2_e, "pho2_e/F" );
  eventTree->Branch( "pho2_pt", &evInfo.pho2_pt, "pho2_pt/F" );
  eventTree->Branch( "pho2_eta", &evInfo.pho2_eta, "pho2_eta/F" );
  eventTree->Branch( "pho2_phi", &evInfo.pho2_phi, "pho2_phi/F" );
  eventTree->Branch( "pho2_idmva", &evInfo.pho2_idmva, "pho2_idmva/F" );
  eventTree->Branch( "pho2_genMatchType", &evInfo.pho2_genMatchType, "pho2_genMatchType/I" );

  eventTree->Branch( "dipho_pt", &evInfo.dipho_pt, "dipho_pt/F" );
  eventTree->Branch( "dipho_m", &evInfo.dipho_m, "dipho_m/F" );
  eventTree->Branch( "dipho_mva", &evInfo.dipho_mva, "dipho_mva/F" );

  eventTree->Branch( "jet_e", &evInfo.jet_e);
  eventTree->Branch( "jet_pt", &evInfo.jet_pt);
  eventTree->Branch( "jet_eta", &evInfo.jet_eta);
  eventTree->Branch( "jet_phi", &evInfo.jet_phi);
  eventTree->Branch( "jet_pujetid", &evInfo.jet_pujetid);
  eventTree->Branch( "jet_bdiscriminant", &evInfo.jet_bdiscriminant);
  eventTree->Branch( "jet_partonFlavour", &evInfo.jet_partonFlavour);
  eventTree->Branch( "jet_hadronFlavour", &evInfo.jet_hadronFlavour);
  eventTree->Branch( "jet_isMatchedToGen", &evInfo.jet_isMatchedToGen);

  eventTree->Branch( "ele_e", &evInfo.ele_e);
  eventTree->Branch( "ele_pt", &evInfo.ele_pt);
  eventTree->Branch( "ele_eta", &evInfo.ele_eta);
  eventTree->Branch( "ele_phi", &evInfo.ele_phi);
  eventTree->Branch( "ele_idmva", &evInfo.ele_idmva);
  eventTree->Branch( "ele_drGsfToPho1", &evInfo.ele_drGsfToPho1);
  eventTree->Branch( "ele_drGsfToPho2", &evInfo.ele_drGsfToPho2);
  eventTree->Branch( "ele_iso", &evInfo.ele_iso);
  eventTree->Branch( "ele_dz", &evInfo.ele_dz);
  eventTree->Branch( "ele_d0", &evInfo.ele_d0);
  eventTree->Branch( "ele_nMissingHits", &evInfo.ele_nMissingHits);
  eventTree->Branch( "ele_passCutBasedIdLoose", &evInfo.ele_passCutBasedIdLoose);
  eventTree->Branch( "ele_isMatchedToGen", &evInfo.ele_isMatchedToGen);

  eventTree->Branch( "mu_pt", &evInfo.mu_pt);
  eventTree->Branch( "mu_eta", &evInfo.mu_eta);
  eventTree->Branch( "mu_phi", &evInfo.mu_phi);
  eventTree->Branch( "mu_iso", &evInfo.mu_iso);
  eventTree->Branch( "mu_isTight", &evInfo.mu_isTight);
  eventTree->Branch( "mu_isMedium", &evInfo.mu_isMedium);
  eventTree->Branch( "mu_isLoose", &evInfo.mu_isLoose);
  eventTree->Branch( "mu_isMatchedToGen", &evInfo.mu_isMatchedToGen);

  eventTree->Branch( "met", &evInfo.met);
  eventTree->Branch( "metx", &evInfo.metx);
  eventTree->Branch( "mety", &evInfo.mety);
  eventTree->Branch( "metphi", &evInfo.metphi);  
  eventTree->Branch( "metSumEt", &evInfo.metSumEt);

}
// ******************************************************************************************


// ******************************************************************************************
void
SimpleTreeMaker::endJob()
{
    cout << "Total nuber of events before preselection = "<< ngen << endl;
    cout << "Number of events after preselection       = "<< npre << endl;
    cout << "Number of events after full preselection  = "<< nfullpre << endl;
 
} // end of endJob
// ******************************************************************************************


// ******************************************************************************************
void
SimpleTreeMaker::initEventStructure()
{
    // per-event tree:
    evInfo.weight = -999.;
    evInfo.puweight = -999.;

    evInfo.run = -999;
    evInfo.lumi = -999.;
    evInfo.event = -999.;

    evInfo.npu = -999;
    evInfo.nvtx = -999;
    evInfo.passHLT = -1;
    
    evInfo.pho1_e  = -999.;
    evInfo.pho1_pt  = -999.;
    evInfo.pho1_eta = -999.;
    evInfo.pho1_phi = -999.;
    evInfo.pho1_idmva = -999.;
    evInfo.pho1_genMatchType = -999.;
    
    evInfo.pho2_e  = -999.;
    evInfo.pho2_pt  = -999.;
    evInfo.pho2_eta = -999.;
    evInfo.pho2_phi = -999.;
    evInfo.pho2_idmva = -999.;
    evInfo.pho2_genMatchType = -999.;

    evInfo.dipho_pt   = -999.;
    evInfo.dipho_m    = -999.;
    evInfo.dipho_mva  = -999.;
    
    evInfo.jet_e .clear();
    evInfo.jet_pt .clear();
    evInfo.jet_eta .clear();
    evInfo.jet_phi .clear();
    evInfo.jet_pujetid .clear();
    evInfo.jet_bdiscriminant .clear();
    evInfo.jet_partonFlavour .clear();
    evInfo.jet_hadronFlavour .clear();
    evInfo.jet_isMatchedToGen .clear();
    
    evInfo.ele_e .clear();
    evInfo.ele_pt .clear();
    evInfo.ele_eta .clear();
    evInfo.ele_phi .clear();
    evInfo.ele_idmva .clear();
    evInfo.ele_drGsfToPho1 .clear();
    evInfo.ele_drGsfToPho2 .clear();
    evInfo.ele_iso .clear();
    evInfo.ele_dz .clear();
    evInfo.ele_d0 .clear();
    evInfo.ele_nMissingHits .clear();
    evInfo.ele_passCutBasedIdLoose .clear();
    evInfo.ele_isMatchedToGen .clear();

    evInfo.mu_pt .clear();
    evInfo.mu_eta .clear();
    evInfo.mu_phi .clear();
    evInfo.mu_iso .clear();
    evInfo.mu_isTight .clear();
    evInfo.mu_isMedium .clear();
    evInfo.mu_isLoose .clear();
    evInfo.mu_isMatchedToGen .clear();

    evInfo.met = -999;
    evInfo.metx = -999;
    evInfo.mety = -999;
    evInfo.metphi = -999;
    evInfo.metSumEt = -999;

}
// ******************************************************************************************


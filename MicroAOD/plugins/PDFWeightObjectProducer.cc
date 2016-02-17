#include <memory>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <cctype>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "flashgg/DataFormats/interface/PDFWeightObject.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "DataFormats/PatCandidates/interface/libminifloat.h"
#include "PhysicsTools/HepMCCandAlgos/interface/PDFWeightsHelper.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace edm;

namespace flashgg {
    
	class PDFWeightProducer : public edm::EDProducer
	{
    public:
        PDFWeightProducer( const edm::ParameterSet & );
        void set_generator_type( vector<string> );
    private:
        void produce( edm::Event &, const edm::EventSetup & );
        EDGetTokenT<LHEEventProduct> LHEEventToken_;
        EDGetTokenT<GenEventInfoProduct> srcTokenGen_;
        string tag_;
        string delimiter_1_;
        string delimiter_2_;
        string delimiter_3_;
        void beginRun( edm::Run const &, edm::EventSetup const &iSetup );
        vector<int> weight_indices;
        vector<int> alpha_indices;
        vector<int> scale_indices;
        string removeSpace( string line );
        PDFWeightsHelper pdfweightshelper_;//tool from HepMCCandAlgos/interface/PDFWeightsHelper
        unsigned int nPdfEigWeights_;
        //		        std::vector<float> pdfeigweights_;
        edm::FileInPath mc2hessianCSV;
        std::vector<double> inpdfweights;
        float gen_weight;
        string pdfset_;
        string pdfid_1;
        string pdfid_2;
        string alphas_id_1;
        string alphas_id_2;
        bool flashgg_flag_;
        bool pdfset_flag_;
        bool alpha_s_flag_;
        bool scale_flag_;
        string generatorType ;
	};
    
	PDFWeightProducer::PDFWeightProducer( const edm::ParameterSet &iConfig ):
		LHEEventToken_( consumes<LHEEventProduct>( iConfig.getParameter<InputTag>( "LHEEventTag" ) ) ),
        srcTokenGen_( consumes<GenEventInfoProduct>( iConfig.getParameter<InputTag>("GenTag") ) )
	{
        consumes<LHERunInfoProduct,edm::InRun> (edm::InputTag("externalLHEProducer"));
        pdfset_flag_ = iConfig.getUntrackedParameter<bool>("pdfset_flag",false);
        flashgg_flag_ = iConfig.getUntrackedParameter<bool>("flashgg_flag",true);       
        alpha_s_flag_ = iConfig.getUntrackedParameter<bool>("alpha_s_flag","true");
        scale_flag_ = iConfig.getUntrackedParameter<bool>("scale_flag","true"); 
		tag_ = iConfig.getUntrackedParameter<string>( "tag", "initrwgt" );
//		pdfid_1 = iConfig.getUntrackedParameter<string>("pdfid_1","0");
//		pdfid_2 = iConfig.getUntrackedParameter<string>("pdfid_2","0");
        pdfset_ = iConfig.getUntrackedParameter<string>("pdfset","NNPDF30_lo_as_0130.LHgrid");
		delimiter_1_ = iConfig.getUntrackedParameter<string>( "delimiter_1", "id=\"" );
		delimiter_2_ = iConfig.getUntrackedParameter<string>( "delimiter_2", "\">" );
        delimiter_3_ = iConfig.getUntrackedParameter<string>( "delimiter_3", "</weightgroup>" );
        nPdfEigWeights_ = iConfig.getParameter<unsigned int>("nPdfEigWeights");
        mc2hessianCSV = iConfig.getParameter<edm::FileInPath>("mc2hessianCSV");

		produces<vector<flashgg::PDFWeightObject> >();
	}

    
    
	void PDFWeightProducer::beginRun( edm::Run const &iRun, edm::EventSetup const &iSetup )
	{

        weight_indices.clear();
        scale_indices.clear();
        alpha_indices.clear();


		Handle<LHERunInfoProduct> run;
		typedef vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
        
        iRun.getByLabel( "externalLHEProducer", run );
		LHERunInfoProduct myLHERunInfoProduct = *( run.product() );
        
        //--- get info from LHERunInfoProduct
		vector<string> weight_lines;
		for( headers_const_iterator iter = myLHERunInfoProduct.headers_begin(); iter != myLHERunInfoProduct.headers_end(); iter++ ) {
            
            vector<string> lines = iter->lines();
            //for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
            //    std::cout << lines.at(iLine);
            //}
            
			if( ( iter->tag() ).compare( tag_ ) == 0 ) {
				//cout << iter->tag() << endl;
				weight_lines = iter->lines();
			}

            

            //if(flashgg_flag_){
            //    PDFWeightProducer::set_generator_type(lines);
            //    //std::cout << " After set_generator_type on " << iter->tag() << " we have pdfids: " << pdfid_1 << " " << pdfid_2 << std::endl;
            //}
            for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
                if (lines[iLine].find("powheg") != std::string::npos){
                    generatorType = "powheg";
                    break;
                }
            }
        }


        // From Josh's slides 13-15: https://indico.cern.ch/event/459797/contribution/2/attachments/1181555/1800214/mcaod-Feb15-2016.pdf
        int pdfidx = 0;
        pdfidx = run->heprup().PDFSUP.first;
        if (pdfidx == -1 && generatorType == "powheg"){
            pdfidx = 260000;    
        }
        
        cout <<" This sample was generated with the following PDFs : " << pdfidx <<endl;
        
        // --- get min and max pdf index for 100 replicas
        pdfid_1 = boost::lexical_cast<std::string>(pdfidx + 1);
        pdfid_2 = boost::lexical_cast<std::string>(pdfidx + 100);

        cout << "PDFs min and max id for MC replicas: " << pdfid_1 << "   " << pdfid_2 <<endl;

        // --- get alphas id
        if (alpha_s_flag_){
         
            if (pdfidx == 292200){
                alphas_id_1 = "292301";
                alphas_id_2 = "292302";
            }
            
            if (pdfidx == 292000){
                alphas_id_1 = "292101";
                alphas_id_2 = "292102";
            }
            
            if (pdfidx == 260000){
                alphas_id_1 = "265000";
                alphas_id_2 = "266000";
            }
            
            if (pdfidx == 260400){
                alphas_id_1 = "265400";
                alphas_id_2 = "266400";
            }
            cout << "alpha_S min and max id             : " << alphas_id_1 << "   " << alphas_id_2 << endl;
        }



        //// Required std::stringstream object
        std::stringstream ss;
        //// Populate
        std::copy(weight_lines.begin(), weight_lines.end(),std::ostream_iterator<std::string>(ss,""));
        cout << ss.str()<<endl;
        boost::property_tree::ptree pt;
        read_xml( ss , pt);
        
        // --- get weight id for pdfs, alpha_s, scale variations
        for( unsigned int iLine = 0; iLine < weight_lines.size(); iLine++ ) {
                        
            // this works for standard centrally produced samples

            // -- pdf and alpha_s variations
            if ( weight_lines[iLine].find("pdfset") != std::string::npos || weight_lines[iLine].find("PDF set") != std::string::npos){

                std::stringstream iss;
                iss.str(weight_lines[iLine]);                
                boost::property_tree::ptree pt;
                read_xml( iss , pt);
                
                string strwid  = pt.get<std::string>("weight.<xmlattr>.id");
                string strw    = pt.get<std::string>("weight");

                int id = stoi(strwid);
                vector<string> strs;
                boost::split(strs, strw, boost::is_any_of("="));
                int pdfindex  = stoi(strs.back());
                                
                if (pdfindex >= stoi(pdfid_1) && pdfindex <= stoi(pdfid_2)){
                    PDFWeightProducer::weight_indices.push_back( id );
                }
                
                if (alpha_s_flag_){
                    if (pdfindex == stoi(alphas_id_1) || pdfindex == stoi(alphas_id_2)){
                        PDFWeightProducer::alpha_indices.push_back( id );
                    }
                }
            }
            // add something similar for non sTaNard samples...
            
            
            // -- scale variations
            if(scale_flag_){

                string m_format = "muR";
                if(pdfset_flag_){
                    m_format= "mur";
                } 
                
                if ( weight_lines[iLine].find(m_format) != std::string::npos ) {
                    std::stringstream iss;
                    iss.str(weight_lines[iLine]);
                    boost::property_tree::ptree pt;
                    read_xml( iss , pt);
                    //std::cout << ">>>>>> weight id = "<< pt.get<std::string>("weight.<xmlattr>.id")<<std::endl;
                    int scaleind = stoi(pt.get<std::string>("weight.<xmlattr>.id"));
                    scale_indices.push_back( scaleind );
                }
            } 
    
        }// end loop over weight lines

    }

	void PDFWeightProducer::produce( Event &evt, const EventSetup & )
	{
		Handle<LHEEventProduct> LHEEventHandle;
		evt.getByToken( LHEEventToken_, LHEEventHandle );

        Handle<GenEventInfoProduct> genInfo;
        evt.getByToken( srcTokenGen_, genInfo );
    
        gen_weight = genInfo->weight();



		std::auto_ptr<vector<flashgg::PDFWeightObject> > PDFWeight( new vector<flashgg::PDFWeightObject> );

		inpdfweights.clear(); 

		flashgg::PDFWeightObject pdfWeight;


        for( unsigned int i = 0; i < LHEEventHandle->weights().size(); i++) {

			int id_i = stoi( LHEEventHandle->weights()[i].id );
			
            //--- get pdf weights
            for( unsigned int j = 0; j < PDFWeightProducer::weight_indices.size(); j++ ){
				int id_j = PDFWeightProducer::weight_indices[j];	
                //cout << " id_i = " << id_i << "   id_j = " << id_j <<endl;
                if( id_i == id_j ){
                    float weight = LHEEventHandle->weights()[i].wgt;
                    inpdfweights.push_back( weight );
				}
			}
            
            //--- get alpha_s weights
            if ( alpha_s_flag_ ){
                for( unsigned int k = 0; k < PDFWeightProducer::alpha_indices.size(); k++ ){
                    int id_k = PDFWeightProducer::alpha_indices[k];
                    if(id_i == id_k ){
                        float alpha = LHEEventHandle->weights()[i].wgt;
                        uint16_t alpha_16 = MiniFloatConverter::float32to16( alpha );
                        pdfWeight.alpha_s_container.push_back(alpha_16);	
                    }
                }
            }
            
            //--- get qcd scale weights
            if ( scale_flag_ ){ 
                for( unsigned int k = 0 ; k < PDFWeightProducer::scale_indices.size() ; k++ ) {
                    int id_k = PDFWeightProducer::scale_indices[k];
                    if ( id_i == id_k ) {
                        float scale = LHEEventHandle->weights()[i].wgt;
                        uint16_t scale_16 = MiniFloatConverter::float32to16( scale );
                        pdfWeight.qcd_scale_container.push_back( scale_16 );
                    }
                }
            }
            
		}// end loop over lhe weights


		cout << "Size of pdf weights    : " << inpdfweights.size() << endl;
		cout << "Size of scale weights  : " << PDFWeightProducer::scale_indices.size() << endl;
		cout << "Size of alpha_s weights: " << PDFWeightProducer::alpha_indices.size() << endl;

		
        // get MCtoHessian pdf weights
        pdfweightshelper_.Init(PDFWeightProducer::weight_indices.size(),nPdfEigWeights_,mc2hessianCSV);
        
        std::vector<double> outpdfweights(nPdfEigWeights_);
        
        double nomlheweight = LHEEventHandle->weights()[0].wgt;
		//cout << "nomlheweight " << nomlheweight << endl;
        
		pdfweightshelper_.DoMC2Hessian(nomlheweight,inpdfweights.data(),outpdfweights.data());
        
        for (unsigned int iwgt=0; iwgt<nPdfEigWeights_; ++iwgt) {
            
            double wgtval = outpdfweights[iwgt];
            
            //cout << "wgtval  " << wgtval << endl;
            
            float real_weight = wgtval*gen_weight/nomlheweight;
            
            uint16_t weight_16 = MiniFloatConverter::float32to16( real_weight );
			
            //cout << "real_weight " << real_weight << endl; 
            //cout << "weight_16 " << weight_16 << endl;
            
            //the is the weight to be used for evaluating uncertainties with hessian weights
            pdfWeight.pdf_weight_container.push_back(weight_16);
            
		}    
        
		PDFWeight->push_back( pdfWeight );

		evt.put( PDFWeight );

        //cout << "FINAL pdf_weight_container size " <<pdfWeight.pdf_weight_container.size() << endl;
        //cout << "FINAL alpha_s_container size " <<pdfWeight.alpha_s_container.size() << endl;
        //cout << "FINAL qcd_scale_container size " <<pdfWeight.qcd_scale_container.size() << endl;
        
        
	}

}

typedef flashgg::PDFWeightProducer FlashggPDFWeightProducer;
DEFINE_FWK_MODULE( FlashggPDFWeightProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

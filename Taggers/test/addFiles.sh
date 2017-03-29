MYDIR=~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170327/

#hadd -f $MYDIR/Sig/output_ggh.root $MYDIR/Sig/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8_*.root
#hadd -f $MYDIR/Sig/output_vbf.root $MYDIR/Sig/VBFHToGG_M125_13TeV_amcatnloFXFX_pythia8_*.root
#hadd -f $MYDIR/Sig/output_vh.root $MYDIR/Sig/VHToGG_M125_13TeV_amcatnloFXFX_pythia8_*.root
#hadd -f $MYDIR/Sig/output_tth.root $MYDIR/Sig/ttHToGG_M125_13TeV_amcatnloFXFX_pythia8_*.root

#hadd -f $MYDIR/Bkg/output_diphoton.root $MYDIR/Bkg/DiPhotonJetsBox_*.root
#hadd -f $MYDIR/Bkg/output_gjet.root $MYDIR/Bkg/GJet*.root
#hadd -f $MYDIR/Bkg/output_qcd.root $MYDIR/Bkg/QCD*.root
#hadd -f $MYDIR/Bkg/output_DYJetsToLL.root $MYDIR/Bkg/DYJetsToLL*.root
#hadd -f $MYDIR/Bkg/output_ZGTo2LG.root $MYDIR/Bkg/ZGTo2LG*root
#hadd -f $MYDIR/Bkg/output_WGToLNuG.root $MYDIR/Bkg/WGToLNuG*.root
#hadd -f $MYDIR/Bkg/output_ZZTo2L2Q.root $MYDIR/Bkg/ZZTo2L2Q*.root
#hadd -f $MYDIR/Bkg/output_WZTo2L2Q.root $MYDIR/Bkg/WZTo2L2Q*.root
#hadd -f $MYDIR/Bkg/output_WW.root $MYDIR/Bkg/WW*.root
#hadd -f $MYDIR/Bkg/output_TTGG_0Jets.root $MYDIR/Bkg/TTGG_0Jets*.root
#hadd -f $MYDIR/Bkg/output_TTGJets.root $MYDIR/Bkg/TTGJets*.root
#hadd -f $MYDIR/Bkg/output_TGJets.root $MYDIR/Bkg/TGJets*.root
#hadd -f $MYDIR/Bkg/output_TTJets.root $MYDIR/Bkg/TTJets*.root


hadd -f $MYDIR/Data/output_data_Run2016B.root $MYDIR/Data/DoubleEG_Run2016B*.root
hadd -f $MYDIR/Data/output_data_Run2016C.root $MYDIR/Data/DoubleEG_Run2016C*.root
hadd -f $MYDIR/Data/output_data_Run2016D.root $MYDIR/Data/DoubleEG_Run2016D*.root
hadd -f $MYDIR/Data/output_data_Run2016E.root $MYDIR/Data/DoubleEG_Run2016E*.root
hadd -f $MYDIR/Data/output_data_Run2016F.root $MYDIR/Data/DoubleEG_Run2016F*.root
hadd -f $MYDIR/Data/output_data_Run2016G.root $MYDIR/Data/DoubleEG_Run2016G*.root
hadd -f $MYDIR/Data/output_data_Run2016H.root $MYDIR/Data/DoubleEG_Run2016H*.root
hadd -f $MYDIR/Data/output_data.root $MYDIR/Data/output_data_Run2016*.root

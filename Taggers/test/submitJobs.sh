#to save on eos
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207_noMuFilter/Sig
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207_noMuFilter/Bkg
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207_noMuFilter/Data

fggRunJobs.py --load jobs_mc_sig_moriond17.json -H -D -P -n 10 -d jobs_vhLepDumper_20170207_noMuFilter_Sig -x cmsRun vhLeptonicTagsDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207_noMuFilter/Sig/

fggRunJobs.py --load jobs_mc_bkg_moriond17.json -H -D -P -n 10 -d jobs_vhLepDumper_20170207_noMuFilter_Bkg -x cmsRun vhLeptonicTagsDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207_noMuFilter/Bkg/

fggRunJobs.py --load jobs_data_moriond17.json -H -D -P -n 10 -d jobs_vhLepDumper_20170207_noMuFilter_Data -x cmsRun vhLeptonicTagsDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207_noMuFilter/Data/



mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207_MuFilter/Sig
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207_MuFilter/Bkg
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207_MuFilter/Data

fggRunJobs.py --load jobs_mc_sig_moriond17.json -H -D -P -n 10 -d jobs_vhLepDumper_20170207_MuFilter_Sig -x cmsRun vhLeptonicTagsDumper_muFilter_cfg.py maxEvents=-1 --no-use-tarball -q 8nh
 --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207_MuFilter/Sig/

fggRunJobs.py --load jobs_mc_bkg_moriond17.json -H -D -P -n 10 -d jobs_vhLepDumper_20170207_MuFilter_Bkg -x cmsRun vhLeptonicTagsDumper_muFilter_cfg.py maxEvents=-1 --no-use-tarball -q 8nh  --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207_MuFilter/Bkg/

fggRunJobs.py --load jobs_data_moriond17.json -H -D -P -n 10 -d jobs_vhLepDumper_20170207_MuFilter_Data -x cmsRun vhLeptonicTagsDumper_muFilter_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207_MuFilter/Data/




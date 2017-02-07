#to save on eos
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207/Sig
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207/Bkg
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207/Data

fggRunJobs.py --load jobs_mc_sig_moriond17.json -H -D -P -n 10 -d jobs_vhLepDumper_20170207_Sig -x cmsRun vhLeptonicTagsDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207/Sig/

fggRunJobs.py --load jobs_mc_bkg_moriond17.json -H -D -P -n 10 -d jobs_vhLepDumper_20170207_Bkg -x cmsRun vhLeptonicTagsDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207/Bkg/

fggRunJobs.py --load jobs_data_moriond17.json -H -D -P -n 10 -d jobs_vhLepDumper_20170207_Data -x cmsRun vhLeptonicTagsDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170207/Data/

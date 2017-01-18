#to save on eos
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170119/Sig
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170119/Bkg
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170119/Data

fggRunJobs.py --load jobs_mc_sig.json -H -D -P -n 10 -d jobs_vhLepDumper_20170119_Sig -x cmsRun vhLeptonicTagsDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170119/Sig/

fggRunJobs.py --load jobs_mc_bkg.json -H -D -P -n 10 -d jobs_vhLepDumper_20170119_Bkg -x cmsRun vhLeptonicTagsDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170119/Bkg/

fggRunJobs.py --load jobs_data.json -H -D -P -n 10 -d jobs_vhLepDumper_20170119_Data -x cmsRun vhLeptonicTagsDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhLepDumper_20170119/Data/

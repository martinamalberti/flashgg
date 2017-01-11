#to save on eos
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhHadDumper_20161223/Sig
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhHadDumper_20161223/Bkg
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhHadDumper_20161223/Data
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhHadDumper_20161223/Data_CS

fggRunJobs.py --load jobs_mc_sig.json -H -D -P -n 10 -d jobs_vhHadDumper_20161223_Sig -x cmsRun vhHadDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhHadDumper_20161223/Sig/

fggRunJobs.py --load jobs_mc_bkg.json -H -D -P -n 10 -d jobs_vhHadDumper_20161223_Bkg -x cmsRun vhHadDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhHadDumper_20161223/Bkg/

fggRunJobs.py --load jobs_data.json -H -D -P -n 10 -d jobs_vhHadDumper_20161223_Data -x cmsRun vhHadDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhHadDumper_20161223/Data/

fggRunJobs.py --load jobs_data.json -H -D -P -n 10 -d jobs_vhHadDumper_20161223_Data_CS -x cmsRun vhHadDumper_CS_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhHadDumper_20161223/Data_CS/

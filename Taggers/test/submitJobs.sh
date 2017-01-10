#to save on eos
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170110/Data
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170110/MC

fggRunJobs.py --load jobsZee_MC.json -H -D -P -n 10 -d jobs_Zee_20161220_MC -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170110/MC

fggRunJobs.py --load jobsZee_SingleElectron.json -H -D -P -n 10 -d jobs_Zee_20161220_Data -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170110/Data
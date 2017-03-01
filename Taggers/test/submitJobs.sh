#to save on eos
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170227/Data
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170227/MC
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170227/MC_LTbinned_v2

fggRunJobs.py --load jobsZee_MC.json -H -D -P -n 10 -d jobs_Zee_20170227_MC -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170227/MC

fggRunJobs.py --load jobsZee_SingleElectron_reMiniAOD.json -H -D -P -n 10 -d jobs_Zee_20170227_Data -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170227/Data

fggRunJobs.py --load jobsZee_MC_LTbinned.json -H -D -P -n 10 -d jobs_Zee_20170227_MC_LTbinned_v2 -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170227/MC_LTbinned_v2

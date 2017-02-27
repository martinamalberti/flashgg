#to save on eos
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170222/Data
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170222/MC
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170222/MC_noSmearings

fggRunJobs.py --load jobsZee_MC.json -H -D -P -n 10 -d jobs_Zee_20170222_MC -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170222/MC

fggRunJobs.py --load jobsZee_SingleElectron_reMiniAOD.json -H -D -P -n 10 -d jobs_Zee_20170222_Data -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170222/Data

fggRunJobs.py --load jobsZee_MC.json -H -D -P -n 10 -d jobs_Zee_20170222_MC_noSmearings -x cmsRun zeeDumper_noSmearings_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170222/MC_noSmearings


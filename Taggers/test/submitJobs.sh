#to save on eos
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170303/Data
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170303/MC
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170303/MC_LTbinned_v2

fggRunJobs.py --load jobsZee_MC.json -H -D -P -n 10 -d jobs_Zee_20170303_MC -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170303/MC

fggRunJobs.py --load jobsZee_SingleElectron_reMiniAOD.json -H -D -P -n 10 -d jobs_Zee_20170303_Data -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170303/Data

fggRunJobs.py --load jobsZee_MC_LTbinned.json -H -D -P -n 10 -d jobs_Zee_20170303_MC_LTbinned_v2 -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170303/MC_LTbinned_v2





mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170303/MC_regrE_v2
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/hgg/hgg_20170303/MC_regrE_v2

fggRunJobs.py --load jobsZee_MC.json -H -D -P -n 10 -d jobs_Zee_20170303_MC_regrE_v2 -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/zee/zee_20170303/MC_regrE_v2

fggRunJobs.py --load jobs_mc_sig_moriond17.json -H -D -P -n 10 -d jobs_Hgg_20170303_MC_regrE_v2 -x cmsRun simpleDiphotonDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/hgg/hgg_20170303/MC_regrE_v2
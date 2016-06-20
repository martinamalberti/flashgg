fggRunJobs.py --load jobs_data.json -H -D -P -n 100 -d jobsData_80X_v0 -x cmsRun tthOptimizationTreeMaker_cfg.py maxEvents=-1 --no-use-tarball -q 8nh
fggRunJobs.py --load jobs_mc_sig.json -H -D -P -n 100 -d jobsSig_80X_v0 -x cmsRun tthOptimizationTreeMaker_cfg.py maxEvents=-1 --no-use-tarball -q 8nh
fggRunJobs.py --load jobs_mc_bkg.json -H -D -P -n 100 -d jobsBkg_80X_v0 -x cmsRun tthOptimizationTreeMaker_cfg.py maxEvents=-1 --no-use-tarball -q 8nh

#to save on eos
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/tth/80X_v1/Data
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/tth/80X_v1/Sig
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/tth/80X_v1/Bkg

fggRunJobs.py --load jobs_data.json -H -D -P -n 50 -d jobsData_80X_v1 -x cmsRun tthOptimizationTreeMaker_cfg.py maxEvents=-1 --no-use-tarball -q cmscaf1nd --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/tth/80X_v1/Data
fggRunJobs.py --load jobs_mc_sig.json -H -D -P -n 10 -d jobsSig_80X_v1 -x cmsRun tthOptimizationTreeMaker_cfg.py maxEvents=-1 --no-use-tarball -q cmscaf1nd --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/tth/80X_v1/Sig
fggRunJobs.py --load jobs_mc_bkg.json -H -D -P -n 100 -d jobsBkg_80X_v1 -x cmsRun tthOptimizationTreeMaker_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/tth/80X_v1/Bkg


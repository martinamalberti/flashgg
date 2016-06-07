fggRunJobs.py --load jobs_data.json -H -D -P -n 100 -d jobsData_80X_v0 -x cmsRun tthOptimizationTreeMaker_cfg.py maxEvents=-1 --no-use-tarball -q 8nh
fggRunJobs.py --load jobs_mc_sig.json -H -D -P -n 100 -d jobsSig_80X_v0 -x cmsRun tthOptimizationTreeMaker_cfg.py maxEvents=-1 --no-use-tarball -q 8nh
fggRunJobs.py --load jobs_mc_bkg.json -H -D -P -n 100 -d jobsBkg_80X_v0 -x cmsRun tthOptimizationTreeMaker_cfg.py maxEvents=-1 --no-use-tarball -q 8nh

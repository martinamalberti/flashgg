#to save on eos
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhHadDumper_20161213

fggRunJobs.py --load jobs_mc_sig.json -H -D -P -n 10 -d jobs_vhHadDumper_20161213 -x cmsRun vhHadDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/vh/vhHadDumper_20161213

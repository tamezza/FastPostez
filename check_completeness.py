#!/usr/bin/env python3
"""
Check local file completeness vs AMI for all zbby DSIDs.
Uses get_dataset_info directly with full dataset names.
"""

import os
import glob
import pyAMI.client
import pyAMI.atlas.api as atlas

LOCAL_BASE = "/data/atlas/tamezza/croland/JetCalibNtuples/samples_bjes_Zbby/bjes_june19"

DATASETS = {
    364418: "mc20_13TeV.364418.MadGraphPythia8EvtGen_ZqqgammaNp0123_DP140_280.deriv.DAOD_PHYS.e5969_s3681_r13144_r13146_p6490",
    364419: "mc20_13TeV.364419.MadGraphPythia8EvtGen_ZqqgammaNp0123_DP280_500.deriv.DAOD_PHYS.e5969_s3681_r13144_r13146_p6490",
    364420: "mc20_13TeV.364420.MadGraphPythia8EvtGen_ZqqgammaNp0123_DP500_1000.deriv.DAOD_PHYS.e5969_s3681_r13144_r13146_p6490",
    364421: "mc20_13TeV.364421.MadGraphPythia8EvtGen_ZqqgammaNp0123_DP1000_2000.deriv.DAOD_PHYS.e5969_s3681_r13144_r13146_p6490",
    364422: "mc20_13TeV.364422.MadGraphPythia8EvtGen_ZqqgammaNp0123_DP2000_inf.deriv.DAOD_PHYS.e5969_s3681_r13144_r13146_p6490",
    364542: "mc20_13TeV.364542.Sherpa_222_NNPDF30NNLO_SinglePhoton_pty_35_70.deriv.DAOD_PHYS.e6788_s3681_r13144_r13146_p6490",
    364543: "mc20_13TeV.364543.Sherpa_222_NNPDF30NNLO_SinglePhoton_pty_70_140.deriv.DAOD_PHYS.e5938_s3681_r13144_r13146_p6490",
    364544: "mc20_13TeV.364544.Sherpa_222_NNPDF30NNLO_SinglePhoton_pty_140_280.deriv.DAOD_PHYS.e5938_s3681_r13144_r13146_p6490",
    364545: "mc20_13TeV.364545.Sherpa_222_NNPDF30NNLO_SinglePhoton_pty_280_500.deriv.DAOD_PHYS.e5938_s3681_r13144_r13146_p6490",
    364546: "mc20_13TeV.364546.Sherpa_222_NNPDF30NNLO_SinglePhoton_pty_500_1000.deriv.DAOD_PHYS.e5938_s3681_r13144_r13146_p6490",
    364547: "mc20_13TeV.364547.Sherpa_222_NNPDF30NNLO_SinglePhoton_pty_1000_E_CMS.deriv.DAOD_PHYS.e6068_s3681_r13144_r13146_p6490",
    500800: "mc20_13TeV.500800.MGPy8EG_tty_yfromdec.deriv.DAOD_PHYS.e8261_s3681_r13144_r13146_p6490",
    504554: "mc20_13TeV.504554.aMCPy8EG_tty_yprod.deriv.DAOD_PHYS.e8261_s3681_r13144_r13146_p6490",
    700352: "mc20_13TeV.700352.Sh_2211_pTZ100_Zqqgamma.deriv.DAOD_PHYS.e8312_s3681_r13144_r13146_p6490",
    700353: "mc20_13TeV.700353.Sh_2211_pTZ100_Zbbgamma.deriv.DAOD_PHYS.e8312_s3681_r13144_r13146_p6490",
    700507: "mc20_13TeV.700507.Sh_2211_pTW140_Wqqgamma.deriv.DAOD_PHYS.e8338_s3681_r13144_r13146_p6490",
    800660: "mc20_13TeV.800660.Py8_gammajet_direct_DP35_50_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800661: "mc20_13TeV.800661.Py8_gammajet_direct_DP50_70_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800662: "mc20_13TeV.800662.Py8_gammajet_direct_DP70_140_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800663: "mc20_13TeV.800663.Py8_gammajet_direct_DP140_280_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800664: "mc20_13TeV.800664.Py8_gammajet_direct_DP280_500_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800665: "mc20_13TeV.800665.Py8_gammajet_direct_DP500_800_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800666: "mc20_13TeV.800666.Py8_gammajet_direct_DP800_1000_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800667: "mc20_13TeV.800667.Py8_gammajet_direct_DP1000_1500_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800668: "mc20_13TeV.800668.Py8_gammajet_direct_DP1500_2000_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800669: "mc20_13TeV.800669.Py8_gammajet_direct_DP2000_2500_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800670: "mc20_13TeV.800670.Py8_gammajet_direct_DP2500_3000_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800671: "mc20_13TeV.800671.Py8_gammajet_direct_DP3000_inf_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800672: "mc20_13TeV.800672.Py8_gammajet_frag_DP8_17_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800673: "mc20_13TeV.800673.Py8_gammajet_frag_DP17_35_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800674: "mc20_13TeV.800674.Py8_gammajet_frag_DP35_50_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800675: "mc20_13TeV.800675.Py8_gammajet_frag_DP50_70_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800676: "mc20_13TeV.800676.Py8_gammajet_frag_DP70_140_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800677: "mc20_13TeV.800677.Py8_gammajet_frag_DP140_280_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800678: "mc20_13TeV.800678.Py8_gammajet_frag_DP280_500_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800679: "mc20_13TeV.800679.Py8_gammajet_frag_DP500_800_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800680: "mc20_13TeV.800680.Py8_gammajet_frag_DP800_1000_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800681: "mc20_13TeV.800681.Py8_gammajet_frag_DP1000_1500_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800682: "mc20_13TeV.800682.Py8_gammajet_frag_DP1500_2000_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
    800683: "mc20_13TeV.800683.Py8_gammajet_frag_DP2000_inf_ordered.deriv.DAOD_PHYS.e8279_s3681_r13144_r13146_p6490",
}

client = pyAMI.client.Client('atlas')

print(f"\n{'DSID':<10} {'Local':>6} {'AMI':>6} {'Status':>12} {'Factor':>8}")
print("-" * 48)

incomplete = []

for dsid, ds_name in DATASETS.items():
    pattern = os.path.join(LOCAL_BASE, f"user.croland.bjes_june19.{dsid}.*", "*.root")
    local_files = [f for f in glob.glob(pattern) if "merged" not in f]
    n_local = len(local_files)

    try:
        info = atlas.get_dataset_info(client, ds_name)
        n_ami = int(info[0].get('nFiles', 0)) if info else 0

        if n_ami == 0:
            status = "AMI=0?"
            factor = "?"
        elif n_local >= n_ami:
            status = "OK"
            factor = "-"
        else:
            status = "INCOMPLETE"
            factor = f"{n_ami/n_local:.1f}x"
            incomplete.append((dsid, n_local, n_ami, ds_name))

    except Exception as e:
        n_ami = -1
        status = "ERR"
        factor = str(e)[:20]

    print(f"{dsid:<10} {n_local:>6} {str(n_ami):>6} {status:>12} {factor:>8}")

print("-" * 48)
print(f"\nSummary: {len(incomplete)} incomplete DSIDs")
for dsid, nl, na, name in incomplete:
    print(f"  DSID {dsid}: {nl}/{na} files  (need {na-nl} more)")
    print(f"    {name}")

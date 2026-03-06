import subprocess
import argparse

4# Run 2 MC20 DSIDs
list_dsid_mc20 = [700398,700399, # Zll+gamma
             345775,345776,345777,345778,345779,345780,345781,345782, # Zll+gamma alternative
             700320, 700321, 700322, 700323, 700324, 700325, # Zll
             700400, 700401, # Z(tau tau, nu nu) + gamma
             700402, 700403, 700404, #W+gamma
             700326, 700327, 700328, 700329, 700330, 700331, 700332, 700333, 700334, #Z (tau, tau) + jets
             366160, 366161, 366162, # VV + gamma
             700352, 700353, 700855, 700849, # Z(qq,bb)+gamma
             364418, 364419, 364420, 364421, 364422, # Zqq+gamma alternative
             364542, 364543, 364544, 364545, 364546, 364547, # gamma+jets
             700507, # Wqq+gamma
             700598, 700599, #Z(qq,bb)+jets
             700855, 700849, #Z(qq,bb)+jets
             811056, 811056, #dijetsjets
             410389, 410087, #tt+gamma
             800658,800659,800660,800661,800662,800663,800664,800665,800666,800667,800668,800669,800670,800671,800672,800673,800674,800675,800676,800677,800678,800679,800680,800681,800682,800683,# alternative gamma+jets,
             700598, 7005989,
             364700, 364701, 364702, 364703, 364704, 364705, 364706, 364707, 364708, 364709, 364710, 364711, 364712, # Dijet
             364686,364687,364688,364689,364690,364691,364692,364693,364694, # Dijet sherpa
             700441, 700843, #Wqq
             410470, 410471, #ttbar
             410644, 410645, 410646, 410647, 410648, 410648, 410658, 410659, # single top
             500800,504554
             ]

# Run 3 MC23 DSIDs (from your screenshot)
list_dsid_mc23 = [700855,  # Z(→bb) + jets (Nominal signal) - Sherpa 2.2.14
                  700849,  # Z(→qq) + jets - Sherpa 2.2.14
                  700843,  # W(→qq) + jets - Sherpa 2.2.14
                  601237,  # Top - PhPy8EG_A14
                  801165, 801166, 801167, 801168, 801169, 801170, 801171, 801172, 801173  # Multijets (Main background) - Pythia8EvtGen_A14
                  ]


def get_metadata(dsid_list, campaign, rtag, ptag, output_file):
    """
    Fetch metadata for a list of DSIDs from AMI
    
    Args:
        dsid_list: List of DSID integers
        campaign: MC campaign (e.g., 'mc20_13TeV', 'mc23_13p6TeV')
        rtag: Reconstruction tag (e.g., 'r13144_r13146', 'r15799')
        ptag: Physics tag (e.g., 'p6490', 'p6697')
        output_file: Output filename
    """
    inDS = '"'
    for dsid in dsid_list:
        inDS += '{}.{}%DAOD_PHYS.%{}%{},'.format(campaign, dsid, rtag, ptag)
    inDS = inDS[:-1] + '"'

    cmd = 'getMetadata.py --fields="dataset_number,crossSection_pb,genFiltEff,kFactor,generator_name" --inDS {} --outFile={}'.format(
        inDS, output_file)

    print(f"Fetching metadata for {campaign}...")
    print(f"Command: {cmd}")
    subprocess.check_call(cmd, shell=True)
    print(f"Metadata saved to {output_file}")

"""
inDS = '"'
for dsid in list_dsid:
    inDS += 'mc20_13TeV.{}%DAOD_PHYS.%r13144_r13146%p6490,'.format(dsid)
inDS = inDS[:-1] + '"'

cmd = 'getMetadata.py --fields="dataset_number,crossSection_pb,genFiltEff,kFactor,generator_name" --inDS {} --outFile=metadata.txt'.format(inDS)

subprocess.check_call(cmd, shell=True)
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fetch MC metadata from AMI')
    parser.add_argument('--campaign', type=str, choices=['mc20', 'mc23', 'both'],
                        default='both', help='MC campaign to fetch')
    parser.add_argument('--output-mc20', type=str, default='metadata_mc20.txt',
                        help='Output file for MC20 metadata')
    parser.add_argument('--output-mc23', type=str, default='metadata_mc23.txt',
                        help='Output file for MC23 metadata')

    args = parser.parse_args()

    if args.campaign in ['mc20', 'both']:
        get_metadata(
            dsid_list=list_dsid_mc20,
            campaign='mc20_13TeV',
            rtag='r13144_r13146',
            ptag='p6490',
            output_file=args.output_mc20
        )

    if args.campaign in ['mc23', 'both']:
        get_metadata(
            dsid_list=list_dsid_mc23,
            campaign='mc23_13p6TeV',
            rtag='r15224_r15225',  # Standard Run 3 rtag
            ptag='p6697',
            output_file=args.output_mc23
        )



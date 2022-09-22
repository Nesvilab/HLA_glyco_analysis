import os
import random

import pandas as pd

# NetMHCIIpan exec file location
NETMHCPAN = ""

# netMHCIIpan-4.1/data/allele.list file
alleles_file = ''
NETMHCPAN_alleles = [x.strip('\n') for x in open(alleles_file).readlines()]


def get_hla_binding(peptide_list, alleles):
    pep_file = "/tmp/pepfile.txt"
    netmhcpan_output = "/tmp/pepfile_netmhcpan"
    stdout = "/tmp/stdout_netmhcpan"
    with open(pep_file, "w") as fh:
        fh.write("\n".join(peptide_list))

    df_list = []
    for allele in alleles.split(','):
        if allele not in NETMHCPAN_alleles:
            print(
                f"Allele: {allele} not regonized by netMHCIIPan 4.1 ==> skipping"
            )
            continue
        os.system(
            f"{NETMHCPAN} -f {pep_file} -a {allele} -inptype 1 -xls -xlsfile {netmhcpan_output} > {stdout}"
        )

        bp_binding = pd.read_csv(netmhcpan_output, sep="\t", header=1)
        bp_binding['Allele'] = allele
        os.system(f"rm {netmhcpan_output} {stdout}")
        df_list.append(bp_binding)
    os.system(f"rm {pep_file}")
    return pd.concat(df_list)

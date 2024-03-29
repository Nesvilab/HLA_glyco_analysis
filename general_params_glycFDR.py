import os

studies = [
    "Bassani_2017_PXD006939", "Bassani_2021_PXD020079",
    "Gerlinger_2019_PXD014017", "Neidert_2021_PXD019643",
    "Neidert_2021_PXD020186", "Purcell_2021_PXD025877",
    "2016_Bassani_PXD004894", "Chong_2020_PXD013649", "Racle_2019_PXD012308"
]
studies = [x + '_glycFDR' for x in studies]
# file_types = ["ion.tsv", "peptide.tsv", "protein.tsv", "psm.tsv"]
file_types = ["peptide.tsv", "protein.tsv", "psm.tsv"]

# root directory
root = ''
out = os.path.join(root, 'reports_glycFDR')

# to retrieve a table path fast
tables = {}
for file_type in file_types:
    key = file_type.rsplit('.', 1)[0]
    tables[key] = os.path.join(out, 'tables', file_type)

# supplementary table  1, sheet 1
tables['annotation'] = ''

# uniprot target decoy database with tag _rev for decoys
database = ''

# MoDec manual inspection table
pwm_annotation_path = ''

# psm table for Purcell_20210_PXD025877
tables[
    'Purcell_20210_PXD025877'] = ''

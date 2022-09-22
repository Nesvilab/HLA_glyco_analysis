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

root = '/storage/data/HLA/_Results'
out = os.path.join(root, 'reports_glycFDR')

# to retrieve a table path fast
tables = {}
for file_type in file_types:
    key = file_type.rsplit('.', 1)[0]
    tables[key] = os.path.join(out, 'tables', file_type)

tables['annotation'] = '/storage/data/HLA/datasets_annotation.tsv'
database = '/storage/dpolasky/projects/HLA_atlas/2019-08-22-td-rev-UP000005640.fas'
pwm_annotation_path = '/storage/data/HLA/_Results/reports_glycFDR/deconvolution/PWM_annotation.tsv'
tables[
    'Purcell_20210_PXD025877'] = '/storage/data/HLA/_Results/Purcell_2021_PXD025877_glycFDR/psm.tsv'

import argparse
import re
import pandas as pd


class Translator():
    def __init__(self, human_cols):
        self._pd_row = pd.DataFrame(data='.', columns=human_cols, index=[0])

    def translate(self, df):
        result_series = df.apply(self.translate_line_to_human, axis=1)
        return pd.concat(result_series.values.tolist(), ignore_index=True)

    def translate_line_to_human(self, line):
        line_w_hum_cols = self._pd_row.merge(line.to_frame().transpose(), how='right').fillna('.')
        semi_split = line_w_hum_cols.at[0, 'INFO'].split(';')
        idx = 1
        for semi in semi_split:
            eq_split = semi.split('=')
            if eq_split[0] in line_w_hum_cols:
                line_w_hum_cols[eq_split[0]] = eq_split[1]
            else:
                comma_split = eq_split[1].split(',')
                for individual_effect in comma_split:
                    if len(line_w_hum_cols.index) - 1 < idx:
                        line_w_hum_cols = pd.concat([line_w_hum_cols, line_w_hum_cols.loc[0].to_frame().transpose()], ignore_index=True)
                    effect_details = re.split(r'\(|\||\)', individual_effect.strip(')'))
                    for ALN_col, ALN_detail in zip(ALNdict[eq_split[0]], effect_details):
                        line_w_hum_cols.loc[idx, ALN_col] = ALN_detail
                    idx = idx + 1
        if len(line_w_hum_cols.index) > 1:
             line_w_hum_cols = line_w_hum_cols.drop([0])
        return line_w_hum_cols


def parse_info_header_to_cols(header, info_id):
    info_line = [line for line in header if 'ID='+info_id in line][0]
    _format = re.search(r'(?<=Format).+', info_line).group()
    details = re.findall(r'\w+', _format)
    col_names = [info_id + '_' + detail for detail in details]
    return col_names


snpEff2HumanVersion = '2.0.0 (2019-12-20)'
print('Hello!')

# Parse args
parser = argparse.ArgumentParser(description='Convert snpEff VCF format to a more human-readable format')
parser.add_argument('inputFile', type=str, help='The snpEff file')
parser.add_argument('-outputFile', type=str, help='The output file')
args = parser.parse_args()
input_file = args.inputFile
output_file = args.outputFile
if output_file is None:
    output_file = input_file.replace('.vcf', '.tsv')

# read header
print('Reading header...')
header = []
with open(input_file, 'r') as f:
    line = f.readline()
    while '#' in line:
        header.append([line])
        line = f.readline()

header_size = len(header)
new_header = header.copy()
new_header[0][0] = new_header[0][0].strip() + '.notReallyButKinda'
new_header.insert(-1, ['##snpEff2HumanVersion=' + snpEff2HumanVersion + ' by Eddie Hujber'])
new_header.insert(-1, ['##snpEff2HumanInputFile=' + input_file])
new_header = [line[0].strip() for line in new_header]

# The header contains lines, each one describes a different subfield within the "INFO" column
info_header = [line for line in new_header if '##INFO=<ID=' in line]
info_cols = [re.search(r'(?<=##INFO=<ID=)\w+', line).group() for line in info_header]
'''
MISSING INFOIDclean LINES
'''

ANNorEFF = 'ANN' if 'ANN' in info_cols else 'EFF'
ALNcols = [ANNorEFF, 'LOF', 'NMD']
ALNdict = {}
for col in ALNcols:
    ALNdict[col] = parse_info_header_to_cols(info_header, col)


print('Translating to human...')
snpEff_DF = pd.read_csv(input_file, sep='\t', header=header_size-1)

human_col_names = []
for col in list(snpEff_DF):
    if 'INFO' not in col:
        human_col_names.append(col)
    else:
        for info in info_cols:
            if info in ALNdict:
                human_col_names += ALNdict[info]
            else:
                human_col_names.append(info)

T = Translator(human_col_names)
human_df = T.translate(snpEff_DF)

print('Writing to file...')
with open(output_file, 'w', newline='') as f:
    f.write('\n'.join(new_header[:-1]) + '\n')
    human_df.to_csv(f, sep='\t', index=False)

print('Done! Enjoy!')

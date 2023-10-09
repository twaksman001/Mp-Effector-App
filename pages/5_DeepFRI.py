##################modules##################

import streamlit as st
import pandas as pd

##################layout##################

st.set_page_config(layout="wide")
st.write('This page creates an interactive spreadsheet from DeepFRI output file for any single effector.')
st.write('If DeepFRI has not been tested for a selected protein, or column filters result in no matching rows, error will occur.')
with st.expander('Column definitions'):
	st.write('salience - list of scores for each consecutive amino acid in the protein primary sequence')

##################load dataframe##################

effectors_structural = ['Mp1', 'Mp2', 'Mp4', 'Mp5', 'Mp6', 'Mp7', 'Mp10', 'Mp11', 'Mp12', 'Mp14', 'Mp15', 'Mp16', 'Mp17', 'Mp19', 'Mp20', 
						'Mp21', 'Mp22', 'Mp23', 'Mp24', 'Mp28', 'Mp29', 'Mp30', 'Mp31', 'Mp32', 'Mp33', 'Mp35', 'Mp36', 'Mp37', 'Mp38', 
						'Mp39', 'Mp40', 'Mp41', 'Mp42', 'Mp43', 'Mp44', 'Mp45', 'Mp46', 'Mp47', 'Mp49', 'Mp50', 'Mp51', 'Mp53', 'Mp54', 
						'Mp55', 'Mp57', 'Mp58', 'Mp60', 'Mp64', 'Mp65', 'Mp66', 'Mp70', 'Mp71', 'Mp72', 'Mp73', 'Mp74', 'Mp76', 'Mp77', 
						'Mp78', 'Mp79', 'Mp81', 'Mp82', 'Mp85', 'Mp90', 'Mp91', 'Mp93', 'Mp94', 'Mp95', 'Mp92a', 'Mp92b', 'MpC002', 'MIF1', 'Mp67']
protein = st.selectbox(label='select effector:', options=effectors_structural, index=0)

def whoami():
	import inspect
	return inspect.stack()[1][3]

def get_file_of_interest(protein, path='Data Files/DeepFRI/*', string='_'):
    
    import glob
    import re
    
    desired_file = None
    
    for file_path in glob.glob(path):
        if re.search(protein+string, file_path):                        
            desired_file = file_path.replace('\\', '/')
    
    if desired_file == None:
        raise ValueError(whoami(), 'fail to get file of interest')
    else:
        return desired_file

def find_DeepFRI_predictions(lines):
    
    import re
    import pandas as pd
    
    GO_terms_ID, GO_terms_name, GO_terms_score, GO_terms_salience = [], [], [], []

    for i in range(len(lines)):
        if 'go_term_score' in lines[i]:
            GO_terms_ID.append(lines[i-2][lines[i-2].find('GO:')+4:lines[i-2].rfind('"')])
            GO_terms_name.append(lines[i-1][lines[i-1].find(': ')+3:lines[i-1].rfind('"')])
            GO_terms_score.append(lines[i][lines[i].find(': ')+2:lines[i].find(': ')+7])
            if 'salience' in lines[i+1]:
                salience = ''
                j = i+2
                while re.search('\d', lines[j]):
                    salience += lines[j]
                    j += 1
                salience = list(map(float, re.split(',\s*', salience[salience.find(',')-1:])))
                GO_terms_salience.append(salience)
            else:
                GO_terms_salience.append('N/A')

    df = pd.DataFrame(data=zip(GO_terms_ID, GO_terms_name, GO_terms_score, GO_terms_salience), columns=['GO_ID', 'GO_name', 'score', 'salience'])

    return df

def DeepFRI_df(protein, path='Data Files/DeepFRI/*', string='_'):

    import pandas as pd
    
    try:
        file = open(get_file_of_interest(protein=protein, path=path, string=string)).read()
    except Exception as e:
        raise ValueError(whoami(), 'fail to get file of interest')
    else:
        dfs = []
        
        lines_top = file[:file.find('contact-map')].split('\n')
        lines_bottom = file[file.find('gcn_bp'):].split('\n')

        for lines in [lines_top, lines_bottom]:
            df = find_DeepFRI_predictions(lines)
            dfs.append(df)

        dfs[0]['DeepFRI_type'], dfs[1]['DeepFRI_type'] = 'sequence', 'structure'
        df = pd.concat(objs=dfs, axis=0).reset_index().drop(columns='index')
        assert df['salience'].str.len().nunique() <= 2

        return df

##################call functions##################	

try:
	df = DeepFRI_df(protein)
except Exception as e:
	st.write('no information available for this protein')

st.dataframe(df)

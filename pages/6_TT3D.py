##################modules##################

import streamlit as st
import pandas as pd
# import h5py

##################layout##################

st.set_page_config(layout="wide")
st.write('This page allows exploration of data from TT3D output files for any single effector.')
# st.write('If TT3D has not been tested for a selected protein, or column filters result in no matching rows, error will occur.')
# with st.expander('Column definitions'):
# 	st.write('salience - list of scores for each consecutive amino acid in the protein primary sequence')

##################load dataframe (top hits)##################

effectors_structural = ['Mp1', 'Mp2', 'Mp4', 'Mp5', 'Mp6', 'Mp7', 'Mp10', 'Mp11', 'Mp12', 'Mp14', 'Mp15', 'Mp16', 'Mp17', 'Mp19', 'Mp20', 
						'Mp21', 'Mp22', 'Mp23', 'Mp24', 'Mp28', 'Mp29', 'Mp30', 'Mp31', 'Mp32', 'Mp33', 'Mp35', 'Mp36', 'Mp37', 'Mp38', 
						'Mp39', 'Mp40', 'Mp41', 'Mp42', 'Mp43', 'Mp44', 'Mp45', 'Mp46', 'Mp47', 'Mp49', 'Mp50', 'Mp51', 'Mp53', 'Mp54', 
						'Mp55', 'Mp57', 'Mp58', 'Mp60', 'Mp64', 'Mp65', 'Mp66', 'Mp70', 'Mp71', 'Mp72', 'Mp73', 'Mp74', 'Mp76', 'Mp77', 
						'Mp78', 'Mp79', 'Mp81', 'Mp82', 'Mp85', 'Mp90', 'Mp91', 'Mp93', 'Mp94', 'Mp95', 'Mp92a', 'Mp92b', 'MpC002', 'MIF1', 'Mp67']

effector_select = st.selectbox(label='choose effector',
								    options=effectors_structural,
								    index=0,
								    key='choose effector'
								    )

threshold = st.slider(label='choose interaction score threshold',
                      min_value=0.0,
                      max_value=1.0,
                      value=0.95,
                      step=0.05,
                      key='interaction score threshold'
                      )

def get_top_hits(protein, threshold=0.95, path='Data Files/TT3D/MpEffectors_AtProteome_All_TT3D_positive.txt'):
    
    import pandas as pd
    
    df = pd.read_csv(filepath_or_buffer=path, sep='\t', names=['effector', 'TAIR AGI', 'score'])
    df = df[(df['effector']==protein) & (df['score']>threshold)]
    
    return df

##################call functions##################	

try:
	df = get_top_hits(protein=effector_select, threshold=threshold)
except Exception as e:
	st.write('error getting dataframe')

st.dataframe(df)

##################load dataframe (contact map)##################

with st.form('select proteins for contact map', clear_on_submit=False):
    target_select = st.selectbox(label='target',
                                    options=df['TAIR AGI'].unique(),
                                    key='target')
    show_contact_map = st.form_submit_button(label='show contact map')

def get_file_LFS(url, output_path):

    import requests

    r = requests.get(url=url)
    open(output_path, mode='w+b').write(r.content)

def get_contact_map(protein1, protein2, path1='Data Files/TT3D/*positive*',
                    path2='Data Files/TT3D/*cmaps*'):
    
    import pandas as pd, os, glob, h5py
        
    found1, found2 = False, False
            
    for file_path in glob.glob(path1):
        file_name, file = file_path.replace('\\', '/')[file_path.replace('\\', '/').rfind('/')+1:], open(file_path).read()

        if all(j in file for j in [protein1, protein2]):
            df = pd.read_csv(filepath_or_buffer=file_path, sep='\t', header=None)

            for i in range(len(df)):
                if set([protein1, protein2]).issubset(set(df.iloc[i, :])):
                    found1 = True
                    break
        else:
            continue
        break
    
    if found1:
        for file_path in glob.glob(path2):
            if file_name[:file_name.find('_TT')+1] in file_path and 'real.h5' not in file_path:
                file_path
                file_name = file_path.replace('\\', '/')[file_path.replace('\\', '/').rfind('/')+1:]
                if not os.path.exists(file_path[:-2]+'real.h5'):
                    get_file_LFS(url='https://github.com/twaksman001/Mp-Effector-App/raw/master/Data%20Files/TT3D/'+file_name,
                                 output_path=file_path[:-2]+'real.h5')
                file = h5py.File(file_path[:-2]+'real.h5', 'r')
                list(file.keys())[:10]
                for key in list(file.keys()):
                    if [protein1, protein2] == key.split('x'):
                        df, found2 = pd.DataFrame(data=file[key]), True
                        break
            else:
                continue
            break
    
    # if found2:
    return df
    # else:
        # return [protein1, protein2, 'score', found1, 'cmaps', found2]

def matrix_heatmap(data, cmap='Greys', vmax=1):#, square=True, cbar=False):
    
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    fig = plt.figure(figsize=(25,25))
    plt.rcParams.update({'font.size': 15})
    sns.heatmap(data=data, cmap=cmap, vmax=vmax, cbar_kws={'shrink':0.2, 'location':'top'})#, square=square, cbar=cbar)
    st.pyplot(fig)

##################call functions##################	

if show_contact_map:
    df1 = get_contact_map(protein1=effector_select, protein2=target_select)
    # if type(df1) == list:
    #     st.write('fail to get contact map dataframe')
    #     st.write(df1)
    # else:
    #     matrix_heatmap(data=df1)
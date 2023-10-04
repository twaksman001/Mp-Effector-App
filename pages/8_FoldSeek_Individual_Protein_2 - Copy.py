# -------------------------------------------------------- PACKAGES --------------------------------------------------------

import pandas as pd  # pip install pandas openpyxl
import plotly.express as px  # pip install plotly-express
import streamlit as st  # pip install streamlit
import math
#import sys

# ------------------------------------------------------------ LAYOUT ------------------------------------------------------------

st.set_page_config(layout='wide')
st.write('This page creates an interactive spreadsheet from FoldSeek output file for any single effector.')
with st.expander('Column definitions'):
	st.write('https://github.com/steineggerlab/foldseek')

# ------------------------------------------------------------ READ EXCEL ------------------------------------------------------------

df = pd.read_csv('Data Files/MpEffectorsBioinfo_Structural_KeyInfo.csv')

# ------------------------------------------------------------ VARIABLES ------------------------------------------------------------

desired_file_path_FoldSeek = 'Data Files/FoldSeek/*.txt'

def get_file_of_interest_FoldSeek(protein):
	
	import glob
	import re
	
	for file_path in glob.glob(desired_file_path_FoldSeek):

		if re.search(protein + '_', file_path) != None:
						
			return file_path.replace('\\', '/')

def FoldSeek_df(protein):
	
	import re
	import pandas as pd
	
	dict_proteinname = {0: 'MpC002', 100: 'MIF1', 92: 'Mp92a'}
	if protein in [0, 100, 92]:
		protein = dict_proteinname[protein]
	else:
		protein = 'Mp' + str(protein)
	
	column_titles = ['query','target','fident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend',
					 'evalue','bits','qlen','tlen','qaln','taln','tca','tseq','taxid','taxon']
	
	if get_file_of_interest_FoldSeek(protein):
	
		file = open(get_file_of_interest_FoldSeek(protein))
		lines = file.read().split('\n')

		columns_dict = {}
		columns_dict['database'] = []

		for i in range(len(column_titles)):

			columns_dict[column_titles[i]] = []
				
		for j in range(len(lines)):

			if 'job' not in lines[j]:

				continue

			line_split = re.split('\t', lines[j][lines[j].find('job.pdb'):])
						
			del line_split[10]

			for i in range(len(column_titles)):

				if i < len(line_split):
				
					columns_dict[column_titles[i]].append(line_split[i])
				
				else:
					
					columns_dict[column_titles[i]].append('N/A')

		database_dict = {'AF-':'AlphaFold', '\d[\da-zA-Z]{3}_':'PDB', 'MGY':'MGnify', 'GMGC':'GMGC'}

		for i in range(len(columns_dict['target'])):

			columns_dict['database'].append('N/A')

			target = columns_dict['target'][i]

			for key in database_dict:

				if re.search(key, target):

					columns_dict['database'][i] = database_dict[key]

		df = pd.DataFrame.from_dict(columns_dict)

		return df
	
	else:
		
		return None
	
def FoldSeek_df_present(protein):
	
	import pandas as pd
	import numpy as np
	
	df = FoldSeek_df(protein)
	
	column_list_strint = ['alnlen','mismatch','gapopen','qstart','qend','tstart','tend','qlen','tlen','bits']
	column_list_strfloat = ['fident','evalue']
	
	if str(df) != 'None':
					
		for column in column_list_strint:
			
			df[column] = df[column].astype(np.int64)
		
		for column in column_list_strfloat:
			
			df[column] = df[column].astype(float)
		
		targetID_list, targetname_list, genus_list, species_list, extra_list = [], [], [], [], []
		
		for target in df['target']:

			target = target + ' '
			
			targetID_list.append(target[:target.find(' ')])
			targetname_list.append(target[target.find(' '):])
		
		for taxon in df['taxon']:
			
			taxon_list = taxon.split(' ')
			
			if len(taxon_list) > 2:
				
				taxon_list[2] = ' '.join(taxon_list[2:])
				del taxon_list[3:]
					   
			while len(taxon_list) < 3:
				
				taxon_list.append('N/A')
			
			genus_list.append(taxon_list[0])
			species_list.append(taxon_list[1])
			extra_list.append(taxon_list[2])
 
		df.insert(3, 'targetID', targetID_list)
		df.insert(4, 'targetname', targetname_list)
		df['genus'], df['species'], df['extra'] = genus_list, species_list, extra_list
		
		df['genus'] = df['genus'] + ' (' + df['genus'].map(df['genus'].value_counts().to_dict()).astype(str) + ')'
		df['species'] = df['species'] + ' (' + df['species'].map(df['species'].value_counts().to_dict()).astype(str) + ')'
		df.insert(24, 'genuscount', df['genus'].map(df['genus'].value_counts().to_dict()))
		df.insert(26, 'speciescount', df['species'].map(df['species'].value_counts().to_dict()))
		
		del (df['query'], df['target'], df['taxon'])
		
		df = df.sort_values(by='bits', axis=0, ascending=False)
		
		return df
	
ID_Number = st.selectbox(label='select Effector:', options=df['ID_No'].unique(), index=0)

if str(FoldSeek_df_present(ID_Number)) != 'None':
	
	df_FoldSeek = FoldSeek_df_present(ID_Number)
	
# ------------------------------------------------------------ SIDEBAR ------------------------------------------------------------
		
def filter_dataframe(df):
	
	# df = df_FoldSeek
	df_filter_dict = {}
	
	columns_stringsearch = {'targetID', 'targetname', 'qaln', 'taln', 'tseq', 'genus', 'species'}
	columns_checkbox = {'evalue'}
	columns_multiselect = {'database'}
	columns_selectslider = {'gapopen', 'qstart', 'qend', 'tstart', 'tend'}
	columns_slider_int = {'alnlen', 'mismatch', 'bits', 'genuscount', 'speciescount'}
	columns_slider_float = {'fident', 'evalue'}
	columns_all_filter = sorted(set.union(columns_stringsearch, columns_checkbox, columns_multiselect, 
									 columns_selectslider, columns_slider_int, columns_slider_float))

	for column in columns_all_filter:
		df_filter_dict[column] = None
	
	with st.sidebar:
		st.header('filters')
		reset_all = st.button('reset all filters')
		st.write('this button shows the unfiltered dataframe, but does not reset filters')
		for column in columns_all_filter:
			with st.expander(column):
				
				if column in columns_stringsearch:
					user_text_input = st.text_input(label=f'text string in {column} (case-sensitive ; if multiple, use only 1 comma to separate)', value='', key=f'text string in {column}')
					if user_text_input:# and reset_all == False:
						df_filter_dict[column] = user_text_input
					# else:
						# df_filter_dict[column] = '.'
				
				if column in columns_checkbox:
					df_filter_dict[f'{column}2'] = None
					threshold_checkbox = st.checkbox(label=f'use threshold for {column}', value=False, key=f'use threshold for {column}')
					if threshold_checkbox:
						user_text_input = st.text_input(label=f'threshold for {column}', value='', key=f'threshold for {column}')
						option = st.selectbox(label='select threshold type:', options=[f'{column} above', f'{column} below'], index=0, key=f'threshold type for {column}')
						set_ = st.checkbox(label=f'activate threshold for {column}', value=False, key=f'activate threshold for {column}')
						df_filter_dict[f'{column}2'] = user_text_input
					# set_ == true
						# if set_ == True and reset_all == False:						
							# if option == 'above':
								# df = df[df[column] >= float(user_text_input)]
							# elif option == 'Below':
								# df = df[df[column] <= float(user_text_input)]					
					# else:
						# user_text_input = None
				
				# if column in columns_multiselect:
					# select_all = st.checkbox(f'select all {column}', value=True)
					# if select_all == True or reset_all == True:				
						# user_cat_input = st.multiselect(
														# f'values for {column}',
														# df[column].unique(),
														# default=list(df[column].unique()),
														# )
					# else:
						# user_cat_input = st.multiselect(
														# f'values for {column}',
														# df[column].unique(),
														# )
					
					# df = df[df[column].isin(user_cat_input)]
				
				if column in columns_multiselect:
					# cat_values = df[column].unique()
					user_cat_input = st.multiselect(
													label=f'values for {column}',
													options=df[column].unique(),
													default=list(df[column].unique()),
													key=f'multiselect values for {column}'
													)
					# st.write(len(user_cat_input), len(cat_values))
					if len(user_cat_input) != df_FoldSeek[column].nunique():
						select_all = st.checkbox(label=f'select all {column}', value=False, key=f'multiselect select all for {column}')
						if select_all == True:# or reset_all == True:				
							st.write(f'all {column} values selected until this box is unchecked')
							user_cat_input = df_FoldSeek[column].unique()
											 # st.multiselect(
															# f'values for {column}',
															# df[column].unique(),
															# )					
					df_filter_dict[column] = user_cat_input
				
				if column in columns_selectslider:	
					# if len(df) != 0:
					if df[column].nunique() > 10:
						_min = df[column].sort_values().min()
						_max = df[column].sort_values().max()
						user_num_cat_input = st.select_slider(
															  label=f'values for {column}',
															  options=df_FoldSeek[column].sort_values(),
															  value=(_min, _max),
															  key=f'selectslider values for {column}'
															  )
						if user_num_cat_input != (_min, _max):
							reset = st.checkbox(label=f'reset {column}', value=False, key=f'reset selectslider values for {column}')
							if reset == True:# or reset_all == True:
								user_num_cat_input = (_min, _max)
								st.write(f'value range for {column} is maximal until this box is unchecked')
						df_filter_dict[column] = user_num_cat_input
					else:
						user_num_cat_input = st.multiselect(
														label=f'values for {column}',
														options=df[column].unique(),
														default=list(df[column].unique()),
														key=f'selectslider multiselect values for {column}'
														)
						# st.write(len(user_cat_input), len(df[column].unique()))
						if len(user_num_cat_input) != df_FoldSeek[column].nunique():
							select_all = st.checkbox(label=f'select all {column}', value=False, key=f'selectslider multiselect select all for {column}')
							if select_all == True:# or reset_all == True:				
								user_num_cat_input = df_FoldSeek[column].unique()
								st.write(f'all {column} values selected until this box is unchecked')
												 # st.multiselect(
																# f'values for {column}',
																# df[column].unique(),
																# )					
					df_filter_dict[column] = user_num_cat_input

				if column in columns_slider_int:
					# st.write(len(df))
					# df.isnull().values.any() == False
					# if df.isnull().values.any() == False:
					# if len(df) != 0:
					_min = int(df[column].min())
					_max = int(df[column].max())
					step = math.ceil((_max - _min)/100)
					if step != 0:
						user_num_input = st.slider(
												   label=f'values for {column}',
												   min_value=_min,
												   max_value=_max,
												   value=(_min, _max),
												   step=step,
												   key=f'slider int values for {column}'
												   )
						# st.write(user_num_input)
						if user_num_input != (_min, _max):
							reset = st.checkbox(label=f'reset {column}', value=False, key=f'reset slider int values for {column}')
							# reset == False
							if reset == True:# or reset_all == True:
								# st.experimental_rerun()
								# break
								# continue
								user_num_input = (_min, _max)
								st.write(f'value range for {column} is maximal until this box is unchecked')
						df_filter_dict[column] = user_num_input
						
				if column in columns_slider_float:
					# if len(df) != 0:
					_min = float(df[column].min())
					_max = float(df[column].max())
					step = (_max - _min)/100
					if step != 0:
						user_num_input = st.slider(
												   label=f'values for {column}',
												   min_value=_min,
												   max_value=_max,
												   value=(_min, _max),
												   step=step,
												   key=f'slider float values for {column}'
												   )
						if user_num_input != (_min, _max):
							reset = st.checkbox(label=f'reset {column}', value=False, key=f'reset slider float values for {column}')
							# reset == False
							if reset == True:# or reset_all == True:
								# st.experimental_rerun()
								# break
								# continue
								user_num_input = (_min, _max)
								st.write(f'value range for {column} is maximal until this box is unchecked')
						df_filter_dict[column] = user_num_input
	
	df_filter_dict
	
	if reset_all:
		# for key in df_filter_dict.keys():
			# df_filter_dict[key] = None
		df = df_FoldSeek
		# for key in st.session_state.keys():
			# del st.session_state[key]
		del st.session_state['slider float values for evalue']
		st.session_state['slider float values for evalue'] = (float(df['evalue'].min()), float(df['evalue'].max()))
	else:
		loc_list = []
		
		for column in columns_all_filter:
			if str(df_filter_dict[column]) == 'None':
				continue
			if column in columns_stringsearch:			
				loc_list.append(set(df[df[column].str.contains(df_filter_dict[column])].index))
			if column in columns_checkbox:
				if threshold_checkbox:
					if set_ == True:# and reset_all == False:						
						if option == f'{column} above':
							loc_list.append(set(df[df[column] >= float(df_filter_dict[f'{column}2'])].index))
						elif option == f'{column} below':
							loc_list.append(set(df[df[column] <= float(df_filter_dict[f'{column}2'])].index))
			if column in columns_multiselect:
				loc_list.append(set(df[df[column].isin(df_filter_dict[column])].index))
			if column in columns_selectslider:
				if 'tuple' in str(type(df_filter_dict[column])):
					loc_list.append(set(df[df[column].between(*df_filter_dict[column])].index))
				if 'list' in str(type(df_filter_dict[column])):
					# st.write(map(int, df_filter_dict[column]))
					# df_filter_dict[column] = list(map(int, df_filter_dict[column]))
					loc_list.append(set(df[df[column].isin(df_filter_dict[column])].index))
			if column in columns_slider_int:
				loc_list.append(set(df[df[column].between(*df_filter_dict[column])].index))
			if column in columns_slider_float:
				loc_list.append(set(df[df[column].between(*df_filter_dict[column])].index))
		
		df = df.loc[list(set.intersection(*loc_list))]
		# st.write(type(loc_list[0]))
		# st.write(loc_list)
		# st.write(loc_list[0].isdisjoint(loc_list[1]))
		# loc_list = []
		
		# if reset_all:
			# for key in df_filter_dict.keys():
				# df_filter_dict[key] = None
			# df = df_FoldSeek
		
		if len(df) == 0:
			st.write('No rows match the given column filters')
			# exit()
	
	st.session_state
	
	return df



df = filter_dataframe(df_FoldSeek)
st.dataframe(df)
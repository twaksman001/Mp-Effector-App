# -------------------------------------------------------- PACKAGES --------------------------------------------------------

import pandas as pd  # pip install pandas openpyxl
import plotly.express as px  # pip install plotly-express
import streamlit as st  # pip install streamlit
import math
#import sys

# ------------------------------------------------------------ LAYOUT ------------------------------------------------------------

st.set_page_config(layout="wide")
st.write('This page creates an interactive spreadsheet from FoldSeek output file for any single effector.')
with st.expander('Column definitions'):
	st.write('https://github.com/steineggerlab/foldseek')

# ------------------------------------------------------------ READ EXCEL ------------------------------------------------------------

df = pd.read_csv('Data Files/MpEffectorsBioinfo_Structural_KeyInfo.csv')

# ------------------------------------------------------------ VARIABLES ------------------------------------------------------------

desired_file_path_FoldSeek = 'Data Files/FoldSeek/*.txt'

def list_unique(lists):
	
	list_unique = []
	
	for list in lists:
		for i in list:
			if i not in list_unique:
				list_unique.append(i)
	
	return list_unique

def df_row_index_list_cond(df, column, list_):
	
	index_list = []
	
	for j in list_:
   
		for i in df.index:
		# for i in range(len(df[column])):
					
			# if any(k == j for k in df[column][i]):
			# if any(k == j for k in df.loc[i, column]):
			if j in df.loc[i, column]:
			
				if i not in index_list:
					
					index_list.append(i)
	
	index_list.sort()
	
	return index_list

def row_string_dynamic_list(df, column):
	
	list_dynamic = []
	
	for string in df[column]:
		if string not in list_dynamic:
			list_dynamic.append(string)
	
	# list_dynamic.sort()

	return list_dynamic

def textinput_row_string_dynamic_list(df, column, textinput):
	
	list_dynamic = []
	
	for string in df[column]:		
		if any(x in string for x in textinput.split(',')):
			list_dynamic.append(string)
	
	return list_dynamic

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
	
ID_Number = st.selectbox("Select Effector:", df["ID_No"].unique())

if str(FoldSeek_df_present(ID_Number)) != 'None':
	
	df_FoldSeek = FoldSeek_df_present(ID_Number)
	
# ------------------------------------------------------------ SIDEBAR ------------------------------------------------------------
		
def filter_dataframe(df):
	
	# df = df_FoldSeek
	df_filter_dict = {}
	
	columns_stringsearch = ['targetID', 'targetname', 'qaln', 'taln', 'tseq', 'genus', 'species']
	columns_checkbox = ['evalue']
	columns_multiselect = ['database']
	columns_selectslider = ['gapopen', 'qstart', 'qend', 'tstart', 'tend']
	columns_slider_int = ['alnlen', 'mismatch', 'bits', 'genuscount', 'speciescount']
	columns_slider_float = ['fident', 'evalue']
	columns_all_filter = list_unique([columns_stringsearch, columns_checkbox, columns_multiselect, 
									 columns_selectslider, columns_slider_int, columns_slider_float])
	
	# df.isnull().values.any() == False
	# st.write(len(df))
	
	for column in columns_all_filter:
		df_filter_dict[column] = None
	
	with st.sidebar:
		st.header('Filters')
		reset_all = st.button('Reset all filters')
		for column in columns_all_filter:
			with st.expander(column):
				
				if column in columns_stringsearch:
					user_text_input = st.text_input(f"Text string in {column} (case-sensitive ; if multiple, use only 1 comma to separate)")
					if user_text_input:# and reset_all == False:
						list_dynamic = textinput_row_string_dynamic_list(df, column, user_text_input)
						df_filter_dict[column] = df_row_index_list_cond(df, column, list_dynamic)
					# else:
						# df_filter_dict[column] = None
				
				if column in columns_checkbox:
					df_filter_dict[f'{column}2'] = None
					threshold_checkbox = st.checkbox(f'Use threshold for {column}')
					if threshold_checkbox:
						user_text_input = st.text_input(f"Threshold for {column}")
						option = st.selectbox("Select threshold type:", [f'{column} above', f'{column} below'])
						set = st.checkbox(f'Activate threshold for {column}')
						df_filter_dict[f'{column}2'] = user_text_input
					# set == True
						# if set == True and reset_all == False:						
							# if option == 'Above':
								# df = df[df[column] >= float(user_text_input)]
							# elif option == 'Below':
								# df = df[df[column] <= float(user_text_input)]					
					# else:
						# user_text_input = None
				
				# if column in columns_multiselect:
					# select_all = st.checkbox(f"Select all {column}", value=True)
					# if select_all == True or reset_all == True:				
						# user_cat_input = st.multiselect(
														# f"Values for {column}",
														# df[column].unique(),
														# default=list(df[column].unique()),
														# )
					# else:
						# user_cat_input = st.multiselect(
														# f"Values for {column}",
														# df[column].unique(),
														# )
					
					# df = df[df[column].isin(user_cat_input)]
				
				if column in columns_multiselect:
					# cat_values = df[column].unique()
					user_cat_input = st.multiselect(
													f"Values for {column}",
													df[column].unique(),
													default=list(df[column].unique()),
													)
					# st.write(len(user_cat_input), len(cat_values))
					if len(user_cat_input) != df_FoldSeek[column].nunique():
						select_all = st.checkbox(f"Select all {column}", value=False)
						if select_all == True:# or reset_all == True:				
							st.write(f'All {column} values selected until this box is unchecked')
							user_cat_input = df_FoldSeek[column].unique()
											 # st.multiselect(
															# f"Values for {column}",
															# df[column].unique(),
															# )					
					df_filter_dict[column] = user_cat_input
				
				if column in columns_selectslider:	
					# if len(df) != 0:
					if df[column].nunique() > 10:
						_min = df[column].sort_values().min()
						_max = df[column].sort_values().max()
						user_num_cat_input = st.select_slider(
															  f"Values for {column}",
															  df_FoldSeek[column].sort_values(),
															  (_min, _max)
															  )
						if user_num_cat_input != (_min, _max):
							reset = st.checkbox(f"Reset {column}", value=False)
							if reset == True:# or reset_all == True:
								user_num_cat_input = (_min, _max)
								st.write(f'Value range for {column} is maximal until this box is unchecked')
						df_filter_dict[column] = user_num_cat_input
					else:
						user_num_cat_input = st.multiselect(
														f"Values for {column}",
														df[column].unique(),
														default=list(df[column].unique()),
														)
						# st.write(len(user_cat_input), len(df[column].unique()))
						if len(user_num_cat_input) != df_FoldSeek[column].nunique():
							select_all = st.checkbox(f"Select all {column}", value=False)
							if select_all == True:# or reset_all == True:				
								user_num_cat_input = df_FoldSeek[column].unique()
								st.write(f'All {column} values selected until this box is unchecked')
												 # st.multiselect(
																# f"Values for {column}",
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
												   f"Values for {column}",
												   min_value=_min,
												   max_value=_max,
												   value=(_min, _max),
												   step=step
												   )
						# st.write(user_num_input)
						if user_num_input != (_min, _max):
							reset = st.checkbox(f"Reset {column}", value=False)
							# reset == False
							if reset == True:# or reset_all == True:
								# st.experimental_rerun()
								# break
								# continue
								user_num_input = (_min, _max)
								st.write(f'Value range for {column} is maximal until this box is unchecked')
						df_filter_dict[column] = user_num_input
						
				if column in columns_slider_float:
					# if len(df) != 0:
					_min = float(df[column].min())
					_max = float(df[column].max())
					step = (_max - _min)/100
					if step != 0:
						user_num_input = st.slider(
												   f"Values for {column}",
												   min_value=_min,
												   max_value=_max,
												   value=(_min, _max),
												   step=step
												   )
						if user_num_input != (_min, _max):
							reset = st.checkbox(f"Reset {column}", value=False)
							# reset == False
							if reset == True:# or reset_all == True:
								# st.experimental_rerun()
								# break
								# continue
								user_num_input = (_min, _max)
								st.write(f'Value range for {column} is maximal until this box is unchecked')
						df_filter_dict[column] = user_num_input
	
	for column in columns_all_filter:
		if str(df_filter_dict[column]) == 'None':
			continue
		if column in columns_stringsearch:			
			df = df.loc[df_filter_dict[column]]
		if column in columns_checkbox:
			if threshold_checkbox:
				if set == True:# and reset_all == False:						
					if option == f'{column} above':
						df = df[df[column] >= float(df_filter_dict[f'{column}2'])]
					elif option == f'{column} below':
						df = df[df[column] <= float(df_filter_dict[f'{column}2'])]
		if column in columns_multiselect:
			df = df[df[column].isin(df_filter_dict[column])]
		if column in columns_selectslider:
			if 'tuple' in str(type(df_filter_dict[column])):
				df = df[df[column].between(*df_filter_dict[column])]
			if 'list' in str(type(df_filter_dict[column])):
				# st.write(map(int, df_filter_dict[column]))
				# df_filter_dict[column] = list(map(int, df_filter_dict[column]))
				df = df[df[column].isin(df_filter_dict[column])]
		if column in columns_slider_int:
			df = df[df[column].between(*df_filter_dict[column])]
		if column in columns_slider_float:
			df = df[df[column].between(*df_filter_dict[column])]

	if reset_all:
		df = df_FoldSeek
	
	if len(df) == 0:
		st.write('No rows match the given column filters')
		# exit()
	
	return df

# ---------------------------------------------------- FILTER DATAFRAME ----------------------------------------------------

df = filter_dataframe(df_FoldSeek)
st.dataframe(df)
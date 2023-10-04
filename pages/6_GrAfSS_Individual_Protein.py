# -------------------------------------------------------- PACKAGES --------------------------------------------------------

import pandas as pd  # pip install pandas openpyxl
import streamlit as st  # pip install streamlit
import math

# ------------------------------------------------------------ LAYOUT ------------------------------------------------------------

st.set_page_config(layout="wide")
st.write('This page creates an interactive spreadsheet from GrAfSS output file for any single effector. '
		 'The spreadsheet is a list of known substructures found in the protein of interest using GrAfSS.')

# ------------------------------------------------------------ READ EXCEL ------------------------------------------------------------

df = pd.read_csv('Data Files/MpEffectorsBioinfo_Structural_KeyInfo.csv')

# ------------------------------------------------------------ VARIABLES and FUNCTIONS ------------------------------------------------------------

desired_file_path_GrAfSS = 'Data Files/GrAfSS/*.txt'

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

def get_file_of_interest_GrAfSS(protein):
	
	import glob
	import re
	
	for file_path in glob.glob(desired_file_path_GrAfSS):

		if re.search(protein + '_', file_path):
						
			return file_path.replace('\\', '/')

def GrAfSS_df(protein):

	import re
	import pandas as pd
	
	dict_proteinname = {0: 'MpC002', 100: 'MIF1', 92: 'Mp92a'}
	if protein in [0, 100, 92]:
		protein = dict_proteinname[protein]
	else:
		protein = 'Mp' + str(protein)	
	
	columns = ['substructure_ID','substructure_name','RMSD','hand']

	if get_file_of_interest_GrAfSS(protein):

		file = open(get_file_of_interest_GrAfSS(protein))
		lines = file.read().split('\n')

		positions = []

		for i in range(len(lines)):
			if 'user.SUM' in lines[i]:
				positions.append(i)

		lines_interest = []

		for i in range(positions[0], len(lines)):
			if ':' in lines[i]:
				line_split = re.split('\s{2,}', lines[i])
				del (line_split[0], line_split[2], line_split[-1])
				if i < positions[1]:
					line_split.append('RH')
				else:
					line_split.append('LH')
				lines_interest.append(line_split)

		df = pd.DataFrame(lines_interest, columns = columns)
	  
		return df
		
	else: 
		return None


ID_Number = st.selectbox("Select Effector:", df["ID_No"].unique())

if str(GrAfSS_df(ID_Number)) != 'None':
	
	df_GrAfSS = GrAfSS_df(ID_Number)

df_GrAfSS['RMSD'] = df_GrAfSS['RMSD'].astype(float)
df_GrAfSS.insert(2, 'substructure_name_count', df_GrAfSS['substructure_name'].map(df_GrAfSS['substructure_name'].value_counts().to_dict()))

columns_stringsearch = ['substructure_name']
columns_checkbox = ['RMSD']
columns_multiselect = ['hand']
columns_selectslider = ['substructure_name_count']
columns_slider_float = ['RMSD']
columns_all_filter = list_unique([columns_stringsearch, columns_checkbox, columns_multiselect, columns_selectslider, columns_slider_float])
	
# ------------------------------------------------------------ SIDEBAR ------------------------------------------------------------
		
def filter_dataframe(df):
	
	df_filter_dict = {}
	
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
				
				if column in columns_checkbox:
					df_filter_dict[f'{column}2'] = None
					threshold_checkbox = st.checkbox(f'Use threshold for {column}')
					if threshold_checkbox:
						user_text_input = st.text_input(f"Threshold for {column}")
						option = st.selectbox("Select threshold type:", [f'{column} above', f'{column} below'])
						set = st.checkbox(f'Activate threshold for {column}')
						df_filter_dict[f'{column}2'] = user_text_input
			
				if column in columns_multiselect:
					user_cat_input = st.multiselect(
													f"Values for {column}",
													df[column].unique(),
													default=list(df[column].unique()),
													)
					if len(user_cat_input) != df_GrAfSS[column].nunique():
						select_all = st.checkbox(f"Select all {column}", value=False)
						if select_all == True:# or reset_all == True:				
							st.write(f'All {column} values selected until this box is unchecked')
							user_cat_input = df_GrAfSS[column].unique()
					df_filter_dict[column] = user_cat_input

				if column in columns_selectslider:	
					# WILL ALWAYS BE UNFILTERED DF
					if df[column].nunique() > 10:
						_min = df[column].sort_values().min()
						_max = df[column].sort_values().max()
						user_num_cat_input = st.select_slider(
															  f"Values for {column}",
															  df_GrAfSS[column].sort_values(),
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
						if len(user_num_cat_input) != df_GrAfSS[column].nunique():
							select_all = st.checkbox(f"Select all {column}", value=False)
							if select_all == True:# or reset_all == True:				
								user_num_cat_input = df_GrAfSS[column].unique()
								st.write(f'All {column} values selected until this box is unchecked')
					df_filter_dict[column] = user_num_cat_input
				
				if column in columns_slider_float:
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
							if reset == True:# or reset_all == True:
								user_num_input = (_min, _max)
								st.write(f'Value range for {column} is maximal until this box is unchecked')
						df_filter_dict[column] = user_num_input
		
	# df_filter_dict

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
				df = df[df[column].isin(df_filter_dict[column])]
		if column in columns_slider_float:
			df = df[df[column].between(*df_filter_dict[column])]

	if reset_all:
		df = df_GrAfSS
	
	if len(df) == 0:
		st.write('No rows match the given column filters')
		# exit()
	
	return df

# ---------------------------------------------------- FILTER DATAFRAME ----------------------------------------------------

df = filter_dataframe(df_GrAfSS)
st.dataframe(df)
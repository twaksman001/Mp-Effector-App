# -------------------------------------------------------- PACKAGES --------------------------------------------------------

import pandas as pd
import streamlit as st
import math

# ------------------------------------------------------------ LAYOUT ------------------------------------------------------------

st.set_page_config(layout="wide")
st.write('This page loads an interactive spreadsheet summarising GrAfSS output for all proteins.')
with st.expander('Column definitions'):
	st.write('subtructure_count: total number of GrAfSS hits for the substructure')
	st.write('effector_count: total number of effectors with at least 1 GrAfSS hit for the substructure')
	st.write('effectors: list of effectors with at least 1 GrAfSS hit for the substructure')
	
# --------------------------------------------------------- READ EXCEL ---------------------------------------------------------

def convert_list_str_int_dfcolumn(list_):
	return list(map(int, list_.split(',')))
	
df_summary = pd.read_csv('Data Files/GrAfSS/GrAfSS_Summary.csv')
df_summary['effectors'] = df_summary['effectors'].apply(convert_list_str_int_dfcolumn)

columns_int = ['substructure_count', 'effector_count']
columns_list = ['effectors']
columns_string = ['substructure']

columns_all_filter = columns_list + columns_string + columns_int

# ------------------------------------------------------ OTHER VARIABLES ------------------------------------------------------

def row_list_dynamic_list(df, column):
	
	list_dynamic = []
	
	for list in df[column]:
		for i in list:
			if i not in list_dynamic:
				list_dynamic.append(i)
	
	list_dynamic.sort()

	return list_dynamic

def row_string_dynamic_list(df, column):
	
	list_dynamic = []
	
	for string in df[column]:
			if string not in list_dynamic:
				list_dynamic.append(string)
	
	# list_dynamic.sort()

	return list_dynamic

def textinput_string_index_list(df, column, textinput):
	
	index_lists = []
	
	for str in textinput.split(','):
		index_lists.append(df[df[column].str.contains(str)].index.values)
	
	index_list_unique = list(set([item for sublist in index_lists for item in sublist]))
		
	return index_list_unique

def textinput_list_index_list_any(df, column, textinput):
	
	index_list = []
	
	for i in range(len(df[column])):		
		counter = 0
		for j in df[column][i]:
			if counter == 0:
				# EFFECTOR IDS ARE EXACT MATCHES
				if any(x == str(j) for x in textinput.split(',')):
					index_list.append(i)
					counter += 1
		
	return index_list

def textinput_list_index_list_all(df, column, textinput):
	
	index_list = []
	
	for i in range(len(df[column])):		
		counter = 0
		for x in textinput.split(','):
			# EFFECTOR IDS ARE EXACT MATCHES
			if any(x == str(j) for j in df[column][i]):
				counter += 1
			else:
				break
		if counter == len(textinput.split(',')):
			index_list.append(i)
	
	return index_list

effector_ID_list = [0, 1, 2, 4, 5, 6, 7, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 
					40, 41, 42, 43, 44, 45, 46, 47, 49, 50, 51, 53, 54, 55, 57, 58, 60, 64, 65, 66, 70, 71, 72, 73, 74, 76, 77, 78, 79,
					81, 82, 85, 90, 91, 92, 93, 94, 95, 100]

# --------------------------------------------------------------- FUNCTIONS ---------------------------------------------------------------

def filter_dataframe(df):
	
	df_filter_dict = {}

	for column in columns_all_filter:
		df_filter_dict[column] = None
	
	with st.sidebar:
		st.header('Filters')
		reset_all = st.button('Reset all filters')
		for column in columns_all_filter:
			with st.expander(column):
				
				if column in columns_list:
					user_text_input = st.text_input(f'Text string(s) in {column} (case-sensitive ; if multiple, use only 1 comma to separate)')
					if user_text_input:# and reset_all == False:
						all_or_any = st.selectbox("All or any", ['All', 'Any'])
						if all_or_any == 'Any':
							df_filter_dict[column] = textinput_list_index_list_any(df, column, user_text_input)
						if all_or_any == 'All':
							df_filter_dict[column] = textinput_list_index_list_all(df, column, user_text_input)
				
				if column in columns_string:
					user_text_input = st.text_input(f'Text string(s) in {column} (case-sensitive ; if multiple, use only 1 comma to separate)')
					if user_text_input:# and reset_all == False:
						df_filter_dict[column] = textinput_string_index_list(df, column, user_text_input)
				
				if column in columns_int:
					_min = int(df[column].min())
					_max = int(df[column].max())
					step = math.ceil((_max - _min)/100)
					if step != 0:
						user_num_input = st.slider(
												   f'Values for {column}',
												   min_value=_min,
												   max_value=_max,
												   value=(_min, _max),
												   step=step
												   )
						if user_num_input != (_min, _max):
							reset = st.checkbox(f'Reset {column}', value=False)
							if reset == True:# or reset_all == True:
								user_num_input = (_min, _max)
								st.write(f'Value range for {column} is maximal until this box is unchecked')
						df_filter_dict[column] = user_num_input

	for column in columns_all_filter:
		if str(df_filter_dict[column]) == 'None':
			continue
		if column in columns_int:
			df = df[df[column].between(*df_filter_dict[column])]
		if column in columns_string:			
			df = df.loc[df_filter_dict[column]]
		if column in columns_list:
			df = df.loc[df_filter_dict[column]]
	
	if reset_all:
		df = df_summary
		for column in columns_all_filter:
			df_filter_dict[column] = None
	
	if len(df) == 0:
		st.write('No rows match the given column filters')
	
	return df

# ---------------------------------------------------- FILTER DATAFRAME ----------------------------------------------------

df = filter_dataframe(df_summary)
st.dataframe(df)
# -------------------------------------------------------- PACKAGES --------------------------------------------------------

import pandas as pd  # pip install pandas openpyxl
import plotly.express as px  # pip install plotly-express
import streamlit as st  # pip install streamlit
import math

# ------------------------------------------------------------ LAYOUT ------------------------------------------------------------

st.set_page_config(layout='wide')

st.write('This page loads an interactive spreadsheet summarising FoldSeek output for all proteins.')
st.write('For each protein, the spreadsheet from the "FoldSeek Individual Protein" page was filtered to evalue < 0.1.')

with st.expander('Column definitions'):
	
	st.subheader('Genera:')
	st.write('Unique genus list for each protein, in the order in which they appear in the filtered spreadsheet.')
	st.write('Taxon information only available for AlphaFold and PDB hits.')
	
	st.subheader('Top hits notes:')
	st.write('Summary information of top hits (max. 100 hits)')
	st.write('outside brackets: name summarizing top hits. Usually InterPro domain/function found in hits, sometimes name from Google UnProtein AI.')
	st.write('inside brackets: summary of statistics of top hits. evalue (max. 2 decimal places), hit number, taxon summary')
	st.write('Usually the vast majority of top hits are in one class.')
	st.write('Did not check very carefully to distinguish arthropod taxons, so small chance of significant presence of other arthropod types if only "insect" is noted.')

# ------------------------------------------------------------ READ EXCEL ------------------------------------------------------------

def convert_list_dfcolumn(df, column):
		
	list_newcolumn = []
	
	for i in range(len(df[column])):
		if df[column][i].startswith('['):
			list_new_cell = []
			for j in df[column][i][1:-1].split(', '):
				list_new_cell.append(j[1:-1])
			list_newcolumn.append(list_new_cell)
		else:
			list_newcolumn.append(['None'])
	
	return list_newcolumn

df_summary = pd.read_csv('Data Files/FoldSeek/FoldSeek_Summary.csv')
df_summary['Genera (not N/A)'] = convert_list_dfcolumn(df_summary, 'Genera (not N/A)')

columns_slider_int = ['Number Hits (e < 0.1)', 'AlphaFold', 'PDB', 'MGnify', 'GMGC']
columns_list = ['Genera (not N/A)']
columns_stringsearch = ['Top Hits Notes']

columns_all_filter = columns_slider_int + columns_list + columns_stringsearch

for column in columns_slider_int:
	df_summary[column] = df_summary[column].astype(int)

# ------------------------------------------------------------ VARIABLES and FUNCTIONS ------------------------------------------------------------

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
				if any(x in j for x in textinput.split(',')):
					index_list.append(i)
					counter += 1
		
	return index_list

def textinput_list_index_list_all(df, column, textinput):
	
	index_list = []
	
	for i in range(len(df[column])):		
		counter = 0
		for x in textinput.split(','):
			if any(x in j for j in df[column][i]):
				counter += 1
			else:
				break
		if counter == len(textinput.split(',')):
			index_list.append(i)
	
	return index_list

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
				
				if column in columns_slider_int:
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
				
				if column in columns_list:
					user_text_input = st.text_input(f'Text string(s) in {column} (case-sensitive ; if multiple, use only 1 comma to separate)')
					if user_text_input:# and reset_all == False:
						all_or_any = st.selectbox("All or any", ['All', 'Any'])
						if all_or_any == 'Any':
							df_filter_dict[column] = textinput_list_index_list_any(df, column, user_text_input)
						if all_or_any == 'All':
							df_filter_dict[column] = textinput_list_index_list_all(df, column, user_text_input)
				
				if column in columns_stringsearch:
					user_text_input = st.text_input(f'Text string(s) in {column} (case-sensitive ; if multiple, use only 1 comma to separate)')
					if user_text_input:# and reset_all == False:
						df_filter_dict[column] = textinput_string_index_list(df, column, user_text_input)
	
	for column in columns_all_filter:
		if str(df_filter_dict[column]) == 'None':
			continue
		if column in columns_slider_int:
			df = df[df[column].between(*df_filter_dict[column])]
		if column in columns_stringsearch:			
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

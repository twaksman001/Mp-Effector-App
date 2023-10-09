####################modules####################

import streamlit as st
import pandas as pd
import math

####################layout####################

st.set_page_config(layout="wide")

st.write('This page loads an interactive spreadsheet of key information for all effectors.')
st.write('Please see specific column definitions below.')

with st.expander('Column definitions'):
	st.write('sequence: primary sequence with predicted signal peptide removed')
	st.write('MW: molecular weight (kDa)')
	st.write('MSA_depth: approximate average (over all amino acids) depth of AlphaFold 2 sequence alignment')
	st.write('pTM: AlphaFold predicted template modelling score')
	st.write('SASA/length: solvent accessible surface area to length ratio of AlphaFold predicted structure')
	st.write('Pearson: Pearson correlation coefficient (r) between AlphaFold 2 '
			 'and OmegaFold per-residue confidence estimates')
	st.write('DALI_Z: DALI Z score comparing AlphaFold 2 and OmegaFold structures')
	st.write('pI: predicted isoelectric point (predicted via XTALPRED)')
	st.write('beta: protein contains beta-sheet secondary structure yes/no')
	st.write('modular: protein contains multiple domains/modules yes/no')
	
####################load dataframe####################

df = pd.read_csv(filepath_or_buffer='Data Files/MpEffectors_KeyInfo.txt', sep='\t')

with st.sidebar:
	st.header('filters')
	reset_all = st.button(label='reset all filters', key='reset_all_button')#, on_click=widgets_initial())

####################dataframe functions####################

def prepare_dataframe(df):

	columns_string = {'protein', 'sequence', 'beta', 'modular'}
	columns_float = {'MSA_Depth_log10', 'pTM', 'SASA/length', 'pI', 'DALI_Z', 'Pearson', 'MW'}
	columns_int = {'ID_Number', 'length', 'MSA_Depth'}
	
	for column in df.columns:
		if column in columns_string:
			df[column] = df[column].astype(str)
			df[column] = df[column].fillna(value = 'n/a')
		elif column in columns_float:
			df[column] = df[column].astype(float)
		elif column in columns_int:
			df[column] = df[column].astype(int)
	
	return df
			
def filter_dataframe(df):
		
	df_filter_dict = {}
	
	columns_stringsearch = {'sequence'}
	columns_checkbox = {}
	columns_multiselect = {'protein'}
	columns_selectslider = {'MSA_Depth', 'DALI_Z'}
	columns_slider_int = {'ID_Number', 'length'}
	columns_slider_float = {'MSA_Depth_log10', 'pTM', 'SASA/length', 'pI', 'Pearson', 'MW'}
	columns_all_filter = sorted(set.union(columns_stringsearch, columns_checkbox, columns_multiselect, 
									 columns_selectslider, columns_slider_int, columns_slider_float))

	for column in columns_all_filter:
		df_filter_dict[column] = None
	
	with st.sidebar:
		for column in columns_all_filter:
			with st.expander(column):	
				
				if column in columns_stringsearch:
					user_text_input = st.text_input(label=f'text string in {column} (case-sensitive ; if multiple, use only 1 comma to separate)', value='', key=f'text string in {column}')
					if user_text_input:
						df_filter_dict[column] = user_text_input
				
				if column in columns_checkbox:
					df_filter_dict[f'{column}2'] = None
					threshold_checkbox = st.checkbox(label=f'use threshold for {column}', value=False, key=f'use threshold for {column}')
					if threshold_checkbox:
						user_text_input = st.text_input(label=f'threshold for {column}', value='', key=f'threshold for {column}')
						option = st.selectbox(label='select threshold type:', options=[f'{column} above', f'{column} below'], index=0, key=f'threshold type for {column}')
						set_ = st.checkbox(label=f'activate threshold for {column}', value=False, key=f'activate threshold for {column}')
						df_filter_dict[f'{column}2'] = user_text_input

				if column in columns_multiselect:
					user_cat_input = st.multiselect(
													label=f'values for {column}',
													options=df[column].unique(),
													default=list(df[column].unique()),
													key=f'multiselect values for {column}'
													)
					if len(user_cat_input) != df[column].nunique():
						select_all = st.checkbox(label=f'select all {column}', value=False, key=f'multiselect select all for {column}')
						if select_all == True:				
							st.write(f'all {column} values selected until this box is unchecked')
							user_cat_input = df[column].unique()							
					df_filter_dict[column] = user_cat_input
				
				if column in columns_selectslider:	
					if df[column].nunique() > 10:
						_min = df[column].sort_values().min()
						_max = df[column].sort_values().max()
						user_num_cat_input = st.select_slider(
															  label=f'values for {column}',
															  options=df[column].sort_values(),
															  value=(_min, _max),
															  key=f'selectslider values for {column}'
															  )
						if user_num_cat_input != (_min, _max):
							reset = st.checkbox(label=f'reset {column}', value=False, key=f'reset selectslider values for {column}')
							if reset == True:
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
						if len(user_num_cat_input) != df[column].nunique():
							select_all = st.checkbox(label=f'select all {column}', value=False, key=f'selectslider multiselect select all for {column}')
							if select_all == True:
								user_num_cat_input = df[column].unique()
								st.write(f'all {column} values selected until this box is unchecked')				
					df_filter_dict[column] = user_num_cat_input

				if column in columns_slider_int:
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
						if user_num_input != (_min, _max):
							reset = st.checkbox(label=f'reset {column}', value=False, key=f'reset slider int values for {column}')
							if reset == True:
								user_num_input = (_min, _max)
								st.write(f'value range for {column} is maximal until this box is unchecked')
						df_filter_dict[column] = user_num_input
						
				if column in columns_slider_float:
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
							if reset == True:
								user_num_input = (_min, _max)
								st.write(f'value range for {column} is maximal until this box is unchecked')
						df_filter_dict[column] = user_num_input
		
	loc_list = []
	
	for column in columns_all_filter:
		if str(df_filter_dict[column]) == 'None':
			continue
		if column in columns_stringsearch:			
			loc_list.append(set(df[df[column].str.contains(df_filter_dict[column])].index))
		if column in columns_checkbox:
			if st.session_state[f'use threshold for {column}']:
				if st.session_state[f'activate threshold for {column}']:						
					if st.session_state[f'threshold type for {column}'] == f'{column} above':
						loc_list.append(set(df[df[column] >= float(df_filter_dict[f'{column}2'])].index))
					elif st.session_state[f'threshold type for {column}'] == f'{column} below':
						loc_list.append(set(df[df[column] <= float(df_filter_dict[f'{column}2'])].index))
		if column in columns_multiselect:
			loc_list.append(set(df[df[column].isin(df_filter_dict[column])].index))
		if column in columns_selectslider:
			if 'tuple' in str(type(df_filter_dict[column])):
				loc_list.append(set(df[df[column].between(*df_filter_dict[column])].index))
			if 'list' in str(type(df_filter_dict[column])):
				loc_list.append(set(df[df[column].isin(df_filter_dict[column])].index))
		if column in columns_slider_int:
			loc_list.append(set(df[df[column].between(*df_filter_dict[column])].index))
		if column in columns_slider_float:
			loc_list.append(set(df[df[column].between(*df_filter_dict[column])].index))
	
	df = df.loc[list(set.intersection(*loc_list))]
	
	if len(df) == 0:
		st.write('No rows match the given column filters')
		# exit()
		
	return df
	
def widgets_initial():
	
	for key in st.session_state.keys():
		del st.session_state[key]
		
	columns_stringsearch = {'sequence'}
	columns_checkbox = {}
	columns_multiselect = {'protein'}
	columns_selectslider = {'MSA_Depth', 'DALI_Z'}
	columns_slider_int = {'length'}
	columns_slider_float = {'MSA_Depth_log10', 'pTM', 'SASA/length', 'pI', 'Pearson', 'MW'}
	columns_all_filter = sorted(set.union(columns_stringsearch, columns_checkbox, columns_multiselect, 
									 columns_selectslider, columns_slider_int, columns_slider_float))
	
	with st.sidebar:
		for column in columns_all_filter:
			with st.expander(column):
				
				if column in columns_stringsearch:
					user_text_input = st.text_input(label=f'text string in {column} (case-sensitive ; if multiple, use only 1 comma to separate)', value='', key=f'text string in {column} -')		
				
				if column in columns_checkbox:
					threshold_checkbox = st.checkbox(label=f'use threshold for {column}', value=False, key=f'use threshold for {column} -')
					if threshold_checkbox:
						user_text_input = st.text_input(label=f'threshold for {column}', value='', key=f'threshold for {column} -')
						option = st.selectbox(label='select threshold type:', options=[f'{column} above', f'{column} below'], index=0, key=f'threshold type for {column} -')
						set_ = st.checkbox(label=f'activate threshold for {column}', value=False, key=f'activate threshold for {column} -')

				if column in columns_multiselect:
					user_cat_input = st.multiselect(
													label=f'values for {column}',
													options=df[column].unique(),
													default=list(df[column].unique()),
													key=f'multiselect values for {column} -'
													)
					if len(user_cat_input) != df[column].nunique():
						select_all = st.checkbox(label=f'select all {column}', value=False, key=f'multiselect select all for {column} -')
						if select_all == True:				
							st.write(f'all {column} values selected until this box is unchecked')
							user_cat_input = df[column].unique()
									
				if column in columns_selectslider:	
					if df[column].nunique() > 10:
						_min = df[column].sort_values().min()
						_max = df[column].sort_values().max()
						user_num_cat_input = st.select_slider(
															  label=f'values for {column}',
															  options=df[column].sort_values(),
															  value=(_min, _max),
															  key=f'selectslider values for {column} -'
															  )
						if user_num_cat_input != (_min, _max):
							reset = st.checkbox(label=f'reset {column}', value=False, key=f'reset selectslider values for {column} -')
							if reset == True:
								user_num_cat_input = (_min, _max)
								st.write(f'value range for {column} is maximal until this box is unchecked')
					else:
						user_num_cat_input = st.multiselect(
														label=f'values for {column}',
														options=df[column].unique(),
														default=list(df[column].unique()),
														key=f'selectslider multiselect values for {column} -'
														)
						if len(user_num_cat_input) != df[column].nunique():
							select_all = st.checkbox(label=f'select all {column}', value=False, key=f'selectslider multiselect select all for {column} -')
							if select_all == True:				
								user_num_cat_input = df[column].unique()
								st.write(f'all {column} values selected until this box is unchecked')
				
				if column in columns_slider_int:
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
												   key=f'slider int values for {column} -'
												   )
						if user_num_input != (_min, _max):
							reset = st.checkbox(label=f'reset {column}', value=False, key=f'reset slider int values for {column} -')
							if reset == True:
								user_num_input = (_min, _max)
								st.write(f'value range for {column} is maximal until this box is unchecked')
						
				if column in columns_slider_float:
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
												   key=f'slider float values for {column} -'
												   )
						if user_num_input != (_min, _max):
							reset = st.checkbox(label=f'reset {column}', value=False, key=f'reset slider float values for {column} -')
							if reset == True:
								user_num_input = (_min, _max)
								st.write(f'value range for {column} is maximal until this box is unchecked')
								
####################call functions####################	

if reset_all:
	widgets_initial()
	df = prepare_dataframe(df)
	st.dataframe(df)
else:
	df = prepare_dataframe(df)
	df = filter_dataframe(df)
	st.dataframe(df)
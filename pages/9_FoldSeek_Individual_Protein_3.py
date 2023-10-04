import streamlit as st  # pip install streamlit
import pandas as pd  # pip install pandas openpyxl
import math
#import sys

st.set_page_config(layout='wide')
st.write('This page creates an interactive spreadsheet from Foldseek output file for any single effector.')
with st.expander('Column definitions'):
	st.write('https://github.com/steineggerlab/Foldseek')

effectors_structural = ['Mp1', 'Mp2', 'Mp4', 'Mp5', 'Mp6', 'Mp7', 'Mp10', 'Mp11', 'Mp12', 'Mp14', 'Mp15', 'Mp16', 'Mp17', 'Mp19', 'Mp20', 
						'Mp21', 'Mp22', 'Mp23', 'Mp24', 'Mp28', 'Mp29', 'Mp30', 'Mp31', 'Mp32', 'Mp33', 'Mp35', 'Mp36', 'Mp37', 'Mp38', 
						'Mp39', 'Mp40', 'Mp41', 'Mp42', 'Mp43', 'Mp44', 'Mp45', 'Mp46', 'Mp47', 'Mp49', 'Mp50', 'Mp51', 'Mp53', 'Mp54', 
						'Mp55', 'Mp57', 'Mp58', 'Mp60', 'Mp64', 'Mp65', 'Mp66', 'Mp70', 'Mp71', 'Mp72', 'Mp73', 'Mp74', 'Mp76', 'Mp77', 
						'Mp78', 'Mp79', 'Mp81', 'Mp82', 'Mp85', 'Mp90', 'Mp91', 'Mp93', 'Mp94', 'Mp95', 'Mp92a', 'Mp92b', 'MpC002', 'MIF1', 'Mp67']
protein = st.selectbox(label='select Effector:', options=effectors_structural, index=0)

try:
	df_Foldseek = pd.read_csv(filepath_or_buffer='Data Files/Foldseek/'+protein+'_Foldseek.txt', sep='\t', index_col=0)
except Exception as e:
	st.write('no information available for this protein')

with st.sidebar:
	st.header('filters')
	reset_all = st.button(label='reset all filters', key='reset_all_button')#, on_click=widgets_initial())

def prepare_dataframe(df):

	columns_string = {'targetID', 'targetname', 'qaln', 'taln', 'tseq', 'genus', 'species', 'db', 'extra', 'tca'}
	columns_float = {'fident', 'evalue', 'prob'}
	columns_int = {'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'alnlen', 'mismatch', 'bits', 'genuscount', 'speciescount', 'qlen', 'tlen'}
	
	for column in df.columns:
		if column in columns_string:
			df[column] = df[column].astype(str)
			df[column] = df[column].fillna(value = 'n/a')
		elif column in columns_float:
			df[column] = df[column].astype(float)
		elif column in columns_int:
			df[column] = df[column].astype(int)
		elif column == 'taxid':
			df[column] = df[column].astype(str)
			df[column] = df[column].str.split('\.').str[0]
			df[column] = df[column].fillna(value = 'n/a')
	
	return df
			
def filter_dataframe(df):
		
	df_filter_dict = {}
	
	columns_stringsearch = {'targetID', 'targetname', 'qaln', 'taln', 'tseq', 'genus', 'species'}
	columns_checkbox = {'evalue', 'prob'}
	columns_multiselect = {'db'}
	columns_selectslider = {'gapopen', 'qstart', 'qend', 'tstart', 'tend'}
	columns_slider_int = {'alnlen', 'mismatch', 'bits', 'genuscount', 'speciescount'}
	columns_slider_float = {'fident', 'evalue', 'prob'}
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
					if len(user_cat_input) != df_Foldseek[column].nunique():
						select_all = st.checkbox(label=f'select all {column}', value=False, key=f'multiselect select all for {column}')
						if select_all == True:				
							st.write(f'all {column} values selected until this box is unchecked')
							user_cat_input = df_Foldseek[column].unique()							
					df_filter_dict[column] = user_cat_input
				
				if column in columns_selectslider:	
					if df[column].nunique() > 10:
						_min = df[column].sort_values().min()
						_max = df[column].sort_values().max()
						user_num_cat_input = st.select_slider(
															  label=f'values for {column}',
															  options=df_Foldseek[column].sort_values(),
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
						if len(user_num_cat_input) != df_Foldseek[column].nunique():
							select_all = st.checkbox(label=f'select all {column}', value=False, key=f'selectslider multiselect select all for {column}')
							if select_all == True:
								user_num_cat_input = df_Foldseek[column].unique()
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
		
	columns_stringsearch = {'targetID', 'targetname', 'qaln', 'taln', 'tseq', 'genus', 'species'}
	columns_checkbox = {'evalue'}
	columns_multiselect = {'db'}
	columns_selectslider = {'gapopen', 'qstart', 'qend', 'tstart', 'tend'}
	columns_slider_int = {'alnlen', 'mismatch', 'bits', 'genuscount', 'speciescount'}
	columns_slider_float = {'fident', 'evalue'}
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
													options=df_Foldseek[column].unique(),
													default=list(df_Foldseek[column].unique()),
													key=f'multiselect values for {column} -'
													)
					if len(user_cat_input) != df_Foldseek[column].nunique():
						select_all = st.checkbox(label=f'select all {column}', value=False, key=f'multiselect select all for {column} -')
						if select_all == True:				
							st.write(f'all {column} values selected until this box is unchecked')
							user_cat_input = df_Foldseek[column].unique()
									
				if column in columns_selectslider:	
					if df_Foldseek[column].nunique() > 10:
						_min = df_Foldseek[column].sort_values().min()
						_max = df_Foldseek[column].sort_values().max()
						user_num_cat_input = st.select_slider(
															  label=f'values for {column}',
															  options=df_Foldseek[column].sort_values(),
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
														options=df_Foldseek[column].unique(),
														default=list(df_Foldseek[column].unique()),
														key=f'selectslider multiselect values for {column} -'
														)
						if len(user_num_cat_input) != df_Foldseek[column].nunique():
							select_all = st.checkbox(label=f'select all {column}', value=False, key=f'selectslider multiselect select all for {column} -')
							if select_all == True:				
								user_num_cat_input = df_Foldseek[column].unique()
								st.write(f'all {column} values selected until this box is unchecked')
				
				if column in columns_slider_int:
					_min = int(df_Foldseek[column].min())
					_max = int(df_Foldseek[column].max())
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
					_min = float(df_Foldseek[column].min())
					_max = float(df_Foldseek[column].max())
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
	
	# st.dataframe(df_Foldseek)	

if reset_all:
	widgets_initial()
	df = prepare_dataframe(df_Foldseek)
	st.dataframe(df)
else:
	df = prepare_dataframe(df_Foldseek)
	df = filter_dataframe(df)
	st.dataframe(df)
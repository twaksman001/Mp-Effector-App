# -------------------------------------------------------- PACKAGES --------------------------------------------------------

import pandas as pd  # pip install pandas openpyxl
import streamlit as st  # pip install streamlit
import math
#import sys

# ------------------------------------------------------------ LAYOUT ------------------------------------------------------------

st.set_page_config(layout="wide")
st.write('This page creates an interactive spreadsheet from HHBlits output file for any single effector.')
st.write('If HHBlits has not been tested for a selected protein, or column filters result in no matching rows, error will occur.')
with st.expander('Column definitions'):
	st.write('targetID - ')
	st.write('probability - ')
	st.write('evalue - ')
	st.write('pvalue - ')
	st.write('score - ')
	st.write('fident - ')
	st.write('similarity - ')
	st.write('sum_probs - ')
	st.write('template_Neff - ')
	st.write('ss - ')
	st.write('alnlen - ')
	st.write('qstart - ')
	st.write('qend - ')
	st.write('tstart - ')
	st.write('tend - ')
	st.write('tlen - ')
	st.write('qaln - ')
	st.write('taln - ')
	st.write('taxID - ')
	st.write('extra - ')
	st.write('RepID - ')

# ------------------------------------------------------------ READ EXCEL ------------------------------------------------------------

df = pd.read_csv('C:/Users/TW43969/OneDrive - University of Dundee/Current/Research/Tom W Research Record/Bioinformatics/Structure/MpEffectorsBioinfo_Structural_KeyInfo_2.csv')

# ------------------------------------------------------------ VARIABLES ------------------------------------------------------------

desired_file_path_HHBlits = 'C:/Users/TW43969/OneDrive - University of Dundee/Current/Research/Tom W Research Record/Bioinformatics/Structure/Sequence Similarity/HHBlits/*.txt'

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

def get_file_of_interest_HHBlits(protein):
	
	import glob
	import re
	
	for file_path in glob.glob(desired_file_path_HHBlits):

		if re.search(protein + '_', file_path) != None:
						
			return file_path.replace('\\', '/')
	
def HHBlits_df(protein):

	import re
	import pandas as pd
		
	dict_proteinname = {0: 'MpC002', 100: 'MIF1', 92: 'Mp92a'}
	if protein in [0, 100, 92]:
		protein = dict_proteinname[protein]
	else:
		protein = 'Mp' + str(protein)
	
	file = get_file_of_interest_HHBlits(protein)
			
	if file:	
		file = open(file).read()
		
		columns = ['targetID','targetname','probability','evalue','pvalue','score','fident','similarity','sum_probs','template_Neff',
				   'ss','alnlen','qstart','qend','tstart','tend','tlen','qaln','taln','taxID','genus','species','extra','RepID']

		columns_dict = {}
		for i in range(len(columns)):
			columns_dict[columns[i]] = []

		hits_list = file[file.find('No Hit'):file.find('\nNo ')].split('\n')[1:-1]
		hits_list_titles = ['probability','evalue','pvalue','score','ss','alnlen','qstart','qend','tstart','tend','tlen']

		alignments = file.split('\n\nNo ')[1:]
		alignments_titles = ['targetID','targetname','genus','species','taxID','extra','RepID','fident','similarity','sum_probs','template_Neff',
							'qaln','taln']	

		for i in range(len(hits_list)):

			######PROCESS HITS LIST######

			hit_info_start = int(str(re.search('\s\d{2,3}\.\d\s{1,}\d', hits_list[i]))[24:26])
			hit_info_end = hits_list[i].rfind('(')
			hit_info_list = re.split('\s{1,}', hits_list[i][hit_info_start+1:hit_info_end])
			hit_info_list = (hit_info_list[:6] + hit_info_list[6].split('-') + hit_info_list[7].split('-') +
							 [hits_list[i][hits_list[i].rfind('(')+1:-1]]) # type == list

			for j in range(len(hits_list_titles)):
				columns_dict[hits_list_titles[j]].append(hit_info_list[j])

			######PROCESS ALIGNMENTS LIST######

			alignment_lines = alignments[i].split('\n')   
			alignment_info_list = []

			# targetID
			alignment_info_list.append(alignment_lines[1][11:alignment_lines[1].find(' ')])
			# targetname
			alignment_info_list.append(alignment_lines[1][alignment_lines[1].find(' '):alignment_lines[1].find('n=')])

			### taxon ###
			taxon_list = alignment_lines[1][alignment_lines[1].find('Tax=')+4:].split(' ')
			# genus
			alignment_info_list.append(taxon_list[0])
			# species
			if re.search('\d', taxon_list[1]):
				alignment_info_list.append('N/A')
			else:
				alignment_info_list.append(taxon_list[1])
			# taxID + extra + RepID
			for j in range(len(taxon_list)):
				if 'TaxID' in taxon_list[j]:
					alignment_info_list.append(taxon_list[j][6:])
					if j >= 3:
						alignment_info_list.append(' '.join(taxon_list[2:j]))
					if j <= 2:
						alignment_info_list.append('N/A')
				if 'RepID' in taxon_list[j]:
					alignment_info_list.append(taxon_list[j][6:])

			# other statistics #
			other_stat_list = re.split('\s{1,}', alignment_lines[2])
			other_stat_list = [int(other_stat_list[4][11:-1])/100, other_stat_list[5][11:], other_stat_list[6][10:], other_stat_list[7][14:]]
			alignment_info_list += other_stat_list

			### sequences ###
			query_sequence_aln = ''
			target_sequence_aln = ''

			for line in alignment_lines:
				if '~' not in line:
					if line.startswith('Q'):
						query_sequence_aln += re.split('\s{1,}', line)[3]
					if line.startswith('T'):
						target_sequence_aln += re.split('\s{1,}', line)[3]

			alignment_info_list += [query_sequence_aln, target_sequence_aln]

			for j in range(len(alignments_titles)):
				columns_dict[alignments_titles[j]].append(alignment_info_list[j])

		df = pd.DataFrame.from_dict(columns_dict)

		return df

	else:
		return None

def HHBlits_df_present(protein):
	
	import pandas as pd
	import numpy as np
	
	df = HHBlits_df(protein)
		
	column_list_strint = ['alnlen','qstart','qend','tstart','tend','tlen']
	column_list_strfloat = ['probability','evalue','pvalue','score','fident','similarity','sum_probs','template_Neff','ss']
	
	if str(df) != 'None':
					
		for column in column_list_strint:
			
			df[column] = df[column].astype(np.int64)
		
		for column in column_list_strfloat:
			
			df[column] = df[column].astype(float)
				
		df['genus'] = df['genus'] + ' (' + df['genus'].map(df['genus'].value_counts().to_dict()).astype(str) + ')'
		df['species'] = df['species'] + ' (' + df['species'].map(df['species'].value_counts().to_dict()).astype(str) + ')'
		df.insert(21, 'genuscount', df['genus'].map(df['genus'].value_counts().to_dict()))
		df.insert(23, 'speciescount', df['species'].map(df['species'].value_counts().to_dict()))
				
		df = df.sort_values(by='probability', axis=0, ascending=False)
		
		return df
	
ID_Number = st.selectbox("Select Effector:", df["ID_No"].unique())

if str(HHBlits_df_present(ID_Number)) != 'None':
	
	df_HHBlits = HHBlits_df_present(ID_Number)
	
# ------------------------------------------------------------ SIDEBAR ------------------------------------------------------------
		
def filter_dataframe(df):
	
	df_filter_dict = {}
	
	columns_stringsearch = ['targetID','targetname','qaln','taln','taxID','genus','species','extra','RepID']
	columns_checkbox = ['probability','evalue','pvalue']
	columns_selectslider = ['qstart','qend','tstart','tend']
	columns_slider_int = ['alnlen','tlen','genuscount','speciescount']
	columns_slider_float = ['probability','evalue','pvalue','score','fident','similarity','sum_probs','template_Neff','ss']
	columns_all_filter = list_unique([columns_stringsearch, columns_checkbox, 
									 columns_selectslider, columns_slider_int, columns_slider_float])
		
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
				
				# if column in columns_multiselect:
					# # cat_values = df[column].unique()
					# user_cat_input = st.multiselect(
													# f"Values for {column}",
													# df[column].unique(),
													# default=list(df[column].unique()),
													# )
					# # st.write(len(user_cat_input), len(cat_values))
					# if len(user_cat_input) != df_HHBlits[column].nunique():
						# select_all = st.checkbox(f"Select all {column}", value=False)
						# if select_all == True:# or reset_all == True:				
							# st.write(f'All {column} values selected until this box is unchecked')
							# user_cat_input = df_HHBlits[column].unique()
											 # # st.multiselect(
															# # f"Values for {column}",
															# # df[column].unique(),
															# # )					
					# df_filter_dict[column] = user_cat_input
				
				if column in columns_selectslider:	
					# if len(df) != 0:
					if df[column].nunique() > 10:
						_min = df[column].sort_values().min()
						_max = df[column].sort_values().max()
						user_num_cat_input = st.select_slider(
															  f"Values for {column}",
															  df_HHBlits[column].sort_values(),
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
						if len(user_num_cat_input) != df_HHBlits[column].nunique():
							select_all = st.checkbox(f"Select all {column}", value=False)
							if select_all == True:# or reset_all == True:				
								user_num_cat_input = df_HHBlits[column].unique()
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
	
	df_filter_dict

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
		# if column in columns_multiselect:
			# df = df[df[column].isin(df_filter_dict[column])]
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
		df = df_HHBlits
	
	if len(df) == 0:
		st.write('No rows match the given column filters')
		# exit()
	
	return df

# ---------------------------------------------------- FILTER DATAFRAME ----------------------------------------------------

df = filter_dataframe(df_HHBlits)
st.dataframe(df)
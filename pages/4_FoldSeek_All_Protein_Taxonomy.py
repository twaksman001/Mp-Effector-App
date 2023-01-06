# -------------------------------------------------------- PACKAGES --------------------------------------------------------

# from pandas.api.types import (
	# is_categorical_dtype,
	# # is_datetime64_any_dtype,
	# is_numeric_dtype,
	# is_object_dtype,
# )
import pandas as pd
import streamlit as st

# ------------------------------------------------------------ LAYOUT ------------------------------------------------------------

st.set_page_config(layout="wide")
st.write('This page loads an interactive spreadsheet of compiled taxonomic information from FoldSeek output files for all effectors.')
with st.expander('Column definitions'):
	st.write('genus_count: total number of FoldSeek hits belonging to genus')
	st.write('effector_count: total number of effectors with at least 1 FoldSeek hit belonging to genus')
	st.write('effectors: list of effectors with at least 1 FoldSeek hit belonging to genus')
	
# --------------------------------------------------------- READ EXCEL ---------------------------------------------------------

def convert_list_str_int_dfcolumn(list_):
	return list(map(int, list_.split(',')))
	
df_taxonomy = pd.read_csv('Data Files/FoldSeek/FoldSeek_Taxonomy.csv')
df_taxonomy = df_taxonomy.drop(columns = 'Unnamed: 0')
df_taxonomy['effectors'] = df_taxonomy['effectors'].apply(convert_list_str_int_dfcolumn)
df_taxonomy['genus'] = df_taxonomy['genus'] + ' (' + df_taxonomy['genus_count'].astype(str) + ')'

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

def textinput_row_string_dynamic_list(df, column, textinput):
	
	list_dynamic = []
	
	for string in df[column]:		
		if any(x in string for x in textinput.split(',')):
			list_dynamic.append(string)
	
	return list_dynamic

effector_ID_list = [0, 1, 2, 4, 5, 6, 7, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 
					40, 41, 42, 43, 44, 45, 46, 47, 49, 50, 51, 53, 54, 55, 57, 58, 60, 64, 65, 66, 70, 71, 72, 73, 74, 76, 77, 78, 79,
					81, 82, 85, 90, 91, 92, 93, 94, 95, 100]

# --------------------------------------------------------------- FUNCTIONS ---------------------------------------------------------------

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

def filter_dataframe(df):

	with st.sidebar:
		for column in df.columns:
			with st.expander(column):
				
				if all(type(k) == list for k in df[column]):
					# select_all = st.checkbox(f"Select all {column}")#, value=True)
					# if select_all:
						# user_list_input = st.multiselect(
														 # f"Values for {column}",
														 # effector_ID_list,
														 # default=effector_ID_list
														 # )
					# else:
						# user_list_input = st.multiselect(
														 # f"Values for {column}",
														 # effector_ID_list_dynamic,
														 # default=effector_ID_list_dynamic
														 # )
					user_list_input = st.multiselect(
													 f"Values for {column}",
													 row_list_dynamic_list(df, column),
													 default=row_list_dynamic_list(df, column)
													 )
					df = df.loc[df_row_index_list_cond(df, column, user_list_input)]
					
					
				# elif df[column].nunique() < 10:
					# user_cat_input = st.multiselect(
													# f"Values for {column}",
													# df[column].unique(),
													# default=list(df[column].unique()),
													# )
					# df = df[df[column].isin(user_cat_input)]
					
				elif all(type(k) == int for k in df[column]):
					_min = int(df[column].min())
					_max = int(df[column].max())
					step = (_max - _min)/100 # not integer
					reset = st.checkbox(f"Reset {column}", value=True)
					if reset:
						user_num_input = st.slider(
												   f"Values for {column}",
												   min_value=_min,
												   max_value=_max,
												   value=(_min, _max),
												   step=1, # not working!
												   )
					else:
						user_num_input = st.slider(
												   f"Values for {column}",
												   min_value=_min,
												   max_value=_max,
												   value=(int(_min + 10*step), int(_max - 10*step)),
												   step=1, # not working!
												   )
					df = df[df[column].between(*user_num_input)]
					
				# elif all(type(k) == float for k in df[column]):
					# _min = float(df[column].min())
					# _max = float(df[column].max())
					# step = (_max - _min) // 100
					# user_num_input = st.slider(
											   # f"Values for {column}",
											   # min_value=_min,
											   # max_value=_max,
											   # value=(_min, _max),
											   # step=step,
											   # )
					# df = df[df[column].between(*user_num_input)]			
				
				elif all(type(k) == str for k in df[column]):
					# select_all = st.checkbox(f"Select all {column}")#, value=True)
					# if select_all:
						# user_cat_input = st.multiselect(
														# f"Values for {column}",
														# df[column].unique(),
														# default=df[column].unique(),
														# )
					# else:
						# user_cat_input = st.multiselect(
														# f"Values for {column}",
														# df[column].unique(),
														# default=df[column].head(21),
														# )
				
					# user_cat_input = st.multiselect(
													 # f"Values for {column}",
													 # row_string_dynamic_list(df, column),
													 # default=row_string_dynamic_list(df, column)
													 # )
					# df = df[df[column].isin(user_cat_input)]
				
				# else:
					user_text_input = st.text_input(
													f"Key letters in {column} (if multiple genera, use only 1 comma to separate)",
													)
					if user_text_input:
						list_dynamic = textinput_row_string_dynamic_list(df, column, user_text_input)
						df = df.loc[df_row_index_list_cond(df, column, list_dynamic)]
							
	return df

# ---------------------------------------------------- FILTER DATAFRAME ----------------------------------------------------

df = filter_dataframe(df_taxonomy)
st.dataframe(df)
# st.dataframe(filter_dataframe(df_taxonomy))

# st.write(type(df_taxonomy['effectors'][1]))
# st.write(type(df_taxonomy['effectors'][1][1]))
# st.write(type(user_list_input))
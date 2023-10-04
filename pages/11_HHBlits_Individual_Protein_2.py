import streamlit as st  # pip install streamlit
import pandas as pd  # pip install pandas openpyxl
import math
#import sys

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

effectors_structural = ['Mp1', 'Mp2', 'Mp4', 'Mp5', 'Mp6', 'Mp7', 'Mp10', 'Mp11', 'Mp12', 'Mp14', 'Mp15', 'Mp16', 'Mp17', 'Mp19', 'Mp20', 
						'Mp21', 'Mp22', 'Mp23', 'Mp24', 'Mp28', 'Mp29', 'Mp30', 'Mp31', 'Mp32', 'Mp33', 'Mp35', 'Mp36', 'Mp37', 'Mp38', 
						'Mp39', 'Mp40', 'Mp41', 'Mp42', 'Mp43', 'Mp44', 'Mp45', 'Mp46', 'Mp47', 'Mp49', 'Mp50', 'Mp51', 'Mp53', 'Mp54', 
						'Mp55', 'Mp57', 'Mp58', 'Mp60', 'Mp64', 'Mp65', 'Mp66', 'Mp70', 'Mp71', 'Mp72', 'Mp73', 'Mp74', 'Mp76', 'Mp77', 
						'Mp78', 'Mp79', 'Mp81', 'Mp82', 'Mp85', 'Mp90', 'Mp91', 'Mp93', 'Mp94', 'Mp95', 'Mp92a', 'Mp92b', 'MpC002', 'MIF1', 'Mp67']
protein = st.selectbox(label='select effector:', options=effectors_structural, index=0)

def get_file_of_interest(protein, path='Data Files/HHBlits/*', string='_HHBlits.txt'):
    
    import glob
    import re
    
    desired_file = None
    
    for file_path in glob.glob(path):
        if re.search(protein+string, file_path):                        
            desired_file = file_path.replace('\\', '/')
    
    if desired_file == None:
        raise ValueError(whoami(), 'fail to get file of interest')
    else:
        return desired_file

def HHBlits_df(protein):

    import re
    import pandas as pd
    
    try:
        file = get_file_of_interest(protein)
    except Exception as e:
        raise ValueError(whoami(), e)
    else:        
        file = open(file).read()
        
        columns = ['targetID','targetname','probability','evalue','pvalue','score','fident','similarity','sum_probs','template_Neff',
                   'ss','alnlen','qstart','qend','tstart','tend','tlen','qaln','taln','taxon','TaxID']

        columns_dict = {}
        for i in range(len(columns)):
            columns_dict[columns[i]] = []

        hits_list = file[file.find('No Hit'):file.find('\nNo ')].split('\n')[1:-1]
        hits_list_titles = ['probability','evalue','pvalue','score','ss','alnlen','qstart','qend','tstart','tend','tlen']

        alignments = file.split('\n\nNo ')[1:]
        alignments_titles = ['targetID','targetname','taxon','TaxID','fident','similarity','sum_probs','template_Neff',
                             'qaln','taln']
        
        assert len(hits_list)==len(alignments)

        for i in range(len(hits_list)):

            """PROCESS HITS LIST"""

            hit_info_start = re.search('\s\d{2,3}\.\d\s{1,}\d', hits_list[i]).span()[0]
            hit_info_end = hits_list[i].rfind('(')
            hit_info_list = re.split('\s{1,}', hits_list[i][hit_info_start+1:hit_info_end])
            hit_info_list = (hit_info_list[:6] + hit_info_list[6].split('-') + hit_info_list[7].split('-') +
                             [hits_list[i][hits_list[i].rfind('(')+1:-1]]) # type == list

            for j in range(len(hits_list_titles)):
                columns_dict[hits_list_titles[j]].append(hit_info_list[j])

            """PROCESS ALIGNMENTS LIST"""

            alignment_lines = alignments[i].split('\n')   

            targetID = lambda x: x[11:alignment_lines[1].find(' ')]
            targetname = lambda x: x[alignment_lines[1].find(' '):alignment_lines[1].find('n=')]
            taxon = lambda x: x[alignment_lines[1].find('Tax=')+4:alignment_lines[1].find(' TaxID')]
            TaxID = lambda x: x[alignment_lines[1].find(' TaxID=')+7:alignment_lines[1].find(' RepID')]

            # other statistics #
            other_stat_list = re.split('\s{1,}', alignment_lines[2])
            other_stat_list = [int(other_stat_list[4][11:-1])/100, other_stat_list[5][11:], other_stat_list[6][10:], other_stat_list[7][14:]]

            # sequences #
            query_sequence_aln = ''
            target_sequence_aln = ''
            for line in alignment_lines:
                if '~' not in line:
                    if line.startswith('Q'):
                        query_sequence_aln += re.split('\s{1,}', line)[3]
                    if line.startswith('T'):
                        target_sequence_aln += re.split('\s{1,}', line)[3]

            alignment_info_list = [targetID(alignment_lines[1]), targetname(alignment_lines[1]), taxon(alignment_lines[1]), TaxID(alignment_lines[1]), 
                                   *other_stat_list, query_sequence_aln, target_sequence_aln]

            for j in range(len(alignments_titles)):
                columns_dict[alignments_titles[j]].append(alignment_info_list[j])

        df = pd.DataFrame.from_dict(data=columns_dict)

        function = lambda x,y: df['taxon'].str.split(' ').str.slice(x,y).str.join(' ')
        df['genus'], df['species'], df['extra'] = function(0,1), function(1,2), function(2,100)
        
        return df

def HHBlits_df_present(protein):
    
    import pandas as pd
    import numpy as np
    
    column_list_strint = ['alnlen','qstart','qend','tstart','tend','tlen']
    column_list_strfloat = ['probability','evalue','pvalue','score','fident','similarity','sum_probs','template_Neff','ss']
    
    try:
        df = HHBlits_df(protein)
    except Exception as e:
        raise ValueError(whoami(), 'fail to get df', e)
    else:     
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

try:
	df_HHBlits = HHBlits_df_present(protein)
except Exception as e:
	st.write('no information available for this protein')
	
with st.sidebar:
	st.header('filters')
	reset_all = st.button(label='reset all filters', key='reset_all_button')#, on_click=widgets_initial())

def filter_dataframe(df):
		
	df_filter_dict = {}
	
	columns_stringsearch = {'targetID','targetname','qaln','taln','taxID','genus','species','extra'}
	columns_checkbox = {'probability','evalue','pvalue'}
	columns_selectslider = {'qstart','qend','tstart','tend'}
	columns_slider_int = {'alnlen','tlen','genuscount','speciescount'}
	columns_slider_float = {'probability','evalue','pvalue','score','fident','similarity','sum_probs','template_Neff','ss'}
	columns_all_filter = sorted(set.union(columns_stringsearch, columns_checkbox, columns_selectslider, columns_slider_int, columns_slider_float))

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

				if column in columns_selectslider:	
					if df[column].nunique() > 10:
						_min = df[column].sort_values().min()
						_max = df[column].sort_values().max()
						user_num_cat_input = st.select_slider(
															  label=f'values for {column}',
															  options=df_HHBlits[column].sort_values(),
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
						if len(user_num_cat_input) != df_HHBlits[column].nunique():
							select_all = st.checkbox(label=f'select all {column}', value=False, key=f'selectslider multiselect select all for {column}')
							if select_all == True:
								user_num_cat_input = df_HHBlits[column].unique()
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
		
	columns_stringsearch = {'targetID','targetname','qaln','taln','taxID','genus','species','extra'}
	columns_checkbox = {'probability','evalue','pvalue'}
	columns_multiselect = {}
	columns_selectslider = {'qstart','qend','tstart','tend'}
	columns_slider_int = {'alnlen','tlen','genuscount','speciescount'}
	columns_slider_float = {'probability','evalue','pvalue','score','fident','similarity','sum_probs','template_Neff','ss'}
	columns_all_filter = sorted(set.union(columns_stringsearch, columns_checkbox, columns_multiselect,
									 columns_selectslider, columns_slider_int, columns_slider_float))
	
	with st.sidebar:
		for column in columns_all_filter:
			with st.expander(column):
				
				if column in columns_stringsearch:
					user_text_input = st.text_input(label=f'text string in {column} (case-sensitive ; if multiple, use only 1 comma to separate)', value='', key=f'text string in {column}')		
				
				if column in columns_checkbox:
					threshold_checkbox = st.checkbox(label=f'use threshold for {column}', value=False, key=f'use threshold for {column}')
					if threshold_checkbox:
						user_text_input = st.text_input(label=f'threshold for {column}', value='', key=f'threshold for {column}')
						option = st.selectbox(label='select threshold type:', options=[f'{column} above', f'{column} below'], index=0, key=f'threshold type for {column}')
						set_ = st.checkbox(label=f'activate threshold for {column}', value=False, key=f'activate threshold for {column}')
						
				if column in columns_selectslider:	
					if df_HHBlits[column].nunique() > 10:
						_min = df_HHBlits[column].sort_values().min()
						_max = df_HHBlits[column].sort_values().max()
						user_num_cat_input = st.select_slider(
															  label=f'values for {column}',
															  options=df_HHBlits[column].sort_values(),
															  value=(_min, _max),
															  key=f'selectslider values for {column}'
															  )
						if user_num_cat_input != (_min, _max):
							reset = st.checkbox(label=f'reset {column}', value=False, key=f'reset selectslider values for {column}')
							if reset == True:
								user_num_cat_input = (_min, _max)
								st.write(f'value range for {column} is maximal until this box is unchecked')
					else:
						user_num_cat_input = st.multiselect(
														label=f'values for {column}',
														options=df_HHBlits[column].unique(),
														default=list(df_HHBlits[column].unique()),
														key=f'selectslider multiselect values for {column}'
														)
						if len(user_num_cat_input) != df_HHBlits[column].nunique():
							select_all = st.checkbox(label=f'select all {column}', value=False, key=f'selectslider multiselect select all for {column}')
							if select_all == True:				
								user_num_cat_input = df_HHBlits[column].unique()
								st.write(f'all {column} values selected until this box is unchecked')
				
				if column in columns_slider_int:
					_min = int(df_HHBlits[column].min())
					_max = int(df_HHBlits[column].max())
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
						
				if column in columns_slider_float:
					_min = float(df_HHBlits[column].min())
					_max = float(df_HHBlits[column].max())
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
	
	st.dataframe(df_HHBlits)	

if reset_all:
	widgets_initial()	
else:
	df = filter_dataframe(df_HHBlits)
	st.dataframe(df)
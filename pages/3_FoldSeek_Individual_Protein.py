# -------------------------------------------------------- PACKAGES --------------------------------------------------------

import pandas as pd  # pip install pandas openpyxl
import plotly.express as px  # pip install plotly-express
import streamlit as st  # pip install streamlit

# ------------------------------------------------------------ LAYOUT ------------------------------------------------------------

st.set_page_config(layout="wide")
st.write('This page creates an interactive spreadsheet from FoldSeek output file for any single effector.')
with st.expander('Column definitions'):
	st.write('https://github.com/steineggerlab/foldseek')

# ------------------------------------------------------------ READ EXCEL ------------------------------------------------------------

df = pd.read_csv('Data Files/MpEffectorsBioinfo_Structural_KeyInfo.csv')
df_taxonomy = pd.read_csv('Data Files/FoldSeek/FoldSeek_Taxonomy.csv')

# ------------------------------------------------------------ VARIABLES ------------------------------------------------------------

desired_file_path_FoldSeek = 'Data Files/FoldSeek/*.txt'
column_list_strint = ['alnlen','mismatch','gapopen','qstart','qend','tstart','tend','qlen','tlen', 'bits']
column_list_strfloat = ['fident', 'evalue']

# ------------------------------------------------------------ FUNCTIONS ------------------------------------------------------------

def get_file_of_interest_FoldSeek(protein):
	
	import glob
	import re
	
	for file_path in glob.glob(desired_file_path_FoldSeek):

		if re.search(protein + '_', file_path):
						
			return file_path.replace('\\', '/')
			
# def FoldSeek_df(protein):
	
	# import re
	# import pandas as pd
	
	# if protein == 0:
		# protein = 'MpC002'
	# elif protein == 100:
		# protein = 'MIF1'
	# else:
		# protein = 'Mp' + str(protein)
	
	# column_titles = ['query','target','fident','alnlen','mismatch','gapopen','qstart','qend','tstart',
					 # 'tend','evalue','bits','qlen','tlen','qaln','taln','tca','tseq','taxid','genus','species']		
	
	# file = open(get_file_of_interest_FoldSeek(protein))
	# lines = file.read().split('\n')
	
	# columns_dict = {}

	# for i in range(len(column_titles)):

		# columns_dict[column_titles[i]] = []

	# columns_dict['extra'] = []

	# for j in range(len(lines)):

		# line_split = re.split('\s{1,}', lines[j][lines[j].find('job.pdb'):])

		# if len(line_split) <= len(column_titles):

			# for i in range(len(column_titles)):

				# if i < len(line_split):

					# columns_dict[column_titles[i]].append(line_split[i])

				# else:

					# columns_dict[column_titles[i]].append('N/A')

			# columns_dict['extra'].append('0')

		# else:

			# for i in range(len(line_split)):

				# if i < len(column_titles):

					# columns_dict[column_titles[i]].append(line_split[i])

				# else:

					# columns_dict['extra'].append(line_split[i:])

					# break

	# df = pd.DataFrame.from_dict(columns_dict)
	
	# return df

def FoldSeek_df(protein):
	
	import re
	import pandas as pd
	
	if protein == 0:
		protein = 'MpC002'
	elif protein == 100:
		protein = 'MIF1'
	else:
		protein = 'Mp' + str(protein)
	
	column_titles = ['query','target','fident','alnlen','mismatch','gapopen','qstart','qend','tstart',
					 'tend','evalue','bits','qlen','tlen','qaln','taln','tca','tseq','taxid','genus','species']	

	if get_file_of_interest_FoldSeek(protein):
		
		file = open(get_file_of_interest_FoldSeek(protein))
		
		lines = file.read().split('\n')
	
		columns_dict = {}

		for i in range(len(column_titles)):

			columns_dict[column_titles[i]] = []

		columns_dict['extra'] = []

		for j in range(len(lines)):

			line_split = re.split('\s{1,}', lines[j][lines[j].find('job.pdb'):])

			if len(line_split) <= len(column_titles):

				for i in range(len(column_titles)):

					if i < len(line_split):

						columns_dict[column_titles[i]].append(line_split[i])

					else:

						columns_dict[column_titles[i]].append('N/A')

				columns_dict['extra'].append('0')

			else:

				for i in range(len(line_split)):

					if i < len(column_titles):

						columns_dict[column_titles[i]].append(line_split[i])

					else:

						columns_dict['extra'].append(line_split[i:])

						break

		df = pd.DataFrame.from_dict(columns_dict)
		
		return df
			
	else:
		
		st.write('No FoldSeek hits for' + protein)

def df_present_streamlit(protein):
	
	df = FoldSeek_df(protein)
	
	df = df.drop(index = df.index[-1], columns = 'query')
	
	import numpy as np
	
	for column in column_list_strint:
		
		df[column] = df[column].astype(np.int64)
	
	for column in column_list_strfloat:
		
		df[column] = df[column].astype(float)
	
	df['genus'] = df['genus'] + ' (' + df['genus'].map(df['genus'].value_counts().to_dict()).astype(str) + ')'
	df['species'] = df['species'] + ' (' + df['species'].map(df['species'].value_counts().to_dict()).astype(str) + ')'
	df.insert(19, 'genuscount', df['genus'].map(df['genus'].value_counts().to_dict()))
	df.insert(21, 'speciescount', df['species'].map(df['species'].value_counts().to_dict()))
	df = df.sort_values(by='bits', axis=0, ascending=False)
	
	return df


ID_Number = st.selectbox("Select Effector:", df["ID_No"].unique())

df_FoldSeek = df_present_streamlit(ID_Number)

# st.dataframe(df_FoldSeek.style.format({'fident': '{:.1f}'}))

# ------------------------------------------------------------ SIDEBAR ------------------------------------------------------------

with st.sidebar:
	
	st.header('Filters')

	# Select fraction of identical matches
	with st.expander('Fraction of identical matches'):
		reset_fident = st.checkbox("Reset fident", value=True)
		if reset_fident:
			min_fident, max_fident = st.slider("fident", float(df_FoldSeek["fident"].min()), float(df_FoldSeek["fident"].max()), (float(df_FoldSeek["fident"].min()), float(df_FoldSeek["fident"].max())), step=0.1)
		else:
			min_fident, max_fident = st.slider("fident", float(df_FoldSeek["fident"].min()), float(df_FoldSeek["fident"].max()), (5.0, 10.0), step=0.1)

	# Select alignment length
	with st.expander('Alignment length'):
		reset_alnlen = st.checkbox("Reset alnlen", value=True)
		if reset_alnlen:
			min_alnlen, max_alnlen = st.slider("alnlen", int(df_FoldSeek["alnlen"].min()), int(df_FoldSeek["alnlen"].max()), (int(df_FoldSeek["alnlen"].min()), int(df_FoldSeek["alnlen"].max())), step=1)
		else:
			min_alnlen, max_alnlen = st.slider("alnlen", int(df_FoldSeek["alnlen"].min()), int(df_FoldSeek["alnlen"].max()), (int(df_FoldSeek["alnlen"].min())+10, int(df_FoldSeek["alnlen"].max())-10), step=1)
	
	# Select number of mismatches
	with st.expander('Number of mismatches'):
		reset_mismatch = st.checkbox("Reset mismatch", value=True)
		if reset_mismatch:
			min_mismatch, max_mismatch = st.slider("mismatch", int(df_FoldSeek["mismatch"].min()), int(df_FoldSeek["mismatch"].max()), (int(df_FoldSeek["mismatch"].min()), int(df_FoldSeek["mismatch"].max())), step=1)
		else:
			min_mismatch, max_mismatch = st.slider("mismatch", int(df_FoldSeek["mismatch"].min()), int(df_FoldSeek["mismatch"].max()), (int(df_FoldSeek["mismatch"].min())+10, int(df_FoldSeek["mismatch"].max())-10), step=1)
	
	# Number of gap open events
	with st.expander('Number of gap open events'):
		reset_gapopen = st.checkbox("Reset gapopen", value=True)
		if reset_gapopen:
			min_gapopen, max_gapopen = st.slider("gapopen", int(df_FoldSeek["gapopen"].min()), int(df_FoldSeek["gapopen"].max()), (int(df_FoldSeek["gapopen"].min()), int(df_FoldSeek["gapopen"].max())), step=1)
		else:
			min_gapopen, max_gapopen = st.slider("gapopen", int(df_FoldSeek["gapopen"].min()), int(df_FoldSeek["gapopen"].max()), (int(df_FoldSeek["gapopen"].min())+1, int(df_FoldSeek["gapopen"].max())-1), step=1)
	
	# Select alignment start position in query sequence
	with st.expander('Alignment start position in query sequence'):
		reset_qstart = st.checkbox("Reset qstart", value=True)
		if reset_qstart:
			min_qstart, max_qstart = st.slider("qstart", int(df_FoldSeek["qstart"].min()), int(df_FoldSeek["qstart"].max()), (int(df_FoldSeek["qstart"].min()), int(df_FoldSeek["qstart"].max())), step=1)
		else:
			min_qstart, max_qstart = st.slider("qstart", int(df_FoldSeek["qstart"].min()), int(df_FoldSeek["qstart"].max()), (int(df_FoldSeek["qstart"].min())+10, int(df_FoldSeek["qstart"].max())-10), step=1)

	# Select alignment end position in query sequence
	with st.expander('Alignment end position in query sequence'):
		reset_qend = st.checkbox("Reset qend", value=True)
		if reset_qend:
			min_qend, max_qend = st.slider("qend", int(df_FoldSeek["qend"].min()), int(df_FoldSeek["qend"].max()), (int(df_FoldSeek["qend"].min()), int(df_FoldSeek["qend"].max())), step=1)
		else:
			min_qend, max_qend = st.slider("qend", int(df_FoldSeek["qend"].min()), int(df_FoldSeek["qend"].max()), (int(df_FoldSeek["qend"].min())+10, int(df_FoldSeek["qend"].max())-10), step=1)

	# Select alignment start position in target sequence
	with st.expander('Alignment start position in target sequence'):
		reset_tstart = st.checkbox("Reset tstart", value=True)
		if reset_tstart:
			min_tstart, max_tstart = st.slider("tstart", int(df_FoldSeek["tstart"].min()), int(df_FoldSeek["tstart"].max()), (int(df_FoldSeek["tstart"].min()), int(df_FoldSeek["tstart"].max())), step=1)
		else:
			min_tstart, max_tstart = st.slider("tstart", int(df_FoldSeek["tstart"].min()), int(df_FoldSeek["tstart"].max()), (int(df_FoldSeek["tstart"].min())+10, int(df_FoldSeek["tstart"].max())-10), step=1)

	# Select alignment end position in target sequence
	with st.expander('Alignment end position in target sequence'):
		reset_tend = st.checkbox("Reset tend", value=True)
		if reset_tend:
			min_tend, max_tend = st.slider("tend", int(df_FoldSeek["tend"].min()), int(df_FoldSeek["tend"].max()), (int(df_FoldSeek["tend"].min()), int(df_FoldSeek["tend"].max())), step=1)
		else:
			min_tend, max_tend = st.slider("tend", int(df_FoldSeek["tend"].min()), int(df_FoldSeek["tend"].max()), (int(df_FoldSeek["tend"].min())+10, int(df_FoldSeek["tend"].max())-10), step=1)

	# Select E value
	with st.expander('E value'):
		reset_evalue = st.checkbox("Reset evalue", value=True)
		if reset_evalue:
			min_evalue, max_evalue = st.select_slider("evalue", df_FoldSeek["evalue"].sort_values(), (df_FoldSeek["evalue"].sort_values().min(), df_FoldSeek["evalue"].sort_values().max()))
		else:
			min_evalue, max_evalue = st.select_slider("evalue", df_FoldSeek["evalue"].sort_values(), (list(df_FoldSeek["evalue"].sort_values())[3], list(df_FoldSeek["evalue"].sort_values())[-3]))

	# Select bit score
	with st.expander('Bit score'):
		reset_bits = st.checkbox("Reset bits", value=True)
		if reset_bits:
			min_bits, max_bits = st.select_slider("bits", df_FoldSeek["bits"].sort_values(), (df_FoldSeek["bits"].sort_values().min(), df_FoldSeek["bits"].sort_values().max()))
		else:
			min_bits, max_bits = st.select_slider("bits", df_FoldSeek["bits"].sort_values(), (list(df_FoldSeek["bits"].sort_values())[3], list(df_FoldSeek["bits"].sort_values())[-3]))
	
	# Select target length
	with st.expander('Target length'):
		reset_tlen = st.checkbox("Reset tlen", value=True)
		if reset_tlen:
			min_tlen, max_tlen = st.slider("tlen", int(df_FoldSeek["tlen"].min()), int(df_FoldSeek["tlen"].max()), (int(df_FoldSeek["tlen"].min()), int(df_FoldSeek["tlen"].max())), step=1)
		else:
			min_tlen, max_tlen = st.slider("tlen", int(df_FoldSeek["tlen"].min()), int(df_FoldSeek["tlen"].max()), (int(df_FoldSeek["tlen"].min())+10, int(df_FoldSeek["tlen"].max())-10), step=1)
	
	# Select genus
	with st.expander('Select genus'):
		all_genus = st.checkbox("Select all genus", value=True)
		if all_genus:
			genus = st.multiselect("genus:",
			df_FoldSeek["genus"].unique(),df_FoldSeek["genus"].unique())
		else:
			genus = st.multiselect("genus:",
			df_FoldSeek["genus"].unique())
	
	# Select genus count
	with st.expander('Select genus count'):
		reset_genuscount = st.checkbox("Reset genus count", value=True)
		if reset_genuscount:
			min_genuscount, max_genuscount = st.select_slider("genuscount", df_FoldSeek["genuscount"].sort_values(), (df_FoldSeek["genuscount"].sort_values().min(), df_FoldSeek["genuscount"].sort_values().max()))
		else:
			min_genuscount, max_genuscount = st.select_slider("genuscount", df_FoldSeek["genuscount"].sort_values(), (list(df_FoldSeek["genuscount"].sort_values())[3], list(df_FoldSeek["genuscount"].sort_values())[-3]))

	# Select species
	with st.expander('Select species'):
		all_species = st.checkbox("Select all species", value=True)
		if all_species:
			species = st.multiselect("species:",
			df_FoldSeek["species"].unique(),df_FoldSeek["species"].unique())
		else:
			species = st.multiselect("species:",
			df_FoldSeek["species"].unique())

# -------------------------------------------------------- FILTER DATAFRAME ---------------------------------------------------------

df_FoldSeek_selection = df_FoldSeek.query(
	"@min_fident <= fident <= @max_fident"
	"& @min_alnlen <= alnlen <= @max_alnlen"
	"& @min_mismatch <= mismatch <= @max_mismatch"
	"& @min_gapopen <= gapopen <= @max_gapopen"
	"& @min_qstart <= qstart <= @max_qstart"
	"& @min_qend <= qend <= @max_qend"
	"& @min_tstart <= tstart <= @max_tstart"
	"& @min_tend <= tend <= @max_tend"
	"& @min_evalue <= evalue <= @max_evalue"
	"& @min_bits <= bits <= @max_bits"
	"& @min_tlen <= tlen <= @max_tlen"
	"& genus == @genus"
	"& @min_genuscount <= genuscount <= @max_genuscount"
	"& species == @species"
	)

st.dataframe(df_FoldSeek_selection)
	
# st.write(type(df['ID_No'][1]))
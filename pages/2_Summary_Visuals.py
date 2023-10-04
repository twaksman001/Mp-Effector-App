import streamlit as st  # pip install streamlit
import pandas as pd  # pip install pandas openpyxl
import plotly.express as px  # pip install plotly-express

st.set_option('deprecation.showPyplotGlobalUse', False)

# ---- LAYOUT ----

st.set_page_config(layout="wide")

st.write("This page creates useful graphs or figures from the summary spreadsheet data. Interactive Plotly scatter charts, "
		 "structure prediction per-residue confidence graphs, or all-against-all charts of multiple variables from summary spreadsheet")

tab1, tab2 = st.tabs(['Plotly', 'Other'])

# ---- READ EXCEL ----

df = pd.read_csv('Data Files/MpEffectorsBioinfo_Structural_KeyInfo.csv')

# st.dataframe(df)

# ---- VARIABLES ----

column_list = list(df.columns)
column_list.append(None)
variables_scatter_hoverdata_default = ['PDBeFold_Q', 'AF2_OF_pLDDT_Pearson', 'AF2_OF_DALI_Z']
variables_PairPlot_default = ['ID_No', 'Length', 'MSA_Depth_Mean_log', 'pLDDT_mean', 'AF2_OF_pLDDT_Pearson']

desired_file_path_AF2 = 'Data Files/AF2/*'
desired_file_path_OF = 'Data Files/OmegaFold/*'

# ---- SIDEBAR ----

with st.sidebar:
	
	st.header('Control Visuals')
	
	# button = st.button('Button')
		
	with st.expander('Scatter 1'):
		x = st.selectbox('x', column_list, 0)
		y = st.selectbox('y', column_list, 7)
		color = st.selectbox('color', column_list, 10)
		symbol = st.selectbox('symbol', column_list, 9)
		size = st.selectbox('size', column_list, 14)
		hover_name = st.selectbox('hover name', column_list, 0)
		hover_info = st.multiselect('hover extra info', df.columns, variables_scatter_hoverdata_default)

		scatter1_display = st.checkbox("Display 1")
	
	with st.expander('Scatter 2'):
		x_2 = st.selectbox('x 2', column_list, 0)
		y_2 = st.selectbox('y 2', column_list, 7)
		color_2 = st.selectbox('color 2', column_list, 10)
		symbol_2 = st.selectbox('symbol 2', column_list, 9)
		size_2 = st.selectbox('size 2', column_list, 14)
		hover_name_2 = st.selectbox('hover name 2', column_list, 0)
		hover_info_2 = st.multiselect('hover extra info 2', df.columns, variables_scatter_hoverdata_default)
		
		scatter2_display = st.checkbox("Display 2")

# ---- VISUALS DEFINITIONS ----

scatter1 = px.scatter(df, x, y, color, symbol, size,
					 hover_name, hover_data=hover_info, trendline='ols', trendline_scope='overall')
scatter1.update_layout(coloraxis_colorbar=dict(yanchor="top", y=0.79, xanchor="left", x=1.01),
					  legend=dict(yanchor="top", y=1, xanchor="left", x=1.02)
					  )
					  
scatter2 = px.scatter(df, x_2, y_2, color_2, symbol_2, size_2,
					  hover_name_2, hover_data=hover_info_2, trendline='ols', trendline_scope='overall')
scatter2.update_layout(coloraxis_colorbar=dict(yanchor="top", y=0.79, xanchor="left", x=1.01),
					  legend=dict(yanchor="top", y=1, xanchor="left", x=1.02)
					  )

# PAE_image = 'C:/Users/TW43969/OneDrive - University of Dundee/Current/Research/Tom W Research Record/Bioinformatics/Structure/Structure Prediction/CoLabFold/Effectors_PAE.png'			  
# pLDDT_image = 'C:/Users/TW43969/OneDrive - University of Dundee/Current/Research/Tom W Research Record/Bioinformatics/Structure/Structure Prediction/CoLabFold/Effectors_pLDDT.png'			  
# Image1 = 'https://cdn-images-1.medium.com/max/1024/1*u9U3YjxT9c9A1FIaDMonHw.png'
# Image2 = 'C:/Users/TW43969/Downloads/Effectors_PAE.png'

# ---- FUNCTIONS ----

def get_file_of_interest_AF2(protein):
	
	import glob
	import re
	
	for file_path in glob.glob(desired_file_path_AF2):

		if re.search(protein + '_.{5}_.*rank_1*', file_path):
			
			return file_path.replace('\\', '/')

def get_PAE_dataframe_from_scores_file(protein):
	
	import re
	import pandas as pd
	
	scores_file = open(get_file_of_interest_AF2(protein)).read()
	
	if scores_file:
	
		scores_file_PAE = scores_file[scores_file.find('"pae"')+8:scores_file.find('plddt')-4]
	
	PAE_df_dict = {}

	counter = 1

	for match in re.finditer('\[[\.\d\s,]*]', scores_file_PAE):

		match_values = match.group()[1:len(match.group())-1]

		PAE_df_dict[counter] = list(map(float,match_values.split(',')))

		counter += 1

	df = pd.DataFrame.from_dict(PAE_df_dict, orient='index', dtype=None, columns=range(1, counter))
	
	return df

def figure_PAE_proteins(proteins, tick_no=5):
	
	import numpy as np
	import seaborn as sea
	import matplotlib.pyplot as plt
	# %matplotlib inline 
	
	if len(proteins) > 30:
		spacing=len(proteins)/160
	else:
		spacing=len(proteins)/60
		
	fig = plt.figure(figsize=(25,25))
	plt.subplots_adjust(hspace=spacing, wspace=spacing)

	n = 1

	for protein in proteins:
		
		if protein == 0:
			protein = 'MpC002'
		elif protein == 92:
			protein = '92a'
		elif protein == 100:
			protein = 'MIF1'
		else:
			protein = 'Mp' + str(protein)
		
		df = get_PAE_dataframe_from_scores_file(protein)
		
		plt.subplot(int(np.ceil(len(proteins)**0.5)), int(np.ceil(len(proteins)**0.5)), n)
		plt.title(protein, fontdict={'weight':'bold'})
		sea.heatmap(df, cmap='rainbow', square=True, linewidths=0,
					xticklabels=len(df)//tick_no, yticklabels=len(df)//tick_no, cbar=False) #, vmin=0, vmax=50)

		n += 1		   

def get_pLDDT_from_scores_file(protein):

	import re
	
	scores_file = open(get_file_of_interest_AF2(protein)).read()
	
	scores_file_pLDDT = scores_file[scores_file.find('"plddt"'):]

	for match in re.finditer('\[[\.\d\s,]*]', scores_file_pLDDT):

		match_values = match.group()[1:len(match.group())-1]

		pLDDT = list(map(float,match_values.split(',')))
	
	return pLDDT

def get_file_of_interest_OF(protein):
	
	import glob
	import re
	
	for file_path in glob.glob(desired_file_path_OF):

		if re.search(protein + '_', file_path):
						
			return file_path.replace('\\', '/')

def pdb_file_bfactor_residue(protein):
	
	import re
	
	bfactors_residue = []
	
	file_lines = open(get_file_of_interest_OF(protein)).read().split('\n')
	
	for i in range(len(file_lines)):
		
		if re.split('\s{1,}', file_lines[i])[5] != re.split('\s{1,}', file_lines[i-1])[5]:
						
			for match in re.finditer('\s.\d\.\d\d\s{5,}[NCO]', file_lines[i]):
			
				bfactors_residue.append(float(match.group()[1:6]))

	return bfactors_residue

def figure_AF2andOF_confidence_proteins(proteins, tick_no=5):

	import numpy as np
	import matplotlib.pyplot as plt
	# %matplotlib inline
	
	if len(proteins) > 30:
		spacing=len(proteins)/160
	else:
		spacing=len(proteins)/100
	
	fig = plt.figure(figsize=(25,25))
	plt.subplots_adjust(hspace=spacing, wspace=spacing)

	n = 1

	for protein in proteins:

		if protein == 0:
			protein = 'MpC002'
		elif protein == 92:
			protein = '92a'
		elif protein == 100:
			protein = 'MIF1'
		else:
			protein = 'Mp' + str(protein)
		
		bfactors = pdb_file_bfactor_residue(protein)
		pLDDT = get_pLDDT_from_scores_file(protein)
		
		plt.subplot(int(np.ceil(len(proteins)**0.5)), int(np.ceil(len(proteins)**0.5)), n)
		plt.title(protein, fontdict={'weight':'bold'})
		if OF_check:
			plt.plot(bfactors, color='blue', label = 'OmegaFold')
		if AF2_check:
			plt.plot(pLDDT, color='orange', label = 'AlphaFold 2')
		
		plt.xticks(ticks=range(0, len(bfactors)+1, len(bfactors)//tick_no))
		plt.yticks(ticks=range(0, 101, 100//tick_no))
		
		# if n==1:# and pLDDT_legend_check:
		plt.legend(loc='lower left')

		n += 1

def figure_corr(proteins,tick_no=5):

	import scipy
	import numpy as np
	import matplotlib.pyplot as plt
	# %matplotlib inline
	
	if len(proteins) > 30:
		spacing=len(proteins)/160
	else:
		spacing=len(proteins)/70
	
	fig = plt.figure(figsize=(25,25))
	plt.subplots_adjust(hspace=spacing, wspace=spacing)

	n = 1
	
	for protein in proteins:
	
		if protein == 0:
			protein = 'MpC002'
		elif protein == 92:
			protein = '92a'
		elif protein == 100:
			protein = 'MIF1'
		else:
			protein = 'Mp' + str(protein)
		
		pLDDT = get_pLDDT_from_scores_file(protein)
		bfactors = pdb_file_bfactor_residue(protein)
		corr1, p = scipy.stats.pearsonr(pLDDT, bfactors)
		corr2, p = scipy.stats.spearmanr(pLDDT, bfactors)

		plt.subplot(int(np.ceil(len(proteins)**0.5)), int(np.ceil(len(proteins)**0.5)), n)
		plt.title(str(protein) + ' ' + str(round(corr1, 2)) + ' ' + str(round(corr2, 2)), fontdict={'weight':'bold'})
		plt.scatter(pLDDT, bfactors, s=1/(len(proteins)/70))
		plt.xlabel('AlphaFold 2')
		plt.ylabel('OmegaFold')
		plt.xticks(ticks=range(0, 101, 100//tick_no))
		plt.yticks(ticks=range(0, 101, 100//tick_no))

		n += 1
	
def subplot_positions_list(number_of_variables):	
	
	positions_all = [1]

	l = number_of_variables
	row = 2
	position = l + 1

	while position <= l**2:

		for i in range(row):

			positions_all.append(position + i)

		position = l*row + 1

		row += 1
	
	return positions_all

def PairPlot(variables, hue, size, style):
	
	import seaborn as sns
	import matplotlib.pyplot as plt
	# %matplotlib inline
	
	sns.set_theme(font_scale = 1.5) #palette = 'bright')
	sns.set_style('white', {'xtick.bottom': True, 'ytick.left': True})

	fig = plt.figure(figsize=(25,25))
	plt.subplots_adjust(hspace=0.05, wspace=0.05)
	
	l = len(variables)
	subplot_positions_all = subplot_positions_list(l)
	n = 1
	
	for variable1 in variables:

		for variable2 in variables:

			if n in subplot_positions_all:

				plt.subplot(l, l, n)

				if variable1 == variable2:   
					plot = sns.kdeplot(
									  data = df, x = variable1, legend = False,
									  hue = hue, palette = 'viridis',
									  )

					plot.set(ylabel = None, yticklabels = [])

					if n != l**2:
						plot.set(
								xlabel = None, xticklabels = [] 
								)
				
				else:
					
					if n == 1+l:
						plot = sns.scatterplot(
											  data = df,
											  x = variable2,
											  y = variable1,
											  legend = True,
											  hue = hue, palette = 'viridis',
											  size = size, style = style,
											  s = 100
											  )
						sns.move_legend(plot, bbox_to_anchor=(l, 2), loc='upper right', borderaxespad=0)
					
					else: 
						plot = sns.scatterplot(
											  data = df,
											  x = variable2,
											  y = variable1,
											  legend = False,
											  hue = hue, palette = 'viridis',
											  size = size, style = style,
											  s = 100
											  )

				if n not in list(range(l**2-l+1, l**2+1)):
					plot.set(
					xlabel = None, xticklabels = [] 
					)

				if n not in list(range(1, l**2, l)):
					plot.set(
					ylabel = None, yticklabels = [] 
					)

			n += 1

def update_plotly_hover_dict(variables):
	
	for variable in variables:
		
		plotly_hover_dict[variable] = True
	
	return plotly_hover_dict

# ---- DISPLAY ----

with tab1:
	
	st.write('Use sidebar to make graphs')
	
	col1, col2 = st.columns(2)
								  
	if scatter1_display:
		col1.plotly_chart(scatter1)
				
	if scatter2_display:
		col2.plotly_chart(scatter2)

with tab2:
	
	tab2_1, tab2_2, tab2_3 = st.tabs(['PAE', 'Per-Residue Confidence', 'PairPlot'])
	
	with tab2_1:	
		
		all_effectors_PAE = st.checkbox("Select all effectors PAE", value=True)
		
		with st.form('Select Effectors PAE', clear_on_submit=False):

			if all_effectors_PAE:
				effectors_PAE = st.multiselect("Effectors PAE:",
				df["ID_No"].unique(),df["ID_No"].unique())
			else:
				effectors_PAE = st.multiselect("Effectors PAE:",
				df["ID_No"].unique())
			
			make_figure_PAE = st.form_submit_button('Make Figure PAE')
					
		if make_figure_PAE:
			
			PAE_Figure = figure_PAE_proteins(effectors_PAE)
			# st.set_option('deprecation.showPyplotGlobalUse', False)
			st.pyplot(PAE_Figure)
			
	with tab2_2:	
		
		all_effectors_pLDDT = st.checkbox("Select all effectors Per-Residue Confidence", value=True)
		
		with st.form('Select Effectors Per-Residue Confidence', clear_on_submit=False):
			
			OF_check = st.checkbox("Omegafold", value=True)
			AF2_check = st.checkbox("AlphaFold 2", value=True)
			
			if all_effectors_pLDDT:
				effectors_pLDDT = st.multiselect("Effectors Per-Residue Confidence:",
				df["ID_No"].unique(),df["ID_No"].unique())
			else:
				effectors_pLDDT = st.multiselect("Effectors Per-Residue Confidence:",
				df["ID_No"].unique())
			
			make_figure_pLDDT_plot = st.form_submit_button('Plot Per-Residue Confidence')
			make_figure_pLDDT_corr = st.form_submit_button('Compare Per-Residue Confidence')
					
		if make_figure_pLDDT_plot:
			
			figure_pLDDT_plot = figure_AF2andOF_confidence_proteins(effectors_pLDDT)
			st.pyplot(figure_pLDDT_plot)
			
		if make_figure_pLDDT_corr:
			
			figure_pLDDT_corr = figure_corr(effectors_pLDDT)
			st.pyplot(figure_pLDDT_corr)
	
	with tab2_3:	
		
		# all_effectors_pLDDT = st.checkbox("Select all effectors Per-Residue Confidence", value=True)
		
		with st.form('Select Variables for PairPlot', clear_on_submit=False):
			
			variables_PairPlot = st.multiselect('Variables for PairPlot', df.columns, variables_PairPlot_default)
			variables_PairPlot_hue = st.selectbox('Select Color Variable', column_list, 10)
			variables_PairPlot_size = st.selectbox('Select Size Variable', column_list, 14)
			variables_PairPlot_style = st.selectbox('Select Marker Variable', column_list, 9)
			
			make_PairPlot = st.form_submit_button('Make PairPlot')
					
		if make_PairPlot:
			
			PairPlot = PairPlot(variables_PairPlot, variables_PairPlot_hue, variables_PairPlot_size, variables_PairPlot_style)
			st.pyplot(PairPlot)
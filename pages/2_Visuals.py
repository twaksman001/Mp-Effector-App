####################modules####################

import streamlit as st
import pandas as pd
import plotly.express as px

####################settings####################

st.set_option(key='deprecation.showPyplotGlobalUse', value=False)

####################layout####################

st.set_page_config(layout='wide')

st.write('This page creates useful graphs or figures from the summary spreadsheet data. Interactive plotly scatter charts, '
		 'structure prediction per-residue confidence graphs, or all-against-all charts of multiple variables from summary spreadsheet.')

tab1, tab2 = st.tabs(tabs=['plotly', 'other'])

####################load dataframe####################

df = pd.read_csv(filepath_or_buffer='Data Files/MpEffectors_KeyInfo.txt', sep='\t')

####################dataframe functions####################

column_list_graphs = list((set(df.columns).union({None})).difference({'protein', 'sequence'}))
variables_pairplot_default = ['length', 'pTM', 'Pearson']

desired_file_path_AF2 = 'Data Files/AF2/*'
desired_file_path_OF = 'Data Files/OmegaFold/*'

def prepare_dataframe(df):

	import numpy as np
	
	columns_string = {'protein', 'sequence', 'beta', 'modular', 'function'}
	columns_float = {'pTM', 'SASA/length', 'flDPNN_AF', 'Pearson', 'DALI_Z', 'RMSD'}
	columns_int = {'length'}
	
	for column in df.columns:
		if column in columns_string:
			df[column] = df[column].astype(str)
		elif column in columns_float:
			df[column] = df[column].astype(float)
		elif column in columns_int:
			df[column] = df[column].astype(int)
	
	# df['RMSD'] = df['RMSD'].fillna('-')
	
	return df

df = prepare_dataframe(df)

####################sidebar####################

with st.sidebar:
	
	st.header(body='control visuals')
			
	with st.expander(label='scatter 1', expanded=False):
		x = st.selectbox(label='x', options=column_list_graphs, index=column_list_graphs.index('pTM'))
		y = st.selectbox(label='y', options=column_list_graphs, index=column_list_graphs.index('SASA/length'))
		color = st.selectbox(label='color', options=column_list_graphs, index=column_list_graphs.index('RMSD'))
		symbol = st.selectbox(label='symbol', options=column_list_graphs, index=column_list_graphs.index('beta'))
		size = st.selectbox(label='size', options=column_list_graphs, index=column_list_graphs.index(None))
		hover_info = st.multiselect(label='hover extra info', options=column_list_graphs, default=column_list_graphs)

		scatter1_display = st.checkbox(label='display 1', key='display 1', value=False)
	
	with st.expander(label='scatter 2', expanded=False):
		x_2 = st.selectbox(label='x 2', options=column_list_graphs, index=column_list_graphs.index('pTM'))
		y_2 = st.selectbox(label='y 2', options=column_list_graphs, index=column_list_graphs.index('SASA/length'))
		color_2 = st.selectbox(label='color 2', options=column_list_graphs, index=column_list_graphs.index('RMSD'))
		symbol_2 = st.selectbox(label='symbol 2', options=column_list_graphs, index=column_list_graphs.index('beta'))
		size_2 = st.selectbox(label='size 2', options=column_list_graphs, index=column_list_graphs.index(None))
		hover_info_2 = st.multiselect(label='hover extra info 2', options=column_list_graphs, default=column_list_graphs)		
		
		scatter2_display = st.checkbox(label='display 2', key='display 2', value=False)

####################visuals definitions####################

scatter1 = px.scatter(data_frame=df, x=x, y=y, color=color, symbol=symbol, size=size, 
					 hover_name='protein', hover_data=hover_info, trendline='ols', trendline_scope='overall')
scatter1.update_layout(coloraxis_colorbar=dict(yanchor='top', y=0.79, xanchor='left', x=1.01),
					  legend=dict(yanchor='top', y=1, xanchor='left', x=1.02))
					  
scatter2 = px.scatter(data_frame=df, x=x_2, y=y_2, color=color_2, symbol=symbol_2, size=size_2, 
					 hover_name='protein', hover_data=hover_info_2, trendline='ols', trendline_scope='overall')
scatter2.update_layout(coloraxis_colorbar=dict(yanchor='top', y=0.79, xanchor='left', x=1.01),
					  legend=dict(yanchor='top', y=1, xanchor='left', x=1.02))

####################functions####################

def whoami():
	import inspect
	return inspect.stack()[1][3]

def get_file(protein, path, string):
	
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

def get_AF2_scores(protein, path, string, type_='pTM'): # type = pTM, pLDDT, PAE
	
	file = get_file(protein=protein, path=path, string=string)
	file = open(file).read()
	
	if type_ == 'pTM':
		pTM = file[file.rfind(' ')+1:-1]
		return pTM
	if type_ == 'pLDDT':
		pLDDT = file[file.rfind('[')+1:file.rfind(']')]
		pLDDT = list(map(float, pLDDT.split(',')))
		return pLDDT
	if type_ == 'PAE':
		import pandas as pd
		PAE = file[file.find('"pae":')+9:file.find('plddt')-5]
		PAE = PAE.split('], [')
		PAE = [i.split(', ') for i in PAE]
		df = pd.DataFrame(data=PAE).astype(float)
		return df

def get_disorder_list(protein):
	
	import pandas as pd
	
	df = pd.read_csv(filepath_or_buffer='Data Files/MpEffectors_fIDPnn_2.txt', sep='\t', index_col=0, header=0)
	list_ = df.loc[protein, 'disorder_score']
	list_ = list_.split(',')
	list_ = [float(i)*100 for i in list_]
	
	return list_

def df_PDBfile(protein, path, string):
	
	import pandas as pd
	
	positions = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26), (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78), (78, 80)]
	columns = ['entity', 'serial', 'atom_aa', 'altloc', 'aa', 'chain', 'res_index', 'icode', 'x', 'y', 'z', 'occ', 'B', 'element', 'charge']
	
	file = get_file(protein=protein, path=path, string=string)
	df = pd.read_fwf(filepath_or_buffer=file, names=columns, colspecs=positions)
	df = df[df['entity']=='ATOM']
	
	return df

def PDB_bfactor_list(protein, path='Data Files/OmegaFold/*', string='_'):
	
	df = df_PDBfile(protein=protein, path=path, string=string)
	b_list = list(df[df['atom_aa'] == 'CA']['B'].astype(float))
	
	return b_list

def figure_perresidue_plot_proteins(proteins, algorithms, tick_no=5):

	import numpy as np
	import matplotlib.pyplot as plt
	
	st.markdown('**legend**: :orange[AlphaFold]  ;  :blue[OmegaFold] ; :green[flDPNN]')

	spacing = (len(proteins)/160 if len(proteins)>=50 else len(proteins)/120 if 30<=len(proteins)<50
				else len(proteins)/80 if 15<=len(proteins)<30 else len(proteins)/40)
	fontsize = (500/len(proteins) if len(proteins)>=50 else 400/len(proteins) if 30<=len(proteins)<50
				else 225/len(proteins) if 15<=len(proteins)<30 else 15)
	fig = plt.figure(figsize=(25,25))
	plt.subplots_adjust(hspace=spacing, wspace=spacing)
	n = 1

	for protein in proteins:

		plt.subplot(int(np.ceil(len(proteins)**0.5)), int(np.ceil(len(proteins)**0.5)), n)
		plt.title(protein, fontdict={'size':fontsize*1.5, 'weight':'bold'})

		if 'OmegaFold' in algorithms:
			OF_conf = PDB_bfactor_list(protein=protein, path='Data Files/OmegaFold/*', string='_')
			plt.plot(OF_conf, color='blue', label='OmegaFold')
		if 'AlphaFold' in algorithms:
			pLDDT = get_AF2_scores(protein=protein, path='Data Files/AF2/*', string='_', type_='pLDDT')
			plt.plot(pLDDT, color='orange', label='AlphaFold')
		if 'flDPNN' in algorithms:
			flDPNN = get_disorder_list(protein=protein)
			plt.plot(flDPNN, color='green', label='flDPNN')

		plt.margins(x=0)
		plt.xlim(0, len(pLDDT))
		plt.xticks(ticks=[int(i) for i in np.linspace(start=0, stop=len(pLDDT), num=tick_no, endpoint=True)], fontsize=fontsize)
		plt.yticks(ticks=range(0, 101, 100//tick_no), fontsize=fontsize)

		n += 1

def figure_corr(proteins, algorithms, tick_no=5):

	import scipy
	import numpy as np
	import matplotlib.pyplot as plt

	spacing = (len(proteins)/120 if len(proteins)>=50 else len(proteins)/80 if 30<=len(proteins)<50
				else len(proteins)/60 if 15<=len(proteins)<30 else len(proteins)/30)
	fontsize = (500/len(proteins) if len(proteins)>=50 else 400/len(proteins) if 30<=len(proteins)<50
				else 225/len(proteins) if 15<=len(proteins)<30 else 15)
	fig = plt.figure(figsize=(25,25))
	plt.subplots_adjust(hspace=spacing, wspace=spacing)

	n = 1
	
	if len(algorithms) > 2:
		st.write('use maximum of two algorithms for scatter plot')
		algorithms = algorithms[:2]

	for protein in proteins:
		
		plot_list = []

		if 'OmegaFold' in algorithms:
			OF_conf = PDB_bfactor_list(protein=protein, path='Data Files/OmegaFold/*', string='_')
			plot_list.append(OF_conf)
		if 'AlphaFold' in algorithms:
			pLDDT = get_AF2_scores(protein=protein, path='Data Files/AF2/*', string='_', type_='pLDDT')
			plot_list.append(pLDDT)
		if 'flDPNN' in algorithms:
			flDPNN = get_disorder_list(protein=protein)
			plot_list.append(flDPNN)
		
		corr1, p = scipy.stats.pearsonr(plot_list[0], plot_list[1])
		corr2, p = scipy.stats.spearmanr(plot_list[0], plot_list[1])

		plt.subplot(int(np.ceil(len(proteins)**0.5)), int(np.ceil(len(proteins)**0.5)), n)
		plt.title(protein, fontdict={'size':fontsize*1.5, 'weight':'bold'}) #+ str(round(corr1, 2)) + ' ' + str(round(corr2, 2)))
		plt.scatter(plot_list[0], plot_list[1], s=70/(len(proteins)))
		plt.xlabel(xlabel=algorithms[0], fontdict={'size':fontsize})
		plt.ylabel(ylabel=algorithms[1], fontdict={'size':fontsize})
		plt.xticks(ticks=range(0, 101, 100//tick_no), fontsize=fontsize)
		plt.yticks(ticks=range(0, 101, 100//tick_no), fontsize=fontsize)

		n += 1
	
def figure_PAE_proteins(proteins, tick_no=5):
	
	import numpy as np
	import seaborn as sea
	import matplotlib.pyplot as plt
	
	spacing = (len(proteins)/160 if len(proteins)>=50 else len(proteins)/120 if 30<=len(proteins)<50
				else len(proteins)/80 if 15<=len(proteins)<30 else len(proteins)/40)
	fontsize = (500/len(proteins) if len(proteins)>=50 else 400/len(proteins) if 30<=len(proteins)<50
				else 225/len(proteins) if 15<=len(proteins)<30 else 15)
	fig = plt.figure(figsize=(25,25))
	plt.subplots_adjust(hspace=spacing, wspace=spacing)

	n = 1

	for protein in proteins:
		
		df = get_AF2_scores(protein=protein, path='Data Files/AF2/*', string='_', type_='PAE')
		ticks = [int(i) for i in np.linspace(start=0, stop=len(df), num=tick_no, endpoint=True)]
		plt.subplot(int(np.ceil(len(proteins)**0.5)), int(np.ceil(len(proteins)**0.5)), n)
		plt.title(label=protein, fontdict={'size':fontsize*1.5, 'weight':'bold'})
		sea.heatmap(data=df, cmap='rainbow', square=True, cbar=False, linewidths=0)

		plt.xticks(ticks=ticks, labels=ticks, fontsize=fontsize, rotation='horizontal')
		plt.yticks(ticks=ticks, labels=ticks, fontsize=fontsize)
		plt.axhline(y=0, linewidth=2, color='black')
		plt.axhline(y=len(df), linewidth=2, color='black')
		plt.axvline(x=0, linewidth=2, color='black')
		plt.axvline(x=len(df), linewidth=2, color='black')
		
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

def pairplot(df, variables, hue, size, style):
	
	import seaborn as sns
	import matplotlib.pyplot as plt
	
	sns.set_theme(font_scale = 1.5) #palette = 'bright')
	sns.set_style('white', {'xtick.bottom': True, 'ytick.left': True})

	fig = plt.figure(figsize=(25,25))
	plt.subplots_adjust(hspace=0.05, wspace=0.05)
	
	l = len(variables)
	subplot_positions_all = subplot_positions_list(l)
	n = 1
	NA_mask = df[list({hue, size, style}.difference({None}))].isna().any(axis=1)

	for variable1 in variables:
		for variable2 in variables:
			if n in subplot_positions_all:

				plt.subplot(l, l, n)

				if variable1 == variable2:   
					plot = sns.kdeplot(data=df, x=variable1, legend=False,
									   hue=hue, palette='viridis')
					if n != 1:
						plot.set(ylabel=None, yticklabels=[])
					else:
						plot.set(ylabel=variable1, yticklabels=[])
					if n != l**2:
						plot.set(xlabel=None, xticklabels=[])
				else:
					if n == 1+l:
						plot = sns.scatterplot(data=df[~NA_mask], x=variable2, y=variable1, legend=True, 
											  hue=hue, palette='viridis', size=size, style=style, s=100)
						plot2 = plt.scatter(df[NA_mask][variable2], df[NA_mask][variable1], 
											color='white', edgecolors='black', label='N/A', marker='s', s=20)
						sns.move_legend(plot, bbox_to_anchor=(l, 2), loc='upper right', borderaxespad=0)
					else:
						plot = sns.scatterplot(data=df[~NA_mask], x=variable2, y=variable1, legend=False,
											  hue=hue, palette='viridis', size=size, style=style, s=100)
						plot2 = plt.scatter(df[NA_mask][variable2], df[NA_mask][variable1], 
											color='white', edgecolors='black', label='N/A', marker='s', s=20)
				if n not in list(range(l**2-l+1, l**2+1)):
					plot.set(xlabel = None, xticklabels = [])
				if n not in list(range(1, l**2, l)):
					plot.set(ylabel = None, yticklabels = [])

			n += 1
	
	# plt.figlegend()#[plot2], ['N/A'])
	
def update_plotly_hover_dict(variables):
	
	for variable in variables:
		plotly_hover_dict[variable] = True
	
	return plotly_hover_dict

####################display####################

with tab1:
	st.write('use sidebar to make graphs')
	col1, col2 = st.columns(spec=[0.5,0.5])
								  
	if scatter1_display:
		col1.plotly_chart(figure_or_data=scatter1, use_container_width=True)
				
	if scatter2_display:
		col2.plotly_chart(figure_or_data=scatter2, use_container_width=True)

with tab2:
	tab2_1, tab2_2, tab2_3 = st.tabs(tabs=['PAE', 'per-residue scores', 'pairplot'])
	
	with tab2_1:	
		all_effectors_PAE = st.checkbox(label='select all effectors PAE', key='select all effectors PAE', value=True)
		
		with st.form(key='select effectors PAE', clear_on_submit=False):
			if all_effectors_PAE:
				effectors_PAE = st.multiselect(label='effectors PAE:', key='effectors PAE:',
											   options=df['protein'].unique(), default=df['protein'].unique())
			else:
				effectors_PAE = st.multiselect(label='effectors PAE:', key='effectors PAE:',
											   options=df['protein'].unique())
			
			make_figure_PAE = st.form_submit_button(label='make Figure PAE', )
	
		if make_figure_PAE:
			PAE_Figure = figure_PAE_proteins(proteins=effectors_PAE, tick_no=5)
			st.pyplot(fig=PAE_Figure, clear_figure=False, use_container_width=True)
			
	with tab2_2:	
		st.write('This section is for comparison of per-residue data from different algorithms. '
		   		 'The data available are AlphaFold pLDDT, OmegaFold confidence (the equivalent of pLDDT), '
				 'and probability of residue occurring in disordered region (flDPNN).')
		all_effectors_pLDDT = st.checkbox(label='select all effectors', 
											key='select all effectors per-residue scores', value=True)
		
		with st.form(key='select effectors per-residue scores', clear_on_submit=False):
			algorithms_pLDDT = st.multiselect(label='select algorithms:', key='algorithms per-residue scores', 
											  options=['AlphaFold','OmegaFold','flDPNN'], default=['AlphaFold','OmegaFold','flDPNN'])

			if all_effectors_pLDDT:
				effectors_pLDDT = st.multiselect(label='select effectors:', key='effectors per-residue scores',
												 options=df['protein'].unique(), default=df['protein'].unique())
			else:
				effectors_pLDDT = st.multiselect(label='select effectors:', key='effectors per-residue scores',
												 options=df['protein'].unique())
			
			make_figure_pLDDT_plot = st.form_submit_button(label='plot per-residue scores')
			make_figure_pLDDT_corr = st.form_submit_button(label='compare per-residue scores')
					
		if make_figure_pLDDT_plot:
			figure_pLDDT_plot = figure_perresidue_plot_proteins(proteins=effectors_pLDDT, algorithms=algorithms_pLDDT, tick_no=5)
			st.pyplot(fig=figure_pLDDT_plot, clear_figure=False, use_container_width=True)
		if make_figure_pLDDT_corr:
			figure_pLDDT_corr = figure_corr(proteins=effectors_pLDDT, algorithms=algorithms_pLDDT, tick_no=5)
			st.pyplot(fig=figure_pLDDT_corr, clear_figure=False, use_container_width=True)
	
	with tab2_3:	
		with st.form(key='form pairplot', clear_on_submit=False):
			variables_pairplot = st.multiselect(label='select variables:', key='variables pairplot', 
												options=column_list_graphs, default=variables_pairplot_default)
			variables_pairplot_hue = st.selectbox(label='select color variable', key='color pairplot', 
												options=column_list_graphs, index=column_list_graphs.index('SASA/length'))
			variables_pairplot_size = st.selectbox(label='select size variable', key='size pairplot', 
												options=column_list_graphs, index=column_list_graphs.index(None))
			variables_pairplot_style = st.selectbox(label='select marker variable', key='marker pairplot', 
												options=column_list_graphs, index=column_list_graphs.index('beta'))
			
			st.write('If any variables selected for color, size or marker contain missing values, the data points are plotted '
					 'as small hollow squares.')
			make_pairplot = st.form_submit_button(label='make pairplot')
					
		if make_pairplot:
			pairplot = pairplot(df=df, variables=variables_pairplot, hue=variables_pairplot_hue, 
					   			size=variables_pairplot_size, style=variables_pairplot_style)
			st.pyplot(fig=pairplot, clear_figure=False, use_container_width=True)
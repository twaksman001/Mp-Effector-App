# ---- PACKAGES ----

import streamlit as st  # pip install streamlit
import pandas as pd  # pip install pandas openpyxl

# ---- LAYOUT ----

st.set_page_config(layout="wide")

st.write('This page loads an interactive spreadsheet of key information for all effectors. ')
st.write('Information is from XTALPRED (columns to the left), '
		 'ColabFold AlphaFold 2 (AF2), or ColabFold OmegaFold (OF). OF was only used in comparison to AF2 protein structures.')
st.write('Please see specific column definitions below.')

with st.expander('Column definitions'):
	st.write('Longest_Disorder_Region_percent: longest contiguous region of primary sequence'
			 'predicted to be disordered, as percentage of total sequence length')
	st.write('MSA_Depth_Mean: approximate average (over all amino acids) depth of AlphaFold 2 sequence alignment')
	st.write('pLDDT_mean: approximate average (over all amino acids) AlphaFold 2 pLDDT')
	st.write('PDBeFold_Q: Q score of PDBeFold top match')
	st.write('Beta_Sheet: protein contains beta-sheet secondary structure yes/no')
	st.write('Foldedness_Rating: rating from 1-5 of how compact/folded the protein structure looks (1 - disordered, 5 - very compact)')
	st.write('Multi_Module: protein contains multiple domains/modules yes/no')
	st.write('AF2_OF_pLDDT_Pearson: Pearson correlation coefficient (r) between AlphaFold 2 '
			 'and OmegaFold per-residue confidence estimates')
	st.write('AF2_OF_DALI_Z: DALI Z score comparing AlphaFold 2 and OmegaFold structures')

# ---- READ EXCEL ----

df = pd.read_csv('Data Files/MpEffectorsBioinfo_Structural_KeyInfo.csv')

# ---- SIDEBAR ----

with st.sidebar:
	
	st.header('Filters')
	
	# Select Effectors
	with st.expander('Effectors'):
		all_Effectors = st.checkbox("Select all Effectors", value=True)
		if all_Effectors:
			ID_Number = st.multiselect("Effectors:",
			df["ID_No"].unique(),df["ID_No"].unique())
		else:
			ID_Number = st.multiselect("Effectors:",
			df["ID_No"].unique())

	with st.expander('ID Number'):
		reset_ID = st.checkbox("Reset ID Number", value=True)
		if reset_ID:
			minID, maxID = st.slider("ID Number", 0, 100, (0, 100), 1)
		else:
			minID, maxID = st.slider("ID Number", 0, 100, (25, 75), 1)

	# Select Length
	with st.expander('Length'):
		reset_length = st.checkbox("Reset Length", value=True)
		if reset_length:
			min_len, max_len = st.slider("Length", int(df["Length"].min()), int(df["Length"].max()), (int(df["Length"].min()), int(df["Length"].max())), step=1)
		else:
			min_len, max_len = st.slider("Length", int(df["Length"].min()), int(df["Length"].max()), (50, 300), step=1)
	
	# Select MW
	with st.expander('MW'):
		reset_MW = st.checkbox("Reset MW", value=True)
		if reset_MW:
			min_MW, max_MW = st.slider("MW (kDa)", int(df["MW"].min()), int(df["MW"].max()), (int(df["MW"].min()), int(df["MW"].max())), step=1)
		else:
			min_MW, max_MW = st.slider("MW (kDa)", int(df["MW"].min()), int(df["MW"].max()), (10, 60), step=1)

	# Select pI
	with st.expander('pI'):
		reset_pI = st.checkbox("Reset pI", value=True)
		if reset_pI:
			min_pI, max_pI = st.slider("pI", 0.0, 14.0, (0.0, 14.0), 0.1)
		else:
			min_pI, max_pI = st.slider("pI", 0.0, 14.0, (4.0, 10.0), 0.1)

	# Select Longest Disorder Region %
	with st.expander('Longest Disorder Region %'):
		reset_disorder = st.checkbox("Reset Longest Disorder Region %", value=True)
		if reset_disorder:
			min_disorder, max_disorder = st.slider("Longest Disorder Region %", 0, 100, (0, 100), 1)
		else:
			min_disorder, max_disorder = st.slider("Longest Disorder Region %", 0, 100, (10, 50), 1)

	# Select MSA Depth Mean
	with st.expander('MSA Depth Mean'):
		reset_MSA_Depth_Mean = st.checkbox("Reset MSA Depth Mean", value=True)
		if reset_MSA_Depth_Mean:
			min_MSA_Depth_Mean, max_MSA_Depth_Mean = st.select_slider("MSA Depth Mean", options=df["MSA_Depth_Mean"].sort_values(), value=(int(df["MSA_Depth_Mean"].min()), int(df["MSA_Depth_Mean"].max())))
		else:
			min_MSA_Depth_Mean, max_MSA_Depth_Mean = st.select_slider("MSA Depth Mean", options=df["MSA_Depth_Mean"].sort_values(), value=(30,100))

	# Select pLDDT mean
	with st.expander('pLDDT Mean'):
		reset_pLDDT_mean = st.checkbox("Reset pLDDT mean", value=True)
		if reset_pLDDT_mean:
			min_pLDDT_mean, max_pLDDT_mean = st.slider("pLDDT mean", 0, 100, (0, 100), 1)
		else:
			min_pLDDT_mean, max_pLDDT_mean = st.slider("pLDDT mean", 0, 100, (50, 100), 1)

	# Select PDBeFold Q
	with st.expander('PDBeFold Q'):
		reset_PDBeFold_Q = st.checkbox("Reset PDBeFold Q", value=True)
		if reset_PDBeFold_Q:
			min_PDBeFold_Q, max_PDBeFold_Q = st.slider("PDBeFold Q", 0.0, 1.0, (0.0, 1.0), 0.01)
		else:
			min_PDBeFold_Q, max_PDBeFold_Q = st.slider("PDBeFold Q", 0.0, 1.0, (0.1, 0.7), 0.01)

	# Select Beta Sheet
	with st.expander('Beta Sheet'):
		all_Beta_Sheet = st.checkbox("Select all Beta Sheet", value=True)	
		if all_Beta_Sheet:
			Beta_Sheet = st.multiselect("Beta Sheet",
										df["Beta_Sheet"].unique(),
										df["Beta_Sheet"].unique())
		else:
			Beta_Sheet =  st.multiselect("Beta Sheet",
										df["Beta_Sheet"].unique())

	# Select Foldedness Rating
	with st.expander('Foldedness Rating'):
		all_Foldedness_Rating = st.checkbox("Select all Foldedness Rating", value=True)	
		if all_Foldedness_Rating:
			Foldedness_Rating = st.multiselect("Foldedness Rating",
				 df["Foldedness_Rating"].unique(),df["Foldedness_Rating"].unique())
		else:
			Foldedness_Rating =  st.multiselect("Foldedness Rating",
				df["Foldedness_Rating"].unique())

	# Select Multi Module
	with st.expander('Multi Module'):
		all_Multi_Module = st.checkbox("Select all Multi Module", value=True)	 
		if all_Multi_Module:
			Multi_Module = st.multiselect("Multi Module",
				 df["Multi_Module"].unique(),df["Multi_Module"].unique())
		else:
			Multi_Module =  st.multiselect("Multi Module",
				df["Multi_Module"].unique())

	# Select Correlation
	with st.expander('AF2-OF pLDDT Pearson'):
		reset_r = st.checkbox("Reset AF2-OF pLDDT Pearson", value=True)
		if reset_r:
			minr, maxr = st.slider("AF2-OF pLDDT Pearson", 0.0, 1.0, (0.0, 1.0), 0.01)
		else:
			minr, maxr = st.slider("AF2-OF pLDDT Pearson", 0.0, 1.0, (0.25, 0.75), 0.01)
			
	# Select DALI Z
	with st.expander('AF2-OF DALI Z'):
		reset_DALI = st.checkbox("Reset AF2-OF DALI Z", value=True)
		if reset_DALI:
			min_DALI, max_DALI = st.slider("AF2-OF DALI Z", 0.0, 100.0, (0.0, 100.0), 0.1)
		else:
			min_DALI, max_DALI = st.slider("AF2-OF DALI Z", 0.0, 100.0, (2.0, 50.0), 0.1)

# ---- FILTER DATAFRAME ----

df_selection = df.query(
	"ID_No == @ID_Number"
	"& @minID <= ID_No <= @maxID"
	"& @min_len <= Length <= @max_len"
	"& @min_MW <= MW <= @max_MW"
	"& @min_pI <= pI <= @max_pI"
	"& @min_disorder <= Longest_Disorder_Region_percent <= @max_disorder"
	"& @min_MSA_Depth_Mean <= MSA_Depth_Mean <= @max_MSA_Depth_Mean"
	"& @min_pLDDT_mean <= pLDDT_mean <= @max_pLDDT_mean"
	"& @min_PDBeFold_Q <= PDBeFold_Q <= @max_PDBeFold_Q"
	"& Beta_Sheet == @Beta_Sheet"
	"& Foldedness_Rating == @Foldedness_Rating"
	"& Multi_Module == @Multi_Module"
	"& @minr <= AF2_OF_pLDDT_Pearson <= @maxr"
	"& @min_DALI <= AF2_OF_DALI_Z <= @max_DALI"
	)

st.dataframe(df_selection)
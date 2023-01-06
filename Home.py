# ---- INSTRUCTIONS ----

# Windows Command Prompt

# pip install pandas openpyxl
# pip install plotly-express
# pip install streamlit

# streamlit run Home.py

import streamlit as st

st.set_page_config(layout="wide")

st.title(":house: Effectors Dashboard")

st.write('This app is for interacting with the effector bioinformatics data. '
		 'For protein structure viewing, PyMOL is best (use grid mode if needed)')

st.markdown('---')

st.subheader('Note on protein names')

st.write('Identification number (ID No) is used to identify each protein. This means that the'
		 '\'Mp\' in the name of the protein has been removed. '
		 'MpC002, Mp92a and MIF1 have the ID No 0, 92 and 100, respectively. Mp92b is excluded.'
		 )

st.markdown('---')

st.subheader('Please navigate to the desired page (sidebar)')

# --- LINKS OF INTEREST ---
# https://discuss.streamlit.io/t/issues-with-nested-button-selectbox/28569/10
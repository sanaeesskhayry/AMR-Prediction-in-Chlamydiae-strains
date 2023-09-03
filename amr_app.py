import streamlit as st
import pandas as pd
from PIL import Image
import pickle
from Bio import SeqIO
from io import StringIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def feature_encoding(fasta_file):
    results = []

    for i, record in enumerate(SeqIO.parse(StringIO(fasta_file), "fasta")):

        analysed_protein = ProteinAnalysis(record.seq)

        # Calculate properties
        charge_at_pH = analysed_protein.charge_at_pH(7)
        # 1-
        molecular_weight = analysed_protein.molecular_weight()
        # 2-
        hydrophobicity = analysed_protein.gravy()
        # 3-
        isoelectric_point = analysed_protein.isoelectric_point()
        # 4-
        aromaticity = analysed_protein.aromaticity()
        # 5-
        instability_index = analysed_protein.instability_index()
        # 6-
        secondary_structure_fraction = analysed_protein.secondary_structure_fraction()
        # 7-
        molar_extinction_coefficient = analysed_protein.molar_extinction_coefficient()
        # 8-
        flexibility = analysed_protein.flexibility()
        # 9-
        get_amino_acids_percent = analysed_protein.get_amino_acids_percent()

        # Store the results
        results.append({
            "Sequence": i+1,
            "charge at pH =7": charge_at_pH,
            "Molecular Weight": molecular_weight,
            "Hydrophobicity": hydrophobicity,
            "Isoelectric Point": isoelectric_point,
            "Aromaticity": aromaticity,
            "Instability Index": instability_index,
            "secondary structure fraction": secondary_structure_fraction,
            "Molar extinction coefficient": molar_extinction_coefficient,
            "Flexibility": sum(flexibility) / len(flexibility),
            "Amino acids percent": get_amino_acids_percent
        })

# Process the data to create separate columns for Amino Acids Percent, Molar Extinction Coefficient, and Secondary Structure Fraction
    for entry in results:
        # Convert Amino Acids Percent dictionary into separate columns
        for aa, percent in entry['Amino acids percent'].items():
            entry[f'{aa} '] = percent

    # Separate Molar Extinction Coefficient into separate columns
        entry['MEC_reduced cysteines'], entry['MEC_disulfid bridges'] = entry['Molar extinction coefficient']

    # Separate Secondary Structure Fraction into separate columns
        entry['SSF_Helix'], entry['SSF_Turn'], entry['SSF_Sheet'] = entry['secondary structure fraction']

# Convert 'results' into a DataFrame
    df = pd.DataFrame(results)

# Drop the original columns that have been separated
    df.drop(['Amino acids percent', 'Molar extinction coefficient',
            'secondary structure fraction'], axis=1, inplace=True)

# Print the DataFrame
    return df

# Model building


def build_model(input_data, sequences):
    # Reads in saved regression model
    with open("amr_predict_model.pkl", "rb") as file:
        loaded_model = pickle.load(file)

    # Apply model to make predictions
    prediction = loaded_model.predict(input_data)

    # Initialize lists to store data
    records = SeqIO.parse(StringIO(sequences), "fasta")
    ids = []
    prediction_list = []

    # Loop through the sequences and extract information
    for seq_record, pred in zip(records, prediction):
        # Extract the sequence ID and description
        sequence_id = seq_record.id
        predict = pred

    # Append the data to lists
        ids.append(sequence_id)
        prediction_list.append(predict)

    # Create a DataFrame
    data = {
        "Sequence": ids,
        "Antimicrobial Profile": prediction_list,
    }

    df = pd.DataFrame(data)

    # Map prediction values to labels (1: resistant, 0: sensitive)
    df["Antimicrobial Profile"] = df["Antimicrobial Profile"].map(
        {1: "Resistant", 0: "Sensitive"})

    st.header('**Prediction output**')
    st.write(df)




# Page title
st.markdown("""
# Antimicrobial Activity Prediction App (Chlamydia Species)
This app leverages protein sequences to anticipate the efficacy of antimicrobial agents against diverse Chlamydia strains, facilitating informed treatment decisions and ultimately enhancing patient care.

---
""")

# Logo image
image = Image.open('logo.png')

st.image(image, use_column_width=True)


# Sidebar
with st.sidebar.header('Upload your protein sequences in a FASTA file'):
    uploaded_file = st.sidebar.file_uploader(
        "Upload your input file", type=['fasta'])
    st.sidebar.markdown("""
[Example input file](https://raw.githubusercontent.com/SanaeEsskhayry/data/master/fasta_file_example.fasta)
""")

# Process the uploaded file
if uploaded_file is not None:
    # Read the contents of the uploaded file as bytes
    content_bytes = uploaded_file.read()

    # Decode the bytes content to a string using UTF-8 encoding
    content_str = content_bytes.decode("utf-8")

    # Parse the FASTA content using BioPython
    records = SeqIO.parse(StringIO(content_str), "fasta")

    # Initialize lists to store data
    ids = []
    descriptions = []
    sequence_data = []

    # Loop through the sequences and extract information
    for seq_record in records:
        # Extract the sequence ID and description
        sequence_id = seq_record.id
        sequence_description = seq_record.description

    # Extract the sequence
        sequence = str(seq_record.seq)

    # Append the data to lists
        ids.append(sequence_id)
        descriptions.append(sequence_description)
        sequence_data.append(sequence)

    # Create a DataFrame
    data = {
        "#": ids,
        "sequence description": descriptions,
        "sequence": sequence_data
    }
    df = pd.DataFrame(data)

    # Display the parsed records and sequences in a single text area
    st.header('**Parsed FASTA records:**')
    st.write(df)

    if st.sidebar.button('Predict'):
        with st.spinner("Predicting Antimicrobial Activity..."):
            input_data = feature_encoding(content_str)
            st.header('**Calculate Protein Properties**')
            st.write(input_data)
            st.write(input_data.shape)

            # Read descriptor list used in previously built model
            st.header('**Subset of descriptors from previously built models**')
            Xlist = list(pd.read_csv('features_list.csv').columns)
            input_data_subset = input_data[Xlist]
            st.write(input_data_subset)
            st.write(input_data_subset.shape)
            # Run the Predictive Model:
            build_model(input_data_subset, content_str)


else:
    st.info('Upload input data in the sidebar to start!')

import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import requests  # To make HTTP requests

# List all available molecular descriptors
descriptor_names = [desc_name[0] for desc_name in Descriptors._descList]
calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)

# Load the trained model and other objects
best_model = joblib.load('best_model_RandomForest.pkl')
label_encoder_smiles = joblib.load('label_encoder_smiles.pkl')
scaler = joblib.load('scaler.pkl')
imputer = joblib.load('imputer.pkl')
feature_columns = joblib.load('feature_columns.pkl')

# Function to calculate all available molecular descriptors for SMILES
def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return calculator.CalcDescriptors(mol)
    else:
        return [np.nan] * len(descriptor_names)

# Function to calculate Morgan fingerprints for SMILES
def calculate_morgan_fingerprints(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        return np.array(fp)
    else:
        return np.array([np.nan] * n_bits)

def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol_wt = round(Descriptors.MolWt(mol), 2)
    log_p = round(Descriptors.MolLogP(mol), 2)
    h_bond_donors = Descriptors.NumHDonors(mol)
    h_bond_acceptors = Descriptors.NumHAcceptors(mol)
    tpsa = round(Descriptors.TPSA(mol), 2)
    return {
        'molecular_weight': mol_wt,
        'log_p': log_p,
        'hydrogen_bond_donors': h_bond_donors,
        'hydrogen_bond_acceptors': h_bond_acceptors,
        'tpsa': tpsa
    }

# Function to generate the best fragment for a given drug SMILES and other properties
def generate_best_fragment(smiles, model, feature_columns):
    fingerprints = calculate_morgan_fingerprints(smiles)
    descriptors = calculate_descriptors(smiles)
    if any(np.isnan(fingerprints)) or any(np.isnan(descriptors)):
        return "Invalid SMILES input"
    
    # Create feature vector directly using numpy
    input_data = np.concatenate([fingerprints, descriptors])
    
    # Ensure the feature vector has all required columns (zero-padding if necessary)
    input_data_full = np.zeros(len(feature_columns))
    feature_length = len(input_data)
    if feature_length <= len(input_data_full):
        input_data_full[:feature_length] = input_data
    
    # Handle missing values
    input_data_imputed = imputer.transform([input_data_full])[0]
    input_data_scaled = scaler.transform([input_data_imputed])[0]
    
    # Predict the fragment
    fragment_label = model.predict([input_data_scaled])[0]
    fragment_smiles = label_encoder_smiles.inverse_transform([int(fragment_label)])[0]
    
    return fragment_smiles

# Function to clean up a molecule
def cleanup_molecule_rdkit(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)

# Function to get 3D structure from SMILES using app.py's endpoint
def get_3d_structure(smiles):
    response = requests.post('https://frag-predict-final.onrender.com/get_3d_structure', json={'smiles': smiles})
    if response.status_code == 200:
        return response.json().get('pdb')
    else:
        print("Error fetching 3D structure:", response.json())
        return None

# Main function to process the SMILES string and get the 3D structure
def main(input_smiles):
    # Generate the best fragment
    best_fragment = generate_best_fragment(input_smiles, best_model, feature_columns)
    cleaned_fragment_smiles = cleanup_molecule_rdkit(best_fragment)
    
    if cleaned_fragment_smiles:
        print("Best fragment SMILES:", cleaned_fragment_smiles)
        pdb_block = get_3d_structure(cleaned_fragment_smiles)
        if pdb_block:
            print("3D Structure (PDB):", pdb_block)
        else:
            print("Failed to get 3D structure.")
    else:
        print("Failed to generate a valid fragment.")

# Example usage
if __name__ == "__main__":
    input_smiles = "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"  # Replace with user input
    main(input_smiles)

import torch
from transformers import RobertaTokenizer, RobertaForMaskedLM
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors
import requests
import re

def is_valid_smiles(smiles):
    """Check if the SMILES is valid."""
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None, mol

def fix_smiles(smiles):
    """Attempt to fix the SMILES string using several strategies."""
    
    # Step 1: Balance parentheses
    print("Attempting to balance parentheses...")
    balanced_smiles = balance_parentheses(smiles)
    if balanced_smiles:
        valid, mol = is_valid_smiles(balanced_smiles)
        if valid:
            return Chem.MolToSmiles(mol, canonical=True)
    
    # Step 2: Identify and fix specific errors
    print("Attempting to identify and fix specific errors...")
    corrected_smiles = correct_smiles(smiles)
    if corrected_smiles:
        valid, mol = is_valid_smiles(corrected_smiles)
        if valid:
            return Chem.MolToSmiles(mol, canonical=True)
    
    print("All attempts to fix the SMILES failed.")
    return None

def balance_parentheses(smiles):
    """Balance the number of opening and closing parentheses."""
    open_paren_count = smiles.count('(')
    close_paren_count = smiles.count(')')
    
    if open_paren_count > close_paren_count:
        # Add missing closing parentheses at the end
        smiles += ')' * (open_paren_count - close_paren_count)
    elif close_paren_count > open_paren_count:
        # Remove extra closing parentheses from the end
        smiles = smiles[::-1].replace(')', '', close_paren_count - open_paren_count)[::-1]
    
    return smiles

def correct_smiles(smiles):
    """Correct specific errors such as extra characters or mismatched ring bonds."""
    smiles_list = list(smiles)
    corrected = False

    # Check for unbalanced parentheses and correct
    paren_stack = []
    for i, char in enumerate(smiles_list):
        if char == '(':
            paren_stack.append(i)
        elif char == ')':
            if paren_stack:
                paren_stack.pop()
            else:
                # Found an unmatched closing parenthesis
                print(f"Removing unmatched closing parenthesis at position {i}")
                smiles_list[i] = ''  # Remove the extra closing parenthesis
                corrected = True
    
    # If there are unmatched opening parentheses, remove them
    while paren_stack:
        index = paren_stack.pop()
        print(f"Removing unmatched opening parenthesis at position {index}")
        smiles_list[index] = ''
        corrected = True
    
    # Check and correct unclosed rings
    ring_numbers = re.findall(r'\d', ''.join(smiles_list))
    if len(ring_numbers) % 2 != 0:
        print("Unclosed ring detected.")
        unclosed_ring_idx = find_unclosed_ring_index(smiles_list, ring_numbers)
        if unclosed_ring_idx is not None:
            print(f"Removing unclosed ring at position {unclosed_ring_idx}")
            smiles_list[unclosed_ring_idx] = ''  # Remove the unclosed ring
            corrected = True

    if corrected:
        return ''.join(smiles_list)
    
    # If no correction was made, return None
    return None

def find_unclosed_ring_index(smiles_list, ring_numbers):
    """Find the index of the unclosed ring number in the SMILES list."""
    ring_dict = {}
    for i, char in enumerate(smiles_list):
        if char.isdigit():
            if char in ring_dict:
                del ring_dict[char]  # Ring is closed, remove from dictionary
            else:
                ring_dict[char] = i  # Unclosed ring found
    if ring_dict:
        # Return the position of the first unclosed ring found
        return list(ring_dict.values())[0]
    return None

# Example usage:
invalid_smiles = "CC=CC1CC)C(=C)O)NCO"
valid_smiles = fix_smiles(invalid_smiles)

if valid_smiles:
    print(f"Closest valid SMILES: {valid_smiles}")
else:
    print("Could not convert to a valid SMILES.")

#======

# Function to clean up a molecule
def cleanup_molecule_rdkit(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)

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

def get_3d_structure(smiles):
    response = requests.post('/get_3d_structure', json={'smiles': smiles})
    if response.status_code == 200:
        return response.json().get('pdb')
    else:
        print("Error fetching 3D structure:", response.json())
        return None
    
# Function to check if a SMILES string is valid
# def is_valid_smiles(smiles):
#     return Chem.MolFromSmiles(smiles) is not None

# # Function to correct common issues in SMILES
# def is_balanced(smiles):
#     return smiles.count('(') == smiles.count(')')

# def remove_invalid_characters(smiles):
#     allowed_chars = set('CNOFPScnBrios1234567890-=#()/\\@[]+%')
#     return ''.join([char for char in smiles if char in allowed_chars])

# def correct_ring_closures(smiles):
#     ring_digits = re.findall(r'\d', smiles)
#     for digit in ring_digits:
#         if smiles.count(digit) % 2 != 0:
#             smiles = smiles[::-1].replace(digit, '', 1)[::-1]
#     return smiles

# def correct_smiles(smiles):
#     if not is_balanced(smiles):
#         open_parens = smiles.count('(')
#         close_parens = smiles.count(')')
#         if open_parens > close_parens:
#             smiles += ')' * (open_parens - close_parens)
#         elif close_parens > open_parens:
#             smiles = smiles.lstrip(')' * (close_parens - open_parens))

#     smiles = remove_invalid_characters(smiles)
#     smiles = correct_ring_closures(smiles)
#     smiles = smiles.strip(')(')

#     return smiles

def predict_fragment_smiles(smiles, protein, max_length=128):

    model = RobertaForMaskedLM.from_pretrained('model-'+str(protein))  # Update with your model path
    tokenizer = RobertaTokenizer.from_pretrained('tokenizer-'+str(protein))  # Update with your tokenizer path
    model.eval()

    inputs = tokenizer(smiles, max_length=max_length, padding='max_length', truncation=True, return_tensors="pt")
    with torch.no_grad():
        outputs = model(input_ids=inputs['input_ids'], attention_mask=inputs['attention_mask'])
    logits = outputs.logits
    predicted_ids = torch.argmax(logits, dim=-1)
    predicted_smiles = tokenizer.decode(predicted_ids[0], skip_special_tokens=True)
    predicted_smiles = fix_smiles(predicted_smiles)
    print("predicted smiles: ", predicted_smiles)
    # If predicted SMILES is invalid, attempt to correct it
    if not is_valid_smiles(predicted_smiles):
        predicted_smiles = correct_smiles(predicted_smiles)
        if not is_valid_smiles(predicted_smiles):
            predicted_smiles = "INVALID"  # Placeholder if correction fails

    return predicted_smiles

def tanimoto_similarity(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 is None or mol2 is None:
        return 0.0
    
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    
    return DataStructs.TanimotoSimilarity(fp1, fp2)

# Example usage
new_drug_smiles = "CC=C(C)C(=O)OC1C(C)=CC23C(=O)C(C=C(COC(C)=O)C(O)C12O)C1C(CC3C)C1(C)C"  # Replace with your input SMILES
predicted_fragment_smiles = predict_fragment_smiles(new_drug_smiles, 'mTOR')

print("Predicted Fragment SMILES:", predicted_fragment_smiles)

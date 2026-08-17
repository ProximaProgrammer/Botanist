#!/usr/bin/env python
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import AllChem,Descriptors
import argparse, subprocess, sys, os
from pathlib import Path
from tkinter import Tk,filedialog
import absorptivity

#reminder: use debug mode and move line by line if needed
'''
Chemical – Name of the chemical
Plant Part – source (leaf, root, seed, etc.)
Compound_ID – Unique identifier for each compound
SMILES – stores keyboard-typable bond/connection topology; 3D structure up to steroisomerism
InChIKey – Standardized chemical identifier
Molecular_Weight – Weight of the compound in Daltons/amu (1.6605e-27 kg)
LogP – Partition coefficient (hydrophobicity measure)
HBD – Number of hydrogen bond donors
HBA – Number of hydrogen bond acceptors
Rotatable_Bond – Number of rotatable bonds

[11 more columns below]
'''
#columns with only one value that don't tell us anything: ["Bioactivity_Label","Lipinski_OK","Has_Reactive"]

Z = {"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,
          "Na":11,"Mg":12,"Al":13,"Si":14,"P":15,"S":16,"Cl":17,"Ar":18,"K":19,"Ca":20,"Sc":21,"Ti":22,"V":23,"Cr":24,"Mn":25,"Fe":26,"Co":27,"Ni":28,"Cu":29,"Zn":30,"Ga":31,"Ge":32,"As":33,"Se":34,"Br":35,"Kr":36,"Rb":37,"Sr":38,"Y":39,"Zr":40,"Nb":41,"Mo":42,"Tc":43,"Ru":44,"Rh":45,"Pd":46,"Ag":47,"Cd":48,"In":49}

def analyze_conjugation(mol):
  assert isinstance(mol, Chem.Mol), "Error: The input of analyze_conjugation should be a mol object." #if accepting SMILES: mol = Chem.MolFromSmiles(smiles)

  # Ensure the molecule is sanitized to compute conjugation flags accurately
  AllChem.SanitizeMol(mol)
  #tag and store conjugation status for every bond
  conjugated_bonds = []
  for bond in mol.GetBonds():
    is_conj = bond.GetIsConjugated()
    bond.SetBoolProp("is_conjugated", is_conj) #adding a custom property is_conjugated on the bond object for downstream use
    if is_conj:
      conjugated_bonds.append(bond)
  
  conjugated_atom_indices = set() #extract unique atom order indices participating in any conjugated bond
  for bond in conjugated_bonds:
    conjugated_atom_indices.add(bond.GetBeginAtomIdx())
    conjugated_atom_indices.add(bond.GetEndAtomIdx())
  atom_pi_electrons = [] #record the pi electrons contributed by these specific atoms
  for idx in conjugated_atom_indices:
    atom = mol.GetAtomWithIdx(idx)
    atom_pi_electrons.append([atom.GetIdx(), Chem.GetNumPiElectrons(atom)])

  output_information = []
  for bond in mol.GetBonds():
    a1 = bond.GetBeginAtom().GetSymbol()
    a2 = bond.GetEndAtom().GetSymbol()
    idx1 = bond.GetBeginAtomIdx()
    idx2 = bond.GetEndAtomIdx()
    is_conj = int(bond.GetBoolProp("is_conjugated"))
    output_information.append([idx1,idx2,Z[a1],Z[a2],is_conj+1]) #atom order/index 1 & 2, element 1 & 2 atomic numbers, bond status (0/1/2) ---- we make all five ints for easy C++ handling (an array can only have one type)
    
  return [mol.GetNumAtoms(),atom_pi_electrons,output_information]

# --- Example Usage
# example_smiles = Chem.MolFromSmiles("C1=CC(=CC=C1C[C@@H](C(=O)O)N)O") #Tyrosine
# print(analyze_conjugation(example_smiles))
# print(f"Total electrons involved in conjugation: {electron_count}\n")
# print("Bond Conjugation Breakdown:")

# for bond in processed_mol.GetBonds():
#   a1 = bond.GetBeginAtom().GetSymbol()
#   a2 = bond.GetEndAtom().GetSymbol()
#   idx1 = bond.GetBeginAtomIdx()
#   idx2 = bond.GetEndAtomIdx()
#   is_conj = bond.GetBoolProp("is_conjugated")
#   print(f"Bond {idx1:>2}({a1})-{idx2:<2}({a2}) | Conjugated: {is_conj}") #WARNINGS: atoms are numbered 0-indexed!

class NotFoundError(Exception):
    def by_name(self):
        sys.stderr.write("Error: Compound not found. Please try alternative names or check the spelling/formatting. If sure your entry was correct, the dataset used may not include it.\nSee https://www.kaggle.com/datasets/yashasvigoswami/phytochemicals")
    def by_id(self):
        sys.stderr.write("Error: Compound not found. Please check the ID entered. If sure your entry was correct, the dataset used may not include it.\nSee https://www.kaggle.com/datasets/yashasvigoswami/phytochemicals")


df = pl.read_csv(Path(__file__).resolve().parent / "prefphytochemicals.csv", has_header=True) #35254 rows, exluding header row, 21 columns
df = df.unique(subset=["Compound_ID"])

def analyze_chemical(compound_id): #make compound name an alternative argument
    try:
        compound_id = "CID_"+str(int(compound_id)) #tests if int was input as required
        n_row = np.argmax(df["Compound_ID"].to_numpy() == compound_id)
        if n_row==0:
            raise NotFoundError
        else:
            SMILES = df.row(n_row)[3] #not using Standard_SMILES since harder to read and have to do extra inferred constructions
            mol = Chem.MolFromSmiles(SMILES)
            mol = Chem.AddHs(mol) #adding implicit hydrogens
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            AllChem.MMFFOptimizeMolecule(mol)

            try:
                total_electrons = sum(atom.GetAtomicNum() for atom in mol.GetAtoms())
                valence_electrons = Descriptors.NumValenceElectrons(mol)
                core_electrons = total_electrons-Descriptors.NumValenceElectrons(mol)
                print("Electrons: total ",total_electrons," | valence ",valence_electrons," | core ",core_electrons)

                #sys.path.append(os.path.abspath(os.path.dirname(__file__)))
                analyze_conj_pair = analyze_conjugation(mol)
                return absorptivity.compute_spectrum(SMILES, analyze_conj_pair[0], analyze_conj_pair[1], analyze_conj_pair[2]) #test
                # Here we call the C++ module and eventually return an absorptivity for the compound, before having one last calculation here (in v1) to get the reflectance spectrum--which should appear continuous with small dips corresponding to absorbance lines with their respective amplitudes.
                # then a graph or two is displayed by analyze_chemical, and if verbosity=Trcue, numeric data is open as well
            except subprocess.CalledProcessError as e:
                print(f"Command failed with exit code {e.returncode}")
                print("Error logs:")
                print(e.stderr)
    except ValueError:
        print("Invalid: must enter a numeric ID (without quotes)") #replace this with an argparse custom error
    except NotFoundError as e:
        e.by_id()

class ARGPARSE_analyze_chemical_id(argparse.Action): #inheriting from the argparse.Action class to run analyze_chemical on the input after the -id flag
   def __call__(self, parser, namespace, input_id, option_string=None):
        result = analyze_chemical(input_id)
        setattr(namespace, self.dest, result) #result saved to args.id
        print("Output of analyze_chemical:",result) #temporary
class ARGPARSE_analyze_chemical_name(argparse.Action):
    def __call__(self, parser, namespace, input_name, option_string=None):
        try:
            name = input_name.replace(" ","-").capitalize() #standardizing the input name to match the dataset
            if name.count(",")>0:
                name = '"' + name + '"' #to not count comma as delimiter, csv enclosed names with commas in quotes
            candidates = df.filter(df.col("Chemical")==name)
            if not candidates.isEmpty():
                input_id = candidates.select("Compound_ID").to_numpy()[0][0]
                result = analyze_chemical(input_id)
                setattr(namespace, self.dest, result)
                print("Output of analyze_chemical:",result) #temporary
            else:
               raise NotFoundError
        except NotFoundError as e:
            e.by_name()


parser = argparse.ArgumentParser(description="Plot absorption/reflection spectrum of phytocompound. Predict plausible concentrations of phytocompounds from observed color of a plant part.") #some primary plant compounds like chlorophyll are excluded from the dataset!

parser.add_argument("--analyze_compound", metavar="CID", nargs="*", type=int, help="Enter CID [digits only]. An approximate spectrum for 120-800 nm is produced and displayed. Multiple CIDs may be input, and all graphs will be displayed.")  
parser.add_argument("-hide", action="store_true", help="Hides graphs. Must use the -v flag so numeric data is output instead.") #return an error if --analyze_compound and --hide are used but not -v
parser.add_argument("-v", "--verbose", action="store_true", help="Lists graph-equivalent numeric data in tabular format. If multiple compounds were input, then each will have their data displayed.")
parser.add_argument("-n", "--name", type=str, action=ARGPARSE_analyze_chemical_name, help="Enter name of compound. If not found, try alternative names.")
parser.add_argument("-id", type=int, action=ARGPARSE_analyze_chemical_id, help="Enter CID/compound ID (digits only): ")

class FileSelectAction(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        # Hide the main Tkinter root window
        root = Tk()
        root.withdraw()

        # Open the file selection dialog
        selected_path = filedialog.askopenfilename(
            title="Select Input File",
            filetypes=[("Image Files", "*.png *.jpg *.jpeg *.bmp"), ("All Files", "*.*")]
        )

        # Store the result as a Path object if a file was chosen, otherwise None
        setattr(namespace, self.dest, Path(selected_path) if selected_path else None)

parser.add_argument(
    "--analyze_plant",
    nargs=0,  # tells argparse not to accept command-line string values
    action=FileSelectAction,
    help="Open a file selection window to choose the input file."
)

args = parser.parse_args()
if args.verbose:
   print("Numeric data below. As it is still cached, you may download it as a csv or txt with the following command: ???")
# args.log.write("argparse log written~")  
# args.log.close() #need to run `chmod +x main.py` and the shebang on the first line to enable terminal execution of file w/ argparse
#run `python main.py` with flags after to see the argparse in action

# Use the safely stored path later in your script (do we need to insert the below block between the first and other two lines immediately above?)
# if args.input_file:
#     with args.input_file.open("r") as f:
#         print(f"Opened file: {args.input_file}") #just a preview of the img, so remove this block of code after testing it, right?

# from PIL import Image, ImageOps #posterization (color discretization) for analysis (if chemicals are unevenly distributed, we can more easily see a subset of them in differently colored regions, hence allowing for better analysis)
# posterized_img = ImageOps.posterize(args.input_file, 2) #second parameter is basically resolution #do we replace the first parameter with 'f' from the block above?
# posterized_img.save(f"./posterized_images/{args.input_file.(name)}.png")

#------------------user search options------------------ (display with argparse --help)
# by name : if there's commas, add quotes to start and end, and capitalize anyway, and replace spaces with '-', before searching by compound name. Say to use alternative names if not found. Issue 'not found' custom warning message
# by id : ask user to 'enter CID [digits only]', then append numeric entry to the string "CID_"




''' User Settings (JSON) -- original:
{
    "workbench.colorTheme": "Abyss",
    "liveServer.settings.donotShowInfoMsg": true,
    "files.autoSave": "afterDelay",
    "liveServer.settings.port": 5500,
    "python.createEnvironment.trigger": "off",
    "workbench.secondarySideBar.defaultVisibility": "hidden",
    "cmake.cmakePath": "/Applications/CMake.app/Contents/bin/cmake",
    "github.copilot.enable": {
        "*": false,
        "plaintext": false,
        "markdown": false,
        "scminput": false
    }
}

also, you added a line to ~/.zshrc
'''
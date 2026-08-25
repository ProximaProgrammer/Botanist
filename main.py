#!/usr/bin/env python
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
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
  Chem.Kekulize(mol, clearAromaticFlags=False) #explicit alternating single/double bonds instead of an ambiguous 1.5 "aromatic" bond order, so C++ gets a real integer bond order without having to re-parse the SMILES string itself; clearAromaticFlags=False keeps GetIsAromatic() available if needed later
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
    order_code = int(round(bond.GetBondTypeAsDouble())) #explicit bond order (1/2/3) now that the molecule is Kekulized -- C++ no longer needs to reconstruct this from the SMILES string
    in_ring = int(bond.IsInRing())
    output_information.append([idx1,idx2,Z[a1],Z[a2],is_conj+1,order_code,in_ring]) #atom order/index 1 & 2, element 1 & 2 atomic numbers, bond status (0/1/2), bond order, bond-in-ring ---- ints for easy C++ handling (an array can only have one type)

  atom_information = [] #idx, atomic number, total (implicit+explicit) H count, formal charge, atom-in-ring -- covers EVERY atom (including unbonded single-atom species like a bare [Zn]), unlike output_information which only covers bonded atoms
  for atom in mol.GetAtoms():
    atom_information.append([atom.GetIdx(), atom.GetAtomicNum(), atom.GetTotalNumHs(), atom.GetFormalCharge(), int(atom.IsInRing())])

  return [mol.GetNumAtoms(),atom_pi_electrons,output_information,atom_information]

def convolute(peaks_df, min_v=int(1.0e7/800), max_v=1+int(1.0e7/120), step=1, HWHM=800, shape="Gaussian"): #HWHM = half-width at half maximum (wavenumber difference of 800 for around 250 nm corresponds to 5 nm), and should be determined by uncertainty (For analyze_chemical, some clever measure of spectrum difference from true data. For analyze_plant, sensor error)
    #get flattened np arrays
    peaks = peaks_df["wavenumber"].to_numpy()
    amplitudes = peaks_df["amplitude"].to_numpy()

    v_grid = np.arange(min_v, max_v, step) #wavenumber axis grid with resolution determined by `step`
    relative_position = v_grid[:, None] - peaks[None, :] #this hack automatically shapes out a matrix with shape (len(x_grid),len(peaks)) so there's an axis/grid column for each peak
    line_shapes = np.exp(-np.log(2) * (relative_position/HWHM)**2)
    total_absorption = np.dot(line_shapes, amplitudes) #adds an amplitude-weighted sum of gaussian intensity grids for each peak, superposing to form the displayed spectrum
    # normalizing_factor = ((np.log(2)/np.pi)**0.5)/HWHM #multiply line_shapes by this to maintain same area, treating original peaks as delta functions with given amplitudes
    
    return pl.DataFrame({"wavenumber": v_grid, "amplitude": total_absorption}) #reminder: this is (still) only *relative* amplitude since we didn't even normalize above
    #print(TEMP.sort("amplitude", descending=True).head(20))

class NotFoundError(Exception):
    def by_name(self):
        sys.exit("Error: Compound not found. Please try alternative names or check the spelling/formatting. If sure your entry was correct, the dataset used may not include it.\nSee https://www.kaggle.com/datasets/yashasvigoswami/phytochemicals")
    def by_id(self):
        sys.exit("Error: Compound not found. Please check the ID entered. If sure your entry was correct, the dataset used may not include it.\n( See https://www.kaggle.com/datasets/yashasvigoswami/phytochemicals )")

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
                return absorptivity.compute_spectrum(SMILES, analyze_conj_pair[0], analyze_conj_pair[1], analyze_conj_pair[2], analyze_conj_pair[3]) #test
                # Here we call the C++ module and eventually return an absorptivity for the compound, before having one last calculation here (in v1) to get the reflectance spectrum--which should appear continuous with small dips corresponding to absorption lines with their respective amplitudes.
                # then a graph or two is displayed by analyze_chemical, and if verbosity=True, numeric data is open as well
            except subprocess.CalledProcessError as e:
                print(f"Command failed with exit code {e.returncode}")
                print("Error logs:")
                print(e.stderr)
    except ValueError:
        print("Invalid: must enter a numeric ID without quotes")
    except NotFoundError as e:
        e.by_id()

def display_single(results):
    x,y = zip(*results)
    results = pl.DataFrame({"wavenumber": x, "amplitude": y})
    results = convolute(results)

    plt.figure()
    plot = plt.plot(1.0e7/results["wavenumber"], results["amplitude"], lw=1)
    plt.xlabel("wavelength (nm)")
    plt.ylabel("(relative) amplitude")
    plt.xlim(120,800) #our wavelength range
    plt.show()

def display_multiple(results):
    print("[this function hasn't been written yet]")
    for result in results:
        print("[add to same plot and add key, then display the same plot for everyone]")

class ARGPARSE_analyze_chemical_id(argparse.Action): #inheriting from the argparse.Action class to run analyze_chemical on the input after the -id flag
   def __call__(self, parser, namespace, input_id, option_string=None):
        result = analyze_chemical(input_id)
        setattr(namespace, self.dest, result) #result saved to args.id
        display_single(result)
class ARGPARSE_analyze_chemical_name(argparse.Action):
    def __call__(self, parser, namespace, input_name, option_string=None):
        try:
            name = input_name[0]
            if name.count("-")>0: #e.g. (+)-DIHYDROCARVEOL     (TRYING TO FIND A (QUOTE-LIKE) CHARACTER THAT DOESN'T TRIGGER TERMINAL COMMAND RESPONSE)
                if name.count("'")==2:
                    name = name.split("'")[1] #removes quotation marks, which were only necessary for command processing
                else:
                    print(name, " quote count: ", name.count("'"))
                    sys.exit("If your compound name has dashes, surround it with dollar signs (`) and try again")
            name = name.replace(" ","-").upper() #standardizing the input name to match the dataset
            if name.count(",")>0:
                name = '"' + name + '"' #to not count comma as delimiter, csv enclosed names with commas in quotes
            candidates = df.filter(pl.col("Chemical")==name)
            if not candidates.is_empty():
                input_id = candidates.select("Compound_ID").to_numpy()[0][0]
                input_id = input_id.split("_")[1] #since the format is CID_12345
                result = analyze_chemical(input_id)
                setattr(namespace, self.dest, result)
                display_single(result)
            else:
                raise NotFoundError
        except NotFoundError as e:
            e.by_name()
class ARGPARSE_multiple_chemical_ids(argparse.Action):
    def __call__(self, parser, namespace, input_ids, option_string=None):
        results = []
        for input_id in input_ids:
            try:
                result = analyze_chemical(input_id)
                results.append(result)
            except NotFoundError as e:
                e.by_id()
        setattr(namespace, self.dest, results) #results saved to args.ids
        display_multiple(results)
class ARGPARSE_multiple_chemical_names(argparse.Action):
    def __call__(self, parser, namespace, input_names, option_string=None):
        results = []
        for name in input_names:
            try:
                if name.count("-")>0:
                    if name.count('"')==2:
                        name = name.split('"')[1]
                    else:
                        sys.exit("If some of your compound names have dashes, surround them with double quotation marks and try again")
                name = name.replace(" ","-").upper()
                if name.count(",")>0:
                    name = '"' + name + '"'
                candidates = df.filter(pl.col("Chemical")==name)
                if not candidates.is_empty():
                    input_id = candidates.select("Compound_ID").to_numpy()[0][0]
                    input_id = input_id.split("_")[1]
                    result = analyze_chemical(input_id)
                    results.append(result)
                else:
                    raise NotFoundError
            except NotFoundError as e:
                e.by_name()
        setattr(namespace, self.dest, results) #results saved to args.names
        display_multiple(results)


parser = argparse.ArgumentParser(description="Plot absorption/reflection spectrum of phytocompound. Predict plausible concentrations of phytocompounds from observed color of a plant part.") #some primary plant compounds like chlorophyll are excluded from the dataset!

parser.add_argument("-analyze_compound", action="store_true", help="An approximate absorption/reflection spectrum for 120-800 nm is computed for a single compound.")  
parser.add_argument("-hide", action="store_true", help="Hides graphs. Must download data in a specified format instead.")
parser.add_argument("-d", "--download", nargs=1, choices=["csv", "txt", "json"], help="Downloads graph-equivalent numeric data in tabular format. After this flag, specify 'csv', 'txt', or 'json'")
parser.add_argument("-n", "--name", type=str, nargs=1, action=ARGPARSE_analyze_chemical_name, help="Enter name of a compound in the dataset. If not found, try alternative names.")
parser.add_argument("--names", type=str, nargs="+", action=ARGPARSE_multiple_chemical_names, help="Enter the names of at least two compounds in the dataset. If not found, try alternative names.")
parser.add_argument("-id", metavar="CID", type=int, nargs=1, action=ARGPARSE_analyze_chemical_id, help="Enter CID/compound ID (digits only): ")
parser.add_argument("-ids", type=int, nargs="+", action=ARGPARSE_multiple_chemical_ids, help="Enter CID/Compound ID (digits only) of at least two compounds in the dataset.")
parser.add_argument("-spectrum_type", choices=["absorption", "reflection"], help="Choose whether to compute molar [absorption] spectrum or [reflection] (observed) spectrum of compound.")
class FileSelectAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        root = Tk()
        root.withdraw() #hide the main Tkinter root window

        # open the file selection dialog
        selected_path_csv = filedialog.askopenfilename(
            title="Select (wavelength,reflectance) .csv file",
            filetypes=[("CSV Files", "*.csv")]
        )
        selected_path_pkl = filedialog.askopenfilename(
            title="Select calibration .pkl file",
            filetypes=[("Pickle Files", "*.pkl")]
        )
        root.destroy()

        #store the results as a Path object if a file was chosen, otherwise None
        setattr(namespace, "csv_path", Path(selected_path_csv) if selected_path_csv else None)
        setattr(namespace, "pkl_path", Path(selected_path_pkl) if selected_path_pkl else None)

        #the colors and locations (and features?) in the image will be used to compute the estimated percentage concentrations of phytocompounds in the plant part, which will be output as a table and/or be graphically displayed.
        #plant_analyzer.py will handle this with a pipeline of imports already built by others
parser.add_argument("--analyze_plant", nargs=0, action=FileSelectAction, help="Opens a file selection window to choose the .csv and .pkl plant spectrometry files. Estimated percentage concentrations of phytocompounds in the plant part will be computed") # nargs=0 tells argparse not to accept command-line string values

args = parser.parse_args()

has_compound = bool(args.analyze_compound or args.id or args.ids or args.name or args.names) #boolean casting to safely handle Nonetype
has_plant = bool(getattr(args, "csv_path", None) and getattr(args, "pkl_path", None))

if (has_compound ^ has_plant):
    print("Plant spectrometry files successfully received.")
else:
    sys.exit("Error: exactly one of --analyze_compound or --analyze_plant must be used.")

if has_plant and args.spectrum_type:
    sys.exit("Error: can currently only specify spectrum type if analyze compound, not when analyzing plant.")
elif has_compound and not bool(args.spectrum_type):
    sys.exit("Error: must specify spectrum type.")

def if_displaying():
    if has_compound:
            if args.spectrum_type=="absorption":
                if args.id or args.name: #only single compound
                    print("[this should display a graph for the single compound specified]")
                elif args.ids or args.names: #multiple compounds
                    print("[this should display a graph with a superposed curve for each compound in the list]")
                else:
                    sys.exit("Error: no compound or plant was specified")
            elif args.spectrum_type=="observed":
                #this will download the reflectance spectrum data instead
                print("[this should graph the observed (reflectance) spectrum data. Remove print statements after debugging and verifying downloads work]")
    elif has_plant: #checking if the attached image path is anything other than None or nonexistent
        SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
        plant_analyzer_path = os.path.join(SCRIPT_DIR, "plant_analyzer.py")
        command = [sys.executable, plant_analyzer_path, "-file_paths", str(args.csv_path), str(args.pkl_path)] #printing rather than downloading for now
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        print(result.stdout)
    else:
        print("[this shouldn't be printing! fix logic branching to return an error earlier]")

def if_downloading():
    if args.download: #NOTE: check that args.download = False when nothing selected?
        print("Downloading numeric data in the specified format:" ,args.download)
        if has_compound:
            if args.spectrum_type=="absorption":
                if args.id or args.name: #only single compound
                    print("[this should download a file for the single compound specified]")
                elif args.ids or args.names: #multiple compounds
                    print("[this should download a file for each compound in the list]")
                else:
                    sys.exit("Error: no compound or plant was specified")
            elif args.spectrum_type=="observed":
                #this will download the reflectance spectrum data instead
                print("[this should download the observed (reflectance) spectrum data. Remove print statements after debugging and verifying downloads work]")
        elif has_plant: #checking if the attached image path is anything other than None or nonexistent
            command = [sys.executable, "plant_analyzer.py", "-file_paths", str(args.csv_path), str(args.pkl_path), "-download_type", args.download]
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            print(result.stdout)
        else:
            print("[this shouldn't be printing! fix logic branching to return an error earlier]")

if not args.hide:
    print("--displaying--")
    if_displaying()
if args.download:
    print("--downloading--")
    if_downloading()
if args.hide and not args.download:
    sys.exit("Error: cannot choose to avoid both graphic display and numeric data download.") #return an error if user neither downloads or wants to display output
        

#print(args.id, args,ids, args.name, args.names) #replace all the temporary print statements with this

# Use the safely stored path later in your script (do we need to insert the below block between the first and other two lines immediately above?)
# if args.input_file:
#     with args.input_file.open("r") as f:
#         print(f"Opened file: {args.input_file}") #just a preview of the img, so remove this block of code after testing it, right?

# from PIL import Image, ImageOps #posterization (color discretization) for analysis (if chemicals are unevenly distributed, we can more easily see a subset of them in differently colored regions, hence allowing for better analysis)
# posterized_img = ImageOps.posterize(args.input_file, 2) #second parameter is basically resolution #do we replace the first parameter with 'f' from the block above?
# posterized_img.save(f"./posterized_images/{args.input_file.(name)}.png")
# OUTPUT OF C++: std::vector<std::pair<int,double>> cache; (wavenumber in cm^-1, amplitude)

# display user search options with argparse --help)
# NOTE: convolute spectrum with appropriate-width kernel?




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
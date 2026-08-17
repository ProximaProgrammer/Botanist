#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> //crucial: makes python lists and c++ vectors compatible

/*
1. build library in c++ that can be used in python (use extern C modifiers and modify CMake settings; see below)
2. filter dataset with polars in python
3. utilize c++ library of computational commands to perform analysis (and output things like spectrography plots?)
4. add an argparse (help and info) so that you can use terminal commands with the project!-------e.g. ```UV_plot -name "Jacaranda mimosifolia" -verbose```

Because C++ supports function overloading, the compiler secretly changes function names, so use the extern "C" modifier to your functions, which tells the compiler to keep the name exactly as written
c.f. https://share.google/aimode/lfeUkixpX63JgVRDW
For plants, add in (a properly weighted) chlorophyll a and b absorption spectrum. c.f.3 https://www.sciencedirect.com/science/article/pii/S0034425708000813
*/

std::vector<std::vector<int>> groupings(const std::vector<std::vector<int>>& adj_matrix, const int num_atoms, const int desired_val) {
    std::vector<std::vector<int>> groups; 
    std::vector<bool> visited(num_atoms, false); //keep outside to preserve state

    for (int i=0; i<num_atoms; i++) {
        if (!visited[i]) {
            std::vector<int> current_group;
            std::vector<int> stack = {i};

            while (!stack.empty()) {
                int atom_idx = stack.back();
                stack.pop_back();

                if (!visited[atom_idx]) {
                    visited[atom_idx] = true;
                    current_group.push_back(atom_idx);

                    //scan neighbors in the adjacency matrix
                    for (int j=0; j<num_atoms; j++) {
                        if (adj_matrix[atom_idx][j] == desired_val && !visited[j]) {
                            stack.push_back(j);
                        }
                    }
                }
            }
            groups.push_back(std::move(current_group)); //Move optimization
        }
    }
    return groups;
}


//only accepting analysis in the wavelength range where individual atom absorptivities are negligible (and other data may be flawed anyway in the vacuum UV range) : 120-800nm (far UV - very near IR)
std::vector<int> cache; //stores memory from last input

std::vector<int> molecule_spectrum(const std::string& SMILES, const size_t& num_atoms, const std::vector<std::pair<int,int>>& pi_electrons_ordered, const std::vector<std::array<int,5>>& bond_info) { //returning a vector of ints since wavenumber is a large number for wavelength << 1 cm
    cache.clear(); //clearing cache and memory for new decomposition
    cache.shrink_to_fit(); 

    //Step 1: Chromophore division
    //let's try the adjacency matrix method for now, otherwise memory-optimized tree method
    std::vector<std::vector<int>> adj_matrix(num_atoms, std::vector<int>(num_atoms, 0)); //nothing = 0, bonds = 1, conj bonds = 2

    //NOTE: atom indices are 0-indexed!!!!!!
    for (const auto& bond : bond_info) {
        int idx1 = bond[0];
        int idx2 = bond[1];
        int bond_status = bond[4]; //bond status (nonexistent=0,exists=1,conjugated=2)
        adj_matrix[idx1][idx2] = bond_status;
        adj_matrix[idx2][idx1] = bond_status; //symmetrical entry since bond is bidirectional
    }

    std::vector<std::vector<int>> chromophores = groupings(adj_matrix, num_atoms, 2); //first we identify chromophores
    std::vector<std::vector<int>> remainder_groups = groupings(adj_matrix, num_atoms, 1); // next we identify remainder groups between the chromophores
    //NOTE: the above are sets of unordered atom indices, but we don't need to sort them (yet)

    //Step 2: Auxochrome identification
    //the auxochromes are identified, and their position relative to the carbonyl group in a chromophore matters
    //NOTE: auxochromes generally AREN'T included in chromophores but only at most share an atom
    /*

    auxochrome list and primary modification algorithm (check that there's both multipliers and shifts)
    https://github.com/CompPhotoChem/ChromoPredict/blob/main/src/chromopredict/woodward_fieser.py

    secondary modification algorithm (structural features identified, incremental changes based on both chromophores and substituents/extensions))
    https://github.com/CompPhotoChem/ChromoPredict/blob/main/src/chromopredict/strucfeatures.py

    [solvent correction values]
    https://github.com/CompPhotoChem/ChromoPredict/blob/main/src/chromopredict/solvent.py
    */


    //Step 3: Modified spectra calculation and superposition
    cache = chromophores[0]; //placeholder
    return cache;
}

// inputs in bond_info: atom order/index 1 & 2, element 1 & 2 atomic numbers, bond status
/* Tests:
BENZENE / C1=CC=CC=C1
WATER,Rhizome / O
TRANS-LINALOOL-OXIDE / C[C@]1(CC[C@H](O1)C(C)(C)O)C=C
ZINC,Root,CID_23994,[Zn]

after testing everything, remove unnecessary input parameters (above, below in the PYBIND11, and in the main.py file)
*/


PYBIND11_MODULE(absorptivity, m) {
    m.doc() = "Absorptivity calculator (v0)";
    m.def("compute_spectrum", [](const std::string SMILES, const int& num_atoms, const std::vector<std::pair<int,int>> pi_electrons_ordered, const std::vector<std::array<int,5>> bond_info) {
        return molecule_spectrum(SMILES,num_atoms,pi_electrons_ordered,bond_info);
    }, "Decomposes a molecule into chromophores and remainder groups. It then identifies auxochromes within each chromophore, calculates the modified spectra, and superposes to obtain the molecule's overall spectrum for 120-800 nm. Line positions and amplitudes are respected but not lineshape.");
}

/*
to remove library: rm [library name]*.so (replace 'rm' with 'find . -name' to check there is only one result first)
to match python version with that of conda (if outside build, remove the ..): cmake -DPYTHON_EXECUTABLE=$(which python3) ..

to update library:

rm absorptivity*.so
clang++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup \
  -arch arm64 -arch x86_64 \
  $(python3-config --includes) \
  $(python3 -m pybind11 --includes) \
  absorptivity.cpp \
  -o absorptivity$(python3-config --extension-suffix)

*/
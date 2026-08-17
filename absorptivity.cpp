#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <set>
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

struct BondGraph {
    std::vector<int> atomic_num;
    std::vector<int> pi_electrons;
    std::vector<std::vector<int>> neighbors;
    std::vector<std::vector<int>> conjugated;
    std::vector<std::vector<int>> order;
    std::vector<bool> atom_in_ring;
    std::vector<std::vector<bool>> bond_in_ring;
};

struct SmilesParse {
    std::vector<int> atoms;
    std::vector<std::array<int,3>> bonds;
};

struct EnoneAssignment {
    int carbonyl = -1;
    int oxygen = -1;
    int alpha = -1;
    int beta = -1;
    std::vector<int> chromophore_atoms;
    std::vector<int> position_by_atom;
    int base_nm = 212;
};

int element_z_from_symbol(std::string symbol) {
    if (symbol.empty()) return 0;
    symbol[0] = static_cast<char>(std::toupper(symbol[0]));
    if (symbol.size() > 1) symbol[1] = static_cast<char>(std::tolower(symbol[1]));

    static const std::map<std::string,int> z = {
        {"H",1},{"B",5},{"C",6},{"N",7},{"O",8},{"F",9},{"P",15},{"S",16},
        {"Cl",17},{"Br",35},{"I",53},{"Zn",30}
    };
    std::map<std::string,int>::const_iterator it = z.find(symbol);
    return it == z.end() ? 0 : it->second;
}

int smi_bond_order(const int pending, const bool left_aromatic, const bool right_aromatic) {
    if (pending > 0) return pending;
    return (left_aromatic && right_aromatic) ? 2 : 1;
}

SmilesParse parse_smiles_bonds(const std::string& smiles) {
    SmilesParse parsed;
    std::vector<bool> aromatic_atoms;
    std::vector<int> branch_stack;
    std::map<std::string, std::pair<int,int>> rings;
    int previous = -1;
    int pending_order = 0;

    for (size_t i = 0; i < smiles.size(); ++i) {
        char ch = smiles[i];
        if (ch == '(') {
            branch_stack.push_back(previous);
            continue;
        }
        if (ch == ')') {
            if (!branch_stack.empty()) {
                previous = branch_stack.back();
                branch_stack.pop_back();
            }
            continue;
        }
        if (ch == '.') {
            previous = -1;
            pending_order = 0;
            continue;
        }
        if (ch == '-' || ch == '/' || ch == '\\') {
            pending_order = 1;
            continue;
        }
        if (ch == '=' || ch == ':') {
            pending_order = 2;
            continue;
        }
        if (ch == '#') {
            pending_order = 3;
            continue;
        }

        std::string ring_label;
        if (ch == '%' && i + 2 < smiles.size() && std::isdigit(smiles[i + 1]) && std::isdigit(smiles[i + 2])) {
            ring_label = smiles.substr(i + 1, 2);
            i += 2;
        } else if (std::isdigit(ch)) {
            ring_label = std::string(1, ch);
        }
        if (!ring_label.empty()) {
            std::map<std::string, std::pair<int,int>>::iterator it = rings.find(ring_label);
            if (it == rings.end()) {
                rings[ring_label] = std::make_pair(previous, pending_order);
            } else if (previous >= 0 && it->second.first >= 0) {
                int stored = it->second.second;
                int order = pending_order > 0 ? pending_order : stored;
                parsed.bonds.push_back({{it->second.first, previous, smi_bond_order(order, aromatic_atoms[it->second.first], aromatic_atoms[previous])}});
                rings.erase(it);
            }
            pending_order = 0;
            continue;
        }

        std::string symbol;
        bool aromatic = false;
        if (ch == '[') {
            ++i;
            while (i < smiles.size() && std::isdigit(smiles[i])) ++i;
            if (i < smiles.size() && std::isalpha(smiles[i])) {
                symbol = std::string(1, smiles[i]);
                aromatic = std::islower(smiles[i]);
                if (i + 1 < smiles.size() && std::islower(smiles[i + 1]) && std::isupper(smiles[i])) {
                    symbol.push_back(smiles[++i]);
                }
            }
            while (i < smiles.size() && smiles[i] != ']') ++i;
        } else if (std::isalpha(ch)) {
            symbol = std::string(1, ch);
            aromatic = std::islower(ch);
            if (i + 1 < smiles.size() && std::islower(smiles[i + 1]) && std::isupper(ch)) {
                symbol.push_back(smiles[++i]);
            }
        }

        int z = element_z_from_symbol(symbol);
        if (z == 0) {
            pending_order = 0;
            continue;
        }

        int current = static_cast<int>(parsed.atoms.size());
        parsed.atoms.push_back(z);
        aromatic_atoms.push_back(aromatic);
        if (previous >= 0) {
            parsed.bonds.push_back({{previous, current, smi_bond_order(pending_order, aromatic_atoms[previous], aromatic_atoms[current])}});
        }
        previous = current;
        pending_order = 0;
    }
    return parsed;
}

bool path_exists_without_edge(const std::vector<std::vector<int>>& neighbors, const int start, const int goal, const int skip_a, const int skip_b) {
    std::vector<bool> seen(neighbors.size(), false);
    std::vector<int> stack(1, start);
    while (!stack.empty()) {
        int cur = stack.back();
        stack.pop_back();
        if (cur == goal) return true;
        if (seen[cur]) continue;
        seen[cur] = true;
        for (int next : neighbors[cur]) {
            if ((cur == skip_a && next == skip_b) || (cur == skip_b && next == skip_a)) continue;
            if (!seen[next]) stack.push_back(next);
        }
    }
    return false;
}

BondGraph build_bond_graph(const size_t& num_atoms, const std::vector<std::pair<int,int>>& pi_electrons_ordered, const std::vector<std::array<int,5>>& bond_info, const std::string& smiles) {
    BondGraph g;
    g.atomic_num.assign(num_atoms, 0);
    g.pi_electrons.assign(num_atoms, 0);
    g.neighbors.assign(num_atoms, std::vector<int>());
    g.conjugated.assign(num_atoms, std::vector<int>(num_atoms, 0));
    g.order.assign(num_atoms, std::vector<int>(num_atoms, 0));
    g.atom_in_ring.assign(num_atoms, false);
    g.bond_in_ring.assign(num_atoms, std::vector<bool>(num_atoms, false));

    for (const auto& pair : pi_electrons_ordered) {
        if (pair.first >= 0 && static_cast<size_t>(pair.first) < num_atoms) g.pi_electrons[pair.first] = pair.second;
    }

    for (const auto& bond : bond_info) {
        int a = bond[0], b = bond[1];
        if (a < 0 || b < 0 || static_cast<size_t>(a) >= num_atoms || static_cast<size_t>(b) >= num_atoms) continue;
        g.atomic_num[a] = bond[2];
        g.atomic_num[b] = bond[3];
        g.neighbors[a].push_back(b);
        g.neighbors[b].push_back(a);
        g.conjugated[a][b] = bond[4];
        g.conjugated[b][a] = bond[4];
    }

    SmilesParse parsed = parse_smiles_bonds(smiles);
    for (size_t i = 0; i < parsed.atoms.size() && i < num_atoms; ++i) {
        if (g.atomic_num[i] == 0) g.atomic_num[i] = parsed.atoms[i];
    }
    for (const auto& bond : parsed.bonds) {
        int a = bond[0], b = bond[1], order = bond[2];
        if (a < 0 || b < 0 || static_cast<size_t>(a) >= num_atoms || static_cast<size_t>(b) >= num_atoms) continue;
        g.order[a][b] = order;
        g.order[b][a] = order;
    }

    for (size_t a = 0; a < g.neighbors.size(); ++a) {
        for (int b : g.neighbors[a]) {
            if (static_cast<int>(a) >= b) continue;
            if (path_exists_without_edge(g.neighbors, static_cast<int>(a), b, static_cast<int>(a), b)) {
                g.atom_in_ring[a] = true;
                g.atom_in_ring[b] = true;
                g.bond_in_ring[a][b] = true;
                g.bond_in_ring[b][a] = true;
            }
        }
    }

    return g;
}

bool is_double_like(const BondGraph& g, const int a, const int b) {
    return g.order[a][b] == 2;
}

int count_h_neighbors(const BondGraph& g, const int atom) {
    int count = 0;
    for (int n : g.neighbors[atom]) if (g.atomic_num[n] == 1) ++count;
    return count;
}

bool has_carbonyl_oxygen(const BondGraph& g, const int carbon, int* oxygen_out = NULL) {
    if (g.atomic_num[carbon] != 6) return false;
    for (int n : g.neighbors[carbon]) {
        if (g.atomic_num[n] == 8 && is_double_like(g, carbon, n)) {
            if (oxygen_out) *oxygen_out = n;
            return true;
        }
    }
    return false;
}

int woodward_base_value(const BondGraph& g, const int carbonyl, const int oxygen) {
    bool single_o = false;
    bool ester = false;
    bool acid = false;
    for (int n : g.neighbors[carbonyl]) {
        if (n == oxygen || g.atomic_num[n] != 8 || is_double_like(g, carbonyl, n)) continue;
        single_o = true;
        for (int nn : g.neighbors[n]) {
            if (nn == carbonyl) continue;
            if (g.atomic_num[nn] == 1) acid = true;
            if (g.atomic_num[nn] == 6) ester = true;
        }
    }
    if (acid) return 196;
    if (ester || single_o) return 195;
    if (count_h_neighbors(g, carbonyl) > 0) return 218;
    return 212;
}

int count_component_double_bonds(const BondGraph& g, const std::vector<int>& component, const bool carbon_only) {
    std::set<int> atoms(component.begin(), component.end());
    int count = 0;
    for (int a : component) {
        for (int b : g.neighbors[a]) {
            if (a >= b || atoms.count(b) == 0 || !is_double_like(g, a, b)) continue;
            if (carbon_only && (g.atomic_num[a] != 6 || g.atomic_num[b] != 6)) continue;
            ++count;
        }
    }
    return count;
}

int count_exocyclic_double_bonds(const BondGraph& g, const std::vector<int>& component) {
    std::set<int> atoms(component.begin(), component.end());
    int count = 0;
    for (int a : component) {
        for (int b : g.neighbors[a]) {
            if (a >= b || atoms.count(b) == 0 || g.atomic_num[a] != 6 || g.atomic_num[b] != 6 || !is_double_like(g, a, b)) continue;
            if (g.atom_in_ring[a] != g.atom_in_ring[b] || (g.atom_in_ring[a] && g.atom_in_ring[b] && !g.bond_in_ring[a][b])) ++count;
        }
    }
    return count;
}

int count_endocyclic_double_bonds(const BondGraph& g, const std::vector<int>& component) {
    std::set<int> atoms(component.begin(), component.end());
    int count = 0;
    for (int a : component) {
        for (int b : g.neighbors[a]) {
            if (a >= b || atoms.count(b) == 0 || g.atomic_num[a] != 6 || g.atomic_num[b] != 6 || !is_double_like(g, a, b)) continue;
            if (g.bond_in_ring[a][b]) ++count;
        }
    }
    return count;
}

std::vector<int> conjugated_component(const BondGraph& g, const int start) {
    std::vector<int> result;
    std::vector<bool> seen(g.atomic_num.size(), false);
    std::vector<int> stack(1, start);
    while (!stack.empty()) {
        int cur = stack.back();
        stack.pop_back();
        if (seen[cur]) continue;
        seen[cur] = true;
        result.push_back(cur);
        for (int next : g.neighbors[cur]) {
            if (g.conjugated[cur][next] == 2 && !seen[next]) stack.push_back(next);
        }
    }
    return result;
}

std::vector<EnoneAssignment> find_enones(const BondGraph& g) {
    std::vector<EnoneAssignment> enones;
    std::set<std::pair<int,int>> seen;
    for (size_t c = 0; c < g.atomic_num.size(); ++c) {
        int oxygen = -1;
        if (!has_carbonyl_oxygen(g, static_cast<int>(c), &oxygen)) continue;
        for (int alpha : g.neighbors[c]) {
            if (g.atomic_num[alpha] != 6 || alpha == oxygen || g.conjugated[c][alpha] != 2) continue;
            for (int beta : g.neighbors[alpha]) {
                if (beta == static_cast<int>(c) || g.atomic_num[beta] != 6 || !is_double_like(g, alpha, beta)) continue;
                if (seen.count(std::make_pair(static_cast<int>(c), beta))) continue;
                seen.insert(std::make_pair(static_cast<int>(c), beta));
                EnoneAssignment item;
                item.carbonyl = static_cast<int>(c);
                item.oxygen = oxygen;
                item.alpha = alpha;
                item.beta = beta;
                item.chromophore_atoms = conjugated_component(g, static_cast<int>(c));
                item.position_by_atom.assign(g.atomic_num.size(), -1);
                item.position_by_atom[c] = 0;
                std::vector<int> stack(1, alpha);
                item.position_by_atom[alpha] = 1;
                while (!stack.empty()) {
                    int cur = stack.back();
                    stack.pop_back();
                    int next_pos = std::min(4, item.position_by_atom[cur] + 1);
                    for (int next : g.neighbors[cur]) {
                        if (next == static_cast<int>(c) || g.atomic_num[next] != 6 || g.conjugated[cur][next] != 2) continue;
                        if (item.position_by_atom[next] == -1 || next_pos < item.position_by_atom[next]) {
                            item.position_by_atom[next] = next_pos;
                            stack.push_back(next);
                        }
                    }
                }
                item.base_nm = woodward_base_value(g, static_cast<int>(c), oxygen);
                enones.push_back(item);
            }
        }
    }
    return enones;
}

int woodward_sub_value(const std::string& kind, const int pos) {
    static const std::map<std::string, std::array<int,4>> values = {
        {"alkoxy",  {{29, 22, 17, 31}}},
        {"hydroxy", {{38, 14, 50,  0}}},
        {"alkyl",   {{11, 19, 18, 18}}},
        {"bromo",   {{38, 33,  0,  0}}},
        {"chloro",  {{28, 22, 12, 12}}},
        {"carboxy", {{12, 12, 12, 12}}}
    };
    std::map<std::string, std::array<int,4>>::const_iterator it = values.find(kind);
    if (it == values.end() || pos < 1) return 0;
    return it->second[std::min(4, pos) - 1];
}

int fieser_sub_value(const std::string& kind) {
    static const std::map<std::string,int> values = {
        {"alkoxy", 6}, {"alkyl", 5}, {"carboxy", 0}, {"dialkylamine", 60},
        {"halogen", 10}, {"phenoxy", 18}, {"thioether", 30}
    };
    std::map<std::string,int>::const_iterator it = values.find(kind);
    return it == values.end() ? 0 : it->second;
}

std::string substituent_kind(const BondGraph& g, const int anchor, const int sub_atom, const std::set<int>& chromophore_atoms) {
    if (chromophore_atoms.count(sub_atom)) return "";
    int z = g.atomic_num[sub_atom];
    if (z == 17) return "chloro";
    if (z == 35) return "bromo";
    if (z == 9 || z == 53) return "halogen";
    if (z == 16) return "thioether";
    if (z == 7) {
        int carbon_neighbors = 0;
        for (int n : g.neighbors[sub_atom]) if (n != anchor && g.atomic_num[n] == 6) ++carbon_neighbors;
        return carbon_neighbors >= 1 ? "dialkylamine" : "";
    }
    if (z == 8) {
        bool attached_carbon = false;
        bool attached_h = false;
        bool attached_aromatic = false;
        for (int n : g.neighbors[sub_atom]) {
            if (n == anchor) continue;
            attached_carbon = attached_carbon || g.atomic_num[n] == 6;
            attached_h = attached_h || g.atomic_num[n] == 1;
            attached_aromatic = attached_aromatic || (g.atomic_num[n] == 6 && g.pi_electrons[n] > 0);
        }
        if (attached_aromatic) return "phenoxy";
        if (attached_carbon) return "alkoxy";
        if (attached_h || g.neighbors[sub_atom].size() == 1) return "hydroxy";
    }
    if (z == 6) {
        int oxygen = -1;
        if (has_carbonyl_oxygen(g, sub_atom, &oxygen)) return "carboxy";
        if (g.pi_electrons[sub_atom] == 0 && g.conjugated[anchor][sub_atom] != 2) return "alkyl";
    }
    return "";
}

int woodward_substituent_increment(const BondGraph& g, const EnoneAssignment& enone) {
    std::set<int> chromophore_atoms(enone.chromophore_atoms.begin(), enone.chromophore_atoms.end());
    std::array<int,4> best = {{0,0,0,0}};
    for (size_t atom = 0; atom < enone.position_by_atom.size(); ++atom) {
        int pos = enone.position_by_atom[atom];
        if (pos <= 0) continue;
        for (int sub : g.neighbors[atom]) {
            std::string kind = substituent_kind(g, static_cast<int>(atom), sub, chromophore_atoms);
            if (kind == "halogen") kind = "chloro";
            best[std::min(4, pos) - 1] = std::max(best[std::min(4, pos) - 1], woodward_sub_value(kind, pos));
        }
    }
    return best[0] + best[1] + best[2] + best[3];
}

int fieser_substituent_increment(const BondGraph& g, const std::vector<int>& component) {
    std::set<int> chromophore_atoms(component.begin(), component.end());
    int total = 0;
    std::set<std::pair<int,int>> counted;
    for (int atom : component) {
        if (g.atomic_num[atom] != 6) continue;
        for (int sub : g.neighbors[atom]) {
            std::string kind = substituent_kind(g, atom, sub, chromophore_atoms);
            if (kind == "chloro" || kind == "bromo") kind = "halogen";
            if (kind.empty()) continue;
            std::pair<int,int> key = std::make_pair(std::min(atom, sub), std::max(atom, sub));
            if (counted.insert(key).second) total += fieser_sub_value(kind);
        }
    }
    return total;
}

int count_alkyl_substituents(const BondGraph& g, const std::vector<int>& component) {
    std::set<int> chromophore_atoms(component.begin(), component.end());
    int count = 0;
    std::set<std::pair<int,int>> counted;
    for (int atom : component) {
        for (int sub : g.neighbors[atom]) {
            if (substituent_kind(g, atom, sub, chromophore_atoms) != "alkyl") continue;
            std::pair<int,int> key = std::make_pair(std::min(atom, sub), std::max(atom, sub));
            if (counted.insert(key).second) ++count;
        }
    }
    return count;
}

bool component_contains_carbonyl(const BondGraph& g, const std::vector<int>& component) {
    for (int atom : component) if (has_carbonyl_oxygen(g, atom)) return true;
    return false;
}

std::vector<int> estimate_rule_based_lambdas(const BondGraph& g, const std::vector<std::vector<int>>& chromophores) {
    std::vector<int> lambdas;
    std::vector<EnoneAssignment> enones = find_enones(g);
    for (const EnoneAssignment& enone : enones) {
        int double_count = count_component_double_bonds(g, enone.chromophore_atoms, false);
        int total = enone.base_nm;
        total += std::max(0, double_count - 2) * 30;
        total += woodward_substituent_increment(g, enone);
        total += 5 * count_exocyclic_double_bonds(g, enone.chromophore_atoms);
        lambdas.push_back(total);
    }

    for (const std::vector<int>& component : chromophores) {
        bool useful = false;
        for (int atom : component) {
            if (g.pi_electrons[atom] > 0 || g.atomic_num[atom] == 6 || g.atomic_num[atom] == 8) {
                useful = true;
                break;
            }
        }
        if (!useful || component_contains_carbonyl(g, component)) continue;

        int double_count = count_component_double_bonds(g, component, true);
        if (double_count >= 4) {
            int m = count_alkyl_substituents(g, component);
            int endo = count_endocyclic_double_bonds(g, component);
            int exo = count_exocyclic_double_bonds(g, component);
            double total = 114.0 + 5.0 * m + double_count * (48.0 - 1.7 * double_count) - 16.5 * endo - 10.0 * exo;
            lambdas.push_back(static_cast<int>(std::round(total)));
        } else if (double_count >= 2) {
            bool cyclic = false;
            for (int atom : component) cyclic = cyclic || g.atom_in_ring[atom];
            int base = cyclic ? 214 : 217;
            int total = base + std::max(0, double_count - 2) * 30;
            total += fieser_substituent_increment(g, component);
            total += 5 * count_exocyclic_double_bonds(g, component);
            lambdas.push_back(total);
        }
    }

    std::sort(lambdas.begin(), lambdas.end());
    lambdas.erase(std::unique(lambdas.begin(), lambdas.end()), lambdas.end());
    return lambdas;
}

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

    BondGraph molecular_graph = build_bond_graph(num_atoms, pi_electrons_ordered, bond_info, SMILES);
    std::vector<int> empirical_lambdas = estimate_rule_based_lambdas(molecular_graph, chromophores);

    //Step 3: Modified spectra calculation and superposition
    cache = empirical_lambdas;
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

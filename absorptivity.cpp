#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <unordered_map>
#include <set>
#include <string>
#include <dlfcn.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> //crucial: makes python lists and c++ vectors compatible

/*
1. build library in c++ that can be used in python (use extern C modifiers and modify CMake settings; see below)
2. filter dataset with polars in python
3. utilize c++ library of computational commands to perform analysis (and output things like spectrography plots?)
4. add an argparse (help and info) so that you can use terminal commands with the project!-------e.g. ```UV_plot -name "Jacaranda mimosifolia" -verbose```

For plants, add in (a properly weighted) chlorophyll a and b absorption spectrum. c.f. https://www.sciencedirect.com/science/article/pii/S0034425708000813
(chlorophyll a/b reference curves are staged in empirical_spectra/ by fetch_empirical_spectra.sh; not yet wired into the chromophore pipeline below -- that's a later step)

Structural data (bond order, aromaticity/ring membership, formal charge, H counts) is now forwarded directly from
RDKit via main.py instead of being reconstructed by re-parsing the SMILES string here. That reconstruction (and the
O(E*(V+E)) per-edge ring-detection DFS it required) has been removed; see build_bond_graph().
*/

// ---------------------------------------------------------------------------
// Sparse per-edge storage: an organic molecule's bond graph has O(n) edges,
// so we key edge attributes (conjugation status / bond order / ring membership)
// by an order-independent pair hash instead of allocating O(n^2) dense matrices.
// ---------------------------------------------------------------------------
inline long long edge_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<long long>(a) << 32) | static_cast<unsigned int>(b);
}
inline int edge_lookup(const std::unordered_map<long long,int>& m, int a, int b) {
    std::unordered_map<long long,int>::const_iterator it = m.find(edge_key(a, b));
    return it == m.end() ? 0 : it->second;
}

// groupings(): connected components restricted to edges whose status == desired_val.
// O(n + m) BFS over the adjacency list (previously O(n^2): a dense adjacency-matrix scan per atom).
std::vector<std::vector<int>> groupings(const std::vector<std::vector<int>>& neighbors,
                                         const std::unordered_map<long long,int>& edge_status,
                                         const int num_atoms, const int desired_val) {
    std::vector<std::vector<int>> groups;
    std::vector<bool> visited(num_atoms, false);

    for (int i = 0; i < num_atoms; i++) {
        if (visited[i]) continue;
        std::vector<int> current_group;
        std::vector<int> stack = {i};
        while (!stack.empty()) {
            int atom_idx = stack.back();
            stack.pop_back();
            if (visited[atom_idx]) continue;
            visited[atom_idx] = true;
            current_group.push_back(atom_idx);
            for (int j : neighbors[atom_idx]) {
                if (!visited[j] && edge_lookup(edge_status, atom_idx, j) == desired_val) stack.push_back(j);
            }
        }
        groups.push_back(std::move(current_group));
    }
    return groups;
}

//only accepting analysis in the wavelength range where individual atom absorptivities are negligible (and other data may be flawed anyway in the vacuum UV range) : 120-800nm (far UV - very near IR)
std::vector<std::pair<int,double>> cache; //stores memory from last input: (wavenumber in cm^-1, amplitude)

struct BondGraph {
    std::vector<int> atomic_num;
    std::vector<int> pi_electrons;
    std::vector<int> total_h;
    std::vector<int> formal_charge;
    std::vector<std::vector<int>> neighbors;       //adjacency list, O(n+m)
    std::unordered_map<long long,int> conjugated;  //edge -> bond status (0 none / 1 bond / 2 conjugated)
    std::unordered_map<long long,int> order;       //edge -> explicit bond order (1/2/3), Kekulized on the python side
    std::unordered_map<long long,int> ring_bond;   //edge -> 1 if the bond is in a ring
    std::vector<bool> atom_in_ring;
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

// A single detected chromophore, its rule-based lambda_max, and enough context to either
// look up a matching empirical template or (for enones) split it when none is found.
struct ChromophoreEstimate {
    std::string base_type;
    int lambda_nm = 0;
    int double_bond_count = 0;
    std::vector<int> atoms;
    bool is_satellite = false; //secondary/extension piece from a split -> reduced amplitude weight
    bool is_enone = false;
    EnoneAssignment enone; //valid only when is_enone
};

// ---------------------------------------------------------------------------
// Graph construction: everything below (bond order, ring membership, formal charge,
// total H count) is now forwarded from RDKit via main.py -- see analyze_conjugation().
// This replaced a SMILES-string re-parser plus an O(E*(V+E)) per-edge "remove edge,
// check connectivity" ring detector (path_exists_without_edge), both now unnecessary.
// ---------------------------------------------------------------------------
BondGraph build_bond_graph(const size_t& num_atoms,
                            const std::vector<std::pair<int,int>>& pi_electrons_ordered,
                            const std::vector<std::array<int,7>>& bond_info,
                            const std::vector<std::array<int,5>>& atom_info) {
    BondGraph g;
    g.atomic_num.assign(num_atoms, 0);
    g.pi_electrons.assign(num_atoms, 0);
    g.total_h.assign(num_atoms, 0);
    g.formal_charge.assign(num_atoms, 0);
    g.neighbors.assign(num_atoms, std::vector<int>());
    g.atom_in_ring.assign(num_atoms, false);

    // atom_info: [idx, atomic_num, total_H, formal_charge, in_ring] -- covers every atom,
    // including unbonded single-atom species (e.g. a bare ion like [Zn]) that have no entry in bond_info at all.
    for (const auto& a : atom_info) {
        int idx = a[0];
        if (idx < 0 || static_cast<size_t>(idx) >= num_atoms) continue;
        g.atomic_num[idx] = a[1];
        g.total_h[idx] = a[2];
        g.formal_charge[idx] = a[3];
        g.atom_in_ring[idx] = a[4] != 0;
    }
    for (const auto& p : pi_electrons_ordered) {
        if (p.first >= 0 && static_cast<size_t>(p.first) < num_atoms) g.pi_electrons[p.first] = p.second;
    }
    // bond_info: [idx1, idx2, Z1, Z2, bond_status, bond_order, bond_in_ring?]
    for (const auto& bond : bond_info) {
        int a = bond[0], b = bond[1];
        if (a < 0 || b < 0 || static_cast<size_t>(a) >= num_atoms || static_cast<size_t>(b) >= num_atoms) continue;
        g.neighbors[a].push_back(b);
        g.neighbors[b].push_back(a);
        long long key = edge_key(a, b);
        g.conjugated[key] = bond[4];
        g.order[key] = bond[5];
        g.ring_bond[key] = bond[6];
    }
    return g;
}

bool is_double_like(const BondGraph& g, const int a, const int b) {
    return edge_lookup(g.order, a, b) == 2;
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
    if (g.total_h[carbonyl] > 0) return 218; //O(1) via forwarded H count, was an O(degree) neighbor scan
    return 212;
}

std::string enone_base_type(int base_nm) {
    switch (base_nm) {
        case 196: return "acid_enone";
        case 195: return "ester_enone";
        case 218: return "aldehyde_enone";
        default:  return "ketone_enone"; //212 and any other fallback
    }
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
            bool bond_in_ring = edge_lookup(g.ring_bond, a, b) == 1;
            if (g.atom_in_ring[a] != g.atom_in_ring[b] || (g.atom_in_ring[a] && g.atom_in_ring[b] && !bond_in_ring)) ++count;
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
            if (edge_lookup(g.ring_bond, a, b) == 1) ++count;
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
            if (edge_lookup(g.conjugated, cur, next) == 2 && !seen[next]) stack.push_back(next);
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
            if (g.atomic_num[alpha] != 6 || alpha == oxygen || edge_lookup(g.conjugated, static_cast<int>(c), alpha) != 2) continue;
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
                        if (next == static_cast<int>(c) || g.atomic_num[next] != 6 || edge_lookup(g.conjugated, cur, next) != 2) continue;
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
    static const std::map<std::string, std::array<int,4>> values = { // refined values, doi: 10.1021/jo01164a003 (via ChromoPredict, github.com/CompPhotoChem/ChromoPredict)
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
    static const std::map<std::string,int> values = { // doi: 10.1021/jo01164a003 (via ChromoPredict)
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
        bool attached_aromatic = false;
        for (int n : g.neighbors[sub_atom]) {
            if (n == anchor) continue;
            attached_carbon = attached_carbon || g.atomic_num[n] == 6;
            attached_aromatic = attached_aromatic || (g.atomic_num[n] == 6 && g.pi_electrons[n] > 0);
        }
        if (attached_aromatic) return "phenoxy";
        if (attached_carbon) return "alkoxy";
        if (g.total_h[sub_atom] > 0 || g.neighbors[sub_atom].size() == 1) return "hydroxy"; //O(1) via forwarded H count
    }
    if (z == 6) {
        int oxygen = -1;
        if (has_carbonyl_oxygen(g, sub_atom, &oxygen)) return "carboxy";
        if (g.pi_electrons[sub_atom] == 0 && edge_lookup(g.conjugated, anchor, sub_atom) != 2) return "alkyl";
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

std::vector<ChromophoreEstimate> estimate_rule_based_lambdas(const BondGraph& g, const std::vector<std::vector<int>>& chromophores) {
    std::vector<ChromophoreEstimate> estimates;
    for (const EnoneAssignment& enone : find_enones(g)) {
        int double_count = count_component_double_bonds(g, enone.chromophore_atoms, false);
        int total = enone.base_nm;
        total += std::max(0, double_count - 2) * 30;
        total += woodward_substituent_increment(g, enone);
        total += 5 * count_exocyclic_double_bonds(g, enone.chromophore_atoms);

        ChromophoreEstimate est;
        est.base_type = enone_base_type(enone.base_nm);
        est.lambda_nm = total;
        est.double_bond_count = double_count;
        est.atoms = enone.chromophore_atoms;
        est.is_enone = true;
        est.enone = enone;
        estimates.push_back(std::move(est));
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

            ChromophoreEstimate est;
            est.base_type = "extended_polyene";
            est.lambda_nm = static_cast<int>(std::round(total));
            est.double_bond_count = double_count;
            est.atoms = component;
            estimates.push_back(std::move(est));
        } else if (double_count >= 2) {
            bool cyclic = false;
            for (int atom : component) cyclic = cyclic || g.atom_in_ring[atom];
            int base = cyclic ? 214 : 217;
            int total = base + std::max(0, double_count - 2) * 30;
            total += fieser_substituent_increment(g, component);
            total += 5 * count_exocyclic_double_bonds(g, component);

            ChromophoreEstimate est;
            est.base_type = cyclic ? "diene_cyclic" : "diene_acyclic";
            est.lambda_nm = total;
            est.double_bond_count = double_count;
            est.atoms = component;
            estimates.push_back(std::move(est));
        }
    }
    return estimates;
}

// ---------------------------------------------------------------------------
// Empirical spectral template library: small digitized/parameterized (wavelength_nm, epsilon)
// reference curves for each base chromophore class, loaded once from empirical_spectra/
// (populated by fetch_empirical_spectra.sh) and cached for the lifetime of the process.
// ---------------------------------------------------------------------------
std::string module_directory() {
    Dl_info info;
    if (dladdr(reinterpret_cast<void*>(&module_directory), &info) && info.dli_fname) {
        std::string path(info.dli_fname);
        size_t slash = path.find_last_of('/');
        return slash == std::string::npos ? "." : path.substr(0, slash);
    }
    return ".";
}

std::map<std::string, std::vector<std::pair<double,double>>> load_empirical_library() {
    std::map<std::string, std::vector<std::pair<double,double>>> library;
    static const char* names[] = {"ketone_enone", "aldehyde_enone", "diene_acyclic", "diene_cyclic", "extended_polyene"};
    std::string dir = module_directory() + "/empirical_spectra/";
    for (const char* name : names) {
        std::ifstream file(dir + name + ".csv");
        if (!file.is_open()) continue;
        std::vector<std::pair<double,double>> points;
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream iss(line);
            double nm, eps;
            if (iss >> nm >> eps) points.push_back(std::make_pair(nm, eps));
        }
        if (!points.empty()) library[name] = std::move(points);
    }
    // No separate heteroanular/homoanular-ring reference was sourced; the acyclic diene
    // lineshape (re-centered/rescaled per fragment below) stands in for the cyclic case too.
    if (!library.count("diene_cyclic") && library.count("diene_acyclic")) {
        library["diene_cyclic"] = library["diene_acyclic"];
    }
    return library;
}

const std::map<std::string, std::vector<std::pair<double,double>>>& empirical_library() {
    static std::map<std::string, std::vector<std::pair<double,double>>> library = load_empirical_library();
    return library;
}

bool has_empirical_spectrum(const std::string& base_type) {
    const auto& lib = empirical_library();
    std::map<std::string, std::vector<std::pair<double,double>>>::const_iterator it = lib.find(base_type);
    return it != lib.end() && !it->second.empty();
}

// Representative peak amplitude for a fragment. Polyenes/dienes use the standard Fieser-Kuhn
// intensity rule (epsilon_max ~= 1.74e4 * n conjugated C=C bonds -- the amplitude counterpart of
// the lambda_max formula already used above; see e.g. Pretsch/Buhlmann/Badertscher, "Structure
// Determination of Organic Compounds"). Enones don't have as clean a closed-form amplitude rule,
// so a representative order-of-magnitude K-band intensity is used instead -- documented as
// approximate, not a per-molecule measurement.
double representative_amplitude(const ChromophoreEstimate& est) {
    double value;
    if (est.base_type == "extended_polyene" || est.base_type.rfind("diene", 0) == 0) {
        value = 1.74e4 * std::max(1, est.double_bond_count);
    } else {
        value = 1.0e4 + 1500.0 * std::max(0, est.double_bond_count - 2);
    }
    return est.is_satellite ? value * 0.35 : value; //split-off extension pieces read as a secondary/shoulder feature, not a second full band
}

// Shift a template curve in WAVELENGTH space (nm), not wavenumber space: Woodward-Fieser/Fieser-Kuhn
// increments are additive in nm, and since wavenumber = 1e7/lambda, a uniform nm shift is NOT a uniform
// cm^-1 shift (it varies roughly as 1/lambda^2). Converting to wavenumber only after shifting keeps every
// line consistent with what the nm-based increment actually means.
std::vector<std::pair<int,double>> build_modified_spectrum(const ChromophoreEstimate& est) {
    const std::vector<std::pair<double,double>>& templ = empirical_library().at(est.base_type);
    double template_peak_nm = templ.front().first, template_peak_eps = templ.front().second;
    for (const auto& pt : templ) {
        if (pt.second > template_peak_eps) { template_peak_eps = pt.second; template_peak_nm = pt.first; }
    }
    double delta_nm = est.lambda_nm - template_peak_nm;
    double target_eps = representative_amplitude(est);
    double scale = template_peak_eps > 0.0 ? target_eps / template_peak_eps : 1.0;

    std::vector<std::pair<int,double>> out;
    out.reserve(templ.size());
    for (const auto& pt : templ) {
        double shifted_nm = pt.first + delta_nm;
        if (shifted_nm < 120.0 || shifted_nm > 800.0) continue; //stay within the analysis window
        int wavenumber = static_cast<int>(std::lround(1.0e7 / shifted_nm)); //rounded wavenumber, cm^-1
        out.push_back(std::make_pair(wavenumber, pt.second * scale));
    }
    return out;
}

// Split an enone lacking a direct template into its minimal Woodward core (carbonyl + O + Cα=Cβ,
// keeping the FULL substituent-corrected lambda_max -- the substituents live on the core) plus a
// secondary "extension" fragment covering any additional conjugated double bonds beyond that core,
// classified and matched like an ordinary diene/polyene. Only meaningful when double_count > 2, and
// only ever applied once (the recursive call below passes a null EnoneAssignment), which keeps the
// process bounded: Woodward-Fieser rules only define a single flat "+30 nm per extra double bond"
// extension term to begin with, so there's no deeper structure to recurse into.
std::vector<ChromophoreEstimate> split_enone(const BondGraph& g, const EnoneAssignment& enone, int total_lambda, int double_count) {
    std::vector<ChromophoreEstimate> parts;
    if (double_count <= 2) return parts; //already minimal; nothing left to peel off

    std::set<int> core_set = {enone.carbonyl, enone.oxygen, enone.alpha, enone.beta};
    std::vector<int> core_atoms, tail_atoms;
    for (int atom : enone.chromophore_atoms) (core_set.count(atom) ? core_atoms : tail_atoms).push_back(atom);
    if (tail_atoms.empty()) return parts;

    ChromophoreEstimate core;
    core.base_type = enone_base_type(enone.base_nm);
    core.lambda_nm = total_lambda; //keep the full, substituent-corrected value: this is where the real band sits
    core.double_bond_count = 2;
    core.atoms = core_atoms;
    parts.push_back(std::move(core));

    bool tail_cyclic = false;
    for (int atom : tail_atoms) tail_cyclic = tail_cyclic || g.atom_in_ring[atom];
    int tail_double_count = double_count - 2;
    ChromophoreEstimate tail;
    tail.base_type = (tail_double_count >= 4) ? "extended_polyene" : (tail_cyclic ? "diene_cyclic" : "diene_acyclic");
    int tail_base = tail_cyclic ? 214 : 217;
    tail.lambda_nm = tail_base + std::max(0, tail_double_count - 2) * 30; //unsubstituted extension estimate: substituents were already scored on the core
    tail.double_bond_count = tail_double_count;
    tail.atoms = tail_atoms;
    tail.is_satellite = true;
    parts.push_back(std::move(tail));
    return parts;
}

std::vector<std::pair<int,double>> resolve_fragment(
    const BondGraph& g,
    const ChromophoreEstimate& est,
    const EnoneAssignment* enone_ref
) {
    // Extended enones must be decomposed before template matching.
    if (enone_ref != NULL && est.double_bond_count > 2) {
        std::vector<ChromophoreEstimate> parts =
            split_enone(g, *enone_ref, est.lambda_nm, est.double_bond_count);

        if (!parts.empty()) {
            std::vector<std::pair<int,double>> out;

            for (const ChromophoreEstimate& part : parts) {
                std::vector<std::pair<int,double>> sub =
                    resolve_fragment(g, part, NULL);
                out.insert(out.end(), sub.begin(), sub.end());
            }

            return out;
        }
    }

    if (has_empirical_spectrum(est.base_type))
        return build_modified_spectrum(est);

    // No empirical template: retain the existing fallback.
    int wavenumber =
        static_cast<int>(std::lround(
            1.0e7 / std::max(1, est.lambda_nm)));

    return {{wavenumber, representative_amplitude(est)}};
}

std::vector<std::pair<int,double>> molecule_spectrum(const std::string& SMILES, const size_t& num_atoms,
                                                       const std::vector<std::pair<int,int>>& pi_electrons_ordered,
                                                       const std::vector<std::array<int,7>>& bond_info,
                                                       const std::vector<std::array<int,5>>& atom_info) { //returns (wavenumber cm^-1, amplitude) pairs
    cache.clear(); //clearing cache and memory for new decomposition
    cache.shrink_to_fit();

    //Step 1: Chromophore division
    //adjacency-list + sparse edge maps: O(n + m), replacing the old O(n^2) dense adjacency-matrix build
    BondGraph molecular_graph = build_bond_graph(num_atoms, pi_electrons_ordered, bond_info, atom_info);
    int n = static_cast<int>(num_atoms);

    //NOTE: atom indices are 0-indexed!!!!!!
    std::vector<std::vector<int>> chromophores = groupings(molecular_graph.neighbors, molecular_graph.conjugated, n, 2); //first we identify chromophores
    std::vector<std::vector<int>> remainder_groups = groupings(molecular_graph.neighbors, molecular_graph.conjugated, n, 1); // next we identify remainder groups between the chromophores
    (void)remainder_groups; //reserved for future auxochrome-context work; not yet consumed downstream
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

    std::vector<ChromophoreEstimate> estimates = estimate_rule_based_lambdas(molecular_graph, chromophores);

    std::cout << "=== ESTIMATES ===\n"; //DIAGNOSTIC
    for (const auto& est : estimates) {
        std::cout << est.base_type << " | lambda=" << est.lambda_nm << " nm | double_bonds=" << est.double_bond_count << " | atoms=" << est.atoms.size() << '\n';
    }
    //Step 3: Modified spectra calculation and superposition
    //Each chromophore is matched against empirical_spectra/ (splitting extended enones once if there's no
    //direct match -- see resolve_fragment/split_enone), shifted in nm-space to its rule-based lambda_max,
    //rescaled to a representative amplitude, and the resulting lines from every fragment are pooled.
    for (const ChromophoreEstimate& est : estimates) {
        std::vector<std::pair<int,double>> piece = resolve_fragment(molecular_graph, est, est.is_enone ? &est.enone : NULL);
        std::cout << "  -> spectral points: " << piece.size() << '\n'; //DIAGNOSTIC
        cache.insert(cache.end(), piece.begin(), piece.end());
    }
    std::sort(cache.begin(), cache.end());
    return cache;
}

// inputs in bond_info: atom order/index 1 & 2, element 1 & 2 atomic numbers, bond status, bond order, bond-in-ring
// inputs in atom_info: atom index, atomic number, total H count, formal charge, atom-in-ring
/* Tests:
BENZENE / C1=CC=CC=C1
WATER,Rhizome / O
TRANS-LINALOOL-OXIDE / C[C@]1(CC[C@H](O1)C(C)(C)O)C=C
ZINC,Root,CID_23994,[Zn]
*/


PYBIND11_MODULE(absorptivity, m) {
    m.doc() = "Absorptivity calculator (v0)";
    m.def("compute_spectrum", [](const std::string SMILES, const int& num_atoms,
                                  const std::vector<std::pair<int,int>> pi_electrons_ordered,
                                  const std::vector<std::array<int,7>> bond_info,
                                  const std::vector<std::array<int,5>> atom_info) {
        return molecule_spectrum(SMILES, static_cast<size_t>(num_atoms), pi_electrons_ordered, bond_info, atom_info);
    }, "Decomposes a molecule into chromophores and remainder groups, identifies auxochromes, matches (or "
       "recursively splits, on a miss) each chromophore against a small empirical UV/Vis template library, and "
       "superposes the shifted/rescaled templates into the molecule's overall spectrum for 120-800 nm as "
       "(wavenumber cm^-1, amplitude) pairs. Amplitudes are per-chromophore relative intensities, not a "
       "concentration/path-length-scaled absorbance.");
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
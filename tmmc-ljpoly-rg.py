# Combined Code for TMMC Polymer Simulation 

import numpy as np
import math
import os
import random
import itertools 
import json 
import matplotlib.pyplot as plt 
from types import SimpleNamespace # For sim_config object

# --- Global Constants & Configuration Placeholders ---
KB = 1.0 # Assuming reduced units for simulation
sim_config = None 

# Create input_config.json Programmatically ---
def create_input_config_json(filename="input_config.json"):
    """
    Defines the default configuration and writes it to a JSON file.
    This function replaces the %%writefile magic for broader compatibility.
    """
    config_content = {
        "conformation_library_file": "/content/my_polymer_library.dat",
        "initial_solvent_box_file": None, # or "/content/pre_equilibrated_solvent.dat"
        "output_dir": "/content/tmmc_output",
        "target_temperature": 1.0,
        "library_temperature": 1.0,
        "random_seed": 12345,
        "num_independent_runs": 1, # In main script, this would mean looping the whole process
        "num_monomers": 20,
        "k_harmonic": 3.8e7, # Ensure units are consistent with LJ params
        "r_eq_harmonic": 1.1225, # In units of sigma
        "eps_pp": 1.0,
        "sigma_pp": 1.0,
        "num_solvent_particles": 467,
        "solvent_density": 0.7, # Not directly used if num_solvent_particles and box_side_length are set
        "eps_ss": 1.0,
        "sigma_ss": 1.0,
        "box_side_length": 8.9, # In units of sigma
        "eps_ms": 1.0, # Polymer-solvent interaction strength
        "sigma_ms": 1.0,
        "cv_type": "Rg",
        "num_cv_bins": 49,
        "cv_min": 1.26, 
        "cv_max": 4.09, # Can be updated by Cell 3 from library
        "num_growth_stages": 31,
        "mc_equilibration_cycles": 10000, # Outer GEE blocks for equilibration
        "mc_production_cycles": 20000,    # Outer GEE blocks for production
        "solvent_sweeps_per_gee_attempt": 10, # Solvent sweeps before each GEE attempt
        "eta_update_frequency": 500,      # In terms of outer GEE blocks
        "eta_damping_factor": 0.2,
        "solvent_max_displacement": 0.15, # In units of sigma
        "num_processors": 1, # For potential future parallelization hints
        "plot_data_log_freq": 200, # Log plotting data every N elementary steps (sweeps)
        "verbose_print_freq": 1000 # Print verbose output every N elementary steps (sweeps)
    }
    try:
        with open(filename, 'w') as f:
            json.dump(config_content, f, indent=4)
        print(f"Cell 1 part: Successfully created configuration file: '{filename}'")
        return filename
    except IOError as e:
        print(f"Cell 1 part ERROR: Could not write configuration file '{filename}': {e}")
        return None

# --- Cell 2 Logic: Load Configuration ---
DEFAULT_CONFIG_PY = {
    "conformation_library_file": None, 
    "initial_solvent_box_file": None,
    "output_dir": "./tmmc_output_default",
    "target_temperature": 1.0, "library_temperature": 1.0, "random_seed": None,
    "num_independent_runs": 1, "num_monomers": 20,
    "k_harmonic": 3.8e7, "r_eq_harmonic": 1.1225,
    "eps_pp": 1.0, "sigma_pp": 1.0,
    "num_solvent_particles": 0, "solvent_density": 0.7, 
    "eps_ss": 1.0, "sigma_ss": 1.0, "box_side_length": 8.9,
    "eps_ms": 1.0, "sigma_ms": 1.0,
    "cv_type": "Rg", "num_cv_bins": 49, "cv_min": None, "cv_max": None,
    "num_growth_stages": 31,
    "mc_equilibration_cycles": 1000, "mc_production_cycles": 2000,
    "solvent_sweeps_per_gee_attempt": 10,
    "eta_update_frequency": 100, "eta_damping_factor": 0.2,
    "solvent_max_displacement": 0.15, "num_processors": 1,
    "plot_data_log_freq": 100, "verbose_print_freq": 500,
    "KB": 1.0 # Reduced Boltzmann constant
}

def load_simulation_config(config_file_path):
    config_data = DEFAULT_CONFIG_PY.copy()
    try:
        with open(config_file_path, 'r') as f:
            file_params = json.load(f)
        config_data.update(file_params)
        print(f"Cell 2 part: Successfully loaded configuration from: {config_file_path}")
    except FileNotFoundError:
        print(f"Cell 2 part WARNING: Config file '{config_file_path}' not found. Using defaults.")
    except json.JSONDecodeError as e:
        print(f"Cell 2 part WARNING: Error decoding JSON from '{config_file_path}': {e}. Using defaults.")
    
    config_namespace = SimpleNamespace(**config_data)
    if config_namespace.conformation_library_file is None:
        print("Cell 2 part ERROR: 'conformation_library_file' must be specified.")
        return None # Or raise error
        
    # Ensure output directory exists
    os.makedirs(config_namespace.output_dir, exist_ok=True)
    
    # Derived/fallback parameters for frequencies if total cycles are very small
    total_outer_cycles = config_namespace.mc_equilibration_cycles + config_namespace.mc_production_cycles
    sweeps_per_block = max(1, config_namespace.solvent_sweeps_per_gee_attempt)
    total_elementary_steps = total_outer_cycles * sweeps_per_block

    if not hasattr(config_namespace, 'plot_data_log_freq') or config_namespace.plot_data_log_freq == 0:
        config_namespace.plot_data_log_freq = max(1, total_elementary_steps // 200 if total_elementary_steps > 0 else 100)
    if not hasattr(config_namespace, 'verbose_print_freq') or config_namespace.verbose_print_freq == 0:
        config_namespace.verbose_print_freq = max(1, total_elementary_steps // 100 if total_elementary_steps > 0 else 500)
    if not hasattr(config_namespace, 'eta_damping_factor'):
        config_namespace.eta_damping_factor = 0.2
    if not hasattr(config_namespace, 'KB'): # Add KB if not in json
        config_namespace.KB = 1.0


    return config_namespace

# --- Cell 3 Logic: Library Processing ---
# Helper functions for Cell 3 (Energies & CV)
def calculate_distance_sq_c3(p1, p2): return np.sum((p1 - p2)**2)
def calculate_distance_c3(p1, p2): return np.sqrt(calculate_distance_sq_c3(p1, p2))
def calculate_bonded_energy_c3(coords, k_harmonic, r_eq_harmonic):
    E = 0.0; N = len(coords)
    if N < 2: return 0.0
    for i in range(N - 1): E += 0.5 * k_harmonic * (calculate_distance_c3(coords[i], coords[i+1]) - r_eq_harmonic)**2
    return E
def calculate_lj_energy_pair_c3(r_sq, epsilon, sigma):
    if r_sq < 1e-12: return 1e10
    sig_sq_over_r_sq = (sigma**2) / r_sq; sr6 = sig_sq_over_r_sq**3; sr12 = sr6**2
    return 4.0 * epsilon * (sr12 - sr6)
def calculate_non_bonded_energy_c3(coords, epsilon_pp, sigma_pp):
    E = 0.0; N = len(coords)
    if N < 3: return 0.0 # For j >= i+2
    for i in range(N):
        for j in range(i + 2, N): E += calculate_lj_energy_pair_c3(calculate_distance_sq_c3(coords[i], coords[j]), epsilon_pp, sigma_pp)
    return E
def calculate_total_internal_energy_c3(coords, config):
    return calculate_bonded_energy_c3(coords, config.k_harmonic, config.r_eq_harmonic) + \
           calculate_non_bonded_energy_c3(coords, config.eps_pp, config.sigma_pp)
def calculate_cv_rg_first_monomer_origin_c3(coords, cv_type="Rg"):
    coords_arr = np.array(coords); N = len(coords_arr)
    if N == 0: return 0.0
    if cv_type.lower() == "rg":
        origin = coords_arr[0]; rg_sq = 0.0
        for i in range(N): rg_sq += np.sum((coords_arr[i] - origin)**2)
        return np.sqrt(rg_sq / N) if N > 0 else 0.0
    elif cv_type.lower() == "end_to_end_distance":
        return calculate_distance_c3(coords_arr[0], coords_arr[-1]) if N >= 2 else 0.0
    raise ValueError(f"Unsupported CV type: {cv_type}")

def yield_polymer_conformations_new_format_c3(library_file_path, num_monomers_per_conf):
    conf_idx = 0; current_conf_coords = []; reading_conf = False; start_ln_for_conf = -1
    try:
        with open(library_file_path, 'r') as f:
            for current_ln, line_content in enumerate(f, 1):
                line = line_content.strip()
                if line.startswith("# Conformation"):
                    if reading_conf and len(current_conf_coords) < num_monomers_per_conf:
                        print(f"Warning: Incomplete conf (started L{start_ln_for_conf}) before new tag L{current_ln}. Skipping.")
                    current_conf_coords = []; reading_conf = True; start_ln_for_conf = current_ln
                    continue
                if reading_conf:
                    if not line: continue
                    parts = line.split()
                    if len(parts) == 3:
                        try:
                            current_conf_coords.append([float(p) for p in parts])
                            if len(current_conf_coords) == num_monomers_per_conf:
                                yield conf_idx, start_ln_for_conf, np.array(current_conf_coords)
                                conf_idx += 1; current_conf_coords = []; reading_conf = False; start_ln_for_conf = -1
                        except ValueError:
                            print(f"Warning: Parse error L{current_ln}: '{line_content.strip()}'. Skipping conf block (L{start_ln_for_conf}).")
                            current_conf_coords = []; reading_conf = False; start_ln_for_conf = -1
                    else:
                        print(f"Warning: Format error L{current_ln}: '{line_content.strip()}'. Skipping conf block (L{start_ln_for_conf}).")
                        current_conf_coords = []; reading_conf = False; start_ln_for_conf = -1
            if reading_conf and current_conf_coords and len(current_conf_coords) < num_monomers_per_conf:
                print(f"Warning: EOF. Last conf (L{start_ln_for_conf}) incomplete. Skipping.")
    except FileNotFoundError: print(f"ERROR: Library file '{library_file_path}' not found."); yield None
    except Exception as e: print(f"ERROR reading library: {e}"); yield None

def process_library_cell3(config_obj): # Renamed sim_config to config_obj
    print("\n--- Cell 3 part: Processing Polymer Library ---")
    if not config_obj: print("ERROR: sim_config not loaded. Cannot run Cell 3 part."); return False
    
    print(f"Using: Library='{config_obj.conformation_library_file}', CV='{config_obj.cv_type}'")
    processed_data = []; processed_count = 0; error_flag_c3 = False
    min_cv_obs = float('inf'); max_cv_obs = float('-inf')

    if not os.path.exists(config_obj.conformation_library_file):
        print(f"ERROR: Library file '{config_obj.conformation_library_file}' does not exist."); return False

    conf_gen = yield_polymer_conformations_new_format_c3(config_obj.conformation_library_file, config_obj.num_monomers)
    for item in conf_gen:
        if item is None: error_flag_c3 = True; break
        original_idx, start_ln, coords = item
        try:
            cv_val = calculate_cv_rg_first_monomer_origin_c3(coords, config_obj.cv_type)
            total_E = calculate_total_internal_energy_c3(coords, config_obj)
        except Exception as e: print(f"Err CV/Energy for conf idx {original_idx} (L{start_ln}): {e}"); continue
        processed_data.append({"original_index": original_idx, "start_line_num": start_ln, "cv": cv_val, "energy": total_E})
        processed_count += 1; min_cv_obs = min(min_cv_obs, cv_val); max_cv_obs = max(max_cv_obs, cv_val)
        if processed_count % 1000 == 0: print(f"  Cell 3: Processed {processed_count} configurations...")
    
    if error_flag_c3: print("Cell 3: Processing stopped due to file reading error."); return False
    print(f"Cell 3: Finished processing. Processed {processed_count} valid configurations.")

    if processed_count > 0:
        processed_data_sorted = sorted(processed_data, key=lambda x: x["cv"])
        print(f"Observed CV range: Min={min_cv_obs:.4f}, Max={max_cv_obs:.4f}")
        updated_cv = False
        if config_obj.cv_min is None or min_cv_obs < config_obj.cv_min: config_obj.cv_min = min_cv_obs; updated_cv=True
        if config_obj.cv_max is None or max_cv_obs > config_obj.cv_max: config_obj.cv_max = max_cv_obs; updated_cv=True
        if updated_cv: print(f"  Updated sim_config.cv_min/max to: {config_obj.cv_min:.4f}/{config_obj.cv_max:.4f}")
        if abs(config_obj.cv_max - config_obj.cv_min) < 1e-9: config_obj.cv_bin_width = 0.0
        elif config_obj.num_cv_bins > 0: config_obj.cv_bin_width = (config_obj.cv_max - config_obj.cv_min) / config_obj.num_cv_bins
        else: config_obj.cv_bin_width = 0.0
        if updated_cv or config_obj.cv_bin_width is None: print(f"  Updated sim_config.cv_bin_width to: {config_obj.cv_bin_width:.6f if config_obj.cv_bin_width is not None else 'N/A'}")

        out_fname = "library_processed_index.txt"
        out_path = os.path.join(config_obj.output_dir, out_fname)
        print(f"\nCell 3: Writing processed library index to: {out_path}")
        try:
            with open(out_path, 'w') as f_out:
                f_out.write(f"# Processed Library Data (Sorted by CV)\n# Source: {config_obj.conformation_library_file}\n")
                f_out.write(f"# CV: {config_obj.cv_type} (Rg origin: first monomer)\n")
                f_out.write(f"# Original_Index Start_Line_of_Comment(1-based) CV_Value Total_Energy(Vacuum_Internal)\n")
                for item_data in processed_data_sorted: f_out.write(f"{item_data['original_index']} {item_data['start_line_num']} {item_data['cv']:.6f} {item_data['energy']:.6f}\n")
            print(f"Cell 3: Successfully wrote index for {len(processed_data_sorted)} configurations.")
        except IOError as e: print(f"Cell 3 ERROR: Could not write output file '{out_path}': {e}"); return False
        return True
    else: print("Cell 3: No valid configurations processed. Output file not written."); return False

# --- Cell 4 Logic: TMMC-GEE Simulation ---
# (Helper functions for Cell 4 are defined globally now, or use Cell 3 versions if names don't clash)
# (GEE_TMMC_Simulator_Cell4 and FreeEnergyCalculator_Cell4 classes as provided in last full Cell 4 response)
# (NpEncoder class as provided)
# (The main if __name__ == "__main__" block for Cell 4 simulation orchestration)

# --- Cell 4: Full GEE_TMMC_Simulator_Cell4 and FreeEnergyCalculator_Cell4 Classes and Main Execution ---
# (This is the content from the previous "Complete Cell 4" response, with minor adjustments for clarity
#  and to ensure it uses the globally defined helper functions or self-contained versions if necessary)

# (Re-defining calculation functions for Cell 4 to avoid potential conflicts if cells run out of order,
#  though it's better if they are defined once globally if identical)

def calculate_lj_energy_pair_modified_c4(r_sq, lambda_val, epsilon_ms, sigma_ms): # Renamed to avoid clash
    if r_sq < 1e-12: return 1e10 if lambda_val > 1e-9 else 0.0
    sigma_sq_over_r_sq = (sigma_ms**2) / r_sq; sr6 = sigma_sq_over_r_sq**3
    denominator_base = 0.5 * ((1.0 - lambda_val)**2) + sr6
    if abs(denominator_base) < 1e-12 : return 1e10 if lambda_val > 1e-9 else 0.0
    return 4.0 * lambda_val * epsilon_ms * (1.0 / (denominator_base**2) - 1.0 / denominator_base)

def calculate_polymer_solvent_interaction_c4(polymer_coords, solvent_coords, lambda_val, eps_ms, sigma_ms):
    E = 0.0; Npoly=len(polymer_coords);Nsol=len(solvent_coords)
    if Nsol == 0: return 0.0
    for i in range(Npoly): E += sum(calculate_lj_energy_pair_modified_c4(calculate_distance_sq(polymer_coords[i], solvent_coords[k]), lambda_val, eps_ms, sigma_ms) for k in range(Nsol))
    return E

def calculate_lj_energy_pair_standard_c4(r_sq, epsilon, sigma):
    if r_sq < 1e-12: return 1e10; sigma_sq = sigma**2; inv_r_sq = 1.0 / r_sq
    sr6 = (sigma_sq * inv_r_sq)**3; sr12 = sr6**2
    return 4.0 * epsilon * (sr12 - sr6)

def calculate_solvent_solvent_energy_c4(solvent_coords, eps_ss, sigma_ss, box_dim_arr, cutoff_sq=None):
    E=0.0;Nsol=len(solvent_coords)
    if Nsol < 2: return 0.0
    for i in range(Nsol):
        for j in range(i+1, Nsol):
            rij=solvent_coords[j]-solvent_coords[i]; rij=rij-box_dim_arr*np.rint(rij/box_dim_arr); dsq=np.sum(rij**2)
            if cutoff_sq is None or dsq<cutoff_sq: E+=calculate_lj_energy_pair_standard_c4(dsq,eps_ss,sigma_ss)
    return E

class GEE_TMMC_Simulator_Cell4:
    def __init__(self, initial_polymer_coords, initial_solvent_config, config_obj_from_sim, 
                 cv_bin_idx, binned_lib_entries_this_bin, lib_file_path, n_monomers_conf): # n_monomers_conf
        self.polymer_coords = initial_polymer_coords 
        self.initial_solvent_config_template = np.copy(initial_solvent_config) 
        self.solvent_coords = np.copy(initial_solvent_config)
        self.sim_config = config_obj_from_sim 
        self.num_growth_stages = self.sim_config.num_growth_stages
        self.max_stage_idx = self.num_growth_stages - 1 
        self.lambda_values = np.linspace(0, 1, self.num_growth_stages)
        self.beta = 1.0 / (KB * self.sim_config.target_temperature)
        self.box_dim_arr = np.array([self.sim_config.box_side_length] * 3)
        self.current_growth_stage_idx = 0
        self.eta = np.ones(self.num_growth_stages)
        self.collection_matrix_C_i_condensed = np.zeros((self.num_growth_stages, 3))
        self.visits_per_stage = np.zeros(self.num_growth_stages)
        self.mc_cycle_history_log = []
        self.total_solvent_moves_attempted = 0; self.total_solvent_moves_accepted = 0
        self.total_gee_moves_attempted = 0; self.total_gee_moves_accepted = 0
        self.gee_attempt_details_log = []
        self.cv_bin_index = cv_bin_idx
        self.binned_library_entries_for_this_bin = binned_lib_entries_this_bin
        self.library_file_path = lib_file_path
        self.num_monomers = n_monomers_conf # Use passed n_monomers
        self.current_polymer_lib_entry = None

    def _select_and_load_new_polymer(self):
        if not self.binned_library_entries_for_this_bin: return False
        selected_entry = random.choice(self.binned_library_entries_for_this_bin)
        new_polymer_coords = load_specific_conformation_new_format(
            self.library_file_path, selected_entry['start_line_num'], self.num_monomers)
        if new_polymer_coords is not None:
            self.polymer_coords = new_polymer_coords; self.current_polymer_lib_entry = selected_entry
            print(f"    Switched to new polymer: Idx {selected_entry['original_index']} (Rg={selected_entry['cv']:.3f}) for CV Bin {self.cv_bin_index}")
            if len(self.initial_solvent_config_template) > 0: self.solvent_coords = np.copy(self.initial_solvent_config_template)
            elif self.sim_config.num_solvent_particles > 0 : self.solvent_coords = np.random.rand(self.sim_config.num_solvent_particles, 3) * self.sim_config.box_side_length
            else: self.solvent_coords = np.array([])
            return True
        print(f"    ERROR: Failed to load new polymer (Idx {selected_entry['original_index']}) for CV Bin {self.cv_bin_index}.")
        return False

    def get_total_system_energy(self, stage_idx):
        if len(self.solvent_coords) == 0 and self.sim_config.num_solvent_particles == 0: return 0.0
        lambda_val = self.lambda_values[stage_idx]
        U_ps = calculate_polymer_solvent_interaction_c4(self.polymer_coords, self.solvent_coords, lambda_val, self.sim_config.eps_ms, self.sim_config.sigma_ms)
        U_ss = calculate_solvent_solvent_energy_c4(self.solvent_coords, self.sim_config.eps_ss, self.sim_config.sigma_ss, self.box_dim_arr)
        return U_ps + U_ss
        
    def run_solvent_mc_sweep(self):
        # ... (same as before)
        if len(self.solvent_coords) == 0: return 0.0
        current_U_system_val = self.get_total_system_energy(self.current_growth_stage_idx)
        accepted_this_sweep = 0; num_particles_to_move = len(self.solvent_coords)
        self.total_solvent_moves_attempted += num_particles_to_move
        for _ in range(num_particles_to_move):
            p_idx = np.random.randint(num_particles_to_move); orig_coord = np.copy(self.solvent_coords[p_idx])
            translate_solvent_molecule_pbc(self.solvent_coords, p_idx, self.sim_config.solvent_max_displacement, self.box_dim_arr)
            new_U = self.get_total_system_energy(self.current_growth_stage_idx); delta_U = new_U - current_U_system_val
            if random.random() < np.exp(-self.beta * delta_U): current_U_system_val = new_U; accepted_this_sweep += 1
            else: self.solvent_coords[p_idx] = orig_coord
        self.total_solvent_moves_accepted += accepted_this_sweep
        return accepted_this_sweep / num_particles_to_move

    def attempt_growth_stage_transition(self, growth_phase_active, current_total_elementary_steps):
        # ... (same as previous version with symmetric log(2.0) for C_i and walk) ...
        self.total_gee_moves_attempted += 1; j = self.current_growth_stage_idx; U_j = self.get_total_system_energy(j)
        k_target_for_biased_walk = -1 
        if growth_phase_active:
            if j < self.max_stage_idx: k_target_for_biased_walk = j + 1
        else: 
            if j > 0: k_target_for_biased_walk = j - 1
        factor_for_Ci_update = 0.0; log2_correction = 0.0
        if k_target_for_biased_walk != -1:
            k_ci = k_target_for_biased_walk; U_k_ci = self.get_total_system_energy(k_ci); delta_U_ci = U_k_ci - U_j                   
            log_boltzmann_for_Ci = -self.beta * delta_U_ci
            if (j == 0 and k_ci == 1) or (j == self.max_stage_idx - 1 and k_ci == self.max_stage_idx): log2_correction = np.log(2.0)
            elif (j == 1 and k_ci == 0) or (j == self.max_stage_idx and k_ci == self.max_stage_idx - 1): log2_correction = -np.log(2.0)
            log_factor_for_Ci = log_boltzmann_for_Ci + log2_correction
            factor_for_Ci_update = np.exp(log_factor_for_Ci) if log_factor_for_Ci < 0 else 1.0
            factor_for_Ci_update = min(1.0, factor_for_Ci_update)
            if k_ci > j: self.collection_matrix_C_i_condensed[j, 2] += factor_for_Ci_update
            else: self.collection_matrix_C_i_condensed[j, 0] += factor_for_Ci_update
            self.collection_matrix_C_i_condensed[j, 1] += (1.0 - factor_for_Ci_update)
        else: self.collection_matrix_C_i_condensed[j, 1] += 1.0
        p_acc_biased_walk = 0.0; log_eta_ratio_walk = 0.0; eta_ratio_val_walk = 1.0; log_boltzmann_walk = 0.0; delta_U_walk = 0.0; k_walk_actual_target = k_target_for_biased_walk
        if k_target_for_biased_walk != -1:
            delta_U_walk = self.get_total_system_energy(k_walk_actual_target) - U_j; log_boltzmann_walk = -self.beta * delta_U_walk
            if abs(self.eta[j]) > 1e-12:
                eta_ratio_val_walk = self.eta[k_walk_actual_target] / self.eta[j]
                if eta_ratio_val_walk > 1e-12: log_eta_ratio_walk = np.log(eta_ratio_val_walk)
                else: log_eta_ratio_walk = -float('inf') 
            elif delta_U_walk * self.beta > 0 : log_eta_ratio_walk = -float('inf')
            log_p_acc_unbounded_walk = log_eta_ratio_walk + log_boltzmann_walk + log2_correction # Use same log2_correction
            p_acc_biased_for_walk = np.exp(log_p_acc_unbounded_walk) if log_p_acc_unbounded_walk < 0 else 1.0
        accepted = random.random() < p_acc_biased_for_walk
        self.gee_attempt_details_log.append({"step":current_total_elementary_steps,"from_j":j,"to_k":k_walk_actual_target,"U_j":U_j,"U_k":self.get_total_system_energy(k_walk_actual_target) if k_walk_actual_target!=-1 else U_j,"delta_U":delta_U_walk,"eta_j":self.eta[j],"eta_k":self.eta[k_walk_actual_target if k_walk_actual_target!=-1 else j],"eta_ratio":eta_ratio_val_walk,"boltzmann_factor":np.exp(log_boltzmann_walk),"log2_corr_applied":log2_correction,"p_acc_biased_walk":p_acc_biased_for_walk,"accepted":accepted,"factor_for_Ci":factor_for_Ci_update})
        if accepted: self.current_growth_stage_idx=k_walk_actual_target; self.total_gee_moves_accepted+=1; return True
        return False

    def update_eta_biasing_factors(self):
        # ... (Same robust eta update logic as before) ...
        temp_F = np.zeros(self.num_growth_stages)
        for j_stage in range(self.num_growth_stages - 1):
            N_j = np.sum(self.collection_matrix_C_i_condensed[j_stage, :])
            N_jplus1 = np.sum(self.collection_matrix_C_i_condensed[j_stage+1, :])
            P_j_to_jplus1 = self.collection_matrix_C_i_condensed[j_stage, 2] / (N_j + 1e-12)
            P_jplus1_to_j = self.collection_matrix_C_i_condensed[j_stage+1, 0] / (N_jplus1 + 1e-12)
            delta_F_step = 0.0
            if P_jplus1_to_j < 1e-12: delta_F_step = 0.0 if P_j_to_jplus1 < 1e-12 else -20.0 
            elif P_j_to_jplus1 < 1e-12: delta_F_step = 20.0
            else:
                ratio = P_j_to_jplus1 / P_jplus1_to_j
                if ratio > 1e-9 : delta_F_step = -(1.0/self.beta) * np.log(ratio)
                elif ratio < -1e-9: delta_F_step = 20.0
            temp_F[j_stage+1] = temp_F[j_stage] + delta_F_step
        if np.any(np.isfinite(temp_F)):
            valid_F_elements = temp_F[np.isfinite(temp_F)]
            if len(valid_F_elements) == 0: self.eta = np.ones(self.num_growth_stages); return
            min_F_val = np.min(valid_F_elements)
            shifted_F = temp_F - min_F_val 
            probabilities_p_j = np.exp(-self.beta * shifted_F); probabilities_p_j[~np.isfinite(probabilities_p_j)] = 0
            sum_p = np.sum(probabilities_p_j)
            if sum_p > 1e-9:
                probabilities_p_j /= sum_p
                new_eta = 1.0 / (probabilities_p_j + 1e-12) 
                valid_new_eta = new_eta[probabilities_p_j > 1e-12]
                min_valid_new_eta = np.min(valid_new_eta) if len(valid_new_eta) > 0 else 1.0
                if min_valid_new_eta < 1e-9: min_valid_new_eta = 1.0
                normalized_new_eta = new_eta / min_valid_new_eta
                damping_factor = self.sim_config.eta_damping_factor 
                self.eta = (1.0 - damping_factor) * self.eta + damping_factor * normalized_new_eta
            else: self.eta = np.ones(self.num_growth_stages)
        else: self.eta = np.ones(self.num_growth_stages)


    def run_simulation_for_cv_bin(self, equilibration_cycles, production_cycles, 
                                  plot_log_freq, verbose_print_freq_val, 
                                  solvent_sweeps_per_gee_attempt_val):
        # ... (Structure with equilibration/production loops and dynamic polymer selection logic is same) ...
        self.mc_cycle_history_log = []; self.gee_attempt_details_log = []
        self.total_solvent_moves_attempted = 0; self.total_solvent_moves_accepted = 0
        self.total_gee_moves_attempted = 0; self.total_gee_moves_accepted = 0
        current_total_elementary_steps = 0; growth_phase = True
        sweeps_per_block = max(1, solvent_sweeps_per_gee_attempt_val)
        num_polymer_selections = 0
        if self.current_polymer_lib_entry: num_polymer_selections = 1 # Count the initial one
        else: # If not set by main block, select one now
            if not self._select_and_load_new_polymer(): print(f"CRITICAL ERROR: Could not load initial polymer for CV Bin {self.cv_bin_index}"); return None,None,None,None,0,0,0,0,[]

        print(f"  Starting Equilibration Phase ({equilibration_cycles} outer GEE blocks)...")
        for outer_cycle_eq in range(equilibration_cycles):
            for _ in range(sweeps_per_block):
                self.run_solvent_mc_sweep(); current_total_elementary_steps += 1
                if current_total_elementary_steps % plot_log_freq == 0: self.mc_cycle_history_log.append((current_total_elementary_steps, self.current_growth_stage_idx))
                if current_total_elementary_steps % verbose_print_freq_val == 0:
                    energy = self.get_total_system_energy(self.current_growth_stage_idx)
                    print(f"    Equil Step {current_total_elementary_steps:8d} | Blk {outer_cycle_eq:6d} | St: {self.current_growth_stage_idx:2d} | E: {energy:10.4f} | Ph: {'Gr' if growth_phase else 'De'}")
            self.attempt_growth_stage_transition(growth_phase, current_total_elementary_steps)
            self.visits_per_stage[self.current_growth_stage_idx] += 1
            if not growth_phase and self.current_growth_stage_idx == 0:
                print(f"    Returned to stage 0 in Equil (Outer Cycle {outer_cycle_eq}). Reselecting polymer.")
                if self._select_and_load_new_polymer(): num_polymer_selections +=1
                else: print("      Failed to load new polymer, continuing with old one.")
                growth_phase = True 
            elif growth_phase and self.current_growth_stage_idx == self.max_stage_idx: growth_phase = False
            elif not growth_phase and self.current_growth_stage_idx == 0: growth_phase = True # Should be caught by above
            eta_upd_outer_freq = max(1, self.sim_config.eta_update_frequency // sweeps_per_block if sweeps_per_block > 0 else self.sim_config.eta_update_frequency)
            if outer_cycle_eq > 0 and outer_cycle_eq % eta_upd_outer_freq == 0 : self.update_eta_biasing_factors()
        print(f"  Equil Complete. Stage: {self.current_growth_stage_idx}, Elm. Steps: {current_total_elementary_steps}")
        print(f"  Starting Production Phase ({production_cycles} outer GEE blocks)...")
        for outer_cycle_prod in range(production_cycles):
            for _ in range(sweeps_per_block):
                self.run_solvent_mc_sweep(); current_total_elementary_steps +=1
                if current_total_elementary_steps % plot_log_freq == 0: self.mc_cycle_history_log.append((current_total_elementary_steps, self.current_growth_stage_idx))
                if current_total_elementary_steps % verbose_print_freq_val == 0:
                    energy = self.get_total_system_energy(self.current_growth_stage_idx)
                    print(f"    Prod. Step {current_total_elementary_steps:8d} | Blk {outer_cycle_prod:6d} | St: {self.current_growth_stage_idx:2d} | E: {energy:10.4f} | Ph: {'Gr' if growth_phase else 'De'}")
            self.attempt_growth_stage_transition(growth_phase, current_total_elementary_steps)
            self.visits_per_stage[self.current_growth_stage_idx] += 1
            if not growth_phase and self.current_growth_stage_idx == 0:
                print(f"    Returned to stage 0 in Prod (Outer Cycle {outer_cycle_prod}). Reselecting polymer.")
                if self._select_and_load_new_polymer(): num_polymer_selections +=1
                else: print("      Failed to load new polymer, continuing with old one.")
                growth_phase = True
            elif growth_phase and self.current_growth_stage_idx == self.max_stage_idx: growth_phase = False
            elif not growth_phase and self.current_growth_stage_idx == 0: growth_phase = True
            eta_upd_outer_freq = max(1, self.sim_config.eta_update_frequency // sweeps_per_block if sweeps_per_block > 0 else self.sim_config.eta_update_frequency)
            if outer_cycle_prod > 0 and outer_cycle_prod % eta_upd_outer_freq == 0 : self.update_eta_biasing_factors()
        print(f"  Prod Complete. Stage: {self.current_growth_stage_idx}, Elm. Steps: {current_total_elementary_steps}")
        print(f"  Total polymer conformations sampled for this CV bin: {num_polymer_selections}")
        mc_cycles_out = [item[0] for item in self.mc_cycle_history_log]
        stage_indices_out = [item[1] for item in self.mc_cycle_history_log]
        return (self.collection_matrix_C_i_condensed, self.visits_per_stage,
                mc_cycles_out, stage_indices_out,
                self.total_solvent_moves_attempted, self.total_solvent_moves_accepted,
                self.total_gee_moves_attempted, self.total_gee_moves_accepted,
                self.gee_attempt_details_log)

# --- FreeEnergyCalculator_Cell4 Class (Definition as in previous complete Cell 4) ---
# ... (This class remains unchanged) ...
class FreeEnergyCalculator_Cell4:
    def __init__(self, config_obj_from_sim): 
        self.sim_config = config_obj_from_sim
        self.beta = 1.0 / (KB * self.sim_config.target_temperature)
    def calculate_F_i_j_minus_F_i_1(self, C_i_condensed_np_arr):
        num_stages = self.sim_config.num_growth_stages; F_diffs = np.zeros(num_stages)
        for j in range(num_stages - 1):
            N_j = np.sum(C_i_condensed_np_arr[j, :]); N_jplus1 = np.sum(C_i_condensed_np_arr[j+1, :])
            P_j_to_jplus1 = C_i_condensed_np_arr[j, 2] / (N_j + 1e-12)
            P_jplus1_to_j = C_i_condensed_np_arr[j+1, 0] / (N_jplus1 + 1e-12)
            delta_F_step = 0.0
            if P_jplus1_to_j < 1e-12: delta_F_step = 0.0 if P_j_to_jplus1 < 1e-12 else -20.0 
            elif P_j_to_jplus1 < 1e-12: delta_F_step = 20.0
            else:
                ratio = P_j_to_jplus1 / P_jplus1_to_j
                if ratio > 1e-9 : delta_F_step = -(1.0/self.beta) * np.log(ratio)
                elif ratio < -1e-9: delta_F_step = 20.0
            F_diffs[j+1] = F_diffs[j] + delta_F_step
        return F_diffs
    def calculate_delta_F_i1_initial(self, ni_counts_per_bin_dict, n1_count_ref_bin):
        beta_eff = 1.0 / (KB * self.sim_config.target_temperature)
        if abs(self.sim_config.library_temperature - self.sim_config.target_temperature) > 1e-3: print(f"Warning: Library/Target temps differ. Eq. 12 uses target T.")
        delta_F_i1_dict = {}
        for bin_idx, ni_count in ni_counts_per_bin_dict.items():
            if ni_count == 0 or n1_count_ref_bin == 0: delta_F_i1_dict[bin_idx] = float('inf')
            else: delta_F_i1_dict[bin_idx] = -(1.0/beta_eff) * np.log(float(ni_count) / float(n1_count_ref_bin))
        return delta_F_i1_dict
    def combine_free_energies(self, F_ij_minus_F_i1_all_bins_dict, delta_F_i1_initial_dict):
        final_delta_F_profiles = {}
        for cv_bin_idx, f_diffs_profile_np_list in F_ij_minus_F_i1_all_bins_dict.items():
            f_diffs_profile_np = np.array(f_diffs_profile_np_list)
            delta_Fi1 = delta_F_i1_initial_dict.get(cv_bin_idx, float('inf'))
            if np.isinf(delta_Fi1): final_delta_F_profiles[cv_bin_idx] = np.full_like(f_diffs_profile_np, float('inf'))
            else: final_delta_F_profiles[cv_bin_idx] = f_diffs_profile_np + delta_Fi1
        return final_delta_F_profiles

# Custom JSON encoder
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer): return int(obj)
        if isinstance(obj, np.floating): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return super(NpEncoder, self).default(obj)

# --- Main Execution Block for Cell 4 ---
def run_simulation_main_logic(config_obj): # Encapsulate main logic
    global sim_config # Make sure sim_config is the one from Cell 2
    sim_config = config_obj

    print(f"TMMC simulations: {sim_config.num_cv_bins} CV bins. Equil: {sim_config.mc_equilibration_cycles} outer blocks, Prod: {sim_config.mc_production_cycles} outer blocks.")
    print(f"Each outer block: {sim_config.solvent_sweeps_per_gee_attempt} solvent sweeps then 1 GEE attempt.")
    print(f"Verbose print freq (elem. steps): {sim_config.verbose_print_freq}, Plot log freq (elem. steps): {sim_config.plot_data_log_freq}")
    print(f"Collection matrix C_i and Biased Walk will include log(2.0) boundary factors for ALL boundaries (0<=>1 and max-1<=>max).")
    
    processed_library_path = os.path.join(sim_config.output_dir, "library_processed_index.txt")
    library_index_data = []; all_cv_values_from_index = []
    try:
        with open(processed_library_path, 'r') as f:
            for line in f:
                if line.startswith("#"): continue; parts = line.split()
                if len(parts) < 4: continue
                library_index_data.append({"original_index": int(parts[0]), "start_line_num": int(parts[1]), "cv": float(parts[2]), "energy": float(parts[3])})
                all_cv_values_from_index.append(float(parts[2]))
        if not library_index_data: raise FileNotFoundError("Processed library index is empty.")
        print(f"Loaded {len(library_index_data)} entries from '{processed_library_path}'.")
    except FileNotFoundError as e: print(f"FATAL ERROR: {e}. Run Cell 3."); raise SystemExit(f"Exiting: {e}")

    # Binning logic
    if sim_config.cv_min is None and all_cv_values_from_index: sim_config.cv_min = min(all_cv_values_from_index)
    elif sim_config.cv_min is None: sim_config.cv_min = 0.0
    if sim_config.cv_max is None and all_cv_values_from_index: sim_config.cv_max = max(all_cv_values_from_index)
    elif sim_config.cv_max is None: sim_config.cv_max = sim_config.cv_min
    if abs(sim_config.cv_max - sim_config.cv_min) < 1e-9: sim_config.cv_bin_width = 0.0
    elif sim_config.num_cv_bins > 0 : sim_config.cv_bin_width = (sim_config.cv_max - sim_config.cv_min) / sim_config.num_cv_bins
    else: sim_config.cv_bin_width = 0.0; print("Warning: num_cv_bins is 0 or negative.")
    print(f"Using CV range for binning: {sim_config.cv_min:.4f} to {sim_config.cv_max:.4f}, width: {sim_config.cv_bin_width:.4f if sim_config.cv_bin_width else 'N/A'}")
    binned_library_indices = [[] for _ in range(sim_config.num_cv_bins)]; ni_counts = np.zeros(sim_config.num_cv_bins, dtype=int)
    for entry in library_index_data:
        cv_val = entry["cv"]; bin_idx = 0
        if sim_config.cv_bin_width > 1e-9 : bin_idx = math.floor((cv_val - sim_config.cv_min) / sim_config.cv_bin_width)
        elif not (sim_config.cv_min -1e-9 <= cv_val <= sim_config.cv_max + 1e-9) : continue
        bin_idx = max(0, min(sim_config.num_cv_bins - 1, bin_idx))
        binned_library_indices[bin_idx].append(entry); ni_counts[bin_idx] += 1
    print(f"Conformations per CV bin (ni_counts): {ni_counts}")

    all_F_ij_minus_F_i1 = {}; all_C_i_condensed_dict = {}; all_visits_dict = {}
    all_mc_cycle_histories_dict = {}; all_stage_histories_dict = {}
    all_gee_attempt_details_log_dict = {}
    cumulative_solvent_attempted = 0; cumulative_solvent_accepted = 0
    cumulative_gee_attempted = 0; cumulative_gee_accepted = 0
    processed_bins_count = 0
    plot_bins_indices = []
    non_empty_bins = [i for i, count in enumerate(ni_counts) if count > 0]
    if non_empty_bins:
        plot_bins_indices.append(non_empty_bins[0])
        if len(non_empty_bins) > 2 : plot_bins_indices.append(non_empty_bins[len(non_empty_bins)//2])
        if len(non_empty_bins) > 1 : plot_bins_indices.append(non_empty_bins[-1])
        plot_bins_indices = sorted(list(set(plot_bins_indices)))
    print(f"Will plot detailed MC history & GEE logs for CV bins: {plot_bins_indices}")

    for i_cv_bin in range(sim_config.num_cv_bins):
        print(f"\n--- Processing CV Bin {i_cv_bin} ---")
        if ni_counts[i_cv_bin] == 0: print(f"CV Bin {i_cv_bin} is empty. Skipping."); continue
        processed_bins_count +=1
        first_lib_entry_for_bin = random.choice(binned_library_indices[i_cv_bin])
        initial_polymer_coords = load_specific_conformation_new_format(
            sim_config.conformation_library_file, first_lib_entry_for_bin['start_line_num'], sim_config.num_monomers)
        if initial_polymer_coords is None: print(f"  ERROR: Could not load initial polymer for bin {i_cv_bin}. Skipping."); continue
        
        initial_solvent_coords_for_sim = np.array([])
        if sim_config.num_solvent_particles > 0:
            initial_solvent_coords_for_sim = np.random.rand(sim_config.num_solvent_particles, 3) * sim_config.box_side_length
            if sim_config.initial_solvent_box_file and os.path.exists(sim_config.initial_solvent_box_file):
                try:
                    loaded_solvent = np.loadtxt(sim_config.initial_solvent_box_file)
                    if loaded_solvent.ndim == 2 and loaded_solvent.shape[1] == 3:
                        if len(loaded_solvent) >= sim_config.num_solvent_particles:
                           initial_solvent_coords_for_sim = loaded_solvent[:sim_config.num_solvent_particles]; print(f"  Loaded initial solvent from {sim_config.initial_solvent_box_file}")
                except Exception as e: print(f"  Warning: Error loading solvent box: {e}. Using random.")
        else: print("  Running with NO SOLVENT MOLECULES.")

        gee_simulator = GEE_TMMC_Simulator_Cell4(
            initial_polymer_coords, initial_solvent_coords_for_sim, sim_config,
            i_cv_bin, binned_library_indices[i_cv_bin], 
            sim_config.conformation_library_file, sim_config.num_monomers)
        gee_simulator.current_polymer_lib_entry = first_lib_entry_for_bin

        sim_results = gee_simulator.run_simulation_for_cv_bin(
                sim_config.mc_equilibration_cycles, sim_config.mc_production_cycles, 
                sim_config.plot_data_log_freq, sim_config.verbose_print_freq,
                sim_config.solvent_sweeps_per_gee_attempt)
        if sim_results is None : # Handle potential critical error in simulator
            print(f"  ERROR: GEE simulation failed for CV Bin {i_cv_bin}. Skipping.")
            continue
        C_i_cond, visits_i, mc_hist, stage_hist, solv_att, solv_acc, gee_att, gee_acc, gee_details = sim_results
        
        all_C_i_condensed_dict[i_cv_bin] = C_i_cond.tolist()
        all_visits_dict[i_cv_bin] = visits_i.tolist()
        cumulative_solvent_attempted += solv_att; cumulative_solvent_accepted += solv_acc
        cumulative_gee_attempted += gee_att; cumulative_gee_accepted += gee_acc
        if i_cv_bin in plot_bins_indices:
            all_mc_cycle_histories_dict[i_cv_bin] = mc_hist
            all_stage_histories_dict[i_cv_bin] = stage_hist
            all_gee_attempt_details_log_dict[i_cv_bin] = gee_details
        fe_calculator = FreeEnergyCalculator_Cell4(sim_config)
        F_diffs_for_bin = fe_calculator.calculate_F_i_j_minus_F_i_1(C_i_cond)
        all_F_ij_minus_F_i1[i_cv_bin] = F_diffs_for_bin
        print(f"  Bin {i_cv_bin}: F_i(target) - F_i(library_state) = {F_diffs_for_bin[-1]:.4f} (kBT_target units)")

    if processed_bins_count > 0:
        # Plotting, Metrics, Collection Matrix Printout, Saving Data (same as before)
        # ... (All the plotting, metrics printing, C_i printing, and JSON saving code from the previous full Cell 4 response) ...
        print("\n--- Generating Plots ---")
        plt.figure(figsize=(10,6)); 
        for i_p, mc_h in all_mc_cycle_histories_dict.items():
            st_h = all_stage_histories_dict[i_p]
            if mc_h and st_h: plt.plot(mc_h,st_h,marker='.',ls='-',ms=1,alpha=0.7,label=f'CV Bin {i_p}')
        plt.xlabel(f"Total Elementary MC Steps"); plt.ylabel("Growth Stage Index")
        plt.title("Growth Stage vs. MC Steps (for selected CV Bins)");
        if all_mc_cycle_histories_dict : plt.legend(markerscale=5); plt.grid(True); plt.tight_layout()
        plt.savefig(os.path.join(sim_config.output_dir, "growth_stage_vs_mc_steps.png")); plt.show()
        
        fe_calc_final = FreeEnergyCalculator_Cell4(sim_config); ref_bin = 0
        while ref_bin < sim_config.num_cv_bins and ni_counts[ref_bin] == 0: ref_bin +=1
        final_FEs = {}; delta_F_i1_plot = {}
        if ref_bin < sim_config.num_cv_bins :
            n1_ref = ni_counts[ref_bin]; ni_dict = {i:c for i,c in enumerate(ni_counts)}
            delta_F_i1_plot = fe_calc_final.calculate_delta_F_i1_initial(ni_dict,n1_ref)
            ref_shift = delta_F_i1_plot.get(ref_bin,0.0)
            if not np.isinf(ref_shift):
                for k_idx in delta_F_i1_plot: 
                    if not np.isinf(delta_F_i1_plot[k_idx]): delta_F_i1_plot[k_idx] -= ref_shift
            final_FEs = fe_calc_final.combine_free_energies(all_F_ij_minus_F_i1, delta_F_i1_plot)
        plt.figure(figsize=(12,7)); lambda_ax = np.linspace(0,1,sim_config.num_growth_stages); plotted_count=0
        for i_p in plot_bins_indices:
             if i_p in all_F_ij_minus_F_i1 and ni_counts[i_p]>0:
                f_prof=all_F_ij_minus_F_i1[i_p]
                if f_prof is not None and not np.all(np.isnan(f_prof)):
                    cv_val_center=sim_config.cv_min+(i_p+0.5)*sim_config.cv_bin_width if sim_config.cv_bin_width else sim_config.cv_min
                    plt.plot(lambda_ax,f_prof,marker='o',ls='-',ms=3,label=f'CV Bin {i_p} ($R_g \sim {cv_val_center:.2f}$)')
                    plotted_count+=1
        plt.xlabel("$\lambda$ (Growth Parameter)"); plt.ylabel("$F_i(\lambda)-F_i(0)$ (kBT units)")
        plt.title("$F_i(j)-F_i(1)$ vs Growth (Selected Bins)")
        if plotted_count>0: plt.legend(fontsize='small',loc='best'); plt.grid(True); plt.tight_layout()
        plt.savefig(os.path.join(sim_config.output_dir, "fe_vs_growth_stages_Fi1.png")); plt.show()
        
        if plot_bins_indices and plot_bins_indices[0] in all_gee_attempt_details_log_dict:
            first_plotted_bin = plot_bins_indices[0]
            print(f"\n--- GEE Move Attempt Details for CV Bin {first_plotted_bin} (last up to 100 attempts) ---")
            details_to_print = all_gee_attempt_details_log_dict[first_plotted_bin][-100:]
            for entry in details_to_print: print(f"  Stp:{entry['step']}, {entry['from_j']}->{entry['to_k']}, dU:{entry['delta_U']:.2f}, eta_r:{entry['eta_ratio']:.2e}, exp(-bdU):{entry['boltzmann_factor']:.2e}, P_acc_walk:{entry['p_acc_biased_walk']:.2e}, log2c_walk:{entry['log2_corr_walk']:.2f}, factor_Ci:{entry.get('factor_for_Ci','N/A'):.2e}, Acc:{entry['accepted']}")
        
        print("\n--- Simulation Performance Metrics ---")
        avg_s_acc = (cumulative_solvent_accepted / cumulative_solvent_attempted) * 100 if cumulative_solvent_attempted > 0 else 0
        avg_g_acc = (cumulative_gee_accepted / cumulative_gee_attempted) * 100 if cumulative_gee_attempted > 0 else 0
        print(f"  Processed {processed_bins_count} non-empty CV Bins.")
        print(f"  Overall Solvent Move Acceptance: {avg_s_acc:.2f}% ({cumulative_solvent_accepted}/{cumulative_solvent_attempted})")
        print(f"  Overall GEE Stage Transition Acceptance: {avg_g_acc:.2f}% ({cumulative_gee_accepted}/{cumulative_gee_attempted})")

        print("\n--- Final Condensed Collection Matrices ---")
        for i_cb, C_list in all_C_i_condensed_dict.items():
            if ni_counts[i_cb] > 0 :
                print(f"CV Bin {i_cb}:"); C_np = np.array(C_list)
                for s_j in range(sim_config.num_growth_stages): print(f"  S {s_j:2d}: C(j-1)={C_np[s_j,0]:<10.2f} C(j)={C_np[s_j,1]:<10.2f} C(j+1)={C_np[s_j,2]:<10.2f}")
                print("-" * 70)
        
        output_summary_path = os.path.join(sim_config.output_dir, "tmmc_summary_data.json")
        try: 
            s_Fij_m_Fi1 = {str(k): (v.tolist() if isinstance(v,np.ndarray) else v) for k,v in all_F_ij_minus_F_i1.items()}
            s_dFi1 = {str(k):v for k,v in delta_F_i1_values_plot.items()}
            s_finalFE = {str(k): (v.tolist() if isinstance(v,np.ndarray) else v) for k,v in final_FEs.items()}
            config_to_save = {k: v for k, v in vars(sim_config).items() if not k.startswith('_') and not callable(v) and k != 'KB' and not isinstance(v, type(plt)) and k not in ['gee_attempt_details_log','mc_cycle_history_log']}
            if 'KB_REDUCED' in config_to_save: del config_to_save['KB_REDUCED']
            summary_data = {"F_ij_minus_F_i1_profiles": s_Fij_m_Fi1, "delta_F_i1_values_vs_F_ref1": s_dFi1, "final_FE_profiles_DeltaFij_vs_F_ref1": s_finalFE,
                            "collection_matrices_condensed": all_C_i_condensed_dict, "visit_histograms": all_visits_dict,
                            "ni_counts_per_bin": ni_counts.tolist(), "performance_metrics": {"proc_bins": processed_bins_count, 
                            "solv_att": cumulative_solvent_attempted, "solv_acc": cumulative_solvent_accepted,
                            "gee_att": cumulative_gee_attempted, "gee_acc": cumulative_gee_accepted},
                            "config_used": config_to_save,
                            "gee_attempt_details_first_plotted_bin": all_gee_attempt_details_log_dict.get(plot_bins_indices[0] if plot_bins_indices else -1, [])[-100:] }
            with open(output_summary_path, 'w') as f: json.dump(summary_data, f, indent=2, cls=NpEncoder)
            print(f"Saved all summary data to {output_summary_path}")
        except Exception as e: print(f"Error saving summary JSON: {e}")
    else: print("No data processed or library index was empty. Skipping further analysis.")
    print("\n--- Cell 4 TMMC Simulation (Complete with Dynamic Polymer & Symmetrical Boundary Corr.) finished. ---")

elif 'sim_config' in globals() and exit_flag_cell4:
    print("Cell 4 execution skipped due to prerequisite errors.")
elif 'sim_config' not in globals():
     print("Cell 4 execution skipped: sim_config not found.")


# --- Main Execution Orchestration ---
# This part would typically be run once at the start of the notebook or script
# after defining all functions and classes.

if __name__ == "__main_runner__": # Use a different guard if pasting all cells together
    print("Starting Full TMMC Simulation Workflow...")
    
    # Cell 1 part: Create config (if not already present)
    config_json_path = "input_config.json"
    if not os.path.exists(config_json_path):
        create_input_config_json(config_json_path)
    else:
        print(f"Using existing '{config_json_path}'")

    # Cell 2 part: Load config
    sim_config = load_simulation_config(config_json_path)
    
    if sim_config:
        # Cell 3 part: Process library
        # Ensure your polymer library file (e.g., /content/my_polymer_library.dat)
        # is uploaded to Colab before running this.
        if not os.path.exists(sim_config.conformation_library_file):
            print(f"FATAL ERROR: Polymer library file '{sim_config.conformation_library_file}' not found for Cell 3 processing!")
            print("Please upload it to Colab at the correct path specified in input_config.json and restart.")
        else:
            library_processing_successful = process_library_cell3(sim_config)
        
            if library_processing_successful:
                # Cell 4 part: Run TMMC simulation
                # (The main logic of Cell 4 is already under `if __name__ == "__main__" ...`)
                # To run it if this combined script is executed, we can call a wrapper
                # or rely on the fact that its `if __name__ == "__main__"` block will run if this
                # script is executed directly.
                # For clarity, let's assume the main block of Cell 4 (pasted above) will run
                # if this whole file is executed as a script and `sim_config` is set.
                # The `exit_flag_cell4` will control if it proceeds.
                
                # Re-check flags after Cell 3
                if 'sim_config' in globals() and sim_config is not None:
                    if not os.path.exists(os.path.join(sim_config.output_dir, "library_processed_index.txt")):
                        print("Cell 4 will be skipped as library_processed_index.txt was not created by Cell 3 part.")
                        exit_flag_cell4 = True # Ensure it's set if cell 3 failed
                    else:
                        exit_flag_cell4 = False # Allow Cell 4 main logic to run
                
                if not exit_flag_cell4:
                    # The main Cell 4 logic defined above will execute due to the
                    # `if __name__ == "__main__" and 'sim_config' in globals() and not exit_flag_cell4:` guard
                    print("Proceeding to Cell 4 main logic execution...")
                    # (No explicit call needed here if the large if __name__ block is part of this script)
                else:
                    print("Cell 4 execution was flagged to be skipped.")
            else:
                print("Cell 3 library processing failed. Skipping Cell 4.")
    else:
        print("Failed to load simulation configuration. Aborting.")

    print("\nFull TMMC Simulation Workflow Script finished.")

# To run this in Colab if pasted as one giant cell, the `if __name__ == "__main_runner__":`
# block would be the entry point. You would remove the separate `if __name__ == "__main__":`
# guard from the Cell 4 main logic and just let it be part of the sequential execution,
# or call a main function that encapsulates Cell 4's core logic.
# For simplicity of direct pasting, the `if __name__ == "__main__" ...` block for Cell 4's
# core logic will execute if this file is run as a script.
# If pasting into one cell, just ensure `sim_config` is populated by the Cell 2 part before Cell 3/4 parts.

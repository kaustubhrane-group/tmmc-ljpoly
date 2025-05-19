# tmmc-ljpoly
# TMMC Simulation for Polymer Free Energy Calculation

This repository contains Python scripts to perform Transition Matrix Monte Carlo (TMMC)
simulations coupled with a Growth Expanded Ensemble (GEE) to calculate the
relative free energies of a polymer in an explicit solvent, as a function of a
collective variable (e.g., Radius of Gyration, Rg). The methodology is based on
the approach described in "Scalable Free Energy Computation of Polymers in Explicit
Solvent Using TMMC and Pre-Generated Conformation Libraries" (tmmc1.pdf).

The simulation is structured to be run in a Jupyter/Google Colab environment,
broken down into logical cells or combined into a single script.

## Simulation Workflow

The simulation proceeds in the following main stages:

1.  **Configuration Setup (Conceptual Cell 1 & 2):**
    * A JSON file (`input_config.json`) is used to define all simulation parameters.
        An example is created if not present.
    * These parameters are loaded into a Python configuration object (`sim_config`).
    * Key parameters include:
        * Paths to input (polymer library) and output files.
        * System details: number of monomers, number of solvent particles, box size.
        * Interaction parameters: Lennard-Jones ($\epsilon, \sigma$) for polymer-polymer,
            solvent-solvent, and polymer-solvent interactions; harmonic bond parameters.
        * TMMC/GEE parameters: number of CV bins, number of growth stages ($\lambda$ stages),
            equilibration and production MC cycles, frequency of GEE moves and $\eta$ updates,
            $\eta$ damping factor.
        * Plotting and logging frequencies.

2.  **Polymer Library Pre-processing (Conceptual Cell 3):**
    * Reads a pre-generated polymer conformation library file. The expected format is:
        ```
        # Conformation 1
        x1 y1 z1
        x2 y2 z2
        ... (N lines for N monomers)
        # Conformation 2
        x1 y1 z1
        ...
        ```
    * For each conformation:
        * Calculates the chosen Collective Variable (CV), e.g., Radius of Gyration ($R_g$)
            with the first monomer as the origin.
        * Calculates the total intramolecular energy (bonded + non-bonded Lennard-Jones
            interactions in vacuum).
    * Sorts all conformations based on their CV values.
    * Updates the `cv_min` and `cv_max` in `sim_config` based on the observed range in the library.
    * Writes an output file (`library_processed_index.txt`) containing:
        `Original_Index Start_Line_of_Comment(1-based) CV_Value Total_Energy`.
        This index file is used by the main simulation to avoid reprocessing the
        large library file repeatedly.

3.  **TMMC-GEE Simulation (Conceptual Cell 4):**
    * Loads the `library_processed_index.txt`.
    * Assigns library conformations to CV bins based on their CV values and `sim_config`.
    * For each non-empty CV bin:
        * **Dynamic Polymer Selection:** An initial polymer conformation belonging to the
            current CV bin is randomly selected from the library index. Its full coordinates
            are loaded from the original large library file using the stored start line number.
            Subsequently, each time the GEE simulation for this CV bin completes a full
            growth and de-growth cycle (returning to $\lambda=0$), a *new* polymer conformation
            is randomly selected from the same CV bin and its coordinates are loaded.
            The collection matrix and other statistics are accumulated over all polymers sampled
            for that CV bin.
        * **Solvent Initialization:** Solvent molecules are placed around the polymer. If an
            `initial_solvent_box_file` is provided in the config, it's used; otherwise,
            random placement occurs (if `num_solvent_particles > 0`).
        * **GEE Simulation:**
            * **Equilibration Phase:** Runs for `sim_config.mc_equilibration_cycles` (outer GEE blocks).
            * **Production Phase:** Runs for `sim_config.mc_production_cycles` (outer GEE blocks).
            * In each outer GEE block:
                * Performs `sim_config.solvent_sweeps_per_gee_attempt` solvent MC sweeps
                    (each sweep is $N_{solvent}$ individual particle translation attempts).
                * Attempts one GEE stage transition ($j \to j \pm 1$).
                    * The **biased walk** (actual MC move) uses acceptance criteria from
                        "tmmc1.pdf" Eqs. (4)-(7), incorporating biasing factors $\eta_j$ and
                        $\log(2.0)$ boundary corrections for transitions involving the
                        first/last stages (factors of 2 or 0.5).
                    * The **collection matrix $C_i$** is updated using a factor derived from
                        the unbiased Boltzmann term $e^{-\beta \Delta U}$ and the same
                        $\log(2.0)$ boundary corrections (factors of 2 or 0.5), but **without**
                        the $\eta_j$ factors, as per user specification for detailed balance in $C_i$.
            * **Biasing Factor Update:** $\eta_j$ values are periodically updated based on the
                current state of the collection matrix $C_i$.
            * Periodic printing of MC step, current growth stage, and system energy.
    * **Free Energy Calculation:**
        * $F_i(j) - F_i(1)$ is calculated for each CV bin $i$ from its final $C_i$ (manuscript Eq. 10-11).
        * $\Delta F_i(1) = F_i(1) - F_1(1)$ is calculated from library counts $n_i$ (manuscript Eq. 12).
        * $\Delta F_i(j) = (F_i(j) - F_i(1)) + \Delta F_i(1)$ gives the final profiles (manuscript Eq. 13).
    * **Output and Analysis:**
        * Plots "Growth stage vs. MC steps" for selected CV bins.
        * Plots "Free energy ($F_i(\lambda) - F_i(0)$) vs. growth stages ($\lambda$)" for selected CV bins.
        * Prints performance metrics (solvent move acceptance, GEE transition acceptance).
        * Prints the final condensed collection matrices for processed CV bins.
        * Saves a summary of all results (FE profiles, $C_i$ matrices, visit histograms,
            metrics, config used) to a JSON file (`tmmc_summary_data.json`).

## Setup and Running

**Prerequisites:**
* Python 3
* NumPy
* Matplotlib
* (If running locally outside Colab, ensure these are installed)

**Configuration (`input_config.json`):**
* Create or modify the `input_config.json` file (or the `config_content` dictionary in the Python script) to set all simulation parameters. Pay close attention to:
    * `conformation_library_file`: Path to your polymer coordinates.
    * `num_monomers`: Must match the library.
    * `num_solvent_particles`: Set to 0 if you want to test the "no solvent" case where $F_i(\lambda)-F_i(0)$ should be 0.
    * Energy parameters (`k_harmonic`, `r_eq_harmonic`, `eps_pp`, `sigma_pp`, `eps_ss`, `sigma_ss`, `eps_ms`, `sigma_ms`). Ensure units are consistent.
    * MC cycle counts (`mc_equilibration_cycles`, `mc_production_cycles`). These are "outer GEE block" counts.
    * `solvent_sweeps_per_gee_attempt`: Number of solvent MC sweeps per GEE attempt.
    * Frequencies for updates and logging (`eta_update_frequency`, `plot_data_log_freq`, `verbose_print_freq`).

**Running in Google Colab:**

1.  **Prepare Files:**
    * Upload your polymer conformation library file (e.g., `my_polymer_library.dat`) to your Colab session (e.g., to `/content/my_polymer_library.dat`). Ensure the format matches the one expected by `yield_polymer_conformations_new_format_c3`.
    * If using an initial solvent box, upload that file too.
2.  **Run the Code:**
    * Paste the entire combined Python code into a single Colab cell.
    * Modify the `config_content` dictionary at the beginning of the script if you haven't created a separate `input_config.json` file. Ensure paths like `conformation_library_file` point to your uploaded files in Colab.
    * Execute the cell.

**Expected Output:**

* Printouts indicating the creation and loading of `input_config.json`.
* Progress messages from Cell 3 (library processing).
* An output file `library_processed_index.txt` in the specified `output_dir`.
* Verbose printouts from Cell 4 during equilibration and production phases, showing MC step, current growth stage, and system energy.
* Plots displayed in Colab and saved as PNG files in `output_dir`:
    * `growth_stage_vs_mc_steps.png`
    * `fe_vs_growth_stages_Fi1.png`
* Printed performance metrics (acceptance rates).
* Printed final condensed collection matrices.
* A `tmmc_summary_data.json` file in `output_dir` containing all key results.

## Notes on the Code and Methodology

* **Units:** The code assumes reduced units where $k_B=1.0$. Ensure your input energy and temperature parameters are consistent.
* **Large Simulations:** For production-quality results, `mc_equilibration_cycles` and `mc_production_cycles` will need to be very large, as suggested by the "tmmc1.pdf" paper. This may require significant runtime.
* **$\eta$ Convergence:** The convergence of biasing factors $\eta_j$ is crucial. Monitor their behavior if possible and ensure sufficient equilibration. The current update mechanism is a common one but might need tuning for very rugged landscapes.
* **Collection Matrix Update:** The current $C_i$ update includes user-specified boundary corrections. The standard TMMC free energy equations (manuscript Eqs. 10-11) assume $C_i$ is built from counts proportional to $\min(1, e^{-\beta \Delta U})$. If the boundary corrections significantly alter the nature of $C_i$ counts, the direct application of those equations might need re-evaluation for strict theoretical consistency, or the interpretation of the resulting "free energies" might need to consider these specific definitions.

## Debugging Low GEE Acceptance

If growth stage transitions are not accepted:
1.  Check the **`gee_attempt_details_log`** printed for selected CV bins. This shows $\Delta U$, $\eta$ ratios, boundary corrections, and the final biased acceptance probability for each GEE attempt.
2.  **Energy Barriers:** Very large positive $\Delta U$ will prevent transitions. This could be due to too few growth stages or issues with initial solvent placement.
3.  **$\eta$ Factor Adaptation:** Requires long equilibration. Ensure $C_i$ is getting populated with off-diagonal elements.
4.  **Simulation Time:** TMMC often needs very long runs.

# General settings
# ================

vpath %.npz data


# Preliminary targets
# ===================

# This target creates an empty data directory
data:
	@echo "Creating the output data directory"
	mkdir data

# Computational targets
# =====================

# Cherenkov radiation and seed solitons
# -------------------------------------

# A really clean soliton with almost no Cherenkov radiation
1.064_seed.npz: data
	@echo "Launching a seed soliton at ω₁=1.064"
	./code/scripts/run_seed_soliton_propagation.py -f1 1.064 ./data/1.064_seed.npz

# A fine-spatial-grid (FSG) sample of the previous solution
1.064_fsg.npz: data 1.064_seed.npz
	@echo "Running a FSG simulation for ω₁=1.064"
	./code/scripts/run_continue_solution_fsg.py ./data/1.064_seed.npz ./data/1.064_fsg.npz

# A solution that produces a nice Cherenkov band due to higher-order
# dispersion effects and FWVM process
1.000_seed.npz: data
	@echo "Launching a seed soliton at ω₁=1.000"
	./code/scripts/run_seed_soliton_propagation.py -f1 1.000 ./data/1.000_seed.npz

# A fine-spatial-grid (FSG) sample of the previous solution
1.000_fsg.npz: data 1.000_seed.npz
	@echo "Running a FSG simulation for ω₁=1.000"
	./code/scripts/run_continue_solution_fsg.py ./data/1.000_seed.npz ./data/1.000_fsg.npz

.PHONY: seed_npz
seed_npz: 1.064_seed.npz 1.000_seed.npz 1.064_fsg.npz 1.000_fsg.npz;

# Scattering process
# ------------------

# A very clean scattering on an ideal soliton
1.064_1.700.npz: data 1.064_seed.npz
	@echo "Scattering a weak DW with ωᵢ=1.700 at a soliton with ω₁=1.064"
	./code/scripts/run_weak_scattering_on_filtered_soliton.py -fi 1.700 ./data/1.064_seed.npz ./data/1.064_1.700.npz

# A fine-spatial-grid (FSG) sample of the previous solution
1.064_1.700_fsg.npz: data 1.064_1.700.npz
	@echo "Running a FSG simulation for ωᵢ=1.700 and ω₁=1.064"
	./code/scripts/run_continue_solution_fsg.py ./data/1.064_1.700.npz ./data/1.064_1.700_fsg.npz

.PHONY: scattering_npz
scattering_npz: 1.064_1.700.npz 1.064_1.700_fsg.npz;

.PHONY: npz
npz: seed_npz scattering_npz

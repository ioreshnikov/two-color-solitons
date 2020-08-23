# General settings
# ================

vpath %.npz data


# Preliminary targets
# ===================

# This target creates an empty data directory
data:
	@echo "Creating the output data directory"
	mkdir data

# This target creates a new virtualenv and then installs all the
# required packages in there
venv:
	@echo "Building a venv"
	python3 -m venv venv
	. venv/bin/activate
	cd code && python3 -m pip install -e .


# Computational targets
# =====================

# Cherenkov radiation and seed solitons
# -------------------------------------

# A really clean soliton with almost no Cherenkov radiation
1.067_seed.npz: data
	@echo "Launching a seed soliton at ω₁=1.067"
	. venv/bin/activate
	./code/scripts/run_seed_soliton_propagation.py -f1 1.067 ./data/1.067_seed.npz

# A fine-spatial-grid (FSG) sample of the previous solution
1.067_fsg.npz: data 1.067_seed.npz
	@echo "Running a FSG simulation for ω₁=1.067"
	. venv/bin/activate
	./code/scripts/run_continue_solution_fsg.py ./data/1.067_seed.npz ./data/1.067_fsg.npz

# A solution that produces a nice Cherenkov band due to higher-order
# dispersion effects and FWVM process
1.000_seed.npz: data
	@echo "Launching a seed soliton at ω₁=1.000"
	. venv/bin/activate
	./code/scripts/run_seed_soliton_propagation.py -f1 1.000 ./data/1.000_seed.npz

# A fine-spatial-grid (FSG) sample of the previous solution
1.000_fsg.npz: data 1.000_seed.npz
	@echo "Running a FSG simulation for ω₁=1.000"
	. venv/bin/activate
	./code/scripts/run_continue_solution_fsg.py ./data/1.000_seed.npz ./data/1.000_fsg.npz

# This soliton we need for a scattering section
1.100_seed.npz: data
	@echo "Launching a seed soliton at ω₁=1.100"
	. venv/bin/activate
	./code/scripts/run_seed_soliton_propagation.py -f1 1.100 ./data/1.100_seed.npz

# A fine-spatial-grid (FSG) sample of the previous solution
1.100_fsg.npz: data 1.100_seed.npz
	@echo "Running a FSG simulation for ω₁=1.100"
	. venv/bin/activate
	./code/scripts/run_continue_solution_fsg.py ./data/1.100_seed.npz ./data/1.100_fsg.npz

# This soliton we need for a scattering section
1.080_seed.npz: data
	@echo "Launching a seed soliton at ω₁=1.080"
	. venv/bin/activate
	./code/scripts/run_seed_soliton_propagation.py -f1 1.080 ./data/1.080_seed.npz

# A fine-spatial-grid (FSG) sample of the previous solution
1.080_fsg.npz: data 1.080_seed.npz
	@echo "Running a FSG simulation for ω₁=1.080"
	. venv/bin/activate
	./code/scripts/run_continue_solution_fsg.py ./data/1.080_seed.npz ./data/1.080_fsg.npz

.PHONY: seed_npz
seed_npz: 1.067_seed.npz 1.000_seed.npz 1.100_seed.npz 1.080_seed.npz 1.067_fsg.npz 1.000_fsg.npz 1.100_fsg.npz 1.080_fsg.npz

# Scattering processes
# --------------------

# A very clean scattering on an ideal soliton
1.067_1.700.npz: data 1.067_seed.npz
	@echo "Scattering a weak DW with ωᵢ=1.700 at a soliton with ω₁=1.067"
	. venv/bin/activate
	./code/scripts/run_weak_scattering_on_filtered_soliton.py -fi 1.700 ./data/1.067_seed.npz ./data/1.067_1.700.npz

# A fine-spatial-grid (FSG) sample of the previous solution
1.067_1.700_fsg.npz: data 1.067_1.700.npz
	@echo "Running a FSG simulation for ωᵢ=1.700 and ω₁=1.067"
	. venv/bin/activate
	./code/scripts/run_continue_solution_fsg.py ./data/1.067_1.700.npz ./data/1.067_1.700_fsg.npz

# A scattering that should show a lot of secondary resonances
1.100_2.720.npz: data 1.100_seed.npz
	@echo "Scattering a weak DW with ωᵢ=2.720 at a soliton with ω₁=1.100"
	. venv/bin/activate
	./code/scripts/run_weak_scattering_on_filtered_soliton.py -fi 2.720 ./data/1.100_seed.npz ./data/1.100_2.720.npz

# A fine-spatial-grid (FSG) sample of the previous solution
1.100_2.720_fsg.npz: data 1.100_2.720.npz
	@echo "Running a FSG simulation for ωᵢ=2.720 and ω₁=1.100"
	. venv/bin/activate
	./code/scripts/run_continue_solution_fsg.py ./data/1.100_2.720.npz ./data/1.100_2.720_fsg.npz

# A scattering that should show a lot of secondary resonances
1.100_3.370.npz: data 1.100_seed.npz
	@echo "Scattering a weak DW with ωᵢ=3.370 at a soliton with ω₁=1.100"
	. venv/bin/activate
	./code/scripts/run_weak_scattering_on_filtered_soliton.py -fi 3.370 ./data/1.100_seed.npz ./data/1.100_3.370.npz

# A fine-spatial-grid (FSG) sample of the previous solution
1.100_3.370_fsg.npz: data 1.100_3.370.npz
	@echo "Running a FSG simulation for ωᵢ=3.370 and ω₁=1.100"
	. venv/bin/activate
	./code/scripts/run_continue_solution_fsg.py ./data/1.100_3.370.npz ./data/1.100_3.370_fsg.npz

# A scattering that should show a lot of secondary resonances
1.080_3.250.npz: data 1.080_seed.npz
	@echo "Scattering a weak DW with ωᵢ=3.250 at a soliton with ω₁=1.080"
	. venv/bin/activate
	./code/scripts/run_weak_scattering_on_filtered_soliton.py -fi 3.250 ./data/1.080_seed.npz ./data/1.080_3.250.npz

# A fine-spatial-grid (FSG) sample of the previous solution
1.080_3.250_fsg.npz: data 1.080_3.250.npz
	@echo "Running a FSG simulation for ωᵢ=3.250 and ω₁=1.080"
	. venv/bin/activate
	./code/scripts/run_continue_solution_fsg.py ./data/1.080_3.250.npz ./data/1.080_3.250_fsg.npz

.PHONY: scattering_npz
scattering_npz: 1.067_1.700.npz 1.100_2.720.npz 1.100_3.370.npz 1.080_3.250.npz 1.067_1.700_fsg.npz 1.100_2.720_fsg.npz 1.100_3.370_fsg.npz 1.080_3.250_fsg.npz

.PHONY: npz
npz: seed_npz scattering_npz

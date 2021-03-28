# General settings
# ================

vpath %.npz data
vpath %.png figs

# Preliminary targets
# ===================

# This target creates an empty data directory
data:
	@echo "Creating the output data directory"
	mkdir data

# This target creates an empty figures directory
figs:
	@echo "Creating the output figure directory"
	mkdir figs

# Computational targets
# =====================

# Define all the seed frequencies to check
SEED_FREQS = $(shell seq -f "1.%03.0f" 0 5 150)

# This macros defines a seed soliton target
define make-seed-computation-target
$(addsuffix _seed.npz, $1): data
	python code/scripts/run_seed_soliton_propagation.py -f1 $(strip $1) data/$(strip $1)_seed.npz
endef

# This macros defines a filtered soliton target
define make-filt-computation-target
$(addsuffix _filt.npz, $1): data $(addsuffix _seed.npz, $1)
	python code/scripts/run_filtered_soliton_propagation.py data/$(strip $1)_seed.npz data/$(strip $1)_filt.npz
endef

# This macros defines a FSG target for filtered result
define make-fsg-computation-target
$(addsuffix _fsg.npz, $1): data $(addsuffix _filt.npz, $1)
	python code/scripts/run_continue_solution_fsg.py data/$(strip $1)_filt.npz data/$(strip $1)_fsg.npz
endef

# This macros plots the propagation results for the seed soliton
define make-seed-plotting-target
$(addsuffix _seed.png, $1): figs $(addsuffix _seed.npz, $1)
	python code/scripts/fig_time_outspectrum.py data/$(strip $1)_seed.npz figs/$(strip $1)_seed.png
endef

# This macros plots the propagation results for the seed soliton
define make-filt-plotting-target
$(addsuffix _filt.png, $1): figs $(addsuffix _filt.npz, $1)
	python code/scripts/fig_time_outspectrum.py data/$(strip $1)_filt.npz figs/$(strip $1)_filt.png
endef

# This macros defines a plotting target for FSG results
define make-fsg-plotting-target
$(addsuffix _fsg.png, $1): figs $(addsuffix _fsg.npz, $1)
	python code/scripts/fig_rescon.py data/$(strip $1)_fsg.npz figs/$(strip $1)_fsg.png
endef

# Generate seed targets for every frequency
$(foreach freq, $(SEED_FREQS), $(eval $(call make-seed-computation-target, $(freq))))
$(foreach freq, $(SEED_FREQS), $(eval $(call make-filt-computation-target, $(freq))))
$(foreach freq, $(SEED_FREQS), $(eval $(call make-fsg-computation-target,  $(freq))))
$(foreach freq, $(SEED_FREQS), $(eval $(call make-seed-plotting-target, $(freq))))
$(foreach freq, $(SEED_FREQS), $(eval $(call make-filt-plotting-target, $(freq))))
$(foreach freq, $(SEED_FREQS), $(eval $(call make-fsg-plotting-target,  $(freq))))

# And one huge phony target combining all of them
.PHONY: seed_npz
seed_npz: $(addsuffix _seed.npz, $(SEED_FREQS))

.PHONY: seed_png
seed_png: $(addsuffix _seed.png, $(SEED_FREQS))

.PHONY: filt_npz
filt_npz: $(addsuffix _filt.npz, $(SEED_FREQS))

.PHONY: filt_png
filt_png: $(addsuffix _filt.png, $(SEED_FREQS))

.PHONY: fsg_npz
fsg_npz: $(addsuffix _fsg.npz, $(SEED_FREQS))

.PHONY: fsg_png
fsg_png: $(addsuffix _fsg.png, $(SEED_FREQS))

.PHONY: npz
npz: seed_npz filt_npz fsg_npz

.PHONY: png
png: seed_png filt_png fsg_png

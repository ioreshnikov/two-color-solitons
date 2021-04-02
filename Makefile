# General settings
# ================

ext = pdf

vpath %.npz npz
vpath %.$(ext) fig


# Brute-force computational targets
# =================================

# The targets in these section are left from the period of active
# research. They do computations and exploratory plotting on a large
# family of parameters, but almost none of it makes it to the final
# paper. We leave them here for posterity.

# Seed propagation, filtered solutions and Cherenkov radiation
# ------------------------------------------------------------

# Define all the seed frequencies to check
SEED_FREQS = $(shell seq -f "1.%03.0f" 0 5 150)

# This macros defines a seed soliton target
define make-seed-computation-target
$(addsuffix _seed.npz, $1):
	python code/scripts/run_seed_soliton_propagation.py -f1 $(strip $1) npz/$(strip $1)_seed.npz
endef

# This macros defines a filtered soliton target
define make-filt-computation-target
$(addsuffix _filt.npz, $1): $(addsuffix _seed.npz, $1)
	python code/scripts/run_filtered_soliton_propagation.py npz/$(strip $1)_seed.npz npz/$(strip $1)_filt.npz
endef

# This macros plots the propagation results for the seed soliton
define make-seed-plotting-target
$(addsuffix _seed.$(ext), $1): $(addsuffix _seed.npz, $1)
	python code/scripts/fig_time_outspectrum.py npz/$(strip $1)_seed.npz fig/$(strip $1)_seed.$(ext)
endef

# This macros plots the propagation results for the seed soliton
define make-filt-plotting-target
$(addsuffix _filt.$(ext), $1): $(addsuffix _filt.npz, $1)
	python code/scripts/fig_time_outspectrum.py npz/$(strip $1)_filt.npz fig/$(strip $1)_filt.$(ext)
endef

$(foreach freq, $(SEED_FREQS), $(eval $(call make-seed-computation-target, $(freq))))
$(foreach freq, $(SEED_FREQS), $(eval $(call make-filt-computation-target, $(freq))))

$(foreach freq, $(SEED_FREQS), $(eval $(call make-seed-plotting-target, $(freq))))
$(foreach freq, $(SEED_FREQS), $(eval $(call make-filt-plotting-target, $(freq))))

.PHONY: seed_npz
seed_npz: $(addsuffix _seed.npz, $(SEED_FREQS))

.PHONY: seed_$(ext)
seed_$(ext): $(addsuffix _seed.$(ext), $(SEED_FREQS))

.PHONY: filt_npz
filt_npz: $(addsuffix _filt.npz, $(SEED_FREQS))

.PHONY: filt_$(ext)
filt_$(ext): $(addsuffix _filt.$(ext), $(SEED_FREQS))

# Weak dispersive waves scattering
# --------------------------------

INCIDENT_FREQS = $(shell seq -f "%1.3f" 1.100 0.100 3.000)

# This macros defines a scattering target for a fixed seed soliton
define make-weak-scattering-computation-target
$(addprefix $(addsuffix _, $1), $(addsuffix .npz, $2)): $(addsuffix _seed, $1).npz
	python code/scripts/run_weak_scattering_on_filtered_soliton.py -fi $2 npz/$(strip $1)_seed.npz npz/$(strip $1)_$(strip $2).npz
endef

# This macros defines a plotting target for the previous scattering problem
define make-weak-scattering-plotting-target
$(addprefix $(addsuffix _, $1), $(addsuffix .$(ext), $2)): $(addprefix $(addsuffix _, $1), $(addsuffix .npz, $2))
	python code/scripts/fig_time_outspectrum.py npz/$(strip $1)_$(strip $2).npz fig/$(strip $1)_$(strip $2).$(ext)
endef

$(foreach freq, $(INCIDENT_FREQS), $(eval $(call make-weak-scattering-computation-target, 1.005, $(freq))))
$(foreach freq, $(INCIDENT_FREQS), $(eval $(call make-weak-scattering-computation-target, 1.070, $(freq))))
$(foreach freq, $(INCIDENT_FREQS), $(eval $(call make-weak-scattering-computation-target, 1.150, $(freq))))
$(foreach freq, $(INCIDENT_FREQS), $(eval $(call make-weak-scattering-plotting-target,    1.005, $(freq))))
$(foreach freq, $(INCIDENT_FREQS), $(eval $(call make-weak-scattering-plotting-target,    1.070, $(freq))))
$(foreach freq, $(INCIDENT_FREQS), $(eval $(call make-weak-scattering-plotting-target,    1.150, $(freq))))

.PHONY: scattering_1.005_npz
scattering_1.005_npz: $(addprefix 1.005_, $(addsuffix .npz, $(INCIDENT_FREQS)))

.PHONY: scattering_1.070_npz
scattering_1.079_npz: $(addprefix 1.070_, $(addsuffix .npz, $(INCIDENT_FREQS)))

.PHONY: scattering_1.150_npz
scattering_1.150_npz: $(addprefix 1.150_, $(addsuffix .npz, $(INCIDENT_FREQS)))

.PHONY: scattering_1.005_$(ext)
scattering_1.005_$(ext): $(addprefix 1.005_, $(addsuffix .$(ext), $(INCIDENT_FREQS)))

.PHONY: scattering_1.070_$(ext)
scattering_1.070_$(ext): $(addprefix 1.070_, $(addsuffix .$(ext), $(INCIDENT_FREQS)))

.PHONY: scattering_1.150_$(ext)
scattering_1.150_$(ext): $(addprefix 1.150_, $(addsuffix .$(ext), $(INCIDENT_FREQS)))


# Figures and computational targets for the paper
# ===============================================

# This target defines a really nice scattering process which is not
# covered by the bruteforce approach.
1.070_1.650.npz: 1.070_seed.npz
	python code/scripts/run_weak_scattering_on_filtered_soliton.py -fi 1.650 npz/1.070_seed.npz npz/1.070_1.650.npz

Fig1.$(ext): 1.000_seed.npz
	python code/scripts/fig_time_outspectrum.py --no-rescon npz/1.000_seed.npz fig/Fig1.$(ext)

Fig2.$(ext): 1.005_filt.npz
	python code/scripts/fig_time_outspectrum.py --vmin 1E-8 npz/1.005_filt.npz fig/Fig2.$(ext)

Fig3.$(ext): 1.070_filt.npz
	python code/scripts/fig_time_outspectrum.py --vmin 1E-8 npz/1.070_filt.npz fig/Fig3.$(ext)

Fig4.$(ext): 1.150_filt.npz
	python code/scripts/fig_time_outspectrum.py --vmin 1E-8 npz/1.150_filt.npz fig/Fig4.$(ext)

Fig5.$(ext): 1.005_1.600.npz
	python code/scripts/fig_time_outspectrum.py npz/1.005_1.600.npz fig/Fig5.$(ext)

Fig6.$(ext): 1.070_1.650.npz
	python code/scripts/fig_time_outspectrum.py npz/1.070_1.650.npz fig/Fig6.$(ext)

.PHONY: fig
fig: Fig1.$(ext) Fig2.$(ext) Fig3.$(ext) Fig4.$(ext) Fig5.$(ext) Fig6.$(ext)

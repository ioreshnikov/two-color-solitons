# General settings
# ================

ext = png

vpath %.npz npz
vpath %.$(ext) fig


# The computational targets in the sections below are left from the
# period of active research. They do computations and exploratory
# plotting on a large family of parameters, but almost none of it
# makes it to the final paper. We leave them here for posterity.


# Brute-force computational targets
# =================================

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
	python code/scripts/fig_time_outspectrum.py --vmin 1E-8 npz/$(strip $1)_filt.npz fig/$(strip $1)_filt.$(ext)
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

# XXX: One wild guess. An uneven seed.
1.085_ue_seed.npz:
	python code/scripts/run_seed_soliton_propagation.py -f1 1.085 -t1 30.000 -t2 10.000  npz/1.085_ue_seed.npz

1.085_ue_filt.npz: 1.085_ue_seed.npz
	python code/scripts/run_filtered_soliton_propagation.py npz/1.085_ue_seed.npz npz/1.085_ue_filt.npz

1.085_ue_seed.$(ext): 1.085_ue_seed.npz
	python code/scripts/fig_time_outspectrum.py npz/1.085_ue_seed.npz fig/1.085_ue_seed.$(ext)

1.085_ue_filt.$(ext): 1.085_ue_filt.npz
	python code/scripts/fig_time_outspectrum.py --vmin 1E-8 npz/1.085_ue_filt.npz fig/1.085_ue_filt.$(ext)

# Weak dispersive waves scattering
# --------------------------------

INCIDENT_FREQS_1 = $(shell seq -f "%1.3f" 1.100 0.100 3.000)
INCIDENT_FREQS_2 = $(shell seq -f "%1.3f" 3.000 0.100 3.500)

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

$(foreach freq, $(INCIDENT_FREQS_1), $(eval $(call make-weak-scattering-computation-target, 1.005, $(freq))))
$(foreach freq, $(INCIDENT_FREQS_1), $(eval $(call make-weak-scattering-plotting-target,    1.005, $(freq))))

$(foreach freq, $(INCIDENT_FREQS_1), $(eval $(call make-weak-scattering-computation-target, 1.070, $(freq))))
$(foreach freq, $(INCIDENT_FREQS_1), $(eval $(call make-weak-scattering-plotting-target,    1.070, $(freq))))

$(foreach freq, $(INCIDENT_FREQS_1), $(eval $(call make-weak-scattering-computation-target, 1.150, $(freq))))
$(foreach freq, $(INCIDENT_FREQS_1), $(eval $(call make-weak-scattering-plotting-target,    1.150, $(freq))))

$(foreach freq, $(INCIDENT_FREQS_2), $(eval $(call make-weak-scattering-computation-target, 1.085_ue, $(freq))))
$(foreach freq, $(INCIDENT_FREQS_2), $(eval $(call make-weak-scattering-plotting-target,    1.085_ue, $(freq))))

.PHONY: scattering_1.005_npz
scattering_1.005_npz: $(addprefix 1.005_, $(addsuffix .npz, $(INCIDENT_FREQS_1)))

.PHONY: scattering_1.070_npz
scattering_1.070_npz: $(addprefix 1.070_, $(addsuffix .npz, $(INCIDENT_FREQS_1)))

.PHONY: scattering_1.085_npz
scattering_1.150_npz: $(addprefix 1.085_ue_, $(addsuffix .npz, $(INCIDENT_FREQS_1)))

.PHONY: scattering_1.150_npz
scattering_1.150_npz: $(addprefix 1.150_, $(addsuffix .npz, $(INCIDENT_FREQS_1)))

.PHONY: scattering_1.005_$(ext)
scattering_1.005_$(ext): $(addprefix 1.005_, $(addsuffix .$(ext), $(INCIDENT_FREQS_1)))

.PHONY: scattering_1.070_$(ext)
scattering_1.070_$(ext): $(addprefix 1.070_, $(addsuffix .$(ext), $(INCIDENT_FREQS_1)))

.PHONY: scattering_1.150_$(ext)
scattering_1.150_$(ext): $(addprefix 1.150_, $(addsuffix .$(ext), $(INCIDENT_FREQS_1)))

.PHONY: scattering_1.085_ue_$(ext)
scattering_1.085_ue_$(ext): $(addprefix 1.085_ue_, $(addsuffix .$(ext), $(INCIDENT_FREQS_2)))


# Figures and computational targets for the paper
# ===============================================

# These targets define really nice scattering processes that are not
# covered by the bruteforce approach.
1.070_1.650.npz: 1.070_seed.npz
	python code/scripts/run_weak_scattering_on_filtered_soliton.py -fi 1.650 npz/1.070_seed.npz npz/1.070_1.650.npz

Fig1.$(ext): 1.000_seed.npz
	python code/scripts/fig_time_outspectrum.py --no-rescon npz/1.000_seed.npz fig/Fig1.$(ext)

Fig2.$(ext): 1.000_filt.npz
	python code/scripts/fig_time_outspectrum.py --vmin 1E-8 npz/1.000_filt.npz fig/Fig2.$(ext)

Fig3.$(ext): 1.070_1.650.npz
	python code/scripts/fig_time_outspectrum.py --vmin 5E-7 npz/1.070_1.650.npz fig/Fig3.$(ext)

Fig4.$(ext): 1.150_1.500.npz
	python code/scripts/fig_time_outspectrum.py --vmin 5E-7 npz/1.150_1.500.npz fig/Fig4.$(ext)

Fig5.$(ext): 1.085_ue_3.500.npz
	python code/scripts/fig_time_outspectrum.py --vmin 5E-7 npz/1.085_ue_3.500.npz fig/Fig5.$(ext)

.PHONY: fig
fig: Fig1.$(ext) Fig2.$(ext) Fig3.$(ext) Fig4.$(ext) Fig5.$(ext)

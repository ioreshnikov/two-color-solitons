# General settings
# ================

# Default image extension (for exploratory plots)
ext = png

vpath %.npz npz
vpath %.pdf fig
vpath %.$(ext) fig

.SUFFIXES: .npz .png .pdf


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
	python code/scripts/plot_timedomain_spectrum_rescon.py npz/$(strip $1)_seed.npz fig/$(strip $1)_seed.$(ext)
endef

# This macros plots the propagation results for the seed soliton
define make-filt-plotting-target
$(addsuffix _filt.$(ext), $1): $(addsuffix _filt.npz, $1)
	python code/scripts/plot_timedomain_spectrum_rescon.py --vmin 1E-8 npz/$(strip $1)_filt.npz fig/$(strip $1)_filt.$(ext)
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
	python code/scripts/plot_timedomain_spectrum_rescon.py npz/1.085_ue_seed.npz fig/1.085_ue_seed.$(ext)

1.085_ue_filt.$(ext): 1.085_ue_filt.npz
	python code/scripts/plot_timedomain_spectrum_rescon.py --vmin 1E-8 npz/1.085_ue_filt.npz fig/1.085_ue_filt.$(ext)

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
	python code/scripts/plot_timedomain_spectrum_rescon.py npz/$(strip $1)_$(strip $2).npz fig/$(strip $1)_$(strip $2).$(ext)
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

# Long-distance propagation of seed solitons
# ------------------------------------------

# XXX: It takes literally days to run the long-distance simulations.

LD_SEED_FREQS = $(shell seq -f "1.%03.0f" 0 10 150)

# This macros defines a seed soliton target (but this time for a long-distance simulation)
define make-ld-seed-computation-target
$(addsuffix _ld_seed.npz, $1):
	python code/scripts/run_ld_seed_soliton_propagation.py -f1 $(strip $1) npz/$(strip $1)_ld_seed.npz
endef

# This macros plots the propagation results for the seed soliton (the long-distance one)
define make-ld-seed-plotting-target
$(addsuffix _ld_seed.$(ext), $1): $(addsuffix _ld_seed.npz, $1)
	python code/scripts/plot_timedomain_freqdomain.py npz/$(strip $1)_ld_seed.npz fig/$(strip $1)_ld_seed.$(ext)
endef

$(foreach freq, $(LD_SEED_FREQS), $(eval $(call make-ld-seed-computation-target, $(freq))))
$(foreach freq, $(LD_SEED_FREQS), $(eval $(call make-ld-seed-plotting-target, $(freq))))

.PHONY: ld_seed_npz
ld_seed_npz: $(addsuffix _ld_seed.npz, $(LD_SEED_FREQS))

.PHONY: ld_seed_$(ext)
ld_seed_$(ext): $(addsuffix _ld_seed.$(ext), $(LD_SEED_FREQS))

define make-seed-soliton-parameters-target
$(addsuffix _ld_seed_sp.npz, $1): $(addsuffix _ld_seed.npz, $1)
	python code/scripts/extract_soliton_parameters.py npz/$(strip $1)_ld_seed.npz npz/$(strip $1)_ld_seed_sp.npz
endef

$(foreach freq, $(LD_SEED_FREQS), $(eval $(call make-ext-soliton-parameters-target, $(freq))))

.PHONY: ld_seed_sp_npz
ld_seed_sp_npz: $(addsuffix _ld_seed_sp.npz, $(LD_SEED_FREQS))

# Itensive dispersive wave scattering
# -----------------------------------

define make-intensive-scattering-computation-target
$(addprefix $(addsuffix _, $1), $(addsuffix _int.npz, $2)): $(addsuffix _seed, $1).npz
	python code/scripts/run_intensive_scattering_on_filtered_soliton.py -fi $2 npz/$(strip $1)_seed.npz npz/$(strip $1)_$(strip $2)_int.npz
endef

define make-intensive-scattering-soliton-parameters-target
$(addprefix $(addsuffix _, $1), $(addsuffix _int_sp.npz, $2)): $(addprefix $(addsuffix _, $1), $(addsuffix _int.npz, $2))
	python code/scripts/extract_soliton_parameters.py npz/$(strip $1)_$(strip $2)_int.npz npz/$(strip $1)_$(strip $2)_int_sp.npz
endef

define make-intensive-scattering-exploratory-plot-target  # just for exploratory needs
$(addprefix $(addsuffix _, $1), $(addsuffix _int.$(ext), $2)): $(addprefix $(addsuffix _, $1), $(addsuffix _int.npz, $2)) $(addprefix $(addsuffix _, $1), $(addsuffix _int_sp.npz, $2))
	python code/scripts/plot_fig_5.py npz/$(strip $1)_$(strip $2)_int.npz npz/$(strip $1)_$(strip $2)_int_sp.npz fig/$(strip $1)_$(strip $2)_int.$(ext)
endef

$(foreach freq, $(INCIDENT_FREQS_1), $(eval $(call make-intensive-scattering-computation-target,        1.010, $(freq))))
$(foreach freq, $(INCIDENT_FREQS_1), $(eval $(call make-intensive-scattering-soliton-parameters-target, 1.010, $(freq))))
$(foreach freq, $(INCIDENT_FREQS_1), $(eval $(call make-intensive-scattering-exploratory-plot-target,   1.010, $(freq))))

.PHONY: scattering_1.010_int_npz
scattering_1.010_int_npz: $(addprefix 1.010_, $(addsuffix _int.npz, $(INCIDENT_FREQS_1)))

.PHONY: scattering_1.010_int_sp
scattering_1.010_sp_npz: $(addprefix 1.010_, $(addsuffix _int_sp.npz, $(INCIDENT_FREQS_1)))

.PHONY: scattering_1.010_int_$(ext)
scattering_1.010_int_$(ext): $(addprefix 1.010_, $(addsuffix _int.$(ext), $(INCIDENT_FREQS_1)))

# Resonance curve around the internal oscillating mode
# ----------------------------------------------------

# This mostly repeats what is done with the larger intensive scattering, but
# with a lower amplitude to try to work around the nonlinearity of the
# frequency oscillations itself.

RC_FREQS = $(shell seq -f "%1.3f" 1.050 0.010 1.200; seq -f "%1.3f" 1.250 0.050 1.350)

define make-rc-computation-target
$(addprefix $(addsuffix _, $1), $(addsuffix _rc10.npz, $2)): $(addsuffix _seed, $1).npz
	python code/scripts/run_intensive_scattering_on_filtered_soliton.py -ai 0.010 -fi $2 npz/$(strip $1)_seed.npz npz/$(strip $1)_$(strip $2)_rc10.npz
endef

define make-rc-soliton-parameters-target
$(addprefix $(addsuffix _, $1), $(addsuffix _rc10_sp.npz, $2)): $(addprefix $(addsuffix _, $1), $(addsuffix _rc10.npz, $2))
	python code/scripts/extract_soliton_parameters.py npz/$(strip $1)_$(strip $2)_rc10.npz npz/$(strip $1)_$(strip $2)_rc10_sp.npz
endef

$(foreach freq, $(RC_FREQS), $(eval $(call make-rc-computation-target,        1.010, $(freq))))
$(foreach freq, $(RC_FREQS), $(eval $(call make-rc-soliton-parameters-target, 1.010, $(freq))))

.PHONY: scattering_1.010_rc10_npz
scattering_1.010_rc10_npz: $(addprefix 1.010_, $(addsuffix _rc10.npz, $(RC_FREQS)))

.PHONY: scattering_1.010_rc10_sp
scattering_1.010_rc10_sp: $(addprefix 1.010_, $(addsuffix _rc10_sp.npz, $(RC_FREQS)))

define make-rc-exploratory-plot-target  # just for exploratory needs
$(addprefix $(addsuffix _, $1), $(addsuffix _rc10_sp.$(ext), $2)): $(addprefix $(addsuffix _, $1), $(addsuffix _rc10.npz, $2)) $(addprefix $(addsuffix _, $1), $(addsuffix _rc10_sp.npz, $2))
	python code/scripts/plot_fig_5.py npz/$(strip $1)_$(strip $2)_rc10.npz npz/$(strip $1)_$(strip $2)_rc10_sp.npz fig/$(strip $1)_$(strip $2)_rc10_sp.png
endef

$(foreach freq, $(RC_FREQS), $(eval $(call make-rc-exploratory-plot-target, 1.010, $(freq))))

.PHONY: scattering_1.010_rc10_$(ext)
scattering_1.010_rc10_$(ext): $(addprefix 1.010_, $(addsuffix _rc10_sp.$(ext), $(RC_FREQS)))


# Special computational targets
# =============================

# This target defines a really nice scattering processes that is not covered
# by the brute force approach. This is the simulation in Fig2.
1.070_1.650.npz: 1.070_seed.npz
	python code/scripts/run_weak_scattering_on_filtered_soliton.py -fi 1.650 npz/1.070_seed.npz npz/1.070_1.650.npz

# These targets are here to explore the vicinity of Fig2.
1.070_1.625.npz: 1.070_seed.npz
	python code/scripts/run_weak_scattering_on_filtered_soliton.py -fi 1.625 npz/1.070_seed.npz npz/1.070_1.625.npz

1.070_1.675.npz: 1.070_seed.npz
	python code/scripts/run_weak_scattering_on_filtered_soliton.py -fi 1.675 npz/1.070_seed.npz npz/1.070_1.675.npz

fig2a.png: 1.070_1.625.npz
	python code/scripts/plot_timedomain_spectrum_rescon.py --vmin 1E-7 npz/1.070_1.625.npz fig/fig2a.png

fig2b.png: 1.070_1.650.npz
	python code/scripts/plot_timedomain_spectrum_rescon.py --vmin 1E-7 npz/1.070_1.650.npz fig/fig2b.png

fig2c.png: 1.070_1.675.npz
	python code/scripts/plot_timedomain_spectrum_rescon.py --vmin 1E-7 npz/1.070_1.675.npz fig/fig2c.png

# Still trying to figure out the scattering we see in Fig3 :) Here we do a
# separate pass at simulating the scattering field only using the linearized
# equation.

# FIG2
1.070_1.650_lin_pot_12.npz: 1.070_seed.npz
	python code/scripts/run_linearized_scattering_on_filtered_soliton.py npz/1.070_seed.npz npz/1.070_1.650_lin_pot_12.npz -fi 1.650 -t pot.1 12 --title "Fig. 2, (12)"

1.070_1.650_lin_pot_12_13.npz: 1.070_seed.npz
	python code/scripts/run_linearized_scattering_on_filtered_soliton.py npz/1.070_seed.npz npz/1.070_1.650_lin_pot_12_13.npz -fi 1.650 -t pot.1 12 13 --title "Fig. 2, + (12) + (13)"

1.070_1.650_lin_pot_12_13_14.npz: 1.070_seed.npz
	python code/scripts/run_linearized_scattering_on_filtered_soliton.py npz/1.070_seed.npz npz/1.070_1.650_lin_pot_12_13_14.npz -fi 1.650 -t pot.1 12 13 14 --title "Fig. 2, + (12) + (13) + (14)"

1.070_1.650_lin_pot_12_13_14_15.npz: 1.070_seed.npz
	python code/scripts/run_linearized_scattering_on_filtered_soliton.py npz/1.070_seed.npz npz/1.070_1.650_lin_pot_12_13_14_15.npz -fi 1.650 -t pot.1 12 13 14 15 --title "Fig. 2, + (12) + (13) + (14) + (15)"

.PHONY: 1.070_1.650_lin.npz
1.070_1.650_lin.npz: 1.070_1.650_lin_pot_12.npz 1.070_1.650_lin_pot_12_13.npz 1.070_1.650_lin_pot_12_13_14.npz 1.070_1.650_lin_pot_12_13_14_15.npz


# FIG3
1.150_1.500_lin_pot_12.npz: 1.070_seed.npz
	python code/scripts/run_linearized_scattering_on_filtered_soliton.py npz/1.150_seed.npz npz/1.070_1.500_lin_pot_12.npz -fi 1.650 -t pot 12 --title "Fig. 3, Pot. + (12)"

1.150_1.500_lin_pot_12_13.npz: 1.070_seed.npz
	python code/scripts/run_linearized_scattering_on_filtered_soliton.py npz/1.150_seed.npz npz/1.070_1.500_lin_pot_12_13.npz -fi 1.650 -t pot 12 13 --title "Fig. 3, Pot. + (12) + (13)"

1.150_1.500_lin_pot_12_13_14.npz: 1.070_seed.npz
	python code/scripts/run_linearized_scattering_on_filtered_soliton.py npz/1.150_seed.npz npz/1.070_1.500_lin_pot_12_13_14.npz -fi 1.650 -t pot 12 13 14 --title "Fig. 3, Pot. + (12) + (13) + (14)"

1.150_1.500_lin_pot_12_13_14_15.npz: 1.070_seed.npz
	python code/scripts/run_linearized_scattering_on_filtered_soliton.py npz/1.150_seed.npz npz/1.070_1.500_lin_pot_12_13_14_15.npz -fi 1.650 -t pot 12 13 14 15 --title "Fig. 3, Pot. + (12) + (13) + (14) + (15)"

.PHONY: 1.150_1.500_lin.npz
1.150_1.500_lin.npz: 1.150_1.500_lin_pot_12.npz 1.150_1.500_lin_pot_12_13.npz 1.150_1.500_lin_pot_12_13_14.npz 1.150_1.500_lin_pot_12_13_14_15.npz

# # FIG4
# 1.085_ue_3.500_lin.npz: 1.085_ue_seed.npz
# 	python code/scripts/run_linearized_scattering_on_filtered_soliton.py npz/1.085_ue_seed.npz npz/1.085_ue_3.500_lin.npz -fi 3.500 -t 12 13


# Figures for the paper
# ===============================================

Fig1.pdf: 1.010_filt.npz
	python code/scripts/plot_timedomain_spectrum_rescon.py --vmin 1E-8 npz/1.010_filt.npz fig/Fig1.pdf

Fig2.pdf: 1.070_1.650.npz
	python code/scripts/plot_timedomain_spectrum_rescon.py --vmin 1E-7 npz/1.070_1.650.npz fig/Fig2.pdf

Fig3.pdf: 1.150_1.500.npz
	python code/scripts/plot_timedomain_spectrum_rescon.py --vmin 1E-7 npz/1.150_1.500.npz fig/Fig3.pdf

Fig4.pdf: 1.085_ue_3.500.npz
	python code/scripts/plot_timedomain_spectrum_rescon.py --vmin 1E-7 npz/1.085_ue_3.500.npz fig/Fig4.pdf

Fig5.pdf: 1.010_2.100_int.npz 1.010_2.100_int_sp.npz
	python code/scripts/plot_fig_5.py --vmin 5E-6 npz/1.010_2.100_int.npz npz/1.010_2.100_int_sp.npz fig/Fig5.pdf

Fig6.pdf: 1.010_1.100_int.npz 1.010_1.100_int_sp.npz
	python code/scripts/plot_fig_6.py --vmin 5E-6 npz/1.010_1.100_int.npz npz/1.010_1.100_int_sp.npz fig/Fig6.pdf

Fig7.pdf: 1.010_seed.npz $(addprefix 1.010_, $(addsuffix _rc10_sp.npz, $(RC_FREQS)))
	python code/scripts/plot_fig_7.py npz/1.010_1.100_int.npz $(addprefix npz/1.010_, $(addsuffix _rc10_sp.npz, $(RC_FREQS))) fig/Fig7.pdf

.PHONY: fig
fig: Fig1.pdf Fig2.pdf Fig3.pdf Fig4.pdf Fig5.pdf Fig6.pdf Fig7.pdf


# Text targets
# ============

.PHONY: draft
draft: fig
	cp fig/Fig*.pdf text/Figures/
	cd text && pdflatex Draft.tex
	cd text && bibtex Draft
	cd text && pdflatex Draft.tex
	cd text && pdflatex Draft.tex


# Paper source archive
# ====================
.PHONY: draft.zip
draft.zip:
	cd text && \
	zip draft.zip Draft.tex Bibliography.bib Figures/*.pdf && \
	mv draft.zip .. && \
	cd -

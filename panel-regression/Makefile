clean:
	rm ./data/*
	rm ./results/*
	rm ./figures/*
	rm ./tables/*

all:
	python3 07a_prepare_data_contemp.py
	python3 01a_prepare_data_distlags.py
	python3 01b_prepare_data_distlags_firstdiff.py
	python3 01c_prepare_data_distlags_secdiff.py
	python3 01d_prepare_data_distlags_bins_firstdiff_L10.py
	python3 01e_create_file_modelselection.py
	Rscript 02a_smooth-models_distlags_L05.R
	Rscript 02b_smooth-models_distlags_L10.R
	python3 02c_plot_figure_01.py
	Rscript 02c_smooth-models_distlags_L20.R
	Rscript 03a_smooth-models_distlags_L10_modelselection.R
	python3 03b_plot_variable_selection_figure_S12.py
	Rscript 04a_binned-models_distlags_L10.R
	python3 04b_plot_bins_figure_S14.py
	python3 05a_extract_dependent.py
	python3 05b_stationarity_figure_S11.py
	python3 05c_multicollinearity.py
	python3 05c_overfitting.py
	python3 05d_make_table_vif.py
	python3 06_plot_distlags_figure13.py
	Rscript 07b_smooth-models_contemp_mobility.R
	Rscript 07c_smooth-models_distlags_mortality.R
	Rscript 07d_smooth-models_distlags_mobility.R
	python3 08_ftest.py


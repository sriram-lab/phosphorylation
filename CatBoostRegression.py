# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 21:46:04 2025

@author: jbren

"""

import pandas as pd
from sklearn.model_selection import GroupShuffleSplit, GroupKFold
import numpy as np
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
import json
import os
import matplotlib.pyplot as plt

# --- 1. Global Constants and Configuration ---
# -- Input Config --
FILE_PATH = r"C:\Users\jbren\Documents\Documents\phosphomutations\tsui_foldx_ss.csv"
EXT_FILE_PATH = r"C:\Users\jbren\Documents\Documents\phosphomutations\phosphoddG_phos.csv"

RUN_ID = 99
TARGET_COL = 'Exp'
GROUP_COL = 'pdb'
FILTER_RESIDUE = None # Set to a character (e.g., 'Y') to filter, or None to disable.

FEATURES_TO_KEEP = [
	#'ResType',
	'SS', 'RelSASA',
	'AltPosition',
	'delta_total_energy',
	# 'delta_Backbone_Hbond',
	# 'delta_Sidechain_Hbond',
	# 'delta_Van_der_Waals',
	# 'delta_Electrostatics',
	# 'delta_Solvation_Polar',
	# 'delta_Solvation_Hydrophobic',
	# 'delta_Van_der_Waals_clashes',
	# 'delta_torsional_clash',
	# 'delta_backbone_clash',
	# 'delta_helix_dipole',
	# 'Hydrophobic-short', 
	# 'total-short',
	# 'Hydrophobic-short'
	# 'Polar-short','Positive-short','Negative-short','total-short'
	'Hydrophobic-long',
	'Polar-long',
    'total-long'
]


# -- Model & Tuning Config -
N_OPTUNA_TRIALS = 20
N_SPLITS_CV = 5
HOLDOUT_TEST_SIZE = 0.20
RANDOM_STATE = 42

# -- 1. Output Path Templates --
MODEL_SAVE_PATH_TPL = "{model_name}_model{run}.{extension}"
CONFIG_SAVE_PATH_TPL = "{model_name}_model_CONFIG{run}.json"
FEATURE_IMPORTANCE_PLOT_PATH_TPL = "{model_name}_feature_importance{run}.png"
FEATURE_IMPORTANCE_CSV_PATH_TPL = "{model_name}_feature_importance{run}.csv"
VALIDATION_PLOT_PATH_TPL = "{model_name}_validation_plot{run}.png"
EXTERNAL_VALIDATION_PLOT_PATH_TPL = "{model_name}_external_validation_plot{run}.png"

# --- Library Availability Checks ---
try:
	from catboost import CatBoostRegressor, Pool
	CATBOOST_AVAILABLE = True
except ImportError:
	CATBOOST_AVAILABLE = False

try:
	from xgboost import XGBRegressor
	XGBOOST_AVAILABLE = True
except ImportError:
	XGBOOST_AVAILABLE = False

try:
	import optuna
	OPTUNA_AVAILABLE = True
except ImportError:
	OPTUNA_AVAILABLE = False

try:
	import matplotlib.pyplot as plt
	PLOT_LIBS_AVAILABLE = True
except ImportError:
	PLOT_LIBS_AVAILABLE = False


# --- Helper Functions ---
def plot_feature_importance(importance_df, model_name, save_path):
	"""Generates and saves a bar plot for feature importances."""
	if not PLOT_LIBS_AVAILABLE: return
	try:
		df_plot = importance_df.head(20).sort_values(by='Importance', ascending=True)
		plt.figure(figsize=(10, 8))
		plt.barh(df_plot['Feature'], df_plot['Importance'], color='skyblue')
		plt.xlabel("Feature Importance")
		plt.ylabel("Feature")
		plt.title(f"Top 20 Feature Importances for {model_name}")
		plt.tight_layout()
		plt.savefig(save_path)
		plt.close()
		print(f"	Feature importance plot saved to {save_path}")
	except Exception as e:
		print(f"	Could not plot feature importance: {e}")

def generate_regression_plot(y_true, y_pred, model_name, plot_title_detail, rmse, pearson_r, plot_filename):
	"""Generates and saves a scatter plot for regression results."""
	if not PLOT_LIBS_AVAILABLE: return
	valid_indices = ~np.isnan(y_true) & ~np.isnan(y_pred)
	y_true_valid, y_pred_valid = y_true[valid_indices], y_pred[valid_indices]
	if len(y_true_valid) == 0: return
	plt.figure(figsize=(8, 6))
	plt.scatter(y_true_valid, y_pred_valid, alpha=0.5)
	min_val, max_val = min(np.nanmin(y_true_valid), np.nanmin(y_pred_valid)), max(np.nanmax(y_true_valid), np.nanmax(y_pred_valid))
	plt.plot([min_val, max_val], [min_val, max_val], 'k--', lw=2)
	plt.xlabel("Experimental 'Exp'")
	plt.ylabel("Predicted 'Exp'")
	plt.title(f"{model_name}: Predicted vs. Exp {plot_title_detail}\nRMSE: {rmse:.4f}, Pearson r: {pearson_r:.4f}")
	try:
		plt.savefig(plot_filename)
		plt.close()
		print(f"	Plot saved to {plot_filename}")
	except Exception as e:
		print(f"	Error saving plot {plot_filename}: {e}")


# *****************************************************************************
#  Core Functions
# *****************************************************************************
def preprocess_data(df):
	"""Applies all necessary preprocessing steps to the dataframe."""
	
	# Create 'res_type' from 'ResCode' if it exists
	if 'ResCode' in df.columns:
		df['res_type'] = df['ResCode'].astype(str).str[0]

	# Impute missing values
	for col in ['RelSASA']:
		if col in df.columns:
			# Note: Imputing with the median of the *current* dataset.
			# For a production system, these medians should be saved from the training set.
			df[col] = df[col].fillna(df[col].median())
			
	if 'SS' in df.columns:
		df['SS'] = df['SS'].fillna(df['SS'].mode()[0])
		# Map SS values to broader categories
		ss_mapping = {
			'-': 'C', 'T': 'C', 'S': 'C', # Coil types -> C
			'E': 'E', 'B': 'E',           # Beta-strand types -> E
			'H': 'H', 'G': 'H', 'I': 'H'  # Helix types -> H
		}
		df['SS'] = df['SS'].map(ss_mapping).fillna(df['SS'])
	
	return df

def load_and_split_data(file_path):
	"""Loads data, checks for features, preprocesses, and splits into train/validation sets."""
	print("--- Loading and Splitting Data ---")
	try:
		df = pd.read_csv(file_path)
		# Filter for a specific residue based on the global switch
		if FILTER_RESIDUE and 'ResCode' in df.columns:
			print(f"Filtering for ResCode starting with '{FILTER_RESIDUE}'...")
			df = df[df['ResCode'].str.startswith(FILTER_RESIDUE)]
	except FileNotFoundError:
		print(f"Error: Data file not found at {file_path}")
		return None, None

	missing_features = [col for col in FEATURES_TO_KEEP if col not in df.columns]
	if missing_features:
		print(f"\nFATAL ERROR: The following required columns are missing from the input file: {missing_features}")
		print("Please add these columns to your CSV file and try again.")
		return None, None
	print("All required features found in the dataframe.")
	
	df = preprocess_data(df.copy())

	gss = GroupShuffleSplit(n_splits=1, test_size=HOLDOUT_TEST_SIZE, random_state=RANDOM_STATE)
	try:
		train_tune_idx, final_val_idx = next(gss.split(df, groups=df[GROUP_COL]))
	except (ValueError, KeyError) as e:
		print(f"Error during GroupShuffleSplit: {e}")
		return None, None

	df_traintune = df.iloc[train_tune_idx].copy()
	df_final_validation = df.iloc[final_val_idx].copy()

	print(f"Train/Tune set size: {len(df_traintune)}, Validation set size: {len(df_final_validation)}")
	return df_traintune, df_final_validation

def train_model(df_traintune, model_name):
	"""Tunes, trains, and saves a single model (CatBoost or XGBoost)."""
	print(f"\n--- Training Model: {model_name} ---")

	base_feature_cols = [col for col in FEATURES_TO_KEEP if col in df_traintune.columns]
		
	X_traintune_base = df_traintune[base_feature_cols]
	y_traintune = df_traintune[TARGET_COL]
	groups = df_traintune[GROUP_COL]

	model_config = {'model_type': model_name.lower(), 'base_feature_cols': base_feature_cols}
	best_params = {}

	if model_name == 'CatBoost':
		if not CATBOOST_AVAILABLE: print("CatBoost not available."); return None, None
		cat_features = [col for col in ['SS', 'res_type'] if col in X_traintune_base.columns]
		cat_indices = [X_traintune_base.columns.get_loc(col) for col in cat_features]
		best_params = {'random_state': RANDOM_STATE, 'verbose': 0, 'cat_features': cat_indices}
		model_config['catboost_categorical_feature_names'] = cat_features
		base_model_class = CatBoostRegressor

	elif model_name == 'XGBoost':
		if not XGBOOST_AVAILABLE: print("XGBoost not available."); return None, None
		cols_to_ohe = [col for col in ['SS', 'res_type'] if col in X_traintune_base.columns]
		# We will handle one-hot encoding inside the Optuna objective
		best_params = {'random_state': RANDOM_STATE, 'objective': 'reg:squarederror'}
		model_config['cols_to_ohe_for_xgb'] = cols_to_ohe
		base_model_class = XGBRegressor
	else:
		raise ValueError(f"Unknown model_name: {model_name}")

	if OPTUNA_AVAILABLE:
		print("	Running Optuna for hyperparameter tuning...")

		def objective(trial):
			"""Optuna objective function to minimize RMSE using GroupKFold CV."""
			# Define Hyperparameter Search Space
			if model_name == 'CatBoost':
				params = {
					'iterations': trial.suggest_int('iterations', 100, 1000),
					'depth': trial.suggest_int('depth', 4, 10),
					'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
					'l2_leaf_reg': trial.suggest_float('l2_leaf_reg', 1e-3, 10.0, log=True),
					'border_count': trial.suggest_int('border_count', 32, 255),
				}
				params.update(best_params) # Add fixed params like random_state, verbose
			
			elif model_name == 'XGBoost':
				params = {
					'n_estimators': trial.suggest_int('n_estimators', 100, 1000),
					'max_depth': trial.suggest_int('max_depth', 3, 9),
					'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
					'subsample': trial.suggest_float('subsample', 0.6, 1.0),
					'colsample_bytree': trial.suggest_float('colsample_bytree', 0.6, 1.0),
					'gamma': trial.suggest_float('gamma', 1e-8, 1.0, log=True),
					'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 1.0, log=True),
					'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 1.0, log=True),
				}
				params.update(best_params)

			# --- Cross-validation with Out-of-Fold Predictions ---
			gkf = GroupKFold(n_splits=N_SPLITS_CV)
			oof_predictions = np.zeros(len(df_traintune))

			for fold, (train_idx, val_idx) in enumerate(gkf.split(X_traintune_base, y_traintune, groups)):
				X_train_fold, X_val_fold = X_traintune_base.iloc[train_idx], X_traintune_base.iloc[val_idx]
				y_train_fold, y_val_fold = y_traintune.iloc[train_idx], y_traintune.iloc[val_idx]

				if model_name == 'XGBoost':
					# Handle OHE for XGBoost inside the fold
					X_train_fold = pd.get_dummies(X_train_fold, columns=cols_to_ohe, prefix=cols_to_ohe, dummy_na=False)
					X_val_fold = pd.get_dummies(X_val_fold, columns=cols_to_ohe, prefix=cols_to_ohe, dummy_na=False)
					# Align columns after OHE
					train_cols = X_train_fold.columns
					val_cols = X_val_fold.columns
					missing_in_val = set(train_cols) - set(val_cols)
					for c in missing_in_val:
						X_val_fold[c] = 0
					missing_in_train = set(val_cols) - set(train_cols)
					for c in missing_in_train:
						X_train_fold[c] = 0
					X_val_fold = X_val_fold[train_cols]
				
				model = base_model_class(**params)
				
				if model_name == 'CatBoost':
					model.fit(X_train_fold, y_train_fold, eval_set=[(X_val_fold, y_val_fold)], early_stopping_rounds=50, verbose=False)
				else: # For XGBoost
					model.fit(X_train_fold, y_train_fold, eval_set=[(X_val_fold, y_val_fold)], verbose=False)

				oof_predictions[val_idx] = model.predict(X_val_fold)

			rmse = np.sqrt(mean_squared_error(y_traintune, oof_predictions))
			return rmse

		study = optuna.create_study(direction='minimize')
		study.optimize(objective, n_trials=N_OPTUNA_TRIALS)

		print(f"	Best trial for {model_name}:")
		print(f"	Value (RMSE): {study.best_value:.5f}")
		print("	Params: ")
		for key, value in study.best_params.items():
			print(f"		{key}: {value}")
		
		# Update best_params with the results from Optuna
		best_params.update(study.best_params)

	print("	Training definitive model on full traintune dataset...")
	# Handle final data prep for XGBoost
	if model_name == 'XGBoost':
		X_train = pd.get_dummies(X_traintune_base, columns=model_config['cols_to_ohe_for_xgb'], prefix=model_config['cols_to_ohe_for_xgb'], dummy_na=False)
		model_config['xgb_feature_columns_after_ohe'] = X_train.columns.tolist()
	else:
		X_train = X_traintune_base

	definitive_model = base_model_class(**best_params)
	definitive_model.fit(X_train, y_traintune)

	print("	Calculating feature importance...")
	importance_df = None
	if model_name == 'CatBoost':
		importance_df = pd.DataFrame({
			'Feature': X_train.columns,
			'Importance': definitive_model.get_feature_importance()
		}).sort_values(by='Importance', ascending=False)
	elif model_name == 'XGBoost':
		importance_df = pd.DataFrame({
			'Feature': X_train.columns,
			'Importance': definitive_model.feature_importances_
		}).sort_values(by='Importance', ascending=False)

	if importance_df is not None:
		plot_path = FEATURE_IMPORTANCE_PLOT_PATH_TPL.format(model_name=model_name.lower(), run=RUN_ID)
		csv_path = FEATURE_IMPORTANCE_CSV_PATH_TPL.format(model_name=model_name.lower(), run=RUN_ID)
		plot_feature_importance(importance_df, model_name, plot_path)
		importance_df.to_csv(csv_path, index=False)

	model_name_lower = model_name.lower()
	extension = "cbm" if model_name == 'CatBoost' else "json"
	model_save_path = MODEL_SAVE_PATH_TPL.format(model_name=model_name_lower, extension=extension, run=RUN_ID)
	config_save_path = CONFIG_SAVE_PATH_TPL.format(model_name=model_name_lower, run=RUN_ID)
	
	# Save final model parameters to config
	# Convert numpy types to native python types for JSON serialization
	for key, val in best_params.items():
		if isinstance(val, (np.int32, np.int64)):
			best_params[key] = int(val)
		if isinstance(val, (np.float32, np.float64)):
			best_params[key] = float(val)
	model_config['final_model_params'] = best_params


	definitive_model.save_model(model_save_path)
	with open(config_save_path, 'w') as f:
		json.dump(model_config, f, indent=4)
	print(f"	Model saved to: {model_save_path}")
	print(f"	Config saved to: {config_save_path}")

	return definitive_model, config_save_path


def evaluate_model(trained_model, df_validation, model_config, eval_type="Final Hold-Out"):
	"""Evaluates a trained model and returns metrics AND predictions."""
	model_name = "XGBoost" if "xgb" in model_config['model_type'] else "CatBoost"
	print(f"\n--- Evaluating Model: {model_name} on {eval_type} Set ---")

	# Use 'exp' for external, 'Exp' for internal, handle missing case
	target_col = TARGET_COL if TARGET_COL in df_validation.columns else 'exp'
	if target_col not in df_validation.columns:
		print(f"Warning: Target column ('{TARGET_COL}' or 'exp') not found. Cannot calculate metrics.")
		return {'model': model_name, 'error': 'Missing target col'}, None

	X_validation_base = df_validation[model_config['base_feature_cols']]
	y_validation = df_validation[target_col]

	# *** CORRECTED: Check if there are enough data points for metric calculation ***
	if len(y_validation) < 2:
		print(f"	Warning: Not enough data points ({len(y_validation)}) in the {eval_type} set to calculate metrics. Skipping evaluation.")
		return {'model': model_name, 'error': f'Not enough data in {eval_type} set'}, None

	if model_name == 'XGBoost':
		X_eval = pd.get_dummies(X_validation_base, columns=model_config['cols_to_ohe_for_xgb'], prefix=model_config['cols_to_ohe_for_xgb'])
		X_eval = X_eval.reindex(columns=model_config['xgb_feature_columns_after_ohe'], fill_value=0)
	else: # CatBoost
		X_eval = X_validation_base

	predictions = trained_model.predict(X_eval)
	rmse = np.sqrt(mean_squared_error(y_validation, predictions))
	pearson_r, _ = pearsonr(y_validation, predictions)

	print(f"	{eval_type} RMSE: {rmse:.4f}")
	print(f"	{eval_type} Pearson r: {pearson_r:.4f}")

	plot_tpl = EXTERNAL_VALIDATION_PLOT_PATH_TPL if eval_type == "External" else VALIDATION_PLOT_PATH_TPL
	plot_filename = plot_tpl.format(model_name=model_name.lower(), run=RUN_ID)
	
	generate_regression_plot(
		y_validation.values, predictions, model_name,
		f"({eval_type} Validation)", rmse, pearson_r,
		plot_filename
	)

	return {'model': model_name, f'{eval_type}_rmse': rmse, f'{eval_type}_pearson': pearson_r}, predictions


def run_training_workflow():
	"""Main orchestrator for training and internal validation."""
	df_traintune, df_validation = load_and_split_data(FILE_PATH)
	if df_traintune is None:
		print("Workflow aborted due to data loading/splitting error.")
		return None, []

	evaluation_results = []
	trained_models_info = [] # To store paths for external validation
	
	for model_name in ['CatBoost', 'XGBoost']:
		trained_model, config_path = train_model(df_traintune, model_name)

		if trained_model:
			with open(config_path, 'r') as f:
				config = json.load(f)

			results, predictions = evaluate_model(trained_model, df_validation, config, eval_type="Final Hold-Out")
			if predictions is not None:
				df_validation[f'{model_name}_predictions'] = predictions
				evaluation_results.append(results)
			
			trained_models_info.append({'name': model_name, 'config_path': config_path, 'model': trained_model})

	print("\n--- Final Workflow Summary ---")
	results_df = pd.DataFrame(evaluation_results)
	if not results_df.empty:
		print(results_df.to_string())
	else:
		print("No models were successfully trained or evaluated.")
	return df_validation, trained_models_info


def run_external_validation(external_data_path, trained_model, model_config):
	"""
	Orchestrates making predictions on an external dataset and then evaluating them.
	"""
	model_name = model_config['model_type'].capitalize()
	print(f"\n--- Running Full External Validation Workflow for {model_name} ---")

	try:
		ext_df = pd.read_csv(external_data_path)
		# Filter for a specific residue based on the global switch
		if FILTER_RESIDUE and 'ResCode' in ext_df.columns:
			print(f"Filtering external data for ResCode starting with '{FILTER_RESIDUE}'...")
			ext_df = ext_df[ext_df['ResCode'].str.startswith(FILTER_RESIDUE)]
	except FileNotFoundError:
		print(f"Error: External data file not found at {external_data_path}")
		return None

	ext_df_processed = preprocess_data(ext_df.copy())
	
	missing_features = [col for col in model_config['base_feature_cols'] if col not in ext_df_processed.columns]
	if missing_features:
		print(f"FATAL ERROR: Required columns missing from external data after processing: {missing_features}")
		return None

	results, predictions = evaluate_model(trained_model, ext_df_processed, model_config, eval_type="External")

	if predictions is not None:
		ext_df_processed[f'{model_name}_predictions'] = predictions
		return ext_df_processed
	return None

# --- Main Execution ---
if __name__ == "__main__":
	# 1. Run the training and internal validation workflow
	internal_results_df, trained_models = run_training_workflow()

	# 2. If successful, run the external validation for each trained model
	if trained_models:
		print("\n\n" + "="*50)
		print("STARTING EXTERNAL VALIDATION")
		print("="*50)
		
		all_external_results = {}
		for model_info in trained_models:
			model_name = model_info['name']
			config_path = model_info['config_path']
			model = model_info['model']
			
			with open(config_path, 'r') as f:
				config = json.load(f)

			external_df_with_preds = run_external_validation(
				external_data_path=EXT_FILE_PATH,
				trained_model=model,
				model_config=config
			)
			
			if external_df_with_preds is not None:
				print(f"External validation for {model_name} completed.")
				all_external_results[model_name] = external_df_with_preds.head()
			else:
				print(f"External validation for {model_name} failed.")
				
		# You can inspect `all_external_results` dictionary to see the head of each result
		# For example: print(all_external_results['CatBoost'])
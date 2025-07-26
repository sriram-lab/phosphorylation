import pandas as pd
from sklearn.metrics import mean_squared_error
import numpy as np
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_regression
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.neighbors import KNeighborsRegressor
import matplotlib.pyplot as plt
import xgboost as xgb

# Load data and preprocess
df = pd.read_csv('tsu_FoldX.csv')
df = df[df['Mut'] == 'Y'].dropna()

# List of feature names
train_features = [
    "Area", "Number of Residues", "EvoEF", "Foldx", "ss", "phi", "psi"]
# Create new columns for the negative values
# for feature in train_features:
#     df[f'neg_{feature}'] = -df[feature]

# # Update the train_features list to include the new columns
# negative_features = [f'neg_{feature}' for feature in train_features]
# train_features.extend(negative_features)

# Extract the features for the training set
X = df[train_features]
y = df['Exp']
scaler1 = StandardScaler()
X_scaled = scaler1.fit_transform(X)

# Feature selection
selector = SelectKBest(score_func=f_regression, k=15)
X_selected = selector.fit_transform(X_scaled, y)
selected_features = [train_features[i] for i in selector.get_support(indices=True)]

# Prepare cross-validation
kf = KFold(n_splits=5, shuffle=True, random_state=42)

# Cross-validation for GLM using sklearn's LinearRegression
glm = LinearRegression()
glm_predictions = cross_val_predict(glm, X_selected, y, cv=kf)

# Fit GLM model on the entire dataset to get coefficients
glm.fit(X_selected, y)

# Cross-validation for GBM
gbm = GradientBoostingRegressor(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=28)
gbm_predictions = cross_val_predict(gbm, X_selected, y, cv=kf)

# Fit GBM model on the entire dataset to get feature importances
gbm.fit(X_selected, y)

# Cross-validation for XGBoost
xgb_model = xgb.XGBRegressor(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=28)
xgb_predictions = cross_val_predict(xgb_model, X_selected, y, cv=kf)

# Fit XGBoost model on the entire dataset to get feature importances
xgb_model.fit(X_selected, y)

# Cross-validation for Nearest Neighbor Regressor using only "Area"
nnr = KNeighborsRegressor(n_neighbors=5)
X_area = df[["Area"]]
scaler2 = StandardScaler()
X_area_scaled = scaler2.fit_transform(X_area)
nnr_predictions = cross_val_predict(nnr, X_area_scaled, y, cv=kf)

# Fit Nearest Neighbor model on the entire dataset
nnr.fit(X_area_scaled, y)

# Evaluate models
def evaluate_model(true, predicted):
    mse = mean_squared_error(true, predicted)
    rmse = np.sqrt(mse)
    pearson_corr, _ = pearsonr(true, predicted)
    spearman_corr, _ = spearmanr(true, predicted)
    return mse, rmse, pearson_corr, spearman_corr

glm_mse, glm_rmse, glm_pearson, glm_spearman = evaluate_model(y, glm_predictions)
gbm_mse, gbm_rmse, gbm_pearson, gbm_spearman = evaluate_model(y, gbm_predictions)
xgb_mse, xgb_rmse, xgb_pearson, xgb_spearman = evaluate_model(y, xgb_predictions)
nnr_mse, nnr_rmse, nnr_pearson, nnr_spearman = evaluate_model(y, nnr_predictions)

print(f"GLM - MSE: {glm_mse}, RMSE: {glm_rmse}, Pearson: {glm_pearson}, Spearman: {glm_spearman}")
print(f"GBM - MSE: {gbm_mse}, RMSE: {gbm_rmse}, Pearson: {gbm_pearson}, Spearman: {gbm_spearman}")
print(f"XGB - MSE: {xgb_mse}, RMSE: {xgb_rmse}, Pearson: {xgb_pearson}, Spearman: {xgb_spearman}")
print(f"NNR - MSE: {nnr_mse}, RMSE: {nnr_rmse}, Pearson: {nnr_pearson}, Spearman: {nnr_spearman}")


# # Add predictions to DataFrame
# df['GLM_Predictions'] = glm_predictions
# df['GBM_Predictions'] = gbm_predictions
# df['XGB_Predictions'] = xgb_predictions
# df['NNR_Predictions'] = nnr_predictions

# # Save updated DataFrame to CSV
# df.to_csv('output_features_with_predictions_Y.csv', index=False)

# print("Updated DataFrame with predictions saved to 'output_features_with_predictions.csv'")

# Evaluation
# glm_mse = mean_squared_error(y, glm_predictions)
# gbm_mse = mean_squared_error(y, gbm_predictions)
# xgb_mse = mean_squared_error(y, xgb_predictions)
# nnr_mse = mean_squared_error(y, nnr_predictions)

# glm_pearson_corr, _ = pearsonr(y, glm_predictions)
# gbm_pearson_corr, _ = pearsonr(y, gbm_predictions)
# xgb_pearson_corr, _ = pearsonr(y, xgb_predictions)
# nnr_pearson_corr, _ = pearsonr(y, nnr_predictions)

# glm_spearman_corr, _ = spearmanr(y, glm_predictions)
# gbm_spearman_corr, _ = spearmanr(y, gbm_predictions)
# xgb_spearman_corr, _ = spearmanr(y, xgb_predictions)
# nnr_spearman_corr, _ = spearmanr(y, nnr_predictions)

# print(f'GLM Cross-validated MSE: {glm_mse:.4f}, Pearson Correlation: {glm_pearson_corr:.4f}, Spearman Correlation: {glm_spearman_corr:.4f}')
# print(f'GBM Cross-validated MSE: {gbm_mse:.4f}, Pearson Correlation: {gbm_pearson_corr:.4f}, Spearman Correlation: {gbm_spearman_corr:.4f}')
# print(f'XGB Cross-validated MSE: {xgb_mse:.4f}, Pearson Correlation: {xgb_pearson_corr:.4f}, Spearman Correlation: {xgb_spearman_corr:.4f}')
# print(f'NNR Cross-validated MSE: {nnr_mse:.4f}, Pearson Correlation: {nnr_pearson_corr:.4f}, Spearman Correlation: {nnr_spearman_corr:.4f}')

# Print feature importances for GBM
print("GBM Feature Importances:")
for feature, importance in zip(train_features, gbm.feature_importances_):
    print(f"{feature}: {importance}")

# Print feature importances for XGBoost
print("XGBoost Feature Importances:")
for feature, importance in zip(train_features, xgb_model.feature_importances_):
    print(f"{feature}: {importance}")

# Print coefficients for GLM
print("GLM Coefficients:")
for feature, coef in zip(['Intercept'] + train_features, np.append(glm.intercept_, glm.coef_)):
    print(f"{feature}: {coef}")

# # Predict on new dataset
# new_df = pd.read_csv('final_output_Y.csv')  # Replace with your actual new dataset file


# # Ensure new_df has a unique identifier column
# if 'ID' not in new_df.columns:
#     new_df['ID'] = range(len(new_df))

# new_df_predictions = new_df.copy().dropna(subset=train_features)

# # Preprocess the data
# X_predict = new_df_predictions[train_features]
# X_predict_scaled = scaler1.transform(X_predict)  # Ensure scaler can handle batch input

# # GLM Predictions
# new_df_predictions['GLM_Predictions'] = glm_new_predictions = glm.predict(X_predict_scaled)

# # GBM Predictions
# new_df_predictions['GBM_Predictions'] = gbm_new_predictions = gbm.predict(X_predict_scaled)

# # XGBoost Predictions
# new_df_predictions['XGB_Predictions'] = xgb_new_predictions = xgb_model.predict(X_predict_scaled)

# # Nearest Neighbor Predictions using only "Area"
# X_predict_area = new_df_predictions[["Area"]]
# X_predict_area_scaled = scaler2.transform(X_predict_area)
# new_df_predictions['NNR_Predictions'] = nnr.predict(X_predict_area_scaled)

# # Merge predictions back to the original new_df using the unique identifier
# new_df = new_df.merge(new_df_predictions[['ID', 'GLM_Predictions', 'GBM_Predictions', 'XGB_Predictions', 'NNR_Predictions']], 
#                       on='ID', how='left')

# # Save the final DataFrame to a CSV file
# new_df.to_csv('final_output_with_predictions.csv', index=False)
# print("Predictions have been merged and the final dataframe is saved to 'final_output_with_predictions.csv'.")

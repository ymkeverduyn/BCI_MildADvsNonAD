# -*- coding: utf-8 -*-
"""BCI_SVM_classifier.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1GCEqU--pLvAzvARUXcZRWw3oNAYwMLtl

## Setup
"""

# Commented out IPython magic to ensure Python compatibility.
#SETUP

COLAB = True
SKTIME_INSTALLED = False

if COLAB:
    from google.colab import drive
    drive.mount('/content/drive')
    # Load the contents of the directory
    !ls
    # Change your working directory to the folder where you stored your files, e.g.
#     %cd /content/drive/MyDrive/BCI project

import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV, cross_val_score, KFold, cross_validate, train_test_split
from sklearn.svm import SVC
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score
from sklearn import preprocessing
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt

"""## Train a SVM on the ERP data

### Load the ERP data
"""

# load from csv
erp_df = pd.read_csv('Kopie van merged_erp_features.csv')

# string labels into numerical categories
erp_df['Condition'] = erp_df['Condition'].map({'MCI': 1, 'NonAD': 0})
print(erp_df['Condition'].value_counts())
display(erp_df)

# exclude participant 9 from ERP data for better comparison
erp_df = erp_df.drop(erp_df[erp_df.Participant == 'sub-009'].index).reset_index(drop=True)

# change the structure of the dataframe so that each participant's data is in either test or train, not both
# and there are both classes in all folds for CV

# Identify the row indices for participants five and eight
participant_5_indices_erp = erp_df[erp_df['Participant'] == 'sub-005'].index
participant_8_indices_erp = erp_df[erp_df['Participant'] == 'sub-008'].index

# Swap the rows for participants five and eight for all columns in Auditory
for column in erp_df.columns:
    erp_df.loc[participant_5_indices_erp, column], erp_df.loc[participant_8_indices_erp, column] = erp_df.loc[participant_8_indices_erp, column].values, erp_df.loc[participant_5_indices_erp, column].values

erp_df = erp_df.reset_index(drop=True)
erp_df.shape

# split into target and features
y_ERP = erp_df['Condition']
X_ERP = erp_df.iloc[:, 1:-3]

ERP_cols = X_ERP.columns

# scale features
x_ERP = X_ERP.values
min_max_scaler = preprocessing.MinMaxScaler()
x_ERP_scaled = min_max_scaler.fit_transform(x_ERP)
X_ERP = pd.DataFrame(x_ERP_scaled, columns=ERP_cols)

print(X_ERP.shape, y_ERP.shape)
display(y_ERP)
display(X_ERP)

# put together features and target
joined_df = X_ERP.join(y_ERP)

display(joined_df)

# export as csv
joined_df.to_csv('ERP_ordered_scaled_labelled.csv')

"""### Train the model with 5-fold nested cross-validation
Find good hyperparameters and test over different folds
"""

# Define the parameter grid
param_grid = {
    'C': [0.1, 1, 10, 100],
    'gamma': [1, 0.1, 0.01, 0.001],
    'kernel': ['rbf', 'linear']
}

# Define the scoring metrics
scoring = {
    'accuracy': make_scorer(accuracy_score),
    'precision': make_scorer(precision_score, zero_division=0),
    'recall': make_scorer(recall_score, zero_division=0),
    'f1': make_scorer(f1_score, zero_division=0)
}

# Initialize the SVM model
svm = SVC()

# Set up GridSearchCV with 5-fold cross-validation
cv = KFold(n_splits=5)
grid_search = GridSearchCV(svm, param_grid, cv=cv, scoring=scoring, refit='accuracy', return_train_score=True)

# Perform the grid search on the data
grid_search.fit(X_ERP, y_ERP)

# Get the best parameters and the corresponding scores
best_C = grid_search.best_params_['C']
best_gamma = grid_search.best_params_['gamma']
best_kernel = grid_search.best_params_['kernel']
best_score = grid_search.best_score_

print(f"Best C ERP: {best_C}")
print(f"Best gamma ERP: {best_gamma}")
print(f"Best kernel ERP: {best_kernel}")
print(f"Best cross-validated accuracy ERP: {best_score}")

# Extract the scores for the best model
best_model = grid_search.best_estimator_
scores = cross_validate(best_model, X_ERP, y_ERP, cv=cv, scoring=scoring, return_train_score=True)

print("Metrics with optimized hyperparameters for ERP:")
print('Test Accuracy:', scores['test_accuracy'], 'Mean:', scores['test_accuracy'].mean())
print('Train Accuracy:', scores['train_accuracy'], 'Mean:', scores['train_accuracy'].mean())
print('Test Precision:', scores['test_precision'], 'Mean:', scores['test_precision'].mean())
print('Test Recall:', scores['test_recall'], 'Mean:', scores['test_recall'].mean())
print('Test F1:', scores['test_f1'], 'Mean:', scores['test_f1'].mean())

# train and test a new model with the parameters found by grid search
X_train, X_test, y_train, y_test = train_test_split(X_ERP, y_ERP, test_size=0.4, shuffle=False)

# create a SVM classifier according to the best hyperparameters found above
svm = SVC(C=0.1, gamma=0.01, kernel='rbf')
svm.fit(X_train, y_train)

# Predict on training and testing data
y_train_pred = svm.predict(X_train)
y_test_pred = svm.predict(X_test)

# Evaluate performance
train_accuracy = accuracy_score(y_train, y_train_pred)
test_accuracy = accuracy_score(y_test, y_test_pred)
test_precision = precision_score(y_test, y_test_pred, average='binary')
test_recall = recall_score(y_test, y_test_pred, average='binary')
test_f1 = f1_score(y_test, y_test_pred, average='binary')

# Print results
print("SVM on ERP data final model:")
print(f"Training Accuracy: {train_accuracy:.4f}")
print(f"Testing Accuracy: {test_accuracy:.4f}")
print(f"Testing Precision: {test_precision:.4f}")
print(f"Testing Recall: {test_recall:.4f}")
print(f"Testing F1 Score: {test_f1:.4f}")

"""### SVM ERP final model confusion matrix"""

# Generate confusion matrices - 1 = AD; 0 = non-AD
cm = confusion_matrix(y_test, y_test_pred)

# Plot confusion matrices
plt.figure(figsize=(12, 5))

# SVM Confusion Matrix
disp_svm = ConfusionMatrixDisplay(confusion_matrix=cm)
disp_svm.plot(cmap=plt.cm.Blues)
plt.title("SVM Confusion Matrix \n (Auditory Dataset)")
plt.show()

"""## Train a SVM on the Resting State Data

### Load the Resting State Data
"""

# load from csv
resting_df = pd.read_csv('Resting State Data /resting_state_complete.csv')

# change the structure of the dataframe so that each participant's data is in either test or train, not both
# and there are both classes in all folds for CV

# Identify the row indices for participants five and eight
participant_5_indices_rest = resting_df[resting_df['Participant'] == 'participant005'].index
participant_8_indices_rest = resting_df[resting_df['Participant'] == 'participant008'].index

# Swap the rows for participants five and eight for all columns in Auditory
for column in resting_df.columns:
    resting_df.loc[participant_5_indices_rest, column], resting_df.loc[participant_8_indices_rest, column] = resting_df.loc[participant_8_indices_rest, column].values, resting_df.loc[participant_5_indices_rest, column].values

resting_df = resting_df.reset_index(drop=True)
resting_df.shape

# string labels into numerical categories
resting_df['Participant'] = resting_df['Participant'].map({'participant001': 0, 'participant002': 1, 'participant003': 0, 'participant004': 1, 'participant005': 1,
                                                           'participant007': 1, 'participant008': 0, 'participant010': 0, 'participant011': 1, 'participant012': 0,})

resting_df = resting_df.rename(columns={'Participant' : 'label'})

print(resting_df.label.value_counts())
display(resting_df)

# split into target and features
y_resting = resting_df['label']
X_resting = resting_df.iloc[:, 1:-1]

resting_cols = X_resting.columns

# scale features
x_resting = X_resting.values
min_max_scaler = preprocessing.MinMaxScaler()
x_resting_scaled = min_max_scaler.fit_transform(x_resting)
X_resting = pd.DataFrame(x_resting_scaled, columns=resting_cols)

print(y_resting.shape, X_resting.shape)
display(X_resting)

"""### Train the model with 5-fold nested cross-validation
Find good hyperparameters and test over differnt folds
"""

# Define the parameter grid
param_grid = {
    'C': [0.1, 1, 10, 100],
    'gamma': [1, 0.1, 0.01, 0.001],
    'kernel': ['rbf', 'linear']
}

# Define the scoring metrics
scoring = {
    'accuracy': make_scorer(accuracy_score),
    'precision': make_scorer(precision_score, zero_division=0),
    'recall': make_scorer(recall_score, zero_division=0),
    'f1': make_scorer(f1_score, zero_division=0)
}

# Initialize the SVM model
svm = SVC()

# Set up GridSearchCV with 5-fold cross-validation
cv = KFold(n_splits=5)
grid_search = GridSearchCV(svm, param_grid, cv=cv, scoring=scoring, refit='accuracy', return_train_score=True)

# Perform the grid search on the data
grid_search.fit(X_resting, y_resting)

# Get the best parameters and the corresponding scores
best_C = grid_search.best_params_['C']
best_gamma = grid_search.best_params_['gamma']
best_kernel = grid_search.best_params_['kernel']
best_score = grid_search.best_score_

print(f"Best C Resting: {best_C}")
print(f"Best gamma Resting: {best_gamma}")
print(f"Best kernel Resting: {best_kernel}")
print(f"Best cross-validated accuracy Resting: {best_score}")

# Extract the scores for the best model
best_model = grid_search.best_estimator_
scores = cross_validate(best_model, X_resting, y_resting, cv=cv, scoring=scoring, return_train_score=True)

print("Metrics with optimized hyperparameters for Resting:")
print('Test Accuracy:', scores['test_accuracy'], 'Mean:', scores['test_accuracy'].mean())
print('Train Accuracy:', scores['train_accuracy'], 'Mean:', scores['train_accuracy'].mean())
print('Test Precision:', scores['test_precision'], 'Mean:', scores['test_precision'].mean())
print('Test Recall:', scores['test_recall'], 'Mean:', scores['test_recall'].mean())
print('Test F1:', scores['test_f1'], 'Mean:', scores['test_f1'].mean())

"""### Perform Feature selection

The classifier is obviously overfitting to the training data in each fold -- we need to do some sort of feature selection
"""

# select the best k features using ANOVA
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif

# loop over some different values for k

k_values = np.arange(1, 51)
test_acc = []
train_acc = []

for k in k_values:

  # define feature selection
  fs = SelectKBest(score_func=f_classif, k=k)
  # apply feature selection
  X_selected = fs.fit_transform(X_resting, y_resting)
  #print(X_selected.shape)

  # train the model with only the k best features & perform 5-fold cross validation

  # Define the parameter grid
  param_grid = {
      'C': [0.1, 1, 10, 100],
      'gamma': [1, 0.1, 0.01, 0.001],
      'kernel': ['rbf', 'linear']
  }

  # Define the scoring metrics
  scoring = {
      'accuracy': make_scorer(accuracy_score),
      'precision': make_scorer(precision_score, zero_division=0),
      'recall': make_scorer(recall_score, zero_division=0),
      'f1': make_scorer(f1_score, zero_division=0)
  }

  # Initialize the SVM model
  svm = SVC()

  # Set up GridSearchCV with 5-fold cross-validation
  cv = KFold(n_splits=5)
  grid_search = GridSearchCV(svm, param_grid, cv=cv, scoring=scoring, refit='accuracy', return_train_score=True)

  # Perform the grid search on the data
  grid_search.fit(X_selected, y_resting)

  # Get the best parameters and the corresponding scores
  best_C = grid_search.best_params_['C']
  best_gamma = grid_search.best_params_['gamma']
  best_kernel = grid_search.best_params_['kernel']
  best_score = grid_search.best_score_

  #print(f"USING THE BEST k = {k} FEATURES:")

  #print(f"Best C Resting: {best_C}")
  #print(f"Best gamma Resting: {best_gamma}")
  #print(f"Best kernel Resting: {best_kernel}")
  #print(f"Best cross-validated accuracy Resting: {best_score}")

  # Extract the scores for the best model
  best_model = grid_search.best_estimator_
  scores = cross_validate(best_model, X_selected, y_resting, cv=cv, scoring=scoring, return_train_score=True)

  #print("Metrics with optimized hyperparameters for Resting:")
  #print('Test Accuracy:', scores['test_accuracy'], 'Mean:', scores['test_accuracy'].mean())
  #print('Train Accuracy:', scores['train_accuracy'], 'Mean:', scores['train_accuracy'].mean())
  #print('Test Precision:', scores['test_precision'], 'Mean:', scores['test_precision'].mean())
  #print('Test Recall:', scores['test_recall'], 'Mean:', scores['test_recall'].mean())
  #print('Test F1:', scores['test_f1'], 'Mean:', scores['test_f1'].mean())

  test_acc.append(scores['test_accuracy'].mean())
  train_acc.append(scores['train_accuracy'].mean())

  #print(test_acc, train_acc)

import matplotlib.pyplot as plt

# Plotting the data
plt.figure(figsize=(10, 6))
plt.plot(k_values, test_acc, label='Test Accuracy', marker='o')
plt.plot(k_values, train_acc, label='Train Accuracy', marker='o')
plt.xlabel('Value for k')
plt.ylabel('Mean Accuracy')
plt.title('SVM Resting: \n Mean Accuracy Over 5 Folds Using k Best Features')
plt.legend()
plt.grid(True)
plt.show()

"""Now that we found the best k for SelectKBest, create a new dataframe with only those columns and run the algorithm on it again to get all performance measures."""

# get the list with the 40 best features

fs = SelectKBest(score_func=f_classif, k=40)
fs.fit_transform(X_resting, y_resting)

mask = fs.get_support()
new_features = [] # list of the K best features

for bool_val, feature in zip(mask, X_resting.columns):
    if bool_val:
        new_features.append(feature)

print(new_features)

X_resting_selected = X_resting[new_features]
display(X_resting_selected)

# export to csv
X_resting_selected.to_csv('SVM_X_resting_selected.csv')

# Define the parameter grid
param_grid = {
    'C': [0.1, 1, 10, 100],
    'gamma': [1, 0.1, 0.01, 0.001],
    'kernel': ['rbf', 'linear']
}

# Define the scoring metrics
scoring = {
    'accuracy': make_scorer(accuracy_score),
    'precision': make_scorer(precision_score, zero_division=0),
    'recall': make_scorer(recall_score, zero_division=0),
    'f1': make_scorer(f1_score, zero_division=0)
}

# Initialize the SVM model
svm = SVC()

# Set up GridSearchCV with 5-fold cross-validation
cv = KFold(n_splits=5)
grid_search = GridSearchCV(svm, param_grid, cv=cv, scoring=scoring, refit='accuracy', return_train_score=True)

# Perform the grid search on the data
grid_search.fit(X_resting_selected, y_resting)

# Get the best parameters and the corresponding scores
best_C = grid_search.best_params_['C']
best_gamma = grid_search.best_params_['gamma']
best_kernel = grid_search.best_params_['kernel']
best_score = grid_search.best_score_

print(f"Best C Resting: {best_C}")
print(f"Best gamma Resting: {best_gamma}")
print(f"Best kernel Resting: {best_kernel}")
print(f"Best cross-validated accuracy Resting: {best_score}")

# Extract the scores for the best model
best_model = grid_search.best_estimator_
scores = cross_validate(best_model, X_resting_selected, y_resting, cv=cv, scoring=scoring, return_train_score=True)

print("Metrics with optimized hyperparameters for Resting selected:")
print('Test Accuracy:', scores['test_accuracy'], 'Mean:', scores['test_accuracy'].mean())
print('Train Accuracy:', scores['train_accuracy'], 'Mean:', scores['train_accuracy'].mean())
print('Test Precision:', scores['test_precision'], 'Mean:', scores['test_precision'].mean())
print('Test Recall:', scores['test_recall'], 'Mean:', scores['test_recall'].mean())
print('Test F1:', scores['test_f1'], 'Mean:', scores['test_f1'].mean())

# train and test a new model with the parameters found by grid search
X_train, X_test, y_train, y_test = train_test_split(X_resting_selected, y_resting, test_size=0.4, shuffle=False)

# create a SVM classifier according to the best hyperparameters found above
svm = SVC(C=100, gamma=1, kernel='linear')
svm.fit(X_train, y_train)

# Predict on training and testing data
y_train_pred = svm.predict(X_train)
y_test_pred = svm.predict(X_test)

# Evaluate performance
train_accuracy = accuracy_score(y_train, y_train_pred)
test_accuracy = accuracy_score(y_test, y_test_pred)
test_precision = precision_score(y_test, y_test_pred, average='binary')
test_recall = recall_score(y_test, y_test_pred, average='binary')
test_f1 = f1_score(y_test, y_test_pred, average='binary')

# Print results
print("SVM on resting data final model:")
print(f"Training Accuracy: {train_accuracy:.4f}")
print(f"Testing Accuracy: {test_accuracy:.4f}")
print(f"Testing Precision: {test_precision:.4f}")
print(f"Testing Recall: {test_recall:.4f}")
print(f"Testing F1 Score: {test_f1:.4f}")

"""### SVM final model resting state confusion matrix"""

# Generate confusion matrices - 1 = AD; 0 = non-AD
cm = confusion_matrix(y_test, y_test_pred)

# Plot confusion matrices
plt.figure(figsize=(12, 5))

# SVM Confusion Matrix
disp_svm = ConfusionMatrixDisplay(confusion_matrix=cm)
disp_svm.plot(cmap=plt.cm.Blues)
plt.title("SVM Confusion Matrix \n (Resting-State Dataset)")
plt.show()
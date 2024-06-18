# -*- coding: utf-8 -*-
"""CombinedFeatures.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1n0dBSGBas7l5ePEj3tX2-PmemXQZvNxf

### Setup
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
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score
from sklearn import preprocessing

"""### Load ERP and Resting state data"""

# combine the ERP data and the resting state date into one df
# load ERP from csv
erp_df = pd.read_csv('Kopie van merged_erp_features.csv')

# string labels into numerical categories
erp_df['Condition'] = erp_df['Condition'].map({'MCI': 1, 'NonAD': 0})
print(erp_df['Condition'].value_counts())
display(erp_df)


# load Resting from csv
resting_df = pd.read_csv('Resting State Data /resting_state_complete.csv')

# string labels into numerical categories
resting_df['label'] = resting_df['Participant'].map({'participant001': 0, 'participant002': 1, 'participant003': 0, 'participant004': 1, 'participant005': 1,
                                                           'participant007': 1, 'participant008': 0, 'participant010': 0, 'participant011': 1, 'participant012': 0,})

print(resting_df.label.value_counts())
display(resting_df)

# remove participant 009 from the ERP dataframe to match the resting state df

erp_df = erp_df.drop(erp_df[erp_df.Participant == 'sub-009'].index).reset_index()
print(erp_df['Participant'].value_counts())
print(resting_df['Participant'].value_counts())
print(erp_df.shape)
display(erp_df)

"""### Merge the Dataframes"""

# get rid of columns we dont need
resting_df = resting_df.drop(columns=['Unnamed: 0', 'Participant'])
print(resting_df.shape)
erp_df = erp_df.drop(columns=['index', 'Epoch', 'Channel', 'Condition'])
print(erp_df.shape)

# merge together
full_df = erp_df.join(resting_df)
display(full_df)

# change the structure of the dataframe so that each participant's data is in either test or train, not both
# and there are both classes in all folds for CV

# Identify the row indices for participants five and eight
participant_5_indices = full_df[full_df['Participant'] == 'sub-005'].index
participant_8_indices = full_df[full_df['Participant'] == 'sub-008'].index

# Swap the rows for participants five and eight for all columns in Auditory
for column in full_df.columns:
    full_df.loc[participant_5_indices, column], full_df.loc[participant_8_indices, column] = full_df.loc[participant_8_indices, column].values, full_df.loc[participant_5_indices, column].values

full_df = full_df.reset_index(drop=True)
full_df.shape

# drop participant column
full_df = full_df.drop(columns=['Participant'])
full_df.shape

print(full_df.label.value_counts())

"""### Train SVM and KNN on merged dataset

Use nested 5-fold cross validation to include hyperparameter optimization
"""

# split into target and features
y = full_df['label']
X = full_df.iloc[:, :-1]

# scale features
x = X.values
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
X = pd.DataFrame(x_scaled)

print(y.shape, X.shape)
display(X)

"""#### SVM"""

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
grid_search.fit(X, y)

# Get the best parameters and the corresponding scores
best_C = grid_search.best_params_['C']
best_gamma = grid_search.best_params_['gamma']
best_kernel = grid_search.best_params_['kernel']
best_score = grid_search.best_score_

print(f"Best C: {best_C}")
print(f"Best gamma: {best_gamma}")
print(f"Best kernel: {best_kernel}")
print(f"Best cross-validated accuracy: {best_score}")

# Extract the scores for the best model
best_model = grid_search.best_estimator_
scores = cross_validate(best_model, X, y, cv=cv, scoring=scoring, return_train_score=True)

print("Metrics with optimized hyperparameters for the SVM:")
print('Test Accuracy:', scores['test_accuracy'], 'Mean:', scores['test_accuracy'].mean())
print('Train Accuracy:', scores['train_accuracy'], 'Mean:', scores['train_accuracy'].mean())
print('Test Precision:', scores['test_precision'], 'Mean:', scores['test_precision'].mean())
print('Test Recall:', scores['test_recall'], 'Mean:', scores['test_recall'].mean())
print('Test F1:', scores['test_f1'], 'Mean:', scores['test_f1'].mean())

"""#### KNN"""

# Define the range of k values to test
param_grid = {'n_neighbors': np.arange(1, 21)}

# Define the scoring metrics with zero_division parameter
scoring = {
    'accuracy': make_scorer(accuracy_score),
    'precision': make_scorer(precision_score, zero_division=0),
    'recall': make_scorer(recall_score, zero_division=0),
    'f1': make_scorer(f1_score, zero_division=0)
}

# Create the KNN classifier
knn = KNeighborsClassifier()

# Set up GridSearchCV with 5-fold cross-validation
cv = KFold(n_splits=5)
grid_search = GridSearchCV(knn, param_grid, cv=cv, scoring=scoring, refit='accuracy', return_train_score=True)

# Perform the grid search on the auditory data
grid_search.fit(X, y)

# Get the best parameters and the corresponding scores
best_k = grid_search.best_params_['n_neighbors']
best_score = grid_search.best_score_

print(f"Best k: {best_k}")
print(f"Best cross-validated accuracy: {best_score}")

# Extract the scores for the best model
best_model = grid_search.best_estimator_
scores = cross_validate(best_model, X, y, cv=cv, scoring=scoring, return_train_score=True)

print("Metrics with optimized k for the KNN:")
print('Test Accuracy:', scores['test_accuracy'], 'Mean:', scores['test_accuracy'].mean())
print('Train Accuracy:', scores['train_accuracy'], 'Mean:', scores['train_accuracy'].mean())
print('Test Precision:', scores['test_precision'], 'Mean:', scores['test_precision'].mean())
print('Test Recall:', scores['test_recall'], 'Mean:', scores['test_recall'].mean())
print('Test F1:', scores['test_f1'], 'Mean:', scores['test_f1'].mean())

"""## Use only selected features on resting and merge with ERP data afterwards

Perform feature selection, since both classifiers are overfitting to the training data.
"""

# read in the ERP data from csv with data scaled and participants 5 & 8 swapped + labels
ERP_df = pd.read_csv('ERP_ordered_scaled_labelled.csv').iloc[:, 1:]
display(ERP_df)

"""### SVM on merged with feature selection"""

# read in csv with selected features, already scaled and with participants 5 & 8 swapped for correctness
SVM_resting_df = pd.read_csv('SVM_X_resting_selected.csv').iloc[:, 1:]
display(SVM_resting_df)

# merge selected features with ERP data

SVM_full = SVM_resting_df.join(ERP_df)
display(SVM_full)

# split into features and target

X_SVM = SVM_full.drop(columns=['Condition'])
y_SVM = SVM_full['Condition']

print(X_SVM.shape, y_SVM.shape)

# perform 5-fold CV and train the model
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
grid_search.fit(X_SVM, y_SVM)

# Get the best parameters and the corresponding scores
best_C = grid_search.best_params_['C']
best_gamma = grid_search.best_params_['gamma']
best_kernel = grid_search.best_params_['kernel']
best_score = grid_search.best_score_

print(f"Best C: {best_C}")
print(f"Best gamma: {best_gamma}")
print(f"Best kernel: {best_kernel}")
print(f"Best cross-validated accuracy: {best_score}")

# Extract the scores for the best model
best_model = grid_search.best_estimator_
scores = cross_validate(best_model, X_SVM, y_SVM, cv=cv, scoring=scoring, return_train_score=True)

print("Metrics with optimized hyperparameters for the SVM on selected features:")
print('Test Accuracy:', scores['test_accuracy'], 'Mean:', scores['test_accuracy'].mean())
print('Train Accuracy:', scores['train_accuracy'], 'Mean:', scores['train_accuracy'].mean())
print('Test Precision:', scores['test_precision'], 'Mean:', scores['test_precision'].mean())
print('Test Recall:', scores['test_recall'], 'Mean:', scores['test_recall'].mean())
print('Test F1:', scores['test_f1'], 'Mean:', scores['test_f1'].mean())

# train and test a new model with the parameters found by grid search
X_train, X_test, y_train, y_test = train_test_split(X_SVM, y_SVM, test_size=0.4, shuffle=False)

# create a SVM classifier according to the best hyperparameters found above
svm = SVC(C=100, gamma=0.1, kernel='rbf')
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
print("SVM on merged data final model:")
print(f"Training Accuracy: {train_accuracy:.4f}")
print(f"Testing Accuracy: {test_accuracy:.4f}")
print(f"Testing Precision: {test_precision:.4f}")
print(f"Testing Recall: {test_recall:.4f}")
print(f"Testing F1 Score: {test_f1:.4f}")

"""### KNN on merged with feature selection"""

# read in csv with selected features, already scaled and with participants 5 & 8 swapped for correctness
KNN_resting_df = pd.read_csv('KNN_X_resting_selected.csv').iloc[:, 1:]
display(KNN_resting_df)

# merge selected features with ERP data

KNN_full = KNN_resting_df.join(ERP_df)
display(KNN_full)

# split into features and target

X_KNN = KNN_full.drop(columns=['Condition'])
y_KNN = KNN_full['Condition']

print(X_KNN.shape, y_SVM.shape)

# perform 5-fold CV and train the model

# Define the range of k values to test
param_grid = {'n_neighbors': np.arange(1, 21)}

# Define the scoring metrics with zero_division parameter
scoring = {
    'accuracy': make_scorer(accuracy_score),
    'precision': make_scorer(precision_score, zero_division=0),
    'recall': make_scorer(recall_score, zero_division=0),
    'f1': make_scorer(f1_score, zero_division=0)
}

# Create the KNN classifier
knn = KNeighborsClassifier()

# Set up GridSearchCV with 5-fold cross-validation
cv = KFold(n_splits=5)
grid_search = GridSearchCV(knn, param_grid, cv=cv, scoring=scoring, refit='accuracy', return_train_score=True)

# Perform the grid search on the auditory data
grid_search.fit(X_KNN, y_KNN)

# Get the best parameters and the corresponding scores
best_k = grid_search.best_params_['n_neighbors']
best_score = grid_search.best_score_

print(f"Best k: {best_k}")
print(f"Best cross-validated accuracy: {best_score}")

# Extract the scores for the best model
best_model = grid_search.best_estimator_
scores = cross_validate(best_model, X_KNN, y_KNN, cv=cv, scoring=scoring, return_train_score=True)

print("Metrics with optimized k for the KNN:")
print('Test Accuracy:', scores['test_accuracy'], 'Mean:', scores['test_accuracy'].mean())
print('Train Accuracy:', scores['train_accuracy'], 'Mean:', scores['train_accuracy'].mean())
print('Test Precision:', scores['test_precision'], 'Mean:', scores['test_precision'].mean())
print('Test Recall:', scores['test_recall'], 'Mean:', scores['test_recall'].mean())
print('Test F1:', scores['test_f1'], 'Mean:', scores['test_f1'].mean())

# train and test a new model with the parameters found by grid search
X_train, X_test, y_train, y_test = train_test_split(X_KNN, y_KNN, test_size=0.4, shuffle=False)

# create a KNN classifier according to the best hyperparameters found above
knn = KNeighborsClassifier(n_neighbors=15)
knn.fit(X_train, y_train)

# Predict on training and testing data
y_train_pred = knn.predict(X_train)
y_test_pred = knn.predict(X_test)

# Evaluate performance
train_accuracy = accuracy_score(y_train, y_train_pred)
test_accuracy = accuracy_score(y_test, y_test_pred)
test_precision = precision_score(y_test, y_test_pred, average='binary')
test_recall = recall_score(y_test, y_test_pred, average='binary')
test_f1 = f1_score(y_test, y_test_pred, average='binary')

# Print results
print("KNN on merged data final model:")
print(f"Training Accuracy: {train_accuracy:.4f}")
print(f"Testing Accuracy: {test_accuracy:.4f}")
print(f"Testing Precision: {test_precision:.4f}")
print(f"Testing Recall: {test_recall:.4f}")
print(f"Testing F1 Score: {test_f1:.4f}")
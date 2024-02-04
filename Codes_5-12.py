### Section 6 CODE 5 ###

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score, roc_auc_score

### Section 6 CODE 6 ###

# Load a marker labeled dataset from step 5
data = pd.read_csv("ExampleMarker1_Processed.tsv", sep='\t')

# Split the dataset into features and target
X = data.drop("LabelName", axis=1)
y = data["LabelName"]

### Section 6 CODE 7 ###

# Split the data into training and testing sets (80% train, 20% test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


### Section 6 CODE 8 ###

# Create and train a Random Forest classifier
clf = RandomForestClassifier(max_depth=4, random_state=42)
clf.fit(X_train, y_train)

### Section 6 CODE 9 ###

# Predict on the test set and calculate evaluation metrics to check effectiveness
y_pred = clf.predict(X_test)
precision = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred)
accuracy = accuracy_score(y_test, y_pred)
auc = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])

### Section 6 CODE 10 ###

# Load a unlabeled sample dataset from step 4
Sample1_Labeled = pd.read_csv("Sample1_Processed.tsv", sep='\t')

# Ensure the features in the unlabeled dataset match the training labeled set
# This includes both the number and the order of features
unlabeled_features = Sample1_Labeled[X.columns]

### Section 6 CODE 11 ###

# Use the trained classifier to make predictions
unlabeled_predictions = clf.predict(unlabeled_features)

# Add these predictions a column to the dataframe and save or process further
Sample1_Labeled['ExampleMarker1_Positive'] = unlabeled_predictions

### Section 6 CODE 12 ###

# You can now save this dataframe or use it for further analysis once all
# random forest classifiers are added to label unlabeled datasets from step 4
Sample1_Labeled.to_csv("Sample1_Labeled.tsv", sep='\t', index=False)

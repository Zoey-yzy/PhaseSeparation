import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, confusion_matrix, classification_report, roc_curve, auc
import matplotlib.pyplot as plt

# Step 1: Load the dataset
# Assuming you have your dataset in a CSV file called 'dataset.csv'
# Make sure to replace 'dataset.csv' with the actual path to your dataset

# df = pd.read_csv("/Users/zoey/Desktop/PS/output/293tAll_2.csv")
df = pd.read_csv("/Users/zoey/Desktop/PS/output/EmbedAll.csv")
# df = df.drop(columns=['Unnamed: 0', 'Names', 'AA_Number', 'Sequences'])

# Step 2: Prepare the data
# Columns 0-5 are features, column 6 is the target
# X = np.concatenate((df.iloc[:, 0:6].values, df.iloc[:, 8:12].values), axis=1)
X = df.iloc[:, 1:1025].values
y = df.iloc[:, 1026].values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
"""
# Split the data into train and test sets
X_class0 = X[y == 0]
X_class1 = X[y == 1]
y_class1 = y[y == 1]
X_0 = pd.DataFrame(X_class0).sample(n=len(X_class1)*2)
X = np.vstack((X_0, X_class1))
y = np.concatenate((np.zeros(len(X_0)),y_class1))
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

train_indices = np.arange(X_train.shape[0])
np.random.shuffle(train_indices)
X_train = X_train[train_indices]
y_train = y_train[train_indices]

test_indices = np.arange(X_test.shape[0])
np.random.shuffle(test_indices)
X_test = X_test[test_indices]
y_test = y_test[test_indices]
"""
# Create the DMatrix for XGBoost
dtrain = xgb.DMatrix(X_train, label=y_train)
dtest = xgb.DMatrix(X_test, label=y_test)

# Step 3: Set up the parameters for the XGBoost model (use the default para. but with binary obj. and log loss)
params = {
    'objective': 'binary:logistic',  # Binary classification
    'eval_metric': 'logloss',  # Evaluation metric
    'max_depth': 6,  # Maximum depth of a tree
    'eta': 0.3,  # Learning rate
    'seed': 0  # Random seed
}

# Train the model
num_rounds = 10
model = xgb.train(params, dtrain, num_rounds)

# Step 4: Make predictions on the test set
y_pred_prob = model.predict(dtest)
y_pred = np.round(y_pred_prob)

# Evaluate the model
accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred)
conf_matrix = confusion_matrix(y_test, y_pred)
class_report = classification_report(y_test, y_pred)

# ROC AUC plot
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
# import scikitplot as skplt
fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
roc_auc = auc(fpr, tpr)

plt.figure()
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:0.2f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Performance on PS/non-PS Proteins with 12 Parameters')
plt.legend(loc="lower right")
plt.show()
##
print("Accuracy:", accuracy)
print("Confusion Matrix:")
print(conf_matrix)
print("Classification Report:")
print(class_report)
## plot the tree
num_trees = model.best_iteration + 1
print("Number of trees:", num_trees)

# Plot all trees
for i in range(num_trees):
    plt.figure(figsize=(20, 10))
    xgb.plot_tree(model, num_trees=i)
    plt.title(f"Tree {i}")
    plt.show()

name_list = [x for x in df.columns.tolist() if x != 'DeepCoil' and x != 'PS']
importance_weight = model.get_score(importance_type='weight')
importance_gain = model.get_score(importance_type='gain')
importance_cover = model.get_score(importance_type='cover')

importance_df_weight = pd.DataFrame(list(importance_weight.items()), columns=['Feature', 'Weight'])
importance_df_gain = pd.DataFrame(list(importance_gain.items()), columns=['Feature', 'Gain'])
importance_df_cover = pd.DataFrame(list(importance_cover.items()), columns=['Feature', 'Cover'])

importance_df_weight['Feature'] = name_list
importance_df_gain['Feature'] = name_list
importance_df_cover['Feature'] = name_list

# Sort DataFrames by importance scores
importance_df_weight = importance_df_weight.sort_values(by='Weight', ascending=False)
importance_df_gain = importance_df_gain.sort_values(by='Gain', ascending=False)
importance_df_cover = importance_df_cover.sort_values(by='Cover', ascending=False)

# Plot feature importance (Weight)
importance_df_weight.plot(kind='bar', x='Feature', y='Weight', legend=False)
plt.title('Feature Importance by Weight')
plt.ylabel('Weight')
plt.show()

# Plot feature importance (Gain)
importance_df_gain.plot(kind='bar', x='Feature', y='Gain', legend=False)
plt.title('Feature Importance by Gain')
plt.ylabel('Gain')
plt.show()

# Plot feature importance (Cover)
importance_df_cover.plot(kind='bar', x='Feature', y='Cover', legend=False)
plt.title('Feature Importance by Cover')
plt.ylabel('Cover')
plt.show()
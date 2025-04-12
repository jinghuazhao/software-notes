# Statistcal learning

## Deep learning

### Autoencoder
* Four types, [https://github.com/nathanhubens/Autoencoders](https://github.com/nathanhubens/Autoencoders).
* [The notMNIST dataset](http://yann.lecun.com/exdb/mnist/), [https://www.datacamp.com/community/tutorials/autoencoder-keras-tutorial](https://www.datacamp.com/community/tutorials/autoencoder-keras-tutorial).
* [Credit fraud example](https://www.kaggle.com/mlg-ulb/creditcardfraud), [https://blogs.rstudio.com/tensorflow/posts/2018-01-24-keras-fraud-autoencoder/](https://blogs.rstudio.com/tensorflow/posts/2018-01-24-keras-fraud-autoencoder/).
* [PCA vs autoencoder](https://www.r-bloggers.com/pca-vs-autoencoders-for-dimensionality-reduction/)

## Machine learning

(as it is called in computer science literature)

* Apress, [https://github.com/apress](https://github.com/apress)
* Packt, [https://github.com/PacktPublishing](https://github.com/PacktPublishing)

## Statistical learning and [ISLR](https://CRAN.R-project.org/package=ISLR)

> Hastie T, Tibshirani R, Friedman J (2009). [The Elements of Statistical Learning: Data Mining, Inference, and Prediction](https://www.springer.com/gb/book/9780387848570). Springer.

> James G, Witten D, Hastie T, Tibshirani R (2013). [An Introduction to Statistical Learning with Applications in R](http://www-bcf.usc.edu/~gareth/ISL/). Springer.

> James G, Witten D, Hastie T, Tibshirani R (2021). [An Introduction to Statistical Learning with Applications in R, 2e](https://www.statlearning.com/), Springer.

> Chollet F, Allaire JJ (2017). [Deep Learning with R](https://livebook.manning.com/book/deep-learning-with-r/), Manning. [Source code](https://www.manning.com/books/deep-learning-with-r), [GitHub](https://github.com/jjallaire/deep-learning-with-r-notebooks).

## Python Code Examples

From: Machine Learning Fundamentals Handbook – Key Concepts, Algorithms, and Python Code Examples, 

Web: <https://www.freecodecamp.org/news/machine-learning-handbook/>

Changes are made a couple of places, e.g., Bagging, where option `base_estimator` is replaced with `estimator` as in
`BaggingRegressor(estimator=base_estimator, n_estimators=10, random_state=42)`, AdaBoost, the missing line is added:
`from sklearn.ensemble import GradientBoostingRegressor`. Now it is also possible to convert the MarkDown document
into a Jupyter notebook (.ipynb).

### Linear Regression

```python
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# Sample Data
stamps_bought = np.array([1, 3, 5, 7, 9]).reshape((-1, 1))  # Reshaping to make it a 2D array
amount_spent = np.array([2, 6, 8, 12, 18])

# Creating a Linear Regression Model
model = LinearRegression()

# Training the Model
model.fit(stamps_bought, amount_spent)

# Predictions
next_month_stamps = 10
predicted_spend = model.predict([[next_month_stamps]])

# Plotting
plt.scatter(stamps_bought, amount_spent, color='blue')
plt.plot(stamps_bought, model.predict(stamps_bought), color='red')
plt.title('Stamps Bought vs Amount Spent')
plt.xlabel('Stamps Bought')
plt.ylabel('Amount Spent ($)')
plt.grid(True)
plt.show()

# Displaying Prediction
print(f"If Alex buys {next_month_stamps} stamps next month, they will likely spend ${predicted_spend[0]:.2f}.")
```

### Logistic Regression

```python
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score

# Sample Data
pages = np.array([100, 150, 200, 250, 300, 350, 400, 450, 500]).reshape(-1, 1)
likes = np.array([0, 1, 1, 1, 0, 0, 0, 0, 0])  # 1: Like, 0: Dislike

# Creating a Logistic Regression Model
model = LogisticRegression()

# Training the Model
model.fit(pages, likes)

# Predictions
predict_book_pages = 260
predicted_like = model.predict([[predict_book_pages]])

# Plotting
plt.scatter(pages, likes, color='forestgreen')
plt.plot(pages, model.predict_proba(pages)[:, 1], color='darkred')
plt.title('Book Pages vs Like/Dislike')
plt.xlabel('Number of Pages')
plt.ylabel('Likelihood of Liking')
plt.axvline(x=predict_book_pages, color='green', linestyle='--')
plt.axhline(y=0.5, color='grey', linestyle='--')
plt.show()

# Displaying Prediction
print(f"Jenny will {'like' if predicted_like[0] == 1 else 'not like'} a book of {predict_book_pages} pages.")
```

### Linear Discrimant Analysis

```python
import numpy as np
import matplotlib.pyplot as plt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

# Sample Data
# [size, sweetness]
fruits_features = np.array([[3, 7], [2, 8], [3, 6], [4, 7], [1, 4], [2, 3], [3, 2], [4, 3]])
fruits_likes = np.array([1, 1, 1, 1, 0, 0, 0, 0])  # 1: Like, 0: Dislike

# Creating an LDA Model
model = LinearDiscriminantAnalysis()

# Training the Model
model.fit(fruits_features, fruits_likes)

# Prediction
new_fruit = np.array([[2.5, 6]])  # [size, sweetness]
predicted_like = model.predict(new_fruit)

# Plotting
plt.scatter(fruits_features[:, 0], fruits_features[:, 1], c=fruits_likes, cmap='viridis', marker='o')
plt.scatter(new_fruit[:, 0], new_fruit[:, 1], color='darkred', marker='x')
plt.title('Fruits Enjoyment Based on Size and Sweetness')
plt.xlabel('Size')
plt.ylabel('Sweetness')
plt.show()

# Displaying Prediction
print(f"Sarah will {'like' if predicted_like[0] == 1 else 'not like'} a fruit of size {new_fruit[0, 0]} and sweetness {new_fruit[0, 1]}.")
```

### Naive Bayes

```python
import numpy as np
import matplotlib.pyplot as plt
from sklearn.naive_bayes import GaussianNB

# Sample Data
# [movie_length, genre_code] (assuming genre is coded as: 0 for Action, 1 for Romance, etc.)
movies_features = np.array([[120, 0], [150, 1], [90, 0], [140, 1], [100, 0], [80, 1], [110, 0], [130, 1]])
movies_likes = np.array([1, 1, 0, 1, 0, 1, 0, 1])  # 1: Like, 0: Dislike

# Creating a Naive Bayes Model
model = GaussianNB()

# Training the Model
model.fit(movies_features, movies_likes)

# Prediction
new_movie = np.array([[100, 1]])  # [movie_length, genre_code]
predicted_like = model.predict(new_movie)

# Plotting
plt.scatter(movies_features[:, 0], movies_features[:, 1], c=movies_likes, cmap='viridis', marker='o')
plt.scatter(new_movie[:, 0], new_movie[:, 1], color='darkred', marker='x')
plt.title('Movie Likes Based on Length and Genre')
plt.xlabel('Movie Length (min)')
plt.ylabel('Genre Code')
plt.show()

# Displaying Prediction
print(f"Tom will {'like' if predicted_like[0] == 1 else 'not like'} a {new_movie[0, 0]}-min long movie of genre code {new_movie[0, 1]}.")
```

### Decision Trees

```python
import numpy as np
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeRegressor, plot_tree

# Sample Data
# [hours_studied]
study_hours = np.array([1, 2, 3, 4, 5, 6, 7, 8]).reshape(-1, 1)
test_scores = np.array([50, 55, 70, 80, 85, 90, 92, 98])

# Creating a Decision Tree Regression Model
model = DecisionTreeRegressor(max_depth=3)

# Training the Model
model.fit(study_hours, test_scores)

# Prediction
new_study_hour = np.array([[5.5]])  # example of hours studied
predicted_score = model.predict(new_study_hour)

# Plotting the Decision Tree
plt.figure(figsize=(12, 8))
plot_tree(model, filled=True, rounded=True, feature_names=["Study Hours"])
plt.title('Decision Tree Regressor Tree')
plt.show()

# Plotting Study Hours vs. Test Scores
plt.scatter(study_hours, test_scores, color='darkred')
plt.plot(np.sort(study_hours, axis=0), model.predict(np.sort(study_hours, axis=0)), color='orange')
plt.scatter(new_study_hour, predicted_score, color='green')
plt.title('Study Hours vs Test Scores')
plt.xlabel('Study Hours')
plt.ylabel('Test Scores')
plt.grid(True)
plt.show()

# Displaying Prediction
print(f"Predicted test score for {new_study_hour[0, 0]} hours of study: {predicted_score[0]:.2f}.")
```

### Bagging

```python
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import BaggingRegressor
from sklearn.tree import DecisionTreeRegressor, plot_tree  # Ensure plot_tree is imported
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

# Sample Data
clients_data = np.array([[2000, 60], [2500, 45], [1800, 75], [2200, 50], [2100, 62], [2300, 70], [1900, 55], [2000, 65]])
weight_loss = np.array([3, 2, 4, 3, 3.5, 4.5, 3.7, 4.2])

# Train-Test Split
X_train, X_test, y_train, y_test = train_test_split(clients_data, weight_loss, test_size=0.25, random_state=42)

# Creating a Bagging Model
base_estimator = DecisionTreeRegressor(max_depth=4)
model = BaggingRegressor(estimator=base_estimator, n_estimators=10, random_state=42)

# Training the Model
model.fit(X_train, y_train)

# Prediction & Evaluation
y_pred = model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)

# Displaying Prediction and Evaluation
print(f"True weight loss: {y_test}")
print(f"Predicted weight loss: {y_pred}")
print(f"Mean Squared Error: {mse:.2f}")

# Visualizing One of the Base Estimators (if desired)
plt.figure(figsize=(12, 8))
tree = model.estimators_[0]
plt.title('One of the Base Decision Trees from Bagging')
plot_tree(tree, filled=True, rounded=True, feature_names=["Calorie Intake", "Workout Duration"])
plt.show()
```

### Random Forest

```python
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report

# Expanded Data
plants_features = np.array([
    [3, 1], [2, 2], [4, 1], [3, 2], [5, 1], [2, 2], [4, 1], [5, 2],
    [3, 1], [4, 2], [5, 1], [3, 2], [2, 1], [4, 2], [3, 1], [4, 2],
    [5, 1], [2, 2], [3, 1], [4, 2], [2, 1], [5, 2], [3, 1], [4, 2]
])
plants_species = np.array([
    0, 1, 0, 1, 0, 1, 0, 1,
    0, 1, 0, 1, 0, 1, 0, 1,
    0, 1, 0, 1, 0, 1, 0, 1
])

# Train-Test Split
X_train, X_test, y_train, y_test = train_test_split(plants_features, plants_species, test_size=0.25, random_state=42)

# Creating a Random Forest Model
model = RandomForestClassifier(n_estimators=10, random_state=42)

# Training the Model
model.fit(X_train, y_train)

# Prediction & Evaluation
y_pred = model.predict(X_test)
classification_rep = classification_report(y_test, y_pred)

# Displaying Prediction and Evaluation
print("Classification Report:")
print(classification_rep)

# Scatter Plot Visualizing Classes
plt.figure(figsize=(8, 4))
for species, marker, color in zip([0, 1], ['o', 's'], ['forestgreen', 'darkred']):
    plt.scatter(plants_features[plants_species == species, 0],
                plants_features[plants_species == species, 1],
                marker=marker, color=color, label=f'Species {species}')

plt.xlabel('Leaf Size')
plt.ylabel('Flower Color (coded)')
plt.title('Scatter Plot of Species')
plt.legend()
plt.tight_layout()
plt.show()

# Visualizing Feature Importances
plt.figure(figsize=(8, 4))
features_importance = model.feature_importances_
features = ["Leaf Size", "Flower Color"]
plt.barh(features, features_importance, color = "darkred")
plt.xlabel('Importance')
plt.ylabel('Feature')
plt.title('Feature Importance')
plt.show()
```

### Boosting

#### AdaBoost

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error

# Seed for reproducibility
np.random.seed(42)

# Generate synthetic data
num_samples = 200
num_rooms = np.random.randint(3, 10, num_samples)
house_age = np.random.randint(1, 100, num_samples)
noise = np.random.normal(0, 50, num_samples)

# Assume a linear relation with price = 50*rooms + 0.5*age + noise
price = 50*num_rooms + 0.5*house_age + noise

# Create DataFrame
data = pd.DataFrame({'num_rooms': num_rooms, 'house_age': house_age, 'price': price})

# Plot
plt.scatter(data['num_rooms'], data['price'], label='Num Rooms vs Price', color = 'forestgreen')
plt.scatter(data['house_age'], data['price'], label='House Age vs Price', color = 'darkred')
plt.xlabel('Feature Value')
plt.ylabel('Price')
plt.legend()
plt.title('Scatter Plots of Features vs Price')
plt.show()

# Splitting data into training and testing sets
X = data[['num_rooms', 'house_age']]
y = data['price']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Initialize and train AdaBoost Regressor model
model_ab = AdaBoostRegressor(n_estimators=100, random_state=42)
model_ab.fit(X_train, y_train)

# Predictions
predictions = model_ab.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, predictions)
rmse = np.sqrt(mse)
print(f"Mean Squared Error: {mse:.2f}")
print(f"Root Mean Squared Error: {rmse:.2f}")

# Visualization: Actual vs Predicted Prices
plt.scatter(y_test, predictions, color = 'darkred')
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=3)
plt.xlabel('Actual Prices')
plt.ylabel('Predicted Prices')
plt.title('Actual vs Predicted House Prices with AdaBoost')
plt.show()
```

#### Gradient Boosting Model

```python
# Initialize and train Gradient Boosting Regressor model
model_gbm = GradientBoostingRegressor(n_estimators=100, learning_rate=0.1, max_depth=1, random_state=42)
model_gbm.fit(X_train, y_train)

# Predictions
predictions = model_gbm.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, predictions)
rmse = np.sqrt(mse)
print(f"Mean Squared Error: {mse:.2f}")
print(f"Root Mean Squared Error: {rmse:.2f}")

# Visualization: Actual vs Predicted Prices
plt.scatter(y_test, predictions, color = 'orange')
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=3)
plt.xlabel('Actual Prices')
plt.ylabel('Predicted Prices')
plt.title('Actual vs Predicted House Prices with GBM')
plt.show()
```

#### Extreme Gradient Boosting (XGBoost)

```python
import xgboost as xgb

# Initialize and train XGBoost model
model_xgb = xgb.XGBRegressor(objective ='reg:squarederror', n_estimators = 100, seed = 42)
model_xgb.fit(X_train, y_train)

# Predictions
predictions = model_xgb.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, predictions)
rmse = np.sqrt(mse)
print(f"Mean Squared Error: {mse:.2f}")
print(f"Root Mean Squared Error: {rmse:.2f}")

# Visualization: Actual vs Predicted Prices
plt.scatter(y_test, predictions, color="forestgreen")
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=3)
plt.xlabel('Actual Prices')
plt.ylabel('Predicted Prices')
plt.title('Actual vs Predicted House Prices with XGBoost')
plt.show()
```

## R code example

### Confusion matrix

```r
library(caret)

# Example data
set.seed(123)
data <- data.frame(
  Actual = sample(c("True", "False"), 100, replace = TRUE),
  Prediction = sample(c("True", "False"), 100, replace = TRUE)
)

# Create confusion matrix
cm <- confusionMatrix(as.factor(data$Prediction), as.factor(data$Actual), positive = "True")
print(cm)
```

# StatFit App

The app is designed to gamify the user experience of using fitness bracelets. It recommends completing different achievements, based on the daily user's activity indicators: wasted calories, steps, and sleep minutes per day. As a result, one can be a part of the leaderboard among other users. Machine Learning algorithms are the driving force behind the recommendation engines in the app.

##  Website Technology Stack

The following tools were used to develop the Website for the app:
- FastAPI
- Uvicorn (server)
- Sqlalchemy
- Postgres

## Machine Learning Technology Stack

Both classification and clusterization tasks were performed as a part of Machine Learning stack under the Python tools. All the supplementary materials can be found at `ml_part` folder (`Stage 1`, `Stage 2`, `Stage 3`, and `Stage 4` Python notebooks).

The clusterization task was conducted in order to distinguish users into different categories, taking into account their daily activity indicators (wasted calories, steps, and sleep minutes). The task itself includes such steps as defining the most appropriate number of clusters, choosing the best clusterization algorithm judging by the set metric, and adjusting the hyperparameters. The choice was made between the following algorithms:
- KMeans
- Agglomerative Clustering algorithm
- Gaussian Mixture model
- DBSCAN

The classification task is aimed at giving a prediction of the user's "fitness group" and, after that, offering an appropriate achievement. Likewise, the classification task includes choosing the best classification algorithm judging by the set metric and adjusting the hyperparameters. The choice was made between the following algorithms:
- KNN
- Decision Tree
- Naive Bayes Classifier
- XGBoost

## How to start the app?

Follow the instructions below:
```
git clone https://github.com/nogirosta/StatFit.git

python -m venv env      #create a virtual environment
.\env\Scripts\activate  #activate your virtual environment

pip install -r .\requirements.txt
uvicorn main:app --reload     #start server
visit  127.0.0.1:8000/
```

## Contact information

If you have any questions or comments about this project, please feel free to contact the authors:

- Rostyslav Fedorchenko
- Anton Tsybulnyk
- Danylo Shmorhun
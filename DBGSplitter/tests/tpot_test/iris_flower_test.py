from tpot import TPOTClassifier
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
import numpy as np

iris = load_iris()
X_train, X_test, y_train, y_test = train_test_split(iris.data.astype(np.float64),
    iris.target.astype(np.float64), train_size=0.75, test_size=0.25)

tpot = TPOTClassifier(generations=5, population_size=20, verbosity=2)
tpot.fit(X_train, y_train)

with open("evaluated_pipelines.detailed", "w") as detailedFile:
	for key,value in tpot.evaluated_individuals_.items():
		detailedFile.write("%s\n%s\n%s\n%s\n\n\n"%(key, 100*"-", str(value), 100*"-"))

with open("evaluated_pipelines.simple", "w") as simpleFile:
	for pipeline, stats in tpot.evaluated_individuals_.items():
		simpleFile.write("%s\t%s\n"%(pipeline, str(stats["internal_cv_score"])))

with open("best_pipeline", "w") as bestFile:
	bestFile.write(tpot.clean_pipeline_string(tpot._optimized_pipeline)+"\n")
	bestFile.write(str(tpot.fitted_pipeline_)+"\n")

with open("best_pipeline_score_test", "w") as bestFileScore:
	bestFileScore.write(str(tpot.score(X_test, y_test)))


tpot.export('tpot_iris_pipeline.py')

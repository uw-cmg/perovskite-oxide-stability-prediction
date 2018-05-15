import numpy as np
import pandas as pd
import sys
from sklearn.kernel_ridge import KernelRidge
from sklearn.ensemble import ExtraTreesClassifier
from feature_vector import gen_test, gen_train
from sklearn.feature_selection import VarianceThreshold
from sklearn import preprocessing

def read_data(cont_name, disc_name):
	cdata = pd.read_csv(cont_name, index_col=0)
	ddata = pd.read_csv(disc_name, index_col=0)
	Xc = cdata[cdata.columns[0:-2]]
	Xd = ddata[ddata.columns[0:-2]]
	y = cdata[cdata.columns[-2]]
	yl = y.copy()
	y2 = cdata[cdata.columns[-1]]
	# label the instance
	# Energy above 40 meV is considered to be unstable
	y[y <= 40] = 1
	y[y > 40] = 0
	return Xc, Xd, y, yl, y2

def data_process(Xc, Xd, Xc_model, Xd_model):
	p = 0.00
	sel = VarianceThreshold()
	sel.fit(Xc_model)
	Xc1 = Xc.loc[:, sel.variances_ > p * (1 - p)]
	sel.fit(Xd_model)
	Xd1 = Xd.loc[:, sel.variances_ > p * (1 - p)]
	# print('Removed {} feature from {} continuous features'.format(Xc.shape[1] - Xc1.shape[1], Xc.shape[1]))
	# print('Removed {} feature from {} discrete features'.format(Xd.shape[1] - Xd1.shape[1], Xd.shape[1]))
	X_total = pd.concat([Xc1, Xd1], axis=1)
	return X_total

def data_scale(X_total, X_total_model):
	scaler = preprocessing.StandardScaler().fit(X_total_model)
	feature_names = list(X_total)
	X_scaled = scaler.transform(X_total)
	return pd.DataFrame(X_scaled,columns=feature_names)

def select_features(index_file, select_n, X_total):
	indices_data = pd.read_csv(index_file, names=['order'])
	indices = np.array(indices_data['order'].tolist())
	selected = indices[:select_n]
	X_features = X_total.ix[:, selected]
	return X_features

def classification(X_scale, X_scale_test, y):
	clf = ExtraTreesClassifier(criterion='entropy', bootstrap=False, max_leaf_nodes=None,
	                           min_impurity_split=0.1, max_features=43, class_weight='balanced',
	                           min_samples_split=5, min_samples_leaf=1, max_depth=18, n_estimators=115)
	X_features1 = select_features('RFE_clf_indices.txt', 70, X_scale)
	X_features_test1 = select_features('RFE_clf_indices.txt', 70, X_scale_test)

	clf.fit(X_features1, y)
	stability_predict = clf.predict(X_features_test1)
	clf_result = pd.DataFrame(stability_predict, columns=['predicted stability'])
	return clf_result

def cut_highEs(X_features, yl, ye):
	# remove outliers
	X_s = X_features.loc[ye < 400]
	yl_s = yl[ye < 400]
	return X_s, yl_s

def reg_EaH(X_scale, X_scale_test, ye):
	reg = KernelRidge(kernel='rbf', alpha=0.007, gamma=0.007)
	X_features2 = select_features('RFE_eah_indices.txt', 70, X_scale)
	X_features_test2 = select_features('RFE_eah_indices.txt', 70, X_scale_test)

	X_s, ye_s = cut_highEs(X_features2, ye, ye)
	reg.fit(X_s, ye_s)
	y_predict = reg.predict(X_features_test2)
	EaH_predict = pd.DataFrame(y_predict, columns=['predicted Energy above hull'])
	return EaH_predict


def reg_FE(X_scale, X_scale_test, yf, ye):
	reg = KernelRidge(kernel='rbf', alpha=0.00464, gamma=0.0215)
	X_features3 = select_features('stability_fe_indices.txt', 20, X_scale)
	X_features_test3 = select_features('stability_fe_indices.txt', 20, X_scale_test)
	X_s, yf_s = cut_highEs(X_features3, yf, ye)
	reg.fit(X_s, yf_s)
	y_predict = reg.predict(X_features_test3)
	FE_predict = pd.DataFrame(y_predict, columns=['predicted Formation Energy'])
	return FE_predict

def write_result(testfile, output, clf_result, EaH_predict, FE_predict):
	test_data = pd.read_excel(testfile)
	raw_composition = test_data[['Material Composition', 'A site #1', 'A site #2',
	                             'A site #3', 'B site #1', 'B site #2', 'B site #3',
	                             'X site', 'Number of elements']]
	result = pd.concat([raw_composition, clf_result, EaH_predict, FE_predict], axis=1)
	result.to_excel(output, index=None)

def wrap_data(trainfile, testfile, id=0):

	gen_train(trainfile, id)
	gen_test(testfile, id)
	ctrain = 'c_{}_train.csv'.format(id)
	ctest = 'c_{}_test.csv'.format(id)
	dtrain = 'd_{}_train.csv'.format(id)
	dtest = 'd_{}_test.csv'.format(id)

	Xc, Xd, y, ye, yf = read_data(ctrain, dtrain)
	X_total = data_process(Xc, Xd, Xc, Xd)
	X_scale = data_scale(X_total, X_total)

	Xc_test, Xd_test, y_test, ye_test, yf_test = read_data(ctest, dtest)
	X_total_test = data_process(Xc_test, Xd_test, Xc, Xd)
	X_scale_test = data_scale(X_total_test, X_total)

	ye = ye.reset_index()['EnergyAboveHull']
	y = y.reset_index()['EnergyAboveHull']
	yf = yf.reset_index()['Formation_energy']

	return X_scale, X_scale_test, y, ye, yf


if __name__ == "__main__":

	trainfile = 'perovskite_DFT_EaH_FormE.xlsx' if len(sys.argv)<=1 else sys.argv[1]
	testfile = 'newCompound.xlsx' if len(sys.argv)<=2 else sys.argv[2]
	id = 0 if len(sys.argv)<=3 else sys.argv[3]

	X_scale, X_scale_test, y, ye, yf = wrap_data(trainfile, testfile, id)
	clf_result = classification(X_scale, X_scale_test, y)
	EaH_predict = reg_EaH(X_scale, X_scale_test, ye)
	FE_predict = reg_FE(X_scale, X_scale_test, yf, ye)
	output = 'prediction_result.xlsx'
	write_result(testfile, output, clf_result, EaH_predict, FE_predict)







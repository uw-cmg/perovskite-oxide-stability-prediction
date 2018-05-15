import numpy as np
import pandas as pd
from xlrd import open_workbook
import csv
import sys

def buildDict(attribute_index, refsheet, ifnulldefault):
	attribute = dict()
	attribute['feature_name'] = refsheet.cell(0, attribute_index).value
	for row_index in range(1, refsheet.nrows):
		element_name = refsheet.cell(row_index, 0).value
		data = refsheet.cell(row_index, attribute_index).value
		if (data == '' or data == ' '):
			attribute[element_name] = ifnulldefault
		else:
			try:
				attribute[element_name] = float(data)
			except:
				print()
				print('error in reading property')
				print(data, row_index, attribute_index)
	return attribute

def radii_shannon(shannon_name):
	dict_radii = pd.read_excel('shannon_perovskite.xlsx')
	A_r = dict(zip(dict_radii['Asite'],dict_radii['Aradii']))
	B_r = dict(zip(dict_radii['Bsite'],dict_radii['Bradii']))
	X_r = dict(zip(dict_radii['Xsite'],dict_radii['Xradii']))
	return A_r, B_r, X_r

A_r, B_r, X_r = radii_shannon('shannon_perovskite.xlsx')

def getComp(row_index, sheetname, col_index, compositionname, dicList, add_shannon=False):
	sites = []
	# each site have at most 3 element
	for index in range(3):
		if (col_index + index >= 8):
			# break when it meets the last column of elements
			break
		if (sheetname.cell(row_index, col_index + index).value != ''):
			name = sheetname.cell(row_index, col_index + index).value
			num = findNum(compositionname, name)
			properties = [name, num]
			if add_shannon:
				if (col_index == 1):
					try:
						properties.append(A_r[name])
					except KeyError:
						print('{} not in A site shannon dictionary'.format(name))
				elif (col_index == 4):
					try:
						properties.append(B_r[name])
					except KeyError:
						print('{} not in B site shannon dictionary'.format(name))
				else:
					try:
						properties.append(X_r[name])
					except KeyError:
						print('{} not in X site shannon dictionary'.format(name))
			for attrDic in dicList:
				properties.append(attrDic[name])
			sites.append(properties)
	return sites


# find the number of element in the composition given element name
def findNum(composition, element):
	if element == '':
		return 0
	pos = composition.find(element) + len(element)
	num = 0
	while (pos < len(composition) and composition[pos].isdigit()):
		num = num * 10 + int(composition[pos])
		pos += 1
	return num

def add_dheader(d_nums, discList):
	dheader = np.array(range(d_nums), dtype=object)
	dheader[0] = ''
	dheader[1] = 'numberofelements'
	site_order = ['Asite', 'Bsite']
	index = 2
	statics_order = ['weighted_avg']
	for statics in statics_order:
		for site in site_order:
			for disc in discList:
				dheader[index] = '{}_{}_{}'.format(site, disc['feature_name'], statics)
				index = index + 1
	for site in site_order:
		for i in range(1):
			dheader[index] = 'num_of_atoms_{}{}'.format(site, i)
			index = index + 1
			for disc in discList:
				dheader[index] = '{}{}_{}'.format(site, i, disc['feature_name'])
				index = index + 1
	statics_order = ['max', 'min', 'range']
	for statics in statics_order:
		for site in site_order:
			for disc in discList:
				dheader[index] = '{}_{}_{}'.format(site, disc['feature_name'], statics)
				index = index + 1
	dheader[index] = 'EnergyAboveHull'
	index = index + 1
	dheader[index] = 'Formation_energy'
	if (index != d_nums - 1):
		print('error in discrete feature name matching')
		print('index={}, d_nums={}'.format(index, d_nums))
	return dheader

def add_cheader(c_nums, contList):
	cheader = np.array(range(c_nums), dtype=object)
	cheader[0] = ''
	cheader[1] = 'goldschmidt_TF'
	cheader[2] = 'goldschmidt_TF_ionic'
	cheader[3] = 'octahedral_factor'
	cheader[4] = 'octahedral_factor_ionic'
	cheader[5] = 'A_O'
	cheader[6] = 'B_O'
	cheader[7] = 'A_B'
	index = 8
	site_order = ['Asite', 'Bsite']
	for site in site_order:
		for i in range(1):
			cheader[index] = 'num_of_atoms_{}{}'.format(site, i)
			index = index + 1
			cheader[index] = 'shannon_radii_{}{}'.format(site, i)
			index = index + 1
			for disc in contList:
				cheader[index] = '{}{}_{}'.format(site, i, disc['feature_name'])
				index = index+1
	AB_order = ['AB_avg', 'AB_diff', 'AB_ratio']
	for statis in AB_order:
		cheader[index] = 'shannon_radii_{}'.format(statis)
		index = index + 1
		for disc in contList:
			cheader[index] = '{}_{}'.format(disc['feature_name'], statis)
			index = index + 1
	statis_order = ['weighted_avg', 'max', 'min', 'range']
	for statis in statis_order:
		for site in site_order:
			cheader[index] = '{}_shannon_radii_{}'.format(site, statis)
			index = index + 1
			for disc in contList:
				cheader[index] = '{}_{}_{}'.format(site, disc['feature_name'], statis)
				index = index + 1
	cheader[index] = 'EnergyAboveHull'
	index = index + 1
	cheader[index] = 'Formation_energy'

	if (index != c_nums - 1):
		print('error in continuous feature name matching')
		print('index={}, c_nums={}'.format(index, c_nums))

	return cheader

def gen_vector(book_name, contList, discList):
	book = open_workbook(book_name)
	sheet = book.sheet_by_index(0)

	dexamples = []
	cexamples = []
	for row_index in range(1, sheet.nrows):  # skip header line
		composition = sheet.cell(row_index, 0).value
		stability = sheet.cell(row_index, 14).value
		formation_E = sheet.cell(row_index, 15).value
		numberofelements = sheet.cell(row_index, 8).value

		Acontsites = getComp(row_index, sheet, 1, composition, contList, True)
		Bcontsites = getComp(row_index, sheet, 4, composition, contList, True)
		Xcontsites = getComp(row_index, sheet, 7, composition[-4:], contList, True)
		Astatics = np.array(Acontsites)[:, 1:len(Acontsites[0])].astype(float)
		Bstatics = np.array(Bcontsites)[:, 1:len(Bcontsites[0])].astype(float)
		Xstatics = np.array(Xcontsites)[:, 1:len(Xcontsites[0])].astype(float)

		# Asingles = [np.zeros((len(Astatics[0]))), np.zeros((len(Astatics[0]))), np.zeros((len(Astatics[0])))]
		# Bsingles = [np.zeros((len(Bstatics[0]))), np.zeros((len(Bstatics[0]))), np.zeros((len(Bstatics[0])))]
		# for i in range(Astatics.shape[0]):
		# 	Asingles[i] = Astatics[i]
		# for i in range(Bstatics.shape[0]):
		# 	Bsingles[i] = Bstatics[i]


		Adiscsites = getComp(row_index, sheet, 1, composition, discList)
		Bdiscsites = getComp(row_index, sheet, 4, composition, discList)
		Adstatics = np.array(Adiscsites)[:, 1:len(Adiscsites[0])].astype(float)
		Bdstatics = np.array(Bdiscsites)[:, 1:len(Bdiscsites[0])].astype(float)

		# Adsingles = [np.zeros((len(Adstatics[0]))), np.zeros((len(Adstatics[0]))), np.zeros((len(Adstatics[0])))]
		# Bdsingles = [np.zeros((len(Bdstatics[0]))), np.zeros((len(Bdstatics[0]))), np.zeros((len(Bdstatics[0])))]
		# for i in range(Adstatics.shape[0]):
		# 	Adsingles[i] = Adstatics[i]
		# for i in range(Bdstatics.shape[0]):
		# 	Bdsingles[i] = Bdstatics[i]

		Atot = np.sum(Astatics[:, 0], axis=0)
		Btot = np.sum(Astatics[:, 0], axis=0)
		Xtot = np.sum(Xstatics[:, 0], axis=0)
		if (Atot != 8 or Btot != 8 or Xtot != 24):
			print(Atot, Btot, Xtot, row_index, composition)
		Acoef = Astatics[:, 0] / Atot if Atot != 0 else 0
		Bcoef = Bstatics[:, 0] / Btot if Btot != 0 else 0

		Asite_weighted = np.dot(Acoef, Astatics)[1:]
		Asite_max = np.amax(Astatics, axis=0)[1:]
		Asite_min = np.amin(Astatics, axis=0)[1:]
		Asite_ptp = np.ptp(Astatics, axis=0)[1:]
		Bsite_weighted = np.dot(Bcoef, Bstatics)[1:]
		Bsite_max = np.amax(Bstatics, axis=0)[1:]
		Bsite_min = np.amin(Bstatics, axis=0)[1:]
		Bsite_ptp = np.ptp(Bstatics, axis=0)[1:]
		Xsite_avg = np.mean(Xstatics, axis=0)[1:]

		Asite_dweighted = np.dot(Acoef, Adstatics)[1:]
		Asite_dmax = np.amax(Adstatics, axis=0)[1:]
		Asite_dmin = np.amin(Adstatics, axis=0)[1:]
		Asite_dptp = np.ptp(Adstatics, axis=0)[1:]
		Bsite_dweighted = np.dot(Bcoef, Bdstatics)[1:]
		Bsite_dmax = np.amax(Bdstatics, axis=0)[1:]
		Bsite_dmin = np.amin(Bdstatics, axis=0)[1:]
		Bsite_dptp = np.ptp(Bdstatics, axis=0)[1:]

		with np.errstate(divide='ignore', invalid='ignore'):
			ABPRatio = np.true_divide(Asite_weighted, Bsite_weighted)
			ABPRatio[ABPRatio == np.inf] = 0
			ABPRatio[ABPRatio == - np.inf] = 0
			ABPRatio = np.nan_to_num(ABPRatio)

		AB_avg = (Asite_weighted + Bsite_weighted) / 2
		AB_diff = (Asite_weighted - Bsite_weighted) / 2
		ra = Asite_weighted[0]
		rb = Bsite_weighted[0]
		rx = Xsite_avg[0]
		goldschmidt_TF = (ra + rx) / (np.sqrt(2) * (rb + rx))
		octahedral = rb/rx
		A_O = ra + rx
		B_O = rb + rx
		A_B = ra + rb
		ra = Asite_weighted[1]
		rb = Bsite_weighted[1]
		rx = Xsite_avg[1]
		goldschmidt_TF_ionic = (ra + rx) / (np.sqrt(2) * (rb + rx))
		octahedral_ionic = rb / rx


		dexample = np.r_[row_index, numberofelements, Asite_dweighted, Bsite_dweighted, Adstatics[0], Bdstatics[0],
		                 #Adsingles[0], Adsingles[1], Adsingles[2], Bdsingles[0], Bdsingles[1], Bdsingles[2],

		                 Asite_dmax, Bsite_dmax, Asite_dmin, Bsite_dmin, Asite_dptp, Bsite_dptp, stability, formation_E]
		dexamples.append(dexample)

		cexample = np.r_[row_index, goldschmidt_TF, goldschmidt_TF_ionic, octahedral, octahedral_ionic, A_O, B_O, A_B,
		                 #Asingles[0], Asingles[1], Asingles[2], Bsingles[0], Bsingles[1], Bsingles[2],
		                 Astatics[0], Bstatics[0],
		                 AB_avg, AB_diff, ABPRatio, Asite_weighted, Bsite_weighted,
		                 Asite_max, Bsite_max, Asite_min, Bsite_min, Asite_ptp, Bsite_ptp,
		                 stability, formation_E]
		cexamples.append(cexample)

	dexamples.insert(0, add_dheader(dexamples[0].size, discList))
	cexamples.insert(0, add_cheader(cexamples[0].size, contList))
	return cexamples, dexamples

def get_dict_list(book_name):
	refbk = open_workbook(book_name)
	refst = refbk.sheet_by_index(0)
	disc_list = []
	for attribute_index in range(1, refst.ncols):
		disc_list.append(buildDict(attribute_index, refst, 0))
	return disc_list

def write_to_csv(output_type, data):
	with open('{}.csv'.format(output_type), 'w', encoding='utf-8') as csv_file:
		writer = csv.writer(csv_file, delimiter=',')
		writer.writerows(data)

def gen_train(filename='training_set3.xlsx', output=0):
	contList = get_dict_list('continousproperty.xlsx')
	discList = get_dict_list('discrictproperty.xlsx')
	cexamples, dexamples = gen_vector(filename, contList, discList)
	write_to_csv('c_{}_train'.format(output), cexamples)
	write_to_csv('d_{}_train'.format(output), dexamples)

def gen_test(filename='Testing_set3.xlsx', output=0):
	contList = get_dict_list('continousproperty.xlsx')
	discList = get_dict_list('discrictproperty.xlsx')
	ctexamples, dtexamples = gen_vector(filename, contList, discList)
	write_to_csv('c_{}_test'.format(output), ctexamples)
	write_to_csv('d_{}_test'.format(output), dtexamples)

if __name__ == "__main__":
	gen_train()
	if len(sys.argv) > 1:
		gen_test()




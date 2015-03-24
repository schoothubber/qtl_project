"""
This script should handle the calculation of a confusion matrix

Based on:

input data from a model
reference data from a database of ARGIS
"""


def identify_true_false_positives(predicted_relations, TG_TF_ref):
	"""
	Take a list of [TG-TF] which are taken from the enriched genelists with
	a [TG-TF] reference list from AtRegNet
	
	Now check which [TG-TF] prediction is True
		and which [TG-TF] prediction is False
		
	predicted_relations = [TG-TF]
	"""
	
		
	true_relation = []
	false_relation = []
		
	for relation in predicted_relations:
		
		if relation in TG_TF_ref:
			
			true_relation.append(relation)
			
		else:
			false_relation.append(relation)

	return true_relation, false_relation
			


def count_false_negatives(total_FN, pred_TP, true_traits):
	"""
	The leftover False Negatives should be the total possible FN
	minus the actual TP
	The total possible FN are those TG-TF that:
	-are in AtRegNet
	-have a TG which is in the dataset with eQTL (true_trait)
	
	
	Also only take those TG-TF references into account that contain a TG
	that was identified as a true trait in get_TGTF_from_genelist()
	"""
	
	
	permitted_FN = []
	for red_dot in total_FN:
		#the FN has to come from a trait that has an eQTL
		#and cannot be correctly predicted by the model
		#because then it is a purple dot
		if red_dot[0] in true_traits and red_dot not in pred_TP:
			permitted_FN.append(red_dot)
			#print len(permitted_FN)
			
	return permitted_FN
	

def calculate_confusion(total_rel, true_pred_rel, false_pred_rel, unpredicted_rel):
	"""
	"""
	total = len(total_rel)
	
	TP = float(len(true_pred_rel))
	FP = float(len(false_pred_rel))
	FN = float(len(unpredicted_rel))
	TN = total - TP - FP - FN
	
	if TP + FN != 0:
		recall = TP / (TP + FN)
	else:
		#pass
		recall = None
		
	if TN + FP != 0:
		specif = TN / (TN + FP)
	else:
		#pass
		specif = None
	
	if TP + FP != 0:
		precis = TP / (TP + FP)
	else:
		#pass
		precis = None
		
	if recall != None and precis != None and precis+recall != 0:
		bal_F1 = 2.0 * (precis * recall) / (precis + recall)
	else:
		bal_F1 = None
		
	
	return TP, FP, FN, TN, recall, specif, precis, bal_F1
	
def print_results(dataset, cutoff, TP, FP, FN, TN, recall, specif, precis, F1):
	"""
	"""
	print "confusion matrix"
	print "dataset: %s"%dataset
	print "cutoff: %s"%cutoff
	print "---------------------------------"
	print "TP %s \t FN %s"%(TP, FN)
	print "FP %s \t TN %s"%(FP, TN)
	print "---------------------------------"
	if recall != None:
		print "recall: %s"%round(recall, 7)
	else:
		print "recall: NA"
	if specif != None:
		print "specificity: %s"%round(specif, 7)
	else:
		print "specificity: NA"
	if precis != None:
		print "precision: %s"%round(precis, 7)
	else:
		print "precision: NA"
	if F1 != None:
		print "F1: %s"%round(F1, 7)
	else:
		print "F1: NA"
	print "#################################"	








	
	

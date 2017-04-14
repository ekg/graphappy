"""
Test genotyping of pedigrees
"""
from nose.tools import raises
from whatshap.core import GenotypeDPTable, ReadSet, Variant, Pedigree, NumericSampleIds, PhredGenotypeLikelihoods
from whatshap.pedigree import centimorgen_to_phred
from .phasingutils import string_to_readset, string_to_readset_pedigree, brute_force_phase


def genotype_pedigree(reads, recombcost, pedigree, expected_genotypes, distrust_genotypes=False, positions=None):
	rs = string_to_readset_pedigree(reads, scaling_quality=10)
	dp_forward_backward = GenotypeDPTable(rs, recombcost, pedigree, positions)
	
	# for each position compare the likeliest genotype to the expected ones
	positions = rs.get_positions()
	for pos in range(len(positions)):
		for individual in range(len(pedigree)):
			likelihoods = dp_forward_backward.get_genotype_likelihoods(individual,pos)
			# find the likeliest genotype
			max_val = -1
			max_index = -1
			for i in range(len(likelihoods)):
				if likelihoods[i] > max_val:
					max_val = likelihoods[i]
					max_index = i
					
			# compare it to the expected genotype
			print(likelihoods,expected_genotypes[individual][pos])
			assert(max_index == expected_genotypes[individual][pos])
		print("\n")							
	
def test_genotyping_empty_trio():
	rs = ReadSet()
	recombcost = []
	pedigree = Pedigree(NumericSampleIds())
	pedigree.add_individual('individual0', [])
	pedigree.add_individual('individual1', [])
	pedigree.add_individual('individual2', [])
	pedigree.add_relationship('individual0', 'individual1', 'individual2')
	dp_forward_backward = GenotypeDPTable(rs, recombcost, pedigree)


def test_genotyping_trio1():
	reads = """
	  A 00
	  A 00
	  B 11
	  B 11
	  C 11
	  C 00
	"""

	expected_genotypes = [[0,0] , [2,2], [1,1]]
	pedigree = Pedigree(NumericSampleIds())
	pedigree.add_individual('individual0',[1,1])
	pedigree.add_individual('individual1',[1,1])
	pedigree.add_individual('individual2',[1,1])
	pedigree.add_relationship('individual0', 'individual1', 'individual2')
	recombcost = [10,10]
	genotype_pedigree(reads, recombcost, pedigree, expected_genotypes)

def test_genotyping_quartet1():
       reads = """
         A 1111
         A 0000
         B 1010
         C 111000
         C 010101
         D 000000
         D 010
         B 0101
         C  1100
         D  10010
         A   0000
         A   1111
         B   1010
         B   0101
       """
       expected_genotypes = [[1,1,1,1,1,1], [1,1,1,1,1,1], [1,2,1,1,0,1], [0,1,0,0,1,0]]
       pedigree = Pedigree(NumericSampleIds())
       pedigree.add_individual('individual0', [0,0,0,0,0,0])
       pedigree.add_individual('individual1', [0,0,0,0,0,0])
       pedigree.add_individual('individual2', [0,0,0,0,0,0])
       pedigree.add_individual('individual3', [0,0,0,0,0,0])
       pedigree.add_relationship('individual0', 'individual1', 'individual2')
       pedigree.add_relationship('individual0', 'individual1', 'individual3')
       recombcost = [3,3,3,4,3,3]
       genotype_pedigree(reads, recombcost, pedigree, expected_genotypes)


#def test_phase_trio2():
#	reads = """
#	  A 00
#	  A 00
#	  B 11
#	  B 11
#	  C 11
#	  C 00
#	"""
#	pedigree = Pedigree(NumericSampleIds())
#	pedigree.add_individual('individual0', [2,2])
#	pedigree.add_individual('individual1', [0,0])
#	pedigree.add_individual('individual2', [1,1])
#	pedigree.add_relationship('individual0', 'individual1', 'individual2')
#	recombcost = [10,10,10]
#	superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
#	assert cost == 8
#	assert len(set(transmission_vector)) == 1
#	all_expected_haplotypes = [
#		('11','11'),
#		('00','00'),
#		('00','11')
#	]
#	assert_haplotypes(superreads_list, all_expected_haplotypes, 2)
#
#
#def test_phase_trio3():
#	reads = """
#	  A 1111
#	  B 1010
#	  C 111000
#	  C 010101
#	  B 0101
#	  A  0000
#	  B  1010
#	  C  1010
#	  C  1100
#	  A   0000
#	  A   1111
#	  B   1010
#	  B    010
#	"""
#	pedigree = Pedigree(NumericSampleIds())
#	pedigree.add_individual('individual0', [1,1,1,1,1,1])
#	pedigree.add_individual('individual1', [1,1,1,1,1,1])
#	pedigree.add_individual('individual2', [1,2,1,1,0,1])
#	pedigree.add_relationship('individual0', 'individual1', 'individual2')
#	recombcost = [3,3,3,4,3,3]
#	superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
#	assert cost == 4
#	assert transmission_vector in ([0,0,0,1,1,1], [1,1,1,0,0,0], [2,2,2,3,3,3], [3,3,3,2,2,2])
#	all_expected_haplotypes = [
#		('111111','000000'),
#		('010101','101010'),
#		('111000','010101')
#	]
#	assert_haplotypes(superreads_list, all_expected_haplotypes, 6)
#
#
#def test_phase_trio4():
#	reads = """
#	  B 101
#	  B 101
#	  B 101
#	  A 111
#	  A 111
#	  A 111
#	  C 111
#	  C 111
#	  C 111
#	"""
#	pedigree = Pedigree(NumericSampleIds())
#	pedigree.add_individual('individual0', [1,1,1])
#	pedigree.add_individual('individual1', [1,1,1])
#	pedigree.add_individual('individual2', [1,1,1])
#	pedigree.add_relationship('individual0', 'individual1', 'individual2')
#	recombcost = [1,1,1]
#	superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
#	assert cost == 2
#	assert transmission_vector in ([0,2,0], [2,0,2], [1,3,1], [3,1,3])
#	all_expected_haplotypes = [
#		('111','000'),
#		('101','010'),
#		('111','000')
#	]
#	assert_haplotypes(superreads_list, all_expected_haplotypes, 3)
#
#
#def test_phase_trio5():
#	reads = """
#	  B 101
#	  B 101
#	  B 101
#	  A 111
#	  A 111
#	  A 111
#	  C 111
#	  C 111
#	  C 111
#	"""
#	pedigree = Pedigree(NumericSampleIds())
#	pedigree.add_individual('individual0', [1,1,1])
#	pedigree.add_individual('individual1', [1,1,1])
#	pedigree.add_individual('individual2', [1,1,1])
#	pedigree.add_relationship('individual0', 'individual1', 'individual2')
#	recombcost = [2,2,2]
#	superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
#	assert cost == 3
#	assert len(set(transmission_vector)) == 1
#	all_expected_haplotypes = [
#		('111','000'),
#		('111','000'),
#		('111','000')
#	]
#	assert_haplotypes(superreads_list, all_expected_haplotypes, 3)
#
#
#def test_phase_trio_pure_genetic():
#	reads = ""
#	pedigree = Pedigree(NumericSampleIds())
#	pedigree.add_individual('individual0', [2,1,1,0])
#	pedigree.add_individual('individual1', [1,2,2,1])
#	pedigree.add_individual('individual2', [1,1,1,0])
#	pedigree.add_relationship('individual0', 'individual1', 'individual2')
#	recombcost = [2,2,2]
#	superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree, positions=[10,20,30,40])
#	assert cost == 0
#	assert len(set(transmission_vector)) == 1
#	all_expected_haplotypes = [
#		('1110','1000'),
#		('1111','0110'),
#		('1000','0110')
#	]
#	assert_haplotypes(superreads_list, all_expected_haplotypes, 4)
#
#
#def test_phase_doubletrio_pure_genetic():
#	reads = ""
#	pedigree = Pedigree(NumericSampleIds())
#	pedigree.add_individual('individualA', [1,2,1,0])
#	pedigree.add_individual('individualB', [1,0,1,1])
#	pedigree.add_individual('individualC', [2,1,1,0])
#	pedigree.add_individual('individualD', [1,2,2,1])
#	pedigree.add_individual('individualE', [1,1,1,0])
#	pedigree.add_relationship('individualA', 'individualB', 'individualC')
#	pedigree.add_relationship('individualC', 'individualD', 'individualE')
#	recombcost = [2,2,2]
#	superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree, positions=[10,20,30,40])
#	assert cost == 0
#	assert len(set(transmission_vector)) == 1
#	all_expected_haplotypes = [
#		('0100','1110'),
#		('0011','1000'),
#		('1110','1000'),
#		('1111','0110'),
#		('1000','0110')
#	]
#	assert_haplotypes(superreads_list, all_expected_haplotypes, 4)
#
#
#def test_phase_quartet1():
#	reads = """
#	  A 111
#	  A 010
#	  A 110
#	  B 001
#	  B 110
#	  B 101
#	  C 001
#	  C 010
#	  C 010
#	  D 001
#	  D 010
#	  D 010
#	"""
#	pedigree = Pedigree(NumericSampleIds())
#	pedigree.add_individual('individual0', [1,2,1])
#	pedigree.add_individual('individual1', [1,1,1])
#	pedigree.add_individual('individual2', [0,1,1])
#	pedigree.add_individual('individual3', [0,1,1])
#	pedigree.add_relationship('individual0', 'individual1', 'individual2')
#	pedigree.add_relationship('individual0', 'individual1', 'individual3')
#	recombcost = [10,10,10]
#	superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
#	assert cost == 2
#	assert len(set(transmission_vector)) == 1
#	all_expected_haplotypes = [
#		('111','010'),
#		('001','110'),
#		('001','010'),
#		('001','010')
#	]
#	assert_haplotypes(superreads_list, all_expected_haplotypes, 3)
#
#
#def test_phase_quartet2():
#	reads = """
#	  A 111111
#	  A 000000
#	  B 010101
#	  B 101010
#	  C 000000
#	  C 010101
#	  D 000000
#	  D 010101
#	"""
#	pedigree = Pedigree(NumericSampleIds())
#	pedigree.add_individual('individual0', [1,1,1,1,1,1])
#	pedigree.add_individual('individual1', [1,1,1,1,1,1])
#	pedigree.add_individual('individual2', [0,1,0,1,0,1])
#	pedigree.add_individual('individual3', [0,1,0,1,0,1])
#	pedigree.add_relationship('individual0', 'individual1', 'individual2')
#	pedigree.add_relationship('individual0', 'individual1', 'individual3')
#	recombcost =[3,3,3,3,3,3]
#
#	superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
#	assert cost == 0
#	assert len(set(transmission_vector)) == 1
#	all_expected_haplotypes = [
#		('111111','000000'),
#		('010101','101010'),
#		('000000','010101'),
#		('000000','010101')
#	]
#	assert_haplotypes(superreads_list, all_expected_haplotypes, 6)
#
#
#def test_phase_quartet3():
#	reads = """
#	  A 1111
#	  A 0000
#	  B 1010
#	  C 111000
#	  C 010101
#	  D 000000
#	  D 010
#	  B 0101
#	  C  1100
#	  D  10010
#	  A   0000
#	  A   1111
#	  B   1010
#	  B   0101
#	"""
#	pedigree = Pedigree(NumericSampleIds())
#	pedigree.add_individual('individual0', [1,1,1,1,1,1])
#	pedigree.add_individual('individual1', [1,1,1,1,1,1])
#	pedigree.add_individual('individual2', [1,2,1,1,0,1])
#	pedigree.add_individual('individual3', [0,1,0,0,1,0])
#	pedigree.add_relationship('individual0', 'individual1', 'individual2')
#	pedigree.add_relationship('individual0', 'individual1', 'individual3')
#	recombcost = [3,3,3,4,3,3]
#	superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
#	print(cost)
#	print(transmission_vector)
#	assert cost == 8
#	# TODO: expect transmission in both trio relations. Update once transmission vectors
#	#       are returned per trio relationship
#	#assert transmission_vector in ([0,0,0,1,1,1], [1,1,1,0,0,0], [2,2,2,3,3,3], [3,3,3,2,2,2])
#	all_expected_haplotypes = [
#		('111111','000000'),
#		('010101','101010'),
#		('111000','010101'),
#		('000000','010010')
#	]
#	assert_haplotypes(superreads_list, all_expected_haplotypes, 6)
#
#
#def test_centimorgen_to_phred():
#	assert round(centimorgen_to_phred(0.10010013353365396)) == 30
#	assert round(centimorgen_to_phred(0.0010000100001343354)) == 50
#	assert round(centimorgen_to_phred(1e-38)) == 400
#
#
#@raises(ValueError)
#def test_centimorgen_to_phred_zero():
#	assert centimorgen_to_phred(0)
#
#
#def test_phase_trio_genotype_likelihoods():
#	reads = """
#	  A 111
#	  A 010
#	  A 110
#	  B 001
#	  B 110
#	  B 101
#	  C 001
#	  C 010
#	  C 010
#	"""
#	pedigree = Pedigree(NumericSampleIds())
#	genotype_likelihoods_mother = [
#		PhredGenotypeLikelihoods(0,0,0),
#		PhredGenotypeLikelihoods(0,0,1),
#		PhredGenotypeLikelihoods(5,0,5)
#	]
#	genotype_likelihoods0 = [PhredGenotypeLikelihoods(0,0,0)] * 3
#	pedigree.add_individual('individual0', [0,0,0], genotype_likelihoods_mother)
#	pedigree.add_individual('individual1', [0,0,0], genotype_likelihoods0)
#	pedigree.add_individual('individual2', [0,0,0], genotype_likelihoods0)
#	pedigree.add_relationship('individual0', 'individual1', 'individual2')
#	recombcost = [10,10,10]
#	superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree, True)
#	assert cost == 3
#	assert len(set(transmission_vector)) == 1
#	all_expected_haplotypes = [
#		('111','010'),
#		('001','110'),
#		('001','010')
#	]
#	assert_haplotypes(superreads_list, all_expected_haplotypes, 3)

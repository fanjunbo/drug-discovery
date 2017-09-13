import math
import numpy as np
import random


class Atom:
	def __init__(self):
		self.atomId = -1
		self.coord = [0.0, 0.0, 0.0]
		self.atomType = None
		self.status = None
	 
	def __init__(self, atomLine):
		tokens = atomLine.strip().split()
		self.atomId = int(tokens[1])
		self.atomType = tokens[2]
		self.coord = [float(tokens[5]), float(tokens[6]), float(tokens[7])]
		self.status = tokens[-1]
	 
	def getAtomId(self):
		return self.atomId
	 
	def getAtomType(self):
		return self.atomType

	def getAtomCoord(self):
		return self.coord

	def getAtomStatus(self):
		return self.status






class Group:
	def __init__(self):
		self.atoms = []
		self.groupId = -1

	def setGroupId(self, groupId):
		self.groupId = groupId

	def getGroupId(self):
		return self.groupId

	def appendAtom(self, atom):
		self.atoms.append(atom)

	def removeAtomById(self, atomId):
		targetIndex = 0
		for i in range(len(self.atoms)):
			if self.atoms[i].getAtomId() == atomId:
				targetIndex
				break
		self.atoms.pop(targetIndex)

	def findAtomById(self, atomId):
		for atom in self.atoms:
			if atom.getAtomId() == atomId:
				return True
		else:
			return False

	def getAtoms(self):
		return self.atoms

	def getGroupSize(self):
		return len(self.atoms)





 
class BranchMatrix:
	def __init__(self):
		self.groupSize = 0
		self.connections = []

	def __init__(self, N):
		self.groupSize = N
		self.connections = [[0 for i in range(N)] for j in range(N)]

	def setGroupSize(self, N):
		self.__init__(N)

	def setBranch(self, source, target):
		self.connections[source][target] = self.connections[target][source] = 1

	def removeBranch(self, source, target):
		self.connections[source][target] = self.connections[target][source] = 0

	def branchStatus(self, source, target):
		return self.connections[source][target]

	def getConnectionsAsGroupId(self):
		res = [ [row, col] for row in range(len(self.connections)) for col in range(len(self.connections[row])) if self.connections[row][col] and row < col ]
		return res

	def getDofVector(self):
		res = [ self.connections[row].count(1) for row in range(len(self.connections)) ]
		return sorted(res, reverse=True)





 
class Ligand:
	def __init__(self, filename):
		# init necessary arguments for a Ligand
		# groups: list of object of Group
		# branch: object of BranchMatrix
		self.torsdof = 0
		self.groups = []
		self.branch = BranchMatrix(0)
		self.ligandAtoms = []
		self.ligandFilename = filename.split('/')[-1]

		self.groupVector = []
		self.dofVector = []
		self.groupTuple = []

		inputFile = open(filename, 'r')
		groupStack, branchConnections = [], []
		# find all groups, all natural branches i.f.o (atomId, atomId)
		for line in inputFile:
			line = line.strip()
			if ('ROOT' in line or 'BRANCH' in line) and ('END' not in line):
				groupStack.append(Group())
			if 'ATOM' in line:
				groupStack[-1].appendAtom(Atom(line))
			if 'ENDROOT' in line or 'ENDBRANCH' in line:
				tempGroup = groupStack[-1]
				self.groups.append(tempGroup)
				self.groups[-1].setGroupId(len(self.groups) - 1)
				groupStack.pop(-1)
				if 'ENDBRANCH' in line:
					branchConnections.append([int(line.split()[1]), int(line.split()[2])])
			if 'TORSDOF' in line:
				self.torsdof = int(line.split()[1])
		# init branch with specific group size
		self.branch.setGroupSize(len(self.groups))
		# find branch connections i.f.o (groupId, groupId) and update branch matrix
		for item in branchConnections:
			sourceGroupId, targetGroupId = 0, 0
			for group in self.groups:
				if group.findAtomById(item[0]):
					sourceGroupId = group.getGroupId()
				if group.findAtomById(item[1]):
					targetGroupId = group.getGroupId()
			self.branch.setBranch(sourceGroupId, targetGroupId)
		# init ligandAtoms 
		for group in self.groups:
			self.ligandAtoms += group.getAtoms()

		self.getGroupVector()
		self.getDofVector()
		self.getGroupTuple()

	def getLigandAtoms(self):
		return self.ligandAtoms

	def getLigandFilename(self):
		return self.ligandFilename

	def getGroupVector(self):
		for group in self.groups:
			self.groupVector.append(group.getGroupSize())
		self.groupVector = sorted(self.groupVector, reverse=True)

	def getDofVector(self):
		self.dofVector = self.branch.getDofVector()

	def getGroupTuple(self):
		res = []
		connections = self.branch.getConnectionsAsGroupId()
		for connect in connections:
			temp = []
			for group in self.groups:
				if group.getGroupId() in connect:
					temp.append(group.getGroupSize())
			temp = sorted(temp, reverse=True)
			res.append(temp)
		self.groupTuple = sorted(res, key=lambda x:(x[0], x[1]), reverse=True)
	
	def getLigandGroups(self):
		return self.groups





 
class USR:
	def __init__(self):
		self.result = []
		self.ligands = []
		self.targetLigand = None
		self.targetVector = []
		self.targets = []

	def setTarget(self, ligand):
		self.targetLigand = ligand
		self.targetVector = self.ligandToVector(ligand)
		self.targets.append(ligand)

	def appendComparativeLigand(self, ligand):
		#(ligand->Ligand)
		self.ligands.append(ligand)

	def manhattonDistance(self, fstPoint, secPoint):
		#(fstPoint->List of Number, secPoint->List of Number)
		distance = 0.0
		for i in range(len(fstPoint)):
			distance += (fstPoint[i] - secPoint[i]) ** 2
		return distance ** 0.5

	def ctdPoint(self, atoms):
		resPoint = [0.0, 0.0, 0.0]
		for atom in atoms:
			atomCoord = atom.getAtomCoord()
			for i in range(len(resPoint)):
				resPoint[i] += atomCoord[i]
		for i in range(len(resPoint)):
			resPoint[i] /= float(len(atoms))
		return resPoint

	def cstPoint(self, atoms):
		ctdPoint = self.ctdPoint(atoms)
		currentMinDistance = self.manhattonDistance(atoms[0].getAtomCoord(), ctdPoint)
		cstPoint = atoms[0].getAtomCoord()
		for i in range(1, len(atoms)):
			if self.manhattonDistance(atoms[i].getAtomCoord(), ctdPoint) < currentMinDistance:
				cstPoint = atoms[i].getAtomCoord()
		return cstPoint

	def fctPoint(self, atoms):
		ctdPoint = self.ctdPoint(atoms)
		currentMaxDistance = self.manhattonDistance(atoms[0].getAtomCoord(), ctdPoint)
		fctPoint = atoms[0].getAtomCoord()
		for i in range(1, len(atoms)):
			if self.manhattonDistance(atoms[i].getAtomCoord(), ctdPoint) > currentMaxDistance:
				fctPoint = atoms[i].getAtomCoord()
		return fctPoint

	def ftfPoint(self ,atoms):
		#(atoms->List of Atoms)
		fctPoint = self.fctPoint(atoms)
		currentMaxDistance = self.manhattonDistance(atoms[0].getAtomCoord(), fctPoint)
		ftfPoint = atoms[0].getAtomCoord()
		for i in range(1, len(atoms)):
			if self.manhattonDistance(atoms[i].getAtomCoord(), fctPoint) > currentMaxDistance:
				ftfPoint = atoms[i].getAtomCoord()
		return ftfPoint

	def getDistancesFromTarget(self, target, atoms):
		#(target->List of Number, atoms->List of Atoms)
		return [self.manhattonDistance(target, atom.getAtomCoord()) for atom in atoms]
		#return [self.manhattonDistance(target.getAtomCoord(), atom.getAtomCoord()) for atom in atoms]

	def kthMoment(self, K, numbers):
		#(K->Integer, numbers->List of Number)
		res = 0.0
		for number in numbers:
			res += number ** K
		res /= float(len(numbers))
		return res ** (1.0/K)


	def ligandToVector(self, ligand):
		#(ligand->Ligand)
		if type(ligand) == type([]):
			atoms = ligand
		else:
			atoms = ligand.getLigandAtoms()
		positions = [self.ctdPoint(atoms), self.cstPoint(atoms), self.fctPoint(atoms), self.ftfPoint(atoms)]
		res = []
		for point in positions:
			for i in range(1, 4):
				distances = self.getDistancesFromTarget(point, atoms)
				res.append(self.kthMoment(i, distances))
		return res


	def usrScore(self, target, source):
		#(target->Ligand, source->Ligand)
		targetVector, sourceVector = np.array(self.ligandToVector(target)), np.array(self.ligandToVector(source))
		error = np.array(map(lambda x: math.fabs(x), targetVector - sourceVector)).sum() / 12.0
		return 1.0 / (1.0 + error)

	def run(self):
		for ligand in self.ligands:
			vector = self.ligandToVector(ligand)
			self.result.append((ligand.getLigandFilename(), self.usrScore(self.targetLigand, ligand)))


	def sort(self):
		self.result = sorted(self.result, key=lambda x:x[1], reverse=True)

	def printTopK(self, K):
		for i in range(K):
			print self.result[i][0], ', ', self.result[i][1]

	def greedyUsrScore(self, source, target):
		sourceGroups, targetGroups = source.getLigandGroups(), target.getLigandGroups()
		sourceGroups = sorted(sourceGroups, key=lambda x:x.getGroupSize(), reverse=True)
		targetGroups = sorted(targetGroups, key=lambda x:x.getGroupSize(), reverse=True)
		if len(targetGroups) < len(sourceGroups):
			temp = targetGroups[:]
			targetGroups = sourceGroups[:]
			sourceGroups = temp
		weightScores = []
		currentMatchedIndex = []
		for fstGroup in sourceGroups:
			currentMaxUsr, currrentAtomsCount = 0.0, 0
			currentIndex = -1
			if fstGroup.getGroupSize() < 3:
				continue
			for i in range(len(targetGroups)):
				secGroup = targetGroups[i]
				if i in currentMatchedIndex:
					continue
				tempUsrScore = self.usrScore(fstGroup.getAtoms(), secGroup.getAtoms())
				if tempUsrScore > currentMaxUsr:
					currentMaxUsr = tempUsrScore
					currrentAtomsCount = fstGroup.getGroupSize() + secGroup.getGroupSize()
					currentIndex = i
			weightScores.append([currentMaxUsr, currrentAtomsCount])
			currentMatchedIndex.append(i)
		#print weightScores
		totalAtomsCount = 0
		for group in sourceGroups:
			totalAtomsCount += group.getGroupSize()
		for group in targetGroups:
			totalAtomsCount += group.getGroupSize()
		greedyUsrScore = 0
		for item in weightScores:
			greedyUsrScore += item[0] * (float(item[1]) / float(totalAtomsCount))
		return greedyUsrScore

	def dofVectorMatchScore(self, fstDofVec, secDofVec):
		if len(fstDofVec) > len(secDofVec):
			fstVec, secVec = secDofVec[:], fstDofVec[:]
		else:
			fstVec, secVec = fstDofVec[:], secDofVec[:]
		res = 0
		for i in range(len(fstVec)):
			for j in range(len(secVec)):
				if math.fabs(fstVec[i] - secVec[j]) <= 1:
					res = res + fstVec[i] + secVec[j]
					secVec[j] = -9999
					break
		total = 0
		for dof in fstDofVec:
			total += dof
		for dof in secDofVec:
			total += dof
		if total:
			return float(res) / float(total)
		else:
			return 0

	def groupVectorMatchScore(self, fstGroupVec, secGroupVec):
		if len(fstGroupVec) > len(secGroupVec):
			fstVec, secVec = secGroupVec[:], fstGroupVec[:]
		else:
			fstVec, secVec = fstGroupVec[:], secGroupVec[:]
		res = 0
		for i in range(len(fstVec)):
			for j in range(len(secVec)):
				if math.fabs(fstVec[i] - secVec[j]) <= 1:
					res += (fstVec[i] + secVec[j])
					secVec[j] = -9999
					break
		total = 0
		for dof in fstGroupVec:
			total += dof
		for dof in secGroupVec:
			total += dof
		return float(res) / float(total)

	def groupTupleMatchScore(self, fstGroupTuple, secGroupTuple):
		if len(fstGroupTuple) < len(secGroupTuple):
			fstTuple, secTuple = fstGroupTuple[:], secGroupTuple[:]
		else:
			fstTuple, secTuple = secGroupTuple[:], fstGroupTuple[:]
		res = total = 0
		matchedIndex = []
		for i in range(len(fstTuple)):
			for j in range(len(secTuple)):
				if fstTuple[i] == secTuple[j] and j not in matchedIndex:
					res += ((fstTuple[i][0] + secTuple[j][1]) * 2)
					matchedIndex.append(j)
					break
		for tup in fstGroupTuple:
			total += (tup[0] + tup[1])
		for tup in secGroupTuple:
			total += (tup[0] + tup[1])
		if total:
			return float(res) / float(total)
		else:
			return 0

	def filterByTarget(self, target):
		res = []
		base  =0.7
		step = 0.1
		while True:
			findEnoughLigands = False
			for ligand in self.ligands:
				#print str(target.dofVector) + '	' + str(ligand.dofVector)
				dofMatchScore = self.dofVectorMatchScore(target.dofVector, ligand.dofVector)
				groupMatchScore = self.groupVectorMatchScore(target.groupVector, ligand.groupVector)
				tupleMatchScore = self.groupTupleMatchScore(target.groupTuple, ligand.groupTuple)
				tempScore = (dofMatchScore + groupMatchScore + tupleMatchScore) / 3.0
				if tempScore > base:
					tempUsrScore = self.usrScore(target, ligand)
					greedyUsrScore = self.greedyUsrScore(target, ligand)
					res.append([ligand.getLigandFilename(), tempScore, tempUsrScore, greedyUsrScore])
			if len(res) < int(32904 * 0.01):
				base -= step
				res = []
			else:
				break
		res = sorted(res, key=lambda x:(x[1], x[3], x[2]), reverse=True)
		return res

	def filter(self):
		res = []
		target = self.targetLigand
		self.filterByTarget(target)

	def getEnrichmentFactor(self, ACTIVE, DECOY):
		NUM = ACTIVE + DECOY
		factors = []
		targetId = 0
		for target in self.targets:
			targetId += 1
			print 'target=', targetId
			results = self.filterByTarget(target)
			print 'target=', targetId, ' len=', len(results)
			if len(results) < int(NUM * 0.01):
				print 'error len=', len(results)
				continue
			topResults = [ results[:int(NUM * 0.01)],results[:int(NUM * 0.005)], results[:int(NUM * 0.0025)] ]
			outputFile = open('./targetResults/target' + str(targetId) + '_result.csv', 'w')
			for temp in topResults[0]:
				outputFile.write(temp[0] + ',' + str(temp[1]) + ',' + str(temp[2]) + ',' + str(temp[3]) + '\n')
			outputFile.close()
			enrichmentfactorResults = []
			for temp in topResults:
				activeCount, count = 0, 0
				for res in temp:
					if 'actives' in res[0]:
						activeCount += 1
					else:
						count += 1
				if count > 0:
					enrichmentfactorResults.append(float(activeCount) / float(count))
				else:
					enrichmentfactorResults.append(1.0)
			factors.append(enrichmentfactorResults)
		res = [0.0, 0.0, 0.0]
		for temp in factors:
			res[0] += temp[0]
			res[1] += temp[1]
			res[2] += temp[2]
		for i in range(len(res)):
			res[i] /= float(len(factors))
		base = float(ACTIVE) / float(DECOY)
		res = np.array(res) / base
		print res
		return res






target = Ligand('actives815.pdbqt')
usr = USR()
#usr.setTarget(target)
#for i in range(1, 842+1):
for i in range(601, 843):
	#seed = random.uniform(1, 50) 
	#if seed >= 1 and seed<= 2:
	usr.setTarget(Ligand('actives' + str(i) + '.pdbqt'))


#usr.setTarget(Ligand('actives641.pdbqt'))
for i in range(1, 842+1):
	usr.appendComparativeLigand(Ligand('actives' + str(i) + '.pdbqt'))
for i in range(1, 32062+1):
 	usr.appendComparativeLigand(Ligand('../decoys/decoys' + str(i) + '.pdbqt'))
print 'init done...'

usr.getEnrichmentFactor(841, 32062)


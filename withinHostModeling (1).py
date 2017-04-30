import random
import sys
import numpy
import math
import pylab

# Creates an empty list to store the fitnesses of virions in a population
pop = []

# Create an empty list that will store ID numbers for virions
ID = []

# Creates an empty list to store longitudinal population fitnesses
graphPop = []

# Creates an empty list to store a generation timeline
genList = []

# Creates an empty list that will keep track of the population's size
popSize = []

# Populates the population list with a starting fitness of 1. Population size is determined by
# the integer in range()
for x in range(10):
	pop.append(1)
	
for x in range(9):
	ID.append(1)

ID.append(2)
#for x in range(2, 12):
	#ID.append(x)

nextID = 3
	
graphPop.append(list(ID))

IDIteration = list(ID)

popSize.append(len(pop))

outFile = open("networks5.txt", 'w')
	
# Creates an iterator for the master while loop below
iter = 0


print("Entering loop\n")

# Master loop that will run the simulation however long we want, limited by iter
# Average infection has 2.6 day replication rate
# Patient 1 infection lasted 2966 days, ~ 1,000 replications
while iter < 12:
	
	for virion in range(len(pop)):
		# Influenza A mutation rate is 4.1e-3
		# V3 window size is 72 bases
		# 4.1e-3 x 72 = 0.2952
		# Probability 1-0.2952=0.7048
		
		# There is a 19% chance that an individual virion acquires a mutation
		if random.random() >= 0.7:
			
			# Create a variable newFit to hold the calculation of the new fitness, which is based on the old fitness
			newFit = pop[virion]
			
			# Generate a random number that will be used to assess whether the mutation adds, subtracts, or does not
			# affect fitness
			mutType = random.random()
			# This if statement series decides whether the mutation is beneficial or detrimental. If it's neither
			# the fitness of the virion will not change
			# If the mutation is beneficial
			if mutType >= 0.95:
				# Inside this if statement we need to assign a fitness increase associated with the beneficial mutation
				newFit+=0.3
				ID[virion] = nextID
				nextID+=1
			
			# If the mutation is detrimental
			if mutType <= 0.9:
				# Inside this elif statement we need to assign a fitness decrease associated with the detrimental mutation
				newFit-=0.8
				ID[virion] = nextID
				nextID+=1
				
			# If the viral sequence reverts to wild type
			else:
				newFit=1
				ID[virion] = 1
				
			# Set the new fitness for the current virion based on beneficial, detrimental or neutral mutation
			pop[virion] = newFit
			
			
	
	weightTotal = sum(pop)
	print(iter)
	
	# Creates an empty list to store calculated weights of each virion's fitness
	weightedPop = []
	
	for v in range(len(pop)):
		weightedPop.append(pop[v]/weightTotal)
	
	# Check to see if the population size has reached the maximum that can be supported
	# If it is smaller than the max of 1000, multiply the current size of the population by 10
	#if len(pop)*10 < 1000:
	
	if iter%2==0:
		newPopSize = len(pop)*3
	print(newPopSize)
	# If it is at or above the max, set the pop size to the max 
	# else:
		# newPopSize = 1000
	
	# Create an empty list that will hold the next generation of the population
	newPop = []
	newID = []
	
	# Randomly choose virions from the current population pop to populate the next generation population based
	# on the weighted scores of each individual virion
	for p in range(newPopSize):
		# Extract a tuple contained in a list which contains the index (which is also the ID)
		f = random.choices(list(enumerate(pop)), weightedPop)
		# Get the first value out of f, which is just a tuple
		g = f[0]
		# Get the second value out of the tuple g, is the fitness value
		value = g[1]
		# Get the first value out of the tuple g, is the index
		idindex = g[0]
		# Set id as the value held in the ID list at the idindex location
		id = ID[idindex]
		
		# Append the value to the newPop list
		newPop.append(value)
		# Append the id to the newID list
		newID.append(id)
		
	# Clear the original pop and ID lists to make room for the next generation
	pop.clear()
	ID.clear()
	

	for x in range(len(newPop)):
		pop.append(newPop[x])
		ID.append(newID[x])
		
	
	graphPop.append(list(ID))
	popSize.append(len(pop))
	
	iter+=1
	
print("Finished the loop\n")	

# Create an empty list that will be used to graph each generation
graphList = []
	

	
maxList = []
countList = []
lastLine = list(graphPop[-1])

# Count each instance of ID in lastLine which is the last line of graphPop
for j in range(len(IDIteration)):
	countList.append((lastLine.count(j), j))

# Place the (value, ID) pair into maxList
maxList.append(max(countList))

# Create maxPair which takes the tuple out of the list
maxPair = maxList[0]

# Create maxID which grabs the ID out of the maxPair pairing
maxID = maxPair[1]

# Loop through each line of graphPop
for line in graphPop:
	
	# Add only the count of the ID stored in maxID, instead of the counts of every ID
	graphList.append(line.count(maxID))

outFile.write("Graph list:\n" + str(graphList) + "\n")

percentList = []
for i in range(len(graphList)):
	percentList.append(graphList[i]/popSize[i])
	
outFile.write("Percent list:\n" + str(percentList) + "\n")
outFile.write("Population Size list:\n" + str(popSize) + "\n")
	
for i in range(13):
	genList.append(i)
	
outFile.write("Generation list:\n" + str(genList) + "\n")
	
pylab.plot(genList, percentList)
pylab.title('Within-Host Viral Modeling of V3')
pylab.xlabel('Generation')
pylab.ylabel('Percent of Viral ID')
pylab.ylim(0, 1)
pylab.show()
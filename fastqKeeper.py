#
# fastqKeeper.py
# Author: Philip Braunstein
#
# Date Created: May 13, 2014
# Last Modified: May 14, 2014
#
# The FastqKeeper class is used to hold the contents of a FASTQ file.
# It also provides general statistics about the FASTQ file such as G/C
# content, number of reads, and so forth.
#
# This class also contains functionality to get the number of sequneces
# that match a nucleotide pattern as well as the number of seqences
# whose average quality score is above a certain threshold.
#

import random
import re

class FastqKeeper:
        ################### Public API ###################
        def getNumReads(self):
                return self.numReads

        # Returns an list of [ID, SEQ, STRAND-SENSE, QSCORE, AVG_QSCORE]
        # it returns a random sequence 
        def getRandomSeq(self):
                return random.choice(self.seqs)    
        
        # Returns the number of reads that have the passed in
        # pattern. Ns in the pattern are interpreted as match any
        # one nucleotide each. 
        # Ns in the source sequence do not have this property.
        def getNumMatches(self, pattern):
                pattern = pattern.upper()

                # Throw exception if invalid characters passed in 
                self.validateNucs(pattern)

                newPattern = ""
                count = 0

                # Replace Ns with . for regex
                for i in range(len(pattern)):
                        if pattern[i] == 'N':
                                newPattern += '.'
                        else:
                                newPattern += pattern[i]

                for seq in self.seqs:
                        if re.search(newPattern, seq[self.SEQ_INDEX]) != None:
                                count += 1

                return count


        # Makes sure all characters in pattern are one 
        # of those listed in valid list in method.
        # Assumes that pattern has already been converted to upper case.
        def validateNucs(self, pattern):
                valid = ['A', 'C', 'G', 'T', 'N']
                for char in pattern:
                        if char not in valid:
                                raise Exception(self.NUC_EXCEPTION.\
                                        format(char))


        
        # Gets the number of sequences with a quality score at or above
        # the character passed in.
        def getNumQualSeqs(self, cutoff):
                if type(cutoff) != str:
                        raise Exception(self.STRING_EXCEPTION)

                if len(cutoff) != 1:
                        raise Exception(self.CHAR_EXCEPTION)

                count = 0
                val = ord(cutoff)

                for seq in self.seqs:
                        if seq[self.AVG_QSCORE_INDEX] >= val:
                                count += 1

                return count
                


        # Getters for nucleotide counts in whole file
        def getA(self):
                return self.A

        def getC(self):
                return self.C

        def getG(self):
                return self.G

        def getT(self):
                return self.T

        def getN(self):
                return self.N

        def getTotal(self):
                return self.total

        # Getters for GC / AT content
        def getGC(self):
                return self.GC

        def getAT(self):
                return self.AT


        # Getters for strand sense percentage
        def getPlus(self):
                return self.plus

        def getMinus(self):
                return self.minus

        
        # Getters for avg length and avg Q score
        def getAvgLen(self):
                return self.avgLen

        def getAvgQScore(self):
                return self.avgQScore



        ################ Private Methods ################
        # INITIALIZER
        # Inputs reads into data structure counts nucleotides
        # calculates avg Q scores for each read and avg length
        # and q scrore for whole file
        #
        # NOTE: A larger, less modularized readIn that called
        # all of the important functions during each line of read in
        # would be faster. Here speed is sacrificed for clarity of 
        # code deliberately.
        def __init__(self, inputFile):
                print "Initializing FastqKeeper:"
                self.initConstants()
                self.seqs = self.readIn(inputFile)
                self.verifyReads()
                self.numReads = len(self.seqs)
                self.countNucs()
                self.assignQScores()
                self.avgLenQScore()
                print "done."
                print


        # Sets constants that will be used by this class
        def initConstants(self):
                # Convenient indexing into fastq sequences
                self.ID_INDEX = 0
                self.SEQ_INDEX = 1
                self.STRAND_INDEX = 2
                self.QSCORE_INDEX = 3
                self.AVG_QSCORE_INDEX = 4

                # Number of lines per read in fastq file
                self.FASTQ_LINES = 4

                # Exceptions
                self.IO_EXCEPTION = "Can't find file \"{}\""
                self.NUC_EXCEPTION = "Unknown nucleotide \"{}\""
                self.BALANCE_EXCEPTION = "Read ID {} has a different" + \
                                " number of qscores and nucleotides"
                self.STRAND_EXCEPTION = "Unknown strand-sense \"{}\" in" + \
                                " read ID {}"
                self.CHAR_EXCEPTION = "Quality score chars must only be" + \
                                " one character long"
                self.STRING_EXCEPTION = "Expected string"
                 
                


        # Reads FASTQ file into list of lists s.t.
        # [[ID, SEQUENCE, STRAND-SENSE, QSCORE]]
        #
        # Raises exception if file doesn't exist
        # Undefined behavior with functions of different file
        # types.
        def readIn(self, inputFile):
                print "Reading in Sequences...."
                seqs = []
                newSeq = None
                count = 0

                try:
                        with open(inputFile, 'r') as filer:
                                for line in filer:
                                        line = line.strip()

                                        # found a new read
                                        if count == 0 or \
                                                count == self.FASTQ_LINES:
                                                if newSeq != None:
                                                        seqs.append(newSeq)

                                                newSeq = [line]

                                                # Reset counter
                                                count = 1
                                        else:
                                                newSeq.append(line)
                                                count += 1

                                # One remaining sequence not yet added
                                seqs.append(newSeq)
                except IOError:
                        raise Exception(self.IO_EXCEPTION.format(inputFile))

                return seqs

        # Checks to make sure that each read has an equal number of
        # quality scores and nucelotides as well as checks for
        # invalid strand-sense
        def verifyReads(self):
                print "Verifying reads...."
                for seq in self.seqs:
                        if len(seq[self.SEQ_INDEX]) != \
                                len(seq[self.QSCORE_INDEX]):
                                raise Exception(self.BALANCE_EXCEPTION.\
                                        format(seq[self.ID_INDEX]))

                        if seq[self.STRAND_INDEX] != '+' and  \
                                seq[self.STRAND_INDEX] != '-':
                                raise Exception(self.STRAND_EXCEPTION.\
                                        format(seq[self.STRAND_INDEX], \
                                                seq[self.ID_INDEX]))


        # Counts the number of nucleotides of each type in whole file
        #
        # Raises an exception if an unknown nucleotide is readIn
        def countNucs(self):
                print "Counting nucleotides...."
                self.A = 0
                self.C = 0
                self.G = 0
                self.T = 0
                self.N = 0

                # Capitalization probably unecessary, but a good idea
                for seq in self.seqs:
                        for char in seq[self.SEQ_INDEX]:
                                if char == 'A' or char == 'a':
                                        self.A += 1
                                elif char == 'C' or char == 'c':
                                        self.C += 1
                                elif char == 'G' or char == 'g':
                                        self.G += 1
                                elif char == 'T' or char == 't':
                                        self.T += 1
                                elif char == 'N' or char == 'n':
                                        self.N += 1
                                else:  # Unknown nucleotide throw exception
                                        raise Exception(self.NUC_EXCEPTION.\
                                                format(char))

                self.total = self.A + self.C + self.G + self.T + \
                                self.N

                self.calculateGC_AT()


        # Calculates GC and AT content to the nearest integer
        # and stores them in the class instance
        def calculateGC_AT(self):
                gc = float(self.G + self.C)
                gc = gc / (self.total - self.N)  # ignore Ns
                gc *= 100

                # to 2 decimal places
                gc = round(gc, 2)

                self.GC = gc
                self.AT = 100 -gc


        # Adds a 5th element to each read in self.seqs
        # that is the average quality score of that read
        def assignQScores(self):
                print "Assessing quality...."
                for seq in self.seqs:
                        seq.append(self.avgQScore(seq))


        # Takes in a read of the form [ID, SEQ, STRAND-SENSE, QSCORE]
        # and returns an int with rounded avg Q-score for that read
        def avgQScore(self, seq):
                qScore = 0
                for letter in seq[self.QSCORE_INDEX]:
                        qScore += ord(letter)  # Get numerical equivalent

                # Get average, round to two decimal places
                return round(qScore / float(len(seq[self.QSCORE_INDEX])), 2)


        def avgLenQScore(self):
                length = 0
                score = 0
                plus = 0

                # Sum all lengths and qscores
                for seq in self.seqs:
                        length += len(seq[self.SEQ_INDEX])
                        score += seq[self.AVG_QSCORE_INDEX]

                        # Only count + strands to get distribution of strand
                        # sense
                        if seq[self.STRAND_INDEX] == '+':
                                plus += 1

                self.avgLen = float(length) / self.getNumReads()
                self.avgQScore = float(score) / self.getNumReads()

                # Round to 2 decimals
                self.plus = round(float(plus) / self.numReads, 2) * 100
                self.minus = 100 - self.plus



        ################### Testing methods #######################
        # This is a testing function
        def twriteOut(self):
                with open ("OUT_TEST", 'w') as filew:
                        for seq in self.seqs:
                                try:
                                        filew.write(seq[0] + "\n")
                                        filew.write(seq[1] + "\n")
                                        filew.write(seq[2] + "\n")
                                        filew.write(seq[3] + "\n")
                                except IndexError:
                                        print seq

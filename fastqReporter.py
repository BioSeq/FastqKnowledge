#
# fastqReporter.py
# Author: Philip Braunstein
#
# Date Created: May 13, 2014
# Last Modified: May 18, 2014
#
# The FastqReporter class inherits from FastqKeeper. It adds functinoality
# to print descriptions of features stored in FastqKeeper in a nicer format.
# This modularity allows a program to call FastqReporter methods to get
# a nice report printed to stdout.
#

from fastqKeeper import FastqKeeper
from math import ceil

class FastqReporter(FastqKeeper):
        ################## Public API #######################
        # Prints the total number of reads
        def printNumReads(self):
                print self.getNumReads(), "total reads"
                print

        # Wrappers for getters of avg length and quality score
        def printAvgReadLength(self):
                print "Average read length:", self.avgLen, "nucleotides"
                print

        def printAvgQualScore(self):
                print "Average quality score:",

                # Convert avg to int then unichar to get ascii
                print unichr(int(round(self.avgQScore)))
                print

        # Prints a randomly chosen sequence in the format that
        # it appeared in the .fastq file
        def printRandomSeq(self):
                seq = self.getRandomSeq()
                seq = seq[:-1]  # Chop off avg q score
                for item in seq:
                        print item
                print


        # Prints a report of the number of nucleotides in the file
        def printNumNucs(self):
                print "NUMBER OF EACH NUCLEOTIDE SEEN IN ALL READS"
                print "A:", self.A
                print "C:", self.C
                print "G:", self.G
                print "T:", self.T
                print "N:", self.N
                print "Total:", self.total
                print


        # Prints G/C A/T content
        def printGC_AT(self):
                print "G/C A/T CONTENT FOR ALL READS"
                print "G/C Content:", str(self.GC) + "%"
                print "A/T Content:", str(self.AT) + "%"
                print


        # Prints percentage of strandedness
        def printStrandedness(self):
                print "STRANDEDNESS FOR ALL READS"
                print "Positive Strand:", str(self.plus) + "%"
                print "Negative Strand:", str(self.minus) + "%"
                print

       
        # Some nice formatting around getNumMatches inherited method 
        def printNumMatches(self, pattern):
                # Call inherited helper function
                num = self.getNumMatches(pattern)
                print "\"" + pattern + "\"", "found in", num, "out of", 
                print self.numReads, "sequences"
                print "(" + str(round(float(num) / self.numReads, 2) \
                        * 100) + "%)"
                print

        # Some nice formatting around the getNumQualSeqs inherited method
        def printNumQualSeqs(self, cutoff):
                num = self.getNumQualSeqs(cutoff)
                print num, "of", self.numReads,
                print "sequences have an average quality score greater than or",
                print "equal to", "\"" +  cutoff + "\""
                print "(" + str(round(float(num) / self.numReads, 2) \
                        * 100) + "%)"
                print


        def getTotalCharsToWrite(self):
                return self.totalCharsToWrite


        # Print the number of pages it would take to write out
        # the FASTQ file in 12 point font.
        def printTotalNumPages(self):
                numPages = int(ceil(self.totalCharsToWrite / float( \
                                        self.charsPerPage)))
                print "It would take", numPages, "pages to write",
                print "this FASTQ file in 12 point font"



        ##################### Private Methods ###########################
        # INITIALIZER
        # Initializer needs to count the total number of chars
        # that the file would take to writte out. First, it calls the
        # super class initializer then runs a method to do this.
        #
        def __init__(self, inputFile):
                FastqKeeper.__init__(self, inputFile)
                print "Initializing FastqReporter:"
                self.charsPerPage = 2812
                self.totalCharsToWrite = self.countTotalCharsToWrite()
                print "done."
                print 


        # Counts the total number of chars that would be required to
        # write out a FASTQ file including the line breaks between the lines. 
        def countTotalCharsToWrite(self):
                print "Calculating total number of characters...."
                count = 0
                for seq in self.seqs:
                        count += 4  # Line breaks after each line
                        count += len(seq[self.ID_INDEX])
                        count += len(seq[self.SEQ_INDEX])
                        count += len(seq[self.STRAND_INDEX])
                        count += len(seq[self.QSCORE_INDEX])

                return count


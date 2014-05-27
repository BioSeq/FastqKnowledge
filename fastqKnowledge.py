#!/usr/bin/env python

#
# fastqKnowledge.py
# Author: Philip Braunstein
#
# Date Created: May 14, 2014
# Last Modified: May 18, 2014
#
# This python file contains a main that takes in a FASTQ file and allows
# the user to query the FastqReporter data model.
#

from sys import argv
from sys import exit

from fastqReporter import FastqReporter

# CONSTANTS
EXTENSION = "fastq"

def main():
        checkArgs()

        sequences = FastqReporter(argv[1])  # Instantiate FastqReporter

        runLoop(sequences)
      


# Checks to make sure appropriate arguments are passed to program
def checkArgs():
        if len(argv) != 2:
                usage()
        if argv[1].split(".")[-1] != EXTENSION:
                print "Please provide a FASTQ file on the command line"
                usage()


# Prints correct usage and exits non-zero
def usage():
        print "USAGE:", argv[0], "FASTQ_FILE.fastq"
        exit(1)


# Runs program and and handles UI until x passed in by user.
def runLoop(sequences):
        running = True
        options()
        while (running):
                choice = raw_input("What would you like to do?\n")
                choice =  choice.strip()  # just to be safe
                choice = choice.lower()  # So users can enter caps or not
                print

                if choice == 'a':
                        sequences.printAvgQualScore()
                elif choice == 'b':
                        eqChar = getASCII() 
                        if eqChar == None:  # Signal that input wasn't valid
                                continue
                        print eqChar
                        print
                elif choice == 'c':
                        eqNum = getNUM()
                        if eqNum == None:
                                continue
                        print eqNum
                        print
                elif choice == 'd':
                        detailedHelp()
                elif choice == 'g':
                        sequences.printGC_AT()
                elif choice == 'l':
                        sequences.printAvgReadLength()
                elif choice == 'm':
                        seq = getSeq()
                        if seq == None:  # Invalid sequence pased in
                                continue
                        sequences.printNumMatches(seq)
                elif choice == 'n':
                        sequences.printNumReads()
                elif choice == 'o':
                        options()
                elif choice == 'p':
                        sequences.printTotalNumPages()
                elif choice == 'q':
                        qScore = getQScore()
                        if qScore == None:  # Signal that it was invalid qscre
                                continue
                        sequences.printNumQualSeqs(qScore)
                elif choice == 'r':
                        sequences.printRandomSeq()
                elif choice == 'u':
                        sequences.printNumNucs()
                elif choice == 'x':
                        print "Goodbye"
                        running = False
                else:
                        print "Unknown option:" + choice
                        print


# Propmts the user for a number to turn into ascii.
# Returns None if character is out of range. This 
# function also handles the error printing to explain
# out of range input.
def getASCII():
        try:
                num = int(raw_input("Which number would you like to convert" +\
                                        " to ASCII?\n"))
        except ValueError:
                print "Please enter a number for this functionality."
                print
                return None

        if num < 33 or num > 126:
                print "Number must be between 33 and 126"
                print
                return None

        return unichr(num)

# Function is capitalized to match style of getASCII
# Prints error message and returns None if input
# is invalid. Otherwise, returns the number ascii code
# for the character the user inputs.
def getNUM():
        chr = raw_input("Which ASCII character would you like to convert" +\
                                " to its number code?\n")

        chr = chr.strip()  # to be safe
        
        if len(chr) != 1:  # invalid input
                print "Please enter a single character"
                print
                return None

        return ord(chr)

# Prints detailed explanation of what this program does
def detailedHelp():
        print "Detailed explanation of the options:"
        print "a -- Finds the average quality score for every nucleotide" +\
                " in every read."
        print "b -- Converts a number to the ASCII character it codes for." +\
                " The number inputted must be between 33 and 126."
        print "c -- Converts an ASCII character to the number that codes" +\
                " it. The character entered by the user cannot be" +\
                " whitespace, and it must only be one character long."
        print "d -- Prints this message."
        print "g -- Get the G/C A/T content for all of the reads. This" +\
                " feature does not include 'N' nucleotides in its" +\
                " calculation."
        print "l -- Get the average number of nucleotides in each read in" +\
                " the sample."
        print "m -- Prints the number of reads that have an exact match of" +\
                " a nucleotide sequence entered by the user. Valid char" +\
                "acters are 'A', 'C', 'G', 'T', 'N'. 'N' acts as a wildcard" +\
                "  and matches any nucleotide."
        print "n -- Prints the total number of reads in the sample."
        print "o -- Prints the options list."
        print "p -- Prints the number of single-sided pages it would take" +\
                " to the whole file on 8.5 x 11 paper in 12-point font."
        print "q -- Prints the number of quality scores at or above a" +\
                " cuttoff that the user enters. Quality scores must be" +\
                " between ! (33) and ~ (126)"
        print "r -- Randomly select and print one reads exactly as it" +\
                " appeared in the FASTQ file."
        print "u -- Print the total number of nucleotides of each type in" +\
                " the sample."

        print "x -- Exit the program."
        print


# User options  how to use system
def options():
        print "Options:"
        print "a -- Average quality score of all reads"
        print "b -- Convert a number to ascii character"
        print "c -- Convert an ascii character to a number"
        print "d -- Detailed help"
        print "g -- G/C and A/T content of all reads"
        print "l -- Average read length of all reads"
        print "m -- Number of reads that match an inputted sequence"
        print "n -- Total number of reads"
        print "o -- options"
        print "p -- Number of pages this file would be in 12-point font"
        print "q -- Number of quality scores at or above an inputted cutoff"
        print "r -- Randomly select and print sequence from this file"
        print "u -- Total number of each nucleotide in the file" 
        print "x -- Exit"
        print



# Gets a nucleotides sequence from the user.
# This function controls for capitalization, and returns None
# if any of the characters in the sequence are not A, C, G, T, or N.
# Otherwise, returns the sequence.
def getSeq():
        seq = raw_input("What nucleotide sequence would you like to search" +\
                        " the reads for?\n")

        seq = seq.upper()  # Contol for capitalization
        if validSeq(seq):
                return seq
        else:  # Invalid characters in seq
                print "Invalid nucleotide sequence."
                print
                return None

# Returns True if a nucleotide sequence (has only A, C, T, G, or N)
# and returns FAlse if it contains any other characters.
# PRECONDITION: Sequence passed in must be in all caps.
def validSeq(seq):
        nucs = ['A', 'C', 'G', 'T', 'N']

        for nuc in seq:
                if nuc not in nucs:
                        return False
        
        return True  # If got to this point, no invalid nucleotides


# Returns a QScores that the user enters. It returns None
# if it is an invalid qscore as well as prints an error message.
def getQScore():
        qScore = raw_input("What quality score cutoff (in ASCII) would" +\
                        " you like to use?\n")

        qScore = qScore.strip()

        if len(qScore) != 1:
                print "Quality Score must be one character long."
                print
                return None

        return qScore




if __name__ == '__main__':
        main()

#!/usr/local/bin/python

import sys
from sys import argv 
import os
import subprocess
import ahocorasick
import timeit

script, original, corrected, dict = argv 


def makeAutomaton(filename):
	"""Build an Aho-Corasick automaton from a dictionary file and return
	it. The lines in the dictionary file must consist of a key and a
	value separated by a space."""
	global automaton
	print "Making automaton..."
	automaton = ahocorasick.Automaton()
	with open(filename) as f:
		for line in f:
			key, value = line.rstrip().split(" ", 1)
			automaton.add_word(key, (key, value))
	automaton.make_automaton()
	

def headerreplace(automaton, original, corrected):
	"""Apply an Aho-Corasick automaton to an input file, replacing the
	first occurrence of a key in each line with the corresponding
	value, and writing the result to the output file."""
	print "Performing search and replace"
	with open(original, 'r') as infile, open(corrected, 'w') as outfile:
		for line in infile:
			for end, (key, value) in automaton.iter(line):
				line = line[:end - len(key) + 1] + value + line[end + 1:]
				break # At most one replacement per line
			outfile.write(line)

def wrapper(func, *args, **kwargs):
	def wrapped():
		return func(*args, **kwargs)
	return wrapped

if __name__ == '__main__':
	wrapped = wrapper(makeAutomaton, dict)
	time = timeit.timeit(wrapped, number = 1)
	print "Automaton made in", time, "s."
	wrapped = wrapper(headerreplace, automaton, original, corrected)
	time2 = timeit.timeit(wrapped, number = 1)
	print "Finished. Search and replace finished in", time2, "s."

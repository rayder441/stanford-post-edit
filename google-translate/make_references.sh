#!/usr/bin/env bash
#
# Corresponds to the order in which source documents
# appear in source.en.txt 
#
cat Flag_of_Japan_Wikipedia.ar Schizophrenia_Wikipedia.ar Infinite_monkey_theorem_Wikipedia.ar 1896_Summer_Olympics_Wikipedia.ar > ar.ref

cat Flag_of_Japan_Wikipedia.de Schizophrenia_Wikipedia.de Infinite_monkey_theorem_Wikipedia.de 1896_Summer_Olympics_Wikipedia.de > de.ref

cat Flag_of_Japan_Wikipedia.fr Schizophrenia_Wikipedia.fr Infinite_monkey_theorem_Wikipedia.fr 1896_Summer_Olympics_Wikipedia.fr > fr.ref

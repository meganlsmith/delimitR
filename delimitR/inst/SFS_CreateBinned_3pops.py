#!/usr/bin/env python

import os
from sys import argv
from sortedcontainers import SortedDict
from collections import Counter
from collections import OrderedDict

script,Traits,file,Threshold,classes,outfile = argv #arguments the user must supply

classes = int(classes)
Threshold=int(Threshold) # make sure threshold is an integer

def setvalues(self, values):
    """
    You can pass in a list of values, which will replace the
    current list. The value list must be the same len as the OrderedDict.
    
    (Or a ``ValueError`` is raised.)
    
    >>> d = OrderedDict(((1, 3), (3, 2), (2, 1)))
    >>> d.setvalues((1, 2, 3))
    >>> d
    OrderedDict([(1, 1), (3, 2), (2, 3)])
    >>> d.setvalues([6])
    Traceback (most recent call last):
    ValueError: Value list is not the same length as the OrderedDict.
    """
    if len(values) != len(self):
        # FIXME: correct error to raise?
        raise ValueError('Value list is not the same length as the '
            'OrderedDict.')
    self.update(zip(self, values))

def pop_association(Traits):
    """Sets the individuals to their respective populations. 
    Also returns the sample counts per population."""
#    print 'Processing traits file...'
    with open(Traits, 'r') as traits:
        Pops = OrderedDict() # initializes a dictionary
        Pop_counts = SortedDict() # initializes a dictionary
        next(traits) # iterates to next item in file, skipping the header
        for line in traits: # for each line
            line = line.strip().split() # strip the line based on the tab and split it into two parts (allele and population)
            allele = '%r' % line[0]
#            Pops[line[0]] = line[1] # create a key in the Pops dictionary named for the allele and assign to it the name of the population to which the allele belongs
            Pops[allele] = line[1]
            if line[1] in Pop_counts: # check to see if the population is already present in the Pop_counts dictionary as a key
                Pop_counts[line[1]] += 1 # if it is, increase the count of alleles for that population
            else: # if it isn't
                Pop_counts[line[1]] = 1 # set the count equal to one
            #print Pops
        return Pops, Pop_counts # return both of these dictionaries for use in other functions

def Empty_AFS(Pops,Pop_counts):
#    print "Creating an empty dictionary for storing the AFS..."
    AFS_Empty = OrderedDict()
    npops = len(Pop_counts)
    for i in range(0,len(Pop_counts)):
         if i == 0:
             count = Pop_counts.items()[i][1] *Threshold/100
             Max_1 = count
         elif i == 1:
             count = Pop_counts.items()[i][1]  *Threshold/100
             Max_2 = count
         elif i == 2:
             count = Pop_counts.items()[i][1]  *Threshold/100
             Max_3 = count
    Max_1 = int(Max_1)
    Max_2 = int(Max_2)
    Max_3 = int(Max_3)

    s = 0
    if npops==3:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    bin = '%r_%r_%r' % (q,r,s)
#                    print bin
                    AFS_Empty.update({bin:0})
#                   AFS_Empty[bin]+=1
    if npops==4:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        bin = '%r_%r_%r_%r' % (q,r,s,t)
#                        print bin
                        AFS_Empty.update({bin:0})
#        print AFS_Empty
#        print len(AFS_Empty)
    if npops==5:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        for u in range(Max_5+1):
                            bin = '%r_%r_%r_%r_%r' % (q,r,s,t,u)
#                           print bin
                            AFS_Empty.update({bin:0})
#                           AFS_Empty[bin]+=1
    if npops==6:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        for u in range(Max_5+1):
                            for v in range(Max_6+1):
                                bin = '%r_%r_%r_%r_%r_%r' % (q,r,s,t,u,v)
#                               print bin
                                AFS_Empty.update({bin:0})
#                               AFS_Empty[bin]+=1
    if npops==7:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        for u in range(Max_5+1):
                            for v in range(Max_6+1):
                                for w in range(Max_7+1):
                                    bin = '%r_%r_%r_%r_%r_%r_%r' % (q,r,s,t,u,v,w)
#                                   print bin
                                    AFS_Empty.update({bin:0})
#                                   AFS_Empty[bin]+=1
    if npops==8:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        for u in range(Max_5+1):
                            for v in range(Max_6+1):
                                for w in range(Max_7+1):
                                    for x in range(Max_8+1):
                                        bin = '%r_%r_%r_%r_%r_%r_%r_%r' % (q,r,s,t,u,v,w,x)
#                                       print bin
                                        AFS_Empty.update({bin:0})
#                                       AFS_Empty[bin]+=1
    if npops==9:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        for u in range(Max_5+1):
                            for v in range(Max_6+1):
                                for w in range(Max_7+1):
                                    for x in range(Max_8+1):
                                        for y in range(Max_9+1):
                                            bin = '%r_%r_%r_%r_%r_%r_%r_%r_%r' % (q,r,s,t,u,v,w,x,y)
#                                           print bin
                                            AFS_Empty.update({bin:0})
#                                           AFS_Empty[bin]+=1
    if npops==10:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        for u in range(Max_5+1):
                            for v in range(Max_6+1):
                                for w in range(Max_7+1):
                                    for x in range(Max_8+1):
                                        for y in range(Max_9+1):
                                            for z in range(Max_10+1):
                                                bin = '%r_%r_%r_%r_%r_%r_%r_%r_%r_%r' % (q,r,s,t,u,v,w,x,y,z)
#                                               print bin
                                                AFS_Empty.update({bin:0})
#                                               AFS_Empty[bin]+=1

#    print len(AFS_Empty)
    return AFS_Empty

def Pooled_AFS(Pops, Pop_counts, nclasses): 
#    print "Making an empty dictionary to hold the pooled AFS..."
    nbins = nclasses ** len(Pop_counts)
    AFS_Pooled = OrderedDict()
    if len(Pop_counts) == 3: 
        for q in range(1,nclasses+1):
            for r in range(1,nclasses+1):
                for s in range(1,nclasses+1):
                    bin = '%r_%r_%r' % (q,r,s)
                    AFS_Pooled.update({bin:0})
    return AFS_Pooled

def Set_Thresholds(Pops,Pop_counts):
#    print Pop_counts['one']
#    print "Setting the thresholds..."
    threshold = 1/float(classes)
    ThreshCounter = SortedDict()
    for pop in range(0,len(Pop_counts)):
        count = Pop_counts.items()[pop][1]
        newcount = count * Threshold/100
#        print '%s is %s percent of %s' % (newcount, Threshold, count)
        dictname = 'Thresh_%s' % Pop_counts.items()[pop][0]
        for i in range(1,classes+1):
            dictname='Thresh_%s_%s' % (Pop_counts.items()[pop][0],i)
            ThreshCounter[dictname] = int(newcount*threshold*i)
    return ThreshCounter

            
def GetPopCounts(AFS,ThreshCounter,Pooled_AFS):
#    print "Pooling my allele counts..."
    Thresholds_pop1 = []
    Thresholds_pop2 = []
    Thresholds_pop3 = []
    for i in range(1,classes+1): 
        lookup_pop1 = 'Thresh_%s_%s' % (Pop_counts.keys()[0], i)
        Thresholds_pop1.append(ThreshCounter[lookup_pop1])
        lookup_pop2 = 'Thresh_%s_%s' % (Pop_counts.keys()[1], i)
        Thresholds_pop2.append(ThreshCounter[lookup_pop2])
        lookup_pop3 = 'Thresh_%s_%s' % (Pop_counts.keys()[2], i)
        Thresholds_pop3.append(ThreshCounter[lookup_pop3])
#        lookup_pop4 = 'Thresh_%s_%s' % (Pop_counts.keys()[3], i)
#        Thresholds_pop4.append(ThreshCounter[lookup_pop4])
    listofkeys = list(AFS.keys())
    for i in range(0, len(AFS)):
        key = listofkeys[i]
        key_split = key.split('_')
        population_1 = key_split[0]
        population_2 = key_split[1]
        population_3 = key_split[2]
#        population_4 = key_split[3]

        if int(population_1) <= Thresholds_pop1[0]:
            pop_1 = 1
        elif classes > 1 and int(population_1) <= Thresholds_pop1[1]:
            pop_1 = 2
        elif classes > 2 and int(population_1) <= Thresholds_pop1[2]:
            pop_1 = 3
        elif classes > 3 and int(population_1) <= Thresholds_pop1[3]:
            pop_1 = 4
        elif classes > 4 and int(population_1) <= Thresholds_pop1[4]:
            pop_1 = 5
        elif classes > 5 and int(population_1) <= Thresholds_pop1[5]:
            pop_1 = 6
        elif classes > 6 and int(population_1) <= Thresholds_pop1[6]:
            pop_1 = 7
        elif classes > 7 and int(population_1) <= Thresholds_pop1[7]:
            pop_1 = 8
        elif classes > 8 and int(population_1) <= Thresholds_pop1[8]:
            pop_1 = 9
        elif classes > 9 and int(population_1) <= Thresholds_pop1[9]:
            pop_1 = 10

        if int(population_2) <= Thresholds_pop2[0]:
            pop_2 = 1
        elif classes > 1 and int(population_2) <= Thresholds_pop2[1]:
            pop_2 = 2
        elif classes > 2 and int(population_2) <= Thresholds_pop2[2]:
            pop_2 = 3
        elif classes > 3 and int(population_2) <= Thresholds_pop2[3]:
            pop_2 = 4
        elif classes > 4 and int(population_2) <= Thresholds_pop2[4]:
            pop_2 = 5
        elif classes > 5 and int(population_2) <= Thresholds_pop2[5]:
            pop_2 = 6
        elif classes > 6 and int(population_2) <= Thresholds_pop2[6]:
            pop_2 = 7
        elif classes > 7 and int(population_2) <= Thresholds_pop2[7]:
            pop_2 = 8
        elif classes > 8 and int(population_2) <= Thresholds_pop2[8]:
            pop_2 = 9
        elif classes > 9 and int(population_2) <= Thresholds_pop2[9]:
            pop_2 = 10

        if int(population_3) <= Thresholds_pop3[0]:
            pop_3 = 1
        elif classes > 1 and int(population_3) <= Thresholds_pop3[1]:
            pop_3 = 2
        elif classes > 2 and int(population_3) <= Thresholds_pop3[2]:
            pop_3 = 3
        elif classes > 3 and int(population_3) <= Thresholds_pop3[3]:
            pop_3 = 4
        elif classes > 4 and int(population_3) <= Thresholds_pop3[4]:
            pop_3 = 5
        elif classes > 5 and int(population_3) <= Thresholds_pop3[5]:
            pop_3 = 6
        elif classes > 6 and int(population_3) <= Thresholds_pop3[6]:
            pop_3 = 7
        elif classes > 7 and int(population_3) <= Thresholds_pop3[7]:
            pop_3 = 8
        elif classes > 8 and int(population_3) <= Thresholds_pop3[8]:
            pop_3 = 9
        elif classes > 9 and int(population_3) <= Thresholds_pop3[9]:
            pop_3 = 10

#        if int(population_4) <= Thresholds_pop4[0]:
#            pop_4 = 1
#        elif classes > 1 and int(population_4) <= Thresholds_pop4[1]:
#            pop_4 = 2
#        elif classes > 2 and int(population_4) <= Thresholds_pop4[2]:
#            pop_4 = 3
#        elif classes > 3 and int(population_4) <= Thresholds_pop4[3]:
#            pop_4 = 4
#        elif classes > 4 and int(population_4) <= Thresholds_pop4[4]:
#            pop_4 = 5
#        elif classes > 5 and int(population_4) <= Thresholds_pop4[5]:
#            pop_4 = 6
#        elif classes > 6 and int(population_4) <= Thresholds_pop4[6]:
#            pop_4 = 7
#        elif classes > 7 and int(population_4) <= Thresholds_pop4[7]:
#            pop_4 = 8
#        elif classes > 8 and int(population_4) <= Thresholds_pop4[8]:
#            pop_4 = 9
#        elif classes > 9 and int(population_4) <= Thresholds_pop4[9]:
#            pop_4 = 10

        keytopop = '%s_%s_%s' % (pop_1,pop_2,pop_3)
        Pooled_AFS[keytopop] += float(AFS[key])
    return Pooled_AFS

def print_AFS(Full_AFS):
#    print "Writing the AFS..."
    for bin in Full_AFS:
        with open(outfile, 'a') as f:
            f.write(str(Full_AFS[bin]))
            f.write('\t')
    with open(outfile, 'a') as f:
        f.write('\n')

def Populate_AFS(file):
#    print "Populating the AFS..."
    with open(file, 'r') as Infile:
        for line in Infile:
            Full_AFS = Empty_AFS(Pops,Pop_counts)
            line = line.strip().split()
            counter = 0
            setvalues(Full_AFS,line)
            tPooled_AFS = Pooled_AFS(Pops, Pop_counts, classes)
            Final_AFS = GetPopCounts(Full_AFS,ThreshCounter,tPooled_AFS)
            print_AFS(Final_AFS)


Pops, Pop_counts = pop_association(Traits)
ThreshCounter = Set_Thresholds(Pops,Pop_counts)
Populate_AFS(file)
#





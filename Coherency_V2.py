# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 14:26:05 2020

@author: paulv
"""


import requests 
import os
import time

os.chdir("D:/Internship AMOLF/Python/random_sample_hd" )

def get_file(filename):
    # Gets the file from the working directory
    file= open(filename, "r" )
    genes=[]
    for line in file:
        line= line.rstrip()
        genes.append(line)
    return genes






def API_request(genes, th_sc):
    th_sc=int(th_sc)
    #Receives response from string-db for the given genes from the modules
    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"
    request_url = "/".join([string_api_url, output_format, method])
    my_genes=genes
    params = {
    "identifiers" : "%0d".join(my_genes), # your protein
    "species" : 6239, # species NCBI identifier 
    "caller_identity" : "Paul van Lent", # your app name
    "required_score":th_sc

    }
    response = requests.post(request_url, data=params)
    return response.text



#output is  table with 13 columns =
#  score   nscore  fscore  pscore  ascore  escore  dscore  tscore
#Combined score, fscore=gene fusion score, gene neighborhood score, 
#pscore phylogenetic score, ascore coexpression score,
# escore experimental score, dscore database score, tscore textming score


def extract_interactions(interactions):
    #This function extracts the proteins involved in an interaction
    prot_interactions=[]
    for line in interactions.strip().split("\n"):
        l = line.strip().split("\t")
        m= (l[2],l[3])
        prot_interactions.append(m)
    prot_interactions= set(prot_interactions)
    return prot_interactions



def coherency_score(protein_interactions, genes):
    #Calculates the average node degree and a score (the number of edges/ total number of edges)
    N=len(genes)
    edges= (N * (N-1))/2
    score= len(protein_interactions)/edges
    avg_degree= (len(protein_interactions)*2)/N
    
    return score, avg_degree
    


def main(th_sc):
    th_sc= int(th_sc)
    dir_list= os.listdir()
    scores=[]
    avg_degrees=[]
    no_known_interaction=0
    for i in dir_list:
        genes=get_file(i)
        interactions= API_request(genes, th_sc)
        if len(interactions)>0:
            print(i)
            prot= extract_interactions(interactions)
            score, avg_degree=coherency_score(prot, genes)
            score= i+"\t"+str(score)
            avg_degree=i+"\t"+str(avg_degree)
            scores.append(score)
            avg_degrees.append(avg_degree)
            time.sleep(1)
        else:
            time.sleep(1)
            no_known_interaction+=1
        edge_file="../edge_score_rs_hd%d.txt" %th_sc
        avg_file="../avg_degree_rs_hd%d.txt" %th_sc
        file= open(avg_file,"w")
        file2=open(edge_file,"w")
        for i in range(0,len(avg_degrees)):
            file.write(str(avg_degrees[i])+"\n")
            file2.write(str(scores[i])+"\n")
        file.close()
        file2.close()

    

if __name__ == '__main__':
    main(150)




import os
import numpy as np
import pandas as pd
import pickle
#from sklearn.externals import joblib
import sys
import joblib

def aminoAcidEmbedding(embeddingFile):
    f=open(embeddingFile,"r")
    lines=f.readlines()
    f.close()
    
    #AA_Emb is a dictionary of (key, value) where key is the amino acid and value is the corresponding Pfam amino acid embedding
    AA_Emb={}
    for line in lines[1:]:
        temp=line[:-1].split()
        AA_Emb[temp[0]]=temp[1:]
    
    return AA_Emb
        

def fastaToFeature(fastaFile):
    f=open(fastaFile,"r")
    lines=f.readlines()
    f.close()
    
    
    #input_sequence is a dictionary of (key, value) where key is the proteinID and value is the amino acid sequence 
    input_sequences={}
    
    fastaSequence=""
    temp=lines[0]
    temp=temp.replace(">sp|","").replace(">","")
    proteinID=temp[:temp.find("|")]
    for line in lines:
        if line.find(">")<0:
            if line[-1]=="\n":
                fastaSequence+=line[:-1]
            else:
                fastaSequence+=line
        else:
            input_sequences[proteinID]=fastaSequence
            temp=line
            temp=temp.replace(">sp|","").replace(">","")
            proteinID=temp[:temp.find("|")]
            fastaSequence=""
    input_sequences[proteinID]=fastaSequence
    
    #segment each amino acid sequence into segment of 15
    #input_segments is a dictionary of (key, value) where key is the proteinID and value is the segments of 15 AA from the amino acid sequences 
    input_segments={}
    for proteinID in input_sequences.keys():
        sequence=input_sequences.get(proteinID)
        sequence="OOOOOOO"+sequence #padding
        sequence=sequence+"OOOOOOO" #padding
            
        segments=[]
        for i in range(0,len(sequence)-14):
            segments.append(sequence[i:i+15])
        input_segments[proteinID]=segments
        
    #create feature vectors
    #input_feature is a dictionary of (key, value) where key is the proteinID and value is the embedding features created from Pfam_emb of segments of 15 amino acids
    aaEmb=aminoAcidEmbedding("Embedding vectors/Pfam.vec")
    
    if not(os.path.exists("tmp")):
        os.mkdir("tmp")
    for proteinID in input_segments.keys():
        f=open("tmp/"+proteinID+".csv","w")
        
        segments=input_segments.get(proteinID)
        for segment in segments:
            for amino_acid in segment:
                emb_values=aaEmb.get(amino_acid)
                for emb_value in emb_values:
                    f.write(str(emb_value)+",")
            f.write("\n")
        f.close()
    return input_sequences


def labelToOneHot(label):# 0--> [1 0], 1 --> [0 1]
    label = label.reshape(len(label), 1)
    label = np.append(label, label, axis = 1)
    label[:,0] = label[:,0] == 0;
    return label

def predict(inputFile, outputFile="Result.txt"):
    input_sequences=fastaToFeature(inputFile)
    
    #loading the model
    try:
        classifier=joblib.load("Model/PfamVecSize8ModelWithRF.sav")
    except (IOError, pickle.UnpicklingError, AssertionError):
        print(pickle.UnpicklingError)
        return True
    
    threshold=0.365 #predefined threshold
    #loop through each protein feature file and write results into output file
    f=open(outputFile,"w")
    f.write("No. of sequence = "+str(len(os.listdir("tmp")))+"\n")
    f.write("Predefined threshold: "+str(threshold)+"\n")
    f.write("------------------------------------------------\n")
    f.write("Position     Residue	  Score	     Prediction\n")
    f.write("------------------------------------------------\n")
    
    
    
    
    for protein_feature_file in os.listdir("tmp"):
        proteinID=protein_feature_file[:-4]
        sequence=input_sequences.get(proteinID)
        length=len(sequence)
        f.write(">"+proteinID+"    Length = "+str(length)+"\n\n")
        f.write(sequence+"\n\n")
        
        dataset = pd.read_csv("tmp/"+protein_feature_file, header=None)
        X_test = dataset.iloc[:, 0:-1].values
        y_pred = classifier.predict_proba(X_test)
        for i in range(len(y_pred)-2):
            if sequence[i]=="N":
                if y_pred[i][1]>=threshold:
                    f.write("{:5.0f}".format(i+1)+"    "+sequence[i:i+3]+"    "+"{:10.4f}".format(y_pred[i][1])+"    Potential glycosylated\n")
                else:
                    f.write("{:5.0f}".format(i+1)+"    "+sequence[i:i+3]+"    "+"{:10.4f}".format(y_pred[i][1])+"    Non-glycosylated\n")
        
        f.write("***********************************\n\n")

	


inputFile= sys.argv[1]
outputFile="Result.txt" 
predict(inputFile, outputFile="Result.txt")


#delete temporary files    
for filename in os.listdir("tmp"):
    os.remove("tmp\\"+filename)


print("Thank you for using NIonPred!!! Please check the prediction results in Result.txt file")



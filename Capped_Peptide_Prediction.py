"""This file is for predicting where capped peptides might by in the secretome."""
secretome = open('/workspaces/Capped-Peptides/Human_Uniprot_Secreted.fasta', 'r') ## Human
#secretome = open('/workspaces/Capped-Peptides/Mouse_Uniprot_Secreted.fasta', 'r') ## Mouse
##Writing the file for output information, predicted peptides and where they come from
output = open('/workspaces/Capped-Peptides/Human_Predicted_Capped_Peptides.csv', 'w') ##Human
#output = open('/workspaces/Capped-Peptides/Mouse_Predicted_Capped_Peptides.csv', 'w') ##Mouse

## Going through the secretome file and making a dictionary of all secreted proteins
index = 0
secreted_genes = {} #storing secreted genes in this dictionary
gene_num = 1
for line in secretome:
    index +=1
    if line[0] == ">":
        if index > 1:
            #print('true')
            secreted_genes[accession] = [name, sequence] #storing name or fasta header and sequence
            gene_num +=1
        line = line.strip()
        end = line.rfind('|')
        name = str(line[12:])
        name = name.replace(',','')
        accession = str(line[4:end])
        sequence = ''
    else:
        line = line.strip()
        sequence += str(line)
secreted_genes[accession] = [name, sequence]


pep_list = [] #list of all peptides as we search through and predict them
pep_num = 0 #count of peptides predicted
gene_list_pep = [] #genes peptides are predicted from

#Now searching through our dictionary of secreted genes for capped peptide predictions
for gene in secreted_genes:
    seq = secreted_genes[gene][1] #sequence of protein
    for position in range(len(seq) - 2): #scanning through looking for the tripeptide motif for amidation
        if (seq[position:position+3] == 'GRR')or (seq[position:position+3] =='GKR') and (position > 1): #checking if we find the amidation motif
            #Now need to look if there is a glutamine upstream within 20 amino acids
            if position < 20 : #condition where there is less than 20 amino acids upstream
                start = 0
                fragment = seq[:position]
            else:
                start = position - 20
                fragment = seq[start:position]

            q_pos = fragment.rfind('Q') #using rfind because looking for glutamine closes to the amidaiton motif
            if q_pos != -1: #the rfind would give -1 if there is no glutamine
                pot_hormone = fragment[q_pos:] #extracting out the sequence of the potential peptide
                q_pos_in_seq = start + q_pos #amino acid position in the starting sequence
                #extracting out the amino acids around both capping sites for future motif analysis
                if q_pos_in_seq > 5:
                    upstream = seq[q_pos_in_seq-4:q_pos_in_seq+5]
                else:
                    upstream = seq[:q_pos_in_seq + 5]
                if len(seq[position+3:]) > 5 :
                    downstream = seq[position-4:position+7]
                else:
                    downstream = seq[position-4:]
                #Now writing a line in the csv if the peptide is longer than a dipeptide
                if len(pot_hormone) > 2:
                    if pot_hormone not in pep_list: #not storing repeats of the same peptide
                        pep_list += [str(pot_hormone)]
                        pep_num +=1
                        num_pep = pep_num
                        if gene not in gene_list_pep:
                            gene_list_pep += [gene] 
                        if num_pep < 10:
                            num_pep = '00' + str(num_pep)
                        elif num_pep <100:
                            num_pep = '0' + str(num_pep)
                        pep_name = 'mouse CAP ' + str(num_pep) #temporary naming
                        output.write(str(pep_name) +',' + str(pot_hormone) + ',' + str(gene) + ',' + str(secreted_genes[gene][0]) + ',' + str(q_pos_in_seq) + ',' + upstream + ',' + downstream + '\n')


output.close()
secretome.close()

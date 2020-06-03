import os
import random


species_one_CDS = []
species_two_CDS = []

#this path represents the location where all the FASTA files
#are located. Then, we list all the files in the directory
#using the listdir function. 
fasta_path = 'C:\\Users\\maddy\\OneDrive\\UW Classes\\Sophmore Year\\Quarter 3\\CSS 383\\Project 2\\FASTA Files\\'

# allSpeciesFiles is an array that stores all 17 Acipenser FASTA files
allSpeciesFiles = os.listdir(fasta_path)

annot_path = 'C:\\Users\\maddy\\OneDrive\\UW Classes\\Sophmore Year\\Quarter 3\\CSS 383\\Project 2\\All 17 Acipenser Annotations\\'

# annotFiles is an array that saves all 17 annotated Acipenser files with CDS
annotFiles = os.listdir(annot_path)

# open a new similarity file which will hold all the similarity percentages for all the species compared. This is so that we don't have to check each species-combo percentage and is for ease of viewing to compare results. 
# change this to the path on your own computer's location where you would like the file to be saved.

#similarity_file = open("C:\\Users\\maddy\\OneDrive\\UW Classes\\Sophmore Year\\Quarter 3\\CSS 383\\Project 2\\similarity_results.txt", "w")


def calculate_CDS_vals(file_name, fasta_file_num, arr_to_fill):
  # open the annotation file
  file_to_parse = open(file_name, 'r')

  # go through each line in the file
  for line in file_to_parse:
    # convert the received line into an array to save all values
      line_arr = line.strip().split("\t")

      # we only want the 13 CDS, so if the line is not CDS, then don't add it to 
      # the array. 
      if("CDS" in line_arr):
          # since one of the CDS lines has a complement, so to get rid of it, we have
          # a replace here to remove excess characters. 
          line_arr[1] = line_arr[1].replace('complement(', '')
          line_arr[1] = line_arr[1].replace(')', '')
          # gets the start and stop index for the CDS so we can retrieve the substring
          # from the FASTA file.                
          start_stop_indexes = line_arr[1].split("..")

          #convert the start and stop to integers to use as indexes
          start = int(start_stop_indexes[0])
          stop = int(start_stop_indexes[1])

          # get the FASTA file and the particular substring with the above start and stop
          # indexes, and save to the array that saves all CDS in the same order everytime.
          fasta_file = open(fasta_path + allSpeciesFiles[fasta_file_num])
          content = fasta_file.read()

          # add the particular substring from the fasta file using the start and stop
          # indexes, into the array that will store all CDS values.
          arr_to_fill.append(content[start:stop])


def do_BLAST_calculation(result_name):
  # run through all 12 CDS for both species and compare each CDS for each species through BLAST and get comparison percentage for those sequences. 
  for k in range(13):
    # make a new file for the first and second species, so that the particular sequence
    # will be stored there and then blast can compare the 2 files, instead of having 
    # to compare strings. 
    file_one = open("C:\\Users\\maddy\\OneDrive\\UW Classes\\Sophmore Year\\Quarter 3\\CSS 383\\Project 2\\file_one.fa", "w")

    file_two = open("C:\\Users\\maddy\\OneDrive\\UW Classes\\Sophmore Year\\Quarter 3\\CSS 383\\Project 2\\file_two.fa", "w")

    # send the particular sequence to the fasta file so BLAST can compare them
    file_one.write(species_one_CDS[k])
    file_two.write(species_two_CDS[k])

    # have the global variable so that it runs on windows/linux
    os.putenv("BLASTDB_LMDB_MAP_SIZE", "1000000")

    # invoke the makeblastdb so that you can run commands through the .exe using BLAST
    firstCommand = "makeblastdb -in {subject_fasta_file} -input_type fasta -dbtype nucl"

    # specify the file name to use in the command statement using .format
    firstCommand = firstCommand.format(subject_fasta_file = species_one_CDS[k])

    # sends the first command, to invoke BLAST and to use its features, to the system
    os.system(firstCommand)

    # command that will be sent to the system to get output results for
    output = "blastn -query {query_fasta_file} -subject {subject_fasta_file} -task blastn -max_hsps 1 -outfmt 6 >> {results_name}"

    # put the appropriate values of the commands above 
    output = output.format(query_fasta_file = "file_one.fa", subject_fasta_file = "file_two.fa", results_name = result_name)

    # execute the statement above to get the actual results
    os.system(output)

    # empty the contents in the file so that every time the new CDS is compared between
    # species, and the previous sequences don't stay there. 
    open('file_one.fa', 'w').close()
    open('file_two.fa', 'w').close()

  # clear the array after each round, so that the sequences from previous comparisons
  # don't stay and are removed to only have 13 new CDS each time. 
  species_one_CDS.clear()
  species_two_CDS.clear()


def main():
  # systematically compare each species to every other species in the same order each time
  # to get all the comparisons done 

  # changed the directory to make sure we are putting the output files in the correct place. May be required based on particular systems the code in run on.
  # this is optional and might not be needed for your computer, but for mine, there were issues with printing to the correct location so i used it. If you get any errors with the line below, comment it so that there will be no errors when it runs. 
  os.chdir("C:\\Users\\maddy\\OneDrive\\UW Classes\\Sophmore Year\\Quarter 3\\CSS 383\\Project 2")

  count = 0

  # range for the j and k decides the number of files that are compared and which particular files, out of the 17 files, are compared 
  for j in range(0,5):           
    for k in range(0,5):
      # made an array to store all elements of the file name, make it easier to keep track when we do analysis of different comparisons.
      species_name_array = annotFiles[j].split("_") 
      species_one_name = ""  
      # if/else statement for including the genuine cases where 2 species breed and there one species x another species, then the size of array is always 7, so checks for that and includes those case, else prints the regular name of the species, removing .txt.
      if(len(species_name_array) == 7):
        species_one_name += (species_name_array[3].rstrip('.txt') + "_" + species_name_array[4].upper() + "_" + species_name_array[6].rstrip('.txt'))
      elif(len(species_name_array) == 4):
        species_one_name += (species_name_array[3].rstrip('.txt'))

      if(j == k):
        continue
      else:
        # if/else statement for including the genuine cases where 2 species breed and there one species x another species, then the size of array is always 7, so checks for that and includes those case, else prints the regular name of the species, removing .txt.
        species_two_array = annotFiles[k].split("_")
        species_two_name = ""
        if(len(species_two_array) == 7):
          species_two_name += (species_two_array[3].rstrip('.txt') + "_" + species_two_array[4].upper() + "_" + species_two_array[6].rstrip('.txt'))
        elif(len(species_two_array) == 4):
          species_two_name += (species_two_array[3].rstrip('.txt'))
      
      calculate_CDS_vals(annot_path + annotFiles[j], j, species_one_CDS)
      calculate_CDS_vals(annot_path + annotFiles[k], k, species_two_CDS)

      results_file = "{species_one}_{species_two}_results.txt"      
      results_file = results_file.format(species_one = species_one_name, species_two = species_two_name) 
      do_BLAST_calculation(results_file)
      count += 1

  print("Total comparisons: ", count)

main()